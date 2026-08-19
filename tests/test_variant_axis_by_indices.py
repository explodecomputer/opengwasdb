"""`VariantAxis.by_indices()` adaptive random-access vs full-scan dispatch
(issue #104 follow-up).

A random-access `by_index()` seek and a sequential `all()` scan return the
same data; the only thing under test here is that `by_indices()` picks the
right strategy for the request size and that both strategies agree.
"""

from __future__ import annotations

from opengwasdb.variants import VariantAxis
from opengwasdb.variants.axis import _BULK_SCAN_FRACTION, write_variant_axis
from opengwasdb.variants.normalise import CanonicalVariant

N_VARIANTS = 1000  # threshold = N_VARIANTS * _BULK_SCAN_FRACTION = 10


def _build_axis(tmp_path):
    variants = [
        CanonicalVariant(chromosome="1", position=1000 + i, effect_allele="A", other_allele="G")
        for i in range(N_VARIANTS)
    ]
    write_variant_axis(tmp_path, variants, {})
    return VariantAxis(tmp_path)


def test_below_threshold_uses_random_access(tmp_path, monkeypatch):
    axis = _build_axis(tmp_path)
    try:
        small = list(range(5))
        assert len(small) <= N_VARIANTS * _BULK_SCAN_FRACTION

        calls = {"by_index": 0, "all": 0}
        real_by_index = axis.by_index

        def spy_by_index(index):
            calls["by_index"] += 1
            return real_by_index(index)

        monkeypatch.setattr(axis, "by_index", spy_by_index)
        monkeypatch.setattr(axis, "all", lambda: (_ for _ in ()).throw(
            AssertionError("full scan should not run below threshold")
        ))

        records = axis.by_indices(small)
        assert sorted(records) == small
        assert calls["by_index"] == len(small)
    finally:
        axis.close()


def test_above_threshold_uses_full_scan(tmp_path, monkeypatch):
    axis = _build_axis(tmp_path)
    try:
        large = list(range(50))
        assert len(large) > N_VARIANTS * _BULK_SCAN_FRACTION

        monkeypatch.setattr(axis, "by_index", lambda index: (_ for _ in ()).throw(
            AssertionError("random access should not run above threshold")
        ))

        records = axis.by_indices(large)
        assert sorted(records) == large
        for i in large:
            assert records[i].variant_index == i
            assert records[i].position == 1000 + i
    finally:
        axis.close()


def test_both_strategies_agree(tmp_path):
    axis = _build_axis(tmp_path)
    try:
        small = list(range(5))
        large = list(range(50))
        via_random_access = axis.by_indices(small)
        via_full_scan = axis.by_indices(large)
        for i in small:
            assert via_random_access[i] == via_full_scan[i]
    finally:
        axis.close()


def test_empty_indices_returns_empty_dict(tmp_path):
    axis = _build_axis(tmp_path)
    try:
        assert axis.by_indices([]) == {}
    finally:
        axis.close()
