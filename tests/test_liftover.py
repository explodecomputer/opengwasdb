"""Tests for opengwasdb.build.liftover."""

from __future__ import annotations

import pytest

from opengwasdb.build.liftover import LiftoverFailureError, build_liftover_lookup


# Known good positions: hg19 1:100000 → hg38 1:100000 (identity in this region)
# hg19 1:1000000 → hg38 1:1064620
KNOWN_HG19_POSITIONS = [
    ("1", 100_000, "A", "G"),
    ("1", 1_000_000, "C", "T"),
]


def test_build_liftover_lookup_returns_dict_for_successful_variants():
    result = build_liftover_lookup(KNOWN_HG19_POSITIONS)

    assert len(result) == 2
    for key in KNOWN_HG19_POSITIONS:
        assert key in result


def test_build_liftover_lookup_output_alids_have_canonical_ordering():
    variants = [("1", 100_000, "G", "A")]  # G > A, so A is A1
    result = build_liftover_lookup(variants)

    assert len(result) == 1
    alid = list(result.values())[0]
    parts = alid.split(":")
    a1, a2 = parts[2], parts[3]
    assert a1 <= a2  # A1 is alphabetically first


def test_build_liftover_lookup_known_coordinate_maps_correctly():
    # hg19 1:100000 → hg38 1:100000 (identity at this known position)
    result = build_liftover_lookup([("1", 100_000, "A", "G")])
    alid = result[("1", 100_000, "A", "G")]
    chrom, pos, a1, a2 = alid.split(":")
    assert chrom == "1"
    assert int(pos) == 100_000  # this position maps identically


def test_build_liftover_lookup_accepts_chr_prefixed_input():
    variants = [("chr1", 100_000, "A", "G")]
    result = build_liftover_lookup(variants)
    assert len(result) == 1
    # Key in dict preserves original form
    assert ("chr1", 100_000, "A", "G") in result


def test_build_liftover_lookup_omits_failed_variants_below_threshold():
    # Positions that don't lift (e.g., in unsequenced gaps) are omitted.
    # hg19 1:200000 fails; 1:100000 succeeds.
    variants = [("1", 200_000, "A", "G"), ("1", 100_000, "C", "T")]
    result = build_liftover_lookup(variants, failure_threshold=0.6)

    assert ("1", 100_000, "C", "T") in result
    assert ("1", 200_000, "A", "G") not in result


def test_build_liftover_lookup_raises_when_failure_rate_exceeds_threshold():
    # Both positions in a known gap so both should fail.
    # hg19 1:200000 and 1:300000 both fail liftover.
    variants = [
        ("1", 200_000, "A", "G"),
        ("1", 300_000, "C", "T"),
    ]
    with pytest.raises(LiftoverFailureError, match="failure rate"):
        build_liftover_lookup(variants, failure_threshold=0.01)


def test_build_liftover_lookup_empty_input_returns_empty_dict():
    result = build_liftover_lookup([])
    assert result == {}


def test_build_liftover_lookup_output_chrom_is_bare():
    result = build_liftover_lookup([("1", 100_000, "A", "G")])
    alid = list(result.values())[0]
    chrom = alid.split(":")[0]
    assert not chrom.startswith("chr")
