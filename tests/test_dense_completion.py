"""Integration tests for Dense reference completion (ADR 0022, ADR 0023)."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pytest
import zarr

from opengwasdb.build.observed import build_dense_observed_from_sources
from opengwasdb.layouts.dense import complete as complete_module
from opengwasdb.layouts.dense.complete import (
    complete_dense_store,
    resume_dense_completion,
)
from opengwasdb.query import query_store
from opengwasdb.validation.validate import validate_store


def _assert_stores_identical(dst_a: Path, dst_b: Path) -> None:
    """Assert two completed stores hold byte-identical z/se/imputed arrays.

    The parity guarantee for completion is stronger than matching cell counts:
    the imputed z/se *values* (and the imputed/on_panel masks) must match too.
    z/se are stored float16, so exact equality (NaN-aware) is the right check.
    """
    root_a = zarr.open_group(str(dst_a / "data.zarr"), mode="r")
    root_b = zarr.open_group(str(dst_b / "data.zarr"), mode="r")
    for name in ("z", "se"):
        arr_a = root_a[name][:]
        arr_b = root_b[name][:]
        assert arr_a.shape == arr_b.shape, f"{name} shape differs"
        assert np.array_equal(arr_a, arr_b, equal_nan=True), f"{name} values differ"
    for name in ("imputed", "on_panel"):
        assert np.array_equal(root_a[name][:], root_b[name][:]), f"{name} mask differs"

SOURCE_HEADER = "\t".join(
    [
        "analysis_id", "phenotype_id", "phenotype_label", "analysis_label",
        "chromosome", "position", "effect_allele", "other_allele",
        "z", "se", "rsid", "stored_effect_scale",
    ]
)

# a1 observes 4 of the 5 chr1 block positions (900000/950000/1000000/1100000), leaving
# 1050000 as a genuine imputation target — enough points for poly_rescale's quality fit.
# a2 observes only chr1:1,000,000 → its block vector has 1 obs point → never imputed.
# chr1:1,200,000 A/G — observed by a1 only; OFF-PANEL (not in any LD block) → a2's
#                       n_missing_off_panel picks this up.
SOURCE_ROWS = [
    "a1\tp1\tHeight\tHeight primary\t1\t900000\tA\tG\t1.0\t0.15\trs0\tsd_units",
    "a1\tp1\tHeight\tHeight primary\t1\t950000\tA\tC\t1.8\t0.15\trs0b\tsd_units",
    "a1\tp1\tHeight\tHeight primary\t1\t1000000\tA\tG\t2.0\t0.15\trs1\tsd_units",
    "a2\tp2\tDisease\tDisease primary\t1\t1000000\tA\tG\t6.0\t0.20\trs1\tlog_or",
    "a1\tp1\tHeight\tHeight primary\t1\t1100000\tA\tC\t3.0\t0.20\trs2\tsd_units",
    "a1\tp1\tHeight\tHeight primary\t1\t1200000\tA\tG\t1.5\t0.30\trs3\tsd_units",
]


def _write_ld_block(
    block_dir: Path,
    block_name: str,
    snps: list[tuple[str, float, int]],
    seed: int = 0,
) -> None:
    """Write one flat-layout LD block: {block_name}.tsv + .unphased.vcor1.gz."""
    block_dir.mkdir(parents=True, exist_ok=True)
    n = len(snps)

    tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
    for alid, eaf, bp in snps:
        chrom, _pos, a1, a2 = alid.split(":")
        tsv_lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t{eaf}\t{bp}")
    (block_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    ld = A @ A.T + np.eye(n) * n * 0.1

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for row in ld:
            gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
    (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())


def _make_ld_panel(tmp_path: Path) -> Path:
    """Two blocks: chr1 (mixes observed + one new panel SNP) and chr2 (pure panel-only)."""
    root = tmp_path / "ld_panel"
    _write_ld_block(
        root / "EUR" / "1", "900000-1300000",
        [
            ("1:900000:A:G", 0.35, 900_000),     # observed by a1
            ("1:950000:A:C", 0.32, 950_000),     # observed by a1
            ("1:1000000:A:G", 0.30, 1_000_000),  # observed by a1, a2
            ("1:1050000:C:T", 0.40, 1_050_000),  # new — imputation target
            ("1:1100000:A:C", 0.45, 1_100_000),  # observed by a1
        ],
        seed=0,
    )
    _write_ld_block(
        root / "EUR" / "2", "1-500000",
        [
            ("2:100000:A:G", 0.20, 100_000),  # new — not observed by anyone
            ("2:200000:C:G", 0.25, 200_000),  # new — not observed by anyone
        ],
        seed=1,
    )
    return root


@pytest.fixture
def source_path(tmp_path: Path) -> Path:
    path = tmp_path / "associations.tsv"
    path.write_text(SOURCE_HEADER + "\n" + "\n".join(SOURCE_ROWS) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def observed_store(tmp_path: Path, source_path: Path) -> Path:
    out = tmp_path / "obs.opengwasdb"
    build_dense_observed_from_sources(
        [source_path], out,
        store_id="test", release_id="obs-v1", reference_assembly="GRCh38",
    )
    return out


@pytest.fixture
def ld_panel(tmp_path: Path) -> Path:
    return _make_ld_panel(tmp_path)


@pytest.fixture
def completed_store(tmp_path: Path, observed_store: Path, ld_panel: Path) -> Path:
    dst = tmp_path / "comp.opengwasdb"
    complete_dense_store(
        observed_store, dst, ld_panel,
        ancestry="EUR", min_cor=0.0, release_id="comp-v1",
    )
    return dst


class TestCompletionFiles:
    def test_creates_store_directory(self, completed_store):
        assert completed_store.exists()
        assert (completed_store / "manifest.json").exists()
        assert (completed_store / "variants.tsv.gz").exists()
        assert (completed_store / "index.sqlite").exists()

    def test_no_leftover_checkpoint_or_work_dir(self, tmp_path, completed_store):
        assert not (tmp_path / f".{completed_store.name}.checkpoint").exists()
        assert not (tmp_path / f".{completed_store.name}.tmp").exists()

    def test_manifest_completion_state(self, completed_store):
        import json
        m = json.loads((completed_store / "manifest.json").read_text())
        assert m["completion_state"] == "reference_completed"
        assert "completion" in m["provenance"]

    def test_variant_table_larger_than_source(self, completed_store, observed_store):
        from opengwasdb.variants.axis import VariantAxis
        obs_ax = VariantAxis(observed_store)
        comp_ax = VariantAxis(completed_store)
        # +2 (chr1 new panel SNP, off-panel chr1:1200000 already in source) +2 (chr2 new panel SNPs)
        assert comp_ax.n_variants == obs_ax.n_variants + 3
        obs_ax.close()
        comp_ax.close()

    def test_imputed_and_on_panel_arrays_present(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r")
        assert "imputed" in root
        assert "on_panel" in root
        assert root["imputed"].shape == root["z"].shape
        assert root["on_panel"].shape[0] == root["z"].shape[0]

    def test_off_panel_variant_never_imputed(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r")
        on_panel = root["on_panel"][:]
        imputed = root["imputed"][:]
        assert not np.any(imputed[on_panel == 0, :] == 1)

    def test_completion_quality_table_exists(self, completed_store):
        import sqlite3
        conn = sqlite3.connect(str(completed_store / "index.sqlite"))
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        conn.close()
        assert "completion_quality" in tables

    def test_analyses_table_has_n_missing_off_panel(self, completed_store):
        import sqlite3
        conn = sqlite3.connect(str(completed_store / "index.sqlite"))
        conn.row_factory = sqlite3.Row
        rows = {
            r["analysis_id"]: r["n_missing_off_panel"]
            for r in conn.execute("SELECT * FROM analyses")
        }
        conn.close()
        # a1 observed the off-panel chr1:1200000 variant; a2 never did.
        assert rows["a1"] == 0
        assert rows["a2"] == 1

    def test_overwrite_raises_without_flag(
        self, tmp_path, observed_store, ld_panel, completed_store
    ):
        with pytest.raises(FileExistsError):
            complete_dense_store(observed_store, completed_store, ld_panel, min_cor=0.0)


class TestRegionCap:
    """QC z-cap: an imputed |z| is clamped to the max observed |z| within
    +/- REGION_CAP_BP (per pleiodb)."""

    def _obs(self):
        # observed variants at 1.0/2.0/5.0 Mb with |z| = 3, 4, 2 (already sorted).
        pos = np.array([1_000_000, 2_000_000, 5_000_000], dtype=np.int64)
        absz = np.array([3.0, 4.0, 2.0], dtype=np.float64)
        return pos, absz

    def test_caps_overshoot_to_window_max(self):
        from opengwasdb.layouts.dense.complete import _region_capped_z

        pos, absz = self._obs()
        # at 1.5 Mb, window +/-1 Mb -> [1.0, 2.0] Mb (|z| 3, 4) -> cap to 4
        assert _region_capped_z(10.0, 1_500_000, pos, absz) == pytest.approx(4.0)
        assert _region_capped_z(-10.0, 1_500_000, pos, absz) == pytest.approx(-4.0)

    def test_within_window_unchanged(self):
        from opengwasdb.layouts.dense.complete import _region_capped_z

        pos, absz = self._obs()
        assert _region_capped_z(2.5, 1_500_000, pos, absz) == pytest.approx(2.5)

    def test_no_observed_in_window_uncapped(self):
        from opengwasdb.layouts.dense.complete import _region_capped_z

        pos, absz = self._obs()
        # 10 Mb is >1 Mb from every observed variant -> no window -> unchanged
        assert _region_capped_z(9.9, 10_000_000, pos, absz) == pytest.approx(9.9)

    def test_window_excludes_distant_large_z(self):
        from opengwasdb.layouts.dense.complete import _region_capped_z

        pos, absz = self._obs()
        # at 5 Mb, window includes only the 5 Mb observed (|z|=2) -> cap to 2,
        # NOT the larger |z|=4 at 2 Mb which is out of the +/-1 Mb window.
        assert _region_capped_z(8.0, 5_000_000, pos, absz) == pytest.approx(2.0)


class TestValidation:
    def test_valid_completed_store_passes(self, completed_store):
        result = validate_store(completed_store)
        assert result.ok, result.errors

    def test_observed_store_still_passes(self, observed_store):
        result = validate_store(observed_store)
        assert result.ok, result.errors

    def test_corrupt_off_panel_imputed_fails(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r+")
        on_panel = root["on_panel"][:]
        off_panel_rows = np.where(on_panel == 0)[0]
        assert len(off_panel_rows) > 0
        imputed = root["imputed"][:]
        imputed[off_panel_rows[0], 0] = 1
        root["imputed"][:] = imputed
        result = validate_store(completed_store)
        assert not result.ok

    def test_corrupt_imputed_nan_z_fails(self, completed_store):
        # An imputed=1 cell whose z is NaN must be caught by the streamed
        # finite-check band pass (issue 045), not just a full-matrix load. Mark
        # an on-panel cell imputed with a NaN z (on-panel so it isolates the
        # finite check from the off-panel-never-imputed check).
        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r+")
        on_panel = root["on_panel"][:]
        on_panel_rows = np.where(on_panel == 1)[0]
        assert len(on_panel_rows) > 0
        r, c = int(on_panel_rows[0]), 0
        imputed = root["imputed"][:]
        imputed[r, c] = 1
        root["imputed"][:] = imputed
        z = root["z"][:]
        z[r, c] = np.nan
        root["z"][:] = z
        result = validate_store(completed_store)
        assert not result.ok
        assert any("NaN z" in e for e in result.errors), result.errors

    def test_corrupt_top_hit_imputed_index_fails(self, completed_store):
        from opengwasdb.layouts.dense.top_hits import threshold_key

        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r+")
        group = root[f"top_hits/{threshold_key(5e-4)}"]
        assert "imputed" in group
        values = group["imputed"][:]
        assert len(values) > 0
        values[0] = 1 - values[0]
        group["imputed"][:] = values

        result = validate_store(completed_store)
        assert not result.ok
        assert any("imputed value inconsistent" in e for e in result.errors), result.errors


class TestSourceFidelity:
    """validate_store(source=...) cross-checks store values against the source.

    The observed store is same-assembly (GRCh38 in, GRCh38 store), so the join
    goes through the canonical-ALID path with no liftover — no bcftools needed.
    """

    def test_fidelity_passes_against_source(self, observed_store, source_path):
        result = validate_store(observed_store, source=source_path)
        assert result.ok, result.errors

    def test_fidelity_detects_sign_flip(self, observed_store, source_path):
        # Flip the sign of one observed, below-threshold cell: internal validation
        # still passes (|z| and the top-hit index are untouched), so only the
        # source-fidelity check can catch it.
        res = query_store(observed_store).lookup(["1:1000000:A:G"], ["a1"])
        assert len(res["z"]) == 1
        row, col = int(res["variant_index"][0]), int(res["analysis_index"][0])

        root = zarr.open_group(str(observed_store / "data.zarr"), mode="r+")
        z = root["z"][:]
        z[row, col] = -z[row, col]
        root["z"][:] = z

        assert validate_store(observed_store).ok  # internal checks still pass
        result = validate_store(observed_store, source=source_path)
        assert not result.ok
        assert any("source-fidelity" in e for e in result.errors), result.errors


class TestQuery:
    def test_analysis_a1_includes_imputed_by_default(self, completed_store):
        q = query_store(completed_store)
        result = q.analysis("a1")
        assert "association_status" in result
        statuses = set(result["association_status"].tolist())
        assert "observed" in statuses
        q.close()

    def test_analysis_observed_only_excludes_imputed(self, completed_store):
        q = query_store(completed_store)
        all_result = q.analysis("a1")
        obs_result = q.analysis("a1", observed_only=True)
        assert len(obs_result["z"]) <= len(all_result["z"])
        assert "imputed" not in set(obs_result["association_status"].tolist())
        q.close()

    def test_association_status_field_always_present(self, completed_store):
        q = query_store(completed_store)
        result = q.analysis("a1")
        assert len(result["association_status"]) == len(result["z"])
        q.close()

    def test_observed_store_has_status_observed_only(self, observed_store):
        q = query_store(observed_store)
        result = q.analysis("a1")
        assert set(result["association_status"].tolist()).issubset({"observed"})
        q.close()

    def test_top_hits_observed_only(self, completed_store):
        q = query_store(completed_store)
        result = q.top_hits(threshold=1.0, observed_only=True)
        assert "imputed" not in set(result["association_status"].tolist())
        q.close()

    def test_top_hit_index_stores_imputed_flags(self, completed_store):
        from opengwasdb.layouts.dense.top_hits import threshold_key

        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r")
        group = root[f"top_hits/{threshold_key(5e-4)}"]
        assert "imputed" in group

        rows = group["variant_index"][:].astype(np.int64)
        cols = group["analysis_index"][:].astype(np.int64)
        indexed = group["imputed"][:].astype(np.uint8)
        gathered = root["imputed"].vindex[rows, cols].astype(np.uint8)
        np.testing.assert_array_equal(indexed, gathered)

    def test_top_hits_uses_indexed_imputed_flags(self, completed_store):
        q = query_store(completed_store)

        def fail_imputed_pairs(rows, cols):  # noqa: ARG001
            raise AssertionError("top_hits should read indexed imputed flags")

        q._imputed_pairs = fail_imputed_pairs
        result = q.top_hits(threshold=5e-4)
        assert len(result["association_status"]) == len(result["z"])
        q.close()

    def test_top_hits_falls_back_without_indexed_imputed_flags(self, completed_store):
        from opengwasdb.layouts.dense.top_hits import threshold_key

        q = query_store(completed_store)
        indexed = q.top_hits(threshold=5e-4)
        q.close()

        root = zarr.open_group(str(completed_store / "data.zarr"), mode="r+")
        del root[f"top_hits/{threshold_key(5e-4)}"]["imputed"]

        q = query_store(completed_store)
        fallback = q.top_hits(threshold=5e-4)
        q.close()

        for name in ("variant_index", "analysis_index", "z", "se", "association_status"):
            np.testing.assert_array_equal(indexed[name], fallback[name])


class TestResume:
    def test_resume_completes_after_simulated_crash(
        self, tmp_path, observed_store, ld_panel, monkeypatch
    ):
        dst = tmp_path / "resumed.opengwasdb"
        calls = {"n": 0}
        real_run_block = complete_module._run_block

        def flaky_run_block(task):
            calls["n"] += 1
            if calls["n"] == 2:
                raise RuntimeError("simulated crash")
            return real_run_block(task)

        monkeypatch.setattr(complete_module, "_run_block", flaky_run_block)
        with pytest.raises(RuntimeError, match="simulated crash"):
            complete_dense_store(observed_store, dst, ld_panel, ancestry="EUR", min_cor=0.0)

        checkpoint_dir = complete_module._checkpoint_dir_for(dst)
        assert checkpoint_dir.exists()
        assert not dst.exists()
        assert not complete_module._work_dir_for(dst).exists()
        checkpointed = list((checkpoint_dir / "blocks").glob("*.npz"))
        assert len(checkpointed) == 1  # exactly one block finished before the crash

        monkeypatch.undo()
        result = resume_dense_completion(checkpoint_dir)
        assert dst.exists()
        assert not checkpoint_dir.exists()
        assert result.n_variants > 0

        assert validate_store(dst).ok

    def test_resume_matches_fresh_run(self, tmp_path, observed_store, ld_panel):
        fresh_dst = tmp_path / "fresh.opengwasdb"
        fresh = complete_dense_store(
            observed_store, fresh_dst, ld_panel, ancestry="EUR", min_cor=0.0
        )

        resumable_dst = tmp_path / "resumable.opengwasdb"
        checkpoint_dir = complete_module._checkpoint_dir_for(resumable_dst)
        checkpoint_dir.mkdir(parents=True)
        (checkpoint_dir / "blocks").mkdir()
        import json
        (checkpoint_dir / "build_params.json").write_text(
            json.dumps({
                "source_path": str(Path(observed_store).resolve()),
                "dest_path": str(resumable_dst.resolve()),
                "ld_dir": str(Path(ld_panel).resolve()),
                "ancestry": "EUR", "min_cor": 0.0, "thresh": 0.9,
                "release_id": None, "ld_panel_id": "eur-hg38-gpm",
            }),
            encoding="utf-8",
        )
        resumed = resume_dense_completion(checkpoint_dir)

        assert resumed.n_variants == fresh.n_variants
        assert resumed.n_imputed == fresh.n_imputed
        assert resumed.n_missing_off_panel == fresh.n_missing_off_panel
        assert resumed.n_missing_imputation_failed == fresh.n_missing_imputation_failed
        _assert_stores_identical(resumable_dst, fresh_dst)


class TestParallel:
    def test_two_workers_matches_serial(self, tmp_path, observed_store, ld_panel):
        serial_dst = tmp_path / "serial.opengwasdb"
        serial = complete_dense_store(
            observed_store, serial_dst, ld_panel, ancestry="EUR", min_cor=0.0, n_workers=1
        )
        parallel_dst = tmp_path / "parallel.opengwasdb"
        parallel = complete_dense_store(
            observed_store, parallel_dst, ld_panel, ancestry="EUR", min_cor=0.0, n_workers=2
        )
        assert parallel.n_variants == serial.n_variants
        assert parallel.n_imputed == serial.n_imputed
        assert parallel.n_missing_off_panel == serial.n_missing_off_panel
        assert parallel.n_missing_imputation_failed == serial.n_missing_imputation_failed
        _assert_stores_identical(parallel_dst, serial_dst)
        assert validate_store(parallel_dst).ok
