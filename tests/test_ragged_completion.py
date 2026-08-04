"""Integration tests for ragged reference completion (issues 039-042)."""

from __future__ import annotations

import gzip
import io
import struct
from pathlib import Path

import numpy as np
import pytest
import zarr

from opengwasdb.layouts.ragged.build_besd import build_ragged_from_besd
from opengwasdb.layouts.ragged.complete import complete_ragged_store
from opengwasdb.query import query_store
from opengwasdb.validation.validate import validate_store


# ── Synthetic BESD fixture (reused from test_ragged_build_besd.py) ─────────

def _write_esi(path: Path, snps: list[dict]) -> None:
    with open(path, "w") as fh:
        for s in snps:
            fh.write(f"{s['chr']}\t{s['snp_id']}\t0\t{s['bp']}\t{s['a1']}\t{s['a2']}\tNA\n")


def _write_epi(path: Path, probes: list[dict]) -> None:
    with open(path, "w") as fh:
        for p in probes:
            fh.write(f"{p['chr']}\t{p['probe_id']}\t0\t{p['bp']}\t{p.get('gene', 'NA')}\t+\n")


def _write_besd_sparse_3f(path: Path, n_probes: int, probe_assocs: list[list[tuple]]) -> None:
    rowid, val, cols = [], [], []
    offset = 0
    for assocs in probe_assocs:
        n = len(assocs)
        cols.append(offset)
        cols.append(offset + n)
        for snp_idx, beta, _ in assocs:
            rowid.append(snp_idx); val.append(beta)
        for _, _, se in assocs:
            rowid.append(0); val.append(se)
        offset += 2 * n
    cols.append(offset)
    col_num = (n_probes << 1) + 1
    with open(path, "wb") as fh:
        fh.write(struct.pack("<I", 0x40400000))
        fh.write(struct.pack("<Q", len(val)))
        fh.write(struct.pack(f"<{col_num}Q", *cols))
        fh.write(struct.pack(f"<{len(val)}I", *rowid))
        fh.write(struct.pack(f"<{len(val)}f", *val))


def _make_besd_fixture(tmp_path: Path) -> Path:
    fixture = tmp_path / "fixture"
    fixture.mkdir()
    snps = [
        {"chr": "1", "snp_id": "rs1001", "bp": 1_000_000, "a1": "A", "a2": "G"},
        {"chr": "1", "snp_id": "rs1002", "bp": 1_100_000, "a1": "C", "a2": "T"},
        {"chr": "1", "snp_id": "rs1003", "bp": 1_200_000, "a1": "A", "a2": "C"},
    ]
    probes = [
        {"chr": "1", "probe_id": "ENSG00000000001", "bp": 1_050_000, "gene": "GENE1"},
        {"chr": "1", "probe_id": "ENSG00000000002", "bp": 1_150_000, "gene": "GENE2"},
    ]
    probe_assocs = [
        [(0, 0.1, 0.02), (1, -0.2, 0.03)],
        [(1, 0.5, 0.05), (2, -0.15, 0.025)],
    ]
    _write_esi(fixture / "test.esi", snps)
    _write_epi(fixture / "test.epi", probes)
    _write_besd_sparse_3f(fixture / "test.besd", len(probes), probe_assocs)
    return fixture / "test"


# ── Synthetic LD panel fixture ──────────────────────────────────────────────

def _make_ld_panel(tmp_path: Path, chrom: str, start: int, end: int) -> Path:
    """Create a tiny synthetic LD panel block directory (flat layout)."""
    block_name = f"{start}-{end}"
    panel_dir = tmp_path / "ld_panel" / "EUR" / chrom
    panel_dir.mkdir(parents=True)

    # Two reference-panel variants, one of which is NOT in the observed store
    ref_snps = [
        ("1:1000000_A_G", 0.3, 1_000_000),   # same as rs1001
        ("1:1050500_C_T", 0.4, 1_050_500),   # new — not in observed store
        ("1:1100000_C_T", 0.45, 1_100_000),  # same as rs1002
    ]
    n = len(ref_snps)

    tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
    for alid, eaf, bp in ref_snps:
        chrom_, rest = alid.split(":")
        bp_str = rest.split("_")[0]
        a1, a2 = rest.split("_")[1], rest.split("_")[2]
        tsv_lines.append(f"{chrom_}\t{alid}\t{a2}\t{a1}\t{eaf}\t{bp_str}")
    (panel_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

    # Random positive-definite LD matrix
    rng = np.random.default_rng(0)
    A = rng.standard_normal((n, n))
    ld = A @ A.T + np.eye(n) * n * 0.1

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for row in ld:
            gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
    (panel_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())

    return tmp_path / "ld_panel"


# ── Tests ───────────────────────────────────────────────────────────────────

@pytest.fixture
def observed_store(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "obs.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="obs-v1", tissue="Blood")
    return out


@pytest.fixture
def ld_panel(tmp_path):
    return _make_ld_panel(tmp_path, "1", 900_000, 1_300_000)


@pytest.fixture
def completed_store(tmp_path, observed_store, ld_panel):
    dst = tmp_path / "comp.opengwasdb"
    complete_ragged_store(
        observed_store, dst, ld_panel,
        ancestry="EUR", cis_window_bp=500_000, min_cor=0.0, release_id="comp-v1",
    )
    return dst


class TestCompletionFiles:
    def test_creates_store_directory(self, completed_store):
        assert completed_store.exists()
        assert (completed_store / "manifest.json").exists()
        assert (completed_store / "variants.tsv.gz").exists()
        assert (completed_store / "traits.tsv.gz").exists()
        assert (completed_store / "index.sqlite").exists()

    def test_manifest_completion_state(self, completed_store):
        import json
        m = json.loads((completed_store / "manifest.json").read_text())
        assert m["completion_state"] == "reference_completed"
        assert "completion" in m["provenance"]

    def test_imputed_array_present(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr" / "ragged"), mode="r")
        assert "imputed" in root

    def test_imputed_array_aligned_with_z(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr" / "ragged"), mode="r")
        assert len(root["imputed"]) == len(root["z"])

    def test_imputed_values_are_0_or_1(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr" / "ragged"), mode="r")
        imp = root["imputed"][:]
        assert np.all((imp == 0) | (imp == 1))

    def test_imputed_1_rows_have_finite_z(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr" / "ragged"), mode="r")
        z = root["z"][:].astype("float32")
        imp = root["imputed"][:]
        assert np.all(np.isfinite(z[imp == 1]))

    def test_variant_table_larger_than_source(self, completed_store, observed_store):
        from opengwasdb.variants.axis import VariantAxis
        obs_ax = VariantAxis(observed_store)
        comp_ax = VariantAxis(completed_store)
        assert comp_ax.n_variants >= obs_ax.n_variants
        obs_ax.close()
        comp_ax.close()

    def test_completion_quality_table_exists(self, completed_store):
        import sqlite3
        conn = sqlite3.connect(str(completed_store / "index.sqlite"))
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        conn.close()
        assert "completion_quality" in tables

    def test_overwrite_raises_without_flag(self, tmp_path, observed_store, ld_panel, completed_store):
        with pytest.raises(FileExistsError):
            complete_ragged_store(observed_store, completed_store, ld_panel, min_cor=0.0)


class TestValidation:
    def test_valid_completed_store_passes(self, completed_store):
        result = validate_store(completed_store)
        assert result.ok, result.errors

    def test_observed_store_still_passes(self, observed_store):
        result = validate_store(observed_store)
        assert result.ok, result.errors

    def test_corrupt_imputed_array_fails(self, completed_store):
        root = zarr.open_group(str(completed_store / "data.zarr" / "ragged"), mode="r+")
        imp = root["imputed"][:]
        # Set imputed=1 where z is NaN
        z = root["z"][:].astype("float32")
        nan_positions = np.where(~np.isfinite(z))[0]
        if len(nan_positions) > 0:
            imp[nan_positions[0]] = 1
            root["imputed"][:] = imp
            result = validate_store(completed_store)
            assert not result.ok


class TestQuery:
    def test_selected_top_hits_match_global_and_filter_imputed(self, completed_store):
        q = query_store(completed_store)
        global_result = q.top_hits(threshold=5e-4)
        selected = q.top_hits(analysis_id="ENSG00000000001", threshold=5e-4)
        observed = q.top_hits(
            analysis_id="ENSG00000000001", threshold=5e-4, observed_only=True
        )
        expected = global_result["analysis_index"] == 0
        for name in ("variant_index", "analysis_index", "z", "se", "association_status"):
            np.testing.assert_array_equal(selected[name], global_result[name][expected])
        assert "imputed" not in set(observed["association_status"].tolist())
        q.close()

    def test_analysis_returns_imputed_by_default(self, completed_store):
        q = query_store(completed_store)
        result = q.analysis("ENSG00000000001")
        assert "association_status" in result
        statuses = set(result["association_status"].tolist())
        # Should have at least "observed"; may also have "imputed" or "missing"
        assert "observed" in statuses
        q.close()

    def test_analysis_observed_only_excludes_imputed(self, completed_store):
        q = query_store(completed_store)
        all_result = q.analysis("ENSG00000000001")
        obs_result = q.analysis("ENSG00000000001", observed_only=True)
        # observed-only must be a subset
        assert len(obs_result["z"]) <= len(all_result["z"])
        statuses = set(obs_result["association_status"].tolist())
        assert "imputed" not in statuses
        q.close()

    def test_association_status_field_always_present(self, completed_store):
        q = query_store(completed_store)
        result = q.analysis("ENSG00000000001")
        assert "association_status" in result
        assert len(result["association_status"]) == len(result["z"])
        q.close()

    def test_observed_store_has_status_observed(self, observed_store):
        q = query_store(observed_store)
        result = q.analysis("ENSG00000000001")
        assert "association_status" in result
        statuses = set(result["association_status"].tolist())
        assert statuses.issubset({"observed"})
        q.close()

    def test_range_by_analysis_observed_only(self, completed_store):
        q = query_store(completed_store)
        result = q.range_by_analysis("1", 900_000, 1_300_000, observed_only=True)
        statuses = set(result["association_status"].tolist())
        assert "imputed" not in statuses
        q.close()
