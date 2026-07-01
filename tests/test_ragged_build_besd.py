"""Integration test for ragged BESD builder (issue 036)."""

import struct
from pathlib import Path

import numpy as np
import pytest
import zarr

from opengwasdb.layouts.ragged.build_besd import build_ragged_from_besd
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.traits.axis import TraitsAxisReader
from opengwasdb.variants.axis import VariantAxis


# ── Synthetic BESD fixture ────────────────────────────────────────────────────

def _write_esi(path: Path, snps: list[dict]) -> None:
    with open(path, "w") as fh:
        for s in snps:
            freq = s.get("freq", "NA")
            fh.write(f"{s['chr']}\t{s['snp_id']}\t0\t{s['bp']}\t{s['a1']}\t{s['a2']}\t{freq}\n")


def _write_epi(path: Path, probes: list[dict]) -> None:
    with open(path, "w") as fh:
        for p in probes:
            gene = p.get("gene", "NA")
            fh.write(f"{p['chr']}\t{p['probe_id']}\t0\t{p['bp']}\t{gene}\t+\n")


def _write_besd_sparse_3f(
    path: Path,
    n_probes: int,
    probe_assocs: list[list[tuple[int, float, float]]],
) -> None:
    """Write a minimal SPARSE_FILE_TYPE_3F BESD file.

    BESD layout: for each probe p, val contains [betas | SEs] and rowid contains
    [snp_indices | zeros] at the same offsets. cols[2p]=beta_start, cols[2p+1]=se_start.
    """
    rowid: list[int] = []
    val: list[float] = []
    cols: list[int] = []
    offset = 0

    for assocs in probe_assocs:
        n = len(assocs)
        cols.append(offset)          # beta_start for this probe
        cols.append(offset + n)      # se_start for this probe
        # Beta positions: meaningful snp_idx and beta values
        for snp_idx, beta, _ in assocs:
            rowid.append(snp_idx)
            val.append(beta)
        # SE positions: snp_idx is don't-care; SE values follow immediately
        for _, _, se in assocs:
            rowid.append(0)
            val.append(se)
        offset += 2 * n

    cols.append(offset)   # final sentinel
    val_num = len(val)
    col_num = (n_probes << 1) + 1

    with open(path, "wb") as fh:
        fh.write(struct.pack("<I", 0x40400000))         # magic 3F
        fh.write(struct.pack("<Q", val_num))             # val_num
        fh.write(struct.pack(f"<{col_num}Q", *cols))    # column offsets
        fh.write(struct.pack(f"<{val_num}I", *rowid))   # row indices
        fh.write(struct.pack(f"<{val_num}f", *val))     # beta/SE values


def _make_besd_fixture(tmp_path: Path) -> Path:
    """Create a 3-probe synthetic BESD dataset in tmp_path/fixture/."""
    fixture = tmp_path / "fixture"
    fixture.mkdir()

    snps = [
        {"chr": "1", "snp_id": "rs1001", "bp": 1_000_000, "a1": "A", "a2": "G"},
        {"chr": "1", "snp_id": "rs1002", "bp": 1_100_000, "a1": "C", "a2": "T"},
        {"chr": "1", "snp_id": "rs1003", "bp": 1_200_000, "a1": "A", "a2": "C"},
        {"chr": "2", "snp_id": "rs2001", "bp": 2_000_000, "a1": "G", "a2": "T"},
        {"chr": "2", "snp_id": "rs2002", "bp": 2_100_000, "a1": "A", "a2": "T"},
    ]
    probes = [
        {"chr": "1", "probe_id": "ENSG00000000001", "bp": 1_050_000, "gene": "GENE1"},
        {"chr": "1", "probe_id": "ENSG00000000002", "bp": 1_150_000, "gene": "GENE2"},
        {"chr": "2", "probe_id": "ENSG00000000003", "bp": 2_050_000, "gene": "GENE3"},
    ]
    # probe 0: SNPs 0,1 (rs1001, rs1002)
    # probe 1: SNP 2 (rs1003)
    # probe 2: SNPs 3,4 (rs2001, rs2002)
    probe_assocs = [
        [(0, 0.1, 0.02), (1, -0.2, 0.03)],
        [(2, 0.5, 0.05)],
        [(3, 0.3, 0.04), (4, -0.15, 0.025)],
    ]

    _write_esi(fixture / "test.esi", snps)
    _write_epi(fixture / "test.epi", probes)
    _write_besd_sparse_3f(fixture / "test.besd", len(probes), probe_assocs)
    return fixture / "test"


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_build_creates_store_files(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"

    result = build_ragged_from_besd(
        prefix, out,
        store_id="test", release_id="v1",
        tissue="Whole_Blood",
    )

    assert out.exists()
    assert (out / "manifest.json").exists()
    assert (out / "variants.tsv.gz").exists()
    assert (out / "traits.tsv.gz").exists()
    assert (out / "index.sqlite").exists()
    assert (out / "data.zarr" / "ragged").exists()

    assert result.n_variants == 5
    assert result.n_analyses == 3


def test_manifest_primary_layout(tmp_path):
    import json
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["primary_layout"] == "ragged"
    assert manifest["completion_state"] == "observed_only"
    assert manifest["reference_assembly"] == "GRCh38"


def test_traits_tsv_contents(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1", tissue="Whole_Blood")

    reader = TraitsAxisReader(out)
    all_traits = list(reader.all())
    assert len(all_traits) == 3
    probe_ids = {r.trait_id for r in all_traits}
    assert probe_ids == {"ENSG00000000001", "ENSG00000000002", "ENSG00000000003"}

    # regional query — chr1 probes
    chr1_probes = reader.range("1", 1_000_000, 1_200_000)
    assert len(chr1_probes) == 2
    assert {r.trait_id for r in chr1_probes} == {"ENSG00000000001", "ENSG00000000002"}


def test_zarr_csr_associations(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    csr = RaggedCSRReader(out)
    assert csr.n_analyses == 3
    assert csr.n_associations == 5  # 2 + 1 + 2

    # Probe 0 has 2 associations
    a0 = csr.get_analysis(0)
    assert len(a0.variant_index) == 2

    # Probe 1 has 1 association
    a1 = csr.get_analysis(1)
    assert len(a1.variant_index) == 1

    # Probe 2 has 2 associations
    a2 = csr.get_analysis(2)
    assert len(a2.variant_index) == 2


def test_z_scores_computed_correctly(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    csr = RaggedCSRReader(out)
    # Probe 1: single assoc beta=0.5, se=0.05 → z=10.0 (stored as float16)
    a1 = csr.get_analysis(1)
    assert len(a1.z) == 1
    assert abs(float(a1.z[0]) - 10.0) < 0.1  # float16 tolerance


def test_top_hit_index_built_inline(tmp_path):
    """build_ragged_from_besd auto-builds the top-hit index."""
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    root = zarr.open_group(str(out / "data.zarr"), mode="r")
    assert "top_hits" in root
    # At least one threshold group must exist
    assert len(list(root["top_hits"].keys())) > 0


def test_top_hit_index_schema(tmp_path):
    """Each threshold group has the expected arrays and attribute."""
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")
    build_ragged_top_hit_indexes(out)  # idempotent rebuild

    root = zarr.open_group(str(out / "data.zarr"), mode="r")
    for key in root["top_hits"]:
        group = root["top_hits"][key]
        for name in ("variant_index", "analysis_index", "abs_z", "z", "se", "p_value"):
            assert name in group, f"missing {name} in {key}"
        n = len(group["variant_index"])
        for name in ("analysis_index", "abs_z", "z", "se", "p_value"):
            assert len(group[name]) == n
        assert "threshold" in group.attrs


def test_top_hit_z_values_match_csr(tmp_path):
    """Top-hit z values round-trip through the float16 CSR correctly."""
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    csr = RaggedCSRReader(out)
    root = zarr.open_group(str(out / "data.zarr"), mode="r")

    # Use the loosest threshold to get all 5 hits
    loosest_key = sorted(root["top_hits"].keys())[-1]
    group = root["top_hits"][loosest_key]
    vis = group["variant_index"][:].astype(int)
    ais = group["analysis_index"][:].astype(int)
    zs = group["z"][:].astype("float32")

    offsets = csr._offsets[:]
    vi_all = csr._variant_index[:]
    z_all = csr._z[:]

    for vi, ai, z_hit in zip(vis, ais, zs):
        start, end = int(offsets[ai]), int(offsets[ai + 1])
        pos = start + int(np.searchsorted(vi_all[start:end], vi))
        assert pos < end and int(vi_all[pos]) == vi
        assert np.isclose(float(z_all[pos]), float(z_hit), rtol=1e-2, atol=1e-2)


def test_top_hits_query_uses_index(tmp_path):
    """RaggedStoreQuery.top_hits reads from the precomputed index."""
    from opengwasdb.query import query_store
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    q = query_store(out)
    # Threshold 5e-4 should return all 5 associations (all have |z| > 3.48)
    result = q.top_hits(threshold=5e-4)
    assert len(result["z"]) == 5
    assert "variant_index" in result and "analysis_index" in result


def test_overwrite_flag(tmp_path):
    prefix = _make_besd_fixture(tmp_path)
    out = tmp_path / "out.opengwasdb"

    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")
    with pytest.raises(FileExistsError):
        build_ragged_from_besd(prefix, out, store_id="test", release_id="v1")

    # With overwrite=True should succeed
    build_ragged_from_besd(prefix, out, store_id="test", release_id="v1", overwrite=True)
    assert out.exists()
