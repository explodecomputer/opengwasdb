"""End-to-end tests for the Hybrid build + unified query (issues 055/056).

Reuses the dense VCF fixture coordinates. The reference panel deliberately holds
only two of the three variants so the third routes to the Ragged Overflow.

Panel (on-panel → Dense Component):
  1:100000:A:G     (HG38_ALID_1)
  1:1564620:A:G    (HG38_ALID_3)
Off-panel (→ Ragged Overflow):
  1:1064620:C:T    (HG38_ALID_2)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import zarr

from opengwasdb.layouts.hybrid.build import build_hybrid_from_vcf_manifest
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.query import query_store
from opengwasdb.validation import validate_store

HG19_POS_1 = 100_000
HG19_POS_2 = 1_000_000
HG19_POS_3 = 1_500_000

HG38_ALID_1 = "1:100000:A:G"
HG38_ALID_2 = "1:1064620:C:T"
HG38_ALID_3 = "1:1564620:A:G"


def _vcf_header(study_type: str = "Continuous") -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size\">\n"
        "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error\">\n"
        "##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"Z-score\">\n"
        f"##SAMPLE=<ID=STUDY1,StudyType={study_type}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1\n"
    )


def _make_vcf(tmp_path: Path, name: str, rows: list[str], study_type: str = "Continuous") -> Path:
    path = tmp_path / f"{name}.vcf"
    path.write_text(_vcf_header(study_type) + "".join(rows), encoding="utf-8")
    return path


def _make_manifest(tmp_path: Path, entries: list[tuple[str, Path, str]]) -> Path:
    manifest = tmp_path / "manifest.tsv"
    lines = ["trait_id\tfile_path\ttrait_name\tn"]
    for trait_id, file_path, trait_name in entries:
        lines.append(f"{trait_id}\t{file_path}\t{trait_name}\t1000")
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def _panel(tmp_path: Path) -> Path:
    panel = tmp_path / "panel.txt"
    panel.write_text(f"{HG38_ALID_1}\n{HG38_ALID_3}\n", encoding="utf-8")
    return panel


@pytest.fixture(params=[1, 2], ids=["serial", "parallel"])
def hybrid_store(tmp_path, request):
    vcf1 = _make_vcf(
        tmp_path,
        "trait_a",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # z -4.0 dense
            f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",  # z -5.0 OVERFLOW
            f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",  # z  3.0 dense
        ],
    )
    vcf2 = _make_vcf(
        tmp_path,
        "trait_b",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t6.0:0.5\n",  # z -12.0 dense
            f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t1.2:0.3\n",  # z   4.0 dense
        ],
    )
    manifest = _make_manifest(
        tmp_path, [("trait_a", vcf1, "Trait A"), ("trait_b", vcf2, "Trait B")]
    )
    store_path = tmp_path / f"store_{request.param}.opengwasdb"
    build_hybrid_from_vcf_manifest(
        manifest,
        store_path,
        reference_panel=_panel(tmp_path),
        store_id="hybrid-test",
        release_id="v1",
        n_workers=request.param,
    )
    return store_path


def test_manifest_is_hybrid(hybrid_store):
    manifest = StoreManifest.load(hybrid_store)
    assert manifest.primary_layout.value == "hybrid"
    assert manifest.association_coverage.value == "full"


def test_store_envelope(hybrid_store):
    assert (hybrid_store / "manifest.json").exists()
    assert (hybrid_store / "variants.tsv.gz").exists()
    assert (hybrid_store / "dense" / "data.zarr").exists()
    assert (hybrid_store / "data.zarr" / "ragged").exists()
    assert (hybrid_store / "dense" / "dense_to_shared.npy").exists()


def test_dense_matrix_is_panel_sized(hybrid_store):
    root = zarr.open_group(str(hybrid_store / "dense" / "data.zarr"), mode="r")
    # 2 panel variants × 2 analyses (off-panel variant is NOT a dense row).
    assert root["z"].shape == (2, 2)


def test_shared_table_is_union(hybrid_store):
    manifest = StoreManifest.load(hybrid_store)
    assert manifest.provenance["hybrid"]["n_panel"] == 2
    assert manifest.provenance["hybrid"]["n_off_panel"] == 1
    assert manifest.provenance["n_variants"] == 3


def test_lookup_on_panel_hits_dense(hybrid_store):
    q = query_store(hybrid_store)
    r = q.lookup([HG38_ALID_1], ["trait_a"])
    assert len(r["z"]) == 1
    assert r["z"][0] == pytest.approx(-4.0, rel=5e-3)
    assert r["association_status"][0] == "observed"


def test_lookup_off_panel_hits_overflow(hybrid_store):
    q = query_store(hybrid_store)
    r = q.lookup([HG38_ALID_2], ["trait_a"])
    assert len(r["z"]) == 1
    assert r["z"][0] == pytest.approx(-5.0, rel=5e-3)


def test_off_panel_variant_absent_from_dense_matrix(hybrid_store):
    """The overflow variant must NOT occupy a dense row (disjoint partition)."""
    q = query_store(hybrid_store)
    variants = q.variants_table()
    off_shared = next(
        k for k, v in variants.items() if v["alid"] == HG38_ALID_2
    )
    dense_to_shared = np.load(hybrid_store / "dense" / "dense_to_shared.npy")
    assert off_shared not in dense_to_shared.tolist()


def test_phewas_off_panel(hybrid_store):
    q = query_store(hybrid_store)
    r = q.phewas(HG38_ALID_2)
    # Only trait_a observed the off-panel variant.
    assert len(r["z"]) == 1
    assert r["z"][0] == pytest.approx(-5.0, rel=5e-3)


def test_analysis_unions_both_components(hybrid_store):
    q = query_store(hybrid_store)
    r = q.analysis("trait_a")
    # trait_a: 2 dense (ALID_1, ALID_3) + 1 overflow (ALID_2) = 3
    assert len(r["z"]) == 3
    assert set(np.round(r["z"], 1).tolist()) == {-4.0, -5.0, 3.0}


def test_range_phewas_unions_components(hybrid_store):
    q = query_store(hybrid_store)
    # Range covering all three variants (1:100000 .. 1:1564620).
    r = q.range_phewas("1", 1, 2_000_000)
    variants = q.variants_table()
    alids = {variants[int(vi)]["alid"] for vi in r["variant_index"]}
    assert {HG38_ALID_1, HG38_ALID_2, HG38_ALID_3} <= alids


def test_top_hits_merges_components(hybrid_store):
    q = query_store(hybrid_store)
    r = q.top_hits(threshold=5e-4)
    variants = q.variants_table()
    alids = {variants[int(vi)]["alid"] for vi in r["variant_index"]}
    # The off-panel overflow hit (z=-5.0) must appear alongside dense hits.
    assert HG38_ALID_2 in alids
    # Ranked by descending |z|: the strongest is trait_b's -12.0 (dense).
    assert abs(r["z"][0]) == pytest.approx(12.0, rel=5e-3)


# ── Validation (issue 059) ───────────────────────────────────────────────────


def test_validate_passes(hybrid_store):
    result = validate_store(hybrid_store)
    assert result.ok, result.errors


def test_validate_catches_imputed_overflow(hybrid_store):
    """An imputed array on the overflow must fail — the overflow is never imputed."""
    import zarr

    root = zarr.open_group(str(hybrid_store / "data.zarr" / "ragged"), mode="a")
    n = int(root["offsets"][:][-1])
    root.create_dataset("imputed", data=np.zeros(max(n, 1), dtype="uint8"))
    result = validate_store(hybrid_store)
    assert not result.ok
    assert any("overflow" in e.lower() and "imputed" in e.lower() for e in result.errors)


def test_validate_catches_disjoint_violation(hybrid_store):
    """Point an overflow entry at an on-panel variant → disjoint partition fails."""
    import zarr

    dense_to_shared = np.load(hybrid_store / "dense" / "dense_to_shared.npy")
    on_panel_shared = int(dense_to_shared[0])
    root = zarr.open_group(str(hybrid_store / "data.zarr" / "ragged"), mode="a")
    vi = root["variant_index"][:]
    if len(vi) == 0:
        pytest.skip("no overflow associations to corrupt")
    vi[0] = on_panel_shared
    root["variant_index"][:] = vi
    result = validate_store(hybrid_store)
    assert not result.ok
    assert any("disjoint" in e.lower() for e in result.errors)
