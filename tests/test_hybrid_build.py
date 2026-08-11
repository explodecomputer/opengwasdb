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


def _make_manifest(
    tmp_path: Path,
    entries: list[tuple[str, Path, str]],
    scales: dict[str, str] | None = None,
    sd_methods: dict[str, str] | None = None,
    sds: dict[str, str] | None = None,
) -> Path:
    """Write the build manifest; see the identical helper in
    test_dense_vcf_build.py for the full rationale -- this builder shares
    its manifest reader with the dense one."""
    scales = scales or {}
    sd_methods = sd_methods or {}
    sds = sds or {}
    manifest = tmp_path / "manifest.tsv"
    lines = [
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale"
        "\toriginal_sd_method\toriginal_sd"
    ]
    for trait_id, file_path, trait_name in entries:
        scale = scales.get(trait_id, "sd")
        sd_method = sd_methods.get(trait_id, "declared_standardised")
        sd = sds.get(trait_id, "")
        lines.append(
            f"{trait_id}\t{file_path}\t{trait_name}\t1000\t{scale}\t{sd_method}\t{sd}"
        )
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


def test_stored_effect_scale_comes_from_manifest_not_header(tmp_path):
    """The ieu-a-7 fix (issue #17), hybrid path: a VCF header says
    ``StudyType=Continuous`` but the manifest declares ``log_or`` -- the
    built store must record the manifest's value."""
    vcf = _make_vcf(
        tmp_path,
        "trait_a",
        [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"],
        study_type="Continuous",
    )
    manifest = _make_manifest(
        tmp_path, [("trait_a", vcf, "Trait A")], scales={"trait_a": "log_or"}
    )
    store_path = tmp_path / "store.opengwasdb"
    build_hybrid_from_vcf_manifest(
        manifest, store_path, reference_panel=_panel(tmp_path),
        store_id="s", release_id="r",
    )

    q = query_store(store_path)
    analyses = q.analyses_table()
    assert analyses[0]["stored_effect_scale"] == "log_or"


def test_missing_required_manifest_field_fails_the_build_loudly(tmp_path):
    """A manifest missing stored_effect_scale must fail the build with a
    clear error before any I/O, not fall back to VCF-header inference or a
    silent default (issue #17)."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text(
        "trait_id\tfile_path\ttrait_name\tn\n"
        f"trait_a\t{vcf}\tTrait A\t1000\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="stored_effect_scale"):
        build_hybrid_from_vcf_manifest(
            manifest_path, tmp_path / "store.opengwasdb",
            reference_panel=_panel(tmp_path), store_id="s", release_id="r",
        )


def test_continuous_trait_rescaled_by_manifest_original_sd(tmp_path):
    """issue #18 AC1, hybrid path: a continuous-trait Analysis with a
    manifest-supplied original_sd != 1 has its se divided by that SD, for
    both the on-panel (Dense Component) and off-panel (overflow) associations
    of the same study."""
    vcf = _make_vcf(
        tmp_path,
        "trait_a",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # on-panel, se=0.5
            f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.4\n",  # off-panel, se=0.4
        ],
    )
    manifest = _make_manifest(
        tmp_path,
        [("trait_a", vcf, "Trait A")],
        sd_methods={"trait_a": "source_provided"},
        sds={"trait_a": "2.0"},
    )
    store_path = tmp_path / "store.opengwasdb"
    build_hybrid_from_vcf_manifest(
        manifest, store_path, reference_panel=_panel(tmp_path), store_id="s", release_id="r",
    )

    q = query_store(store_path)
    on_panel = q.lookup([HG38_ALID_1], ["trait_a"])
    off_panel = q.lookup([HG38_ALID_2], ["trait_a"])
    assert on_panel["se"][0] == pytest.approx(0.25, rel=5e-3)  # 0.5 / 2.0
    assert off_panel["se"][0] == pytest.approx(0.2, rel=5e-3)  # 0.4 / 2.0


def test_binary_trait_never_rescaled(tmp_path):
    """issue #18 AC2, hybrid path: original_sd_method=binary_trait is never
    rescaled, for either the on-panel or off-panel associations of the study."""
    vcf = _make_vcf(
        tmp_path,
        "trait_b",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # on-panel, se=0.5
            f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.4\n",  # off-panel, se=0.4
        ],
    )
    manifest = _make_manifest(
        tmp_path,
        [("trait_b", vcf, "Trait B")],
        scales={"trait_b": "log_or"},
        sd_methods={"trait_b": "binary_trait"},
    )
    store_path = tmp_path / "store.opengwasdb"
    build_hybrid_from_vcf_manifest(
        manifest, store_path, reference_panel=_panel(tmp_path), store_id="s", release_id="r",
    )

    q = query_store(store_path)
    on_panel = q.lookup([HG38_ALID_1], ["trait_b"])
    off_panel = q.lookup([HG38_ALID_2], ["trait_b"])
    assert on_panel["se"][0] == pytest.approx(0.5, rel=5e-3)
    assert off_panel["se"][0] == pytest.approx(0.4, rel=5e-3)


def test_original_sd_method_unavailable_fails_the_build_loudly(tmp_path):
    """issue #18 AC3, hybrid path: original_sd_method=unavailable fails the
    build rather than assuming sd=1."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest = _make_manifest(
        tmp_path, [("trait_a", vcf, "Trait A")], sd_methods={"trait_a": "unavailable"}
    )

    with pytest.raises(ValueError, match="unavailable"):
        build_hybrid_from_vcf_manifest(
            manifest, tmp_path / "store.opengwasdb",
            reference_panel=_panel(tmp_path), store_id="s", release_id="r",
        )


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
    assert list(zip(r["analysis_index"], r["variant_index"], strict=True)) == sorted(
        zip(r["analysis_index"], r["variant_index"], strict=True)
    )


def test_top_hits_selects_and_merges_one_analysis(hybrid_store):
    q = query_store(hybrid_store)
    global_result = q.top_hits(threshold=5e-4)
    selected = q.top_hits(analysis_id="trait_a", threshold=5e-4)
    expected = global_result["analysis_index"] == 0

    for name in ("variant_index", "analysis_index", "z", "se", "association_status"):
        np.testing.assert_array_equal(selected[name], global_result[name][expected])
    assert len(selected["z"]) == 2  # one Dense hit and one overflow hit
    assert q.top_hits(analysis_id="unknown", threshold=5e-4)["z"].size == 0


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


def test_validate_catches_overflow_top_hit_offsets(hybrid_store):
    root = zarr.open_group(str(hybrid_store / "data.zarr"), mode="r+")
    offsets = root["top_hits/p_5e_04/analysis_offsets"][:]
    offsets[-1] -= 1
    root["top_hits/p_5e_04/analysis_offsets"][:] = offsets

    result = validate_store(hybrid_store)
    assert not result.ok
    assert any("invalid analysis offsets" in error for error in result.errors)
