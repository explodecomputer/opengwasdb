"""End-to-end tests for the two-pass VCF dense build pipeline.

All fixtures use synthetic GWAS-VCF files written to tmp_path.
Real pyliftover is used with known hg19 positions that map successfully to hg38.

Known positions:
  hg19 1:100000  → hg38 1:100000   (REF=A, ALT=G → ALID 1:100000:A:G, flip=True  → z=-z)
  hg19 1:1000000 → hg38 1:1064620  (REF=C, ALT=T → ALID 1:1064620:C:T, flip=True  → z=-z)
  hg19 1:1500000 → hg38 1:1564620  (REF=G, ALT=A → ALID 1:1564620:A:G, flip=False → z unchanged)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.query import query_store
from opengwasdb.validation import validate_store

# hg19 positions used in fixtures and their expected hg38 positions
HG19_POS_1 = 100_000   # → hg38 100000  REF=A ALT=G  (ALT>REF → flip, stored z = -z)
HG19_POS_2 = 1_000_000  # → hg38 1064620 REF=C ALT=T  (ALT>REF → flip, stored z = -z)
HG19_POS_3 = 1_500_000  # → hg38 1564620 REF=G ALT=A  (ALT<REF → no flip)

HG38_ALID_1 = "1:100000:A:G"
HG38_ALID_2 = "1:1064620:C:T"
HG38_ALID_3 = "1:1564620:A:G"


def _vcf_header(study_type: str = "Continuous") -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
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
    """Write the build manifest: trait_id/file_path/trait_name/n (also the
    Analysis Catalogue's BUILD_COLUMNS) plus the required stored_effect_scale
    (issue #17) and original_sd_method/original_sd (issue #18). `scales`/
    `sd_methods`/`sds` override the value per trait_id (defaults `"sd"` and
    `"declared_standardised"` -- i.e. no rescaling, `original_sd` blank),
    letting a test declare values that disagree with its VCF's own
    ``##SAMPLE`` header -- the manifest always wins.
    """
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


@pytest.fixture
def two_trait_store(tmp_path):
    """Store built from two VCF fixtures with three variants each."""
    vcf1 = _make_vcf(
        tmp_path,
        "trait_a",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # z=4.0, flip→-4.0
            f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",  # z=5.0, flip→-5.0
            f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",  # z=3.0, no flip
        ],
    )
    vcf2 = _make_vcf(
        tmp_path,
        "trait_b",
        [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t6.0:0.5\n",  # z=12.0, flip→-12.0
            f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t1.2:0.3\n",  # z=4.0, no flip
        ],
        # The ieu-a-7 scenario: header says Continuous, but this is really a
        # case-control trait -- the manifest (below) declares log_or, and
        # that value must win (issue #17), not this header.
        study_type="Continuous",
    )
    manifest = _make_manifest(
        tmp_path,
        [("trait_a", vcf1, "Trait A"), ("trait_b", vcf2, "Trait B")],
        scales={"trait_b": "log_or"},
    )
    store_path = tmp_path / "store.opengwasdb"
    build_dense_from_vcf_manifest(
        manifest,
        store_path,
        store_id="test-store",
        release_id="v1",
    )
    return store_path


def test_build_creates_standard_store_envelope(two_trait_store):
    assert (two_trait_store / "manifest.json").exists()
    assert (two_trait_store / "index.sqlite").exists()
    assert (two_trait_store / "data.zarr").exists()
    assert (two_trait_store / "variants.tsv.gz").exists()
    assert (two_trait_store / "variant_offsets.npy").exists()


def test_validate_store_passes(two_trait_store):
    result = validate_store(two_trait_store)
    assert result.ok, result.errors


def test_manifest_json_has_grch38_assembly(two_trait_store):
    manifest = json.loads((two_trait_store / "manifest.json").read_text())
    assert manifest["reference_assembly"] == "GRCh38"
    assert manifest["completion_state"] == "observed_only"


def test_store_has_correct_dimensions(two_trait_store):
    import zarr

    root = zarr.open_group(str(two_trait_store / "data.zarr"), mode="r")
    assert root["z"].shape == (3, 2)
    assert root["se"].shape == (3, 2)


def test_allele_flip_z_negated_when_alt_not_a1(two_trait_store):
    """Variants where ALT > REF (A1=REF) should have z negated."""
    query = query_store(two_trait_store)
    # Use lookup to get trait_a's z for ALID_1 directly
    result = query.lookup([HG38_ALID_1], ["trait_a"])
    assert len(result["z"]) == 1
    # ALT=G > REF=A → z was negated; ES=2.0/SE=0.5=4.0 → stored z=-4.0
    assert result["z"][0] == pytest.approx(-4.0, rel=5e-3)


def test_z_not_negated_when_alt_is_a1(two_trait_store):
    """Variants where ALT < REF (A1=ALT) should preserve z sign."""
    query = query_store(two_trait_store)
    result = query.lookup([HG38_ALID_3], ["trait_a"])
    assert len(result["z"]) == 1
    # ALT=A < REF=G → A is A1, no flip; ES=0.6/SE=0.2=3.0 → stored z=3.0
    assert result["z"][0] == pytest.approx(3.0, rel=5e-3)


def test_missing_cells_are_absent(two_trait_store):
    """trait_b does not have a value for HG38_ALID_2; only one analysis returned."""
    query = query_store(two_trait_store)
    result = query.phewas(HG38_ALID_2)
    # Only trait_a has data for variant at HG38_ALID_2
    assert len(result["z"]) == 1
    analyses = query.analyses_table()
    trait_b_idx = next(k for k, v in analyses.items() if v["analysis_id"] == "trait_b")
    assert trait_b_idx not in result["analysis_index"].tolist()


def test_range_query_returns_expected_variants(two_trait_store):
    query = query_store(two_trait_store)
    result = query.range_phewas("1", 50_000, 200_000)
    variants = query.variants_table()
    alids = {variants[int(vi)]["alid"] for vi in result["variant_index"]}
    assert HG38_ALID_1 in alids


def test_analysis_query_returns_all_variants_for_trait(two_trait_store):
    query = query_store(two_trait_store)
    result = query.analysis("trait_a")
    assert len(result["z"]) == 3
    assert all(np.isfinite(result["z"]))


def test_stored_effect_scale_comes_from_manifest_not_header(two_trait_store):
    """The ieu-a-7 fix (issue #17): trait_b's VCF header says
    ``StudyType=Continuous``, but its manifest row declares
    ``stored_effect_scale=log_or`` -- the built store must record the
    manifest's value, not the header's."""
    query = query_store(two_trait_store)
    analyses = query.analyses_table()
    by_id = {v["analysis_id"]: v for v in analyses.values()}
    assert by_id["trait_a"]["stored_effect_scale"] == "sd"
    assert by_id["trait_b"]["stored_effect_scale"] == "log_or"


def test_missing_required_manifest_field_fails_the_build_loudly(tmp_path):
    """A manifest missing stored_effect_scale must fail the build with a
    clear error before any I/O, not fall back to VCF-header inference or a
    silent default (issue #17)."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest_path = tmp_path / "manifest.tsv"
    # No stored_effect_scale column at all.
    manifest_path.write_text(
        "trait_id\tfile_path\ttrait_name\tn\n"
        f"trait_a\t{vcf}\tTrait A\t1000\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="stored_effect_scale"):
        build_dense_from_vcf_manifest(
            manifest_path, tmp_path / "store.opengwasdb", store_id="s", release_id="r"
        )


def test_missing_original_sd_method_fails_the_build_loudly(tmp_path):
    """A manifest missing original_sd_method must fail the build the same way
    a missing stored_effect_scale does (issue #18)."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale\n"
        f"trait_a\t{vcf}\tTrait A\t1000\tsd\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="original_sd_method"):
        build_dense_from_vcf_manifest(
            manifest_path, tmp_path / "store.opengwasdb", store_id="s", release_id="r"
        )


def test_continuous_trait_rescaled_by_manifest_original_sd(tmp_path):
    """issue #18 AC1: a continuous-trait Analysis with a manifest-supplied
    original_sd != 1 has its se divided by that SD in the built store; z is
    unchanged (z = beta/se is invariant to dividing both by the same
    constant)."""
    vcf = _make_vcf(
        tmp_path,
        "trait_a",
        [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"],  # z=4.0, flip→-4.0, se=0.5
    )
    manifest = _make_manifest(
        tmp_path,
        [("trait_a", vcf, "Trait A")],
        sd_methods={"trait_a": "source_provided"},
        sds={"trait_a": "2.0"},
    )
    store_path = tmp_path / "store.opengwasdb"
    build_dense_from_vcf_manifest(manifest, store_path, store_id="s", release_id="r")

    result = query_store(store_path).analysis("trait_a")
    assert result["z"][0] == pytest.approx(-4.0, rel=5e-3)
    assert result["se"][0] == pytest.approx(0.25, rel=5e-3)  # 0.5 / 2.0


def test_binary_trait_never_rescaled(tmp_path):
    """issue #18 AC2: original_sd_method=binary_trait is never rescaled by an
    inapplicable SD scalar, regardless of stored_effect_scale."""
    vcf = _make_vcf(
        tmp_path,
        "trait_b",
        [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"],
        study_type="CaseControl",
    )
    manifest = _make_manifest(
        tmp_path,
        [("trait_b", vcf, "Trait B")],
        scales={"trait_b": "log_or"},
        sd_methods={"trait_b": "binary_trait"},
    )
    store_path = tmp_path / "store.opengwasdb"
    build_dense_from_vcf_manifest(manifest, store_path, store_id="s", release_id="r")

    result = query_store(store_path).analysis("trait_b")
    assert result["se"][0] == pytest.approx(0.5, rel=5e-3)


def test_original_sd_method_unavailable_fails_the_build_loudly(tmp_path):
    """issue #18 AC3: original_sd_method=unavailable is handled the same way
    #17 handles any other missing required field -- flagged/failed, not
    silently assumed to be 1."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest = _make_manifest(
        tmp_path, [("trait_a", vcf, "Trait A")], sd_methods={"trait_a": "unavailable"}
    )

    with pytest.raises(ValueError, match="unavailable"):
        build_dense_from_vcf_manifest(
            manifest, tmp_path / "store.opengwasdb", store_id="s", release_id="r"
        )


def test_sd_rescale_method_without_original_sd_fails_the_build_loudly(tmp_path):
    """A method that carries an SD magnitude (e.g. source_provided) but no
    usable original_sd value must fail loudly rather than silently skip
    rescaling (issue #18)."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest = _make_manifest(
        tmp_path, [("trait_a", vcf, "Trait A")], sd_methods={"trait_a": "source_provided"}
    )  # original_sd left blank

    with pytest.raises(ValueError, match="original_sd"):
        build_dense_from_vcf_manifest(
            manifest, tmp_path / "store.opengwasdb", store_id="s", release_id="r"
        )


def test_stray_original_sd_without_a_rescale_method_fails_the_build_loudly(tmp_path):
    """A tier that carries no SD magnitude (declared_standardised, binary_trait)
    must reject a stray original_sd value rather than silently ignoring it --
    a manifest declaring both is self-contradictory (issue #18)."""
    rows = [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"]
    vcf = _make_vcf(tmp_path, "trait_a", rows)
    manifest = _make_manifest(
        tmp_path,
        [("trait_a", vcf, "Trait A")],
        sd_methods={"trait_a": "declared_standardised"},
        sds={"trait_a": "1.5"},
    )

    with pytest.raises(ValueError, match="original_sd"):
        build_dense_from_vcf_manifest(
            manifest, tmp_path / "store.opengwasdb", store_id="s", release_id="r"
        )


def test_liftover_failure_above_threshold_raises(tmp_path):
    """A manifest where all VCF positions fail liftover raises LiftoverFailureError."""
    from opengwasdb.build.liftover import LiftoverFailureError

    vcf = _make_vcf(
        tmp_path,
        "bad_trait",
        [
            "1\t200000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n",
            "1\t300000\t.\tC\tT\t.\tPASS\t.\tES:SE\t0.5:0.2\n",
        ],
    )
    manifest = _make_manifest(tmp_path, [("bad_trait", vcf, "Bad Trait")])

    with pytest.raises(LiftoverFailureError):
        build_dense_from_vcf_manifest(
            manifest,
            tmp_path / "store.opengwasdb",
            store_id="s",
            release_id="r",
            liftover_failure_threshold=0.01,
        )


class TestParallel:
    def test_two_workers_matches_serial(self, tmp_path):
        vcf1 = _make_vcf(
            tmp_path,
            "trait_a",
            [
                f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",
                f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",
                f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",
            ],
        )
        vcf2 = _make_vcf(
            tmp_path,
            "trait_b",
            [
                f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t6.0:0.5\n",
                f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t1.2:0.3\n",
            ],
            study_type="CaseControl",
        )
        manifest = _make_manifest(
            tmp_path,
            [("trait_a", vcf1, "Trait A"), ("trait_b", vcf2, "Trait B")],
        )

        serial_path = tmp_path / "serial.opengwasdb"
        build_dense_from_vcf_manifest(
            manifest, serial_path, store_id="s", release_id="r", n_workers=1
        )
        parallel_path = tmp_path / "parallel.opengwasdb"
        build_dense_from_vcf_manifest(
            manifest, parallel_path, store_id="s", release_id="r", n_workers=2
        )

        assert validate_store(parallel_path).ok

        import zarr

        serial_root = zarr.open_group(str(serial_path / "data.zarr"), mode="r")
        parallel_root = zarr.open_group(str(parallel_path / "data.zarr"), mode="r")
        serial_z = serial_root["z"][:]
        parallel_z = parallel_root["z"][:]
        serial_se = serial_root["se"][:]
        parallel_se = parallel_root["se"][:]

        assert serial_z.shape == parallel_z.shape
        np.testing.assert_array_equal(np.isnan(serial_z), np.isnan(parallel_z))
        np.testing.assert_allclose(
            serial_z[~np.isnan(serial_z)], parallel_z[~np.isnan(parallel_z)]
        )
        np.testing.assert_allclose(
            serial_se[~np.isnan(serial_se)], parallel_se[~np.isnan(parallel_se)]
        )

        # Inline-harvested top hits must match between serial and parallel.
        for key in ("p_5e_04", "p_5e_06"):
            s = serial_root[f"top_hits/{key}"]
            p = parallel_root[f"top_hits/{key}"]
            np.testing.assert_array_equal(s["variant_index"][:], p["variant_index"][:])
            np.testing.assert_array_equal(s["analysis_index"][:], p["analysis_index"][:])
            np.testing.assert_array_equal(s["z"][:], p["z"][:])


class TestTopHitHarvest:
    def test_harvest_matches_full_scan(self, tmp_path):
        """Top hits harvested during Pass 2 must equal a full-matrix rescan."""
        import zarr

        from opengwasdb.layouts.dense.top_hits import (
            build_top_hit_indexes,
            threshold_key,
        )

        # trait z-scores: 4.0 and 5.0 clear the loosest tier; 3.0 does not.
        vcf1 = _make_vcf(
            tmp_path,
            "trait_a",
            [
                f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",   # z=4.0
                f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",   # z=5.0
                f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",   # z=3.0
            ],
        )
        vcf2 = _make_vcf(
            tmp_path,
            "trait_b",
            [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t6.0:0.5\n"],     # z=12.0
        )
        manifest = _make_manifest(
            tmp_path, [("trait_a", vcf1, "Trait A"), ("trait_b", vcf2, "Trait B")]
        )
        store = tmp_path / "store.opengwasdb"
        build_dense_from_vcf_manifest(manifest, store, store_id="s", release_id="r", n_workers=2)

        root = zarr.open_group(str(store / "data.zarr"), mode="r")
        harvested = {
            t: root[f"top_hits/{threshold_key(t)}"]["z"][:]
            for t in (5e-4, 5e-6, 5e-8)
        }

        # Rebuild the same index by full-matrix scan and compare.
        build_top_hit_indexes(store)
        for t in (5e-4, 5e-6, 5e-8):
            rescanned = root[f"top_hits/{threshold_key(t)}"]["z"][:]
            np.testing.assert_array_equal(harvested[t], rescanned)

        # Sanity: the loosest tier caught the three |z|>=3.4808 cells
        # (z = -12, -5, -4); the z=3.0 cell is below the 3.4808 cutoff.
        assert sorted(harvested[5e-4].tolist()) == [-12.0, -5.0, -4.0]

    def test_index_z_equals_stored_matrix(self, tmp_path):
        """Issue 046: the top-hit index z must equal the stored (float16) matrix
        value exactly, so the index agrees with what a query reads from `z`, and
        the store validates cleanly."""
        import zarr

        from opengwasdb.layouts.dense.top_hits import threshold_key

        vcf = _make_vcf(
            tmp_path,
            "trait_a",
            [
                f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # z=4.0
                f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",  # z=5.0
                f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",  # z=3.0
            ],
        )
        manifest = _make_manifest(tmp_path, [("trait_a", vcf, "Trait A")])
        store = tmp_path / "store.opengwasdb"
        build_dense_from_vcf_manifest(manifest, store, store_id="s", release_id="r", n_workers=2)

        assert validate_store(store).ok

        root = zarr.open_group(str(store / "data.zarr"), mode="r")
        z_matrix = root["z"][:]
        for t in (5e-4, 5e-6, 5e-8):
            g = root[f"top_hits/{threshold_key(t)}"]
            assert "imputed" not in g
            rows = g["variant_index"][:]
            cols = g["analysis_index"][:]
            index_z = g["z"][:]
            gathered = z_matrix[rows, cols].astype("float32")
            np.testing.assert_array_equal(index_z, gathered)


class TestBandStreaming:
    def _three_trait_manifest(self, tmp_path):
        rows = [
            f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",
            f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",
            f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",
        ]
        entries = []
        for k in range(3):
            vcf = _make_vcf(tmp_path, f"trait_{k}", rows)
            entries.append((f"trait_{k}", vcf, f"Trait {k}"))
        return _make_manifest(tmp_path, entries)

    def test_short_final_band_matches_single_band(self, tmp_path):
        """A 2-wide analysis chunk over 3 analyses (bands [0:2],[2:3]) must equal
        a single-band build byte-for-byte — exercises the short final band."""
        import zarr

        manifest = self._three_trait_manifest(tmp_path)

        single = tmp_path / "single.opengwasdb"
        build_dense_from_vcf_manifest(
            manifest, single, store_id="s", release_id="r", n_workers=2,
            chunk_shape=(1000, 1000),
        )
        banded = tmp_path / "banded.opengwasdb"
        build_dense_from_vcf_manifest(
            manifest, banded, store_id="s", release_id="r", n_workers=2,
            chunk_shape=(1000, 2),
        )

        assert validate_store(banded).ok
        rs = zarr.open_group(str(single / "data.zarr"), mode="r")
        rb = zarr.open_group(str(banded / "data.zarr"), mode="r")
        for name in ("z", "se"):
            a, b = rs[name][:], rb[name][:]
            np.testing.assert_array_equal(np.isnan(a), np.isnan(b))
            np.testing.assert_array_equal(a[~np.isnan(a)], b[~np.isnan(b)])
        # banded store really used a 2-wide analysis chunk
        assert rb["z"].chunks[1] == 2


class TestForkSafeLookup:
    def test_last_wins_dedup_on_collision(self, tmp_path):
        """Two source variants mapping to the same row (as a liftover collision
        would) resolve to one row, keeping the last stream occurrence."""
        from opengwasdb.layouts.dense.build_vcf import (
            _build_variant_key_index,
            _resolve_column,
        )

        # Both hg19 keys map to the same hg38 ALID → same row 0.
        hg19_lookup = {("1", 100, "A", "G"): "1:100:A:G", ("1", 200, "A", "G"): "1:100:A:G"}
        keys, rows = _build_variant_key_index(hg19_lookup, {"1:100:A:G": 0})

        vcf = _make_vcf(
            tmp_path,
            "t",
            [
                "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n",  # z=2.0, flip → -2.0
                "1\t200\t.\tA\tG\t.\tPASS\t.\tES:SE\t3.0:0.5\n",  # z=6.0, flip → -6.0 (later)
            ],
        )
        r, z, _se = _resolve_column(str(vcf), keys, rows)
        assert r.tolist() == [0]
        assert z[0] == pytest.approx(-6.0, rel=5e-3)  # last occurrence wins

    def test_absent_variant_not_mismapped(self, tmp_path):
        """A variant not in the panel must be dropped, not snapped to a neighbour."""
        from opengwasdb.layouts.dense.build_vcf import (
            _build_variant_key_index,
            _resolve_column,
        )

        hg19_lookup = {("1", 100, "A", "G"): "1:100:A:G", ("1", 300, "A", "G"): "1:300:A:G"}
        keys, rows = _build_variant_key_index(
            hg19_lookup, {"1:100:A:G": 0, "1:300:A:G": 1}
        )
        vcf = _make_vcf(
            tmp_path,
            "t",
            [
                "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n",  # in panel → row 0
                "1\t200\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # ABSENT → dropped
            ],
        )
        r, _z, _se = _resolve_column(str(vcf), keys, rows)
        assert r.tolist() == [0]  # only the in-panel variant, no mis-map to row 1

    def test_batched_matches_whole_file(self, tmp_path, monkeypatch):
        """A tiny batch size (many batches + a cross-batch collision) resolves to
        the same column as processing the whole file at once."""
        import opengwasdb.layouts.dense.build_vcf as bv

        # pos 100 and 200 collide onto row 0; pos 300 -> row 1.
        hg19_lookup = {
            ("1", 100, "A", "G"): "1:100:A:G",
            ("1", 200, "A", "G"): "1:100:A:G",
            ("1", 300, "C", "T"): "1:300:C:T",
        }
        keys, rows = bv._build_variant_key_index(
            hg19_lookup, {"1:100:A:G": 0, "1:300:C:T": 1}
        )
        vcf = _make_vcf(
            tmp_path,
            "t",
            [
                "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n",  # row 0
                "1\t300\t.\tC\tT\t.\tPASS\t.\tES:SE\t2.0:0.5\n",  # row 1
                "1\t200\t.\tA\tG\t.\tPASS\t.\tES:SE\t3.0:0.5\n",  # row 0 again, last wins
            ],
        )
        whole = bv._resolve_column(str(vcf), keys, rows)
        monkeypatch.setattr(bv, "_RESOLVE_BATCH", 1)  # one association per batch
        batched = bv._resolve_column(str(vcf), keys, rows)

        for a, b in zip(whole, batched, strict=True):
            np.testing.assert_array_equal(a, b)
        r, z, _se = batched
        assert sorted(r.tolist()) == [0, 1]
        # row 0 kept the later (pos 200) occurrence: z=6.0 flipped to -6.0
        assert z[r.tolist().index(0)] == pytest.approx(-6.0, rel=5e-3)


def test_ez_preferred_over_es_se(tmp_path):
    """When EZ is present and finite, it is used instead of ES/SE."""
    vcf = _make_vcf(
        tmp_path,
        "ez_trait",
        [f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tEZ:ES:SE\t7.5:2.0:0.5\n"],
    )
    manifest = _make_manifest(tmp_path, [("ez_trait", vcf, "EZ Trait")])
    store_path = tmp_path / "store.opengwasdb"
    build_dense_from_vcf_manifest(manifest, store_path, store_id="s", release_id="r")

    query = query_store(store_path)
    result = query.analysis("ez_trait")
    assert len(result["z"]) == 1
    # EZ=7.5, ALT>REF → flip → stored z = -7.5
    assert result["z"][0] == pytest.approx(-7.5, rel=5e-3)


def test_build_honours_source_reader_capability_column(tmp_path, monkeypatch):
    """A non-GWAS-VCF, non-bcftools reader can drive a build end-to-end when a
    manifest row declares a different source_reader_capability (issue #20) --
    the builder never assumes GWAS-VCF, it resolves whatever the manifest names."""
    from opengwasdb.readers import registry as readers_registry
    from opengwasdb.readers.fake import FakeReader
    from opengwasdb.readers.interface import ReaderAssociation

    fake_capability = "opengwasdb.test-fake"

    def _fake_factory(path, stored_effect_scale):
        return FakeReader(
            associations=[
                ReaderAssociation(
                    chromosome="1",
                    position=HG19_POS_1,
                    ref="A",
                    alt="G",
                    z=-2.0,
                    se=0.5,
                    stored_effect_scale=stored_effect_scale,
                )
            ]
        )

    monkeypatch.setitem(readers_registry._READERS, fake_capability, _fake_factory)

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale"
        "\toriginal_sd_method\toriginal_sd\tsource_reader_capability\n"
        f"fake_trait\tunused-placeholder-path\tFake Trait\t1000\tsd"
        f"\tdeclared_standardised\t\t{fake_capability}\n",
        encoding="utf-8",
    )
    store_path = tmp_path / "store.opengwasdb"
    build_dense_from_vcf_manifest(manifest, store_path, store_id="s", release_id="r")

    query = query_store(store_path)
    result = query.analysis("fake_trait")
    assert len(result["z"]) == 1
    assert result["z"][0] == pytest.approx(-2.0, rel=5e-3)
