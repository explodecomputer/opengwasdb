"""Dense/Hybrid manifest Analytical Metadata passthrough (issue #86)."""

from __future__ import annotations

from pathlib import Path

import pytest

from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.layouts.hybrid.build import build_hybrid_from_vcf_manifest
from opengwasdb.model.analyses import read_analyses


def _source_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "binary.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size">\n'
        '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1\n"
        "1\t100000\t.\tA\tG\t.\tPASS\t.\tES:SE\t0.6:0.3\n",
        encoding="utf-8",
    )
    return path


def _manifest(tmp_path: Path, source: Path) -> Path:
    path = tmp_path / "manifest.tsv"
    path.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale"
        "\toriginal_sd_method\toriginal_sd\tassigned_ancestry"
        "\tancestry_assignment_method\tsample_size_kind\tsample_size_scope"
        "\tn_cases\tn_controls\toriginal_effect_scale\n"
        f"binary\t{source}\tBinary trait\t10000\tlog_or\tbinary_trait\t\tEUR"
        "\tsource_trusted_no_af\tcase_control\tanalysis_level\t1000\t9000\tlog_or\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("layout", ["dense", "hybrid"])
def test_manifest_analytical_metadata_passes_through_unchanged(tmp_path: Path, layout: str):
    source = _source_vcf(tmp_path)
    manifest = _manifest(tmp_path, source)
    store = tmp_path / f"{layout}.opengwasdb"

    if layout == "dense":
        build_dense_from_vcf_manifest(
            manifest, store, store_id="metadata-test", release_id="r1"
        )
        analyses_paths = [store / "analyses.tsv"]
    else:
        panel = tmp_path / "panel.txt"
        panel.write_text("1:100000:A:G\n", encoding="utf-8")
        build_hybrid_from_vcf_manifest(
            manifest,
            store,
            reference_panel=panel,
            store_id="metadata-test",
            release_id="r1",
        )
        analyses_paths = [store / "analyses.tsv", store / "dense" / "analyses.tsv"]

    expected = {
        "sample_size_kind": "case_control",
        "sample_size_scope": "analysis_level",
        "sample_size": "10000",
        "n_cases": "1000",
        "n_controls": "9000",
        "original_effect_scale": "log_or",
        "ancestry_assignment_method": "source_trusted_no_af",
    }
    for analyses_path in analyses_paths:
        row = read_analyses(analyses_path).rows[0]
        assert {column: row[column] for column in expected} == expected
