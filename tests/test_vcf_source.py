"""Tests for opengwasdb.build.vcf_source."""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import patch

import pytest

from opengwasdb.build.vcf_source import (
    read_vcf_study_type,
    stream_vcf_associations,
    stream_vcf_variants,
)
from opengwasdb.model.enums import StoredEffectScale


def _write_vcf(path: Path, body: str, study_type: str = "Continuous") -> None:
    """Write a minimal GWAS-VCF fixture to ``path``."""
    header = f"""\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">
##SAMPLE=<ID=STUDY1,StudyType={study_type}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1
"""
    path.write_text(header + body, encoding="utf-8")


def test_stream_vcf_variants_yields_biallelic_records(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(
        vcf,
        "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t0.5:0.1\n"
        "1\t200\t.\tC\tT\t.\tPASS\t.\tES:SE\t0.3:0.15\n",
    )

    variants = list(stream_vcf_variants(vcf))

    assert variants == [("1", 100, "A", "G"), ("1", 200, "C", "T")]


def test_stream_vcf_variants_skips_multi_allelic(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(
        vcf,
        "1\t100\t.\tA\tG,T\t.\tPASS\t.\tES:SE\t0.5:0.1\n"
        "1\t200\t.\tC\tT\t.\tPASS\t.\tES:SE\t0.3:0.15\n",
    )

    variants = list(stream_vcf_variants(vcf))

    assert len(variants) == 1
    assert variants[0] == ("1", 200, "C", "T")


def test_stream_vcf_variants_normalises_chr_prefix(tmp_path):
    vcf = tmp_path / "test.vcf"
    header = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">
##SAMPLE=<ID=STUDY1,StudyType=Continuous>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1
"""
    vcf.write_text(
        header + "chr1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t0.5:0.1\n",
        encoding="utf-8",
    )

    variants = list(stream_vcf_variants(vcf))

    assert variants[0][0] == "1"


def test_stream_vcf_associations_uses_ez_when_present(tmp_path):
    vcf = tmp_path / "test.vcf"
    # REF=A, ALT=G: ALT>REF → flip; EZ=5.0 → stored z=-5.0
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tEZ:ES:SE\t5.0:0.5:0.1\n")

    assocs = list(stream_vcf_associations(vcf))

    assert len(assocs) == 1
    chrom, pos, ref, alt, z, se, scale = assocs[0]
    assert z == pytest.approx(-5.0, rel=1e-4)


def test_stream_vcf_associations_falls_back_to_es_over_se(tmp_path):
    vcf = tmp_path / "test.vcf"
    # EZ absent, fall back to ES/SE: ES=1.0/SE=0.5=2.0; REF=A,ALT=G→flip→-2.0
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n")

    assocs = list(stream_vcf_associations(vcf))

    assert len(assocs) == 1
    z = assocs[0][4]
    assert z == pytest.approx(-2.0, rel=1e-4)


def test_stream_vcf_associations_negates_z_when_alt_is_not_a1(tmp_path):
    # ALT=G > REF=A → A is A1, G is A2 (effect allele not A1) → z negated
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n")

    assocs = list(stream_vcf_associations(vcf))
    z = assocs[0][4]
    # ALT=G > REF=A, so A is A1. Effect allele is G (A2), so z is negated.
    assert z == pytest.approx(-2.0, rel=1e-4)


def test_stream_vcf_associations_no_flip_when_alt_is_a1(tmp_path):
    # ALT=A < REF=G → A is A1 (ALT is A1) → no flip
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tG\tA\t.\tPASS\t.\tES:SE\t1.0:0.5\n")

    assocs = list(stream_vcf_associations(vcf))
    z = assocs[0][4]
    assert z == pytest.approx(2.0, rel=1e-4)


def test_stream_vcf_associations_se_unchanged_after_flip(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.25\n")

    assocs = list(stream_vcf_associations(vcf))
    se = assocs[0][5]
    assert se == pytest.approx(0.25, rel=1e-4)


def test_stream_vcf_associations_skips_zero_se(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(
        vcf,
        "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.0\n"
        "1\t200\t.\tC\tT\t.\tPASS\t.\tES:SE\t0.5:0.1\n",
    )

    assocs = list(stream_vcf_associations(vcf))
    assert len(assocs) == 1
    assert assocs[0][1] == 200


def test_stream_vcf_associations_continuous_study_type(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n", study_type="Continuous")

    assocs = list(stream_vcf_associations(vcf))
    assert assocs[0][6] == StoredEffectScale.SD_UNITS


def test_stream_vcf_associations_case_control_study_type(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n", study_type="CaseControl")

    assocs = list(stream_vcf_associations(vcf))
    assert assocs[0][6] == StoredEffectScale.LOG_OR


def test_stream_vcf_associations_raises_on_unknown_study_type(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n", study_type="Unknown")

    with pytest.raises(ValueError, match="StudyType"):
        list(stream_vcf_associations(vcf))


def test_stream_vcf_associations_raises_on_missing_study_type(tmp_path):
    vcf = tmp_path / "no_sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size\">\n"
        "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1\n"
        "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="StudyType"):
        list(stream_vcf_associations(vcf))


def test_read_vcf_study_type_returns_scale(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "", study_type="CaseControl")

    assert read_vcf_study_type(vcf) == StoredEffectScale.LOG_OR


def test_all_three_functions_raise_when_bcftools_not_on_path(tmp_path):
    vcf = tmp_path / "test.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n")

    with patch("shutil.which", return_value=None):
        with pytest.raises(RuntimeError, match="bcftools"):
            list(stream_vcf_variants(vcf))

    with patch("shutil.which", return_value=None):
        with pytest.raises(RuntimeError, match="bcftools"):
            list(stream_vcf_associations(vcf))

    with patch("shutil.which", return_value=None):
        with pytest.raises(RuntimeError, match="bcftools"):
            read_vcf_study_type(vcf)
