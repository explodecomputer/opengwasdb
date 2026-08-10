"""Tests for opengwasdb.readers (issue #19): capability resolution, the
GWAS-VCF reader, the in-memory fake, and the conformance suite the two share.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers import (
    GWAS_VCF_CAPABILITY,
    FakeReader,
    GwasVcfReader,
    ReaderAssociation,
    SiteMetrics,
    resolve_reader,
)


def _write_vcf(path: Path, body: str, study_type: str = "Continuous") -> None:
    """Write a minimal GWAS-VCF fixture to ``path``, including AF alongside
    the effect/precision fields ``tests/test_vcf_source.py`` already uses."""
    header = f"""\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size relative to ALT">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size">
##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score, if used to derive EFFECT/SE">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">
##SAMPLE=<ID=STUDY1,StudyType={study_type}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1
"""
    path.write_text(header + body, encoding="utf-8")


def _bgzip_index(vcf: Path) -> Path:
    """Compress + tabix-index a plain VCF (what extract_at_sites's `-R` needs),
    matching ``tests/test_ancestry_assignment.py``'s identical helper for
    ``extract_af_at_sites`` -- both readers hit the same bcftools -R
    requirement."""
    out = vcf.with_suffix(".vcf.gz")
    subprocess.run(
        ["bcftools", "view", str(vcf), "-Oz", "-o", str(out), "--write-index=tbi"],
        check=True,
        capture_output=True,
    )
    return out


# --- resolve_reader(): the single, documented resolution point ---


def test_resolve_reader_returns_gwas_vcf_reader_for_its_capability(tmp_path):
    vcf = tmp_path / "study.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n")

    reader = resolve_reader(GWAS_VCF_CAPABILITY, vcf)

    assert isinstance(reader, GwasVcfReader)


def test_resolve_reader_rejects_unknown_capability(tmp_path):
    with pytest.raises(ValueError, match="unknown source reader capability"):
        resolve_reader("opengwasdb.some-future-format", tmp_path / "study.vcf")


# --- Shared conformance suite: same assertions against GwasVcfReader and FakeReader ---


def _gwas_vcf_reader(tmp_path: Path) -> GwasVcfReader:
    vcf = tmp_path / "study.vcf"
    _write_vcf(
        vcf,
        # REF=A, ALT=G -> A1=A, effect allele (ALT=G) is A2 -> z negated, AF flipped.
        "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE:AF\t1.0:0.5:0.30\n"
        # REF=G, ALT=A -> A1=A, effect allele (ALT=A) is A1 -> no flip.
        "1\t200\t.\tG\tA\t.\tPASS\t.\tES:SE:AF\t1.0:0.5:0.70\n",
    )
    # Bgzip+index: stream_associations works on a plain VCF, but
    # extract_at_sites's bcftools -R region lookup requires an index.
    return GwasVcfReader(_bgzip_index(vcf))


def _fake_reader() -> FakeReader:
    return FakeReader(
        associations=[
            ReaderAssociation(
                chromosome="1",
                position=100,
                ref="A",
                alt="G",
                z=-2.0,
                se=0.5,
                stored_effect_scale=StoredEffectScale.SD,
            ),
            ReaderAssociation(
                chromosome="1",
                position=200,
                ref="G",
                alt="A",
                z=2.0,
                se=0.5,
                stored_effect_scale=StoredEffectScale.SD,
            ),
        ],
        sites={
            "1:100:A:G": SiteMetrics(af=0.70, se=0.5),
            "1:200:A:G": SiteMetrics(af=0.70, se=0.5),
        },
    )


@pytest.fixture(params=["gwas_vcf", "fake"])
def reader(request, tmp_path):
    if request.param == "gwas_vcf":
        return _gwas_vcf_reader(tmp_path)
    return _fake_reader()


def _malformed_gwas_vcf_reader(tmp_path: Path) -> GwasVcfReader:
    """An unrecognised StudyType is malformed input GWAS-VCF must reject --
    unlike a non-positive SE, which the underlying bcftools stream already
    filters out silently (existing, unchanged behaviour; see
    ``test_stream_vcf_associations_skips_zero_se`` in test_vcf_source.py).
    Construction itself does not raise: GwasVcfReader is a thin path wrapper,
    so the malformed header is only discovered once streaming is attempted."""
    vcf = tmp_path / "malformed.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n", study_type="Unknown")
    return GwasVcfReader(vcf)


def _malformed_fake_reader() -> FakeReader:
    """A negative SE violates ReaderAssociation's own invariant -- raised at
    construction, which is what lets "rejection of malformed input" be a
    genuinely shared, reader-agnostic conformance check rather than one only
    a file-parsing reader can exercise."""
    return FakeReader(
        associations=[
            ReaderAssociation(
                chromosome="1",
                position=100,
                ref="A",
                alt="G",
                z=1.0,
                se=-0.5,
                stored_effect_scale=StoredEffectScale.SD,
            )
        ]
    )


@pytest.fixture(params=["gwas_vcf", "fake"])
def malformed_reader_factory(request, tmp_path):
    """A zero-arg callable that builds (and, for the fake, only then raises
    on) a reader over malformed input -- a factory, not an already-built
    reader, so the ValueError each param raises at a different point
    (construction for the fake, first-iteration for GWAS-VCF) is caught by
    the test's own `pytest.raises`, not by fixture setup."""
    if request.param == "gwas_vcf":
        return lambda: _malformed_gwas_vcf_reader(tmp_path)
    return _malformed_fake_reader


def test_reader_conformance_rejects_malformed_input(malformed_reader_factory):
    with pytest.raises(ValueError):
        list(malformed_reader_factory().stream_associations())


def test_reader_conformance_orients_associations_to_a1(reader):
    associations = list(reader.stream_associations())

    assert len(associations) == 2
    by_position = {a.position: a for a in associations}
    assert by_position[100].z == pytest.approx(-2.0)
    assert by_position[200].z == pytest.approx(2.0)


def test_reader_conformance_se_is_non_negative(reader):
    associations = list(reader.stream_associations())

    assert associations
    assert all(a.se >= 0 for a in associations)


def test_reader_conformance_extract_at_sites_returns_a1_oriented_af(reader):
    sites = reader.extract_at_sites(["1:100:A:G", "1:200:A:G"])

    assert sites["1:100:A:G"].af == pytest.approx(0.70)
    assert sites["1:200:A:G"].af == pytest.approx(0.70)
    assert all(metrics.se >= 0 for metrics in sites.values())


def test_reader_conformance_extract_at_sites_ignores_unrequested_alids(reader):
    sites = reader.extract_at_sites(["1:100:A:G"])

    assert set(sites) == {"1:100:A:G"}


# --- GWAS-VCF-specific malformed-input rejection (existing vcf_source behaviour) ---


def test_gwas_vcf_reader_rejects_unknown_study_type(tmp_path):
    vcf = tmp_path / "study.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.5\n", study_type="Unknown")

    with pytest.raises(ValueError, match="StudyType"):
        list(GwasVcfReader(vcf).stream_associations())


def test_gwas_vcf_reader_extract_at_sites_with_no_sites_returns_empty(tmp_path):
    vcf = tmp_path / "study.vcf"
    _write_vcf(vcf, "1\t100\t.\tA\tG\t.\tPASS\t.\tES:SE:AF\t1.0:0.5:0.30\n")

    assert GwasVcfReader(vcf).extract_at_sites([]) == {}


# --- FakeReader: usable with no bcftools/fixture VCF dependency ---


def test_fake_reader_needs_no_bcftools_or_fixture_vcf():
    reader = FakeReader(
        associations=[
            ReaderAssociation(
                chromosome="1",
                position=100,
                ref="A",
                alt="G",
                z=1.0,
                se=0.1,
                stored_effect_scale=StoredEffectScale.SD,
            )
        ]
    )

    associations = list(reader.stream_associations())

    assert len(associations) == 1
    assert associations[0].z == 1.0
