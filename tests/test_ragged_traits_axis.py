"""Tests for traits.tsv.gz writer and reader (issue 034)."""

import pytest

from opengwasdb.traits.axis import TraitRecord, TraitsAxisReader, write_traits_axis


def _make_records() -> list[TraitRecord]:
    return [
        TraitRecord(
            analysis_index=0, analysis_id="ENSG00000000003::Whole_Blood",
            probe_id="ENSG00000000003", n=500,
            probe_chr="1", probe_bp=1_500_000,
            gene_id="ENSG00000000003", gene_name="TSPAN6",
            tissue="Whole_Blood", context=None,
        ),
        TraitRecord(
            analysis_index=1, analysis_id="ENSG00000000005::Whole_Blood",
            probe_id="ENSG00000000005", n=500,
            probe_chr="1", probe_bp=1_600_000,
            gene_id="ENSG00000000005", gene_name="TNMD",
            tissue="Whole_Blood", context=None,
        ),
        TraitRecord(
            analysis_index=2, analysis_id="ENSG00000000419::Whole_Blood",
            probe_id="ENSG00000000419", n=500,
            probe_chr="2", probe_bp=3_000_000,
            gene_id="ENSG00000000419", gene_name="DPM1",
            tissue="Whole_Blood", context=None,
        ),
        TraitRecord(
            analysis_index=3, analysis_id="ENSG00000000457::Adipose",
            probe_id="ENSG00000000457", n=350,
            probe_chr="2", probe_bp=4_000_000,
            gene_id="ENSG00000000457", gene_name="SCYL3",
            tissue="Adipose_Subcutaneous", context="sex_combined",
        ),
        TraitRecord(
            analysis_index=4, analysis_id="cg00000029::Whole_Blood",
            probe_id="cg00000029", n=800,
            probe_chr="2", probe_bp=5_000_000,
            gene_id=None, gene_name=None,
            tissue="Whole_Blood", context=None,
        ),
    ]


def test_write_and_tabix_range(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    assert (tmp_path / "traits.tsv.gz").exists()
    assert (tmp_path / "traits.tsv.gz.tbi").exists()

    reader = TraitsAxisReader(tmp_path)
    results = reader.range("1", 1_000_000, 2_000_000)
    assert len(results) == 2
    assert {r.probe_id for r in results} == {"ENSG00000000003", "ENSG00000000005"}


def test_range_second_chromosome(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    results = reader.range("2", 2_500_000, 5_500_000)
    assert len(results) == 3
    assert {r.probe_id for r in results} == {"ENSG00000000457", "cg00000029", "ENSG00000000419"}


def test_range_empty(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    results = reader.range("3", 1, 1_000_000)
    assert results == []


def test_by_probe_id(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    results = reader.by_probe_id("ENSG00000000419")
    assert len(results) == 1
    assert results[0].analysis_index == 2
    assert results[0].gene_name == "DPM1"
    assert results[0].tissue == "Whole_Blood"


def test_optional_fields_roundtrip(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    results = reader.by_probe_id("cg00000029")
    assert len(results) == 1
    r = results[0]
    assert r.gene_id is None
    assert r.gene_name is None
    assert r.context is None
    assert r.n == 800


def test_context_field_roundtrip(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    results = reader.by_probe_id("ENSG00000000457")
    assert results[0].context == "sex_combined"


def test_all_iterator(tmp_path):
    records = _make_records()
    write_traits_axis(tmp_path, records)

    reader = TraitsAxisReader(tmp_path)
    all_records = list(reader.all())
    assert len(all_records) == 5
    # Sorted by chr/bp
    assert all_records[0].probe_chr == "1"
    assert all_records[2].probe_chr == "2"


def test_no_tabix_when_positions_missing(tmp_path):
    records = [
        TraitRecord(
            analysis_index=0, analysis_id="ukb-b-1234",
            probe_id="ukb-b-1234", n=5000,
            probe_chr=None, probe_bp=None,
            gene_id=None, gene_name=None, tissue=None, context=None,
        ),
    ]
    write_traits_axis(tmp_path, records)

    assert (tmp_path / "traits.tsv.gz").exists()
    assert not (tmp_path / "traits.tsv.gz.tbi").exists()
