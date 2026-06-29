from __future__ import annotations

import pytest

from opengwasdb.build.source import SourceRowError, read_normalised_associations


def test_source_reader_normalises_alid_and_flips_signed_statistics(source_path):
    records = read_normalised_associations(source_path)

    assert records[0].variant.alid == "1:100:A:G"
    assert records[0].z == 2.0
    assert records[0].se == 0.1
    assert records[1].variant.alid == "1:100:A:G"
    assert records[1].z == 6.0
    assert records[2].variant.alid == "1:200:C:T"
    assert records[2].z == -3.0
    assert records[2].se == 0.2


def test_source_reader_derives_z_from_beta_when_needed(tmp_path):
    source = tmp_path / "beta.tsv"
    source.write_text(
        "\t".join(
            [
                "analysis_id",
                "chromosome",
                "position",
                "effect_allele",
                "other_allele",
                "beta",
                "se",
            ]
        )
        + "\n"
        + "a1\t1\t10\tG\tA\t0.6\t0.2\n",
        encoding="utf-8",
    )

    [record] = read_normalised_associations(source)

    assert record.variant.alid == "1:10:A:G"
    assert record.z == pytest.approx(-3.0)
    assert record.se == 0.2


def test_source_reader_supports_gzipped_tsv(tmp_path):
    import gzip

    source = tmp_path / "source.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "analysis_id",
                    "chromosome",
                    "position",
                    "effect_allele",
                    "other_allele",
                    "z",
                    "se",
                ]
            )
            + "\n"
            + "a1\t1\t10\tA\tG\t2\t0.5\n"
        )

    [record] = read_normalised_associations(source)

    assert record.variant.alid == "1:10:A:G"
    assert record.z == 2.0


def test_source_reader_rejects_invalid_rows(tmp_path):
    source = tmp_path / "invalid.tsv"
    source.write_text(
        "\t".join(
            [
                "analysis_id",
                "chromosome",
                "position",
                "effect_allele",
                "other_allele",
                "z",
                "se",
            ]
        )
        + "\n"
        + "a1\t1\t10\tA\tA\t1\t0.1\n",
        encoding="utf-8",
    )

    with pytest.raises(SourceRowError, match="identical"):
        read_normalised_associations(source)
