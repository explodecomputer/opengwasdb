from __future__ import annotations

import json
import math

from typer.testing import CliRunner

from opengwasdb.cli.main import _format_p, app


def test_format_p_underflow_and_nan():
    # Ordinary p, e.g. p=0.05 -> log10(p) ~= -1.301.
    assert _format_p(math.log10(0.05)) == "0.05"
    # Past float64's own representable range (issue #104: FADS1/FADS2 reaches
    # |z|=47.8), an explicit sentinel replaces a silent 0.
    assert _format_p(-320.0) == "<1e-300"
    assert _format_p(math.nan) == "NA"


def test_cli_build_validate_info_and_query_workflow(tmp_path, source_path):
    runner = CliRunner()
    store_path = tmp_path / "cli-store.opengwasdb"

    build = runner.invoke(
        app,
        [
            "build-dense",
            str(source_path),
            str(store_path),
            "--store-id",
            "cli-fixture",
            "--release-id",
            "observed-v1",
        ],
    )
    assert build.exit_code == 0, build.output
    assert json.loads(build.output)["n_variants"] == 3

    validate = runner.invoke(app, ["validate", str(store_path)])
    assert validate.exit_code == 0, validate.output
    assert validate.output.strip() == "valid"

    info = runner.invoke(app, ["info", str(store_path)])
    assert info.exit_code == 0, info.output
    assert "store_id: cli-fixture" in info.output
    assert "primary_layout: dense" in info.output

    phewas = runner.invoke(app, ["query-phewas", str(store_path), "rs1", "--format", "json"])
    assert phewas.exit_code == 0, phewas.output
    phewas_rows = json.loads(phewas.output)
    assert sorted(r["analysis_index"] for r in phewas_rows) == [0, 1]
    # json format is unchanged -- no resolved/human-readable fields.
    assert set(phewas_rows[0]) == {"variant_index", "analysis_index", "z", "se"}

    range_query = runner.invoke(app, ["query-range-phewas", str(store_path), "1", "150", "350",
                                       "--format", "json"])
    assert range_query.exit_code == 0, range_query.output
    range_rows = json.loads(range_query.output)
    assert len(range_rows) == 2

    analysis = runner.invoke(app, ["query-analysis", str(store_path), "a1", "--format", "json"])
    assert analysis.exit_code == 0, analysis.output
    assert len(json.loads(analysis.output)) == 2

    top_hits = runner.invoke(app, ["query-top-hits", str(store_path), "--format", "json"])
    assert top_hits.exit_code == 0, top_hits.output
    assert [row["z"] for row in json.loads(top_hits.output)] == [6.0, 6.0]


def test_cli_query_defaults_to_resolved_tsv(tmp_path, source_path):
    runner = CliRunner()
    store_path = tmp_path / "cli-store.opengwasdb"
    build = runner.invoke(
        app,
        [
            "build-dense", str(source_path), str(store_path),
            "--store-id", "cli-fixture", "--release-id", "observed-v1",
        ],
    )
    assert build.exit_code == 0, build.output

    phewas = runner.invoke(app, ["query-phewas", str(store_path), "rs1"])
    assert phewas.exit_code == 0, phewas.output
    lines = phewas.output.strip("\n").split("\n")
    header, *rows = [line.split("\t") for line in lines]
    # rsid is opt-in (--variant-info), not part of the default columns
    # (issue #104 follow-up): it's the one identity field that still needs
    # a variants.tsv.gz lookup, so it's the one thing a caller pays for
    # only when they ask for it.
    assert header == [
        "analysis_id", "analysis_label", "chromosome", "position", "alid",
        "effect_allele", "other_allele", "z", "se", "p", "association_status",
    ]
    assert len(rows) == 2
    by_analysis = {row[0]: row for row in rows}
    assert set(by_analysis) == {"a1", "a2"}
    a1_row = by_analysis["a1"]
    assert a1_row[1] == "Height primary"  # analysis_label
    assert a1_row[2] == "1"  # chromosome
    assert a1_row[3] == "100"  # position
    assert a1_row[7] == "2"  # z
    assert a1_row[10] == "observed"  # association_status
    assert float(a1_row[9]) < 0.05  # p, parseable and plausible

    # An explicit --format tsv is equivalent to the default.
    explicit = runner.invoke(app, ["query-phewas", str(store_path), "rs1", "--format", "tsv"])
    assert explicit.output == phewas.output

    range_query = runner.invoke(app, ["query-range-phewas", str(store_path), "1", "1", "500"])
    assert range_query.exit_code == 0, range_query.output
    range_lines = range_query.output.strip("\n").split("\n")
    assert all(len(line.split("\t")) == 11 for line in range_lines)


def test_cli_query_variant_info_flag_adds_rsid_and_eaf_columns(tmp_path, source_path):
    runner = CliRunner()
    store_path = tmp_path / "cli-store.opengwasdb"
    build = runner.invoke(
        app,
        [
            "build-dense", str(source_path), str(store_path),
            "--store-id", "cli-fixture", "--release-id", "observed-v1",
        ],
    )
    assert build.exit_code == 0, build.output

    without_it = runner.invoke(app, ["query-phewas", str(store_path), "rs1"])
    assert without_it.exit_code == 0, without_it.output
    assert "rsid" not in without_it.output.splitlines()[0].split("\t")

    with_it = runner.invoke(app, ["query-phewas", str(store_path), "rs1", "--variant-info"])
    assert with_it.exit_code == 0, with_it.output
    header, *rows = [line.split("\t") for line in with_it.output.strip("\n").split("\n")]
    assert header == [
        "analysis_id", "analysis_label", "rsid", "chromosome", "position", "alid",
        "effect_allele", "other_allele", "z", "se", "p", "eaf", "association_status",
    ]
    assert all(row[2] == "rs1" for row in rows)
    # This fixture's source has no allele frequency, so eaf is the store's own
    # missing marker rather than a fabricated value (ADR 0036).
    assert all(row[11] == "." for row in rows)


def test_cli_regenerate_overview_rewrites_from_persisted_data_only(tmp_path, source_path):
    runner = CliRunner()
    store_path = tmp_path / "cli-store.opengwasdb"
    build = runner.invoke(
        app,
        [
            "build-dense",
            str(source_path),
            str(store_path),
            "--store-id",
            "cli-fixture",
            "--release-id",
            "observed-v1",
        ],
    )
    assert build.exit_code == 0, build.output

    index_sqlite_before = (store_path / "index.sqlite").read_bytes()
    analyses_tsv_before = (store_path / "analyses.tsv").read_text(encoding="utf-8")

    # Simulate a stale/corrupted overview.html the command must overwrite.
    (store_path / "overview.html").write_text("stale", encoding="utf-8")

    result = runner.invoke(app, ["regenerate-overview", str(store_path)])
    assert result.exit_code == 0, result.output
    assert "wrote" in result.output

    content = (store_path / "overview.html").read_text(encoding="utf-8")
    assert "cli-fixture" in content  # header reads manifest.json fresh
    assert "stale" not in content

    # No rebuild -- analyses.tsv and index.sqlite are untouched (issue #23 AC3).
    assert (store_path / "index.sqlite").read_bytes() == index_sqlite_before
    assert (store_path / "analyses.tsv").read_text(encoding="utf-8") == analyses_tsv_before


def _hybrid_vcf(tmp_path, name, rows):
    header = (
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size\">\n"
        "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error\">\n"
        "##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"Z-score\">\n"
        "##SAMPLE=<ID=S,StudyType=Continuous>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    p = tmp_path / f"{name}.vcf"
    p.write_text(header + "".join(rows), encoding="utf-8")
    return p


def test_cli_build_hybrid_validate_and_query(tmp_path):
    runner = CliRunner()
    vcf = _hybrid_vcf(
        tmp_path,
        "trait_a",
        [
            "1\t100000\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",   # on-panel
            "1\t1000000\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n",  # OFF-panel -> overflow
            "1\t1500000\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n",  # on-panel
        ],
    )
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale\toriginal_sd_method\n"
        f"trait_a\t{vcf}\tTrait A\t1000\tsd\tdeclared_standardised\n",
        encoding="utf-8",
    )
    panel = tmp_path / "panel.txt"
    panel.write_text("1:100000:A:G\n1:1564620:A:G\n", encoding="utf-8")

    store = tmp_path / "hybrid-cli.opengwasdb"
    build = runner.invoke(
        app,
        [
            "build-hybrid", str(manifest), str(store),
            "--reference-panel", str(panel),
            "--store-id", "hyb-cli", "--release-id", "v1",
        ],
    )
    assert build.exit_code == 0, build.output
    out = json.loads(build.output.strip().splitlines()[-1])
    assert out["n_panel"] == 2 and out["n_off_panel"] == 1

    validate = runner.invoke(app, ["validate", str(store)])
    assert validate.exit_code == 0, validate.output

    info = runner.invoke(app, ["info", str(store)])
    assert "primary_layout: hybrid" in info.output

    # Off-panel variant is served from the overflow.
    lookup = runner.invoke(
        app, ["query-lookup", str(store), "1:1064620:C:T", "trait_a", "--format", "json"]
    )
    assert lookup.exit_code == 0, lookup.output
    assert len(json.loads(lookup.output.strip().splitlines()[-1])) == 1
