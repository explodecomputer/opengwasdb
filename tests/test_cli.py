from __future__ import annotations

import json

from typer.testing import CliRunner

from opengwasdb.cli.main import app


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

    phewas = runner.invoke(app, ["query-phewas", str(store_path), "rs1"])
    assert phewas.exit_code == 0, phewas.output
    phewas_rows = json.loads(phewas.output)
    assert sorted(r["analysis_index"] for r in phewas_rows) == [0, 1]

    range_query = runner.invoke(app, ["query-range-phewas", str(store_path), "1", "150", "350"])
    assert range_query.exit_code == 0, range_query.output
    range_rows = json.loads(range_query.output)
    assert len(range_rows) == 2

    analysis = runner.invoke(app, ["query-analysis", str(store_path), "a1"])
    assert analysis.exit_code == 0, analysis.output
    assert len(json.loads(analysis.output)) == 2

    top_hits = runner.invoke(app, ["query-top-hits", str(store_path)])
    assert top_hits.exit_code == 0, top_hits.output
    assert [row["z"] for row in json.loads(top_hits.output)] == [6.0, 6.0]


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
        app, ["query-lookup", str(store), "1:1064620:C:T", "trait_a"]
    )
    assert lookup.exit_code == 0, lookup.output
    assert len(json.loads(lookup.output.strip().splitlines()[-1])) == 1
