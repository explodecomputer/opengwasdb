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

    variant = runner.invoke(app, ["query-variant", str(store_path), "rs1"])
    assert variant.exit_code == 0, variant.output
    variant_rows = json.loads(variant.output)
    assert [row["analysis_id"] for row in variant_rows] == ["a1", "a2"]

    range_query = runner.invoke(app, ["query-range", str(store_path), "1", "150", "350"])
    assert range_query.exit_code == 0, range_query.output
    assert [row["alid"] for row in json.loads(range_query.output)] == ["1:200:C:T", "1:300:A:G"]

    analysis = runner.invoke(app, ["query-analysis", str(store_path), "a1"])
    assert analysis.exit_code == 0, analysis.output
    assert len(json.loads(analysis.output)) == 2

    phewas = runner.invoke(app, ["query-phewas", str(store_path), "1:100:A:G"])
    assert phewas.exit_code == 0, phewas.output
    assert len(json.loads(phewas.output)) == 2

    top_hits = runner.invoke(app, ["query-top-hits", str(store_path)])
    assert top_hits.exit_code == 0, top_hits.output
    assert [row["z"] for row in json.loads(top_hits.output)] == [6.0, 6.0]
