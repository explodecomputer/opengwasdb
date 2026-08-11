"""Tests for opengwasdb/layouts/dense/overview.py's house-style redesign (issue #36)."""

from __future__ import annotations

import json
from pathlib import Path

from opengwasdb.layouts.dense.overview import write_overview_html
from opengwasdb.model.analyses import AnalysesTable


def _table() -> AnalysesTable:
    return AnalysesTable(
        fieldnames=("analysis_index", "analysis_id", "phenotype_label", "sample_size"),
        rows=(
            {
                "analysis_index": "0",
                "analysis_id": "a1",
                "phenotype_label": "Height",
                "sample_size": "1000",
            },
            {
                "analysis_index": "1",
                "analysis_id": "a2",
                "phenotype_label": "",
                "sample_size": "",
            },
        ),
    )


def _write_manifest(store_path: Path) -> None:
    manifest = {
        "store_id": "fixture-store",
        "release_id": "v1",
        "completion_state": "observed_only",
    }
    (store_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def test_overview_html_uses_house_style_palette(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    # Palette custom properties lifted from docs/opengwasdb-storage-format.html.
    assert "--paper: #f1f2ef" in content
    assert "--accent: #1f5fd1" in content


def test_overview_html_header_summarises_manifest(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    assert "fixture-store" in content
    assert "release v1" in content
    assert "observed_only" in content
    assert "2 Analyses" in content


def test_overview_html_degrades_gracefully_without_manifest(tmp_path):
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    assert "2 Analyses" in content


def test_overview_html_analysis_id_is_sticky_and_first(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    header_row = content.split("<thead>")[1].split("</thead>")[0]
    # analysis_id is displayed first (sticky identity column), ahead of
    # analysis_index even though analysis_index leads the on-disk column
    # order (store-format spec §7a) -- the .tsv order itself is untouched.
    assert header_row.index("analysis_id") < header_row.index("analysis_index")
    assert 'class="sticky-col"' in header_row


def test_overview_html_blank_cells_render_as_muted_em_dash(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    assert '<td class="blank">—</td>' in content


def test_overview_html_preserves_search_and_sort_script(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _table())
    content = out.read_text(encoding="utf-8")
    assert 'id="search"' in content
    assert "sorted-asc" in content
    assert "addEventListener" in content
