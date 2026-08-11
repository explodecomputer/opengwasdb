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


def _guide_section(content: str) -> str:
    return content.split('id="tab-guide"')[1]


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
    assert 'id="search-analyses"' in content
    assert "sorted-asc" in content
    assert "addEventListener" in content


def _ancestry_table() -> AnalysesTable:
    return AnalysesTable(
        fieldnames=(
            "analysis_index",
            "analysis_id",
            "assigned_ancestry",
            "ancestry_prop_AFR",
            "ancestry_prop_EUR",
        ),
        rows=(
            {
                "analysis_index": "0",
                "analysis_id": "a1",
                "assigned_ancestry": "EUR",
                "ancestry_prop_AFR": "0.02",
                "ancestry_prop_EUR": "0.95",
            },
            {
                "analysis_index": "1",
                "analysis_id": "a2",
                "assigned_ancestry": "",
                "ancestry_prop_AFR": "",
                "ancestry_prop_EUR": "",
            },
        ),
    )


def test_overview_html_has_analyses_and_ancestry_tabs(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _ancestry_table())
    content = out.read_text(encoding="utf-8")
    assert '<button data-tab="analyses" class="active">Analyses</button>' in content
    assert '<button data-tab="ancestry">Ancestry</button>' in content
    assert 'id="tab-analyses" class="tab-panel"' in content
    assert 'id="tab-ancestry" class="tab-panel hidden"' in content
    # Client-side tab switching only -- no reload, no fetch/XHR.
    assert "fetch(" not in content
    assert "location.href" not in content


def test_ancestry_tab_scoped_to_ancestry_columns_with_inline_bars(tmp_path):
    _write_manifest(tmp_path)
    out = write_overview_html(tmp_path, _ancestry_table())
    content = out.read_text(encoding="utf-8")
    ancestry_section = content.split('id="tab-ancestry"')[1]
    header_row = ancestry_section.split("<thead>")[1].split("</thead>")[0]
    # Scoped to analysis_id/assigned_ancestry/ancestry_prop_* -- not every
    # analyses.tsv column (e.g. analysis_index is absent from this tab).
    assert "analysis_id" in header_row
    assert "assigned_ancestry" in header_row
    assert "ancestry_prop_AFR" in header_row
    assert "ancestry_prop_EUR" in header_row
    assert "analysis_index" not in header_row
    # a1's proportions render as an inline width-bar, not a bare number.
    assert '<span class="bar-wrap"><span class="bar-fill" style="width:95.0%">' in ancestry_section
    assert '<span class="bar-wrap"><span class="bar-fill" style="width:2.0%">' in ancestry_section
    # a2 has no composition data -- still the muted blank, not a zero-width bar.
    assert '<td class="blank">—</td>' in ancestry_section


def test_guide_tab_lists_files_actually_present(tmp_path):
    _write_manifest(tmp_path)
    (tmp_path / "index.sqlite").write_bytes(b"")
    (tmp_path / "data.zarr").mkdir()
    (tmp_path / "analyses.tsv").write_text("", encoding="utf-8")

    out = write_overview_html(tmp_path, _table())
    guide = _guide_section(out.read_text(encoding="utf-8"))

    assert 'class="fname">manifest.json<span class="badge human">human</span>' in guide
    assert 'class="fname">index.sqlite<span class="badge internal">internal</span>' in guide
    assert 'class="fname">data.zarr<span class="badge internal">internal</span>' in guide
    # overview.html describes itself even though it doesn't exist on disk yet
    # at scan time (it's being written by this very call) -- it must not be
    # silently missing from its own Guide tab.
    assert 'class="fname">overview.html<span class="badge human">human</span>' in guide


def test_guide_tab_falls_back_for_unrecognised_files(tmp_path):
    _write_manifest(tmp_path)
    (tmp_path / "some_future_file.bin").write_bytes(b"")

    out = write_overview_html(tmp_path, _table())
    guide = _guide_section(out.read_text(encoding="utf-8"))

    assert "some_future_file.bin" in guide
    assert 'class="fname">some_future_file.bin<span class="badge internal">internal</span>' in guide
    assert "Internal store file." in guide


def test_guide_tab_differs_between_dense_and_hybrid_directory_contents(tmp_path):
    dense_dir = tmp_path / "dense"
    dense_dir.mkdir()
    _write_manifest(dense_dir)
    dense_out = write_overview_html(dense_dir, _table())
    dense_guide = _guide_section(dense_out.read_text(encoding="utf-8"))
    assert 'class="fname">dense<' not in dense_guide

    hybrid_dir = tmp_path / "hybrid"
    hybrid_dir.mkdir()
    (hybrid_dir / "dense").mkdir()
    _write_manifest(hybrid_dir)
    hybrid_out = write_overview_html(hybrid_dir, _table())
    hybrid_guide = _guide_section(hybrid_out.read_text(encoding="utf-8"))
    assert 'class="fname">dense<span class="badge internal">internal</span>' in hybrid_guide
