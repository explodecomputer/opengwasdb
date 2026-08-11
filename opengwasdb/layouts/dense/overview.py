"""`overview.html` — a generated, human-browsable rendering of a Store
Release's Analytical Metadata (store-format spec §1/§7a, ADR 0030/0032).

Build-time-derived, like `analyses.tsv` itself: written once by the builder
(or the migration script, issue #24) and never independently authored or
edited. Self-contained (no external assets, no network requests) so it opens
correctly from a downloaded store directory with no network access. Visual
identity (palette, typography) matches `docs/opengwasdb-storage-format.html`,
so OpenGWASDB's generated HTML artifacts read as one family (ADR 0032).

A single file, tab-switched client-side rather than split across sibling
files (ADR 0032): AC1/AC2 (store-format spec §20, ADR 0030) both name
`overview.html` specifically as the file a user browses Analyses in, and
tabs keep the "no external assets, no network" self-containment property
trivially true for the whole page at once (issue #37).
"""

from __future__ import annotations

import html
import json
from pathlib import Path

from opengwasdb.model.analyses import ANCESTRY_PROP_PREFIX, AnalysesTable

_STYLE = """
:root {
  --paper: #f1f2ef;
  --surface: #fbfbf9;
  --ink: #14181c;
  --ink-soft: #4b5459;
  --muted: #8b9095;
  --hairline: #dcdcd6;
  --hairline-strong: #c3c4be;
  --accent: #1f5fd1;
  --accent-soft: #dce6fa;
  --accent-wash: #eef3fc;
  --mono: ui-monospace, "SF Mono", "Cascadia Code", "Roboto Mono", Consolas, "Liberation Mono", monospace;
  --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--paper);
  color: var(--ink);
  font-family: var(--sans);
  font-size: 15px;
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
}
header {
  padding: 28px 32px 20px;
  border-bottom: 1px solid var(--hairline);
  background: var(--surface);
}
.eyebrow {
  font-family: var(--mono); font-size: 11px; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--accent); margin-bottom: 8px; display: flex; align-items: center; gap: 8px;
}
.eyebrow::before {
  content: ""; width: 7px; height: 7px; background: var(--accent); border-radius: 1px;
  display: inline-block;
}
header h1 { margin: 0; font-size: 22px; font-weight: 650; letter-spacing: -0.01em; }
header .meta { margin-top: 6px; font-size: 13px; color: var(--ink-soft); font-family: var(--mono); }
.tabs {
  display: flex; gap: 4px; padding: 0 32px; background: var(--surface);
  border-bottom: 1px solid var(--hairline);
}
.tabs button {
  font: inherit; font-size: 13px; padding: 10px 14px; border: none; background: none; cursor: pointer;
  color: var(--ink-soft); border-bottom: 2px solid transparent; margin-bottom: -1px;
}
.tabs button:hover { color: var(--ink); }
.tabs button.active { color: var(--accent); border-bottom-color: var(--accent); font-weight: 600; }
main { padding: 20px 32px 60px; }
.tab-panel.hidden { display: none; }
input[type="text"] {
  font: inherit; padding: 0.5rem 0.75rem; width: 24rem; max-width: 100%; margin-bottom: 0.9rem;
  border: 1px solid var(--hairline-strong); border-radius: 6px; background: var(--surface); color: var(--ink);
}
input[type="text"]:focus { outline: 2px solid var(--accent-soft); border-color: var(--accent); }
.table-scroll {
  overflow-x: auto; border: 1px solid var(--hairline); border-radius: 8px; background: var(--surface);
}
table { border-collapse: collapse; font-size: 12.5px; width: max-content; min-width: 100%; }
th, td {
  padding: 6px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--hairline);
}
th {
  cursor: pointer; background: var(--surface); position: sticky; top: 0; z-index: 2;
  color: var(--muted); font-weight: 600; font-size: 11px; text-transform: uppercase; letter-spacing: 0.03em;
  border-bottom: 1px solid var(--hairline-strong);
}
th.sorted-asc::after { content: " \\25B2"; color: var(--accent); }
th.sorted-desc::after { content: " \\25BC"; color: var(--accent); }
tbody tr:nth-child(even) { background: var(--paper); }
td.blank { color: var(--muted); }
th.sticky-col, td.sticky-col {
  position: sticky; left: 0; z-index: 1; background: var(--surface); box-shadow: 1px 0 0 var(--hairline-strong);
}
th.sticky-col { z-index: 3; }
tbody tr:nth-child(even) td.sticky-col { background: var(--paper); }
.bar-wrap {
  display: inline-block; width: 60px; height: 0.65em; background: var(--hairline); border-radius: 2px;
  overflow: hidden; vertical-align: middle; margin-right: 6px;
}
.bar-fill { display: block; height: 100%; background: var(--accent); }
"""

_SCRIPT = """
function initTable(table) {
  const search = document.getElementById(table.dataset.searchable);
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);

  search.addEventListener("input", () => {
    const q = search.value.toLowerCase();
    for (const row of rows) {
      row.style.display = row.textContent.toLowerCase().includes(q) ? "" : "none";
    }
  });

  table.tHead.querySelectorAll("th").forEach((th, colIndex) => {
    let ascending = true;
    th.addEventListener("click", () => {
      table.tHead.querySelectorAll("th").forEach(
        (h) => h.classList.remove("sorted-asc", "sorted-desc")
      );
      const sorted = rows.slice().sort((a, b) => {
        const av = a.cells[colIndex].textContent;
        const bv = b.cells[colIndex].textContent;
        const an = parseFloat(av), bn = parseFloat(bv);
        const cmp = (!isNaN(an) && !isNaN(bn)) ? an - bn : av.localeCompare(bv);
        return ascending ? cmp : -cmp;
      });
      sorted.forEach((row) => tbody.appendChild(row));
      th.classList.add(ascending ? "sorted-asc" : "sorted-desc");
      ascending = !ascending;
    });
  });
}
document.querySelectorAll("table[data-searchable]").forEach(initTable);

document.querySelectorAll(".tabs button").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tabs button").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".tab-panel").forEach((p) => p.classList.add("hidden"));
    btn.classList.add("active");
    document.getElementById("tab-" + btn.dataset.tab).classList.remove("hidden");
  });
});
"""

_TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{style}</style>
</head>
<body>
<header>
  <div class="eyebrow">OpenGWASDB Store Release</div>
  <h1>{title}</h1>
  <div class="meta">{meta}</div>
</header>
<nav class="tabs">
  <button data-tab="analyses" class="active">Analyses</button>
  <button data-tab="ancestry">Ancestry</button>
</nav>
<main>
<section id="tab-analyses" class="tab-panel">
{analyses_section}
</section>
<section id="tab-ancestry" class="tab-panel hidden">
{ancestry_section}
</section>
</main>
<script>{script}</script>
</body>
</html>
"""


def _read_manifest_summary(output_path: Path) -> tuple[str, str, str]:
    """`(store_id, release_id, completion_state)`, best-effort from
    `manifest.json`.

    `overview.html` is regenerable from a store's already-persisted data
    (issue #23 AC3), so a missing or unreadable manifest degrades to blank
    header fields rather than failing the whole render.
    """
    try:
        manifest = json.loads((output_path / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "", "", ""
    return (
        str(manifest.get("store_id", "")),
        str(manifest.get("release_id", "")),
        str(manifest.get("completion_state", "")),
    )


def _display_fieldnames(table: AnalysesTable) -> tuple[str, ...]:
    """`analysis_id` first (the sticky identity column), then every other
    column in `analyses.tsv`'s own order (store-format spec §7a) -- a
    display-only reordering; the on-disk column order is unchanged.
    """
    if "analysis_id" not in table.fieldnames:
        return table.fieldnames
    rest = [name for name in table.fieldnames if name != "analysis_id"]
    return ("analysis_id", *rest)


def _ancestry_prop_columns(fieldnames: tuple[str, ...]) -> list[str]:
    return sorted(name for name in fieldnames if name.startswith(ANCESTRY_PROP_PREFIX))


def _ancestry_fieldnames(table: AnalysesTable) -> tuple[str, ...]:
    """`analysis_id`, `assigned_ancestry`, then every `ancestry_prop_<population>`
    column present -- the same Ancestry Composition data already in
    `analyses.tsv`, just a narrower column set than the Analyses tab
    (issue #37)."""
    return ("analysis_id", "assigned_ancestry", *_ancestry_prop_columns(table.fieldnames))


def _bar_html(value: str) -> str:
    """An inline width-bar for a 0..1 proportion, followed by its value --
    an at-a-glance composition view with no charting library (issue #37)."""
    try:
        fraction = float(value)
    except ValueError:
        return html.escape(value)
    pct = max(0.0, min(1.0, fraction)) * 100
    return (
        f'<span class="bar-wrap"><span class="bar-fill" style="width:{pct:.1f}%"></span></span>'
        f"{html.escape(value)}"
    )


def _render_row(
    row: dict[str, str],
    fieldnames: tuple[str, ...],
    *,
    sticky_column: str | None,
    bar_columns: frozenset[str],
) -> str:
    cells = []
    for name in fieldnames:
        value = row.get(name, "")
        classes = " ".join(
            c
            for c in ("sticky-col" if name == sticky_column else "", "" if value else "blank")
            if c
        )
        class_attr = f' class="{classes}"' if classes else ""
        if value and name in bar_columns:
            body = _bar_html(value)
        else:
            body = html.escape(value) if value else "—"
        cells.append(f"<td{class_attr}>{body}</td>")
    return "<tr>" + "".join(cells) + "</tr>"


def _render_table_section(
    *,
    table_id: str,
    search_id: str,
    search_placeholder: str,
    fieldnames: tuple[str, ...],
    rows: tuple[dict[str, str], ...],
    sticky_column: str | None = None,
    bar_columns: frozenset[str] = frozenset(),
) -> str:
    header_cells = "".join(
        f'<th class="sticky-col">{html.escape(name)}</th>'
        if name == sticky_column
        else f"<th>{html.escape(name)}</th>"
        for name in fieldnames
    )
    body_rows = "\n".join(
        _render_row(row, fieldnames, sticky_column=sticky_column, bar_columns=bar_columns)
        for row in rows
    )
    return (
        f'<input id="{search_id}" type="text" placeholder="{html.escape(search_placeholder)}">\n'
        '<div class="table-scroll">\n'
        f'<table id="{table_id}" data-searchable="{search_id}">\n'
        f"<thead><tr>{header_cells}</tr></thead>\n"
        f"<tbody>\n{body_rows}\n</tbody>\n"
        "</table>\n"
        "</div>"
    )


def write_overview_html(
    output_path: str | Path, table: AnalysesTable, *, title: str = "OpenGWASDB Analyses"
) -> Path:
    """Render `table` (an `analyses.tsv` read) as `overview.html` at `output_path`."""
    output_path = Path(output_path)
    store_id, release_id, completion_state = _read_manifest_summary(output_path)
    meta_parts = [
        part
        for part in (store_id, f"release {release_id}" if release_id else "", completion_state)
        if part
    ]
    meta_parts.append(f"{len(table.rows)} Analyses")

    analyses_section = _render_table_section(
        table_id="analyses",
        search_id="search-analyses",
        search_placeholder="Filter rows...",
        fieldnames=_display_fieldnames(table),
        rows=table.rows,
        sticky_column="analysis_id",
    )
    ancestry_fieldnames = _ancestry_fieldnames(table)
    ancestry_section = _render_table_section(
        table_id="ancestry",
        search_id="search-ancestry",
        search_placeholder="Filter analyses...",
        fieldnames=ancestry_fieldnames,
        rows=table.rows,
        sticky_column="analysis_id",
        bar_columns=frozenset(_ancestry_prop_columns(table.fieldnames)),
    )

    out_path = output_path / "overview.html"
    out_path.write_text(
        _TEMPLATE.format(
            title=html.escape(title),
            style=_STYLE,
            meta=html.escape(" · ".join(meta_parts)),
            analyses_section=analyses_section,
            ancestry_section=ancestry_section,
            script=_SCRIPT,
        ),
        encoding="utf-8",
    )
    return out_path
