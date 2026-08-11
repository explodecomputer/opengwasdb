"""`overview.html` — a generated, human-browsable rendering of `analyses.tsv`
(store-format spec §1/§7a, ADR 0030).

Build-time-derived, like `analyses.tsv` itself: written once by the builder
(or the migration script, issue #24) and never independently authored or
edited. Self-contained (no external assets) so it opens correctly from a
downloaded store directory with no network access.
"""

from __future__ import annotations

import html
from pathlib import Path

from opengwasdb.model.analyses import AnalysesTable

_TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 1.5rem; }}
  #search {{ padding: 0.4rem; width: 24rem; max-width: 100%; margin-bottom: 0.75rem; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85rem; }}
  th, td {{
    border: 1px solid #ccc; padding: 0.3rem 0.5rem; text-align: left; white-space: nowrap;
  }}
  th {{ cursor: pointer; background: #f0f0f0; position: sticky; top: 0; }}
  th.sorted-asc::after {{ content: " \\25B2"; }}
  th.sorted-desc::after {{ content: " \\25BC"; }}
  tr:nth-child(even) {{ background: #fafafa; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p>{n_rows} Analyses.</p>
<input id="search" type="text" placeholder="Filter rows...">
<table id="analyses">
<thead><tr>{header_cells}</tr></thead>
<tbody>
{body_rows}
</tbody>
</table>
<script>
const search = document.getElementById("search");
const table = document.getElementById("analyses");
const tbody = table.tBodies[0];
const rows = Array.from(tbody.rows);

search.addEventListener("input", () => {{
  const q = search.value.toLowerCase();
  for (const row of rows) {{
    row.style.display = row.textContent.toLowerCase().includes(q) ? "" : "none";
  }}
}});

table.tHead.querySelectorAll("th").forEach((th, colIndex) => {{
  let ascending = true;
  th.addEventListener("click", () => {{
    table.tHead.querySelectorAll("th").forEach(
      (h) => h.classList.remove("sorted-asc", "sorted-desc")
    );
    const sorted = rows.slice().sort((a, b) => {{
      const av = a.cells[colIndex].textContent;
      const bv = b.cells[colIndex].textContent;
      const an = parseFloat(av), bn = parseFloat(bv);
      const cmp = (!isNaN(an) && !isNaN(bn)) ? an - bn : av.localeCompare(bv);
      return ascending ? cmp : -cmp;
    }});
    sorted.forEach((row) => tbody.appendChild(row));
    th.classList.add(ascending ? "sorted-asc" : "sorted-desc");
    ascending = !ascending;
  }});
}});
</script>
</body>
</html>
"""


def write_overview_html(
    output_path: str | Path, table: AnalysesTable, *, title: str = "OpenGWASDB Analyses"
) -> Path:
    """Render `table` (an `analyses.tsv` read) as `overview.html` at `output_path`."""
    header_cells = "".join(f"<th>{html.escape(name)}</th>" for name in table.fieldnames)
    body_rows = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(row.get(name, ''))}</td>" for name in table.fieldnames)
        + "</tr>"
        for row in table.rows
    )
    out_path = Path(output_path) / "overview.html"
    out_path.write_text(
        _TEMPLATE.format(
            title=html.escape(title),
            n_rows=len(table.rows),
            header_cells=header_cells,
            body_rows=body_rows,
        ),
        encoding="utf-8",
    )
    return out_path
