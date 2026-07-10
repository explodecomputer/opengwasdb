"""Analysis Catalogue reader/writer (ADR 0027).

The Catalogue is a versioned TSV with one row per Analysis: the four build-manifest
columns (``trait_id``, ``file_path``, ``trait_name``, ``n``) followed by the
ancestry annotations. Because it is a **superset** of the build manifest and the
builder reads only its four columns via ``DictReader``, a build is a pure row
filter over the Catalogue with no format translation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.ancestry.mixture import AncestryAssignment

# The build-manifest columns, first and in this order (superset invariant).
BUILD_COLUMNS = ["trait_id", "file_path", "trait_name", "n"]

# Annotation columns written after the build columns and before the per-super-pop
# composition columns and the version stamps.
_ANNOTATION_COLUMNS = [
    "assigned_ancestry",
    "reported_population",
    "af_overlap",
    "nnls_residual",
    "dominant_superpop",
    "dominant_proportion",
    "runner_up_margin",
    "gate_reason",
]
_VERSION_COLUMNS = ["catalogue_version", "ancestry_reference_version"]

_UNASSIGNED = "Unassigned"


@dataclass(frozen=True)
class CatalogueRow:
    """One annotated Analysis: build columns + ancestry assignment."""

    trait_id: str
    file_path: str
    trait_name: str
    n: int
    reported_population: str
    assignment: AncestryAssignment


def _prop_column(superpop: str) -> str:
    return f"ancestry_prop_{superpop}"


def catalogue_fieldnames(superpops: list[str]) -> list[str]:
    """Full Catalogue header for the given super-population set."""
    return [
        *BUILD_COLUMNS,
        *_ANNOTATION_COLUMNS,
        *[_prop_column(sp) for sp in superpops],
        *_VERSION_COLUMNS,
    ]


def _row_to_dict(
    row: CatalogueRow,
    superpops: list[str],
    *,
    catalogue_version: str,
    ancestry_reference_version: str,
) -> dict[str, str]:
    a = row.assignment
    out: dict[str, str] = {
        "trait_id": row.trait_id,
        "file_path": row.file_path,
        "trait_name": row.trait_name,
        "n": str(row.n),
        "assigned_ancestry": a.assigned_ancestry or _UNASSIGNED,
        "reported_population": row.reported_population,
        "af_overlap": str(a.af_overlap),
        "nnls_residual": _fmt(a.residual),
        "dominant_superpop": a.dominant_superpop or "",
        "dominant_proportion": _fmt(a.dominant_proportion),
        "runner_up_margin": _fmt(a.runner_up_margin),
        "gate_reason": a.gate_reason,
        "catalogue_version": catalogue_version,
        "ancestry_reference_version": ancestry_reference_version,
    }
    for sp in superpops:
        out[_prop_column(sp)] = _fmt(a.superpop_composition.get(sp, 0.0))
    return out


def write_catalogue(
    path: str | Path,
    rows: list[CatalogueRow],
    superpops: list[str],
    *,
    catalogue_version: str,
    ancestry_reference_version: str,
) -> Path:
    """Write the Catalogue TSV (superset of the build manifest)."""
    path = Path(path)
    fieldnames = catalogue_fieldnames(superpops)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                _row_to_dict(
                    row,
                    superpops,
                    catalogue_version=catalogue_version,
                    ancestry_reference_version=ancestry_reference_version,
                )
            )
    return path


def read_catalogue(path: str | Path) -> list[dict[str, str]]:
    """Read a Catalogue TSV back as a list of raw string dict rows."""
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fmt(value: float) -> str:
    if value != value:  # NaN
        return "nan"
    return f"{value:.6g}"
