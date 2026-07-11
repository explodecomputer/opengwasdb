"""Per-Analysis ancestry provenance carried alongside a store (ADR 0027/0028).

A store inherits its Analyses' **Assigned Ancestry** from the Catalogue subset it
was built from. Rather than change the traits-axis schema (so ``build-hybrid`` /
``complete-hybrid`` stay unchanged), that provenance rides in two sidecar files in
the store directory:

  * ``ancestry.tsv`` — one row per Analysis: ``trait_id``, ``assigned_ancestry``,
    ``completed_against`` (the panel ancestry it was imputed against, or empty).
  * ``ancestry_provenance.json`` — store-level provenance: the Catalogue version,
    the subset filter that selected these Analyses, and the reference version.

Ancestry-Matched Completion (issue 066) reads ``assigned_ancestry`` to decide
which Analyses to impute and writes back ``completed_against``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ANCESTRY_SIDECAR = "ancestry.tsv"
ANCESTRY_PROVENANCE = "ancestry_provenance.json"

_SIDECAR_FIELDS = ["trait_id", "assigned_ancestry", "completed_against"]


def sidecar_path(store_path: str | Path) -> Path:
    return Path(store_path) / ANCESTRY_SIDECAR


def provenance_path(store_path: str | Path) -> Path:
    return Path(store_path) / ANCESTRY_PROVENANCE


def write_ancestry_sidecar(
    store_path: str | Path,
    analyses: list[tuple[str, str]],
    *,
    completed_against: dict[str, str] | None = None,
) -> Path:
    """Write ``ancestry.tsv`` from ``(trait_id, assigned_ancestry)`` pairs."""
    completed_against = completed_against or {}
    path = sidecar_path(store_path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_SIDECAR_FIELDS, delimiter="\t")
        writer.writeheader()
        for trait_id, ancestry in analyses:
            writer.writerow(
                {
                    "trait_id": trait_id,
                    "assigned_ancestry": ancestry,
                    "completed_against": completed_against.get(trait_id, ""),
                }
            )
    return path


def read_ancestry_sidecar(store_path: str | Path) -> list[dict[str, str]]:
    """Read ``ancestry.tsv`` rows; empty list when the sidecar is absent."""
    path = sidecar_path(store_path)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_ancestry_map(store_path: str | Path) -> dict[str, str]:
    """Return ``{trait_id: assigned_ancestry}`` from the sidecar."""
    return {r["trait_id"]: r["assigned_ancestry"] for r in read_ancestry_sidecar(store_path)}


def update_completed_against(store_path: str | Path, completed_against: dict[str, str]) -> Path:
    """Rewrite the sidecar's ``completed_against`` column for the given Analyses."""
    rows = read_ancestry_sidecar(store_path)
    if not rows:
        raise FileNotFoundError(f"no ancestry sidecar to update in {store_path}")
    for row in rows:
        if row["trait_id"] in completed_against:
            row["completed_against"] = completed_against[row["trait_id"]]
    path = sidecar_path(store_path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_SIDECAR_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_ancestry_provenance(
    store_path: str | Path,
    *,
    catalogue_version: str,
    subset_filter: str,
    ancestry_reference_version: str,
    n_analyses: int,
) -> Path:
    """Write ``ancestry_provenance.json`` linking the store to its Catalogue."""
    path = provenance_path(store_path)
    payload = {
        "catalogue_version": catalogue_version,
        "subset_filter": subset_filter,
        "ancestry_reference_version": ancestry_reference_version,
        "n_analyses": n_analyses,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_ancestry_provenance(store_path: str | Path) -> dict[str, object]:
    """Read ``ancestry_provenance.json``; empty dict when absent."""
    path = provenance_path(store_path)
    if not path.exists():
        return {}
    data: dict[str, object] = json.loads(path.read_text(encoding="utf-8"))
    return data
