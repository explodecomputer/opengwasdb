"""Resolve each Analysis's final store routing (ADR 0027/0028).

Combines the AF-recovered **Assigned Ancestry** with two fallbacks decided for the
ieu-a/ieu-b/bbj-a build:

1. **Ancestry** — if AF assignment succeeded, use it; otherwise fall back to the
   **Reported Population** mapped to a super-population. Analyses whose reported
   label maps to no super-population (e.g. ``"Mixed"``) get no routing ancestry.
2. **Coverage** — a genome-wide-store eligibility gate independent of ancestry: at
   least ``min_variants`` variants, present on all autosomes, and not concentrated
   on a single chromosome. Low-coverage Analyses keep their ancestry but are marked
   ``store_eligible = False`` (flagged, not dropped) — an array-based study images
   poorly against a genome-wide panel, but is recoverable later.

``store_eligible`` = has a routing ancestry **and** passes coverage; a store build
is then the filter ``routing_ancestry == EUR and store_eligible``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

# Free-text Reported Population → super-population. Substring-matched, lower-cased.
# Mixed / unknown map to None (no routing ancestry). Aligned to the reference's
# super-pops (AFR, AMR, EAS, EUR, MID, NAF, SAS).
_REPORTED_PATTERNS: list[tuple[str, str | None]] = [
    ("european", "EUR"),
    ("europe", "EUR"),
    ("sardinian", "EUR"),
    ("east asian", "EAS"),
    ("japanese", "EAS"),
    ("chinese", "EAS"),
    ("korean", "EAS"),
    ("african", "AFR"),
    ("hispanic", "AMR"),
    ("latin", "AMR"),
    ("south asian", "SAS"),
    ("indian", "SAS"),
    ("iranian", "MID"),
    ("middle east", "MID"),
    ("north africa", "NAF"),
    ("mixed", None),
    ("trans-ethnic", None),
    ("admixed", None),
]

DEFAULT_MIN_VARIANTS = 500_000
DEFAULT_MIN_AUTOSOMES = 22
DEFAULT_MAX_FRAC_LARGEST = 0.20

ROUTING_COLUMNS = [
    "routing_ancestry",
    "routing_source",
    "total_variants",
    "genome_distributed",
    "low_coverage",
    "store_eligible",
]


def reported_to_superpop(reported: str) -> str | None:
    """Map a free-text Reported Population to a super-population, or None."""
    text = (reported or "").strip().lower()
    if not text:
        return None
    for needle, superpop in _REPORTED_PATTERNS:
        if needle in text:
            return superpop
    return None


@dataclass(frozen=True)
class RoutingResult:
    routing_ancestry: str  # super-pop, or "" if unroutable
    routing_source: str  # "af-assigned" | "reported-fallback" | ""
    total_variants: int
    genome_distributed: bool
    low_coverage: bool
    store_eligible: bool


def resolve_routing(
    assigned_ancestry: str,
    reported_population: str,
    *,
    total_variants: int,
    n_autosomes: int,
    frac_largest_chrom: float,
    min_variants: int = DEFAULT_MIN_VARIANTS,
    min_autosomes: int = DEFAULT_MIN_AUTOSOMES,
    max_frac_largest: float = DEFAULT_MAX_FRAC_LARGEST,
) -> RoutingResult:
    """Resolve final routing ancestry + coverage eligibility for one Analysis."""
    if assigned_ancestry and assigned_ancestry != "Unassigned":
        ancestry: str | None = assigned_ancestry
        source = "af-assigned"
    else:
        ancestry = reported_to_superpop(reported_population)
        source = "reported-fallback" if ancestry else ""

    genome_distributed = n_autosomes >= min_autosomes and frac_largest_chrom <= max_frac_largest
    low_coverage = not (total_variants >= min_variants and genome_distributed)
    store_eligible = ancestry is not None and not low_coverage

    return RoutingResult(
        routing_ancestry=ancestry or "",
        routing_source=source,
        total_variants=total_variants,
        genome_distributed=genome_distributed,
        low_coverage=low_coverage,
        store_eligible=store_eligible,
    )


def read_coverage(coverage_path: str | Path) -> dict[str, tuple[int, int, float]]:
    """Read ``trait_id → (total_variants, n_autosomes, frac_largest_chrom)``."""
    out: dict[str, tuple[int, int, float]] = {}
    with open(coverage_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["trait_id"]] = (
                int(r["total_variants"]),
                int(r["n_autosomes"]),
                float(r["frac_largest_chrom"]),
            )
    return out


def finalize_catalogue(
    catalogue_path: str | Path,
    coverage_path: str | Path,
    out_path: str | Path,
    *,
    min_variants: int = DEFAULT_MIN_VARIANTS,
) -> dict[str, int]:
    """Augment an assignment Catalogue with routing + coverage columns.

    Reads the assignment Catalogue and a coverage TSV (``trait_id``,
    ``total_variants``, ``n_autosomes``, ``frac_largest_chrom``), writes a new
    Catalogue with the original columns plus ``ROUTING_COLUMNS``, and returns a
    small tally. Rows without coverage data are treated as low-coverage.
    """
    coverage = read_coverage(coverage_path)
    with open(catalogue_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        base_fields = list(reader.fieldnames or [])
        rows = list(reader)

    fieldnames = base_fields + [c for c in ROUTING_COLUMNS if c not in base_fields]
    tally = {"kept": 0, "dropped_ancestry": 0, "dropped_coverage": 0, "reported_fallback": 0}

    out_path = Path(out_path)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            tot, n_auto, frac = coverage.get(row["trait_id"], (0, 0, 1.0))
            r = resolve_routing(
                row.get("assigned_ancestry", ""),
                row.get("reported_population", ""),
                total_variants=tot,
                n_autosomes=n_auto,
                frac_largest_chrom=frac,
                min_variants=min_variants,
            )
            row.update(
                routing_ancestry=r.routing_ancestry,
                routing_source=r.routing_source,
                total_variants=str(r.total_variants),
                genome_distributed=str(r.genome_distributed),
                low_coverage=str(r.low_coverage),
                store_eligible=str(r.store_eligible),
            )
            writer.writerow(row)
            if r.store_eligible:
                tally["kept"] += 1
                if r.routing_source == "reported-fallback":
                    tally["reported_fallback"] += 1
            elif not r.routing_ancestry:
                tally["dropped_ancestry"] += 1
            else:
                tally["dropped_coverage"] += 1
    return tally
