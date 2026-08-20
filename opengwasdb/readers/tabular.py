"""Shared mechanics for canonical tabular summary-statistics readers."""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers.gwas_vcf import is_palindromic
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics, SourceVariant
from opengwasdb.stats import parse_af

_MISSING = {"", ".", "NA", "NaN", "nan", "None"}


def parse_finite_float(value: str | None) -> float | None:
    if value is None or value.strip() in _MISSING:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def parse_positive_float(value: str | None) -> float | None:
    parsed = parse_finite_float(value)
    return parsed if parsed is not None and parsed > 0.0 else None


__all__ = [
    "TabularRow",
    "extract_at_sites",
    # Re-exported: the readers built on this module import their column
    # parsing from here, but the allele-frequency rule itself is shared with
    # every other source path (ADR 0036), so it lives in opengwasdb.stats.
    "parse_af",
    "parse_finite_float",
    "parse_positive_float",
    "stream_associations",
    "stream_variants",
]


@dataclass(frozen=True)
class TabularRow:
    """Format-neutral row after source-specific column mapping."""

    chromosome: str
    position: int
    ref: str
    alt: str
    alid: str
    flipped: bool
    beta: float | None
    se: float | None
    af_alt: float | None
    rsid: str = ""  # the source's own identifier for this row; "" when it names none


def stream_associations(
    rows: Iterable[TabularRow], stored_effect_scale: StoredEffectScale
) -> Iterator[ReaderAssociation]:
    for row in rows:
        if row.beta is None or row.se is None:
            continue
        z = row.beta / row.se
        eaf = row.af_alt
        if row.flipped:
            z = -z
            # `af_alt` is the frequency of the source's alt allele; flipping
            # swapped which allele is stored as the effect one (ADR 0036).
            eaf = None if eaf is None else 1.0 - eaf
        yield ReaderAssociation(
            chromosome=row.chromosome,
            position=row.position,
            ref=row.ref,
            alt=row.alt,
            z=z,
            se=row.se,
            stored_effect_scale=stored_effect_scale,
            eaf=eaf,
        )


def stream_variants(rows: Iterable[TabularRow]) -> Iterator[SourceVariant]:
    for row in rows:
        yield SourceVariant(
            chromosome=row.chromosome,
            position=row.position,
            ref=row.ref,
            alt=row.alt,
            rsid=row.rsid,
        )


def extract_at_sites(
    rows: Iterable[TabularRow], alids: Iterable[str]
) -> dict[str, SiteMetrics]:
    wanted = alids if isinstance(alids, (set, frozenset, dict)) else set(alids)
    if not wanted:
        return {}
    result: dict[str, SiteMetrics] = {}
    for row in rows:
        if row.alid not in wanted or row.af_alt is None or row.se is None:
            continue
        if is_palindromic(row.ref, row.alt):
            continue
        result[row.alid] = SiteMetrics(
            af=(1.0 - row.af_alt) if row.flipped else row.af_alt,
            se=row.se,
        )
    return result
