"""Shared mechanics for canonical tabular summary-statistics readers."""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers.gwas_vcf import is_palindromic
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics

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


def parse_af(value: str | None) -> float | None:
    parsed = parse_finite_float(value)
    return parsed if parsed is not None and 0.0 <= parsed <= 1.0 else None


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


def stream_associations(
    rows: Iterable[TabularRow], stored_effect_scale: StoredEffectScale
) -> Iterator[ReaderAssociation]:
    for row in rows:
        if row.beta is None or row.se is None:
            continue
        z = row.beta / row.se
        if row.flipped:
            z = -z
        yield ReaderAssociation(
            chromosome=row.chromosome,
            position=row.position,
            ref=row.ref,
            alt=row.alt,
            z=z,
            se=row.se,
            stored_effect_scale=stored_effect_scale,
        )


def stream_variants(rows: Iterable[TabularRow]) -> Iterator[tuple[str, int, str, str]]:
    for row in rows:
        yield row.chromosome, row.position, row.ref, row.alt


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
