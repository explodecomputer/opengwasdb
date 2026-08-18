"""FinnGen R13 tabular summary-statistics SourceReader.

FinnGen publishes one tab-delimited, bgzip-compressed file per endpoint on
GRCh38.  ``alt`` is the effect allele; the reader therefore converts
``beta / sebeta`` to the package-wide canonical A1 orientation while leaving
the source REF/ALT labels intact on streamed records.
"""

from __future__ import annotations

import csv
import gzip
import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers.gwas_vcf import is_palindromic
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics
from opengwasdb.variants.normalise import VariantNormalisationError, orient_to_canonical

FINNGEN_R13_CAPABILITY = "opengwasdb.finngen-r13"

_MISSING = {"", ".", "NA", "NaN", "nan", "None"}


def _finite(value: str | None) -> float | None:
    if value is None or value.strip() in _MISSING:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _positive(value: str | None) -> float | None:
    parsed = _finite(value)
    return parsed if parsed is not None and parsed > 0.0 else None


def _af(value: str | None) -> float | None:
    parsed = _finite(value)
    return parsed if parsed is not None and 0.0 <= parsed <= 1.0 else None


@dataclass(frozen=True)
class _FinnGenRow:
    chromosome: str
    position: int
    ref: str
    alt: str
    alid: str
    flipped: bool
    beta: float | None
    se: float | None
    af_alt: float | None


def _iter_rows(path: str | Path) -> Iterator[_FinnGenRow]:
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # R13 uses #chrom. Accepting chrom as well preserves compatibility
            # with older captured FinnGen releases without changing semantics.
            chromosome = row.get("#chrom", row.get("chrom", ""))
            ref = row.get("ref")
            alt = row.get("alt")
            if ref is None or alt is None:
                continue
            try:
                orientation = orient_to_canonical(chromosome, row.get("pos", ""), alt, ref)
            except VariantNormalisationError:
                continue
            yield _FinnGenRow(
                chromosome=orientation.variant.chromosome,
                position=orientation.variant.position,
                ref=ref,
                alt=alt,
                alid=orientation.variant.alid,
                flipped=orientation.flipped,
                beta=_finite(row.get("beta")),
                se=_positive(row.get("sebeta")),
                af_alt=_af(row.get("af_alt")),
            )


@dataclass(frozen=True)
class FinnGenR13Reader:
    """Reader for one FinnGen R13 endpoint summary-statistics file."""

    path: str | Path
    stored_effect_scale: StoredEffectScale = StoredEffectScale.LOG_OR

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        for row in _iter_rows(self.path):
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
                stored_effect_scale=self.stored_effect_scale,
            )

    def stream_variants(self) -> Iterator[tuple[str, int, str, str]]:
        for row in _iter_rows(self.path):
            yield row.chromosome, row.position, row.ref, row.alt

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        wanted = alids if isinstance(alids, (set, frozenset, dict)) else set(alids)
        if not wanted:
            return {}
        result: dict[str, SiteMetrics] = {}
        # A single sequential scan provides both requested metrics.
        for row in _iter_rows(self.path):
            if row.alid not in wanted or row.af_alt is None or row.se is None:
                continue
            if is_palindromic(row.ref, row.alt):
                continue
            result[row.alid] = SiteMetrics(
                af=(1.0 - row.af_alt) if row.flipped else row.af_alt,
                se=row.se,
            )
        return result
