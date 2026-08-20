"""FinnGen R13 tabular summary-statistics SourceReader.

FinnGen publishes one tab-delimited, bgzip-compressed file per endpoint on
GRCh38.  ``alt`` is the effect allele; the reader therefore converts
``beta / sebeta`` to the package-wide canonical A1 orientation while leaving
the source REF/ALT labels intact on streamed records.
"""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics, SourceVariant
from opengwasdb.readers.tabular import (
    TabularRow,
    extract_at_sites,
    parse_af,
    parse_finite_float,
    parse_positive_float,
    stream_associations,
    stream_variants,
)
from opengwasdb.variants.normalise import VariantNormalisationError, orient_to_canonical

FINNGEN_R13_CAPABILITY = "opengwasdb.finngen-r13"


def _first_rsid(value: str | None) -> str:
    """FinnGen's `rsids` column is comma-separated where dbSNP names one
    position more than once. The Store Variant Table has one rsid per row, so
    take the first and leave the rest unrecorded rather than inventing a
    multi-value convention no reader or query path understands (issue #109).
    """
    if not value:
        return ""
    first = value.split(",")[0].strip()
    return first if first.startswith("rs") else ""

def _iter_rows(path: str | Path) -> Iterator[TabularRow]:
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # R13 uses #chrom. Accepting chrom as well preserves compatibility
            # with older captured FinnGen releases without changing semantics.
            chromosome = row.get("#chrom", row.get("chrom", ""))
            # FinnGen's chromosome vocabulary is 1-23, where 23 is chromosome X.
            chromosome = "X" if chromosome.strip() == "23" else chromosome
            ref = row.get("ref")
            alt = row.get("alt")
            if ref is None or alt is None:
                continue
            try:
                orientation = orient_to_canonical(chromosome, row.get("pos", ""), alt, ref)
            except VariantNormalisationError:
                continue
            yield TabularRow(
                chromosome=orientation.variant.chromosome,
                position=orientation.variant.position,
                ref=ref,
                alt=alt,
                alid=orientation.variant.alid,
                flipped=orientation.flipped,
                beta=parse_finite_float(row.get("beta")),
                se=parse_positive_float(row.get("sebeta")),
                af_alt=parse_af(row.get("af_alt")),
                rsid=_first_rsid(row.get("rsids")),
            )


@dataclass(frozen=True)
class FinnGenR13Reader:
    """Reader for one FinnGen R13 endpoint summary-statistics file."""

    path: str | Path
    stored_effect_scale: StoredEffectScale = StoredEffectScale.LOG_OR

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        yield from stream_associations(_iter_rows(self.path), self.stored_effect_scale)

    def stream_variants(self) -> Iterator[SourceVariant]:
        yield from stream_variants(_iter_rows(self.path))

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        return extract_at_sites(_iter_rows(self.path), alids)
