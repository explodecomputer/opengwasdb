"""GWAS-SSF SourceReader (issue #84).

Adapts the orientation and beta/se -> z parsing already proven by
`opengwasdb.layouts.ragged.build_ssf._read_filtered` to the `SourceReader`
interface (issue #19), so a filtered/harmonised GWAS-Catalog-SSF file can
route through `opengwasdb.readers.registry.resolve_reader` into the Dense
and Hybrid builders (issue #20), not only the Ragged-only path that module
serves. `stream_associations`/`stream_variants` share one row parser
(`_iter_rows`) with `extract_at_sites`, so orientation and column handling
live in exactly one place. Unlike that Ragged path, `ReaderAssociation` has
no rsid field -- rsid is not part of the `SourceReader` interface any more
than it is for `GwasVcfReader` -- so this reader does not parse `rsid`/
`variant_id` at all.

`ref`/`alt` on each `ReaderAssociation`/`stream_variants` tuple are the
source's own `other_allele`/`effect_allele` labelling (mirroring GWAS-VCF's
REF/ALT, where ALT is likewise the effect allele) -- not reordered to
canonical A1/A2, per the interface's contract.

`extract_at_sites` has no GWAS-VCF/bcftools equivalent to call into (issue
#21 built that combined AF+SE lookup around bcftools -R specifically): it
scans the file once, reading `effect_allele_frequency` where the file
carries that column and dropping requested sites the file has no AF for
(SiteMetrics never fabricates an AF), oriented to canonical A1 and excluding
palindromic (A/T, C/G) variants like `GwasVcfReader.extract_at_sites` does,
since neither this reader nor its callers have strand information to
resolve them.
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

GWAS_SSF_CAPABILITY = "opengwasdb.gwas-ssf"

_MISSING = {"", ".", "NA", "NaN", "nan", "None"}


def _opt(value: str | None) -> str | None:
    if value is None or value.strip() in _MISSING:
        return None
    return value.strip()


def _parse_finite_float(value: str | None) -> float | None:
    text = _opt(value)
    if text is None:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _parse_positive_float(value: str | None) -> float | None:
    parsed = _parse_finite_float(value)
    return parsed if parsed is not None and parsed > 0.0 else None


def _parse_af(value: str | None) -> float | None:
    parsed = _parse_finite_float(value)
    if parsed is None or not (0.0 <= parsed <= 1.0):
        return None
    return parsed


@dataclass(frozen=True)
class _ParsedRow:
    """One GWAS-SSF row's identity, source-oriented alleles, and statistics."""

    chromosome: str
    position: int
    alid: str
    ref: str
    alt: str
    flipped: bool
    beta: float | None
    se: float | None
    af: float | None


def _iter_rows(path: str | Path) -> Iterator[_ParsedRow]:
    """Parse each row of a filtered/harmonised GWAS-SSF file once.

    A row with an unparseable chromosome, position, or allele pair cannot be
    represented as a variant at all and is dropped from every stream; a row
    with a valid identity but an unusable `beta`/`standard_error` still
    yields a `_ParsedRow` (its `beta`/`se` are `None`) so `stream_variants`
    can still see it per the interface's superset contract.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            effect_allele = row.get("effect_allele")
            other_allele = row.get("other_allele")
            if effect_allele is None or other_allele is None:
                continue
            try:
                ori = orient_to_canonical(
                    row.get("chromosome", ""),
                    row.get("base_pair_location", ""),
                    effect_allele,
                    other_allele,
                )
            except VariantNormalisationError:
                continue
            se = _parse_positive_float(row.get("standard_error"))
            beta = _parse_finite_float(row.get("beta"))
            yield _ParsedRow(
                chromosome=ori.variant.chromosome,
                position=ori.variant.position,
                alid=ori.variant.alid,
                ref=other_allele,
                alt=effect_allele,
                flipped=ori.flipped,
                beta=beta,
                se=se,
                af=_parse_af(row.get("effect_allele_frequency")),
            )


@dataclass(frozen=True)
class GwasSsfReader:
    """SourceReader for one filtered/harmonised GWAS-Catalog-SSF file.

    `stored_effect_scale` is Analytical Metadata for this file's Analysis,
    resolved by the caller from the build manifest (issue #16's schema) --
    the file's own columns carry no effect-scale concept, so it is never
    derived from `path` itself.
    """

    path: str | Path
    stored_effect_scale: StoredEffectScale = StoredEffectScale.SD

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
            yield (row.chromosome, row.position, row.ref, row.alt)

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        wanted = alids if isinstance(alids, (set, frozenset, dict)) else set(alids)
        if not wanted:
            return {}
        out: dict[str, SiteMetrics] = {}
        for row in _iter_rows(self.path):
            if row.alid not in wanted or row.af is None or row.se is None:
                continue
            if is_palindromic(row.ref, row.alt):
                continue
            out[row.alid] = SiteMetrics(
                af=(1.0 - row.af) if row.flipped else row.af,
                se=row.se,
            )
        return out
