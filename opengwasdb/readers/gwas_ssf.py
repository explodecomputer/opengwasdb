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
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics
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

GWAS_SSF_CAPABILITY = "opengwasdb.gwas-ssf"

def _iter_rows(path: str | Path) -> Iterator[TabularRow]:
    """Parse each row of a filtered/harmonised GWAS-SSF file once.

    A row with an unparseable chromosome, position, or allele pair cannot be
    represented as a variant at all and is dropped from every stream; a row
    with a valid identity but an unusable `beta`/`standard_error` still
    yields a `TabularRow` (its `beta`/`se` are `None`) so `stream_variants`
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
            se = parse_positive_float(row.get("standard_error"))
            beta = parse_finite_float(row.get("beta"))
            yield TabularRow(
                chromosome=ori.variant.chromosome,
                position=ori.variant.position,
                alid=ori.variant.alid,
                ref=other_allele,
                alt=effect_allele,
                flipped=ori.flipped,
                beta=beta,
                se=se,
                af_alt=parse_af(row.get("effect_allele_frequency")),
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
        yield from stream_associations(_iter_rows(self.path), self.stored_effect_scale)

    def stream_variants(self) -> Iterator[tuple[str, int, str, str]]:
        yield from stream_variants(_iter_rows(self.path))

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        return extract_at_sites(_iter_rows(self.path), alids)
