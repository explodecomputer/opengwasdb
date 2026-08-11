"""The Source Reader interface (issue #19; extended by issues #20 and #21).

`opengwasdb-stores` declares one `source_reader_capability` string per Source
Collection (ADR-0009, e.g. `"opengwasdb.gwas-vcf"`). This module defines the
interface that string resolves to (see `opengwasdb.readers.registry`): the
things the pipeline needs from source data -- streaming every variant
(issue #20, a builder's union-variant pass), streaming associations for the
build, and extracting allele frequency and standard error at a requested set
of sites for annotation (ancestry assignment, phenotype-SD estimation --
issue #21).

The interface is structural (`typing.Protocol`), not an ABC -- this package
has no abstract-base-class precedent elsewhere, and Protocol lets
`GwasVcfReader` and `FakeReader` satisfy it without a shared base class,
consistent with the codebase's existing preference for dataclasses and duck
typing over inheritance.

The dense and hybrid builders (issue #20) and ancestry assignment / phenotype-SD
estimation (issue #21) now resolve a reader through this interface instead of
importing a source module directly.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Protocol

import numpy as np

from opengwasdb.model.enums import StoredEffectScale


@dataclass(frozen=True)
class ReaderAssociation:
    """One source association's position, effect, and precision.

    `ref`/`alt` are the source's own allele labelling, not reordered to
    canonical A1/A2 -- `z` is already sign-corrected to the A1 = min(ref, alt)
    convention every reader in this package follows. There is no
    `analysis_id`: a source file may cover one Analysis (GWAS-VCF) or many (a
    multi-analysis tabular file), so identity assignment stays the caller's
    responsibility. `stored_effect_scale` is likewise never derived from the
    source file itself (issue #17): it is Analytical Metadata a reader
    receives from its caller (ultimately the build manifest, validated
    against issue #16's schema) and attaches to every association it yields.
    """

    chromosome: str
    position: int
    ref: str
    alt: str
    z: float
    se: float
    stored_effect_scale: StoredEffectScale

    def __post_init__(self) -> None:
        if self.se < 0:
            raise ValueError(f"se must be non-negative, got {self.se!r}")


@dataclass(frozen=True)
class SiteMetrics:
    """A1-oriented allele frequency and standard error at one canonical site."""

    af: float
    se: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.af <= 1.0:
            raise ValueError(f"af must be in [0, 1], got {self.af!r}")
        if self.se < 0:
            raise ValueError(f"se must be non-negative, got {self.se!r}")


class SourceReader(Protocol):
    """One Source Format's reader: association streaming plus per-site AF/SE
    lookup for annotation."""

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        """Yield every association in the source, oriented per ReaderAssociation."""
        ...

    def stream_variants(self) -> Iterator[tuple[str, int, str, str]]:
        """Yield ``(chromosome, position, ref, alt)`` for every biallelic variant
        in the source, independent of whether it carries a usable association.

        A superset of what :meth:`stream_associations` yields positions for --
        a record dropped there for an invalid/missing effect size or SE may
        still belong on a builder's union-variant axis (issue #20), so callers
        needing full variant coverage use this rather than filtering
        :meth:`stream_associations` themselves.
        """
        ...

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        """Return `{canonical_alid: SiteMetrics}` for the requested sites found
        in the source. A requested alid absent from the result was not found
        or was dropped (palindromic, unparseable AF/SE, ...) -- callers must
        not assume every requested alid comes back.
        """
        ...


def site_metrics_arrays(sites: dict[str, SiteMetrics]) -> tuple[np.ndarray, np.ndarray]:
    """``(se, af)`` arrays from an :meth:`SourceReader.extract_at_sites` result.

    The shape :func:`opengwasdb.build.phenotype_sd.estimate_phenotype_sd` needs
    (issue #21) -- the one adapter between the reader interface and the
    estimator, so ancestry assignment (`{alid: af}`) and SD estimation
    (`se`/`af` arrays) share the same underlying extraction call rather than
    each reshaping it independently.
    """
    se = np.fromiter((m.se for m in sites.values()), dtype=np.float64, count=len(sites))
    af = np.fromiter((m.af for m in sites.values()), dtype=np.float64, count=len(sites))
    return se, af


def af_only(sites: dict[str, SiteMetrics]) -> dict[str, float]:
    """``{alid: af}`` from an :meth:`SourceReader.extract_at_sites` result --
    the shape :func:`opengwasdb.ancestry.mixture.assign_ancestry` needs.
    Ancestry assignment has no use for `se`, but the extraction it drives is
    the same combined AF+SE call SD estimation uses (issue #21)."""
    return {alid: metrics.af for alid, metrics in sites.items()}
