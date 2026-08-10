"""The Source Reader interface (issue #19).

`opengwasdb-stores` declares one `source_reader_capability` string per Source
Collection (ADR-0009, e.g. `"opengwasdb.gwas-vcf"`). This module defines the
interface that string resolves to (see `opengwasdb.readers.registry`): the two
things the pipeline needs from source data, streaming associations for the
build and extracting allele frequency and standard error at a requested set
of sites for annotation (ancestry assignment, phenotype-SD estimation).

The interface is structural (`typing.Protocol`), not an ABC -- this package
has no abstract-base-class precedent elsewhere, and Protocol lets
`GwasVcfReader` and `FakeReader` satisfy it without a shared base class,
consistent with the codebase's existing preference for dataclasses and duck
typing over inheritance.

Wiring existing builders to resolve a reader through this interface instead
of importing a source module directly (`opengwasdb.build.vcf_source`,
`opengwasdb.build.source`) is out of scope here -- this ticket only
introduces the seam; no existing call site changes.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Protocol

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

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        """Return `{canonical_alid: SiteMetrics}` for the requested sites found
        in the source. A requested alid absent from the result was not found
        or was dropped (palindromic, unparseable AF/SE, ...) -- callers must
        not assume every requested alid comes back.
        """
        ...
