"""In-memory SourceReader fake (issue #19; extended for #20).

Lets annotator/consumer tests exercise the `SourceReader` interface with no
bcftools subprocess or fixture VCF file required.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics


class FakeReader:
    """A SourceReader backed by in-memory data supplied at construction.

    `variants` defaults to the `(chromosome, position, ref, alt)` of each
    supplied association when omitted -- the common case where nothing needs
    to exercise `stream_variants`'s superset relationship to
    `stream_associations` (issue #20). Pass it explicitly to simulate a
    variant present in the source but dropped from the association stream
    (e.g. invalid SE).
    """

    def __init__(
        self,
        associations: Iterable[ReaderAssociation] = (),
        sites: dict[str, SiteMetrics] | None = None,
        variants: Iterable[tuple[str, int, str, str]] | None = None,
    ) -> None:
        self._associations = list(associations)
        self._sites = dict(sites) if sites is not None else {}
        self._variants = (
            list(variants)
            if variants is not None
            else [(a.chromosome, a.position, a.ref, a.alt) for a in self._associations]
        )

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        yield from self._associations

    def stream_variants(self) -> Iterator[tuple[str, int, str, str]]:
        yield from self._variants

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        wanted = set(alids)
        return {alid: metrics for alid, metrics in self._sites.items() if alid in wanted}
