"""In-memory SourceReader fake (issue #19).

Lets annotator/consumer tests exercise the `SourceReader` interface with no
bcftools subprocess or fixture VCF file required.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics


class FakeReader:
    """A SourceReader backed by in-memory data supplied at construction."""

    def __init__(
        self,
        associations: Iterable[ReaderAssociation] = (),
        sites: dict[str, SiteMetrics] | None = None,
    ) -> None:
        self._associations = list(associations)
        self._sites = dict(sites) if sites is not None else {}

    def stream_associations(self) -> Iterator[ReaderAssociation]:
        yield from self._associations

    def extract_at_sites(self, alids: Iterable[str]) -> dict[str, SiteMetrics]:
        wanted = set(alids)
        return {alid: metrics for alid, metrics in self._sites.items() if alid in wanted}
