"""Source Reader Capability seam (issue #19): resolves a Source Collection's
`source_reader_capability` string to a concrete `SourceReader`.
"""

from __future__ import annotations

from opengwasdb.readers.fake import FakeReader
from opengwasdb.readers.gwas_vcf import GWAS_VCF_CAPABILITY, GwasVcfReader
from opengwasdb.readers.interface import ReaderAssociation, SiteMetrics, SourceReader
from opengwasdb.readers.registry import resolve_reader

__all__ = [
    "GWAS_VCF_CAPABILITY",
    "FakeReader",
    "GwasVcfReader",
    "ReaderAssociation",
    "SiteMetrics",
    "SourceReader",
    "resolve_reader",
]
