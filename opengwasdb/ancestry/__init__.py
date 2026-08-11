"""Ancestry assignment from summary-statistic allele frequencies (ADR 0028).

Recovers each Analysis's ancestry from its GWAS-VCF ``FORMAT/AF`` by fitting a
non-negative sum-to-one mixture (NNLS) of a fine **Ancestry Reference Panel**
(Privé 2022's UK Biobank 21-group reference), aggregating the proportions to
super-populations, and admitting a single **Assigned Ancestry** only when it
clears a multi-gate rule (proportion, margin, overlap, residual). The output is
written into the **Analysis Catalogue** (ADR 0027).

Site-level AF extraction lives in `opengwasdb.readers.gwas_vcf` (issue #21),
not here -- this package consumes it through the `SourceReader` interface
rather than owning any bcftools-specific code of its own.
"""

from __future__ import annotations

from opengwasdb.ancestry.catalogue import (
    BUILD_COLUMNS,
    CatalogueRow,
    catalogue_fieldnames,
    read_catalogue,
    write_catalogue,
)
from opengwasdb.ancestry.mixture import (
    AncestryAssignment,
    Gates,
    assign_ancestry,
    assign_from_vcf,
)
from opengwasdb.ancestry.pipeline import (
    SourceRow,
    annotate_catalogue,
    read_source_manifest,
)
from opengwasdb.ancestry.reference import AncestryReference, load_reference

__all__ = [
    "AncestryAssignment",
    "AncestryReference",
    "BUILD_COLUMNS",
    "CatalogueRow",
    "Gates",
    "SourceRow",
    "annotate_catalogue",
    "assign_ancestry",
    "assign_from_vcf",
    "catalogue_fieldnames",
    "load_reference",
    "read_catalogue",
    "read_source_manifest",
    "write_catalogue",
]
