"""Shared `Analysis` construction for the Ragged builders (BESD + SSF, issue #69).

Both ingestion paths collect the same molecular/context fields (gene, tissue,
context, genomic position, sample size) from their own differently-shaped
per-analysis record before handing them to the shared `Analysis` model --
one conversion instead of two independently typed-out copies.
"""
from __future__ import annotations

from opengwasdb.model.analyses import Analysis


def molecular_analysis(
    analysis_id: str,
    *,
    gene_id: str | None,
    gene_name: str | None,
    tissue: str | None,
    context: str | None,
    trait_chr: str | None,
    trait_bp: int | None,
    n: int | None,
    stored_effect_scale: str = "",
    assigned_ancestry: str = "",
) -> Analysis:
    return Analysis(
        analysis_id=analysis_id,
        gene_id=gene_id or "",
        gene_name=gene_name or "",
        tissue=tissue or "",
        context=context or "",
        trait_chr=trait_chr or "",
        trait_bp=str(trait_bp) if trait_bp is not None else "",
        sample_size=str(n) if n is not None else "",
        stored_effect_scale=stored_effect_scale,
        assigned_ancestry=assigned_ancestry,
    )
