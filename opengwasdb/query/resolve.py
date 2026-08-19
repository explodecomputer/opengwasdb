"""Resolve index-keyed query results to human-readable rows (issue #104).

Every association-returning method on `StoreQuery`/`RaggedStoreQuery`/
`HybridStoreQuery` (`opengwasdb.query.facade`) returns index-keyed parallel
arrays (ADR-0020) -- `variant_index`/`analysis_index` rather than variant/
analysis identity. `resolve_rows()` joins those indices back to identity
using the store's already-open `AnalysesIndex`/`VariantAxis`: `analyses.all()`
is already fully loaded at store-open (ADR 0030, no extra I/O), and each
variant is resolved with a single bounded `variant_axis.by_index()` seek,
cached per unique index, the first time that variant is seen -- never a scan
of the whole (potentially tens-of-millions-of-rows) `variants.tsv.gz` Store
Variant Table the way an external post-processor would have to.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from opengwasdb.index import AnalysesIndex
from opengwasdb.stats import log10_p_two_sided
from opengwasdb.variants import VariantAxis, VariantRecord


def resolve_rows(
    analyses: AnalysesIndex,
    variant_axis: VariantAxis,
    result: dict[str, np.ndarray],
) -> Iterator[dict[str, object]]:
    """Yield one human-readable row per association in `result`.

    Each row carries analysis_id/analysis_label, variant identity
    (rsid/chromosome/position/alid/alleles -- rsid is "." when the store has
    none for that variant, matching the Store Variant Table's own missing
    marker), z, se, log10_p, and association_status. A generator that
    resolves lazily, one variant at a time as rows are consumed, rather than
    materialising the whole resolved table up front -- query-analysis on a
    dense store can return millions of rows.
    """
    analysis_rows = analyses.all()
    log10_p = log10_p_two_sided(result["z"])
    variant_cache: dict[int, VariantRecord | None] = {}
    for i in range(len(result["variant_index"])):
        vi = int(result["variant_index"][i])
        if vi not in variant_cache:
            variant_cache[vi] = variant_axis.by_index(vi)
        variant = variant_cache[vi]
        assert variant is not None, f"variant_index {vi} not found in Store Variant Table"
        analysis = analysis_rows[int(result["analysis_index"][i])]
        yield {
            "analysis_id": analysis["analysis_id"],
            "analysis_label": analysis.get("analysis_label", ""),
            "rsid": variant.rsid or ".",
            "chromosome": variant.chromosome,
            "position": variant.position,
            "alid": variant.alid,
            "effect_allele": variant.effect_allele,
            "other_allele": variant.other_allele,
            "z": float(result["z"][i]),
            "se": float(result["se"][i]),
            "log10_p": float(log10_p[i]),
            "association_status": str(result["association_status"][i]),
        }
