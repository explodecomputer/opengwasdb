"""Resolve index-keyed query results to human-readable rows (issue #104).

Every association-returning method on `StoreQuery`/`RaggedStoreQuery`/
`HybridStoreQuery` (`opengwasdb.query.facade`) returns index-keyed parallel
arrays (ADR-0020) -- `variant_index`/`analysis_index` rather than variant/
analysis identity. `resolve_rows()` is the join back to identity: it builds
the variant side as one lookup table via `variant_axis.by_indices()` (which
itself adaptively picks a bounded random-access resolution or, once the
result touches a large-enough share of the Store Variant Table -- e.g. a
`query-analysis` pull of most of a dense column -- one sequential scan
instead; see that method's docstring) and the analysis side via
`analyses.all()` (already fully loaded at store-open, ADR 0030, no extra
I/O), then does the row-by-row join as a cheap in-memory dict lookup on
both sides. Rows are still yielded lazily, so the CLI's CSV writer streams
rather than materialising the whole resolved table.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from opengwasdb.index import AnalysesIndex
from opengwasdb.stats import log10_p_two_sided
from opengwasdb.variants import VariantAxis


def resolve_rows(
    analyses: AnalysesIndex,
    variant_axis: VariantAxis,
    result: dict[str, np.ndarray],
) -> Iterator[dict[str, object]]:
    """Yield one human-readable row per association in `result`.

    Each row carries analysis_id/analysis_label, variant identity
    (rsid/chromosome/position/alid/alleles -- rsid is "." when the store has
    none for that variant, matching the Store Variant Table's own missing
    marker), z, se, log10_p, and association_status.
    """
    analysis_rows = analyses.all()
    variant_rows = variant_axis.by_indices(result["variant_index"])
    log10_p = log10_p_two_sided(result["z"])
    for i in range(len(result["variant_index"])):
        variant = variant_rows[int(result["variant_index"][i])]
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
