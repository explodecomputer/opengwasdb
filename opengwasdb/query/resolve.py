"""Resolve index-keyed query results to human-readable rows (issue #104).

Every association-returning method on `StoreQuery`/`RaggedStoreQuery`/
`HybridStoreQuery` (`opengwasdb.query.facade`) returns index-keyed parallel
arrays (ADR-0020) -- `variant_index`/`analysis_index` rather than variant/
analysis identity. `resolve_rows()` is the join back to identity.

chromosome/position/effect_allele/other_allele/alid come from
`variant_axis.identity_by_indices()` -- an ALID *is*
`chromosome:position:effect_allele:other_allele`, and that method recovers
it from the store's already-loaded mmap'd ALID search index with zero
`variants.tsv.gz` I/O. rsid is the one identity field that index can't
produce today (it isn't part of the alid), so it's resolved via
`variant_axis.by_indices()` -- which itself adaptively picks a bounded
random-access resolution or, once the request touches a large-enough share
of the Store Variant Table, one sequential scan instead (see that method's
docstring) -- and only when `include_variant_info=True`, since that's the
one part of this join that can cost anywhere from milliseconds to tens of
seconds on a large dense query. Named for what it gates rather than for
rsid specifically: eaf is expected to join this same flag once #106 lands
(store-level effect allele frequency), since it will need the same
variants.tsv.gz-backed (or similarly non-free) lookup rsid does today. The
analysis side comes from `analyses.all()`, already fully loaded at
store-open (ADR 0030, no extra I/O). Rows are yielded lazily, so the CLI's
CSV writer streams rather than materialising the whole resolved table.
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
    *,
    include_variant_info: bool = False,
) -> Iterator[dict[str, object]]:
    """Yield one human-readable row per association in `result`.

    Each row carries analysis_id/analysis_label, variant identity
    (chromosome/position/alid/alleles, plus rsid when
    `include_variant_info=True` -- "." when the store has none for that
    variant, matching the Store Variant Table's own missing marker), z, se,
    log10_p, and association_status. `rsid` is omitted from the row
    entirely when `include_variant_info=False` (the default) rather than
    filled with a placeholder, since resolving it is the expensive part of
    this join at scale and skipping it is the point of the flag.
    """
    analysis_rows = analyses.all()
    variant_index = result["variant_index"]
    log10_p = log10_p_two_sided(result["z"])

    identity = variant_axis.identity_by_indices(variant_index)
    variant_rows = (
        variant_axis.by_indices(variant_index)
        if include_variant_info or identity is None
        else None
    )

    for i in range(len(variant_index)):
        vi = int(variant_index[i])
        if identity is not None:
            chromosome = identity["chromosome"][i]
            position = int(identity["position"][i])
            effect_allele = identity["effect_allele"][i]
            other_allele = identity["other_allele"][i]
            alid = identity["alid"][i]
        else:
            assert variant_rows is not None
            variant = variant_rows[vi]
            chromosome = variant.chromosome
            position = variant.position
            effect_allele = variant.effect_allele
            other_allele = variant.other_allele
            alid = variant.alid

        analysis = analysis_rows[int(result["analysis_index"][i])]
        row: dict[str, object] = {
            "analysis_id": analysis["analysis_id"],
            "analysis_label": analysis.get("analysis_label", ""),
            "chromosome": chromosome,
            "position": position,
            "alid": alid,
            "effect_allele": effect_allele,
            "other_allele": other_allele,
            "z": float(result["z"][i]),
            "se": float(result["se"][i]),
            "log10_p": float(log10_p[i]),
            "association_status": str(result["association_status"][i]),
        }
        if include_variant_info:
            assert variant_rows is not None
            rsid_record = variant_rows.get(vi)
            row["rsid"] = (rsid_record.rsid or ".") if rsid_record is not None else "."
        yield row
