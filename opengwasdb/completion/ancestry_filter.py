"""The per-Analysis ancestry-match filter for Reference Completion (ADR 0028).

An Analysis carrying an ``assigned_ancestry`` (from Catalogue-driven ancestry
assignment) should only be imputed against a matching-ancestry LD panel;
others are carried through observed-only. Store releases with no ancestry
information anywhere impute every Analysis -- the same behaviour as before
ADR 0028 existed.

That fallback is backward compatibility, not a judgement that the store is
single-ancestry: nothing here can tell "every Analysis really is EUR" apart
from "nobody ran ancestry assignment on a multi-cohort manifest". Since
completion imputes from the panel's own LD structure and EAF
(``opengwasdb.completion.block``), the second case silently applies one
ancestry's LD/EAF to Analyses it does not describe. The fallback therefore
warns rather than only logging at info level (issue #108).

Shared by Dense, Ragged, and Hybrid completion so the filter applies
uniformly regardless of which layout's entry point is called directly,
rather than only when the caller happens to be ``complete_hybrid_store``.
"""
from __future__ import annotations

import logging
from typing import Any

log = logging.getLogger(__name__)


def derive_impute_analysis_ids(
    analyses_rows: Any,
    ancestry: str,
) -> set[str] | None:
    """Return the Analysis ids to impute, or ``None`` to impute all.

    ``None`` when no row carries ``assigned_ancestry`` -- the source has no
    ancestry information to filter on. Otherwise the set of
    ``analysis_id``\\ s whose ``assigned_ancestry`` matches *ancestry*.

    Returning an empty set is meaningfully different from returning ``None``:
    it means ancestry *is* known and nothing matches this panel, so nothing
    should be imputed.
    """
    rows = list(analyses_rows)
    if not any(row.get("assigned_ancestry") for row in rows):
        log.warning(
            "No assigned_ancestry on any of %d analyses: imputing all of them against "
            "the %s panel regardless of their true ancestry. If this store spans "
            "multiple ancestries, its imputed cells will carry %s LD and allele "
            "frequencies. Run ancestry assignment (ADR 0028) to impute only "
            "matching analyses.",
            len(rows), ancestry, ancestry,
        )
        return None
    matched = {row["analysis_id"] for row in rows if row.get("assigned_ancestry") == ancestry}
    log.info(
        "Ancestry-matched completion: %d/%d analyses match panel ancestry %s",
        len(matched), len(rows), ancestry,
    )
    return matched
