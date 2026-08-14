"""The ``provenance["completion"]`` sub-dict shared by Dense and Ragged
completed-manifest construction -- the parameters every Reference Completion
run records, regardless of layout. Each layout still builds its own
``StoreManifest`` (the rest of the manifest, and which extra completion
counters it reports, differ) and merges its own layout-specific fields on
top via ``extra``.

Hybrid's ``provenance["completion"]`` key is structurally different -- it
records *which component* was completed, not imputation parameters -- so it
is not built from this helper.
"""
from __future__ import annotations

from typing import Any

from opengwasdb.completion.block import COMPLETION_METHOD


def build_completion_provenance(
    *,
    ld_panel_id: str,
    ancestry: str,
    min_cor: float,
    thresh: float,
    n_variants_total: int,
    n_variants_new: int,
    method: str = COMPLETION_METHOD,
    **extra: Any,
) -> dict[str, Any]:
    return {
        "method": method,
        "ld_panel_id": ld_panel_id,
        "ancestry": ancestry,
        "min_cor": min_cor,
        "pca_thresh": thresh,
        "n_variants_total": n_variants_total,
        "n_variants_new": n_variants_new,
        **extra,
    }
