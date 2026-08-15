"""Ragged top-hit index builder: harvest candidate cells from the CSR store.

The Top-Hit Index format itself -- write/read/counts/validation -- lives in
``opengwasdb.top_hits``. This module owns only what is genuinely Ragged-specific:
deriving each association's Analysis from the CSR offsets.

Ragged stores do not yet persist Top-Hit Counts (unlike Dense/Hybrid's
``opengwasdb.layouts.dense.build.add_hit_counts``, which writes them to
``analyses.tsv``): Ragged has no ``analyses.tsv`` of its own yet, still relying
on a divergent SQLite ``analyses`` table (issue #63, a gap ADR-0030 flagged but
never resolved for Ragged). Adding Top-Hit Count columns to that SQLite table
was tried and reverted (issues #53/#61 discussion) -- it would have deepened a
schema the project has already decided to retire. Revisit once #63 lands.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.top_hits.format import TOP_HIT_THRESHOLDS, threshold_key, z_critical
from opengwasdb.top_hits.writer import write


def build_ragged_top_hit_indexes(
    store_path: str | Path,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> None:
    """Build ranked top-hit arrays for each configured p-value threshold.

    Writes to data.zarr/top_hits/<key>/ using the same schema as the dense
    builder so the query facade and validator can share one code path.
    """
    store_path = Path(store_path)
    csr = RaggedCSRReader(store_path)

    offsets = csr._offsets[:]
    vi_all = csr._variant_index[:].astype(np.int32)
    z_all = csr._z[:].astype(np.float32)
    se_all = csr._se[:].astype(np.float32)
    n_analyses = len(offsets) - 1
    imputed_all = (
        csr._root["imputed"][:].astype(np.uint8) if "imputed" in csr._root else None
    )

    # Derive analysis_index for every association via searchsorted on CSR offsets.
    # offsets[i+1] is the exclusive end of analysis i → searchsorted(offsets[1:], pos) gives i.
    positions = np.arange(len(vi_all), dtype=np.int64)
    analysis_indices = np.searchsorted(offsets[1:], positions, side="right").astype(np.int32)

    # write() re-filters per threshold, so the full (already-sparse) CSR can be
    # passed directly -- no candidate pre-filtering needed the way Dense's
    # band scan needs one to avoid holding a whole matrix in memory.
    write(
        store_path, vi_all, analysis_indices, z_all, se_all, n_analyses,
        thresholds, imputed=imputed_all,
    )

    abs_z = np.abs(z_all)
    for threshold in thresholds:
        n_hits = int(np.count_nonzero(abs_z >= z_critical(threshold)))
        print(f"  {threshold_key(threshold)}: {n_hits:,} hits")
