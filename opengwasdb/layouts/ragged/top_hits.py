"""Ragged top-hit index builder: harvest candidate cells from the CSR store.

The Top-Hit Index format itself -- write/read/counts/validation -- lives in
``opengwasdb.top_hits``. This module owns only what is genuinely Ragged-specific:
deriving each association's Analysis from the CSR offsets, and (unlike Dense
and Hybrid, whose per-Analysis metadata lives in ``analyses.tsv``) persisting
Top-Hit Counts onto Ragged's own SQLite ``analyses`` table (ADR 0030).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.model.analyses import TOP_HIT_COUNT_COLUMNS
from opengwasdb.top_hits.format import TOP_HIT_THRESHOLDS, threshold_key, z_critical
from opengwasdb.top_hits.reader import counts
from opengwasdb.top_hits.writer import write

if TYPE_CHECKING:
    from opengwasdb.store.open import OpenGWASDBStore, StagedRelease


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


def apply_hit_counts(store: OpenGWASDBStore | StagedRelease, n_analyses: int) -> None:
    """Write per-Analysis Top-Hit Counts from ``store``'s already-built
    top-hit index onto its SQLite ``analyses`` table (ADR 0032).

    The Ragged counterpart of ``opengwasdb.layouts.dense.build.add_hit_counts``:
    Ragged Analytical Metadata lives in SQLite, not ``analyses.tsv``
    (ADR 0030), so counts land in a column ``UPDATE`` rather than a rebuilt
    ``AnalysisMetadata`` list. Call after ``build_ragged_top_hit_indexes()``.
    """
    hit_counts = counts(store.path, n_analyses)
    conn = store.index_connection()
    try:
        set_clause = ", ".join(f"{column} = ?" for column in TOP_HIT_COUNT_COLUMNS)
        for i in range(n_analyses):
            conn.execute(
                f"UPDATE analyses SET {set_clause} WHERE analysis_index = ?",
                (*(hit_counts[column][i] for column in TOP_HIT_COUNT_COLUMNS), i),
            )
        conn.commit()
    finally:
        conn.close()
