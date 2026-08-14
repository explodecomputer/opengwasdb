"""Dense top-hit index builder: harvest candidate cells from the dense matrix.

The Top-Hit Index format itself -- write/read/counts/validation -- lives in
``opengwasdb.top_hits``. This module owns only what is genuinely Dense-specific:
scanning the dense ``z``/``se`` grid in row-bands to collect candidate cells.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import zarr

from opengwasdb.top_hits.format import TOP_HIT_THRESHOLDS, z_critical
from opengwasdb.top_hits.writer import write


def build_top_hit_indexes(
    store_path: str | Path,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> None:
    """(Re)build ranked top-hit arrays by scanning the stored dense matrix.

    Used by build paths that do not harvest hits inline, and to rebuild the index
    on an existing store. Scans ``z`` in row-bands (never the full matrix in RAM)
    and thresholds on the **stored** values, so the index matches exactly what a
    query reads back from ``z`` (issue 046). Collects only candidate cells
    (``|z| >= z_critical(loosest)``).
    """

    root = zarr.open_group(str(Path(store_path) / "data.zarr"), mode="r")
    z_arr = root["z"]
    se_arr = root["se"]
    imputed_arr = root["imputed"] if "imputed" in root else None
    n_variants, n_analyses = int(z_arr.shape[0]), int(z_arr.shape[1])
    loosest = z_critical(max(thresholds))
    band_rows = max(int(z_arr.chunks[0]), 250_000)

    rows_parts: list[np.ndarray] = []
    cols_parts: list[np.ndarray] = []
    z_parts: list[np.ndarray] = []
    se_parts: list[np.ndarray] = []
    imputed_parts: list[np.ndarray] = []
    for r0 in range(0, n_variants, band_rows):
        r1 = min(r0 + band_rows, n_variants)
        z_band = z_arr[r0:r1]
        mask = np.abs(z_band.astype("float32")) >= loosest  # NaN compares False
        br, bc = np.where(mask)
        if len(br):
            se_band = se_arr[r0:r1]
            rows_parts.append(br.astype(np.int64) + r0)
            cols_parts.append(bc.astype(np.int64))
            z_parts.append(z_band[br, bc].astype("float32"))
            se_parts.append(se_band[br, bc].astype("float32"))
            if imputed_arr is not None:
                imputed_band = imputed_arr[r0:r1]
                imputed_parts.append(imputed_band[br, bc].astype("uint8"))

    if rows_parts:
        rows = np.concatenate(rows_parts)
        cols = np.concatenate(cols_parts)
        z = np.concatenate(z_parts)
        se = np.concatenate(se_parts)
        imputed = np.concatenate(imputed_parts) if imputed_parts else None
    else:
        rows = np.empty(0, dtype=np.int64)
        cols = np.empty(0, dtype=np.int64)
        z = np.empty(0, dtype=np.float32)
        se = np.empty(0, dtype=np.float32)
        imputed = np.empty(0, dtype=np.uint8) if imputed_arr is not None else None
    write(store_path, rows, cols, z, se, n_analyses, thresholds, imputed=imputed)
