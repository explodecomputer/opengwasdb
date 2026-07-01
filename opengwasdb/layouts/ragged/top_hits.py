"""Ragged top-hit index builder — mirrors opengwasdb/layouts/dense/top_hits.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.layouts.dense.constants import TOP_HIT_THRESHOLDS
from opengwasdb.layouts.dense.top_hits import threshold_key
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader


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

    # Derive analysis_index for every association via searchsorted on CSR offsets.
    # offsets[i+1] is the exclusive end of analysis i → searchsorted(offsets[1:], pos) gives i.
    positions = np.arange(len(vi_all), dtype=np.int64)
    analysis_indices = np.searchsorted(offsets[1:], positions, side="right").astype(np.int32)

    abs_z = np.abs(z_all)
    # Vectorised p-value: p = erfc(|z|/sqrt2).  Compute the z threshold once per
    # threshold value and do a numpy comparison instead of per-element Python calls.
    sqrt2 = math.sqrt(2.0)

    root = zarr.open_group(str(store_path / "data.zarr"), mode="a")
    top = root.require_group("top_hits")
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

    for threshold in thresholds:
        # Binary-search for z threshold equivalent to the p-value cutoff
        lo, hi, mid = 0.0, 40.0, 0.0
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if math.erfc(mid / sqrt2) > threshold:
                lo = mid
            else:
                hi = mid
        z_thresh = float(mid)

        keep = abs_z >= z_thresh
        if not keep.any():
            continue

        kept_vi = vi_all[keep]
        kept_ai = analysis_indices[keep]
        kept_abs = abs_z[keep]
        kept_z = z_all[keep]
        kept_se = se_all[keep]
        # Compute float64 p-values only for the survivors
        kept_p = np.array(
            [math.erfc(float(v) / sqrt2) for v in kept_abs.tolist()],
            dtype=np.float64,
        )

        # Sort by descending |z|, tie-break by analysis_index then variant_index
        order = np.lexsort((kept_ai, kept_vi, -kept_abs))
        kept_vi = kept_vi[order]
        kept_ai = kept_ai[order]
        kept_abs = kept_abs[order]
        kept_z = kept_z[order]
        kept_se = kept_se[order]
        kept_p = kept_p[order]

        key = threshold_key(threshold)
        if key in top:
            del top[key]
        group = top.create_group(key)
        chunk = max(1, min(len(kept_vi), 100_000))

        for name, data, dtype in [
            ("variant_index", kept_vi, "uint32"),
            ("analysis_index", kept_ai, "uint32"),
            ("abs_z", kept_abs, "float32"),
            ("z", kept_z, "float32"),
            ("se", kept_se, "float32"),
            ("p_value", kept_p, "float64"),
        ]:
            group.create_dataset(
                name,
                data=data.astype(dtype),
                chunks=(chunk,),
                compressor=compressor,
                dtype=dtype,
            )
        group.attrs["threshold"] = threshold
        print(f"  {key}: {len(kept_vi):,} hits")

    top.attrs["thresholds"] = list(thresholds)
