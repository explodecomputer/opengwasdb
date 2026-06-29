"""Dense top-hit index builder."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.layouts.dense.constants import TOP_HIT_THRESHOLDS
from opengwasdb.stats import p_value_from_z


def threshold_key(threshold: float) -> str:
    """Stable Zarr group key for a p-value threshold."""

    return f"p_{threshold:.0e}".replace("-", "_").replace("+", "")


def build_top_hit_indexes(
    store_path: str | Path,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> None:
    """Build ranked dense top-hit arrays for each configured p-value threshold."""

    root = zarr.open_group(str(Path(store_path) / "data.zarr"), mode="a")
    z = root["z"][:].astype("float32")
    top = root.require_group("top_hits")
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

    finite = np.isfinite(z)
    rows, cols = np.where(finite)
    values = z[rows, cols]
    p_values = np.array([p_value_from_z(float(value)) for value in values], dtype="float64")
    abs_z = np.abs(values).astype("float32")

    for threshold in thresholds:
        key = threshold_key(threshold)
        if key in top:
            del top[key]
        group = top.create_group(key)
        keep = p_values <= threshold
        kept_rows = rows[keep].astype("uint32")
        kept_cols = cols[keep].astype("uint32")
        kept_abs_z = abs_z[keep].astype("float32")
        kept_z = values[keep].astype("float32")
        kept_p = p_values[keep].astype("float64")
        order = np.lexsort((kept_cols, kept_rows, -kept_abs_z))
        kept_rows = kept_rows[order]
        kept_cols = kept_cols[order]
        kept_abs_z = kept_abs_z[order]
        kept_z = kept_z[order]
        kept_p = kept_p[order]
        chunk = max(1, min(len(kept_rows), 100_000))
        group.create_dataset(
            "variant_index",
            data=kept_rows,
            chunks=(chunk,),
            compressor=compressor,
            dtype="uint32",
        )
        group.create_dataset(
            "analysis_index",
            data=kept_cols,
            chunks=(chunk,),
            compressor=compressor,
            dtype="uint32",
        )
        group.create_dataset(
            "abs_z",
            data=kept_abs_z,
            chunks=(chunk,),
            compressor=compressor,
            dtype="float32",
        )
        group.create_dataset(
            "z",
            data=kept_z,
            chunks=(chunk,),
            compressor=compressor,
            dtype="float32",
        )
        group.create_dataset(
            "p_value",
            data=kept_p,
            chunks=(chunk,),
            compressor=compressor,
            dtype="float64",
        )
        group.attrs["threshold"] = threshold
    top.attrs["thresholds"] = list(thresholds)
