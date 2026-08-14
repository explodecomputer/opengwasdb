"""Write the Top-Hit Index from pre-collected candidate cells."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc
from scipy.special import erfc  # type: ignore[import-untyped]

from opengwasdb.top_hits.format import (
    TOP_HIT_CHUNK_SIZE,
    TOP_HIT_THRESHOLDS,
    threshold_key,
    z_critical,
)


def write(
    store_path: str | Path,
    rows: np.ndarray,
    cols: np.ndarray,
    z: np.ndarray,
    se: np.ndarray,
    n_analyses: int,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
    imputed: np.ndarray | None = None,
    chunk_size: int = TOP_HIT_CHUNK_SIZE,
) -> None:
    """Write ranked top-hit groups from pre-collected candidate cells.

    ``rows``/``cols``/``z``/``se`` describe candidate cells — every cell with
    ``|z| >= z_critical(max(thresholds))`` (the loosest tier). Extra cells below
    every threshold are harmless; each tier re-filters by its own ``z_critical``.
    A caller with a naturally sparse source (e.g. a Ragged CSR) may pass every
    association without pre-filtering -- the per-threshold filter below still
    applies, just over a larger input. Cells are ordered by analysis index and
    canonical genomic position. When present, ``imputed`` is written in the
    same order so completed-store queries can label top hits without random
    reads back into the association arrays.

    ``n_analyses`` sizes ``analysis_offsets`` -- the caller supplies it rather
    than this function inferring it from physical storage shape, since that
    shape differs by layout (a Dense matrix's column count; a Ragged CSR's
    offsets length).
    """

    rows = np.asarray(rows, dtype="uint32")
    cols = np.asarray(cols, dtype="uint32")
    z = np.asarray(z, dtype="float32")
    se = np.asarray(se, dtype="float32")
    imputed_values = None if imputed is None else np.asarray(imputed, dtype="uint8")
    abs_z = np.abs(z).astype("float32")

    root = zarr.open_group(str(Path(store_path) / "data.zarr"), mode="a")
    top = root.require_group("top_hits")
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

    for threshold in thresholds:
        key = threshold_key(threshold)
        if key in top:
            del top[key]
        group = top.create_group(key)
        keep = abs_z >= z_critical(threshold)
        kept_rows = rows[keep]
        kept_cols = cols[keep]
        kept_abs_z = abs_z[keep]
        kept_z = z[keep]
        kept_se = se[keep]
        kept_imputed = None if imputed_values is None else imputed_values[keep]
        kept_p = erfc(kept_abs_z.astype("float64") / math.sqrt(2.0))
        order = np.lexsort((kept_rows, kept_cols))
        kept_rows = kept_rows[order]
        kept_cols = kept_cols[order]
        kept_abs_z = kept_abs_z[order]
        kept_z = kept_z[order]
        kept_se = kept_se[order]
        kept_imputed = None if kept_imputed is None else kept_imputed[order]
        kept_p = kept_p[order]
        offsets = np.empty(n_analyses + 1, dtype="uint64")
        offsets[0] = 0
        np.cumsum(
            np.bincount(kept_cols, minlength=n_analyses),
            dtype=np.uint64,
            out=offsets[1:],
        )
        chunk = max(1, min(len(kept_rows), chunk_size))
        group.create_dataset(
            "analysis_offsets", data=offsets, chunks=(len(offsets),),
            compressor=compressor, dtype="uint64"
        )
        group.create_dataset(
            "variant_index", data=kept_rows, chunks=(chunk,), compressor=compressor, dtype="uint32"
        )
        group.create_dataset(
            "analysis_index", data=kept_cols, chunks=(chunk,), compressor=compressor, dtype="uint32"
        )
        group.create_dataset(
            "abs_z", data=kept_abs_z, chunks=(chunk,), compressor=compressor, dtype="float32"
        )
        group.create_dataset(
            "z", data=kept_z, chunks=(chunk,), compressor=compressor, dtype="float32"
        )
        group.create_dataset(
            "se", data=kept_se, chunks=(chunk,), compressor=compressor, dtype="float32"
        )
        group.create_dataset(
            "p_value", data=kept_p, chunks=(chunk,), compressor=compressor, dtype="float64"
        )
        if kept_imputed is not None:
            group.create_dataset(
                "imputed",
                data=kept_imputed,
                chunks=(chunk,),
                compressor=compressor,
                dtype="uint8",
            )
        group.attrs["threshold"] = threshold
        group.attrs["order"] = "analysis_index,variant_index"
    top.attrs["thresholds"] = list(thresholds)
