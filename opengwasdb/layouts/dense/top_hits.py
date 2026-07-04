"""Dense top-hit index builder."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc
from scipy.special import erfc, erfcinv

from opengwasdb.layouts.dense.constants import TOP_HIT_THRESHOLDS


def threshold_key(threshold: float) -> str:
    """Stable Zarr group key for a p-value threshold."""

    return f"p_{threshold:.0e}".replace("-", "_").replace("+", "")


def z_critical(threshold: float) -> float:
    """|z| cutoff equivalent to two-sided p <= threshold.

    The stored p-value is ``p = erfc(|z| / sqrt(2))`` (see the query facade and
    the old scalar ``p_value_from_z``), so ``p <= threshold`` is exactly
    ``|z| >= sqrt(2) * erfcinv(threshold)``. Thresholding on ``|z|`` lets the
    hot path avoid computing a p-value for every cell.
    """

    return math.sqrt(2.0) * float(erfcinv(threshold))


def write_top_hit_indexes(
    store_path: str | Path,
    rows: np.ndarray,
    cols: np.ndarray,
    z: np.ndarray,
    se: np.ndarray,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> None:
    """Write ranked top-hit groups from pre-collected candidate cells.

    ``rows``/``cols``/``z``/``se`` describe candidate cells — every cell with
    ``|z| >= z_critical(max(thresholds))`` (the loosest tier). Extra cells below
    every threshold are harmless; each tier re-filters by its own ``z_critical``.
    Because candidates are a tiny fraction of a dense matrix, only the
    significant cells are ever held in memory — there is no full-matrix scan
    here. Cells are ranked within each analysis by descending ``|z|``.
    """

    rows = np.asarray(rows, dtype="uint32")
    cols = np.asarray(cols, dtype="uint32")
    z = np.asarray(z, dtype="float32")
    se = np.asarray(se, dtype="float32")
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
        kept_p = erfc(kept_abs_z.astype("float64") / math.sqrt(2.0))
        order = np.lexsort((kept_cols, kept_rows, -kept_abs_z))
        kept_rows = kept_rows[order]
        kept_cols = kept_cols[order]
        kept_abs_z = kept_abs_z[order]
        kept_z = kept_z[order]
        kept_se = kept_se[order]
        kept_p = kept_p[order]
        chunk = max(1, min(len(kept_rows), 100_000))
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
        group.attrs["threshold"] = threshold
    top.attrs["thresholds"] = list(thresholds)


def build_top_hit_indexes(
    store_path: str | Path,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> None:
    """Build ranked top-hit arrays by scanning the full dense matrix.

    Used by build paths that do not harvest hits inline (the VCF builder
    harvests during its fill pass instead — see ``build_dense_from_vcf_manifest``
    and ``write_top_hit_indexes``). Collects only the candidate cells
    (``|z| >= z_critical(loosest)``) so it never materialises a p-value for
    every finite cell.
    """

    root = zarr.open_group(str(Path(store_path) / "data.zarr"), mode="r")
    z = root["z"][:]
    se = root["se"][:]
    loosest = z_critical(max(thresholds))
    mask = np.isfinite(z) & (np.abs(z) >= loosest)
    rows, cols = np.where(mask)
    write_top_hit_indexes(store_path, rows, cols, z[rows, cols], se[rows, cols], thresholds)
