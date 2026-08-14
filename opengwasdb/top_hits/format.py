"""Top-Hit Index format constants and threshold math, shared by every layout."""

from __future__ import annotations

import math

from scipy.special import erfcinv  # type: ignore[import-untyped]

#: p-value tiers persisted as separate threshold groups (store-format spec
#: §7a, ADR 0032). Positionally paired with
#: ``opengwasdb.model.analyses.TOP_HIT_COUNT_COLUMNS``.
TOP_HIT_THRESHOLDS: tuple[float, ...] = (5e-8, 5e-6, 5e-4)

#: Zarr chunk length for every array in a threshold group.
TOP_HIT_CHUNK_SIZE = 16_384


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
