"""Read the Top-Hit Index back: one threshold group's slice, or per-Analysis counts."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.model.analyses import TOP_HIT_COUNT_COLUMNS
from opengwasdb.top_hits.format import TOP_HIT_THRESHOLDS, threshold_key

#: Positional pairing of TOP_HIT_THRESHOLDS with the analyses.tsv/SQLite
#: column each tier persists to (opengwasdb.model.analyses.TOP_HIT_COUNT_COLUMNS).
#: A dict, not a zip over the caller's own `thresholds`, so counts() accepts
#: any subset of TOP_HIT_THRESHOLDS rather than silently requiring exactly three.
_THRESHOLD_COLUMNS = dict(zip(TOP_HIT_THRESHOLDS, TOP_HIT_COUNT_COLUMNS, strict=True))

#: Fallback source for a column value at (variant_index, analysis_index)
#: pairs, used when a threshold group predates that column being written
#: directly into the index.
CellFallback = Callable[[np.ndarray, np.ndarray], np.ndarray]


class TopHitReader:
    """Read one threshold tier's group without exposing its physical arrays.

    ``slice()`` is the one entry point callers need: it resolves the bounds
    for an (optional) Analysis, reads every column the query facade's result
    shape requires, and applies ``se_fallback``/``imputed_fallback`` when a
    group predates having its own ``se``/``imputed`` column -- the same
    fallback every layout's ``top_hits()`` otherwise re-implemented at each
    call site.
    """

    def __init__(self, group: zarr.Group):
        self.group = group

    def bounds(self, analysis_index: int | None) -> tuple[int, int]:
        if analysis_index is None:
            return 0, int(self.group["z"].shape[0])
        offsets = self.group["analysis_offsets"]
        if analysis_index < 0 or analysis_index + 1 >= int(offsets.shape[0]):
            return 0, 0
        pair = offsets[analysis_index : analysis_index + 2]
        return int(pair[0]), int(pair[1])

    def slice(
        self,
        analysis_index: int | None,
        *,
        se_fallback: CellFallback | None = None,
        imputed_fallback: CellFallback | None = None,
    ) -> dict[str, np.ndarray]:
        """Return ``{variant_index, analysis_index, z, se, imputed}`` for one
        Analysis (or the whole group when ``analysis_index`` is ``None``).

        ``imputed`` defaults to all-zeros when the group has no ``imputed``
        column and no ``imputed_fallback`` is given -- an Observed-Only store
        has nothing to impute, so that is the correct answer, not a missing one.
        """
        bounds = self.bounds(analysis_index)
        variant_index = self._read("variant_index", bounds, "int32")
        analysis_index_arr = self._read("analysis_index", bounds, "int32")
        z = self._read("z", bounds, "float32")

        if "se" in self.group:
            se = self._read("se", bounds, "float32")
        elif se_fallback is not None:
            se = se_fallback(variant_index, analysis_index_arr)
        else:
            raise KeyError("top-hit group has no 'se' column and no se_fallback was given")

        if "imputed" in self.group:
            imputed = self._read("imputed", bounds, "uint8")
        elif imputed_fallback is not None:
            imputed = imputed_fallback(variant_index, analysis_index_arr)
        else:
            imputed = np.zeros(len(variant_index), dtype=np.uint8)

        return {
            "variant_index": variant_index,
            "analysis_index": analysis_index_arr,
            "z": z,
            "se": se,
            "imputed": imputed,
        }

    def _read(self, name: str, bounds: tuple[int, int], dtype: str) -> np.ndarray:
        start, stop = bounds
        return np.asarray(self.group[name][start:stop], dtype=dtype)


def counts(
    store_path: str | Path,
    n_analyses: int,
    thresholds: tuple[float, ...] = TOP_HIT_THRESHOLDS,
) -> dict[str, list[int]]:
    """Per-Analysis hit counts for each threshold tier, from an already-built
    top-hit index (``writer.write()`` -- every layout's builder goes through
    it). Keyed by the persisted column each tier corresponds to (ADR 0032),
    in ``analysis_index`` order.

    ``n_analyses`` is needed for the fallback path below, not just for
    sizing a zero-filled result: real pre-issue-#22-era stores (e.g. ukb-b)
    have a top-hit index that predates the ``analysis_offsets`` array this
    function otherwise reads directly, so it falls back to counting the
    flat ``analysis_index`` array by hand for those.
    """
    root = zarr.open_group(str(Path(store_path) / "data.zarr"), mode="r")
    top = root["top_hits"]
    result: dict[str, list[int]] = {}
    for threshold in thresholds:
        column = _THRESHOLD_COLUMNS[threshold]
        group = top[threshold_key(threshold)]
        if "analysis_offsets" in group:
            offsets = np.asarray(group["analysis_offsets"], dtype=np.int64)
            result[column] = (offsets[1:] - offsets[:-1]).tolist()
        else:
            analysis_index = np.asarray(group["analysis_index"], dtype=np.int64)
            result[column] = np.bincount(analysis_index, minlength=n_analyses)[:n_analyses].tolist()
    return result
