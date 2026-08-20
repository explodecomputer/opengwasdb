"""Zarr-backed Compressed Sparse Row storage for Ragged Layout associations."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np
import zarr
from numcodecs import Blosc

RAGGED_ZARR_PATH = "data.zarr/ragged"
_COMPRESSOR = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
# Chunk size for the flat association arrays (~400 KB per chunk at float16).
_ASSOC_CHUNK = 200_000
_OFFSET_CHUNK = 10_000


class AnalysisAssociations(NamedTuple):
    variant_index: np.ndarray  # int32
    z: np.ndarray              # float16
    se: np.ndarray             # float16
    eaf: np.ndarray            # float32, all-NaN when the store carries no EAF


class RaggedCSRWriter:
    """Accumulate per-analysis associations and flush to zarr CSR arrays."""

    def __init__(self) -> None:
        self._variant_indices: list[np.ndarray] = []
        self._zscores: list[np.ndarray] = []
        self._ses: list[np.ndarray] = []
        self._eafs: list[np.ndarray] = []
        self._any_eaf = False
        self._offsets: list[int] = [0]

    def add_analysis(
        self,
        variant_index: np.ndarray,
        z: np.ndarray,
        se: np.ndarray,
        eaf: np.ndarray | None = None,
    ) -> None:
        """Append one analysis. Arrays must be parallel and the same length.

        `eaf` is optional (ADR 0036): an Analysis whose source reports no
        frequency passes None and contributes all-NaN, so the flat array stays
        aligned with `z`/`se` whatever mix of Analyses a build spans. If *no*
        Analysis supplies one, `flush` writes no `eaf` array at all and the
        store looks exactly as it did before EAF existed.
        """
        n = len(variant_index)
        self._variant_indices.append(np.asarray(variant_index, dtype=np.int32))
        self._zscores.append(np.asarray(z, dtype=np.float16))
        self._ses.append(np.asarray(se, dtype=np.float16))
        if eaf is None:
            self._eafs.append(np.full(n, np.nan, dtype=np.float32))
        else:
            self._eafs.append(np.asarray(eaf, dtype=np.float32))
            self._any_eaf = True
        self._offsets.append(self._offsets[-1] + n)

    @property
    def n_analyses(self) -> int:
        return len(self._offsets) - 1

    @property
    def n_associations(self) -> int:
        return self._offsets[-1]

    def flush(self, store_path: str | Path) -> None:
        """Write CSR arrays to data.zarr/ragged/ inside store_path."""
        out = Path(store_path) / RAGGED_ZARR_PATH
        root = zarr.open_group(str(out), mode="w")

        offsets_arr = np.asarray(self._offsets, dtype=np.int64)

        if self.n_associations > 0:
            vi_arr = np.concatenate(self._variant_indices).astype(np.int32)
            z_arr = np.concatenate(self._zscores).astype(np.float16)
            se_arr = np.concatenate(self._ses).astype(np.float16)
            eaf_arr = np.concatenate(self._eafs).astype(np.float32)
        else:
            vi_arr = np.empty(0, dtype=np.int32)
            z_arr = np.empty(0, dtype=np.float16)
            se_arr = np.empty(0, dtype=np.float16)
            eaf_arr = np.empty(0, dtype=np.float32)

        root.create_dataset(
            "offsets", data=offsets_arr,
            chunks=(_OFFSET_CHUNK,), compressor=_COMPRESSOR, dtype=np.int64,
        )
        root.create_dataset(
            "variant_index", data=vi_arr,
            chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.int32,
        )
        root.create_dataset(
            "z", data=z_arr,
            chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.float16,
        )
        root.create_dataset(
            "se", data=se_arr,
            chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.float16,
        )
        if self._any_eaf:
            root.create_dataset(
                "eaf", data=eaf_arr,
                chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.float32,
            )
        root.attrs["layout"] = "ragged"
        root.attrs["completion_state"] = "observed_only"
        root.attrs["n_analyses"] = self.n_analyses
        root.attrs["n_associations"] = self.n_associations


class RaggedCSRReader:
    """Read per-analysis associations from zarr CSR arrays."""

    def __init__(self, store_path: str | Path):
        path = Path(store_path) / RAGGED_ZARR_PATH
        self._root = zarr.open_group(str(path), mode="r")
        self._offsets: zarr.Array = self._root["offsets"]
        self._variant_index: zarr.Array = self._root["variant_index"]
        self._z: zarr.Array = self._root["z"]
        self._se: zarr.Array = self._root["se"]
        # Absent on stores built before ADR 0036, and on stores whose sources
        # report no frequency at all -- both read back as all-NaN.
        self._eaf: zarr.Array | None = self._root["eaf"] if "eaf" in self._root else None

    @property
    def n_analyses(self) -> int:
        return int(self._root.attrs.get("n_analyses", len(self._offsets) - 1))

    @property
    def n_associations(self) -> int:
        return int(self._root.attrs.get("n_associations", len(self._variant_index)))

    def get_analysis(self, analysis_index: int) -> AnalysisAssociations:
        """Return (variant_index, z, se, eaf) arrays for one analysis. O(1) zarr reads."""
        offsets = self._offsets[analysis_index: analysis_index + 2]
        start, end = int(offsets[0]), int(offsets[1])
        if start == end:
            return AnalysisAssociations(
                variant_index=np.empty(0, dtype=np.int32),
                z=np.empty(0, dtype=np.float16),
                se=np.empty(0, dtype=np.float16),
                eaf=np.empty(0, dtype=np.float32),
            )
        return AnalysisAssociations(
            variant_index=self._variant_index[start:end],
            z=self._z[start:end],
            se=self._se[start:end],
            eaf=self.eaf_slice(start, end),
        )

    def eaf_slice(self, start: int, end: int) -> np.ndarray:
        """`eaf[start:end]`, or all-NaN when this store carries no EAF array."""
        if self._eaf is None:
            return np.full(end - start, np.nan, dtype=np.float32)
        return np.asarray(self._eaf[start:end], dtype=np.float32)

    def get_analyses(self, analysis_indices: list[int]) -> list[AnalysisAssociations]:
        """Return associations for multiple analyses."""
        return [self.get_analysis(i) for i in analysis_indices]
