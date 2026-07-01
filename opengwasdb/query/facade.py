"""Layout-independent query facade."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.index import analysis_by_id, connect
from opengwasdb.layouts.dense.top_hits import threshold_key
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.model.enums import PrimaryStorageLayout
from opengwasdb.stats import p_value_from_z
from opengwasdb.store.open import OpenGWASDBStore, open_store
from opengwasdb.traits.axis import TraitsAxisReader
from opengwasdb.variants import VariantAxis


def _empty_result() -> dict[str, np.ndarray]:
    return {
        "variant_index": np.empty(0, dtype="int32"),
        "analysis_index": np.empty(0, dtype="int32"),
        "z": np.empty(0, dtype="float32"),
        "se": np.empty(0, dtype="float32"),
    }


class StoreQuery:
    """Public query object that hides the physical store layout — Dense stores."""

    def __init__(self, store: OpenGWASDBStore):
        self.store = store
        self._connection = connect(store.index_path)
        self._root = zarr.open_group(str(store.data_path), mode="r")
        self._variant_axis = VariantAxis(store.path, self._connection)

    def close(self) -> None:
        self._variant_axis.close()
        self._connection.close()

    def variants_table(self) -> dict[int, dict]:
        """Return all variants keyed by variant_index."""
        return {
            r.variant_index: {
                "alid": r.alid,
                "chromosome": r.chromosome,
                "position": r.position,
                "effect_allele": r.effect_allele,
                "other_allele": r.other_allele,
                "rsid": r.rsid,
            }
            for r in self._variant_axis.all()
        }

    def analyses_table(self) -> dict[int, dict]:
        """Return all analyses keyed by analysis_index."""
        rows = self._connection.execute(
            "SELECT * FROM analyses ORDER BY analysis_index"
        ).fetchall()
        return {
            int(row["analysis_index"]): {
                "analysis_id": str(row["analysis_id"]),
                "phenotype_id": row["phenotype_id"],
                "phenotype_label": row["phenotype_label"],
                "analysis_label": row["analysis_label"],
                "stored_effect_scale": str(row["stored_effect_scale"]),
            }
            for row in rows
        }

    def analysis(self, analysis_id: str) -> dict[str, np.ndarray]:
        """Return all finite associations for one analysis."""
        analysis = analysis_by_id(self._connection, analysis_id)
        if analysis is None:
            return _empty_result()
        col = int(analysis["analysis_index"])
        z_col = self._root["z"][:, col].astype("float32")
        se_col = self._root["se"][:, col].astype("float32")
        mask = np.isfinite(z_col) & np.isfinite(se_col)
        rows = np.where(mask)[0].astype("int32")
        return {
            "variant_index": rows,
            "analysis_index": np.full(len(rows), col, dtype="int32"),
            "z": z_col[mask],
            "se": se_col[mask],
        }

    def phewas(self, identifier: str) -> dict[str, np.ndarray]:
        """Return one variant across all analyses."""
        variant = self._variant_axis.by_identifier(identifier)
        if variant is None:
            return _empty_result()
        row = variant.variant_index
        z_row = self._root["z"][row, :].astype("float32")
        se_row = self._root["se"][row, :].astype("float32")
        mask = np.isfinite(z_row) & np.isfinite(se_row)
        cols = np.where(mask)[0].astype("int32")
        return {
            "variant_index": np.full(len(cols), row, dtype="int32"),
            "analysis_index": cols,
            "z": z_row[mask],
            "se": se_row[mask],
        }

    def range(self, chromosome: str, start: int, end: int) -> dict[str, np.ndarray]:
        """Return finite associations in a genomic range."""
        row_indices = self._variant_axis.range_indices(chromosome, start, end)
        if len(row_indices) == 0:
            return _empty_result()
        z_block = self._root["z"].oindex[row_indices, :].astype("float32")
        se_block = self._root["se"].oindex[row_indices, :].astype("float32")
        mask = np.isfinite(z_block) & np.isfinite(se_block)
        rows_rel, cols = np.where(mask)
        return {
            "variant_index": row_indices[rows_rel],
            "analysis_index": cols.astype("int32"),
            "z": z_block[mask],
            "se": se_block[mask],
        }

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
    ) -> dict[str, np.ndarray]:
        """Return finite associations for a specific variant × analysis set."""
        variants = [
            v
            for id_ in identifiers
            if (v := self._variant_axis.by_identifier(id_)) is not None
        ]
        analyses = [
            a
            for aid in analysis_ids
            if (a := analysis_by_id(self._connection, aid)) is not None
        ]
        if not variants or not analyses:
            return _empty_result()
        row_indices = [v.variant_index for v in variants]
        col_indices = [int(a["analysis_index"]) for a in analyses]
        z_block = self._root["z"].oindex[row_indices, :][:, col_indices].astype("float32")
        se_block = self._root["se"].oindex[row_indices, :][:, col_indices].astype("float32")
        mask = np.isfinite(z_block) & np.isfinite(se_block)
        rows_rel, cols_rel = np.where(mask)
        return {
            "variant_index": np.array([row_indices[r] for r in rows_rel], dtype="int32"),
            "analysis_index": np.array([col_indices[c] for c in cols_rel], dtype="int32"),
            "z": z_block[mask],
            "se": se_block[mask],
        }

    def top_hits(
        self,
        *,
        threshold: float = 5e-8,
        limit: int | None = None,
    ) -> dict[str, np.ndarray]:
        """Return ranked top-hit associations using the dense top-hit index."""
        key = threshold_key(threshold)
        path = f"top_hits/{key}"
        if path not in self._root:
            return _empty_result()
        group = self._root[path]
        variant_indices = group["variant_index"][:].astype("int32")
        analysis_indices = group["analysis_index"][:].astype("int32")
        z_values = group["z"][:].astype("float32")
        if "se" in group:
            se_values = group["se"][:].astype("float32")
        else:
            unique_rows = np.unique(variant_indices).astype("int64")
            row_offsets = {int(r): i for i, r in enumerate(unique_rows)}
            se_block = self._root["se"].oindex[unique_rows, :]
            se_values = np.array(
                [
                    float(se_block[row_offsets[int(r)], int(c)])
                    for r, c in zip(variant_indices, analysis_indices, strict=True)
                ],
                dtype="float32",
            )
        if limit is not None:
            variant_indices = variant_indices[:limit]
            analysis_indices = analysis_indices[:limit]
            z_values = z_values[:limit]
            se_values = se_values[:limit]
        return {
            "variant_index": variant_indices,
            "analysis_index": analysis_indices,
            "z": z_values,
            "se": se_values,
        }


class RaggedStoreQuery:
    """Public query object that hides the physical store layout — Ragged stores."""

    def __init__(self, store: OpenGWASDBStore):
        self.store = store
        self._csr = RaggedCSRReader(store.path)
        self._variant_axis = VariantAxis(store.path)
        self._traits_reader = TraitsAxisReader(store.path)
        self._db: sqlite3.Connection = sqlite3.connect(
            str(store.path / "index.sqlite")
        )
        self._db.row_factory = sqlite3.Row

    def close(self) -> None:
        self._variant_axis.close()
        self._traits_reader.close()
        self._db.close()

    def __enter__(self) -> RaggedStoreQuery:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _resolve_analysis_id(self, analysis_id: str) -> int | None:
        row = self._db.execute(
            "SELECT analysis_index FROM analyses WHERE probe_id = ? OR analysis_id = ? LIMIT 1",
            (analysis_id, analysis_id),
        ).fetchone()
        return None if row is None else int(row["analysis_index"])

    def analysis(self, analysis_id: str) -> dict[str, np.ndarray]:
        """All associations for one analysis (probe_id or analysis_id lookup)."""
        idx = self._resolve_analysis_id(analysis_id)
        if idx is None:
            return _empty_result()
        assoc = self._csr.get_analysis(idx)
        return {
            "variant_index": assoc.variant_index.astype("int32"),
            "analysis_index": np.full(len(assoc.z), idx, dtype="int32"),
            "z": assoc.z.astype("float32"),
            "se": assoc.se.astype("float32"),
        }

    def range(self, chromosome: str, start: int, end: int) -> dict[str, np.ndarray]:
        """All associations where the variant falls in [start, end]."""
        variant_set = set(
            self._variant_axis.range_indices(chromosome, start, end).tolist()
        )
        if not variant_set:
            return _empty_result()

        offsets = self._csr._offsets[:]
        vi_all = self._csr._variant_index[:]
        z_all = self._csr._z[:]
        se_all = self._csr._se[:]

        mask = np.isin(vi_all, np.array(sorted(variant_set), dtype=np.int32))
        hit_positions = np.where(mask)[0]
        if len(hit_positions) == 0:
            return _empty_result()

        # Derive analysis_index for each hit via searchsorted on offsets
        analysis_indices = np.searchsorted(offsets[1:], hit_positions, side="right").astype("int32")

        return {
            "variant_index": vi_all[hit_positions].astype("int32"),
            "analysis_index": analysis_indices,
            "z": z_all[hit_positions].astype("float32"),
            "se": se_all[hit_positions].astype("float32"),
        }

    def range_by_probe(self, chromosome: str, start: int, end: int) -> dict[str, np.ndarray]:
        """All associations for analyses whose probe/TSS falls in [start, end]."""
        trait_records = self._traits_reader.range(chromosome, start, end)
        if not trait_records:
            return _empty_result()

        all_vi: list[np.ndarray] = []
        all_ai: list[np.ndarray] = []
        all_z: list[np.ndarray] = []
        all_se: list[np.ndarray] = []

        for rec in trait_records:
            assoc = self._csr.get_analysis(rec.analysis_index)
            if len(assoc.z) == 0:
                continue
            all_vi.append(assoc.variant_index.astype("int32"))
            all_ai.append(np.full(len(assoc.z), rec.analysis_index, dtype="int32"))
            all_z.append(assoc.z.astype("float32"))
            all_se.append(assoc.se.astype("float32"))

        if not all_vi:
            return _empty_result()
        return {
            "variant_index": np.concatenate(all_vi),
            "analysis_index": np.concatenate(all_ai),
            "z": np.concatenate(all_z),
            "se": np.concatenate(all_se),
        }

    def phewas(self, identifier: str) -> dict[str, np.ndarray]:
        """All analyses that have an association for a given variant identifier.

        O(n_total_associations) scan — acceptable for exploratory use; add a
        variant-centric CSR index (issue deferred) for production phewas.
        """
        variant = self._variant_axis.by_identifier(identifier)
        if variant is None:
            return _empty_result()
        target_vi = np.int32(variant.variant_index)

        offsets = self._csr._offsets[:]
        vi_all = self._csr._variant_index[:]
        z_all = self._csr._z[:]
        se_all = self._csr._se[:]

        hit_positions = np.where(vi_all == target_vi)[0]
        if len(hit_positions) == 0:
            return _empty_result()

        analysis_indices = np.searchsorted(offsets[1:], hit_positions, side="right").astype("int32")
        return {
            "variant_index": np.full(len(hit_positions), target_vi, dtype="int32"),
            "analysis_index": analysis_indices,
            "z": z_all[hit_positions].astype("float32"),
            "se": se_all[hit_positions].astype("float32"),
        }

    def top_hits(
        self,
        *,
        threshold: float = 5e-8,
        limit: int | None = None,
    ) -> dict[str, np.ndarray]:
        """Associations passing a significance threshold (full scan)."""
        offsets = self._csr._offsets[:]
        vi_all = self._csr._variant_index[:]
        z_all = self._csr._z[:]
        se_all = self._csr._se[:]

        import math
        # p = erfc(|z|/sqrt2) ≤ threshold  →  |z| ≥ z_thresh
        # Compute z_thresh once via binary search on erfc (avoids per-element Python loop).
        sqrt2 = math.sqrt(2.0)
        lo, hi = 0.0, 40.0
        mid = 0.0
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if math.erfc(mid / sqrt2) > threshold:
                lo = mid
            else:
                hi = mid
        z_thresh = float(mid)
        z_f32 = z_all.astype("float32")
        mask = np.abs(z_f32) >= z_thresh
        hit_positions = np.where(mask)[0]

        if len(hit_positions) == 0:
            return _empty_result()

        # Sort by descending |z|
        order = np.argsort(-np.abs(z_f32[hit_positions]))
        hit_positions = hit_positions[order]
        if limit is not None:
            hit_positions = hit_positions[:limit]

        analysis_indices = np.searchsorted(offsets[1:], hit_positions, side="right").astype("int32")
        return {
            "variant_index": vi_all[hit_positions].astype("int32"),
            "analysis_index": analysis_indices,
            "z": z_f32[hit_positions],
            "se": se_all[hit_positions].astype("float32"),
        }

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
    ) -> dict[str, np.ndarray]:
        """Associations for a specific variant × analysis set."""
        variants = [
            v
            for id_ in identifiers
            if (v := self._variant_axis.by_identifier(id_)) is not None
        ]
        if not variants:
            return _empty_result()

        target_vi = {v.variant_index for v in variants}
        all_vi, all_ai, all_z, all_se = [], [], [], []

        for aid in analysis_ids:
            idx = self._resolve_analysis_id(aid)
            if idx is None:
                continue
            assoc = self._csr.get_analysis(idx)
            if len(assoc.variant_index) == 0:
                continue
            sub_mask = np.isin(assoc.variant_index, np.array(sorted(target_vi), dtype=np.int32))
            if not sub_mask.any():
                continue
            all_vi.append(assoc.variant_index[sub_mask].astype("int32"))
            all_ai.append(np.full(sub_mask.sum(), idx, dtype="int32"))
            all_z.append(assoc.z[sub_mask].astype("float32"))
            all_se.append(assoc.se[sub_mask].astype("float32"))

        if not all_vi:
            return _empty_result()
        return {
            "variant_index": np.concatenate(all_vi),
            "analysis_index": np.concatenate(all_ai),
            "z": np.concatenate(all_z),
            "se": np.concatenate(all_se),
        }


def query_store(path: str | Path) -> StoreQuery | RaggedStoreQuery:
    """Open a store and return the layout-independent query facade."""
    store = open_store(path)
    if store.manifest.primary_layout is PrimaryStorageLayout.RAGGED:
        return RaggedStoreQuery(store)
    return StoreQuery(store)
