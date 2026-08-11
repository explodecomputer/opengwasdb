"""Layout-independent query facade."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.index import AnalysesIndex, connect
from opengwasdb.layouts.dense.top_hits import DenseTopHitReader, threshold_key
from opengwasdb.layouts.hybrid.layout import dense_component_path, dense_to_shared_path
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.model.enums import CompletionState, PrimaryStorageLayout
from opengwasdb.store.open import OpenGWASDBStore, open_store
from opengwasdb.traits.axis import TraitsAxisReader
from opengwasdb.variants import VariantAxis


def _empty_result() -> dict[str, np.ndarray]:
    return {
        "variant_index": np.empty(0, dtype="int32"),
        "analysis_index": np.empty(0, dtype="int32"),
        "z": np.empty(0, dtype="float32"),
        "se": np.empty(0, dtype="float32"),
        "association_status": np.empty(0, dtype=object),
    }


def _status_array(imputed_flags: np.ndarray, z_vals: np.ndarray) -> np.ndarray:
    """Derive association_status strings from imputed mask and z values."""
    out = np.where(imputed_flags == 1, "imputed", "observed").astype(object)
    out[~np.isfinite(z_vals)] = "missing"
    return out


class StoreQuery:
    """Public query object that hides the physical store layout — Dense stores."""

    def __init__(self, store: OpenGWASDBStore):
        self.store = store
        self._connection = connect(store.index_path)
        self._analyses = AnalysesIndex(store.path)
        self._root = zarr.open_group(str(store.data_path), mode="r")
        self._variant_axis = VariantAxis(store.path, self._connection)
        self._is_completed = (
            store.manifest.completion_state is CompletionState.REFERENCE_COMPLETED
        )
        self._imputed: zarr.Array | None = (
            self._root["imputed"] if self._is_completed and "imputed" in self._root else None
        )

    def _imputed_pairs(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """Imputed flags for elementwise (row, col) pairs; all-zeros when not a completed store."""
        if self._imputed is None or len(rows) == 0:
            return np.zeros(len(rows), dtype=np.uint8)
        return self._imputed.vindex[rows, cols].astype(np.uint8)

    @staticmethod
    def _contiguous_row_slice(row_indices: np.ndarray) -> slice | None:
        if len(row_indices) == 0:
            return None
        start = int(row_indices[0])
        stop = int(row_indices[-1]) + 1
        if stop - start != len(row_indices):
            return None
        if not np.array_equal(row_indices, np.arange(start, stop, dtype=row_indices.dtype)):
            return None
        return slice(start, stop)

    @classmethod
    def _read_row_block(
        cls, array: zarr.Array, row_indices: np.ndarray, dtype: str | np.dtype
    ) -> np.ndarray:
        row_slice = cls._contiguous_row_slice(row_indices)
        if row_slice is not None:
            return array[row_slice, :].astype(dtype)
        return array.oindex[row_indices, :].astype(dtype)

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
        return {
            index: {
                "analysis_id": row["analysis_id"],
                "phenotype_id": row.get("phenotype_id") or None,
                "phenotype_label": row.get("phenotype_label") or None,
                "analysis_label": row.get("analysis_label") or None,
                "stored_effect_scale": row["stored_effect_scale"],
            }
            for index, row in self._analyses.all().items()
        }

    def analysis(self, analysis_id: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
        """Return all finite associations for one analysis."""
        analysis = self._analyses.by_id(analysis_id)
        if analysis is None:
            return _empty_result()
        col = int(analysis["analysis_index"])
        z_col = self._root["z"][:, col].astype("float32")
        se_col = self._root["se"][:, col].astype("float32")
        mask = np.isfinite(z_col) & np.isfinite(se_col)
        rows = np.where(mask)[0].astype("int32")
        cols = np.full(len(rows), col, dtype="int32")
        z_vals = z_col[mask]
        se_vals = se_col[mask]
        imp = self._imputed_pairs(rows, cols)
        if observed_only:
            keep = imp == 0
            rows, cols, z_vals, se_vals, imp = (
                rows[keep], cols[keep], z_vals[keep], se_vals[keep], imp[keep]
            )
        return {
            "variant_index": rows,
            "analysis_index": cols,
            "z": z_vals,
            "se": se_vals,
            "association_status": _status_array(imp, z_vals),
        }

    def phewas(self, identifier: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
        """Return one variant across all analyses."""
        variant = self._variant_axis.by_identifier(identifier)
        if variant is None:
            return _empty_result()
        row = variant.variant_index
        z_row = self._root["z"][row, :].astype("float32")
        se_row = self._root["se"][row, :].astype("float32")
        mask = np.isfinite(z_row) & np.isfinite(se_row)
        cols = np.where(mask)[0].astype("int32")
        rows = np.full(len(cols), row, dtype="int32")
        z_vals = z_row[mask]
        se_vals = se_row[mask]
        imp = self._imputed_pairs(rows, cols)
        if observed_only:
            keep = imp == 0
            rows, cols, z_vals, se_vals, imp = (
                rows[keep], cols[keep], z_vals[keep], se_vals[keep], imp[keep]
            )
        return {
            "variant_index": rows,
            "analysis_index": cols,
            "z": z_vals,
            "se": se_vals,
            "association_status": _status_array(imp, z_vals),
        }

    def range_phewas(
        self, chromosome: str, start: int, end: int, *, observed_only: bool = False
    ) -> dict[str, np.ndarray]:
        """Return finite associations for all variants in a genomic range (regional PheWAS)."""
        row_indices = self._variant_axis.range_indices(chromosome, start, end)
        if len(row_indices) == 0:
            return _empty_result()
        z_block = self._read_row_block(self._root["z"], row_indices, "float32")
        se_block = self._read_row_block(self._root["se"], row_indices, "float32")
        mask = np.isfinite(z_block) & np.isfinite(se_block)
        rows_rel, cols = np.where(mask)
        rows = row_indices[rows_rel].astype("int32")
        cols = cols.astype("int32")
        z_vals = z_block[mask]
        se_vals = se_block[mask]
        if self._imputed is None:
            imp = np.zeros(len(z_vals), dtype=np.uint8)
        else:
            imp_block = self._read_row_block(self._imputed, row_indices, np.uint8)
            imp = imp_block[mask]
        if observed_only:
            keep = imp == 0
            rows, cols, z_vals, se_vals, imp = (
                rows[keep], cols[keep], z_vals[keep], se_vals[keep], imp[keep]
            )
        return {
            "variant_index": rows,
            "analysis_index": cols,
            "z": z_vals,
            "se": se_vals,
            "association_status": _status_array(imp, z_vals),
        }

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
        *,
        observed_only: bool = False,
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
            if (a := self._analyses.by_id(aid)) is not None
        ]
        if not variants or not analyses:
            return _empty_result()
        row_indices = [v.variant_index for v in variants]
        col_indices = [int(a["analysis_index"]) for a in analyses]
        # Surgical orthogonal read: fetch only the chunks intersecting the
        # requested rows × cols, not the full analysis width per row. Under a
        # narrow analysis chunk this reads far fewer chunks (issue 052).
        z_block = self._root["z"].oindex[row_indices, col_indices].astype("float32")
        se_block = self._root["se"].oindex[row_indices, col_indices].astype("float32")
        mask = np.isfinite(z_block) & np.isfinite(se_block)
        rows_rel, cols_rel = np.where(mask)
        rows = np.array([row_indices[r] for r in rows_rel], dtype="int32")
        cols = np.array([col_indices[c] for c in cols_rel], dtype="int32")
        z_vals = z_block[mask]
        se_vals = se_block[mask]
        imp = self._imputed_pairs(rows, cols)
        if observed_only:
            keep = imp == 0
            rows, cols, z_vals, se_vals, imp = (
                rows[keep], cols[keep], z_vals[keep], se_vals[keep], imp[keep]
            )
        return {
            "variant_index": rows,
            "analysis_index": cols,
            "z": z_vals,
            "se": se_vals,
            "association_status": _status_array(imp, z_vals),
        }

    def top_hits(
        self,
        *,
        analysis_id: str | None = None,
        threshold: float = 5e-8,
        limit: int | None = None,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        """Return genomic-order top hits, optionally for one analysis."""
        key = threshold_key(threshold)
        path = f"top_hits/{key}"
        if path not in self._root:
            return _empty_result()
        group = self._root[path]
        analysis_index: int | None = None
        if analysis_id is not None:
            analysis = self._analyses.by_id(analysis_id)
            if analysis is None or "analysis_offsets" not in group:
                return _empty_result()
            analysis_index = int(analysis["analysis_index"])
        reader = DenseTopHitReader(group)
        bounds = reader.bounds(analysis_index)
        variant_indices = reader.read("variant_index", bounds, "int32")
        analysis_indices = reader.read("analysis_index", bounds, "int32")
        z_values = reader.read("z", bounds, "float32")
        if "se" in group:
            se_values = reader.read("se", bounds, "float32")
        else:
            # Pointwise (coordinate) read: fetch se at exactly the index cells,
            # not the full analysis width per row (issue 052).
            se_values = self._root["se"].vindex[
                variant_indices.astype("int64"), analysis_indices.astype("int64")
            ].astype("float32")
        if "imputed" in group:
            imp = reader.read("imputed", bounds, "uint8")
        else:
            imp = self._imputed_pairs(variant_indices, analysis_indices)
        if observed_only:
            keep = imp == 0
            variant_indices, analysis_indices, z_values, se_values, imp = (
                variant_indices[keep], analysis_indices[keep],
                z_values[keep], se_values[keep], imp[keep],
            )
        if limit is not None:
            variant_indices = variant_indices[:limit]
            analysis_indices = analysis_indices[:limit]
            z_values = z_values[:limit]
            se_values = se_values[:limit]
            imp = imp[:limit]
        return {
            "variant_index": variant_indices,
            "analysis_index": analysis_indices,
            "z": z_values,
            "se": se_values,
            "association_status": _status_array(imp, z_values),
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
        self._is_completed = (
            store.manifest.completion_state is CompletionState.REFERENCE_COMPLETED
        )
        # Load imputed mask when present (reference-completed stores).
        ragged_path = store.path / "data.zarr" / "ragged"
        self._imputed: zarr.Array | None = None
        if self._is_completed:
            try:
                _root = zarr.open_group(str(ragged_path), mode="r")
                if "imputed" in _root:
                    self._imputed = _root["imputed"]
            except Exception:  # noqa: BLE001
                pass

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
            "SELECT analysis_index FROM analyses WHERE trait_id = ? OR analysis_id = ? LIMIT 1",
            (analysis_id, analysis_id),
        ).fetchone()
        return None if row is None else int(row["analysis_index"])

    def _get_imputed_slice(self, start: int, end: int) -> np.ndarray:
        """Return imputed mask slice [start:end]; all-zeros if not a completed store."""
        if self._imputed is not None:
            return self._imputed[start:end].astype(np.uint8)
        return np.zeros(end - start, dtype=np.uint8)

    def analysis(self, analysis_id: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
        """All associations for one analysis (probe_id or analysis_id lookup)."""
        idx = self._resolve_analysis_id(analysis_id)
        if idx is None:
            return _empty_result()
        offsets = self._csr._offsets[idx: idx + 2]
        start, end = int(offsets[0]), int(offsets[1])
        if start == end:
            return _empty_result()
        vi = self._csr._variant_index[start:end].astype("int32")
        z = self._csr._z[start:end].astype("float32")
        se = self._csr._se[start:end].astype("float32")
        imp = self._get_imputed_slice(start, end)
        if observed_only:
            mask = imp == 0
            vi, z, se, imp = vi[mask], z[mask], se[mask], imp[mask]
        status = _status_array(imp, z)
        return {
            "variant_index": vi,
            "analysis_index": np.full(len(z), idx, dtype="int32"),
            "z": z,
            "se": se,
            "association_status": status,
        }

    def range_phewas(
        self,
        chromosome: str,
        start: int,
        end: int,
        *,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        """All associations where the variant falls in [start, end] (regional PheWAS)."""
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

        analysis_indices = np.searchsorted(offsets[1:], hit_positions, side="right").astype("int32")
        imp = (
            self._imputed[hit_positions].astype(np.uint8)
            if self._imputed is not None
            else np.zeros(len(hit_positions), dtype=np.uint8)
        )
        if observed_only:
            keep = imp == 0
            hit_positions = hit_positions[keep]
            analysis_indices = analysis_indices[keep]
            imp = imp[keep]

        z_out = z_all[hit_positions].astype("float32")
        se_out = se_all[hit_positions].astype("float32")
        return {
            "variant_index": vi_all[hit_positions].astype("int32"),
            "analysis_index": analysis_indices,
            "z": z_out,
            "se": se_out,
            "association_status": _status_array(imp, z_out),
        }

    def range_by_analysis(
        self,
        chromosome: str,
        start: int,
        end: int,
        *,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        """All associations for analyses whose probe/TSS falls in [start, end]."""
        trait_records = self._traits_reader.range(chromosome, start, end)
        if not trait_records:
            return _empty_result()

        all_vi: list[np.ndarray] = []
        all_ai: list[np.ndarray] = []
        all_z: list[np.ndarray] = []
        all_se: list[np.ndarray] = []
        all_status: list[np.ndarray] = []

        for rec in trait_records:
            ai = rec.analysis_index
            offsets = self._csr._offsets[ai: ai + 2]
            s, e = int(offsets[0]), int(offsets[1])
            if s == e:
                continue
            vi = self._csr._variant_index[s:e].astype("int32")
            z = self._csr._z[s:e].astype("float32")
            se = self._csr._se[s:e].astype("float32")
            imp = self._get_imputed_slice(s, e)
            if observed_only:
                keep = imp == 0
                vi, z, se, imp = vi[keep], z[keep], se[keep], imp[keep]
            if len(z) == 0:
                continue
            all_vi.append(vi)
            all_ai.append(np.full(len(z), ai, dtype="int32"))
            all_z.append(z)
            all_se.append(se)
            all_status.append(_status_array(imp, z))

        if not all_vi:
            return _empty_result()
        return {
            "variant_index": np.concatenate(all_vi),
            "analysis_index": np.concatenate(all_ai),
            "z": np.concatenate(all_z),
            "se": np.concatenate(all_se),
            "association_status": np.concatenate(all_status),
        }

    def phewas(self, identifier: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
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
        imp = (
            self._imputed[hit_positions].astype(np.uint8)
            if self._imputed is not None
            else np.zeros(len(hit_positions), dtype=np.uint8)
        )
        if observed_only:
            keep = imp == 0
            hit_positions = hit_positions[keep]
            analysis_indices = analysis_indices[keep]
            imp = imp[keep]

        z_out = z_all[hit_positions].astype("float32")
        se_out = se_all[hit_positions].astype("float32")
        return {
            "variant_index": np.full(len(hit_positions), target_vi, dtype="int32"),
            "analysis_index": analysis_indices,
            "z": z_out,
            "se": se_out,
            "association_status": _status_array(imp, z_out),
        }

    def top_hits(
        self,
        *,
        analysis_id: str | None = None,
        threshold: float = 5e-8,
        limit: int | None = None,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        """Associations passing a significance threshold.

        Uses the precomputed top-hit index when available (fast path);
        falls back to a full CSR scan otherwise.
        observed_only=True forces a full scan to exclude imputed associations.
        """
        key = threshold_key(threshold)
        root = zarr.open_group(str(self.store.path / "data.zarr"), mode="r")
        path = f"top_hits/{key}"
        if path in root:
            group = root[path]
            analysis_index = None
            if analysis_id is not None:
                analysis_index = self._resolve_analysis_id(analysis_id)
                if analysis_index is None or "analysis_offsets" not in group:
                    return _empty_result()
            reader = DenseTopHitReader(group)
            bounds = reader.bounds(analysis_index)
            vi = reader.read("variant_index", bounds, "int32")
            ai = reader.read("analysis_index", bounds, "int32")
            z = reader.read("z", bounds, "float32")
            se = reader.read("se", bounds, "float32")
            imp = (
                reader.read("imputed", bounds, "uint8")
                if "imputed" in group else np.zeros(len(vi), dtype=np.uint8)
            )
            if observed_only:
                keep = imp == 0
                vi, ai, z, se, imp = vi[keep], ai[keep], z[keep], se[keep], imp[keep]
            if limit is not None:
                vi, ai, z, se = vi[:limit], ai[:limit], z[:limit], se[:limit]
                imp = imp[:limit]
            return {
                "variant_index": vi, "analysis_index": ai, "z": z, "se": se,
                "association_status": _status_array(imp, z),
            }

        # Fallback: full scan
        import math
        offsets = self._csr._offsets[:]
        vi_all = self._csr._variant_index[:]
        z_all = self._csr._z[:]
        se_all = self._csr._se[:]

        sqrt2 = math.sqrt(2.0)
        lo, hi, mid = 0.0, 40.0, 0.0
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

        order = np.argsort(-np.abs(z_f32[hit_positions]))
        hit_positions = hit_positions[order]
        if limit is not None:
            hit_positions = hit_positions[:limit]

        analysis_indices = np.searchsorted(offsets[1:], hit_positions, side="right").astype("int32")
        imp_all = (
            self._imputed[:].astype(np.uint8) if self._imputed is not None
            else np.zeros(len(vi_all), dtype=np.uint8)
        )
        imp_hits = imp_all[hit_positions]
        if observed_only:
            keep = imp_hits == 0
            hit_positions = hit_positions[keep]
            analysis_indices = analysis_indices[keep]
            imp_hits = imp_hits[keep]

        z_out = z_f32[hit_positions]
        se_out = se_all[hit_positions].astype("float32")
        return {
            "variant_index": vi_all[hit_positions].astype("int32"),
            "analysis_index": analysis_indices,
            "z": z_out,
            "se": se_out,
            "association_status": _status_array(imp_hits, z_out),
        }

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
        *,
        observed_only: bool = False,
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
        all_vi, all_ai, all_z, all_se, all_status = [], [], [], [], []

        for aid in analysis_ids:
            idx = self._resolve_analysis_id(aid)
            if idx is None:
                continue
            offsets_pair = self._csr._offsets[idx: idx + 2]
            s, e = int(offsets_pair[0]), int(offsets_pair[1])
            if s == e:
                continue
            vi = self._csr._variant_index[s:e].astype("int32")
            z = self._csr._z[s:e].astype("float32")
            se = self._csr._se[s:e].astype("float32")
            imp = self._get_imputed_slice(s, e)
            sub_mask = np.isin(vi, np.array(sorted(target_vi), dtype=np.int32))
            if not sub_mask.any():
                continue
            vi, z, se, imp = vi[sub_mask], z[sub_mask], se[sub_mask], imp[sub_mask]
            if observed_only:
                keep = imp == 0
                vi, z, se, imp = vi[keep], z[keep], se[keep], imp[keep]
            if len(z) == 0:
                continue
            all_vi.append(vi)
            all_ai.append(np.full(len(z), idx, dtype="int32"))
            all_z.append(z)
            all_se.append(se)
            all_status.append(_status_array(imp, z))

        if not all_vi:
            return _empty_result()
        return {
            "variant_index": np.concatenate(all_vi),
            "analysis_index": np.concatenate(all_ai),
            "z": np.concatenate(all_z),
            "se": np.concatenate(all_se),
            "association_status": np.concatenate(all_status),
        }


def _concat_results(parts: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    """Concatenate query-result dicts (each a plain component read)."""
    parts = [p for p in parts if len(p["z"])]
    if not parts:
        return _empty_result()
    return {
        "variant_index": np.concatenate([p["variant_index"] for p in parts]).astype("int32"),
        "analysis_index": np.concatenate([p["analysis_index"] for p in parts]).astype("int32"),
        "z": np.concatenate([p["z"] for p in parts]).astype("float32"),
        "se": np.concatenate([p["se"] for p in parts]).astype("float32"),
        "association_status": np.concatenate([p["association_status"] for p in parts]),
    }


class HybridStoreQuery:
    """Query facade for a Hybrid store (ADR 0026).

    A thin integration layer: it dispatches to the nested Dense Component's
    ``StoreQuery`` (on-panel variants) and to the Ragged Overflow CSR (off-panel
    variants), remaps the Dense Component's panel-local ``variant_index`` onto the
    shared variant index space, and concatenates. A variant is in exactly one
    component (on-panel xor off-panel), so results are a plain union with no dedup.
    """

    def __init__(self, store: OpenGWASDBStore):
        self.store = store
        self._dense_store = open_store(dense_component_path(store.path))
        self._dense = StoreQuery(self._dense_store)
        self._dense_to_shared = np.load(dense_to_shared_path(store.path)).astype("int32")
        self._csr = RaggedCSRReader(store.path)  # overflow at store/data.zarr/ragged
        self._connection = connect(store.index_path)
        self._analyses = AnalysesIndex(store.path)  # shared analyses.tsv
        self._variant_axis = VariantAxis(store.path, self._connection)  # shared union table

    def close(self) -> None:
        self._dense.close()
        self._variant_axis.close()
        self._connection.close()

    def __enter__(self) -> HybridStoreQuery:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    # ── shared-table tables ──────────────────────────────────────────────────
    def analyses_table(self) -> dict[int, dict]:
        return self._dense.analyses_table()

    def variants_table(self) -> dict[int, dict]:
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

    # ── dispatch helpers ─────────────────────────────────────────────────────
    def _remap_dense(self, result: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Translate a Dense Component result's panel-local variant_index to shared."""
        if len(result["variant_index"]):
            result["variant_index"] = self._dense_to_shared[
                result["variant_index"].astype("int64")
            ].astype("int32")
        return result

    def _shared_is_on_panel(self, shared_idx: int) -> bool:
        pos = int(np.searchsorted(self._dense_to_shared, shared_idx))
        return pos < len(self._dense_to_shared) and int(self._dense_to_shared[pos]) == shared_idx

    def _overflow_for_analysis(self, col: int) -> dict[str, np.ndarray]:
        offsets = self._csr._offsets[col: col + 2]
        s, e = int(offsets[0]), int(offsets[1])
        if s == e:
            return _empty_result()
        vi = self._csr._variant_index[s:e].astype("int32")
        z = self._csr._z[s:e].astype("float32")
        se = self._csr._se[s:e].astype("float32")
        return {
            "variant_index": vi,
            "analysis_index": np.full(len(z), col, dtype="int32"),
            "z": z,
            "se": se,
            "association_status": _status_array(np.zeros(len(z), dtype=np.uint8), z),
        }

    def _overflow_by_variants(self, shared_indices: set[int]) -> dict[str, np.ndarray]:
        """All overflow associations whose (off-panel) variant is in the set."""
        if not shared_indices:
            return _empty_result()
        offsets = self._csr._offsets[:]
        vi_all = self._csr._variant_index[:]
        wanted = np.fromiter(shared_indices, dtype=np.int32, count=len(shared_indices))
        mask = np.isin(vi_all, wanted)
        hits = np.where(mask)[0]
        if len(hits) == 0:
            return _empty_result()
        analysis_indices = np.searchsorted(offsets[1:], hits, side="right").astype("int32")
        z = self._csr._z[:][hits].astype("float32")
        se = self._csr._se[:][hits].astype("float32")
        return {
            "variant_index": vi_all[hits].astype("int32"),
            "analysis_index": analysis_indices,
            "z": z,
            "se": se,
            "association_status": _status_array(np.zeros(len(hits), dtype=np.uint8), z),
        }

    # ── public query surface ─────────────────────────────────────────────────
    def analysis(self, analysis_id: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
        dense = self._remap_dense(self._dense.analysis(analysis_id, observed_only=observed_only))
        analysis = self._analyses.by_id(analysis_id)
        overflow = (
            self._overflow_for_analysis(int(analysis["analysis_index"]))
            if analysis is not None else _empty_result()
        )
        return _concat_results([dense, overflow])

    def phewas(self, identifier: str, *, observed_only: bool = False) -> dict[str, np.ndarray]:
        variant = self._variant_axis.by_identifier(identifier)
        if variant is None:
            return _empty_result()
        if self._shared_is_on_panel(variant.variant_index):
            return self._remap_dense(
                self._dense.phewas(variant.alid, observed_only=observed_only)
            )
        return self._overflow_by_variants({int(variant.variant_index)})

    def range_phewas(
        self, chromosome: str, start: int, end: int, *, observed_only: bool = False
    ) -> dict[str, np.ndarray]:
        dense = self._remap_dense(
            self._dense.range_phewas(chromosome, start, end, observed_only=observed_only)
        )
        shared_idx = self._variant_axis.range_indices(chromosome, start, end)
        off_panel = {
            int(i) for i in shared_idx.tolist() if not self._shared_is_on_panel(int(i))
        }
        overflow = self._overflow_by_variants(off_panel)
        return _concat_results([dense, overflow])

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
        *,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        dense = self._remap_dense(
            self._dense.lookup(identifiers, analysis_ids, observed_only=observed_only)
        )
        # Off-panel identifiers: resolve on the shared table, keep off-panel ones.
        off_shared: set[int] = set()
        for id_ in identifiers:
            rec = self._variant_axis.by_identifier(id_)
            if rec is not None and not self._shared_is_on_panel(rec.variant_index):
                off_shared.add(int(rec.variant_index))
        wanted_cols = {
            int(a["analysis_index"])
            for aid in analysis_ids
            if (a := self._analyses.by_id(aid)) is not None
        }
        overflow = self._overflow_by_variants(off_shared)
        if len(overflow["z"]) and wanted_cols:
            keep = np.isin(overflow["analysis_index"], np.fromiter(
                wanted_cols, dtype="int32", count=len(wanted_cols)))
            overflow = {k: v[keep] for k, v in overflow.items()}
        elif not wanted_cols:
            overflow = _empty_result()
        return _concat_results([dense, overflow])

    def _overflow_top_hits(
        self, threshold: float, analysis_index: int | None = None
    ) -> dict[str, np.ndarray]:
        key = threshold_key(threshold)
        root = zarr.open_group(str(self.store.path / "data.zarr"), mode="r")
        path = f"top_hits/{key}"
        if path not in root:
            return _empty_result()
        group = root[path]
        if analysis_index is not None and "analysis_offsets" not in group:
            return _empty_result()
        reader = DenseTopHitReader(group)
        bounds = reader.bounds(analysis_index)
        vi = reader.read("variant_index", bounds, "int32")
        ai = reader.read("analysis_index", bounds, "int32")
        z = reader.read("z", bounds, "float32")
        se = reader.read("se", bounds, "float32")
        return {
            "variant_index": vi,
            "analysis_index": ai,
            "z": z,
            "se": se,
            "association_status": _status_array(np.zeros(len(z), dtype=np.uint8), z),
        }

    def top_hits(
        self,
        *,
        analysis_id: str | None = None,
        threshold: float = 5e-8,
        limit: int | None = None,
        observed_only: bool = False,
    ) -> dict[str, np.ndarray]:
        analysis_index = None
        if analysis_id is not None:
            analysis = self._analyses.by_id(analysis_id)
            if analysis is None:
                return _empty_result()
            analysis_index = int(analysis["analysis_index"])
        dense = self._remap_dense(self._dense.top_hits(
            analysis_id=analysis_id, threshold=threshold, observed_only=observed_only
        ))
        overflow = self._overflow_top_hits(
            threshold, analysis_index
        )  # overflow is always observed
        merged = _concat_results([dense, overflow])
        if len(merged["z"]):
            order = np.lexsort((merged["variant_index"], merged["analysis_index"]))
            merged = {k: v[order] for k, v in merged.items()}
        if limit is not None:
            merged = {k: v[:limit] for k, v in merged.items()}
        return merged


def query_store(path: str | Path) -> StoreQuery | RaggedStoreQuery | HybridStoreQuery:
    """Open a store and return the layout-independent query facade."""
    store = open_store(path)
    if store.manifest.primary_layout is PrimaryStorageLayout.RAGGED:
        return RaggedStoreQuery(store)
    if store.manifest.primary_layout is PrimaryStorageLayout.HYBRID:
        return HybridStoreQuery(store)
    return StoreQuery(store)
