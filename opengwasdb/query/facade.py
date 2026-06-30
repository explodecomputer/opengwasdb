"""Layout-independent query facade."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import zarr

from opengwasdb.index import analysis_by_id, connect
from opengwasdb.layouts.dense.top_hits import threshold_key
from opengwasdb.model.enums import PrimaryStorageLayout
from opengwasdb.store.open import OpenGWASDBStore, open_store
from opengwasdb.variants import VariantAxis


def _empty_result() -> dict[str, np.ndarray]:
    return {
        "variant_index": np.empty(0, dtype="int32"),
        "analysis_index": np.empty(0, dtype="int32"),
        "z": np.empty(0, dtype="float32"),
        "se": np.empty(0, dtype="float32"),
    }


class StoreQuery:
    """Public query object that hides the physical store layout."""

    def __init__(self, store: OpenGWASDBStore):
        self.store = store
        if store.manifest.primary_layout is not PrimaryStorageLayout.DENSE:
            raise NotImplementedError("only dense layout is implemented in v0.1")
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
        variants = self._variant_axis.range(chromosome, start, end)
        if not variants:
            return _empty_result()
        row_indices = [v.variant_index for v in variants]
        z_block = self._root["z"].oindex[row_indices, :].astype("float32")
        se_block = self._root["se"].oindex[row_indices, :].astype("float32")
        mask = np.isfinite(z_block) & np.isfinite(se_block)
        rows_rel, cols = np.where(mask)
        return {
            "variant_index": np.array([row_indices[r] for r in rows_rel], dtype="int32"),
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


def query_store(path: str | Path) -> StoreQuery:
    """Open a store and return the layout-independent query facade."""
    return StoreQuery(open_store(path))
