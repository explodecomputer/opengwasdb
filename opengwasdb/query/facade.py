"""Layout-independent query facade."""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.index import analysis_by_id, connect
from opengwasdb.layouts.dense.top_hits import threshold_key
from opengwasdb.model.enums import PrimaryStorageLayout
from opengwasdb.query.results import AssociationResult
from opengwasdb.store.open import OpenGWASDBStore, open_store
from opengwasdb.variants import VariantAxis, VariantRecord


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

    def variant(self, identifier: str) -> list[AssociationResult]:
        """Return all finite associations for a canonical ALID or alias."""

        variant = self._variant_axis.by_identifier(identifier)
        if variant is None:
            return []
        return self._results_for_variant(variant)

    def range(
        self,
        chromosome: str,
        start: int,
        end: int,
        *,
        analysis_id: str | None = None,
    ) -> list[AssociationResult]:
        """Return finite associations in a genomic range."""

        variants = self._variant_axis.range(chromosome, start, end)
        if analysis_id is None:
            return self._results_for_variants(variants)
        analysis = analysis_by_id(self._connection, analysis_id)
        if analysis is None:
            return []
        col = int(analysis["analysis_index"])
        results: list[AssociationResult] = []
        for variant in variants:
            z = float(self._root["z"][variant.variant_index, col])
            se = float(self._root["se"][variant.variant_index, col])
            if math.isfinite(z) and math.isfinite(se):
                results.append(AssociationResult.from_rows(variant, analysis, z, se))
        return results

    def analysis(self, analysis_id: str) -> list[AssociationResult]:
        """Return all finite associations for one analysis."""

        analysis = analysis_by_id(self._connection, analysis_id)
        if analysis is None:
            return []
        col = int(analysis["analysis_index"])
        z_col = self._root["z"][:, col]
        se_col = self._root["se"][:, col]
        finite_rows = [
            row
            for row, (z, se) in enumerate(zip(z_col, se_col, strict=True))
            if math.isfinite(float(z)) and math.isfinite(float(se))
        ]
        variants = self._variant_axis.by_indices(finite_rows)
        return [
            AssociationResult.from_rows(
                variants[row],
                analysis,
                float(z_col[row]),
                float(se_col[row]),
            )
            for row in finite_rows
            if row in variants
        ]

    def analysis_arrays(self, analysis_id: str) -> dict[str, object] | None:
        """Return dense statistic arrays for one analysis without row materialisation.

        This is the efficient full-study extraction path for dense stores. The
        row-materialised `analysis()` method is useful for small extracts and
        tests, but millions of associations should use array/streaming results.
        """

        analysis = analysis_by_id(self._connection, analysis_id)
        if analysis is None:
            return None
        col = int(analysis["analysis_index"])
        return {
            "analysis_id": str(analysis["analysis_id"]),
            "analysis_index": col,
            "z": self._root["z"][:, col],
            "se": self._root["se"][:, col],
        }

    def phewas(self, identifier: str) -> list[AssociationResult]:
        """Return one variant across all analyses."""

        return self.variant(identifier)

    def lookup(
        self,
        identifiers: list[str],
        analysis_ids: list[str],
    ) -> list[AssociationResult]:
        """Return finite associations for a specific variant x analysis set."""

        variants = [
            variant
            for identifier in identifiers
            if (variant := self._variant_axis.by_identifier(identifier)) is not None
        ]
        analyses = [
            analysis
            for analysis_id in analysis_ids
            if (analysis := analysis_by_id(self._connection, analysis_id)) is not None
        ]
        if not variants or not analyses:
            return []
        row_indices = [variant.variant_index for variant in variants]
        col_indices = [int(analysis["analysis_index"]) for analysis in analyses]
        z_block = self._root["z"].oindex[row_indices, :][:, col_indices]
        se_block = self._root["se"].oindex[row_indices, :][:, col_indices]
        results: list[AssociationResult] = []
        for row_offset, variant in enumerate(variants):
            for col_offset, analysis in enumerate(analyses):
                z = float(z_block[row_offset, col_offset])
                se = float(se_block[row_offset, col_offset])
                if math.isfinite(z) and math.isfinite(se):
                    results.append(AssociationResult.from_rows(variant, analysis, z, se))
        return results

    def top_hits(
        self,
        *,
        threshold: float = 5e-8,
        limit: int | None = None,
    ) -> list[AssociationResult]:
        """Return ranked top-hit associations using the dense top-hit index."""

        key = threshold_key(threshold)
        path = f"top_hits/{key}"
        if path not in self._root:
            return []
        group = self._root[path]
        variant_indices = group["variant_index"][:]
        analysis_indices = group["analysis_index"][:]
        z_values = group["z"][:]
        se_values = group["se"][:] if "se" in group else None
        if limit is not None:
            variant_indices = variant_indices[:limit]
            analysis_indices = analysis_indices[:limit]
            z_values = z_values[:limit]
            if se_values is not None:
                se_values = se_values[:limit]
        variants = self._variant_axis.by_indices([int(row) for row in variant_indices])
        analyses = self._analyses_for_indices([int(col) for col in analysis_indices])
        if se_values is None:
            unique_rows = np.unique(variant_indices).astype("int64")
            row_offsets = {int(row): i for i, row in enumerate(unique_rows)}
            se_block = self._root["se"].oindex[unique_rows, :]
            se_values = np.array(
                [
                    float(se_block[row_offsets[int(row)], int(col)])
                    for row, col in zip(variant_indices, analysis_indices, strict=True)
                ],
                dtype="float32",
            )
        return [
            AssociationResult.from_rows(
                variants[int(row)],
                analyses[int(col)],
                float(z),
                float(se),
            )
            for row, col, z, se in zip(
                variant_indices,
                analysis_indices,
                z_values,
                se_values,
                strict=True,
            )
            if int(row) in variants and int(col) in analyses
        ]

    def _results_for_variant(self, variant: VariantRecord) -> list[AssociationResult]:
        z_row = self._root["z"][variant.variant_index, :]
        se_row = self._root["se"][variant.variant_index, :]
        analyses = self._analyses_by_index()
        return [
            AssociationResult.from_rows(variant, analyses[col], float(z), float(se))
            for col, (z, se) in enumerate(zip(z_row, se_row, strict=True))
            if col in analyses and math.isfinite(float(z)) and math.isfinite(float(se))
        ]

    def _results_for_variants(self, variants: list[VariantRecord]) -> list[AssociationResult]:
        if not variants:
            return []
        rows = [variant.variant_index for variant in variants]
        z_block = self._root["z"].oindex[rows, :]
        se_block = self._root["se"].oindex[rows, :]
        analyses = self._analyses_by_index()
        results: list[AssociationResult] = []
        for row_offset, variant in enumerate(variants):
            for col, analysis in analyses.items():
                z = float(z_block[row_offset, col])
                se = float(se_block[row_offset, col])
                if math.isfinite(z) and math.isfinite(se):
                    results.append(AssociationResult.from_rows(variant, analysis, z, se))
        return results

    def _analyses_by_index(self) -> dict[int, sqlite3.Row]:
        rows = self._connection.execute("SELECT * FROM analyses ORDER BY analysis_index").fetchall()
        return {int(row["analysis_index"]): row for row in rows}

    def _analyses_for_indices(self, indices: list[int]) -> dict[int, sqlite3.Row]:
        unique = sorted(set(indices))
        if not unique:
            return {}
        placeholders = ",".join("?" for _ in unique)
        rows = self._connection.execute(
            f"SELECT * FROM analyses WHERE analysis_index IN ({placeholders})",
            unique,
        ).fetchall()
        return {int(row["analysis_index"]): row for row in rows}


def query_store(path: str | Path) -> StoreQuery:
    """Open a store and return the layout-independent query facade."""

    return StoreQuery(open_store(path))
