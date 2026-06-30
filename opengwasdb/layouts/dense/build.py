"""Dense Observed-Only writer."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.build.source import NormalisedAssociation
from opengwasdb.index import connect, initialise_schema, set_metadata
from opengwasdb.layouts.dense.constants import (
    DEFAULT_CHUNK_SHAPE,
    DEFAULT_COMPRESSOR,
    DEFAULT_DTYPE,
)
from opengwasdb.layouts.dense.top_hits import build_top_hit_indexes
from opengwasdb.model.enums import AssociationCoverage, CompletionState, PrimaryStorageLayout
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.variants import (
    VARIANT_AXIS_FORMAT,
    VARIANT_OFFSETS_FILENAME,
    VARIANT_TABIX_FILENAME,
    VARIANT_TABLE_FILENAME,
    CanonicalVariant,
    chromosome_sort_key,
    write_variant_axis,
)


@dataclass(frozen=True)
class DenseBuildResult:
    """Paths and dimensions for a built Dense Observed-Only store."""

    output_path: Path
    n_variants: int
    n_analyses: int


@dataclass(frozen=True)
class AnalysisMetadata:
    analysis_id: str
    phenotype_id: str | None
    phenotype_label: str | None
    analysis_label: str | None
    stored_effect_scale: str


def build_dense_observed_store(
    records: list[NormalisedAssociation],
    output_path: str | Path,
    *,
    store_id: str,
    release_id: str,
    reference_assembly: str,
    chunk_shape: tuple[int, int] = DEFAULT_CHUNK_SHAPE,
    dtype: str = DEFAULT_DTYPE,
    overwrite: bool = False,
) -> DenseBuildResult:
    """Write a Dense Observed-Only Store Release from normalised associations."""

    if not records:
        raise ValueError("cannot build a store with no association records")

    out = Path(output_path)
    if out.exists():
        if not overwrite:
            raise FileExistsError(f"output path already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    work = out.with_name(f".{out.name}.tmp")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    try:
        variants = _collect_variants(records)
        analyses = _collect_analyses(records)
        variant_index = {variant.alid: i for i, variant in enumerate(variants)}
        analysis_index = {analysis.analysis_id: i for i, analysis in enumerate(analyses)}

        z = np.full((len(variants), len(analyses)), np.nan, dtype=dtype)
        se = np.full((len(variants), len(analyses)), np.nan, dtype=dtype)

        seen_cells: set[tuple[int, int]] = set()
        for record in records:
            row = variant_index[record.variant.alid]
            col = analysis_index[record.analysis_id]
            cell = (row, col)
            if cell in seen_cells:
                raise ValueError(
                    f"duplicate association for variant {record.variant.alid} "
                    f"and analysis {record.analysis_id}"
                )
            seen_cells.add(cell)
            z[row, col] = record.z
            se[row, col] = record.se

        rsid_by_alid = _first_rsids_by_alid(records)
        _write_manifest(work, store_id, release_id, reference_assembly, records, chunk_shape, dtype)
        write_variant_axis(work, variants, rsid_by_alid)
        _write_index(work, variants, analyses, records, chunk_shape, dtype)
        _write_zarr(work, z, se, chunk_shape, dtype)
        build_top_hit_indexes(work)
        if out.exists():
            shutil.rmtree(out)
        work.rename(out)
        return DenseBuildResult(output_path=out, n_variants=len(variants), n_analyses=len(analyses))
    except Exception:
        shutil.rmtree(work, ignore_errors=True)
        raise


def _collect_variants(records: list[NormalisedAssociation]) -> list[CanonicalVariant]:
    by_alid = {record.variant.alid: record.variant for record in records}
    return sorted(
        by_alid.values(),
        key=lambda variant: (
            chromosome_sort_key(variant.chromosome),
            variant.position,
            variant.effect_allele,
            variant.other_allele,
        ),
    )


def _collect_analyses(records: list[NormalisedAssociation]) -> list[AnalysisMetadata]:
    by_id: dict[str, AnalysisMetadata] = {}
    for record in records:
        existing = by_id.get(record.analysis_id)
        current = AnalysisMetadata(
            analysis_id=record.analysis_id,
            phenotype_id=record.phenotype_id,
            phenotype_label=record.phenotype_label,
            analysis_label=record.analysis_label,
            stored_effect_scale=record.stored_effect_scale.value,
        )
        if existing is None:
            by_id[record.analysis_id] = current
            continue
        if existing.stored_effect_scale != current.stored_effect_scale:
            raise ValueError(f"analysis {record.analysis_id} has mixed stored_effect_scale values")
    return [by_id[key] for key in sorted(by_id)]


def _write_manifest(
    output_path: Path,
    store_id: str,
    release_id: str,
    reference_assembly: str,
    records: list[NormalisedAssociation],
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    manifest = StoreManifest(
        store_id=store_id,
        release_id=release_id,
        format_version="0.1",
        primary_layout=PrimaryStorageLayout.DENSE,
        association_coverage=AssociationCoverage.FULL,
        completion_state=CompletionState.OBSERVED_ONLY,
        reference_assembly=reference_assembly,
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            "builder": "opengwasdb.v0.1_dense_observed",
            "source_record_count": len(records),
            "dense": {
                "statistic_arrays": ["z", "se"],
                "dtype": dtype,
                "chunk_shape": list(chunk_shape),
                "compressor": DEFAULT_COMPRESSOR,
                "top_hit_thresholds": [5e-8, 5e-6, 5e-4],
                "variant_axis": {
                    "format": VARIANT_AXIS_FORMAT,
                    "table": VARIANT_TABLE_FILENAME,
                    "tabix_index": VARIANT_TABIX_FILENAME,
                    "row_offsets": VARIANT_OFFSETS_FILENAME,
                },
            },
        },
    )
    (output_path / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_index(
    output_path: Path,
    variants: list[CanonicalVariant],
    analyses: list[AnalysisMetadata],
    records: list[NormalisedAssociation],
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    with connect(output_path / "index.sqlite") as connection:
        initialise_schema(connection)
        set_metadata(connection, "schema_version", 2)
        set_metadata(connection, "n_variants", len(variants))
        set_metadata(connection, "n_analyses", len(analyses))
        set_metadata(
            connection,
            "dense",
            {
                "dtype": dtype,
                "chunk_shape": list(chunk_shape),
                "compressor": DEFAULT_COMPRESSOR,
                "variant_axis": {
                    "format": VARIANT_AXIS_FORMAT,
                    "table": VARIANT_TABLE_FILENAME,
                    "tabix_index": VARIANT_TABIX_FILENAME,
                    "row_offsets": VARIANT_OFFSETS_FILENAME,
                },
            },
        )
        variant_indices = {variant.alid: i for i, variant in enumerate(variants)}
        aliases: set[tuple[str, int]] = set()
        for record in records:
            if record.rsid:
                aliases.add((record.rsid, variant_indices[record.variant.alid]))
        connection.executemany(
            "INSERT OR IGNORE INTO variant_aliases(alias, variant_index) VALUES (?, ?)",
            sorted(aliases),
        )
        connection.executemany(
            """
            INSERT INTO analyses(
                analysis_index, analysis_id, phenotype_id, phenotype_label,
                analysis_label, stored_effect_scale
            )
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    i,
                    analysis.analysis_id,
                    analysis.phenotype_id,
                    analysis.phenotype_label,
                    analysis.analysis_label,
                    analysis.stored_effect_scale,
                )
                for i, analysis in enumerate(analyses)
            ],
        )
        connection.commit()


def _first_rsids_by_alid(records: list[NormalisedAssociation]) -> dict[str, str]:
    rsid_by_alid: dict[str, str] = {}
    for record in records:
        if record.rsid:
            rsid_by_alid.setdefault(record.variant.alid, record.rsid)
    return rsid_by_alid


def _write_zarr(
    output_path: Path,
    z: np.ndarray,
    se: np.ndarray,
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    # Clip chunk shape to array dimensions so zarr's declared shape matches what
    # is physically stored — oversized chunks cause zarr to allocate a large
    # decompression buffer even when the array is narrower than chunk_shape[1].
    effective_chunks = (min(chunk_shape[0], z.shape[0]), min(chunk_shape[1], z.shape[1]))
    root = zarr.open_group(str(output_path / "data.zarr"), mode="w")
    root.create_dataset("z", data=z, chunks=effective_chunks, compressor=compressor, dtype=dtype)
    root.create_dataset("se", data=se, chunks=effective_chunks, compressor=compressor, dtype=dtype)
    root.attrs["layout"] = "dense"
    root.attrs["completion_state"] = "observed_only"
    root.attrs["compressor"] = DEFAULT_COMPRESSOR
    root.attrs["chunk_shape"] = list(effective_chunks)
