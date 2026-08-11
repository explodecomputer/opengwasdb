"""Dense Observed-Only writer."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.build.source import NormalisedAssociation
from opengwasdb.index import connect, create_lookup_indexes, initialise_schema, set_metadata
from opengwasdb.layouts.dense.constants import (
    DEFAULT_CHUNK_SHAPE,
    DEFAULT_COMPRESSOR,
    DEFAULT_DTYPE,
)
from opengwasdb.layouts.dense.overview import write_overview_html
from opengwasdb.layouts.dense.top_hits import build_top_hit_indexes, read_top_hit_counts
from opengwasdb.model.analyses import (
    ANCESTRY_PROP_PREFIX,
    SHARED_CORE_COLUMNS,
    STORE_ONLY_COLUMNS,
    AnalysesTable,
    write_analyses,
)
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
    """One Analysis's Analytical Metadata (store-format spec §7a, ADR 0030).

    Every field beyond the original five is optional and blank by default:
    this dataclass is shared by every build path, and most of them (the
    manifest-free tiny builder, a bare VCF manifest with none of the optional
    columns) genuinely have no source for sample size, ancestry, or SD
    provenance -- writing "" for an unresolved column is honest, not a
    fabrication, and matches ``opengwasdb.model.analyses``'s own tolerance for
    blank shared-core columns.
    """

    analysis_id: str
    phenotype_id: str | None
    phenotype_label: str | None
    analysis_label: str | None
    stored_effect_scale: str
    assigned_ancestry: str = ""
    ancestry_assignment_method: str = ""
    ancestry_prop: dict[str, str] = field(default_factory=dict)  # {population: proportion}
    sample_size: str = ""
    original_sd_method: str = ""
    original_sd: str = ""
    completed_against: str = ""
    completion_median_pearson_r: str = ""
    completion_n_imputed_total: str = ""
    completion_n_missing_total: str = ""
    n_hits_5e8: str = ""
    n_hits_5e6: str = ""
    n_hits_5e4: str = ""


# SHARED_CORE_COLUMNS/STORE_ONLY_COLUMNS come from opengwasdb.model.analyses
# (issue #16) rather than being retyped here, so the two never drift.
# ancestry_prop_<population> is a dynamic column set, not part of either fixed
# tuple (opengwasdb.model.analyses.classify_column matches it by prefix) --
# _fieldnames_for below splices in whichever populations the analyses being
# written actually carry, at the position store-format spec §7a declares.
_SAMPLE_SIZE_KIND_INDEX = SHARED_CORE_COLUMNS.index("sample_size_kind")


def _fieldnames_for(analyses: list[AnalysisMetadata]) -> tuple[str, ...]:
    prop_populations = sorted({pop for a in analyses for pop in a.ancestry_prop})
    return (
        *SHARED_CORE_COLUMNS[:_SAMPLE_SIZE_KIND_INDEX],
        *(f"{ANCESTRY_PROP_PREFIX}{pop}" for pop in prop_populations),
        *SHARED_CORE_COLUMNS[_SAMPLE_SIZE_KIND_INDEX:],
        *STORE_ONLY_COLUMNS,
    )


def _analysis_metadata_row(
    index: int, a: AnalysisMetadata, fieldnames: tuple[str, ...]
) -> dict[str, str]:
    row = {
        "analysis_index": str(index),
        "analysis_id": a.analysis_id,
        "phenotype_id": a.phenotype_id or "",
        "phenotype_label": a.phenotype_label or "",
        "analysis_label": a.analysis_label or "",
        "stored_effect_scale": a.stored_effect_scale,
        "assigned_ancestry": a.assigned_ancestry,
        "ancestry_assignment_method": a.ancestry_assignment_method,
        "sample_size_kind": "",
        "sample_size_scope": "analysis_level" if a.sample_size else "",
        "sample_size": a.sample_size,
        "n_cases": "",
        "n_controls": "",
        "original_effect_scale": "",
        "original_sd": a.original_sd,
        "original_sd_method": a.original_sd_method,
        "original_sd_dispersion": "",
        "completed_against": a.completed_against,
        "completion_median_pearson_r": a.completion_median_pearson_r,
        "completion_n_imputed_total": a.completion_n_imputed_total,
        "completion_n_missing_total": a.completion_n_missing_total,
        "n_hits_5e8": a.n_hits_5e8,
        "n_hits_5e6": a.n_hits_5e6,
        "n_hits_5e4": a.n_hits_5e4,
    }
    for column in fieldnames:
        if column.startswith(ANCESTRY_PROP_PREFIX):
            row[column] = a.ancestry_prop.get(column[len(ANCESTRY_PROP_PREFIX):], "")
    return row


def analysis_metadata_from_row(row: dict[str, str]) -> AnalysisMetadata:
    """The inverse of :func:`_analysis_metadata_row`: reconstruct an
    ``AnalysisMetadata`` from one already-written ``analyses.tsv`` row.

    Shared by the dense and hybrid completion pipelines (issue #22), which
    both need to carry a source store's Analytical Metadata forward into a
    freshly written destination ``analyses.tsv`` -- callers that only need to
    override a few completion-specific fields (``completed_against``, the
    ``completion_*`` rollups) do so with :func:`dataclasses.replace` on the
    result rather than reconstructing the whole row by hand.
    """
    return AnalysisMetadata(
        analysis_id=row["analysis_id"],
        phenotype_id=row.get("phenotype_id") or None,
        phenotype_label=row.get("phenotype_label") or None,
        analysis_label=row.get("analysis_label") or None,
        stored_effect_scale=row["stored_effect_scale"],
        assigned_ancestry=row.get("assigned_ancestry", ""),
        ancestry_assignment_method=row.get("ancestry_assignment_method", ""),
        ancestry_prop={
            key[len(ANCESTRY_PROP_PREFIX):]: value
            for key, value in row.items()
            if key.startswith(ANCESTRY_PROP_PREFIX) and value
        },
        sample_size=row.get("sample_size", ""),
        original_sd_method=row.get("original_sd_method", ""),
        original_sd=row.get("original_sd", ""),
        completed_against=row.get("completed_against", ""),
        completion_median_pearson_r=row.get("completion_median_pearson_r", ""),
        completion_n_imputed_total=row.get("completion_n_imputed_total", ""),
        completion_n_missing_total=row.get("completion_n_missing_total", ""),
        n_hits_5e8=row.get("n_hits_5e8", ""),
        n_hits_5e6=row.get("n_hits_5e6", ""),
        n_hits_5e4=row.get("n_hits_5e4", ""),
    )


def add_hit_counts(
    store_path: str | Path, analyses: list[AnalysisMetadata]
) -> list[AnalysisMetadata]:
    """Add per-Analysis Top-Hit Counts from ``store_path``'s already-built
    top-hit index onto whatever each ``AnalysisMetadata`` already carries
    (ADR 0032).

    Blank ``n_hits_*`` fields act as zero, so this both sets counts from
    scratch (a fresh build's analyses start blank) and accumulates a Hybrid
    store's Dense Component and Ragged Overflow Component contributions onto
    the same Analysis when called once per component -- the two partition an
    Analysis's associations disjointly (CONTEXT.md, Query Component).
    """
    counts = read_top_hit_counts(store_path, len(analyses))
    return [
        replace(
            a,
            n_hits_5e8=str(int(a.n_hits_5e8 or 0) + counts["n_hits_5e8"][i]),
            n_hits_5e6=str(int(a.n_hits_5e6 or 0) + counts["n_hits_5e6"][i]),
            n_hits_5e4=str(int(a.n_hits_5e4 or 0) + counts["n_hits_5e4"][i]),
        )
        for i, a in enumerate(analyses)
    ]


def build_analyses_table(analyses: list[AnalysisMetadata]) -> AnalysesTable:
    fieldnames = _fieldnames_for(analyses)
    return AnalysesTable(
        fieldnames=fieldnames,
        rows=tuple(_analysis_metadata_row(i, a, fieldnames) for i, a in enumerate(analyses)),
    )


def write_analyses_tsv(output_path: Path, analyses: list[AnalysisMetadata]) -> Path:
    """Write ``analyses.tsv`` + ``overview.html`` -- the sole source of truth
    for Analytical Metadata (ADR 0030, issue #22), replacing the old SQLite
    ``analyses`` table and its own generated human-browsable rendering.
    """
    output_path = Path(output_path)
    table = build_analyses_table(analyses)
    path = write_analyses(output_path / "analyses.tsv", table)
    write_overview_html(output_path, table)
    return path


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
        write_analyses_tsv(work, add_hit_counts(work, analyses))
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
