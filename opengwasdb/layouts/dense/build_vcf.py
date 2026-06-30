"""Two-pass Dense Observed-Only writer from GWAS-VCF manifests with inline liftover."""

from __future__ import annotations

import csv
import json
import logging
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.build.liftover import LiftoverFailureError, build_liftover_lookup
from opengwasdb.build.vcf_source import (
    read_vcf_study_type,
    stream_vcf_associations,
    stream_vcf_variants,
)
from opengwasdb.index import connect, initialise_schema, set_metadata
from opengwasdb.layouts.dense.build import AnalysisMetadata, DenseBuildResult
from opengwasdb.layouts.dense.constants import (
    DEFAULT_CHUNK_SHAPE,
    DEFAULT_COMPRESSOR,
    DEFAULT_DTYPE,
)
from opengwasdb.layouts.dense.top_hits import build_top_hit_indexes
from opengwasdb.model.enums import AssociationCoverage, CompletionState, PrimaryStorageLayout
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.variants import CanonicalVariant, write_variant_axis
from opengwasdb.variants.normalise import chromosome_sort_key

log = logging.getLogger(__name__)

__all__ = ["build_dense_from_vcf_manifest", "LiftoverFailureError"]


@dataclass(frozen=True)
class _ManifestRow:
    trait_id: str
    file_path: str
    trait_name: str
    n: int


def build_dense_from_vcf_manifest(
    manifest_path: str | Path,
    output_path: str | Path,
    *,
    chain_file: str | Path | None = None,
    store_id: str,
    release_id: str,
    liftover_failure_threshold: float = 0.01,
    chunk_shape: tuple[int, int] = DEFAULT_CHUNK_SHAPE,
    dtype: str = DEFAULT_DTYPE,
    overwrite: bool = False,
) -> DenseBuildResult:
    """Build a Dense Observed-Only Store from a manifest of GWAS-VCF files.

    VCF files are assumed to be in GRCh37/hg19 coordinates.  All variant
    positions are lifted to GRCh38/hg38 inline; the output store uses hg38
    coordinates.

    Two-pass streaming: Pass 1 collects the union variant set and runs liftover
    once.  Pass 2 fills zarr columns one analysis at a time.  The full
    association list is never materialised in memory.

    Parameters
    ----------
    manifest_path:
        TSV with columns ``trait_id``, ``file_path``, ``trait_name``, ``n``.
    output_path:
        Destination directory for the store.
    chain_file:
        Optional path to a pyliftover chain file.  When None, pyliftover
        downloads the hg19→hg38 chain automatically.
    store_id / release_id:
        Identifiers written to ``manifest.json``.
    liftover_failure_threshold:
        Maximum fraction of variants allowed to fail liftover (default 0.01).
        Raises ``LiftoverFailureError`` if exceeded.
    """
    manifest_rows = _read_manifest(manifest_path)
    if not manifest_rows:
        raise ValueError(f"manifest {manifest_path} contains no rows")

    out = Path(output_path)
    if out.exists():
        if not overwrite:
            raise FileExistsError(f"output path already exists: {out}")
        shutil.rmtree(out)
    out.mkdir(parents=True)

    # ------------------------------------------------------------------
    # Pass 1: collect union variant set across all VCFs
    # ------------------------------------------------------------------
    log.info("Pass 1: collecting union variant set from %d VCFs", len(manifest_rows))
    hg19_tuples: set[tuple[str, int, str, str]] = set()
    for row in manifest_rows:
        for variant in stream_vcf_variants(row.file_path):
            hg19_tuples.add(variant)
    log.info("Pass 1 complete: %d unique hg19 variants", len(hg19_tuples))

    # ------------------------------------------------------------------
    # Liftover: hg19 → hg38 (single LiftOver object for entire batch)
    # ------------------------------------------------------------------
    log.info("Running liftover hg19 → hg38 (%d variants)", len(hg19_tuples))
    hg19_lookup = build_liftover_lookup(
        hg19_tuples,
        from_build="hg19",
        to_build="hg38",
        failure_threshold=liftover_failure_threshold,
        chain_file=chain_file,
    )
    log.info("Liftover complete: %d variants mapped", len(hg19_lookup))

    # Sort hg38 ALIDs by (chromosome, position, a1, a2)
    hg38_alids = sorted(set(hg19_lookup.values()), key=_alid_sort_key)
    n_variants = len(hg38_alids)
    n_analyses = len(manifest_rows)
    variant_index: dict[str, int] = {alid: i for i, alid in enumerate(hg38_alids)}
    analysis_index: dict[str, int] = {row.trait_id: i for i, row in enumerate(manifest_rows)}

    # ------------------------------------------------------------------
    # Read study types (lightweight header scan per VCF)
    # ------------------------------------------------------------------
    analyses: list[AnalysisMetadata] = []
    for row in manifest_rows:
        stored_effect_scale = read_vcf_study_type(row.file_path)
        analyses.append(
            AnalysisMetadata(
                analysis_id=row.trait_id,
                phenotype_id=row.trait_id,
                phenotype_label=row.trait_name,
                analysis_label=row.trait_id,
                stored_effect_scale=stored_effect_scale.value,
            )
        )

    # ------------------------------------------------------------------
    # Write SQLite index + tabix variant axis
    # ------------------------------------------------------------------
    _write_index(out, hg38_alids, analyses, chunk_shape, dtype)
    canonical_variants = [
        CanonicalVariant(
            chromosome=chrom,
            position=int(pos_str),
            effect_allele=a1,
            other_allele=a2,
        )
        for alid in hg38_alids
        for chrom, pos_str, a1, a2 in [alid.split(":")]
    ]
    write_variant_axis(out, canonical_variants, {})

    # ------------------------------------------------------------------
    # Allocate output arrays (O(n_variants × n_analyses) peak memory)
    # ------------------------------------------------------------------
    z_mat = np.full((n_variants, n_analyses), np.nan, dtype=dtype)
    se_mat = np.full((n_variants, n_analyses), np.nan, dtype=dtype)

    # ------------------------------------------------------------------
    # Pass 2: fill zarr columns one analysis at a time
    # ------------------------------------------------------------------
    log.info("Pass 2: filling %d × %d association matrix", n_variants, n_analyses)
    for row in manifest_rows:
        col_idx = analysis_index[row.trait_id]
        for chrom, pos, ref, alt, z, se, _ in stream_vcf_associations(row.file_path):
            hg38_alid = hg19_lookup.get((chrom, pos, ref, alt))
            if hg38_alid is None:
                continue
            row_idx = variant_index.get(hg38_alid)
            if row_idx is None:
                continue
            z_mat[row_idx, col_idx] = z
            se_mat[row_idx, col_idx] = se
        log.debug("Pass 2: filled column %d (%s)", col_idx, row.trait_id)

    # ------------------------------------------------------------------
    # Write zarr + manifest + top-hit indexes
    # ------------------------------------------------------------------
    _write_zarr(out, z_mat, se_mat, chunk_shape, dtype)
    _write_manifest(
        out, store_id, release_id, n_variants, n_analyses, chain_file, chunk_shape, dtype
    )
    build_top_hit_indexes(out)

    return DenseBuildResult(output_path=out, n_variants=n_variants, n_analyses=n_analyses)


def _read_manifest(manifest_path: str | Path) -> list[_ManifestRow]:
    with open(manifest_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [
            _ManifestRow(
                trait_id=row["trait_id"],
                file_path=row["file_path"],
                trait_name=row.get("trait_name", row["trait_id"]),
                n=int(row.get("n", 0) or 0),
            )
            for row in reader
        ]


def _alid_sort_key(alid: str) -> tuple:
    chrom, pos_str, a1, a2 = alid.split(":")
    return (chromosome_sort_key(chrom), int(pos_str), a1, a2)


def _write_index(
    output_path: Path,
    hg38_alids: list[str],
    analyses: list[AnalysisMetadata],
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    with connect(output_path / "index.sqlite") as connection:
        initialise_schema(connection)
        set_metadata(connection, "schema_version", 1)
        set_metadata(connection, "n_variants", len(hg38_alids))
        set_metadata(connection, "n_analyses", len(analyses))
        set_metadata(
            connection,
            "dense",
            {"dtype": dtype, "chunk_shape": list(chunk_shape), "compressor": DEFAULT_COMPRESSOR},
        )
        connection.executemany(
            """
            INSERT INTO variants(
                variant_index, alid, chromosome, position, effect_allele, other_allele
            )
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [
                (i, alid, *_parse_alid(alid))
                for i, alid in enumerate(hg38_alids)
            ],
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
                    a.analysis_id,
                    a.phenotype_id,
                    a.phenotype_label,
                    a.analysis_label,
                    a.stored_effect_scale,
                )
                for i, a in enumerate(analyses)
            ],
        )
        connection.commit()


def _parse_alid(alid: str) -> tuple[str, int, str, str]:
    parts = alid.split(":")
    return parts[0], int(parts[1]), parts[2], parts[3]


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


def _write_manifest(
    output_path: Path,
    store_id: str,
    release_id: str,
    n_variants: int,
    n_analyses: int,
    chain_file: str | Path | None,
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
        reference_assembly="GRCh38",
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            "builder": "opengwasdb.v0.1_dense_vcf_two_pass",
            "chain_file": str(chain_file) if chain_file else "pyliftover_builtin_hg19_hg38",
            "n_variants": n_variants,
            "n_analyses": n_analyses,
            "dense": {
                "statistic_arrays": ["z", "se"],
                "dtype": dtype,
                "chunk_shape": list(chunk_shape),
                "compressor": DEFAULT_COMPRESSOR,
                "top_hit_thresholds": [5e-8, 5e-6, 5e-4],
            },
        },
    )
    (output_path / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
