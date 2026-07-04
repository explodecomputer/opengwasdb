"""Two-pass Dense Observed-Only writer from GWAS-VCF manifests with inline liftover."""

from __future__ import annotations

import csv
import json
import logging
import multiprocessing
import shutil
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    TOP_HIT_THRESHOLDS,
)
from opengwasdb.layouts.dense.top_hits import write_top_hit_indexes, z_critical
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


# Pass 2 read-only lookups + spill dir, set in the parent process immediately
# before the process pool is created. Forked workers (fork start method — see
# _fork_pool()) inherit these via copy-on-write, so the ~n_variants-entry dicts
# are never pickled or sent over IPC to any of the n_workers processes.
#
# Why disk-spill and not return arrays: an earlier design returned each file's
# result over IPC. At genome-wide scale (~9.85M rows/file × thousands of files)
# that pipe traffic deadlocked the pool. Workers now write a compact per-file
# .npz to _pass2_spill_dir and return only the tiny path string, so no large
# object ever crosses the process boundary.
_pass2_hg19_lookup: dict[tuple[str, int, str, str], str] | None = None
_pass2_variant_index: dict[str, int] | None = None
_pass2_spill_dir: Path | None = None

# Top hits are harvested inline during Pass 2 rather than by a post-hoc scan of
# the full matrix (which had to reload ~200 GB of float32 and compute a p-value
# for every finite cell). Each worker emits only the cells clearing the loosest
# threshold's |z| cutoff; the parent accumulates them and writes the index once.
_TOP_HIT_Z_CRIT = z_critical(max(TOP_HIT_THRESHOLDS))


def _fork_pool(n_workers: int) -> ProcessPoolExecutor:
    """A ProcessPoolExecutor pinned to fork start — required for _pass2_worker
    to see the read-only lookups without re-pickling them per task. Only correct
    on platforms with fork (Linux); would silently hand workers empty lookups
    under spawn/forkserver."""
    fork_ctx = multiprocessing.get_context("fork")
    return ProcessPoolExecutor(max_workers=n_workers, mp_context=fork_ctx)


def _pass2_worker(task: tuple[int, str]) -> int:
    """Stream one VCF's associations, resolve them against the (forked,
    read-only) hg19->hg38 lookup and variant index, and spill the resolved
    (row_index, z, se) triples to ``{spill_dir}/{col_idx}.npz``. Returns
    ``col_idx`` only — the compact result stays on disk, never in a pipe.

    Deduplicates by row_idx, keeping the last occurrence — matches the
    overwrite semantics of the sequential loop, and keeps the caller's
    ``z_mat[rows, col_idx] = ...`` assignment free of duplicate indices
    (numpy doesn't guarantee a winner order for those).
    """
    assert _pass2_hg19_lookup is not None
    assert _pass2_variant_index is not None
    assert _pass2_spill_dir is not None
    col_idx, file_path = task
    last_by_row: dict[int, tuple[float, float]] = {}
    for chrom, pos, ref, alt, z, se, _ in stream_vcf_associations(file_path):
        hg38_alid = _pass2_hg19_lookup.get((chrom, pos, ref, alt))
        if hg38_alid is None:
            continue
        row_idx = _pass2_variant_index.get(hg38_alid)
        if row_idx is None:
            continue
        last_by_row[row_idx] = (z, se)

    n = len(last_by_row)
    rows = np.fromiter(last_by_row.keys(), dtype=np.int64, count=n)
    zs = np.fromiter((v[0] for v in last_by_row.values()), dtype=np.float32, count=n)
    ses = np.fromiter((v[1] for v in last_by_row.values()), dtype=np.float32, count=n)

    # Harvest this analysis's top hits inline: only the cells clearing the
    # loosest threshold's |z| cutoff. Typically a tiny fraction of the column.
    hit = np.abs(zs) >= _TOP_HIT_Z_CRIT

    # Atomic spill: write to a temp path then rename, so a crashed worker never
    # leaves a half-written .npz that the parent would try to load. Both names
    # end in .npz because np.savez appends that suffix unless it is already
    # present — a .tmp suffix would be silently rewritten to .tmp.npz.
    final = _pass2_spill_dir / f"{col_idx}.npz"
    tmp = _pass2_spill_dir / f"{col_idx}.tmp.npz"
    np.savez(
        tmp,
        rows=rows, z=zs, se=ses,
        hit_rows=rows[hit], hit_z=zs[hit], hit_se=ses[hit],
    )
    tmp.replace(final)
    return col_idx


def _log_progress(
    label: str, completed: int, total: int, start_time: float, extra: str, every: int
) -> None:
    if completed % every != 0 and completed != total:
        return
    elapsed = time.monotonic() - start_time
    eta = (elapsed / completed) * (total - completed) if completed else 0.0
    log.info(
        "%s: %d/%d done (%s) — elapsed %s, ETA %s",
        label, completed, total, extra, _fmt_duration(elapsed), _fmt_duration(eta),
    )


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
    n_workers: int = 1,
) -> DenseBuildResult:
    """Build a Dense Observed-Only Store from a manifest of GWAS-VCF files.

    VCF files are assumed to be in GRCh37/hg19 coordinates.  All variant
    positions are lifted to GRCh38/hg38 inline; the output store uses hg38
    coordinates.

    Two-pass streaming: Pass 1 collects the union variant set and runs liftover
    once.  Pass 2 fills zarr columns one analysis at a time.  The full
    association list is never materialised in memory.

    Both passes process one file per analysis and are independent across
    files, so n_workers > 1 parallelises with a fork-based process pool —
    each analysis column is disjoint, so results merge back with no
    coordination beyond the final array assignment.

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
    n_workers:
        Process pool size for Pass 1 and Pass 2. 1 (default) runs both passes
        as a plain sequential loop. Requires the fork start method (Linux).
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
    # Pass 1 is intentionally serial. It streams each VCF's variants into one
    # growing union set; parallelising it would force each worker to ship its
    # whole variant set back over IPC, and since same-cohort VCFs share nearly
    # identical variant lists the union converges almost immediately — so the
    # parallel version pays a large IPC cost for no real speedup (and deadlocked
    # at genome-wide scale). The expensive, parallelised work is Pass 2.
    log.info(
        "Pass 1: collecting union variant set from %d VCFs (serial)",
        len(manifest_rows),
    )
    hg19_tuples: set[tuple[str, int, str, str]] = set()
    pass1_start = time.monotonic()
    n_rows = len(manifest_rows)
    for i, row in enumerate(manifest_rows):
        hg19_tuples.update(stream_vcf_variants(row.file_path))
        _log_progress(
            "Pass 1", i + 1, n_rows, pass1_start,
            f"{len(hg19_tuples)} unique variants so far", every=250,
        )
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
    # Allocate output arrays (O(n_variants × n_analyses)).
    #
    # Backed by MAP_SHARED memory-mapped files rather than anonymous RAM. The
    # Pass 2 process pool forks workers that inherit this address space; an
    # anonymous array allocated before the fork COW-doubles as the parent
    # scatters writes into it (children pin the untouched originals — issue 043).
    # MAP_SHARED pages are never copied on write, so the matrix stays a single
    # physical copy in page cache regardless of fork timing. Workers never touch
    # these arrays. The files are removed in the finally below.
    # ------------------------------------------------------------------
    mat_dir = Path(tempfile.mkdtemp(prefix=f".{out.name}.mat.", dir=out.parent))
    z_mat = np.memmap(mat_dir / "z.dat", dtype=dtype, mode="w+", shape=(n_variants, n_analyses))
    se_mat = np.memmap(mat_dir / "se.dat", dtype=dtype, mode="w+", shape=(n_variants, n_analyses))
    z_mat[:] = np.nan
    se_mat[:] = np.nan

    # ------------------------------------------------------------------
    # Pass 2: fill zarr columns one analysis at a time
    # ------------------------------------------------------------------
    log.info(
        "Pass 2: filling %d × %d association matrix (n_workers=%d)",
        n_variants, n_analyses, n_workers,
    )
    pass2_start = time.monotonic()
    # Top-hit candidates harvested inline (cells clearing the loosest |z| cutoff).
    hit_rows_parts: list[np.ndarray] = []
    hit_cols_parts: list[np.ndarray] = []
    hit_z_parts: list[np.ndarray] = []
    hit_se_parts: list[np.ndarray] = []

    def _collect_hits(col_idx: int, rr: np.ndarray, zz: np.ndarray, ss: np.ndarray) -> None:
        hit = np.abs(zz) >= _TOP_HIT_Z_CRIT
        if hit.any():
            hit_rows_parts.append(rr[hit])
            hit_cols_parts.append(np.full(int(hit.sum()), col_idx, dtype=np.int64))
            hit_z_parts.append(zz[hit])
            hit_se_parts.append(ss[hit])

    if n_workers <= 1:
        for i, row in enumerate(manifest_rows):
            col_idx = analysis_index[row.trait_id]
            # Same dedup-then-vectorise shape as _pass2_worker, so serial and
            # parallel builds produce byte-identical matrices and hit sets.
            last_by_row: dict[int, tuple[float, float]] = {}
            for chrom, pos, ref, alt, z, se, _ in stream_vcf_associations(row.file_path):
                hg38_alid = hg19_lookup.get((chrom, pos, ref, alt))
                if hg38_alid is None:
                    continue
                row_idx = variant_index.get(hg38_alid)
                if row_idx is None:
                    continue
                last_by_row[row_idx] = (z, se)
            if last_by_row:
                n = len(last_by_row)
                rr = np.fromiter(last_by_row.keys(), dtype=np.int64, count=n)
                zz = np.fromiter((v[0] for v in last_by_row.values()), dtype=np.float32, count=n)
                ss = np.fromiter((v[1] for v in last_by_row.values()), dtype=np.float32, count=n)
                z_mat[rr, col_idx] = zz
                se_mat[rr, col_idx] = ss
                _collect_hits(col_idx, rr, zz, ss)
            _log_progress(
                "Pass 2", i + 1, n_analyses, pass2_start, f"last: {row.trait_id}", every=25
            )
    else:
        global _pass2_hg19_lookup, _pass2_variant_index, _pass2_spill_dir
        _pass2_hg19_lookup = hg19_lookup
        _pass2_variant_index = variant_index
        spill_dir = Path(tempfile.mkdtemp(prefix=f".{out.name}.pass2spill.", dir=out.parent))
        _pass2_spill_dir = spill_dir
        id_by_col = {analysis_index[row.trait_id]: row.trait_id for row in manifest_rows}
        try:
            with _fork_pool(n_workers) as pool:
                tasks = [
                    (analysis_index[row.trait_id], row.file_path) for row in manifest_rows
                ]
                futures = [pool.submit(_pass2_worker, t) for t in tasks]
                for i, fut in enumerate(as_completed(futures)):
                    col_idx = fut.result()
                    spill = spill_dir / f"{col_idx}.npz"
                    with np.load(spill) as data:
                        z_mat[data["rows"], col_idx] = data["z"]
                        se_mat[data["rows"], col_idx] = data["se"]
                        hr = data["hit_rows"]
                        if len(hr):
                            hit_rows_parts.append(hr)
                            hit_cols_parts.append(np.full(len(hr), col_idx, dtype=np.int64))
                            hit_z_parts.append(data["hit_z"])
                            hit_se_parts.append(data["hit_se"])
                    spill.unlink()
                    _log_progress(
                        "Pass 2", i + 1, n_analyses, pass2_start,
                        f"last: {id_by_col[col_idx]}", every=25,
                    )
        finally:
            _pass2_hg19_lookup = None
            _pass2_variant_index = None
            _pass2_spill_dir = None
            shutil.rmtree(spill_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Write zarr + manifest + top-hit indexes (from inline-harvested hits)
    # ------------------------------------------------------------------
    if hit_rows_parts:
        all_rows = np.concatenate(hit_rows_parts)
        all_cols = np.concatenate(hit_cols_parts)
        all_z = np.concatenate(hit_z_parts)
        all_se = np.concatenate(hit_se_parts)
    else:
        all_rows = np.empty(0, dtype=np.int64)
        all_cols = np.empty(0, dtype=np.int64)
        all_z = np.empty(0, dtype=np.float32)
        all_se = np.empty(0, dtype=np.float32)

    z_mat.flush()
    se_mat.flush()
    _write_zarr(out, z_mat, se_mat, chunk_shape, dtype)
    _write_manifest(
        out, store_id, release_id, n_variants, n_analyses, chain_file, chunk_shape, dtype
    )
    log.info("Writing top-hit index from %d harvested candidate cells", len(all_rows))
    write_top_hit_indexes(out, all_rows, all_cols, all_z, all_se)

    # Release the memory-mapped matrices and delete their backing files.
    del z_mat, se_mat
    shutil.rmtree(mat_dir, ignore_errors=True)
    log.info("Build complete: %d variants × %d analyses", n_variants, n_analyses)

    return DenseBuildResult(output_path=out, n_variants=n_variants, n_analyses=n_analyses)


def _fmt_duration(seconds: float) -> str:
    total = int(seconds)
    h, rem = divmod(total, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


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
