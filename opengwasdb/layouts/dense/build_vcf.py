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


# Pass 2 fork-safe lookup + spill dir, set in the parent immediately before the
# process pool is created. Forked workers (fork start method — see _fork_pool())
# inherit these; because the lookup is a pair of numpy arrays (one contiguous C
# buffer with a single object refcount) rather than Python dicts, a worker reading
# it via searchsorted never touches per-element refcounts, so it does not
# refcount-COW ~n_variants dict pages per worker (issue 043 item 2).
#
# Why disk-spill and not return arrays: an earlier design returned each file's
# result over IPC. At genome-wide scale (~9.85M rows/file × thousands of files)
# that pipe traffic deadlocked the pool. Workers write a compact per-file .npz to
# _pass2_spill_dir and return only col_idx, so no large object crosses the pipe.
_pass2_keys_sorted: np.ndarray | None = None  # sorted 'S' byte keys
_pass2_rows_sorted: np.ndarray | None = None  # int32 row per key, same order
_pass2_spill_dir: Path | None = None

# Top hits are harvested inline during Pass 2 rather than by a post-hoc scan of
# the full matrix (which had to reload ~200 GB of float32 and compute a p-value
# for every finite cell). Each column emits only the cells clearing the loosest
# threshold's |z| cutoff; the band-write phase accumulates them for the index.
_TOP_HIT_Z_CRIT = z_critical(max(TOP_HIT_THRESHOLDS))


def _fork_pool(n_workers: int) -> ProcessPoolExecutor:
    """A ProcessPoolExecutor pinned to fork start — required for _pass2_worker to
    inherit the numpy lookup arrays without re-pickling them per task. Only
    correct on platforms with fork (Linux)."""
    fork_ctx = multiprocessing.get_context("fork")
    return ProcessPoolExecutor(max_workers=n_workers, mp_context=fork_ctx)


def _encode_variant_keys(
    chrom: object, pos: object, ref: object, alt: object
) -> np.ndarray:
    """Encode variant fields as ``chrom:pos:ref:alt`` byte-string keys, vectorised.

    Used identically for the parent's sorted key table and each worker's query
    keys, so the two encodings match exactly. ``pos`` int->bytes via numpy astype
    (``np.int64(10).astype('S') == b'10'``) avoids per-element Python string work.
    """
    chrom_b = np.asarray(chrom, dtype="S")
    pos_b = np.asarray(pos, dtype=np.int64).astype("S")
    ref_b = np.asarray(ref, dtype="S")
    alt_b = np.asarray(alt, dtype="S")
    key = np.char.add(chrom_b, b":")
    key = np.char.add(key, pos_b)
    key = np.char.add(key, b":")
    key = np.char.add(key, ref_b)
    key = np.char.add(key, b":")
    key = np.char.add(key, alt_b)
    return key


def _build_variant_key_index(
    hg19_lookup: dict[tuple[str, int, str, str], str],
    variant_index: dict[str, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Compose the two Pass 2 dicts into a fork-safe sorted numpy lookup.

    Returns ``(keys_sorted, rows_sorted)``: a sorted array of byte keys for every
    hg19 variant that maps to a stored row, and the matching int32 row indices.
    Workers binary-search this instead of chaining two Python dicts.
    """
    chroms: list[str] = []
    poss: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    rows: list[int] = []
    for (chrom, pos, ref, alt), alid in hg19_lookup.items():
        row = variant_index.get(alid)
        if row is None:
            continue
        chroms.append(chrom)
        poss.append(pos)
        refs.append(ref)
        alts.append(alt)
        rows.append(row)
    keys = _encode_variant_keys(chroms, poss, refs, alts)
    rows_arr = np.array(rows, dtype=np.int32)
    order = np.argsort(keys, kind="stable")
    return keys[order], rows_arr[order]


# Streamed in fixed-size batches so a worker never materialises a whole
# (genome-wide) VCF as Python lists at once — peak per-worker memory is one batch
# of columns plus the accumulated matched rows, not the entire file (issue 043).
_RESOLVE_BATCH = 250_000


def _match_batch(
    chroms: list[str],
    poss: list[int],
    refs: list[str],
    alts: list[str],
    zs: list[float],
    ses: list[float],
    keys_sorted: np.ndarray,
    rows_sorted: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resolve one batch of associations to (rows, z, se) for cells whose variant
    is present in the panel (exact key match). Order-preserving."""
    query = _encode_variant_keys(chroms, poss, refs, alts)
    idx = np.searchsorted(keys_sorted, query)
    idx_clip = np.minimum(idx, len(keys_sorted) - 1)
    matched = keys_sorted[idx_clip] == query
    rows = rows_sorted[idx_clip[matched]].astype(np.int64)
    z_arr = np.array(zs, dtype=np.float32)[matched]
    se_arr = np.array(ses, dtype=np.float32)[matched]
    return rows, z_arr, se_arr


def _resolve_column(
    file_path: str, keys_sorted: np.ndarray, rows_sorted: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stream one VCF, resolve each association's variant to a row via the sorted
    key lookup, and return deduped ``(rows int64, z f32, se f32)`` for the column.

    The stream is processed in ``_RESOLVE_BATCH``-sized batches (vectorised
    searchsorted per batch), so worker peak memory is bounded by one batch rather
    than the whole file. Batches are matched in stream order and concatenated, so
    the final last-wins dedup by row is identical to processing the whole file at
    once: when two source variants lift to the same row, the later stream
    occurrence wins, making the scattered matrix cell deterministic.
    """
    rows_parts: list[np.ndarray] = []
    z_parts: list[np.ndarray] = []
    se_parts: list[np.ndarray] = []
    chroms: list[str] = []
    poss: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    zs: list[float] = []
    ses: list[float] = []

    def _flush() -> None:
        if not zs:
            return
        r, z, se = _match_batch(chroms, poss, refs, alts, zs, ses, keys_sorted, rows_sorted)
        rows_parts.append(r)
        z_parts.append(z)
        se_parts.append(se)
        chroms.clear()
        poss.clear()
        refs.clear()
        alts.clear()
        zs.clear()
        ses.clear()

    for chrom, pos, ref, alt, z, se, _ in stream_vcf_associations(file_path):
        chroms.append(chrom)
        poss.append(pos)
        refs.append(ref)
        alts.append(alt)
        zs.append(z)
        ses.append(se)
        if len(zs) >= _RESOLVE_BATCH:
            _flush()
    _flush()

    if not rows_parts:
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float32)
        return empty_i, empty_f, empty_f

    rows = np.concatenate(rows_parts)
    z_arr = np.concatenate(z_parts)
    se_arr = np.concatenate(se_parts)

    if len(rows):
        # last-wins dedup by row: unique on the reversed rows returns the first
        # index in reversed order == the last occurrence in original order.
        _, first_in_rev = np.unique(rows[::-1], return_index=True)
        keep = np.sort(len(rows) - 1 - first_in_rev)
        rows, z_arr, se_arr = rows[keep], z_arr[keep], se_arr[keep]
    return rows, z_arr, se_arr


def _spill_column(
    spill_dir: Path, col_idx: int, rows: np.ndarray, z: np.ndarray, se: np.ndarray
) -> None:
    """Atomically spill one resolved column (+ its harvested top hits) to
    ``{spill_dir}/{col_idx}.npz`` (temp-then-rename; both names end in .npz because
    np.savez appends that suffix unless already present)."""
    hit = np.abs(z) >= _TOP_HIT_Z_CRIT
    final = spill_dir / f"{col_idx}.npz"
    tmp = spill_dir / f"{col_idx}.tmp.npz"
    np.savez(tmp, rows=rows, z=z, se=se, hit_rows=rows[hit], hit_z=z[hit], hit_se=se[hit])
    tmp.replace(final)


def _pass2_worker(task: tuple[int, str]) -> int:
    """Resolve one VCF column against the fork-inherited numpy lookup and spill it.
    Returns col_idx only — the compact result stays on disk, never in a pipe."""
    assert _pass2_keys_sorted is not None
    assert _pass2_rows_sorted is not None
    assert _pass2_spill_dir is not None
    col_idx, file_path = task
    rows, z, se = _resolve_column(file_path, _pass2_keys_sorted, _pass2_rows_sorted)
    _spill_column(_pass2_spill_dir, col_idx, rows, z, se)
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
    # Fork-safe Pass 2 lookup: compose the two dicts into sorted numpy arrays so
    # workers binary-search them instead of chaining Python dicts — no per-worker
    # refcount-COW of ~n_variants dict pages (issue 043 item 2).
    # ------------------------------------------------------------------
    keys_sorted, rows_sorted = _build_variant_key_index(hg19_lookup, variant_index)
    max_key_len = keys_sorted.dtype.itemsize if len(keys_sorted) else 0
    log.info(
        "Pass 2 lookup: %d variant keys, max key length %d bytes",
        len(keys_sorted), max_key_len,
    )
    del hg19_lookup, variant_index  # free the parent-side dicts before Pass 2

    # ------------------------------------------------------------------
    # Create the empty z/se zarr datasets (NaN fill). Pass 2 streams each analysis
    # column to disk; the band-write phase then fills chunk-column bands without
    # ever holding the full (n_variants × n_analyses) matrix in memory (issue 043).
    # ------------------------------------------------------------------
    effective_chunks = _create_dense_zarr(out, n_variants, n_analyses, chunk_shape, dtype)

    # ------------------------------------------------------------------
    # Pass 2 fill: resolve each analysis column and spill it to disk. No matrix is
    # resident here — the parent only waits for completion.
    # ------------------------------------------------------------------
    log.info(
        "Pass 2: resolving %d analyses × %d variants (n_workers=%d)",
        n_analyses, n_variants, n_workers,
    )
    pass2_start = time.monotonic()
    spill_dir = Path(tempfile.mkdtemp(prefix=f".{out.name}.pass2spill.", dir=out.parent))
    try:
        if n_workers <= 1:
            for i, row in enumerate(manifest_rows):
                col_idx = analysis_index[row.trait_id]
                rows, z, se = _resolve_column(row.file_path, keys_sorted, rows_sorted)
                _spill_column(spill_dir, col_idx, rows, z, se)
                _log_progress(
                    "Pass 2", i + 1, n_analyses, pass2_start, f"last: {row.trait_id}", every=25
                )
        else:
            global _pass2_keys_sorted, _pass2_rows_sorted, _pass2_spill_dir
            _pass2_keys_sorted = keys_sorted
            _pass2_rows_sorted = rows_sorted
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
                        _log_progress(
                            "Pass 2", i + 1, n_analyses, pass2_start,
                            f"last: {id_by_col[col_idx]}", every=25,
                        )
            finally:
                _pass2_keys_sorted = None
                _pass2_rows_sorted = None
                _pass2_spill_dir = None

        # --------------------------------------------------------------
        # Band-write phase: stream the retained column spills into the zarr in
        # chunk-column bands, harvesting top-hit candidates as we go. Peak memory
        # is one band (n_variants × chunk-analysis-width), not the full matrix.
        # --------------------------------------------------------------
        all_rows, all_cols, all_z, all_se = _write_dense_bands(
            out, spill_dir, n_variants, n_analyses, effective_chunks, dtype, pass2_start
        )
    finally:
        shutil.rmtree(spill_dir, ignore_errors=True)

    _write_manifest(
        out, store_id, release_id, n_variants, n_analyses, chain_file, chunk_shape, dtype
    )
    log.info("Writing top-hit index from %d harvested candidate cells", len(all_rows))
    write_top_hit_indexes(out, all_rows, all_cols, all_z, all_se)
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


def _create_dense_zarr(
    output_path: Path,
    n_variants: int,
    n_analyses: int,
    chunk_shape: tuple[int, int],
    dtype: str,
) -> tuple[int, int]:
    """Create the empty NaN-filled z/se datasets and return the effective chunks.

    Chunk shape is clipped to the array dimensions (ADR-0021) so zarr's declared
    chunks match what is physically stored. The arrays are created without data;
    ``_write_dense_bands`` fills them by chunk-column band afterwards.
    """
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    effective_chunks = (min(chunk_shape[0], n_variants), min(chunk_shape[1], n_analyses))
    root = zarr.open_group(str(output_path / "data.zarr"), mode="w")
    for name in ("z", "se"):
        root.create_dataset(
            name,
            shape=(n_variants, n_analyses),
            chunks=effective_chunks,
            compressor=compressor,
            dtype=dtype,
            fill_value=float("nan"),
        )
    root.attrs["layout"] = "dense"
    root.attrs["completion_state"] = "observed_only"
    root.attrs["compressor"] = DEFAULT_COMPRESSOR
    root.attrs["chunk_shape"] = list(effective_chunks)
    return effective_chunks


def _write_dense_bands(
    output_path: Path,
    spill_dir: Path,
    n_variants: int,
    n_analyses: int,
    effective_chunks: tuple[int, int],
    dtype: str,
    pass2_start: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Stream the retained per-column spills into the zarr in chunk-column bands.

    For each band of ``band_cols`` (= chunk analysis-width) analyses: reset a
    reusable ``(n_variants × band_cols)`` NaN buffer, scatter each column's spill
    into it (and accumulate that column's harvested top-hit candidates), then
    write the band to ``z``/``se`` as whole chunk-columns (no read-modify-write).
    Peak memory is one band, never the full matrix. Returns the concatenated
    top-hit candidate arrays ``(rows, cols, z, se)`` for the index build.
    """
    root = zarr.open_group(str(output_path / "data.zarr"), mode="a")
    z_arr = root["z"]
    se_arr = root["se"]
    band_cols = effective_chunks[1]
    z_band = np.empty((n_variants, band_cols), dtype=dtype)
    se_band = np.empty((n_variants, band_cols), dtype=dtype)

    hit_rows_parts: list[np.ndarray] = []
    hit_cols_parts: list[np.ndarray] = []
    hit_z_parts: list[np.ndarray] = []
    hit_se_parts: list[np.ndarray] = []

    log.info("Band-write: %d analyses to zarr in bands of %d", n_analyses, band_cols)
    for c0 in range(0, n_analyses, band_cols):
        c1 = min(c0 + band_cols, n_analyses)
        w = c1 - c0
        z_band[:, :w] = np.nan
        se_band[:, :w] = np.nan
        for c in range(c0, c1):
            local = c - c0
            with np.load(spill_dir / f"{c}.npz") as data:
                rows = data["rows"]
                z_band[rows, local] = data["z"]
                se_band[rows, local] = data["se"]
                hr = data["hit_rows"]
                if len(hr):
                    hit_rows_parts.append(hr)
                    hit_cols_parts.append(np.full(len(hr), c, dtype=np.int64))
                    hit_z_parts.append(data["hit_z"])
                    hit_se_parts.append(data["hit_se"])
            (spill_dir / f"{c}.npz").unlink()
        z_arr[:, c0:c1] = z_band[:, :w]
        se_arr[:, c0:c1] = se_band[:, :w]
        _log_progress(
            "Band-write", c1, n_analyses, pass2_start, f"cols {c0}:{c1}", every=band_cols
        )

    if hit_rows_parts:
        return (
            np.concatenate(hit_rows_parts),
            np.concatenate(hit_cols_parts),
            np.concatenate(hit_z_parts),
            np.concatenate(hit_se_parts),
        )
    return (
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.float32),
        np.empty(0, dtype=np.float32),
    )


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
