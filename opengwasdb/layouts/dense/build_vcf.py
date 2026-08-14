"""Two-pass Dense Observed-Only writer from GWAS-VCF manifests with inline liftover.

Association streaming and the union-variant pass go through a ``SourceReader``
resolved from each row's ``source_reader_capability`` (issue #20) rather than
importing ``opengwasdb.build.vcf_source`` directly -- GWAS-VCF remains the only
implementation, so nothing about the build changes, but no bcftools-specific
assumption is left in this module.
"""

from __future__ import annotations

import csv
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
from numcodecs import Blosc

from opengwasdb.build.liftover import LiftoverFailureError, build_liftover_lookup
from opengwasdb.index import initialise_schema, set_metadata
from opengwasdb.layouts.dense.build import (
    AnalysisMetadata,
    DenseBuildResult,
    add_hit_counts,
    write_analyses_tsv,
)
from opengwasdb.layouts.dense.constants import (
    DEFAULT_CHUNK_SHAPE,
    DEFAULT_COMPRESSOR,
    DEFAULT_DTYPE,
    TOP_HIT_THRESHOLDS,
)
from opengwasdb.layouts.dense.top_hits import write_top_hit_indexes, z_critical
from opengwasdb.model.analyses import ANCESTRY_PROP_PREFIX
from opengwasdb.model.enums import (
    AncestryAssignmentMethod,
    AssociationCoverage,
    CompletionState,
    OriginalSdMethod,
    PrimaryStorageLayout,
    StoredEffectScale,
)
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.readers.gwas_vcf import GWAS_VCF_CAPABILITY
from opengwasdb.readers.registry import resolve_reader
from opengwasdb.store.open import CURRENT_FORMAT_VERSION, OpenGWASDBStore, StagedRelease
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
    stored_effect_scale: str
    se_divisor: float  # divides original_se to standardise to SD units (issue #18); 1.0 = no-op
    source_reader_capability: str  # resolves to a SourceReader (issue #20); GWAS_VCF_CAPABILITY
    # when the manifest omits the column -- the only format this builder has ever supported.
    original_sd_method: str  # raw manifest value, carried into analyses.tsv (issue #22)
    original_sd: str  # raw manifest value ("" when the sd_method tier carries no magnitude)
    assigned_ancestry: str  # optional manifest column (issue #22); "" when the manifest omits it
    # or carries no assignment for this row -- e.g. a Catalogue subset's kept rows (already
    # filtered to one target ancestry) stamp this in verbatim; a bare manifest has no column.
    ancestry_prop: dict[str, str]  # {population: proportion}, from any ancestry_prop_<pop>
    # column(s) the manifest carries -- a Catalogue subset's kept rows already have these
    # (the Catalogue writes one per reference super-population); empty for a bare manifest.


def _manifest_row_to_analysis_metadata(row: _ManifestRow) -> AnalysisMetadata:
    """A manifest row's Analytical Metadata (issue #22), shared by the dense
    and hybrid manifest builders.

    ``ancestry_assignment_method`` is derived, not guessed: every
    ``assigned_ancestry`` this builder ever sees comes from a Catalogue's
    AF-based NNLS mixture fit (`opengwasdb.ancestry.mixture.assign_ancestry`),
    so naming that method for a row that carries an assignment is reporting
    fact, not inference. A row with no ``assigned_ancestry`` gets no method
    either -- ancestry was never attempted for it, which is a different state
    from ``AncestryAssignmentMethod.UNASSIGNED`` (attempted, gates failed).
    """
    return AnalysisMetadata(
        analysis_id=row.trait_id,
        phenotype_id=row.trait_id,
        phenotype_label=row.trait_name,
        analysis_label=row.trait_id,
        stored_effect_scale=row.stored_effect_scale,
        assigned_ancestry=row.assigned_ancestry,
        ancestry_assignment_method=(
            AncestryAssignmentMethod.AF_ASSIGNED.value if row.assigned_ancestry else ""
        ),
        ancestry_prop=row.ancestry_prop,
        sample_size=str(row.n) if row.n else "",
        original_sd_method=row.original_sd_method,
        original_sd=row.original_sd,
    )


# original_sd_method tiers that carry an actual phenotype-SD magnitude to divide
# by (ADR-0029 methods 2-5). declared_standardised implies sd=1 (no magnitude
# recorded); binary_trait is not SD-scale at all -- neither rescales.
_SD_RESCALE_METHODS = frozenset(
    {
        OriginalSdMethod.SOURCE_PROVIDED,
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF,
        OriginalSdMethod.ESTIMATED_FROM_REFERENCE_MAF,
        OriginalSdMethod.ESTIMATED_FROM_BETA_DISTRIBUTION,
    }
)


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


def _apply_se_divisor(se: np.ndarray, se_divisor: float) -> np.ndarray:
    """Continuous-trait phenotype-SD standardisation (issue #18):
    ``stored_se = original_se / sd``. Shared by the dense and hybrid builders,
    both of which resolve an association stream to ``(index, z, se)`` and need
    the same no-op-when-1.0 divide applied to ``se`` before spilling."""
    if se_divisor == 1.0:
        return se
    return se / np.float32(se_divisor)


def _resolve_column(
    file_path: str,
    keys_sorted: np.ndarray,
    rows_sorted: np.ndarray,
    se_divisor: float = 1.0,
    *,
    capability: str = GWAS_VCF_CAPABILITY,
    stored_effect_scale: str = StoredEffectScale.SD.value,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stream one source file, resolve each association's variant to a row via
    the sorted key lookup, and return deduped ``(rows int64, z f32, se f32)``
    for the column.

    The stream is processed in ``_RESOLVE_BATCH``-sized batches (vectorised
    searchsorted per batch), so worker peak memory is bounded by one batch rather
    than the whole file. Batches are matched in stream order and concatenated, so
    the final last-wins dedup by row is identical to processing the whole file at
    once: when two source variants lift to the same row, the later stream
    occurrence wins, making the scattered matrix cell deterministic.

    ``capability`` resolves a ``SourceReader`` (issue #20) rather than this
    module streaming a VCF itself -- ``stored_effect_scale`` is required to
    construct one but unused past construction here (it is a Pass-2 concern
    only for readers that attach it to each yielded association).

    ``se_divisor`` divides the returned ``se`` (continuous-trait phenotype-SD
    standardisation, issue #18: ``stored_se = original_se / sd``). ``z`` is left
    untouched -- ``z = beta/se`` is invariant to dividing both by the same
    constant, so only ``se`` needs rescaling. Defaults to 1.0 (no-op) for
    binary-trait analyses and callers that pre-date issue #18.
    """
    reader = resolve_reader(capability, file_path, StoredEffectScale(stored_effect_scale))
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

    for assoc in reader.stream_associations():
        chroms.append(assoc.chromosome)
        poss.append(assoc.position)
        refs.append(assoc.ref)
        alts.append(assoc.alt)
        zs.append(assoc.z)
        ses.append(assoc.se)
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
    se_arr = _apply_se_divisor(se_arr, se_divisor)
    return rows, z_arr, se_arr


def _spill_column(
    spill_dir: Path, col_idx: int, rows: np.ndarray, z: np.ndarray, se: np.ndarray
) -> None:
    """Atomically spill one resolved column to ``{spill_dir}/{col_idx}.npz``
    (temp-then-rename; both names end in .npz because np.savez appends that suffix
    unless already present).

    Top hits are NOT harvested here — they are harvested during the band-write
    phase from the *stored* float16 values, so the index matches exactly what a
    query reads back from the ``z`` array (issue 046)."""
    final = spill_dir / f"{col_idx}.npz"
    tmp = spill_dir / f"{col_idx}.tmp.npz"
    np.savez(tmp, rows=rows, z=z, se=se)
    tmp.replace(final)


def _pass2_worker(task: tuple[int, str, float, str, str]) -> int:
    """Resolve one column against the fork-inherited numpy lookup and spill it.
    Returns col_idx only — the compact result stays on disk, never in a pipe."""
    assert _pass2_keys_sorted is not None
    assert _pass2_rows_sorted is not None
    assert _pass2_spill_dir is not None
    col_idx, file_path, se_divisor, capability, stored_effect_scale = task
    rows, z, se = _resolve_column(
        file_path,
        _pass2_keys_sorted,
        _pass2_rows_sorted,
        se_divisor,
        capability=capability,
        stored_effect_scale=stored_effect_scale,
    )
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
        TSV with columns ``trait_id``, ``file_path``, ``trait_name``, ``n``
        (also the Analysis Catalogue's ``BUILD_COLUMNS`` --
        `opengwasdb.ancestry.catalogue` -- so a Catalogue-annotated file
        remains readable here), plus a required ``stored_effect_scale``
        (issue #17): Analytical Metadata the manifest supplies, never
        inferred from the VCF header. A missing column or out-of-vocabulary
        value fails the build before any I/O.
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
    with OpenGWASDBStore.staging(out, overwrite=overwrite) as staged:
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
            reader = resolve_reader(
                row.source_reader_capability, row.file_path, StoredEffectScale(row.stored_effect_scale)
            )
            hg19_tuples.update(reader.stream_variants())
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
        # Analytical Metadata: stored_effect_scale comes from the manifest, not
        # the VCF header (issue #17 -- the ieu-a-7 fix: the source header is not
        # authoritative for effect scale).
        # ------------------------------------------------------------------
        analyses: list[AnalysisMetadata] = [
            _manifest_row_to_analysis_metadata(row) for row in manifest_rows
        ]

        # ------------------------------------------------------------------
        # Write SQLite index + analyses.tsv + tabix variant axis
        # ------------------------------------------------------------------
        _write_index(staged, hg38_alids, analyses, chunk_shape, dtype)
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
        # Provenance: record the source (hg19) canonical ALID each row was lifted
        # from. A single hg38 ALID can be the liftover target of several hg19
        # variants (a collision); those rows are ambiguous, so leave them blank.
        hg38_to_source: dict[str, str | None] = {}
        for (chrom, pos, ref, alt), hg38_alid in hg19_lookup.items():
            a1, a2 = sorted((ref, alt))
            origin = f"{chrom}:{pos}:{a1}:{a2}"
            if hg38_alid in hg38_to_source:
                if hg38_to_source[hg38_alid] != origin:
                    hg38_to_source[hg38_alid] = None  # collision → ambiguous
            else:
                hg38_to_source[hg38_alid] = origin
        source_alids = [hg38_to_source.get(alid) for alid in hg38_alids]
        write_variant_axis(staged.path, canonical_variants, {}, source_alids)

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
        effective_chunks = _create_dense_zarr(staged, n_variants, n_analyses, chunk_shape, dtype)

        # ------------------------------------------------------------------
        # Pass 2 fill: resolve each analysis column and spill it to disk. No matrix is
        # resident here — the parent only waits for completion.
        # ------------------------------------------------------------------
        log.info(
            "Pass 2: resolving %d analyses × %d variants (n_workers=%d)",
            n_analyses, n_variants, n_workers,
        )
        pass2_start = time.monotonic()
        spill_dir = Path(
            tempfile.mkdtemp(prefix=f".{out.name}.pass2spill.", dir=staged.path.parent)
        )
        try:
            if n_workers <= 1:
                for i, row in enumerate(manifest_rows):
                    col_idx = analysis_index[row.trait_id]
                    rows, z, se = _resolve_column(
                        row.file_path,
                        keys_sorted,
                        rows_sorted,
                        row.se_divisor,
                        capability=row.source_reader_capability,
                        stored_effect_scale=row.stored_effect_scale,
                    )
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
                            (
                                analysis_index[row.trait_id],
                                row.file_path,
                                row.se_divisor,
                                row.source_reader_capability,
                                row.stored_effect_scale,
                            )
                            for row in manifest_rows
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
                staged, spill_dir, n_variants, n_analyses, effective_chunks, dtype, pass2_start
            )
        finally:
            shutil.rmtree(spill_dir, ignore_errors=True)

        _write_manifest(
            staged, store_id, release_id, n_variants, n_analyses, chain_file, chunk_shape, dtype
        )
        log.info("Writing top-hit index from %d harvested candidate cells", len(all_rows))
        write_top_hit_indexes(staged.path, all_rows, all_cols, all_z, all_se)
        write_analyses_tsv(staged.path, add_hit_counts(staged.path, analyses))
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
    """Read the build manifest: the existing ``trait_id``/``file_path``/
    ``trait_name``/``n`` columns -- also the Analysis Catalogue's
    ``BUILD_COLUMNS`` (`opengwasdb.ancestry.catalogue`), so a
    Catalogue-annotated file remains readable here for those four columns --
    plus required ``stored_effect_scale`` (issue #17) and ``original_sd_method``
    (issue #18) columns, and an ``original_sd`` column required only for the
    ``original_sd_method`` tiers that carry an actual SD magnitude.

    Neither column is part of the Catalogue/ancestry manifest shape: ancestry
    assignment runs on allele frequencies alone and may run before a study's
    effect scale or phenotype SD is even resolved, so these stay independent
    build inputs rather than part of one combined schema. A manifest missing a
    required column, carrying an out-of-vocabulary value, declaring
    ``original_sd_method=unavailable``, omitting ``original_sd`` for a tier
    that needs it, or supplying a stray ``original_sd`` for a tier that carries
    no SD magnitude (``declared_standardised``, ``binary_trait``), fails the
    build loudly rather than falling back to the old VCF-header inference or a
    silently assumed ``sd=1`` (issue #18 AC3).

    An optional ``source_reader_capability`` column (issue #20) selects the
    ``SourceReader`` each row is built through; a manifest that omits it (every
    manifest before this change) defaults every row to ``GWAS_VCF_CAPABILITY``,
    the only format this builder has ever supported.

    An optional ``assigned_ancestry`` column (issue #22) carries a row's
    Assigned Ancestry straight into the built store's ``analyses.tsv`` --
    a Catalogue subset's kept rows already have this column (the Catalogue
    is a superset of the build manifest), so ``opengwasdb.ancestry.subset``
    no longer needs a separate post-build sidecar write to record it. Blank
    or absent for a manifest with no ancestry annotation.
    """
    with open(manifest_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    for column in ("stored_effect_scale", "original_sd_method"):
        if column not in fieldnames:
            raise ValueError(f"manifest {manifest_path} is missing required column: {column!r}")
    result = []
    for row in rows:
        trait_id = row["trait_id"]
        scale = row["stored_effect_scale"]
        try:
            StoredEffectScale(scale)
        except ValueError as exc:
            raise ValueError(
                f"manifest {manifest_path}: analysis {trait_id!r} has invalid "
                f"stored_effect_scale {scale!r}"
            ) from exc
        sd_method_raw = row["original_sd_method"]
        try:
            sd_method = OriginalSdMethod(sd_method_raw)
        except ValueError as exc:
            raise ValueError(
                f"manifest {manifest_path}: analysis {trait_id!r} has invalid "
                f"original_sd_method {sd_method_raw!r}"
            ) from exc
        if sd_method is OriginalSdMethod.UNAVAILABLE:
            raise ValueError(
                f"manifest {manifest_path}: analysis {trait_id!r} has "
                "original_sd_method='unavailable' -- its phenotype SD could not be "
                "established upstream, so the build cannot standardise its effects (issue #18)"
            )
        se_divisor = 1.0
        original_sd_raw = row.get("original_sd", "")
        if sd_method in _SD_RESCALE_METHODS:
            try:
                se_divisor = float(original_sd_raw)
            except ValueError as exc:
                raise ValueError(
                    f"manifest {manifest_path}: analysis {trait_id!r} has "
                    f"original_sd_method={sd_method.value!r} but original_sd "
                    f"{original_sd_raw!r} is not a valid number"
                ) from exc
            if not se_divisor > 0:
                raise ValueError(
                    f"manifest {manifest_path}: analysis {trait_id!r} has "
                    f"non-positive original_sd {original_sd_raw!r}"
                )
        elif original_sd_raw:
            raise ValueError(
                f"manifest {manifest_path}: analysis {trait_id!r} has "
                f"original_sd_method={sd_method.value!r}, which carries no SD "
                f"magnitude, but original_sd={original_sd_raw!r} was supplied"
            )
        result.append(
            _ManifestRow(
                trait_id=trait_id,
                file_path=row["file_path"],
                trait_name=row.get("trait_name", trait_id),
                n=int(row.get("n", 0) or 0),
                stored_effect_scale=scale,
                se_divisor=se_divisor,
                source_reader_capability=row.get("source_reader_capability") or GWAS_VCF_CAPABILITY,
                original_sd_method=sd_method_raw,
                original_sd=original_sd_raw,
                assigned_ancestry=row.get("assigned_ancestry") or "",
                ancestry_prop={
                    key[len(ANCESTRY_PROP_PREFIX):]: value
                    for key, value in row.items()
                    if key.startswith(ANCESTRY_PROP_PREFIX) and value
                },
            )
        )
    return result


def _alid_sort_key(alid: str) -> tuple:
    chrom, pos_str, a1, a2 = alid.split(":")
    return (chromosome_sort_key(chrom), int(pos_str), a1, a2)


def _write_index(
    staged: StagedRelease,
    hg38_alids: list[str],
    analyses: list[AnalysisMetadata],
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    with staged.index_connection() as connection:
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
        connection.commit()


def _parse_alid(alid: str) -> tuple[str, int, str, str]:
    parts = alid.split(":")
    return parts[0], int(parts[1]), parts[2], parts[3]


def _create_dense_zarr(
    staged: StagedRelease,
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
    root = staged.arrays(mode="w")
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
    staged: StagedRelease,
    spill_dir: Path,
    n_variants: int,
    n_analyses: int,
    effective_chunks: tuple[int, int],
    dtype: str,
    pass2_start: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Stream the retained per-column spills into the zarr in chunk-column bands.

    ``z`` and ``se`` are written in two separate passes over a **single reused**
    ``(n_variants × band_cols)`` buffer, so only one band is ever resident —
    halving the band-write peak versus holding z and se together. The top-hit
    harvest runs in the z-pass and reads each hit's ``se`` straight from the spill
    (rounded to the stored dtype), so the se-band is never needed for it. Spills
    are retained until the se-pass consumes them. Returns the concatenated top-hit
    candidate arrays ``(rows, cols, z, se)`` for the index build.
    """
    root = staged.arrays(mode="a")
    z_arr = root["z"]
    se_arr = root["se"]
    band_cols = effective_chunks[1]
    band = np.empty((n_variants, band_cols), dtype=dtype)  # reused for z then se

    hit_rows_parts: list[np.ndarray] = []
    hit_cols_parts: list[np.ndarray] = []
    hit_z_parts: list[np.ndarray] = []
    hit_se_parts: list[np.ndarray] = []

    # Pass 1 — z. Fill, write z, and harvest top hits (se pulled from the spill,
    # rounded to the stored dtype so the index matches what queries read from z).
    log.info("Band-write z: %d analyses in bands of %d", n_analyses, band_cols)
    for c0 in range(0, n_analyses, band_cols):
        c1 = min(c0 + band_cols, n_analyses)
        w = c1 - c0
        band[:, :w] = np.nan
        for c in range(c0, c1):
            local = c - c0
            with np.load(spill_dir / f"{c}.npz") as data:
                rows = data["rows"]
                band[rows, local] = data["z"]
                zc = band[rows, local]  # stored float16 z
                hit = np.abs(zc.astype(np.float32)) >= _TOP_HIT_Z_CRIT
                if np.any(hit):
                    hit_rows_parts.append(rows[hit])
                    hit_cols_parts.append(np.full(int(np.count_nonzero(hit)), c, dtype=np.int64))
                    hit_z_parts.append(zc[hit].astype(np.float32))
                    hit_se_parts.append(data["se"][hit].astype(dtype).astype(np.float32))
        z_arr[:, c0:c1] = band[:, :w]
        _log_progress(
            "Band-write z", c1, n_analyses, pass2_start, f"cols {c0}:{c1}", every=band_cols
        )

    # Pass 2 — se. Reuse the buffer; fill, write se, and delete each spill.
    for c0 in range(0, n_analyses, band_cols):
        c1 = min(c0 + band_cols, n_analyses)
        w = c1 - c0
        band[:, :w] = np.nan
        for c in range(c0, c1):
            local = c - c0
            with np.load(spill_dir / f"{c}.npz") as data:
                band[data["rows"], local] = data["se"]
            (spill_dir / f"{c}.npz").unlink()
        se_arr[:, c0:c1] = band[:, :w]
        _log_progress(
            "Band-write se", c1, n_analyses, pass2_start, f"cols {c0}:{c1}", every=band_cols
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
    staged: StagedRelease,
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
        format_version=CURRENT_FORMAT_VERSION,
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
    staged.write_manifest(manifest)
