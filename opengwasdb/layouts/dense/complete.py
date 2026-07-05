"""Dense Reference Completion — enhancement pipeline.

Builds a Dense Reference-Completed Store Release from a Dense Observed-Only
Full Coverage source, per ADR 0022 (dense axis = source ∪ reference panel)
and ADR 0023 (LD-block process-pool parallelism with checkpointed resume).

Pipeline shape:
  Phase 1 (sequential): enumerate the genome-wide LD block set, build the
    union variant axis, seed z/se from the source, compute the per-Analysis
    n_missing_off_panel scalar.
  Phase 2 (parallel, n_workers processes over LD blocks): each worker opens
    the source store and LD panel itself, imputes serially per Analysis
    within its block, and writes its own checkpoint file.
  Phase 3 (sequential): merge all block results into the seeded z/se arrays,
    write the final zarr, completion_quality rows, top-hit indexes, and
    manifest.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import zarr
from numcodecs import Blosc
from threadpoolctl import threadpool_limits

from opengwasdb.completion.impute import impute_z_block, scalar_n_se
from opengwasdb.completion.ld_panel import (
    list_all_blocks,
    list_chromosomes,
    load_block,
    load_ld_eigenvectors,
)
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
    CanonicalVariant,
    VariantAxis,
    VariantNormalisationError,
    chromosome_sort_key,
    orient_to_canonical,
    parse_canonical_alid,
    write_variant_axis,
)

log = logging.getLogger(__name__)

_COMPRESSOR = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
_LD_PANEL_ID = "eur-hg38-gpm"
_COMPLETION_METHOD = "elastic_net_eigenvectors_v2_regioncap"

# Region-based imputation z-cap (QC): an imputed |z| may not exceed the largest
# observed |z| within +/- this many bp of it (same chromosome / LD block). Caps
# spurious large imputed z-scores from LD extrapolation. See pleiodb.
REGION_CAP_BP = 1_000_000


@dataclass(frozen=True)
class CompletionResult:
    output_path: Path
    n_variants: int
    n_analyses: int
    n_imputed: int
    n_missing_off_panel: int
    n_missing_imputation_failed: int


def _bare_alid(snp_id: str) -> str:
    return snp_id[3:] if snp_id.startswith("chr") else snp_id


def _canonical_panel_alid(snp_id: str) -> str | None:
    """Normalise an LD-panel SNP id to a canonical store ALID (``chr:pos:a1:a2``).

    Handles both the panel's ``[chr]CHR:POS_REF_ALT`` form (underscores) and an
    already-colon-delimited ``CHR:POS:A1:A2`` form, orienting alleles to the
    canonical A1 = min(ref, alt). Returns None for ids that don't parse.
    """
    s = _bare_alid(snp_id)
    parts = s.split(":")
    if len(parts) == 4:
        chrom, pos_s, ref, alt = parts
    elif len(parts) == 2 and parts[1].count("_") == 2:
        chrom, rest = parts
        pos_s, ref, alt = rest.split("_")
    else:
        return None
    try:
        return orient_to_canonical(chrom, int(pos_s), ref, alt).variant.alid
    except (VariantNormalisationError, ValueError):
        return None


def _snp_position(snp_id: str) -> int:
    """Extract the base-pair position from a panel SNP id (either format).
    Returns -1 if it cannot be parsed."""
    parts = _bare_alid(snp_id).split(":")
    field = parts[1] if len(parts) >= 2 else parts[0]
    try:
        return int(field.split("_")[0])
    except ValueError:
        return -1


def _sanitize_block_id(block_id: str) -> str:
    return block_id.replace("/", "__")


def _checkpoint_dir_for(dest_path: Path) -> Path:
    dest_path = Path(dest_path)
    return dest_path.parent / f".{dest_path.name}.checkpoint"


def _work_dir_for(dest_path: Path) -> Path:
    dest_path = Path(dest_path)
    return dest_path.parent / f".{dest_path.name}.tmp"


# ── Phase 2: per-block worker ───────────────────────────────────────────────


@dataclass(frozen=True)
class _BlockTask:
    tsv_path: Path
    source_path: Path
    min_cor: float
    thresh: float
    checkpoint_path: Path


@dataclass(frozen=True)
class BlockCompletionResult:
    block_id: str
    # (analysis_index, pearson_r | None, n_imputed, n_missing_imputation_failed)
    quality_rows: list[tuple[int, float | None, int, int]]
    # (alid, analysis_index, z, se)
    fills: list[tuple[str, int, float, float]]


def _write_checkpoint(path: Path, result: BlockCompletionResult) -> None:
    q_ai = np.array([r[0] for r in result.quality_rows], dtype=np.int32)
    q_pearson = np.array(
        [r[1] if r[1] is not None else np.nan for r in result.quality_rows], dtype=np.float64
    )
    q_nimp = np.array([r[2] for r in result.quality_rows], dtype=np.int32)
    q_nmiss = np.array([r[3] for r in result.quality_rows], dtype=np.int32)
    f_alid = np.array([r[0] for r in result.fills], dtype="U64")
    f_ai = np.array([r[1] for r in result.fills], dtype=np.int32)
    f_z = np.array([r[2] for r in result.fills], dtype=np.float64)
    f_se = np.array([r[3] for r in result.fills], dtype=np.float64)

    tmp_path = path.with_name(path.name + ".tmp")
    with open(tmp_path, "wb") as fh:
        np.savez(
            fh,
            block_id=np.array([result.block_id]),
            q_ai=q_ai, q_pearson=q_pearson, q_nimp=q_nimp, q_nmiss=q_nmiss,
            f_alid=f_alid, f_ai=f_ai, f_z=f_z, f_se=f_se,
        )
    os.replace(tmp_path, path)


# Kept alive for the worker's lifetime so the BLAS thread cap persists (a bare
# threadpool_limits() call would reset on garbage collection).
_worker_thread_limiter: Any = None


def _init_block_worker() -> None:
    """Cap each pool worker's BLAS (OpenBLAS/MKL) to one thread. numpy's linear
    algebra and sklearn otherwise spawn one thread per core *inside every worker*,
    so n_workers processes each with ~n_core threads massively oversubscribe the
    CPU (~1000 threads on 256 cores). One BLAS thread per worker gives clean
    process-level parallelism — scale with n_workers, not threads."""
    global _worker_thread_limiter
    _worker_thread_limiter = threadpool_limits(limits=1)


def _region_capped_z(
    zv: float,
    pos: int,
    obs_pos_sorted: np.ndarray,
    obs_absz_sorted: np.ndarray,
    window_bp: int = REGION_CAP_BP,
) -> float:
    """Clamp an imputed z-score so its magnitude does not exceed the largest
    observed |z| within +/- ``window_bp`` of ``pos`` (QC, per pleiodb). Returns
    ``zv`` unchanged when there is no observed variant in the window.

    ``obs_pos_sorted`` must be ascending, with ``obs_absz_sorted`` the matching
    |z| of the observed variants.
    """
    lo = int(np.searchsorted(obs_pos_sorted, pos - window_bp, side="left"))
    hi = int(np.searchsorted(obs_pos_sorted, pos + window_bp, side="right"))
    if hi <= lo:
        return zv
    z_region_max = float(obs_absz_sorted[lo:hi].max())
    if abs(zv) <= z_region_max:
        return zv
    return z_region_max if zv >= 0 else -z_region_max


def _run_block(task: _BlockTask) -> BlockCompletionResult | None:
    """Complete one LD block for every Analysis. Runs inside a worker process.

    Opens the (read-only, immutable) source store and LD panel itself, so no
    payload beyond a lightweight block descriptor needs to be pickled in.
    """
    block = load_block(task.tsv_path)
    if block is None:
        return None

    try:
        eigenvalues, eigenvectors = load_ld_eigenvectors(block, task.thresh)
    except Exception as exc:  # noqa: BLE001
        log.warning("Block %s: cannot load eigenvectors (%s) — skipping", block.block_id, exc)
        result = BlockCompletionResult(block_id=block.block_id, quality_rows=[], fills=[])
        _write_checkpoint(task.checkpoint_path, result)
        return result

    # Canonical ALID per block variant (None if it doesn't parse); the LD
    # eigenvectors are indexed by block position, so every block SNP keeps its
    # slot even when unmatched.
    canonical_alids = [_canonical_panel_alid(s) for s in block.snp_ids]

    src_axis = VariantAxis(task.source_path)
    try:
        src_root = zarr.open_group(str(Path(task.source_path) / "data.zarr"), mode="r")
        n_analyses = int(src_root["z"].shape[1])

        src_rows: list[int | None] = []
        for alid in canonical_alids:
            parsed = parse_canonical_alid(alid) if alid is not None else None
            rec = src_axis.by_alid(parsed) if parsed is not None else None
            src_rows.append(rec.variant_index if rec is not None else None)

        matched_local = [i for i, r in enumerate(src_rows) if r is not None]
        matched_src = [src_rows[i] for i in matched_local]

        z_obs = np.full((len(canonical_alids), n_analyses), np.nan, dtype=np.float64)
        se_obs = np.full((len(canonical_alids), n_analyses), np.nan, dtype=np.float64)
        if matched_local:
            z_obs[matched_local, :] = src_root["z"].oindex[matched_src, :].astype(np.float64)
            se_obs[matched_local, :] = src_root["se"].oindex[matched_src, :].astype(np.float64)
    finally:
        src_axis.close()

    eaf = block.eaf
    # Base-pair position per block variant (all same chromosome within a block),
    # for the region-based imputation z-cap below.
    positions = np.array([_snp_position(s) for s in block.snp_ids], dtype=np.int64)
    quality_rows: list[tuple[int, float | None, int, int]] = []
    fills: list[tuple[str, int, float, float]] = []

    for ai in range(n_analyses):
        z_dense = z_obs[:, ai]
        n_obs = int(np.isfinite(z_dense).sum())
        n_miss_block = int((~np.isfinite(z_dense)).sum())
        if n_obs < 2:
            quality_rows.append((ai, None, 0, n_miss_block))
            continue

        z_imp_arr, corr = impute_z_block(z_dense, eigenvectors, eigenvalues, min_cor=task.min_cor)
        if z_imp_arr is None:
            quality_rows.append((ai, float(corr) if np.isfinite(corr) else None, 0, n_miss_block))
            continue

        obs_mask = np.isfinite(z_dense)
        se_all = scalar_n_se(se_obs[obs_mask, ai], eaf[obs_mask], eaf)

        # Region z-cap (QC, per pleiodb): an imputed |z| must not exceed the
        # largest observed |z| among variants within +/- REGION_CAP_BP of it —
        # imputation can extrapolate spurious large z-scores. Observed positions
        # are sorted once so each missing variant's window is a searchsorted range.
        obs_idx = np.where(obs_mask)[0]
        o_order = np.argsort(positions[obs_idx])
        obs_pos_s = positions[obs_idx][o_order]
        obs_absz_s = np.abs(z_dense[obs_idx])[o_order]

        missing_mask = ~obs_mask
        n_filled = 0
        for i in np.where(missing_mask)[0]:
            alid = canonical_alids[i]
            zv, sev = z_imp_arr[i], se_all[i]
            if alid is None or not (np.isfinite(zv) and np.isfinite(sev)):
                continue
            zv = _region_capped_z(float(zv), int(positions[i]), obs_pos_s, obs_absz_s)
            fills.append((alid, ai, zv, float(sev)))
            n_filled += 1

        quality_rows.append((ai, float(corr), n_filled, n_miss_block - n_filled))

    result = BlockCompletionResult(block_id=block.block_id, quality_rows=quality_rows, fills=fills)
    _write_checkpoint(task.checkpoint_path, result)
    # Fills stay on disk (the checkpoint); the parent reads them back from the
    # checkpoints in Phase 3. Returning them here would push potentially millions
    # of tuples per block through the pool result queue and accumulate them in the
    # parent (issue 044 follow-up), so return an empty-fills marker.
    return BlockCompletionResult(block_id=block.block_id, quality_rows=[], fills=[])


# ── Public entry points ─────────────────────────────────────────────────────


def complete_dense_store(
    source_path: str | Path,
    dest_path: str | Path,
    ld_dir: str | Path,
    *,
    ancestry: str = "EUR",
    min_cor: float = 0.7,
    thresh: float = 0.9,
    release_id: str | None = None,
    ld_panel_id: str = _LD_PANEL_ID,
    n_workers: int = 1,
    overwrite: bool = False,
) -> CompletionResult:
    """Produce a Dense Reference-Completed Store Release from a Full Coverage
    Dense Observed-Only source.

    source_path: existing Dense Observed-Only, Full Coverage store.
    dest_path:   new store directory to create.
    ld_dir:      root of LD panel; blocks at ld_dir/{ancestry}/{chr}/{block}.*
    """
    dst = Path(dest_path)
    checkpoint_dir = _checkpoint_dir_for(dst)

    if dst.exists():
        if not overwrite:
            raise FileExistsError(f"Destination already exists: {dst}. Use overwrite=True.")
        shutil.rmtree(dst)

    if checkpoint_dir.exists():
        if not overwrite:
            raise FileExistsError(
                f"A checkpoint directory already exists at {checkpoint_dir}. "
                "Use resume_dense_completion() to continue it, or overwrite=True to discard it."
            )
        shutil.rmtree(checkpoint_dir)

    (checkpoint_dir / "blocks").mkdir(parents=True)
    build_params = {
        "source_path": str(Path(source_path).resolve()),
        "dest_path": str(dst.resolve()),
        "ld_dir": str(Path(ld_dir).resolve()),
        "ancestry": ancestry,
        "min_cor": min_cor,
        "thresh": thresh,
        "release_id": release_id,
        "ld_panel_id": ld_panel_id,
    }
    (checkpoint_dir / "build_params.json").write_text(
        json.dumps(build_params, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    result = _run_completion(
        Path(source_path), dst, Path(ld_dir),
        ancestry=ancestry, min_cor=min_cor, thresh=thresh,
        release_id=release_id, ld_panel_id=ld_panel_id,
        n_workers=n_workers, checkpoint_dir=checkpoint_dir,
    )
    shutil.rmtree(checkpoint_dir)
    return result


def resume_dense_completion(
    checkpoint_dir: str | Path,
    *,
    n_workers: int = 1,
) -> CompletionResult:
    """Resume an interrupted complete_dense_store() run.

    Takes only the checkpoint directory path — all other build parameters are
    loaded from the build_params.json written on the first run, so a resumed
    run can never silently apply a different parameter set than the one its
    existing per-block checkpoints were computed under.
    """
    checkpoint_dir = Path(checkpoint_dir)
    params = json.loads((checkpoint_dir / "build_params.json").read_text(encoding="utf-8"))

    result = _run_completion(
        Path(params["source_path"]), Path(params["dest_path"]), Path(params["ld_dir"]),
        ancestry=params["ancestry"], min_cor=params["min_cor"], thresh=params["thresh"],
        release_id=params["release_id"], ld_panel_id=params["ld_panel_id"],
        n_workers=n_workers, checkpoint_dir=checkpoint_dir,
    )
    shutil.rmtree(checkpoint_dir)
    return result


# ── Shared pipeline core ────────────────────────────────────────────────────


def _run_completion(
    source_path: Path,
    dest_path: Path,
    ld_dir: Path,
    *,
    ancestry: str,
    min_cor: float,
    thresh: float,
    release_id: str | None,
    ld_panel_id: str,
    n_workers: int,
    checkpoint_dir: Path,
) -> CompletionResult:
    src = Path(source_path)
    dst = Path(dest_path)
    work = _work_dir_for(dst)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    try:
        manifest = StoreManifest.load(src)
        if manifest.primary_layout is not PrimaryStorageLayout.DENSE:
            raise ValueError(
                f"source store is not Dense (primary_layout={manifest.primary_layout})"
            )
        if manifest.completion_state is not CompletionState.OBSERVED_ONLY:
            raise ValueError(
                f"source store is not Observed-Only (completion_state={manifest.completion_state})"
            )
        if manifest.association_coverage is not AssociationCoverage.FULL:
            raise ValueError(
                "Dense reference completion only supports Full Coverage sources "
                f"(association_coverage={manifest.association_coverage})"
            )
        print(f"Source store: {manifest.store_id} / {manifest.release_id}")

        # ── Phase 1: union variant axis + seeded z/se ───────────────────────
        src_variant_axis = VariantAxis(src)
        src_variants = src_variant_axis.all()
        src_variant_axis.close()
        src_alid_to_idx = {v.alid: v.variant_index for v in src_variants}

        with connect(src / "index.sqlite") as src_conn:
            src_analyses = [dict(r) for r in src_conn.execute(
                "SELECT * FROM analyses ORDER BY analysis_index"
            ).fetchall()]
        n_analyses = len(src_analyses)
        print(f"Source: {len(src_variants):,} variants, {n_analyses:,} analyses")

        print("Enumerating genome-wide LD blocks...")
        tsv_paths: list[Path] = []
        panel_alids: set[str] = set()
        for chrom in list_chromosomes(ld_dir, ancestry):
            for block in list_all_blocks(ld_dir, ancestry, chrom):
                tsv_paths.append(block.tsv_path)
                for snp_id in block.snp_ids:
                    ca = _canonical_panel_alid(snp_id)
                    if ca is not None:
                        panel_alids.add(ca)
        print(f"LD panel: {len(tsv_paths):,} blocks, {len(panel_alids):,} panel variants")

        new_canonical: list[CanonicalVariant] = []
        seen_new: set[str] = set()
        for alid in panel_alids:
            if alid in src_alid_to_idx:
                continue
            parts = alid.split(":")
            if len(parts) != 4:
                continue
            chrom, pos_str, a1, a2 = parts
            try:
                cv_result = orient_to_canonical(chrom, int(pos_str), a1, a2)
            except (VariantNormalisationError, ValueError):
                continue
            if cv_result.variant.alid in src_alid_to_idx or cv_result.variant.alid in seen_new:
                continue
            seen_new.add(cv_result.variant.alid)
            new_canonical.append(cv_result.variant)

        merged_variants: list[CanonicalVariant] = [
            CanonicalVariant(v.chromosome, v.position, v.effect_allele, v.other_allele)
            for v in src_variants
        ] + new_canonical
        merged_variants.sort(
            key=lambda v: (
                chromosome_sort_key(v.chromosome), v.position, v.effect_allele, v.other_allele
            )
        )
        new_alid_to_idx: dict[str, int] = {v.alid: i for i, v in enumerate(merged_variants)}
        n_variants = len(merged_variants)
        print(
            f"Union variant axis: {n_variants:,} variants "
            f"({len(new_canonical):,} new panel variants)"
        )

        on_panel = np.zeros(n_variants, dtype=bool)
        for alid in panel_alids:
            idx = new_alid_to_idx.get(alid)
            if idx is not None:
                on_panel[idx] = True

        rsid_by_alid = {v.alid: v.rsid for v in src_variants if v.rsid}
        print("Writing variants.tsv.gz...")
        write_variant_axis(work, merged_variants, rsid_by_alid)

        # Inverse map output_row -> source variant_index (-1 for panel-only rows).
        # z/se are seeded from the source band-by-band during the write (issue 044),
        # so the full source matrix is never loaded.
        src_root = zarr.open_group(str(src / "data.zarr"), mode="r")
        out_to_src = np.full(n_variants, -1, dtype=np.int64)
        for v in src_variants:
            out_to_src[new_alid_to_idx[v.alid]] = v.variant_index

        print("Writing index.sqlite...")
        with connect(work / "index.sqlite") as dst_db:
            initialise_schema(dst_db)
            dst_db.execute(
                "ALTER TABLE analyses ADD COLUMN n_missing_off_panel INTEGER NOT NULL DEFAULT 0"
            )
            dst_db.execute(
                """
                CREATE TABLE completion_quality (
                    analysis_index INTEGER NOT NULL,
                    block_id       TEXT    NOT NULL,
                    pearson_r      REAL,
                    n_imputed      INTEGER NOT NULL,
                    n_missing      INTEGER NOT NULL,
                    PRIMARY KEY (analysis_index, block_id)
                )
                """
            )
            set_metadata(dst_db, "schema_version", 2)
            set_metadata(dst_db, "n_variants", n_variants)
            set_metadata(dst_db, "n_analyses", n_analyses)
            dst_db.commit()
            # analyses rows are inserted after the band write, once
            # n_missing_off_panel is known (issue 044).

        # ── Phase 2: parallel LD-block completion ───────────────────────────
        print(f"Running reference completion across {len(tsv_paths):,} LD blocks "
              f"(n_workers={n_workers})...")
        blocks_dir = checkpoint_dir / "blocks"
        blocks_dir.mkdir(parents=True, exist_ok=True)

        pending: list[_BlockTask] = []
        n_existing = 0
        for tsv_path in tsv_paths:
            block_id = f"{tsv_path.parent.name}/{tsv_path.stem}"
            ckpt_path = blocks_dir / f"{_sanitize_block_id(block_id)}.npz"
            if ckpt_path.exists():
                n_existing += 1  # fills read from the checkpoint in Phase 3
            else:
                pending.append(_BlockTask(
                    tsv_path=tsv_path, source_path=src,
                    min_cor=min_cor, thresh=thresh, checkpoint_path=ckpt_path,
                ))

        if pending:
            print(f"  {n_existing:,} blocks already checkpointed, {len(pending):,} remaining")
        # Each block writes its own checkpoint; the parent keeps nothing per block.
        if n_workers <= 1:
            for i, task in enumerate(pending):
                _run_block(task)
                if (i + 1) % 200 == 0:
                    print(f"  {i + 1:,} / {len(pending):,} blocks")
        else:
            with ProcessPoolExecutor(
                max_workers=n_workers, initializer=_init_block_worker
            ) as pool:
                futures = [pool.submit(_run_block, task) for task in pending]
                for i, fut in enumerate(as_completed(futures)):
                    fut.result()  # propagate worker errors; result is on disk
                    if (i + 1) % 200 == 0:
                        print(f"  {i + 1:,} / {len(pending):,} blocks")

        # ── Phase 3: stream fills from checkpoints, band-write, finalise ────
        # Read each block's checkpoint as compact arrays and resolve its fill
        # ALIDs to union rows with one vectorised searchsorted — never Python
        # tuples/lists, so the parent's fill memory is packed arrays (~12 B/cell)
        # not ~O(n_imputed) tuples (issue 044 follow-up).
        print("Merging block results from checkpoints...")
        union_alids = np.fromiter(
            new_alid_to_idx.keys(), dtype="U64", count=len(new_alid_to_idx)
        )
        union_rows = np.fromiter(
            new_alid_to_idx.values(), dtype=np.int32, count=len(new_alid_to_idx)
        )
        o = np.argsort(union_alids)
        union_alids_s = union_alids[o]
        union_rows_s = union_rows[o]

        quality_rows: list[tuple[int, str, float | None, int, int]] = []
        fr_parts: list[np.ndarray] = []
        fai_parts: list[np.ndarray] = []
        fz_parts: list[np.ndarray] = []
        fse_parts: list[np.ndarray] = []
        for ckpt in sorted(blocks_dir.glob("*.npz")):
            with np.load(ckpt, allow_pickle=False) as d:
                bid = str(d["block_id"][0])
                for ai, p, ni, nm in zip(
                    d["q_ai"], d["q_pearson"], d["q_nimp"], d["q_nmiss"], strict=True
                ):
                    quality_rows.append(
                        (int(ai), bid, None if not np.isfinite(p) else float(p), int(ni), int(nm))
                    )
                f_alid = d["f_alid"]
                if len(f_alid):
                    idx = np.minimum(
                        np.searchsorted(union_alids_s, f_alid), len(union_alids_s) - 1
                    )
                    matched = union_alids_s[idx] == f_alid
                    if matched.any():
                        fr_parts.append(union_rows_s[idx[matched]])
                        fai_parts.append(d["f_ai"][matched].astype(np.int32))
                        fz_parts.append(d["f_z"][matched].astype(np.float32))
                        fse_parts.append(d["f_se"][matched].astype(np.float32))

        if fr_parts:
            fill_row = np.concatenate(fr_parts)
            fill_ai = np.concatenate(fai_parts)
            fill_z = np.concatenate(fz_parts)
            fill_se = np.concatenate(fse_parts)
        else:
            fill_row = np.empty(0, dtype=np.int32)
            fill_ai = np.empty(0, dtype=np.int32)
            fill_z = np.empty(0, dtype=np.float32)
            fill_se = np.empty(0, dtype=np.float32)
        order = np.argsort(fill_row, kind="stable")
        fill_row, fill_ai, fill_z, fill_se = (
            fill_row[order], fill_ai[order], fill_z[order], fill_se[order]
        )

        print("Writing data.zarr (band-streamed)...")
        effective_chunks = _create_completed_zarr(
            work, n_variants, n_analyses, on_panel, DEFAULT_CHUNK_SHAPE, DEFAULT_DTYPE
        )
        n_missing_off_panel, n_missing_imputation_failed, total_imputed = _write_completed_bands(
            work, src_root, out_to_src, on_panel,
            fill_row, fill_ai, fill_z, fill_se,
            effective_chunks, n_variants, n_analyses,
        )
        n_missing_off_panel_total = int(n_missing_off_panel.sum())
        print(
            f"Completion done: {total_imputed:,} imputed, "
            f"{n_missing_imputation_failed:,} imputation-failed, "
            f"{n_missing_off_panel_total:,} off-panel missing"
        )

        print("Writing analyses table...")
        with connect(work / "index.sqlite") as dst_db:
            dst_db.executemany(
                """
                INSERT INTO analyses (
                    analysis_index, analysis_id, phenotype_id, phenotype_label,
                    analysis_label, stored_effect_scale, n_missing_off_panel
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        int(row["analysis_index"]), row["analysis_id"], row["phenotype_id"],
                        row["phenotype_label"], row["analysis_label"], row["stored_effect_scale"],
                        int(n_missing_off_panel[i]),
                    )
                    for i, row in enumerate(src_analyses)
                ],
            )
            dst_db.commit()

        print(f"Writing {len(quality_rows):,} completion quality rows...")
        with connect(work / "index.sqlite") as dst_db:
            dst_db.executemany(
                "INSERT INTO completion_quality "
                "(analysis_index, block_id, pearson_r, n_imputed, n_missing) "
                "VALUES (?, ?, ?, ?, ?)",
                quality_rows,
            )
            dst_db.commit()

        print("Building top-hit indexes...")
        build_top_hit_indexes(work)

        new_release_id = release_id or f"{manifest.release_id}-completed"
        completed_manifest = StoreManifest(
            store_id=manifest.store_id,
            release_id=new_release_id,
            format_version=manifest.format_version,
            primary_layout=manifest.primary_layout,
            association_coverage=manifest.association_coverage,
            completion_state=CompletionState.REFERENCE_COMPLETED,
            reference_assembly=manifest.reference_assembly,
            created_at=datetime.now(UTC).isoformat(),
            provenance={
                **manifest.provenance,
                "source_release_id": manifest.release_id,
                "completion": {
                    "method": _COMPLETION_METHOD,
                    "ld_panel_id": ld_panel_id,
                    "ancestry": ancestry,
                    "min_cor": min_cor,
                    "pca_thresh": thresh,
                    "n_variants_total": n_variants,
                    "n_variants_new": len(new_canonical),
                    "n_imputed": total_imputed,
                    "n_missing_off_panel": n_missing_off_panel_total,
                    "n_missing_imputation_failed": n_missing_imputation_failed,
                },
            },
        )
        (work / "manifest.json").write_text(
            json.dumps(completed_manifest.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        if dst.exists():
            shutil.rmtree(dst)
        work.rename(dst)

        result = CompletionResult(
            output_path=dst,
            n_variants=n_variants,
            n_analyses=n_analyses,
            n_imputed=total_imputed,
            n_missing_off_panel=n_missing_off_panel_total,
            n_missing_imputation_failed=n_missing_imputation_failed,
        )
        print(
            f"Reference completion complete: {result.n_variants:,} variants, "
            f"{result.n_analyses:,} analyses ({result.n_imputed:,} imputed, "
            f"{result.n_missing_off_panel:,} off-panel missing, "
            f"{result.n_missing_imputation_failed:,} imputation-failed)"
        )
        return result
    except Exception:
        shutil.rmtree(work, ignore_errors=True)
        raise


# Row-band height for streaming the completed matrix — the seed/fill/write pass
# never holds the full (n_variants × n_analyses) matrices in RAM (issue 044).
_BAND_ROWS = 250_000


def _create_completed_zarr(
    output_path: Path,
    n_variants: int,
    n_analyses: int,
    on_panel: np.ndarray,
    chunk_shape: tuple[int, int],
    dtype: str,
) -> tuple[int, int]:
    """Create empty z/se (NaN), imputed (0), and the 1-D on_panel datasets. The
    matrices are filled by ``_write_completed_bands``; on_panel is small enough to
    write in one shot."""
    effective_chunks = (min(chunk_shape[0], n_variants), min(chunk_shape[1], n_analyses))
    root = zarr.open_group(str(output_path / "data.zarr"), mode="w")
    for name in ("z", "se"):
        root.create_dataset(
            name, shape=(n_variants, n_analyses), chunks=effective_chunks,
            compressor=_COMPRESSOR, dtype=dtype, fill_value=float("nan"),
        )
    root.create_dataset(
        "imputed", shape=(n_variants, n_analyses), chunks=effective_chunks,
        compressor=_COMPRESSOR, dtype="uint8", fill_value=0,
    )
    root.create_dataset(
        "on_panel", data=on_panel.astype(np.uint8),
        chunks=(effective_chunks[0],), compressor=_COMPRESSOR, dtype="uint8",
    )
    root.attrs["layout"] = "dense"
    root.attrs["completion_state"] = "reference_completed"
    root.attrs["compressor"] = DEFAULT_COMPRESSOR
    root.attrs["chunk_shape"] = list(effective_chunks)
    return effective_chunks


def _write_completed_bands(
    output_path: Path,
    src_root: Any,
    out_to_src: np.ndarray,
    on_panel: np.ndarray,
    fill_row: np.ndarray,
    fill_ai: np.ndarray,
    fill_z: np.ndarray,
    fill_se: np.ndarray,
    effective_chunks: tuple[int, int],
    n_variants: int,
    n_analyses: int,
) -> tuple[np.ndarray, int, int]:
    """Seed z/se from the source, apply the imputed fills, and write z/se/imputed
    one row-band at a time. ``z`` and ``se`` are written in two passes over a
    **single reused** float32 band buffer (plus a uint8 imputed band in the z-pass),
    so peak memory is ~one band rather than z + se + imputed held together. Both
    passes fill the same cells because source z/se missingness is consistent (a
    validated store invariant) — each pass reads only its own source array once, so
    there is no extra I/O. Returns ``(n_missing_off_panel[n_analyses],
    n_missing_imputation_failed, total_imputed)``. ``fill_*`` must be sorted by
    ``fill_row``.
    """
    root = zarr.open_group(str(output_path / "data.zarr"), mode="a")
    z_arr, se_arr, imp_arr = root["z"], root["se"], root["imputed"]
    src_z, src_se = src_root["z"], src_root["se"]
    band_rows = max(int(effective_chunks[0]), _BAND_ROWS)

    n_missing_off_panel = np.zeros(n_analyses, dtype=np.int64)
    n_missing_imputation_failed = 0
    total_imputed = 0
    band = np.empty((band_rows, n_analyses), dtype=np.float32)  # reused for z then se

    # Pass 1 — z + imputed mask + missingness counts.
    for r0 in range(0, n_variants, band_rows):
        r1 = min(r0 + band_rows, n_variants)
        br = r1 - r0
        zb = band[:br]
        zb[:] = np.nan
        imp_band = np.zeros((br, n_analyses), dtype=np.uint8)

        valid = np.where(out_to_src[r0:r1] >= 0)[0]
        if len(valid):
            srows = out_to_src[r0:r1][valid]
            zb[valid, :] = src_z.oindex[srows, :].astype(np.float32)

        off_local = np.where(on_panel[r0:r1] == 0)[0]
        if len(off_local):
            n_missing_off_panel += np.isnan(zb[off_local, :]).sum(axis=0).astype(np.int64)

        lo = int(np.searchsorted(fill_row, r0, side="left"))
        hi = int(np.searchsorted(fill_row, r1, side="left"))
        if hi > lo:
            lr = fill_row[lo:hi] - r0
            ai = fill_ai[lo:hi]
            fillable = ~np.isfinite(zb[lr, ai])
            if fillable.any():
                lrm, aim = lr[fillable], ai[fillable]
                zb[lrm, aim] = fill_z[lo:hi][fillable]
                imp_band[lrm, aim] = 1
                total_imputed += int(fillable.sum())

        on_local = np.where(on_panel[r0:r1] == 1)[0]
        if len(on_local):
            n_missing_imputation_failed += int(np.isnan(zb[on_local, :]).sum())

        z_arr[r0:r1] = zb
        imp_arr[r0:r1] = imp_band

    # Pass 2 — se (same cells filled, by the missingness-consistency invariant).
    for r0 in range(0, n_variants, band_rows):
        r1 = min(r0 + band_rows, n_variants)
        br = r1 - r0
        sb = band[:br]
        sb[:] = np.nan

        valid = np.where(out_to_src[r0:r1] >= 0)[0]
        if len(valid):
            srows = out_to_src[r0:r1][valid]
            sb[valid, :] = src_se.oindex[srows, :].astype(np.float32)

        lo = int(np.searchsorted(fill_row, r0, side="left"))
        hi = int(np.searchsorted(fill_row, r1, side="left"))
        if hi > lo:
            lr = fill_row[lo:hi] - r0
            ai = fill_ai[lo:hi]
            fillable = ~np.isfinite(sb[lr, ai])
            if fillable.any():
                sb[lr[fillable], ai[fillable]] = fill_se[lo:hi][fillable]

        se_arr[r0:r1] = sb

    return n_missing_off_panel, n_missing_imputation_failed, total_imputed
