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

import numpy as np
import zarr
from numcodecs import Blosc

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
_COMPLETION_METHOD = "elastic_net_eigenvectors_v1"


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


def _load_checkpoint(path: Path) -> BlockCompletionResult:
    data = np.load(path, allow_pickle=False)
    block_id = str(data["block_id"][0])
    quality_rows = [
        (int(ai), None if not np.isfinite(p) else float(p), int(ni), int(nm))
        for ai, p, ni, nm in zip(
            data["q_ai"], data["q_pearson"], data["q_nimp"], data["q_nmiss"], strict=True
        )
    ]
    fills = [
        (str(alid), int(ai), float(z), float(se))
        for alid, ai, z, se in zip(
            data["f_alid"], data["f_ai"], data["f_z"], data["f_se"], strict=True
        )
    ]
    return BlockCompletionResult(block_id=block_id, quality_rows=quality_rows, fills=fills)


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

    bare_alids = [_bare_alid(s) for s in block.snp_ids]

    src_axis = VariantAxis(task.source_path)
    try:
        src_root = zarr.open_group(str(Path(task.source_path) / "data.zarr"), mode="r")
        n_analyses = int(src_root["z"].shape[1])

        src_rows: list[int | None] = []
        for alid in bare_alids:
            parsed = parse_canonical_alid(alid)
            rec = src_axis.by_alid(parsed) if parsed is not None else None
            src_rows.append(rec.variant_index if rec is not None else None)

        matched_local = [i for i, r in enumerate(src_rows) if r is not None]
        matched_src = [src_rows[i] for i in matched_local]

        z_obs = np.full((len(bare_alids), n_analyses), np.nan, dtype=np.float64)
        se_obs = np.full((len(bare_alids), n_analyses), np.nan, dtype=np.float64)
        if matched_local:
            z_obs[matched_local, :] = src_root["z"].oindex[matched_src, :].astype(np.float64)
            se_obs[matched_local, :] = src_root["se"].oindex[matched_src, :].astype(np.float64)
    finally:
        src_axis.close()

    eaf = block.eaf
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

        missing_mask = ~obs_mask
        n_filled = 0
        for i in np.where(missing_mask)[0]:
            zv, sev = z_imp_arr[i], se_all[i]
            if np.isfinite(zv) and np.isfinite(sev):
                fills.append((bare_alids[i], ai, float(zv), float(sev)))
                n_filled += 1

        quality_rows.append((ai, float(corr), n_filled, n_miss_block - n_filled))

    result = BlockCompletionResult(block_id=block.block_id, quality_rows=quality_rows, fills=fills)
    _write_checkpoint(task.checkpoint_path, result)
    return result


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
                    panel_alids.add(_bare_alid(snp_id))
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

        print("Seeding z/se from source...")
        src_root = zarr.open_group(str(src / "data.zarr"), mode="r")
        src_z = src_root["z"][:].astype(np.float32)
        src_se = src_root["se"][:].astype(np.float32)

        z = np.full((n_variants, n_analyses), np.nan, dtype=np.float32)
        se = np.full((n_variants, n_analyses), np.nan, dtype=np.float32)
        for v in src_variants:
            new_row = new_alid_to_idx[v.alid]
            z[new_row, :] = src_z[v.variant_index, :]
            se[new_row, :] = src_se[v.variant_index, :]

        off_panel_mask = ~on_panel
        if off_panel_mask.any():
            n_missing_off_panel = np.isnan(z[off_panel_mask, :]).sum(axis=0).astype(np.int64)
        else:
            n_missing_off_panel = np.zeros(n_analyses, dtype=np.int64)

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

        # ── Phase 2: parallel LD-block completion ───────────────────────────
        print(f"Running reference completion across {len(tsv_paths):,} LD blocks "
              f"(n_workers={n_workers})...")
        blocks_dir = checkpoint_dir / "blocks"
        blocks_dir.mkdir(parents=True, exist_ok=True)

        results: list[BlockCompletionResult] = []
        pending: list[_BlockTask] = []
        for tsv_path in tsv_paths:
            block_id = f"{tsv_path.parent.name}/{tsv_path.stem}"
            ckpt_path = blocks_dir / f"{_sanitize_block_id(block_id)}.npz"
            if ckpt_path.exists():
                results.append(_load_checkpoint(ckpt_path))
            else:
                pending.append(_BlockTask(
                    tsv_path=tsv_path, source_path=src,
                    min_cor=min_cor, thresh=thresh, checkpoint_path=ckpt_path,
                ))

        if pending:
            print(f"  {len(tsv_paths) - len(pending):,} blocks already checkpointed, "
                  f"{len(pending):,} remaining")
        if n_workers <= 1:
            for i, task in enumerate(pending):
                r = _run_block(task)
                if r is not None:
                    results.append(r)
                if (i + 1) % 200 == 0:
                    print(f"  {i + 1:,} / {len(pending):,} blocks")
        else:
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_block, task): task for task in pending}
                for i, fut in enumerate(as_completed(futures)):
                    r = fut.result()
                    if r is not None:
                        results.append(r)
                    if (i + 1) % 200 == 0:
                        print(f"  {i + 1:,} / {len(pending):,} blocks")

        # ── Phase 3: merge, write, finalise ─────────────────────────────────
        print("Merging block results...")
        imputed = np.zeros((n_variants, n_analyses), dtype=np.uint8)
        quality_rows: list[tuple[int, str, float | None, int, int]] = []
        total_imputed = 0
        for r in results:
            for ai, pearson, nimp, nmiss in r.quality_rows:
                quality_rows.append((ai, r.block_id, pearson, nimp, nmiss))
            for alid, ai, zv, sev in r.fills:
                new_idx = new_alid_to_idx.get(alid)
                if new_idx is None:
                    continue
                if not np.isfinite(z[new_idx, ai]):
                    z[new_idx, ai] = zv
                    se[new_idx, ai] = sev
                    imputed[new_idx, ai] = 1
                    total_imputed += 1

        n_missing_imputation_failed = int(np.isnan(z[on_panel, :]).sum()) if on_panel.any() else 0
        n_missing_off_panel_total = int(n_missing_off_panel.sum())
        print(
            f"Completion done: {total_imputed:,} imputed, "
            f"{n_missing_imputation_failed:,} imputation-failed, "
            f"{n_missing_off_panel_total:,} off-panel missing"
        )

        print("Writing data.zarr...")
        _write_dense_zarr(work, z, se, imputed, on_panel, DEFAULT_CHUNK_SHAPE, DEFAULT_DTYPE)

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


def _write_dense_zarr(
    output_path: Path,
    z: np.ndarray,
    se: np.ndarray,
    imputed: np.ndarray,
    on_panel: np.ndarray,
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    effective_chunks = (min(chunk_shape[0], z.shape[0]), min(chunk_shape[1], z.shape[1]))
    root = zarr.open_group(str(output_path / "data.zarr"), mode="w")
    root.create_dataset(
        "z", data=z.astype(dtype), chunks=effective_chunks, compressor=_COMPRESSOR, dtype=dtype
    )
    root.create_dataset(
        "se", data=se.astype(dtype), chunks=effective_chunks, compressor=_COMPRESSOR, dtype=dtype
    )
    root.create_dataset(
        "imputed", data=imputed, chunks=effective_chunks, compressor=_COMPRESSOR, dtype="uint8"
    )
    root.create_dataset(
        "on_panel", data=on_panel.astype(np.uint8),
        chunks=(effective_chunks[0],), compressor=_COMPRESSOR, dtype="uint8",
    )
    root.attrs["layout"] = "dense"
    root.attrs["completion_state"] = "reference_completed"
    root.attrs["compressor"] = DEFAULT_COMPRESSOR
    root.attrs["chunk_shape"] = list(effective_chunks)
