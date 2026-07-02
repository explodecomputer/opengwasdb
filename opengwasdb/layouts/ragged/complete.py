"""Ragged Reference Completion — enhancement pipeline.

Reads an observed-only ragged store and writes a new Reference-Completed
Store Release by imputing z-scores and SE for all LD reference panel
variants within each analysis's declared cis window.
"""
from __future__ import annotations

import json
import logging
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

from opengwasdb.completion.impute import impute_z_block, scalar_n_se
from opengwasdb.completion.ld_panel import (
    LDBlock,
    find_blocks,
    load_ld_eigenvectors,
    match_variants,
)
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RAGGED_ZARR_PATH, RaggedCSRReader
from opengwasdb.model.enums import CompletionState
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.traits.axis import TraitsAxisReader
from opengwasdb.variants.axis import (
    VariantAxis,
    write_variant_axis,
)
from opengwasdb.variants.normalise import (
    CanonicalVariant,
    VariantNormalisationError,
    chromosome_sort_key,
    normalise_chromosome,
    orient_to_canonical,
)

log = logging.getLogger(__name__)

_COMPRESSOR = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
_ASSOC_CHUNK = 200_000
_OFFSET_CHUNK = 10_000

_LD_PANEL_ID = "eur-hg38-gpm"
_COMPLETION_METHOD = "elastic_net_eigenvectors_v1"


@dataclass(frozen=True)
class CompletionResult:
    output_path: Path
    n_variants: int
    n_analyses: int
    n_associations: int
    n_imputed: int
    n_missing: int


def complete_ragged_store(
    source_path: str | Path,
    dest_path: str | Path,
    ld_dir: str | Path,
    *,
    ancestry: str = "EUR",
    cis_window_bp: int = 1_000_000,
    min_cor: float = 0.7,
    thresh: float = 0.9,
    release_id: str | None = None,
    ld_panel_id: str = _LD_PANEL_ID,
    overwrite: bool = False,
) -> CompletionResult:
    """Produce a Reference-Completed Store Release from an observed-only ragged store.

    source_path: existing observed-only ragged store.
    dest_path:   new store directory to create.
    ld_dir:      root of LD panel; blocks at ld_dir/{ancestry}/{chr}/{block}.*
    """
    src = Path(source_path)
    dst = Path(dest_path)

    if dst.exists() and not overwrite:
        raise FileExistsError(f"Destination already exists: {dst}. Use overwrite=True.")
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    manifest = StoreManifest.load(src)
    print(f"Source store: {manifest.store_id} / {manifest.release_id}")
    print(f"  completion_state: {manifest.completion_state}")

    # ── 1. Load source variant axis and traits ────────────────────────────────
    src_variant_axis = VariantAxis(src)
    src_variants = src_variant_axis.all()
    src_alids: list[str] = [v.alid for v in src_variants]
    alid_to_src_idx: dict[str, int] = {v.alid: v.variant_index for v in src_variants}

    traits_reader = TraitsAxisReader(src)
    trait_records = list(traits_reader.all())
    n_analyses = len(trait_records)

    print(f"Source: {len(src_alids):,} variants, {n_analyses:,} analyses")

    # ── 2. First pass: collect all LD panel variants across all cis windows ───
    print("Collecting LD panel variants across all cis windows...")
    new_alids: set[str] = set()  # LD panel variants not in source store
    new_variant_eaf: dict[str, float] = {}

    for i, rec in enumerate(trait_records):
        if rec.trait_chr is None or rec.trait_bp is None:
            continue
        chrom = rec.trait_chr
        start = max(1, rec.trait_bp - cis_window_bp)
        end = rec.trait_bp + cis_window_bp

        blocks = find_blocks(ld_dir, ancestry, chrom, start, end)
        for block in blocks:
            for ld_row, snp_id in enumerate(block.snp_ids):
                bare = snp_id[3:] if snp_id.startswith("chr") else snp_id
                if bare not in alid_to_src_idx:
                    new_alids.add(bare)
                    if np.isfinite(block.eaf[ld_row]):
                        new_variant_eaf[bare] = float(block.eaf[ld_row])

        if (i + 1) % 1000 == 0:
            print(f"  Scanned {i + 1:,} / {n_analyses:,} analyses, "
                  f"{len(new_alids):,} new LD panel variants found")

    print(f"New LD panel variants: {len(new_alids):,}")

    # ── 3. Build merged variant table ────────────────────────────────────────
    # Parse and normalise the new LD panel variants; merge with source variants.
    new_canonical: list[CanonicalVariant] = []
    new_alid_to_cv: dict[str, CanonicalVariant] = {}
    for alid in new_alids:
        try:
            parts = alid.split(":")
            if len(parts) != 4:
                continue
            chrom, pos_str, a1, a2 = parts
            cv_result = orient_to_canonical(chrom, int(pos_str), a1, a2)
            if cv_result.variant.alid not in alid_to_src_idx:
                new_canonical.append(cv_result.variant)
                new_alid_to_cv[cv_result.variant.alid] = cv_result.variant
        except (VariantNormalisationError, ValueError):
            continue

    # Sort new variants by genomic position
    new_canonical.sort(
        key=lambda v: (chromosome_sort_key(v.chromosome), v.position, v.effect_allele, v.other_allele)
    )

    # Merge: source variants (keep indices) + new (sorted by position)
    # Re-sort ALL variants together for a coherent sorted table.
    all_src = [(chromosome_sort_key(v.chromosome), v.position, v.effect_allele, v.other_allele, v)
               for v in src_variants]
    all_new = [(chromosome_sort_key(v.chromosome), v.position, v.effect_allele, v.other_allele, v)
               for v in new_canonical]
    all_variants_sorted = sorted(all_src + all_new, key=lambda x: x[:4])

    merged_variants: list[CanonicalVariant] = [x[4] for x in all_variants_sorted]
    new_alid_to_idx: dict[str, int] = {v.alid: i for i, v in enumerate(merged_variants)}

    # Build remapping: old source variant_index → new variant_index
    src_idx_remap: list[int] = [0] * len(src_alids)
    for v in src_variants:
        src_idx_remap[v.variant_index] = new_alid_to_idx[v.alid]

    print(f"Merged variant table: {len(merged_variants):,} variants")

    print("Writing variants.tsv.gz...")
    write_variant_axis(dst, merged_variants, {})

    # ── 4. Copy traits axis ───────────────────────────────────────────────────
    shutil.copy2(src / "traits.tsv.gz", dst / "traits.tsv.gz")
    tbi_src = src / "traits.tsv.gz.tbi"
    if tbi_src.exists():
        shutil.copy2(tbi_src, dst / "traits.tsv.gz.tbi")

    # ── 5. Build new SQLite: copy analyses, add completion_quality table ──────
    print("Writing index.sqlite...")
    src_db = sqlite3.connect(str(src / "index.sqlite"))
    src_db.row_factory = sqlite3.Row
    dst_db = sqlite3.connect(str(dst / "index.sqlite"))
    dst_db.row_factory = sqlite3.Row

    # Copy analyses table schema and data
    src_schema = src_db.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='analyses'"
    ).fetchone()
    if src_schema:
        dst_db.execute(src_schema["sql"])
    for row in src_db.execute("SELECT * FROM analyses ORDER BY analysis_index"):
        cols = list(row.keys())
        placeholders = ",".join("?" * len(cols))
        dst_db.execute(
            f"INSERT INTO analyses ({','.join(cols)}) VALUES ({placeholders})",
            list(row),
        )
    dst_db.execute("CREATE INDEX IF NOT EXISTS idx_analyses_trait_id  ON analyses(trait_id)")
    dst_db.execute("CREATE INDEX IF NOT EXISTS idx_analyses_gene_id   ON analyses(gene_id)")
    dst_db.execute("CREATE INDEX IF NOT EXISTS idx_analyses_trait_loc ON analyses(trait_chr, trait_bp)")

    dst_db.execute("""
        CREATE TABLE completion_quality (
            analysis_index INTEGER NOT NULL,
            block_id       TEXT    NOT NULL,
            pearson_r      REAL,
            n_imputed      INTEGER NOT NULL,
            n_missing      INTEGER NOT NULL,
            PRIMARY KEY (analysis_index, block_id)
        )
    """)
    dst_db.commit()
    src_db.close()

    # ── 6. Build completion CSR with imputed mask ─────────────────────────────
    print("Running reference completion per analysis...")
    src_csr = RaggedCSRReader(src)

    all_vi: list[np.ndarray] = []
    all_z: list[np.ndarray] = []
    all_se: list[np.ndarray] = []
    all_imp: list[np.ndarray] = []
    offsets: list[int] = [0]

    quality_rows: list[tuple] = []
    total_imputed = 0
    total_missing = 0

    for ai, rec in enumerate(trait_records):
        obs = src_csr.get_analysis(ai)

        # Remap observed variant indices to new indexing scheme
        obs_vi_new = np.array([src_idx_remap[vi] for vi in obs.variant_index.tolist()], dtype=np.int32)
        obs_z = obs.z.astype(np.float32)
        obs_se = obs.se.astype(np.float32)

        if rec.trait_chr is None or rec.trait_bp is None:
            # No cis region — write observed associations as-is, no completion.
            all_vi.append(obs_vi_new)
            all_z.append(obs_z)
            all_se.append(obs_se)
            all_imp.append(np.zeros(len(obs_vi_new), dtype=np.uint8))
            offsets.append(offsets[-1] + len(obs_vi_new))
            continue

        chrom = rec.trait_chr
        cis_start = max(1, rec.trait_bp - cis_window_bp)
        cis_end = rec.trait_bp + cis_window_bp

        blocks = find_blocks(ld_dir, ancestry, chrom, cis_start, cis_end)
        if not blocks:
            all_vi.append(obs_vi_new)
            all_z.append(obs_z)
            all_se.append(obs_se)
            all_imp.append(np.zeros(len(obs_vi_new), dtype=np.uint8))
            offsets.append(offsets[-1] + len(obs_vi_new))
            continue

        # Build per-analysis ALID → z/se lookup from observed data
        obs_alid_to_z: dict[str, float] = {}
        obs_alid_to_se: dict[str, float] = {}
        for vi_old, z_val, se_val in zip(
            obs.variant_index.tolist(), obs.z.tolist(), obs.se.tolist()
        ):
            alid = src_alids[vi_old]
            obs_alid_to_z[alid] = float(z_val)
            obs_alid_to_se[alid] = float(se_val)

        # Collect all reference panel variants in cis window
        ref_alids_in_cis: list[str] = []
        for block in blocks:
            for snp_id in block.snp_ids:
                bare = snp_id[3:] if snp_id.startswith("chr") else snp_id
                ref_alids_in_cis.append(bare)

        # Deduplicate (variants can appear in multiple blocks)
        seen: set[str] = set()
        unique_ref_alids: list[str] = []
        for alid in ref_alids_in_cis:
            if alid not in seen:
                seen.add(alid)
                unique_ref_alids.append(alid)

        # Process each block: run imputation, record quality
        alid_imputed_z: dict[str, float] = {}
        alid_imputed_se: dict[str, float] = {}

        for block in blocks:
            bl_alids = [
                (snp_id[3:] if snp_id.startswith("chr") else snp_id)
                for snp_id in block.snp_ids
            ]

            # Build dense z-score vector for this block
            z_dense = np.array(
                [obs_alid_to_z.get(a, float("nan")) for a in bl_alids],
                dtype=np.float64,
            )
            n_obs_in_block = int(np.isfinite(z_dense).sum())

            try:
                eigenvalues, eigenvectors = load_ld_eigenvectors(block, thresh)
            except Exception as exc:  # noqa: BLE001
                log.warning("Block %s: cannot load eigenvectors (%s) — skipping", block.block_id, exc)
                quality_rows.append((ai, block.block_id, None, 0, int((~np.isfinite(z_dense)).sum())))
                continue

            if n_obs_in_block < 2:
                n_miss = int((~np.isfinite(z_dense)).sum())
                quality_rows.append((ai, block.block_id, None, 0, n_miss))
                continue

            z_imp_arr, corr = impute_z_block(z_dense, eigenvectors, eigenvalues, min_cor=min_cor)

            n_miss_block = int((~np.isfinite(z_dense)).sum())
            if z_imp_arr is None:
                quality_rows.append((ai, block.block_id, float(corr) if np.isfinite(corr) else None, 0, n_miss_block))
                continue

            # SE via scalar-N model using LD panel EAF
            eaf_obs_in_block = np.array(
                [float(block.eaf[i]) for i in range(block.n_ld_snps)],
                dtype=np.float64,
            )
            se_obs_in_block = np.array(
                [obs_alid_to_se.get(a, float("nan")) for a in bl_alids],
                dtype=np.float64,
            )
            se_all = scalar_n_se(
                se_obs_in_block[np.isfinite(z_dense)],
                eaf_obs_in_block[np.isfinite(z_dense)],
                eaf_obs_in_block,
            )

            # Fill missing positions
            missing_mask = ~np.isfinite(z_dense)
            n_filled = 0
            for i, (alid, is_miss, z_v, se_v) in enumerate(
                zip(bl_alids, missing_mask, z_imp_arr, se_all)
            ):
                if is_miss and np.isfinite(z_v) and np.isfinite(se_v):
                    alid_imputed_z[alid] = float(z_v)
                    alid_imputed_se[alid] = float(se_v)
                    n_filled += 1

            quality_rows.append((ai, block.block_id, float(corr), n_filled, n_miss_block - n_filled))
            total_imputed += n_filled

        # Assemble merged association sequence for this analysis
        # All reference panel variants in cis window + any observed off-panel
        ref_vi: list[int] = []
        ref_z: list[float] = []
        ref_se: list[float] = []
        ref_imp: list[int] = []

        seen_alids: set[str] = set()
        # Reference panel variants (observed or imputed or missing)
        for alid in unique_ref_alids:
            vi = new_alid_to_idx.get(alid)
            if vi is None:
                continue
            seen_alids.add(alid)
            if alid in obs_alid_to_z:
                ref_vi.append(vi)
                ref_z.append(obs_alid_to_z[alid])
                ref_se.append(obs_alid_to_se[alid])
                ref_imp.append(0)
            elif alid in alid_imputed_z:
                ref_vi.append(vi)
                ref_z.append(alid_imputed_z[alid])
                ref_se.append(alid_imputed_se[alid])
                ref_imp.append(1)
            else:
                # Missing reference variant
                ref_vi.append(vi)
                ref_z.append(float("nan"))
                ref_se.append(float("nan"))
                ref_imp.append(0)
                total_missing += 1

        # Off-panel observed variants (observed but not in reference panel)
        for vi_old, z_val, se_val in zip(
            obs.variant_index.tolist(), obs.z.tolist(), obs.se.tolist()
        ):
            alid = src_alids[vi_old]
            if alid not in seen_alids:
                ref_vi.append(new_alid_to_idx[alid])
                ref_z.append(float(z_val))
                ref_se.append(float(se_val))
                ref_imp.append(0)

        # Sort by variant_index
        order = np.argsort(ref_vi)
        vi_arr = np.array(ref_vi, dtype=np.int32)[order]
        z_arr = np.array(ref_z, dtype=np.float16)[order]
        se_arr = np.array(ref_se, dtype=np.float16)[order]
        imp_arr = np.array(ref_imp, dtype=np.uint8)[order]

        all_vi.append(vi_arr)
        all_z.append(z_arr)
        all_se.append(se_arr)
        all_imp.append(imp_arr)
        offsets.append(offsets[-1] + len(vi_arr))

        if (ai + 1) % 500 == 0:
            print(f"  {ai + 1:,} / {n_analyses:,} analyses | "
                  f"imputed {total_imputed:,} | missing {total_missing:,}")

    print(f"Completion done: {total_imputed:,} imputed, {total_missing:,} missing")

    # ── 7. Write zarr CSR with imputed mask ───────────────────────────────────
    print("Writing zarr CSR...")
    n_total = offsets[-1]
    offsets_arr = np.array(offsets, dtype=np.int64)
    vi_all = np.concatenate(all_vi) if all_vi else np.empty(0, dtype=np.int32)
    z_all = np.concatenate(all_z) if all_z else np.empty(0, dtype=np.float16)
    se_all = np.concatenate(all_se) if all_se else np.empty(0, dtype=np.float16)
    imp_all = np.concatenate(all_imp) if all_imp else np.empty(0, dtype=np.uint8)

    ragged_path = dst / RAGGED_ZARR_PATH
    ragged_path.mkdir(parents=True, exist_ok=True)
    root = zarr.open_group(str(ragged_path), mode="w")
    root.create_dataset("offsets", data=offsets_arr, chunks=(_OFFSET_CHUNK,), compressor=_COMPRESSOR, dtype=np.int64)
    root.create_dataset("variant_index", data=vi_all, chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.int32)
    root.create_dataset("z", data=z_all, chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.float16)
    root.create_dataset("se", data=se_all, chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.float16)
    root.create_dataset("imputed", data=imp_all, chunks=(_ASSOC_CHUNK,), compressor=_COMPRESSOR, dtype=np.uint8)
    root.attrs["layout"] = "ragged"
    root.attrs["completion_state"] = "reference_completed"
    root.attrs["n_analyses"] = n_analyses
    root.attrs["n_associations"] = n_total

    # ── 8. Write completion quality rows ──────────────────────────────────────
    print(f"Writing {len(quality_rows):,} completion quality rows...")
    dst_db.executemany(
        "INSERT INTO completion_quality (analysis_index, block_id, pearson_r, n_imputed, n_missing) "
        "VALUES (?, ?, ?, ?, ?)",
        quality_rows,
    )
    dst_db.commit()
    dst_db.close()

    # ── 9. Build top-hit indexes ──────────────────────────────────────────────
    print("Building top-hit indexes...")
    build_ragged_top_hit_indexes(dst)

    # ── 10. Write manifest ────────────────────────────────────────────────────
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
                "cis_window_bp": cis_window_bp,
                "min_cor": min_cor,
                "pca_thresh": thresh,
                "n_variants_total": len(merged_variants),
                "n_variants_new": len(new_canonical),
                "n_imputed": total_imputed,
                "n_missing": total_missing,
            },
        },
    )
    (dst / "manifest.json").write_text(
        json.dumps(completed_manifest.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    src_variant_axis.close()
    traits_reader.close()

    result = CompletionResult(
        output_path=dst,
        n_variants=len(merged_variants),
        n_analyses=n_analyses,
        n_associations=n_total,
        n_imputed=total_imputed,
        n_missing=total_missing,
    )
    print(
        f"Reference completion complete: {result.n_variants:,} variants, "
        f"{result.n_analyses:,} analyses, {result.n_associations:,} associations "
        f"({result.n_imputed:,} imputed, {result.n_missing:,} missing)"
    )
    return result
