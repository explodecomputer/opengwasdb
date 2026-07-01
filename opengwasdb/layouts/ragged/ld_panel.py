"""Read-only interface for the LD reference panel used in Ragged Reference Completion.

Panel layout (flat/production):
    ld_dir/{ancestry}/{chr}/{block_name}.tsv
    ld_dir/{ancestry}/{chr}/{block_name}.unphased.vcor1.gz
    ld_dir/{ancestry}/{chr}/{block_name}.ldeig.npz   (optional cache)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from opengwasdb.layouts.ragged.impute import ld_pca

log = logging.getLogger(__name__)


@dataclass
class LDBlock:
    block_id: str       # "{chr}/{block_name}"
    chrom: str
    start_bp: int
    end_bp: int
    tsv_path: Path
    ld_path: Path
    ldeig_npz_path: Path | None
    n_ld_snps: int
    snp_ids: list[str]  # ALID strings as stored in the TSV (may have chr prefix)
    eaf: np.ndarray     # float64, length n_ld_snps; NaN where absent


def _strip_chr(s: str) -> str:
    return s[3:] if s.startswith("chr") else s


def _read_block_tsv(tsv_path: Path) -> tuple[list[str], np.ndarray, int, int]:
    """Parse a block TSV and return (snp_ids, eaf, min_bp, max_bp)."""
    import csv

    snp_ids: list[str] = []
    eafs: list[float] = []
    bps: list[int] = []

    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            snp_ids.append(row["SNP"])
            try:
                eafs.append(float(row.get("EAF") or "nan"))
            except ValueError:
                eafs.append(float("nan"))
            try:
                bps.append(int(row["BP"]))
            except (KeyError, ValueError):
                bps.append(0)

    eaf_arr = np.array(eafs, dtype=np.float64)
    min_bp = int(min(bps)) if bps else 0
    max_bp = int(max(bps)) if bps else 0
    return snp_ids, eaf_arr, min_bp, max_bp


def find_blocks(
    ld_dir: str | Path,
    ancestry: str,
    chrom: str,
    start_bp: int,
    end_bp: int,
) -> list[LDBlock]:
    """Return all LD blocks whose genomic extent intersects [start_bp, end_bp]."""
    panel_dir = Path(ld_dir) / ancestry / chrom
    if not panel_dir.is_dir():
        return []

    blocks: list[LDBlock] = []
    for tsv_path in sorted(panel_dir.glob("*.tsv")):
        block_name = tsv_path.stem
        ld_path = tsv_path.with_suffix(".unphased.vcor1.gz")
        if not ld_path.exists():
            continue

        try:
            snp_ids, eaf, blk_start, blk_end = _read_block_tsv(tsv_path)
        except Exception as exc:  # noqa: BLE001
            log.warning("Could not read %s: %s — skipping", tsv_path, exc)
            continue

        if blk_end < start_bp or blk_start > end_bp:
            continue

        npz_path = tsv_path.with_suffix(".ldeig.npz")
        blocks.append(LDBlock(
            block_id=f"{chrom}/{block_name}",
            chrom=chrom,
            start_bp=blk_start,
            end_bp=blk_end,
            tsv_path=tsv_path,
            ld_path=ld_path,
            ldeig_npz_path=npz_path if npz_path.exists() else None,
            n_ld_snps=len(snp_ids),
            snp_ids=snp_ids,
            eaf=eaf,
        ))

    return blocks


def match_variants(
    block: LDBlock,
    store_alids: list[str],
) -> tuple[list[int], list[int]]:
    """Match store ALIDs to LD panel row indices.

    Returns (store_variant_positions, ld_row_indices) — parallel index arrays.
    Returns empty lists when fewer than 2 variants match.
    """
    alid_to_store: dict[str, int] = {alid: i for i, alid in enumerate(store_alids)}

    store_positions: list[int] = []
    ld_row_indices: list[int] = []

    for ld_row, snp_id in enumerate(block.snp_ids):
        bare = _strip_chr(snp_id)
        if bare in alid_to_store:
            store_positions.append(alid_to_store[bare])
            ld_row_indices.append(ld_row)

    if len(store_positions) < 2:
        return [], []
    return store_positions, ld_row_indices


def load_ld_eigenvectors(
    block: LDBlock,
    thresh: float = 0.9,
) -> tuple[np.ndarray, np.ndarray]:
    """Load (eigenvalues, eigenvectors) for a block.

    Uses the .ldeig.npz cache when available; falls back to loading the full
    LD matrix and running ld_pca.
    """
    if block.ldeig_npz_path is not None:
        try:
            data = np.load(str(block.ldeig_npz_path))
            vals: np.ndarray = data["values"].astype(np.float64)
            vecs: np.ndarray = data["vectors"].astype(np.float64)
            total = float(np.maximum(vals, 0).sum()) or 1.0
            cumvar = np.cumsum(np.maximum(vals, 0)) / total
            k = int(np.searchsorted(cumvar, thresh)) + 1
            k = min(k, vecs.shape[1])
            return vals[:k], vecs[:, :k]
        except Exception as exc:  # noqa: BLE001
            log.warning("Could not load %s: %s — falling back to LD matrix", block.ldeig_npz_path, exc)

    return _load_from_ld_matrix(block, thresh)


def _load_from_ld_matrix(block: LDBlock, thresh: float) -> tuple[np.ndarray, np.ndarray]:
    import pandas as pd
    try:
        full = pd.read_csv(
            block.ld_path, header=None, delimiter="\t"
        ).values.astype(np.float64)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Cannot load LD matrix {block.ld_path}: {exc}") from exc
    return ld_pca(full, thresh)


def load_ld_matrix(block: LDBlock, row_indices: list[int]) -> np.ndarray:
    """Load LD submatrix for the given row indices."""
    import pandas as pd
    full = pd.read_csv(
        block.ld_path, header=None, delimiter="\t"
    ).values.astype(np.float64)
    idx = np.array(row_indices)
    return full[np.ix_(idx, idx)]
