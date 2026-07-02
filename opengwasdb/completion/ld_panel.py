"""Read-only interface for the LD reference panel used in Reference Completion.

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

from opengwasdb.completion.impute import ld_pca

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


def load_block(tsv_path: Path) -> LDBlock | None:
    """Load one LD block from its TSV path. Returns None if unreadable/incomplete.

    chrom is taken from the parent directory name, matching the panel's
    ld_dir/{ancestry}/{chr}/{block_name}.tsv layout.
    """
    tsv_path = Path(tsv_path)
    ld_path = tsv_path.with_suffix(".unphased.vcor1.gz")
    if not ld_path.exists():
        return None

    try:
        snp_ids, eaf, blk_start, blk_end = _read_block_tsv(tsv_path)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not read %s: %s — skipping", tsv_path, exc)
        return None

    chrom = tsv_path.parent.name
    npz_path = tsv_path.with_suffix(".ldeig.npz")
    return LDBlock(
        block_id=f"{chrom}/{tsv_path.stem}",
        chrom=chrom,
        start_bp=blk_start,
        end_bp=blk_end,
        tsv_path=tsv_path,
        ld_path=ld_path,
        ldeig_npz_path=npz_path if npz_path.exists() else None,
        n_ld_snps=len(snp_ids),
        snp_ids=snp_ids,
        eaf=eaf,
    )


def find_blocks(
    ld_dir: str | Path,
    ancestry: str,
    chrom: str,
    start_bp: int,
    end_bp: int,
) -> list[LDBlock]:
    """Return all LD blocks whose genomic extent intersects [start_bp, end_bp]."""
    blocks: list[LDBlock] = []
    for block in list_all_blocks(ld_dir, ancestry, chrom):
        if block.end_bp < start_bp or block.start_bp > end_bp:
            continue
        blocks.append(block)
    return blocks


def list_chromosomes(ld_dir: str | Path, ancestry: str) -> list[str]:
    """Return the chromosome names present in the panel for *ancestry*."""
    ancestry_dir = Path(ld_dir) / ancestry
    if not ancestry_dir.is_dir():
        return []
    return sorted(
        (p.name for p in ancestry_dir.iterdir() if p.is_dir()),
        key=lambda c: (len(c), c),
    )


def list_all_blocks(ld_dir: str | Path, ancestry: str, chrom: str) -> list[LDBlock]:
    """Return every LD block for one chromosome, genome-wide (no region filter).

    Used for Full Coverage Dense completion, where the block set is the
    reference panel's entire block list rather than a per-analysis cis window.
    """
    panel_dir = Path(ld_dir) / ancestry / chrom
    if not panel_dir.is_dir():
        return []

    blocks: list[LDBlock] = []
    for tsv_path in sorted(panel_dir.glob("*.tsv")):
        block = load_block(tsv_path)
        if block is not None:
            blocks.append(block)
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
