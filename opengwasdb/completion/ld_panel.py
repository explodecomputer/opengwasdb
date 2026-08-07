"""Read-only interface for the LD reference panel used in Reference Completion.

Panel layout (flat/production):
    ld_dir/{ancestry}/{chr}/{block_name}.tsv                  (variant table, required)
  and at least one of:
    ld_dir/{ancestry}/{chr}/{block_name}.ldeig.npz             (eigendecomposition; preferred)
    ld_dir/{ancestry}/{chr}/{block_name}.unphased.vcor1.gz     (LD matrix; optional/legacy,
                                                                 used to derive eigenvectors
                                                                 only when the npz is absent)

A block with neither the npz nor the matrix cannot yield eigenvectors and is
skipped. See ADR 0031 and store-format spec §13.1.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from opengwasdb.completion.impute import ld_pca
from opengwasdb.variants import VariantNormalisationError, orient_to_canonical

log = logging.getLogger(__name__)


@dataclass
class LDBlock:
    block_id: str       # "{chr}/{block_name}"
    chrom: str
    start_bp: int
    end_bp: int
    tsv_path: Path
    ld_path: Path | None       # LD matrix; optional, may be absent (eigendecomposition-only panel)
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

    A block is admitted when it can produce eigenvectors by either route — a
    stored eigendecomposition or a raw LD matrix — and skipped only when
    neither artifact is present.
    """
    tsv_path = Path(tsv_path)
    ld_path_candidate = tsv_path.with_suffix(".unphased.vcor1.gz")
    ld_path = ld_path_candidate if ld_path_candidate.exists() else None
    npz_path_candidate = tsv_path.with_suffix(".ldeig.npz")
    ldeig_npz_path = npz_path_candidate if npz_path_candidate.exists() else None
    if ld_path is None and ldeig_npz_path is None:
        return None

    try:
        snp_ids, eaf, blk_start, blk_end = _read_block_tsv(tsv_path)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not read %s: %s — skipping", tsv_path, exc)
        return None

    chrom = tsv_path.parent.name
    return LDBlock(
        block_id=f"{chrom}/{tsv_path.stem}",
        chrom=chrom,
        start_bp=blk_start,
        end_bp=blk_end,
        tsv_path=tsv_path,
        ld_path=ld_path,
        ldeig_npz_path=ldeig_npz_path,
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


def canonical_panel_alid(snp_id: str) -> str | None:
    """Normalise an LD-panel SNP id to a canonical store ALID (``chr:pos:a1:a2``).

    Handles both the panel's ``[chr]CHR:POS_REF_ALT`` form (underscores) and an
    already-colon-delimited ``CHR:POS:A1:A2`` form, orienting alleles to the
    canonical A1 = min(ref, alt). Returns None — after logging a warning
    identifying the offending id — when it doesn't parse under either
    convention, so a panel-wide format mistake surfaces as warnings rather
    than as silently reduced completion coverage.
    """
    s = _strip_chr(snp_id)
    parts = s.split(":")
    if len(parts) == 4:
        chrom, pos_s, ref, alt = parts
    elif len(parts) == 2 and parts[1].count("_") == 2:
        chrom, rest = parts
        pos_s, ref, alt = rest.split("_")
    else:
        log.warning("Could not parse LD panel SNP id %r — unrecognised format", snp_id)
        return None
    try:
        return orient_to_canonical(chrom, int(pos_s), ref, alt).variant.alid
    except (VariantNormalisationError, ValueError) as exc:
        log.warning("Could not canonicalise LD panel SNP id %r: %s", snp_id, exc)
        return None


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
            stored_k = vecs.shape[1]
            if k > stored_k:
                achieved = float(cumvar[stored_k - 1]) if stored_k > 0 else 0.0
                log.warning(
                    "Block %s: stored eigendecomposition has only %d component(s), "
                    "needs %d to reach thresh=%.3f cumulative variance (achieved %.4f) "
                    "— using what's stored",
                    block.block_id, stored_k, k, thresh, achieved,
                )
                k = stored_k
            return vals[:k], vecs[:, :k]
        except Exception as exc:  # noqa: BLE001
            log.warning(
                "Could not load %s: %s — falling back to LD matrix", block.ldeig_npz_path, exc
            )

    return _load_from_ld_matrix(block, thresh)


def _load_from_ld_matrix(block: LDBlock, thresh: float) -> tuple[np.ndarray, np.ndarray]:
    if block.ld_path is None:
        raise RuntimeError(
            f"Block {block.block_id}: no LD matrix available and no usable "
            ".ldeig.npz eigendecomposition"
        )
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
    if block.ld_path is None:
        raise RuntimeError(f"Block {block.block_id} has no LD matrix (ld_path is None)")
    import pandas as pd
    full = pd.read_csv(
        block.ld_path, header=None, delimiter="\t"
    ).values.astype(np.float64)
    idx = np.array(row_indices)
    return full[np.ix_(idx, idx)]
