"""Per-LD-block completion checkpoints -- the on-disk record one process-pool
worker writes for one block, read back by the parent in Phase 3. Shared by
dense and ragged completion so a worker crash or a resumed run behaves
identically in both: each block's result lives at its own path, written
atomically, independent of every other block's.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ALID strings are short and fixed-format (chr:pos:a1:a2); a bounded numpy
# byte-string dtype avoids the overhead of an object array in a checkpoint
# that may hold millions of fill rows.
ALID_DTYPE = "S64"


def checkpoint_dir_for(dest_path: Path) -> Path:
    dest_path = Path(dest_path)
    return dest_path.parent / f".{dest_path.name}.checkpoint"


def sanitize_block_id(block_id: str) -> str:
    return block_id.replace("/", "__")


@dataclass(frozen=True)
class BlockCompletionResult:
    block_id: str
    # (analysis_index, pearson_r | None, n_imputed, n_missing)
    quality_rows: list[tuple[int, float | None, int, int]]
    # (alid, analysis_index, z, se)
    fills: list[tuple[str, int, float, float]]


def write_block_checkpoint(path: Path, result: BlockCompletionResult) -> None:
    q_ai = np.array([r[0] for r in result.quality_rows], dtype=np.int32)
    q_pearson = np.array(
        [r[1] if r[1] is not None else np.nan for r in result.quality_rows], dtype=np.float64
    )
    q_nimp = np.array([r[2] for r in result.quality_rows], dtype=np.int32)
    q_nmiss = np.array([r[3] for r in result.quality_rows], dtype=np.int32)
    f_alid = np.array([r[0].encode("ascii") for r in result.fills], dtype=ALID_DTYPE)
    f_ai = np.array([r[1] for r in result.fills], dtype=np.int32)
    f_z = np.array([r[2] for r in result.fills], dtype=np.float32)
    f_se = np.array([r[3] for r in result.fills], dtype=np.float32)

    tmp_path = path.with_name(path.name + ".tmp")
    with open(tmp_path, "wb") as fh:
        np.savez(
            fh,
            block_id=np.array([result.block_id]),
            q_ai=q_ai, q_pearson=q_pearson, q_nimp=q_nimp, q_nmiss=q_nmiss,
            f_alid=f_alid, f_ai=f_ai, f_z=f_z, f_se=f_se,
        )
    os.replace(tmp_path, path)


def read_block_checkpoint(path: Path) -> BlockCompletionResult:
    with np.load(path, allow_pickle=False) as d:
        block_id = str(d["block_id"][0])
        quality_rows = [
            (int(ai), None if not np.isfinite(p) else float(p), int(ni), int(nm))
            for ai, p, ni, nm in zip(
                d["q_ai"], d["q_pearson"], d["q_nimp"], d["q_nmiss"], strict=True
            )
        ]
        f_alid = d["f_alid"]
        alids = [a.decode("ascii") if isinstance(a, bytes) else str(a) for a in f_alid]
        fills = [
            (alid, int(ai), float(z), float(se))
            for alid, ai, z, se in zip(alids, d["f_ai"], d["f_z"], d["f_se"], strict=True)
        ]
    return BlockCompletionResult(block_id=block_id, quality_rows=quality_rows, fills=fills)
