"""The ``completion_quality`` SQLite schema shared by Dense and Ragged Reference
Completion (issue #22) and checked by the validator -- one definition instead
of three independently typed-out copies.
"""
from __future__ import annotations

from typing import Any

import numpy as np

COMPLETION_QUALITY_COLUMNS = frozenset(
    {"analysis_index", "block_id", "pearson_r", "n_imputed", "n_missing"}
)

CREATE_COMPLETION_QUALITY_SQL = """
    CREATE TABLE completion_quality (
        analysis_index INTEGER NOT NULL,
        block_id       TEXT    NOT NULL,
        pearson_r      REAL,
        n_imputed      INTEGER NOT NULL,
        n_missing      INTEGER NOT NULL,
        PRIMARY KEY (analysis_index, block_id)
    )
"""


def create_completion_quality_table(connection: Any) -> None:
    connection.execute(CREATE_COMPLETION_QUALITY_SQL)


def completion_quality_rollup(connection: Any, n_analyses: int) -> dict[int, tuple[str, str, str]]:
    """Per-Analysis ``(completion_median_pearson_r, completion_n_imputed_total,
    completion_n_missing_total)`` strings, rolled up from the LD-block-by-Analysis
    ``completion_quality`` table (issue #22) -- the large, fine-grained table
    stays SQLite-only per ADR 0030; only this rollup travels into
    ``analyses.tsv``. Shared by Dense and Ragged completion (issue #70): Dense
    computes its own ``completion_n_missing_total`` from off-panel missingness
    instead (a different quantity -- variants never attempted, not variants
    attempted-and-failed) and ignores this tuple's third element; Ragged has
    no separate off-panel concept, so it uses this rollup's value directly.
    """
    per_analysis_r: dict[int, list[float]] = {}
    per_analysis_imputed: dict[int, int] = {}
    per_analysis_missing: dict[int, int] = {}
    for row in connection.execute(
        "SELECT analysis_index, pearson_r, n_imputed, n_missing FROM completion_quality"
    ):
        idx = int(row["analysis_index"])
        if row["pearson_r"] is not None:
            per_analysis_r.setdefault(idx, []).append(float(row["pearson_r"]))
        per_analysis_imputed[idx] = per_analysis_imputed.get(idx, 0) + int(row["n_imputed"])
        per_analysis_missing[idx] = per_analysis_missing.get(idx, 0) + int(row["n_missing"])
    rollup: dict[int, tuple[str, str, str]] = {}
    for idx in range(n_analyses):
        r_values = per_analysis_r.get(idx, [])
        median_r = f"{float(np.median(r_values)):.6g}" if r_values else ""
        rollup[idx] = (
            median_r,
            str(per_analysis_imputed.get(idx, 0)),
            str(per_analysis_missing.get(idx, 0)),
        )
    return rollup
