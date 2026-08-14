"""The ``completion_quality`` SQLite schema shared by Dense and Ragged Reference
Completion (issue #22) and checked by the validator -- one definition instead
of three independently typed-out copies.
"""
from __future__ import annotations

from typing import Any

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
