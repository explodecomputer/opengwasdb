"""The ``analyses`` SQLite table schema shared by both ragged builders
(``build_ragged_from_besd``, ``build_ragged_from_ssf``) -- one definition
instead of two independently typed-out copies, which had already started to
drift once ``assigned_ancestry`` (ADR 0028) needed adding to both.
"""
from __future__ import annotations

import sqlite3

CREATE_ANALYSES_TABLE_SQL = """
    CREATE TABLE analyses (
        analysis_index    INTEGER PRIMARY KEY,
        analysis_id       TEXT NOT NULL,
        trait_id          TEXT NOT NULL,
        gene_id           TEXT,
        gene_name         TEXT,
        tissue            TEXT,
        context           TEXT,
        trait_chr         TEXT,
        trait_bp          INTEGER,
        n                 INTEGER,
        assigned_ancestry TEXT
    )
"""

_INSERT_ANALYSIS_SQL = """
    INSERT INTO analyses
      (analysis_index, analysis_id, trait_id, gene_id, gene_name,
       tissue, context, trait_chr, trait_bp, n, assigned_ancestry)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def create_analyses_table(conn: sqlite3.Connection) -> None:
    conn.execute(CREATE_ANALYSES_TABLE_SQL)
    conn.execute("CREATE INDEX idx_analyses_trait_id  ON analyses(trait_id)")
    conn.execute("CREATE INDEX idx_analyses_gene_id   ON analyses(gene_id)")
    conn.execute("CREATE INDEX idx_analyses_trait_loc ON analyses(trait_chr, trait_bp)")


def insert_analysis_row(
    conn: sqlite3.Connection,
    *,
    analysis_index: int,
    analysis_id: str,
    trait_id: str,
    gene_id: str | None,
    gene_name: str | None,
    tissue: str | None,
    context: str | None,
    trait_chr: str | None,
    trait_bp: int | None,
    n: int | None,
    assigned_ancestry: str | None,
) -> None:
    conn.execute(
        _INSERT_ANALYSIS_SQL,
        (analysis_index, analysis_id, trait_id, gene_id, gene_name,
         tissue, context, trait_chr, trait_bp, n, assigned_ancestry or None),
    )
