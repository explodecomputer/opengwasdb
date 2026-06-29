"""SQLite schema and lookup helpers for Store Releases."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, cast

from opengwasdb.variants import VariantNormalisationError
from opengwasdb.variants.normalise import normalise_allele, normalise_chromosome


def connect(path: str | Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    return connection


def initialise_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS variants (
            variant_index INTEGER PRIMARY KEY,
            alid TEXT NOT NULL,
            chromosome TEXT NOT NULL,
            position INTEGER NOT NULL,
            effect_allele TEXT NOT NULL,
            other_allele TEXT NOT NULL,
            rsid TEXT
        );

        CREATE TABLE IF NOT EXISTS variant_aliases (
            alias TEXT PRIMARY KEY,
            variant_index INTEGER NOT NULL REFERENCES variants(variant_index)
        );

        CREATE TABLE IF NOT EXISTS analyses (
            analysis_index INTEGER PRIMARY KEY,
            analysis_id TEXT NOT NULL UNIQUE,
            phenotype_id TEXT,
            phenotype_label TEXT,
            analysis_label TEXT,
            stored_effect_scale TEXT NOT NULL
        );
        """
    )


def create_lookup_indexes(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE INDEX IF NOT EXISTS idx_variants_range
            ON variants(chromosome, position, variant_index);
        """
    )


def set_metadata(connection: sqlite3.Connection, key: str, value: Any) -> None:
    payload = json.dumps(value, sort_keys=True)
    connection.execute(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        (key, payload),
    )


def get_metadata(connection: sqlite3.Connection, key: str, default: Any = None) -> Any:
    row = connection.execute("SELECT value FROM metadata WHERE key = ?", (key,)).fetchone()
    if row is None:
        return default
    return json.loads(row["value"])


def count_rows(connection: sqlite3.Connection, table: str) -> int:
    row = connection.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()
    return int(row["n"])


def variant_by_identifier(connection: sqlite3.Connection, identifier: str) -> sqlite3.Row | None:
    parsed = _parse_canonical_alid(identifier)
    if parsed is not None:
        chromosome, position, effect_allele, other_allele = parsed
        row = connection.execute(
            """
            SELECT *
            FROM variants
            WHERE chromosome = ?
              AND position = ?
              AND effect_allele = ?
              AND other_allele = ?
            """,
            (chromosome, position, effect_allele, other_allele),
        ).fetchone()
        if row is not None:
            return cast(sqlite3.Row, row)
    alias = connection.execute(
        """
        SELECT v.*
        FROM variant_aliases a
        JOIN variants v ON v.variant_index = a.variant_index
        WHERE a.alias = ?
        """,
        (identifier,),
    ).fetchone()
    return cast(sqlite3.Row | None, alias)


def _parse_canonical_alid(identifier: str) -> tuple[str, int, str, str] | None:
    parts = identifier.split(":")
    if len(parts) != 4:
        return None
    chromosome, position_text, effect_allele, other_allele = parts
    try:
        position = int(position_text)
        if position <= 0:
            return None
        return (
            normalise_chromosome(chromosome),
            position,
            normalise_allele(effect_allele),
            normalise_allele(other_allele),
        )
    except (TypeError, ValueError, VariantNormalisationError):
        return None


def analysis_by_id(connection: sqlite3.Connection, analysis_id: str) -> sqlite3.Row | None:
    return cast(sqlite3.Row | None, connection.execute(
        "SELECT * FROM analyses WHERE analysis_id = ?",
        (analysis_id,),
    ).fetchone())
