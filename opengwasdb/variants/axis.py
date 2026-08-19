"""Tabix-backed Store Variant Axis."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pysam

from opengwasdb.variants.normalise import (
    CanonicalVariant,
    VariantNormalisationError,
    normalise_allele,
    normalise_chromosome,
)

VARIANT_TABLE_FILENAME = "variants.tsv.gz"
VARIANT_TABIX_FILENAME = "variants.tsv.gz.tbi"
VARIANT_OFFSETS_FILENAME = "variant_offsets.npy"
VARIANT_ALID_BYTES_FILENAME = "variant_alid_bytes.npy"
VARIANT_ALID_ROWS_FILENAME = "variant_alid_rows.npy"
VARIANT_AXIS_FORMAT = "tabix_tsv_v1"
VARIANT_HEADER = (
    "#chromosome\tposition\tvariant_index\teffect_allele\tother_allele\talid\trsid"
    "\tsource_alid\n"
)
# Fixed byte-width for ALID encoding in the mmap'd search index.
# Supports chromosomes up to 3 chars, positions up to 9 digits, alleles up to ~20 chars.
_ALID_DTYPE = "|S64"

# `by_indices()` switchover point between random-access and full-scan
# resolution -- see that method's docstring for the measured costs behind it.
_BULK_SCAN_FRACTION = 0.01


@dataclass(frozen=True)
class VariantRecord(Mapping[str, Any]):
    """One Store Variant Table row."""

    variant_index: int
    alid: str
    chromosome: str
    position: int
    effect_allele: str
    other_allele: str
    rsid: str | None
    # Canonical ALID of the source association on the *build* (pre-liftover)
    # assembly that occupies this row — the store's provenance link back to the
    # original dataset. None for same-assembly builds (no liftover) and for rows
    # where several source variants collided onto one lifted ALID (ambiguous).
    source_alid: str | None = None

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __iter__(self) -> Iterator[str]:
        return iter(
            (
                "variant_index",
                "alid",
                "chromosome",
                "position",
                "effect_allele",
                "other_allele",
                "rsid",
                "source_alid",
            )
        )

    def __len__(self) -> int:
        return 8


@dataclass(frozen=True)
class ParsedAlid:
    chromosome: str
    position: int
    effect_allele: str
    other_allele: str


def variant_table_path(store_path: str | Path) -> Path:
    return Path(store_path) / VARIANT_TABLE_FILENAME


def variant_tabix_path(store_path: str | Path) -> Path:
    return Path(store_path) / VARIANT_TABIX_FILENAME


def variant_offsets_path(store_path: str | Path) -> Path:
    return Path(store_path) / VARIANT_OFFSETS_FILENAME


def variant_alid_bytes_path(store_path: str | Path) -> Path:
    return Path(store_path) / VARIANT_ALID_BYTES_FILENAME


def variant_alid_rows_path(store_path: str | Path) -> Path:
    return Path(store_path) / VARIANT_ALID_ROWS_FILENAME


def parse_canonical_alid(identifier: str) -> ParsedAlid | None:
    """Parse a canonical ALID-like identifier, returning None for aliases."""

    parts = str(identifier).split(":")
    if len(parts) != 4:
        return None
    chromosome, position_text, effect_allele, other_allele = parts
    try:
        position = int(position_text)
        if position <= 0:
            return None
        return ParsedAlid(
            chromosome=normalise_chromosome(chromosome),
            position=position,
            effect_allele=normalise_allele(effect_allele),
            other_allele=normalise_allele(other_allele),
        )
    except (TypeError, ValueError, VariantNormalisationError):
        return None


def write_variant_axis(
    store_path: str | Path,
    variants: list[CanonicalVariant],
    rsid_by_alid: Mapping[str, str],
    source_alids: Sequence[str | None] | None = None,
) -> None:
    """Write the Store Variant Table, row-offset sidecar, and tabix index.

    ``source_alids`` (optional, parallel to ``variants``) records the source-
    assembly canonical ALID that each row was lifted from — provenance back to
    the original dataset. Pass ``None`` (or per-row ``None``/``""``) for
    same-assembly builds or ambiguous (collided) rows; the column is written as
    ``"."`` there.
    """

    if source_alids is not None and len(source_alids) != len(variants):
        raise ValueError(
            f"source_alids has {len(source_alids)} entries but there are {len(variants)} variants"
        )
    store = Path(store_path)
    table_path = variant_table_path(store)
    offsets: list[int] = []
    with pysam.BGZFile(str(table_path), "w") as handle:  # type: ignore[call-arg]
        handle.write(VARIANT_HEADER.encode("utf-8"))
        for variant_index, variant in enumerate(variants):
            offsets.append(int(handle.tell()))
            rsid = rsid_by_alid.get(variant.alid) or "."
            source_alid = (source_alids[variant_index] if source_alids is not None else None) or "."
            line = (
                f"{variant.chromosome}\t{variant.position}\t{variant_index}\t"
                f"{variant.effect_allele}\t{variant.other_allele}\t{variant.alid}\t{rsid}\t"
                f"{source_alid}\n"
            )
            handle.write(line.encode("utf-8"))
    np.save(variant_offsets_path(store), np.asarray(offsets, dtype=np.uint64))
    pysam.tabix_index(
        str(table_path),
        seq_col=0,
        start_col=1,
        end_col=1,
        meta_char="#",
        zerobased=False,
        force=True,
    )
    # Write mmap'd ALID search index: two parallel arrays sorted by ALID bytes.
    alid_bytes = np.array([v.alid for v in variants], dtype=_ALID_DTYPE)
    row_indices = np.arange(len(variants), dtype="int32")
    sort_order = np.argsort(alid_bytes)
    np.save(variant_alid_bytes_path(store), alid_bytes[sort_order])
    np.save(variant_alid_rows_path(store), row_indices[sort_order])


class VariantAxis:
    """Read variants from a tabix-indexed Store Variant Table."""

    def __init__(self, store_path: str | Path, aliases: sqlite3.Connection | None = None):
        self.store_path = Path(store_path)
        self.aliases = aliases
        self.table_path = variant_table_path(self.store_path)
        self.tabix_path = variant_tabix_path(self.store_path)
        self.offsets_path = variant_offsets_path(self.store_path)
        self._offsets = np.load(self.offsets_path, mmap_mode="r")
        self._tabix = pysam.TabixFile(str(self.table_path))
        # mmap'd ALID search index — present on stores built after issue 029.
        # Falls back gracefully to tabix-per-call if absent (older stores).
        _bytes_path = variant_alid_bytes_path(self.store_path)
        _rows_path = variant_alid_rows_path(self.store_path)
        if _bytes_path.exists() and _rows_path.exists():
            self._alid_bytes: np.ndarray | None = np.load(_bytes_path, mmap_mode="r")
            self._alid_rows: np.ndarray | None = np.load(_rows_path, mmap_mode="r")
        else:
            self._alid_bytes = None
            self._alid_rows = None

    @property
    def n_variants(self) -> int:
        return int(len(self._offsets))

    def close(self) -> None:
        self._tabix.close()

    def by_identifier(self, identifier: str) -> VariantRecord | None:
        parsed = parse_canonical_alid(identifier)
        if parsed is not None:
            return self.by_alid(parsed)
        if self.aliases is None:
            return None
        row = self.aliases.execute(
            """
            SELECT variant_index
            FROM variant_aliases
            WHERE alias = ?
            ORDER BY variant_index
            LIMIT 1
            """,
            (identifier,),
        ).fetchone()
        if row is None:
            return None
        return self.by_index(int(row["variant_index"]))

    def indices_by_identifiers(self, identifiers: Sequence[str]) -> np.ndarray:
        """Resolve many identifiers to Store-local Variant Indices in one
        batched pass, without materialising a `VariantRecord` (and its
        per-call `variants.tsv.gz` BGZF open/seek, see `by_index`) for each
        one -- the fast path dense `lookup()` needs (issue #3).

        Canonical ALIDs are resolved with a single vectorised `searchsorted`
        over the mmap'd ALID index. Aliases/rsids have no vectorised index and
        fall back to `by_identifier()`'s per-identifier SQLite lookup, same as
        stores that predate the mmap ALID index (`self._alid_bytes is None`).
        Identifiers that don't resolve are dropped, matching `by_identifier()`'s
        None-on-miss semantics. The returned order is not guaranteed to match
        the input order.
        """
        parsed_pairs: list[tuple[str, ParsedAlid]] = []
        alias_identifiers: list[str] = []
        for identifier in identifiers:
            parsed = parse_canonical_alid(identifier)
            if parsed is None:
                alias_identifiers.append(identifier)
            else:
                query = (
                    f"{parsed.chromosome}:{parsed.position}:"
                    f"{parsed.effect_allele}:{parsed.other_allele}"
                )
                parsed_pairs.append((query, parsed))

        indices: list[int] = []
        if parsed_pairs:
            alid_bytes, alid_rows = self._alid_bytes, self._alid_rows
            if alid_bytes is not None and alid_rows is not None:
                queries = np.array([q for q, _ in parsed_pairs], dtype=_ALID_DTYPE)
                positions = np.searchsorted(alid_bytes, queries)
                in_bounds = positions < len(alid_bytes)
                hit = np.zeros(len(queries), dtype=bool)
                hit[in_bounds] = alid_bytes[positions[in_bounds]] == queries[in_bounds]
                indices.extend(int(row) for row in alid_rows[positions[hit]])
            else:
                for _, parsed in parsed_pairs:
                    record = self.by_alid(parsed)
                    if record is not None:
                        indices.append(record.variant_index)

        for identifier in alias_identifiers:
            record = self.by_identifier(identifier)
            if record is not None:
                indices.append(record.variant_index)

        return np.array(indices, dtype="int32")

    def by_alid(self, parsed: ParsedAlid) -> VariantRecord | None:
        if self._alid_bytes is not None:
            query = np.array(
                [f"{parsed.chromosome}:{parsed.position}:{parsed.effect_allele}:{parsed.other_allele}"],
                dtype=_ALID_DTYPE,
            )
            idx = int(np.searchsorted(self._alid_bytes, query[0]))
            if idx < len(self._alid_bytes) and self._alid_bytes[idx] == query[0]:
                return self.by_index(int(self._alid_rows[idx]))
            return None
        for record in self.range(parsed.chromosome, parsed.position, parsed.position):
            if (
                record.effect_allele == parsed.effect_allele
                and record.other_allele == parsed.other_allele
            ):
                return record
        return None

    def range(self, chromosome: str, start: int, end: int) -> list[VariantRecord]:
        chrom = normalise_chromosome(chromosome)
        try:
            lines = self._tabix.fetch(chrom, max(0, int(start) - 1), int(end))
        except ValueError:
            return []
        return [_parse_variant_line(line) for line in lines]

    def range_indices(self, chromosome: str, start: int, end: int) -> np.ndarray:
        """Return variant indices for a genomic range as int32 array, without object allocation."""
        chrom = normalise_chromosome(chromosome)
        try:
            lines = self._tabix.fetch(chrom, max(0, int(start) - 1), int(end))
        except ValueError:
            return np.empty(0, dtype="int32")
        indices = []
        for line in lines:
            # variant_index is the third tab-separated field (index 2)
            fields = line.split("\t", 3)
            indices.append(int(fields[2]))
        return np.array(indices, dtype="int32")

    def by_index(self, variant_index: int) -> VariantRecord | None:
        if variant_index < 0 or variant_index >= self.n_variants:
            return None
        with pysam.BGZFile(str(self.table_path), "r") as handle:  # type: ignore[call-arg]
            handle.seek(int(self._offsets[int(variant_index)]))
            line = handle.readline().decode("utf-8")
        record = _parse_variant_line(line)
        if record.variant_index != variant_index:
            raise ValueError(
                f"variant offset for row {variant_index} points to row {record.variant_index}"
            )
        return record

    def by_indices(self, indices: Iterable[int]) -> dict[int, VariantRecord]:
        """Resolve many Store-local Variant Indices to `VariantRecord`s.

        Adaptively picks between two strategies (issue #104 follow-up,
        measured on the UKB pilot store's ~9.85M-row Store Variant Table):
        a random-access `by_index()` seek costs ~300us/row regardless of
        whether the underlying BGZF handle is reused (the cost is the
        seek-and-decompress-from-the-nearest-block itself, not file-open
        overhead), while one sequential `all()` scan costs ~1-2us/row. The
        two cross over well under 1% of the table, so once a caller wants
        more than `_BULK_SCAN_FRACTION` of it, one full scan beats N random
        seeks by orders of magnitude -- the difference between a
        sub-second dense "query-analysis" resolution and one that never
        finishes in practice.
        """
        wanted = sorted(set(int(item) for item in indices))
        if not wanted:
            return {}
        if len(wanted) > self.n_variants * _BULK_SCAN_FRACTION:
            wanted_set = set(wanted)
            return {
                record.variant_index: record
                for record in self.all()
                if record.variant_index in wanted_set
            }
        records: dict[int, VariantRecord] = {}
        for index in wanted:
            if (record := self.by_index(index)) is not None:
                records[index] = record
        return records

    def all(self) -> list[VariantRecord]:
        with pysam.BGZFile(str(self.table_path), "r") as handle:  # type: ignore[call-arg]
            records: list[VariantRecord] = []
            for raw_line in handle:
                line = raw_line.decode("utf-8")
                if line.startswith("#"):
                    continue
                records.append(_parse_variant_line(line))
        return records


def _parse_variant_line(line: str) -> VariantRecord:
    fields = line.rstrip("\n").split("\t")
    # 7 fields = pre-source_alid stores; 8 = current (trailing source_alid column).
    if len(fields) == 7:
        chromosome, position, variant_index, effect_allele, other_allele, alid, rsid = fields
        source_alid = "."
    elif len(fields) == 8:
        (
            chromosome, position, variant_index, effect_allele, other_allele,
            alid, rsid, source_alid,
        ) = fields
    else:
        raise ValueError(f"variant row has {len(fields)} fields, expected 7 or 8")
    return VariantRecord(
        variant_index=int(variant_index),
        alid=alid,
        chromosome=chromosome,
        position=int(position),
        effect_allele=effect_allele,
        other_allele=other_allele,
        rsid=None if rsid in {"", "."} else rsid,
        source_alid=None if source_alid in {"", "."} else source_alid,
    )
