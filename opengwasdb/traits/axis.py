"""Tabix-backed Store Traits Axis."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

TRAITS_TABLE_FILENAME = "traits.tsv.gz"
TRAITS_TABIX_FILENAME = "traits.tsv.gz.tbi"
TRAITS_AXIS_FORMAT = "tabix_tsv_v1"
TRAITS_HEADER = (
    "#probe_chr\tprobe_bp\tanalysis_index\tanalysis_id\t"
    "probe_id\tn\tgene_id\tgene_name\ttissue\tcontext\n"
)


@dataclass(frozen=True)
class TraitRecord(Mapping[str, Any]):
    """One Store Traits Table row."""

    analysis_index: int
    analysis_id: str
    probe_id: str
    n: int | None
    probe_chr: str | None
    probe_bp: int | None
    gene_id: str | None
    gene_name: str | None
    tissue: str | None
    context: str | None

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __iter__(self) -> Iterator[str]:
        return iter((
            "analysis_index", "analysis_id", "probe_id", "n",
            "probe_chr", "probe_bp", "gene_id", "gene_name", "tissue", "context",
        ))

    def __len__(self) -> int:
        return 10

    def has_position(self) -> bool:
        return self.probe_chr is not None and self.probe_bp is not None


def traits_table_path(store_path: str | Path) -> Path:
    return Path(store_path) / TRAITS_TABLE_FILENAME


def traits_tabix_path(store_path: str | Path) -> Path:
    return Path(store_path) / TRAITS_TABIX_FILENAME


def _fmt(value: Any) -> str:
    return "NA" if value is None else str(value)


def _parse_optional_str(s: str) -> str | None:
    return None if s == "NA" else s


def _parse_optional_int(s: str) -> int | None:
    return None if s == "NA" else int(s)


def write_traits_axis(
    store_path: str | Path,
    records: list[TraitRecord],
) -> None:
    """Write traits.tsv.gz and tabix index.

    Records with probe_chr/probe_bp are sorted by position and tabix-indexed.
    Records without position are appended at the end (no tabix index created
    for the whole file if any record lacks coordinates).
    """
    store = Path(store_path)
    table_path = traits_table_path(store)

    has_positions = all(r.has_position() for r in records)

    if has_positions and records:
        sorted_records = sorted(records, key=lambda r: (_chr_sort_key(r.probe_chr), r.probe_bp))
    else:
        sorted_records = list(records)

    with pysam.BGZFile(str(table_path), "w") as handle:  # type: ignore[call-arg]
        handle.write(TRAITS_HEADER.encode("utf-8"))
        for rec in sorted_records:
            line = (
                f"{_fmt(rec.probe_chr)}\t{_fmt(rec.probe_bp)}\t{rec.analysis_index}\t"
                f"{rec.analysis_id}\t{rec.probe_id}\t{_fmt(rec.n)}\t"
                f"{_fmt(rec.gene_id)}\t{_fmt(rec.gene_name)}\t"
                f"{_fmt(rec.tissue)}\t{_fmt(rec.context)}\n"
            )
            handle.write(line.encode("utf-8"))

    if has_positions and records:
        pysam.tabix_index(
            str(table_path),
            seq_col=0,
            start_col=1,
            end_col=1,
            meta_char="#",
            zerobased=False,
            force=True,
        )


def _chr_sort_key(chrom: str | None) -> tuple[int, str]:
    if chrom is None:
        return (999, "")
    c = chrom.lstrip("chr")
    try:
        return (int(c), "")
    except ValueError:
        return (998, c)


def _parse_record(line: str) -> TraitRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 10:
        raise ValueError(f"traits row has {len(fields)} fields, expected 10: {line!r}")
    probe_chr, probe_bp, analysis_index, analysis_id, probe_id, n, gene_id, gene_name, tissue, context = fields
    return TraitRecord(
        probe_chr=_parse_optional_str(probe_chr),
        probe_bp=_parse_optional_int(probe_bp),
        analysis_index=int(analysis_index),
        analysis_id=analysis_id,
        probe_id=probe_id,
        n=_parse_optional_int(n),
        gene_id=_parse_optional_str(gene_id),
        gene_name=_parse_optional_str(gene_name),
        tissue=_parse_optional_str(tissue),
        context=_parse_optional_str(context),
    )


class TraitsAxisReader:
    """Read analyses from a tabix-indexed Store Traits Table."""

    def __init__(self, store_path: str | Path):
        self.store_path = Path(store_path)
        self._table_path = traits_table_path(self.store_path)
        self._tabix_path = traits_tabix_path(self.store_path)
        self._tabix: pysam.TabixFile | None = None
        if self._tabix_path.exists():
            self._tabix = pysam.TabixFile(str(self._table_path))

    def close(self) -> None:
        if self._tabix is not None:
            self._tabix.close()

    def __enter__(self) -> TraitsAxisReader:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def range(self, chromosome: str, start: int, end: int) -> list[TraitRecord]:
        """Return all analyses with probe position in [start, end] on chromosome."""
        if self._tabix is None:
            return []
        try:
            lines = self._tabix.fetch(chromosome, max(0, start - 1), end)
            return [_parse_record(line) for line in lines]
        except ValueError:
            return []

    def by_probe_id(self, probe_id: str) -> list[TraitRecord]:
        """Return all analyses matching probe_id (linear scan)."""
        results: list[TraitRecord] = []
        with pysam.BGZFile(str(self._table_path), "r") as handle:  # type: ignore[call-arg]
            for raw in handle:
                line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
                if line.startswith("#"):
                    continue
                rec = _parse_record(line)
                if rec.probe_id == probe_id:
                    results.append(rec)
        return results

    def all(self) -> Iterator[TraitRecord]:
        """Iterate over all analysis records."""
        with pysam.BGZFile(str(self._table_path), "r") as handle:  # type: ignore[call-arg]
            for raw in handle:
                line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
                if line.startswith("#"):
                    continue
                yield _parse_record(line)
