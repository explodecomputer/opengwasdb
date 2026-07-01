"""BESD format reader — port of besdq.besd_reader for OpenGWASDB."""

from __future__ import annotations

import struct
from pathlib import Path
from typing import NamedTuple

import numpy as np


class SnpRecord(NamedTuple):
    row_idx: int
    chromosome: str
    snp_id: str
    bp: int
    a1: str | None
    a2: str | None
    freq: float | None


class ProbeRecord(NamedTuple):
    row_idx: int
    chromosome: str
    probe_id: str
    probe_bp: int
    gene: str | None
    orientation: str | None


def read_esi(esi_path: str | Path) -> list[SnpRecord]:
    snps: list[SnpRecord] = []
    with open(esi_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                snps.append(SnpRecord(
                    row_idx=len(snps),
                    chromosome=parts[0],
                    snp_id=parts[1],
                    bp=int(parts[3]),
                    a1=parts[4] if len(parts) > 4 else None,
                    a2=parts[5] if len(parts) > 5 else None,
                    freq=float(parts[6]) if len(parts) > 6 and parts[6] != "NA" else None,
                ))
            except (ValueError, IndexError):
                continue
    return snps


def read_epi(epi_path: str | Path) -> list[ProbeRecord]:
    probes: list[ProbeRecord] = []
    with open(epi_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                probes.append(ProbeRecord(
                    row_idx=len(probes),
                    chromosome=parts[0],
                    probe_id=parts[1],
                    probe_bp=int(parts[3]),
                    gene=parts[4] if len(parts) > 4 else None,
                    orientation=parts[5] if len(parts) > 5 else None,
                ))
            except (ValueError, IndexError):
                continue
    return probes


_MAGIC_SPARSE_3F = 0x40400000
_MAGIC_SPARSE_3 = 3
_RESERVED_UNITS = 16


class BESDReader:
    """Read SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F BESD format files.

    Stores the entire file in numpy arrays for fast per-probe slicing.
    For very large datasets (>100M associations) consider streaming — see issue 036.
    """

    def __init__(self, besd_path: str | Path, n_probes: int):
        self._besd_path = Path(besd_path)
        self._n_probes = n_probes
        self.format_type: str = ""
        self._cols: np.ndarray | None = None
        self._rowid: np.ndarray | None = None
        self._val: np.ndarray | None = None
        self._val_num: int = 0
        self._load()

    def _load(self) -> None:
        with open(self._besd_path, "rb") as fh:
            magic = struct.unpack("<I", fh.read(4))[0]
            if magic == _MAGIC_SPARSE_3F:
                self.format_type = "3F"
                self._parse(fh, skip_reserved=False)
            elif magic == _MAGIC_SPARSE_3:
                self.format_type = "3"
                self._parse(fh, skip_reserved=True)
            else:
                raise ValueError(f"Unsupported BESD magic: 0x{magic:08x}")

    def _parse(self, fh, skip_reserved: bool) -> None:
        if skip_reserved:
            fh.read((_RESERVED_UNITS - 1) * 4)
        self._val_num = struct.unpack("<Q", fh.read(8))[0]
        col_num = (self._n_probes << 1) + 1
        self._cols = np.frombuffer(fh.read(col_num * 8), dtype=np.int64).copy()
        self._rowid = np.frombuffer(fh.read(self._val_num * 4), dtype=np.uint32).copy()
        self._val = np.frombuffer(fh.read(self._val_num * 4), dtype=np.float32).copy()

    def get_probe_associations(
        self, probe_idx: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (snp_indices, betas, ses) as numpy arrays for one probe.

        snp_indices: uint32 ESI row indices
        betas: float32
        ses: float32
        """
        if (
            probe_idx >= self._n_probes
            or self._cols is None
            or self._rowid is None
            or self._val is None
        ):
            empty = np.empty(0, dtype=np.uint32)
            return empty, np.empty(0, dtype=np.float32), np.empty(0, dtype=np.float32)

        beta_start = int(self._cols[probe_idx << 1])
        se_start = int(self._cols[(probe_idx << 1) + 1])
        n = se_start - beta_start
        if n <= 0:
            empty = np.empty(0, dtype=np.uint32)
            return empty, np.empty(0, dtype=np.float32), np.empty(0, dtype=np.float32)

        snp_idx = self._rowid[beta_start:se_start]
        betas = self._val[beta_start:se_start]
        ses = self._val[se_start: se_start + n]
        return snp_idx, betas, ses
