"""Ancestry Reference Panel loader (ADR 0028; fit-fine, label-coarse).

Reads the static Ancestry Reference Panel produced by
``scripts/build_ancestry_reference.py`` — per-fine-group allele frequencies over
~5.8M variants, keyed to canonical ALID and oriented to A1 — plus the
fine→super-population grouping map. The panel is fit at fine granularity and the
mixture proportions are aggregated to super-populations for the routing label.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Fixed leading columns of ``ref_freqs.hg38.tsv.gz`` before the per-group columns.
_LEAD_COLUMNS = ("alid", "chromosome", "position", "effect_allele", "other_allele", "rsid")


@dataclass(frozen=True)
class AncestryReference:
    """Fine-group allele-frequency reference with a super-population map.

    ``freqs`` is ``(n_variants, n_groups)`` — the frequency of each variant's
    canonical A1 in each fine group. ``group_superpop_index`` gives, per group
    column, the row of ``superpops`` it aggregates to (fit-fine, label-coarse).
    """

    alids: np.ndarray  # (V,) str, canonical ALIDs
    freqs: np.ndarray  # (V, G) float64, freq of canonical A1 per fine group
    groups: list[str]  # (G,) fine-group names, column order of ``freqs``
    superpops: list[str]  # (S,) sorted unique super-populations
    group_superpop_index: np.ndarray  # (G,) int, group col -> superpop row
    group_to_superpop: dict[str, str]
    index: dict[str, int]  # alid -> row in ``freqs``

    @property
    def n_variants(self) -> int:
        return int(self.freqs.shape[0])

    def aggregate(self, fine_proportions: np.ndarray) -> np.ndarray:
        """Sum fine-group proportions into super-population proportions."""
        out = np.zeros(len(self.superpops), dtype=np.float64)
        np.add.at(out, self.group_superpop_index, fine_proportions)
        return out


def load_group_map(groups_path: str | Path) -> dict[str, str]:
    """Read the ``group\tsuper_pop`` fine→super-population map."""
    mapping: dict[str, str] = {}
    with open(groups_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["group"]] = row["super_pop"]
    if not mapping:
        raise ValueError(f"empty group map: {groups_path}")
    return mapping


def load_reference(
    freqs_path: str | Path,
    groups_path: str | Path,
    *,
    maf_floor: float = 0.01,
) -> AncestryReference:
    """Load the Ancestry Reference Panel and its fine→super-population map.

    ``freqs_path`` is the (optionally gzipped) TSV written by
    ``scripts/build_ancestry_reference.py``. Variants whose mean reference MAF is
    below ``maf_floor`` are dropped as uninformative for the mixture fit (set
    ``maf_floor=0`` to keep everything, e.g. on tiny synthetic panels).
    """
    group_to_superpop = load_group_map(groups_path)

    opener = gzip.open if str(freqs_path).endswith(".gz") else open
    alids: list[str] = []
    freq_rows: list[list[float]] = []
    with opener(freqs_path, "rt", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        for name in _LEAD_COLUMNS:
            if name not in header:
                raise ValueError(f"reference missing column {name!r}: {freqs_path}")
        i_alid = header.index("alid")
        group_cols = [(i, name) for i, name in enumerate(header) if name not in _LEAD_COLUMNS]
        if not group_cols:
            raise ValueError(f"reference has no group frequency columns: {freqs_path}")
        groups = [name for _i, name in group_cols]
        missing = [g for g in groups if g not in group_to_superpop]
        if missing:
            raise ValueError(f"reference groups absent from group map: {missing}")

        col_idx = [i for i, _name in group_cols]
        for row in reader:
            alids.append(row[i_alid])
            freq_rows.append([_to_float(row[i]) for i in col_idx])

    freqs = np.asarray(freq_rows, dtype=np.float64)
    alids_arr = np.asarray(alids, dtype=object)

    if maf_floor > 0 and freqs.size:
        mean_f = np.nanmean(freqs, axis=1)
        maf = np.minimum(mean_f, 1.0 - mean_f)
        keep = maf >= maf_floor
        freqs = freqs[keep]
        alids_arr = alids_arr[keep]

    superpops = sorted({group_to_superpop[g] for g in groups})
    sp_row = {sp: i for i, sp in enumerate(superpops)}
    group_superpop_index = np.asarray(
        [sp_row[group_to_superpop[g]] for g in groups], dtype=np.int64
    )
    index = {alid: i for i, alid in enumerate(alids_arr.tolist())}

    return AncestryReference(
        alids=alids_arr,
        freqs=freqs,
        groups=groups,
        superpops=superpops,
        group_superpop_index=group_superpop_index,
        group_to_superpop=group_to_superpop,
        index=index,
    )


def _to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")
