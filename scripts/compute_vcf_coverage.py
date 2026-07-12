#!/usr/bin/env python3
"""Compute per-study genome coverage from tabix indexes (no VCF scan).

For every study in a manifest, reads the ``.tbi`` index via ``bcftools index -s``
(per-contig record counts) and writes a coverage TSV: ``trait_id``,
``total_variants``, ``n_autosomes`` (of 1–22 present), ``frac_largest_chrom``
(share of variants on the most-populated autosome — high values flag single-region
studies). Feeds the genome-wide-store eligibility gate in
``opengwasdb.ancestry.routing``.

Usage:
  uv run python scripts/compute_vcf_coverage.py \
      --manifest /local-scratch/data/opengwas/ancestry_reference/ieu_manifest.tsv \
      --out      /local-scratch/data/opengwas/ancestry_reference/coverage.tsv \
      --workers 32
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def _coverage(trait_id: str, file_path: str) -> tuple[str, int, int, float]:
    try:
        out = subprocess.run(
            ["bcftools", "index", "-s", file_path],
            capture_output=True, text=True, timeout=180,
        ).stdout
    except (subprocess.SubprocessError, OSError):
        return (trait_id, 0, 0, 1.0)
    per: dict[str, int] = {}
    for line in out.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        contig = parts[0][3:] if parts[0].lower().startswith("chr") else parts[0]
        try:
            per[contig] = per.get(contig, 0) + int(parts[2])
        except ValueError:
            continue
    autosomes = [per.get(str(i), 0) for i in range(1, 23)]
    total = sum(per.values())
    n_auto = sum(1 for v in autosomes if v > 0)
    frac_largest = (max(autosomes) / total) if total else 1.0
    return (trait_id, total, n_auto, frac_largest)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    with open(args.manifest, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        results = list(ex.map(lambda r: _coverage(r["trait_id"], r["file_path"]), rows))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait_id", "total_variants", "n_autosomes", "frac_largest_chrom"])
        for tid, tot, na, fl in results:
            w.writerow([tid, tot, na, f"{fl:.4f}"])
    print(f"Wrote coverage for {len(results)} studies → {args.out}")


if __name__ == "__main__":
    main()
