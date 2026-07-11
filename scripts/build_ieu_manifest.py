#!/usr/bin/env python3
"""Build a source manifest for the ieu-a / ieu-b / bbj-a collections.

Enumerates every study directory under the IGD tree for the given prefixes that
has a ``{id}.vcf.gz``, and reads its OpenGWAS ``{id}.json`` metadata for the
Reported Population, trait, and sample size. Emits a manifest TSV that the
``assign-ancestry`` pipeline consumes (build columns + ``reported_population``).

Usage:
  uv run python scripts/build_ieu_manifest.py \
      --igd /local-scratch/data/opengwas/igd \
      --prefixes ieu-a ieu-b bbj-a \
      --out /local-scratch/data/opengwas/ancestry_reference/ieu_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def _natural_key(study: str) -> tuple:
    m = re.match(r"([a-z]+-[a-z]+)-(\d+)", study)
    return (m.group(1), int(m.group(2))) if m else (study, 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--igd", type=Path, default=Path("/local-scratch/data/opengwas/igd"))
    ap.add_argument("--prefixes", nargs="+", default=["ieu-a", "ieu-b", "bbj-a"])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    studies: list[str] = []
    for prefix in args.prefixes:
        studies += [d.name for d in args.igd.glob(f"{prefix}-*") if d.is_dir()]
    studies = sorted(set(studies), key=_natural_key)

    rows = []
    n_no_vcf = n_no_json = 0
    for study in studies:
        vcf = args.igd / study / f"{study}.vcf.gz"
        if not vcf.exists():
            n_no_vcf += 1
            continue
        meta = args.igd / study / f"{study}.json"
        population = trait = ""
        n = 0
        if meta.exists():
            try:
                d = json.loads(meta.read_text())
                population = str(d.get("population", "") or "")
                trait = str(d.get("trait", study) or study)
                ss = d.get("sample_size")
                n = int(ss) if ss not in (None, "", "NULL") else 0
            except (json.JSONDecodeError, ValueError):
                n_no_json += 1
        else:
            n_no_json += 1
            trait = study
        rows.append((study, str(vcf), trait, n, population))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait_id", "file_path", "trait_name", "n", "reported_population"])
        w.writerows(rows)

    print(f"Wrote {len(rows)} studies → {args.out}")
    print(f"  (skipped {n_no_vcf} without a VCF; {n_no_json} without usable JSON metadata)")
    # Reported-population tally for a quick sanity check.
    from collections import Counter
    tally = Counter(r[4] or "(none)" for r in rows)
    for pop, c in tally.most_common():
        print(f"  {c:5d}  {pop}")


if __name__ == "__main__":
    main()
