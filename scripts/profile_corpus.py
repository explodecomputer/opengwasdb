#!/usr/bin/env python3
"""Profile a GWAS-VCF corpus into an ingestion-scenario inventory.

Every GWAS-VCF declares the same FORMAT block, so the scenarios come from what is
actually *populated* plus the header labels. For each study this samples a few
thousand records to detect which FORMAT fields carry data (ES, SE, EZ, LP, AF, SS,
NC) and reads the ##SAMPLE header (StudyType, TotalCases/TotalControls) and contig
build. Emits one row per study with the raw signals; downstream code derives the
scenario axes (effect-stat completeness, scale, AF state, sample-size kind, build).

Usage:
  uv run python scripts/profile_corpus.py \
      --igd /local-scratch/data/opengwas/igd \
      --prefixes ieu-a ieu-b bbj-a ebi-a \
      --out /local-scratch/data/opengwas/ancestry_reference/corpus_profile.tsv \
      --workers 48 --sample 3000
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_STUDYTYPE = re.compile(r"StudyType=([^,>\s]+)")
_TOTCASES = re.compile(r"TotalCases=([0-9.]+)")
_TOTCTRL = re.compile(r"TotalControls=([0-9.]+)")
_CONTIG1 = re.compile(r"##contig=<ID=(?:chr)?1,length=(\d+)")
FORMAT_FIELDS = ["ES", "SE", "EZ", "LP", "AF", "SS", "NC"]
# GRCh37 vs GRCh38 chr1 length — distinguishes build when no assembly tag is present.
_HG19_CHR1, _HG38_CHR1 = 249250621, 248956422


def _profile(trait_id: str, vcf: str, sample: int) -> dict:
    row: dict = {"trait_id": trait_id, "study_type": "", "total_cases": "",
                 "total_controls": "", "build": "unknown"}
    for f in FORMAT_FIELDS:
        row[f"has_{f}"] = False
    try:
        hdr = subprocess.run(["bcftools", "view", "-h", vcf],
                             capture_output=True, text=True, timeout=120).stdout
    except (subprocess.SubprocessError, OSError):
        return row
    if (m := _STUDYTYPE.search(hdr)):
        row["study_type"] = m.group(1)
    if (m := _TOTCASES.search(hdr)):
        row["total_cases"] = m.group(1)
    if (m := _TOTCTRL.search(hdr)):
        row["total_controls"] = m.group(1)
    if (m := _CONTIG1.search(hdr)):
        length = int(m.group(1))
        row["build"] = "hg19" if length == _HG19_CHR1 else ("hg38" if length == _HG38_CHR1 else "other")

    fmt = "\t".join(f"[%{f}]" for f in FORMAT_FIELDS)
    try:
        proc = subprocess.Popen(["bcftools", "query", "-f", fmt + "\n", vcf],
                                stdout=subprocess.PIPE, text=True)
        seen = {f: False for f in FORMAT_FIELDS}
        n = 0
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == len(FORMAT_FIELDS):
                for f, v in zip(FORMAT_FIELDS, parts, strict=True):
                    if not seen[f] and v not in (".", "", "nan"):
                        seen[f] = True
            n += 1
            if n >= sample or all(seen.values()):
                break
        proc.stdout.close()
        proc.terminate()
        for f in FORMAT_FIELDS:
            row[f"has_{f}"] = seen[f]
    except (subprocess.SubprocessError, OSError):
        pass
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--igd", type=Path, default=Path("/local-scratch/data/opengwas/igd"))
    ap.add_argument("--prefixes", nargs="+", default=["ieu-a", "ieu-b", "bbj-a", "ebi-a"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--sample", type=int, default=3000)
    args = ap.parse_args()

    studies = []
    for prefix in args.prefixes:
        for d in args.igd.glob(f"{prefix}-*"):
            vcf = d / f"{d.name}.vcf.gz"
            if vcf.exists():
                studies.append((d.name, str(vcf)))
    studies.sort()
    print(f"Profiling {len(studies)} studies…")

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(lambda s: _profile(s[0], s[1], args.sample), studies))

    fields = ["trait_id", "collection", "build", "study_type", "total_cases",
              "total_controls"] + [f"has_{f}" for f in FORMAT_FIELDS]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            r["collection"] = r["trait_id"].rsplit("-", 1)[0]
            w.writerow({k: r.get(k, "") for k in fields})
    print(f"Wrote {len(rows)} rows → {args.out}")


if __name__ == "__main__":
    main()
