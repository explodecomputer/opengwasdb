#!/usr/bin/env python3
"""Ancestry assignment on chr1 of real ieu-a GWAS-VCFs (development data).

Runs the real assignment pipeline (AF extraction with hg19→hg38 liftover → NNLS
mixture → multi-gate rule) on chromosome 1 of a handful of ieu-a studies, and
writes an Analysis Catalogue TSV for the QMD report. Uses chr1 only so the whole
sweep runs in seconds while still overlapping ~10^5 reference sites per study —
enough to clear the overlap gate. Reported Population is read from each study's
OpenGWAS ``{id}.json`` (it never routes; it only labels the τ/δ scatter).

Usage:
  uv run python benchmarks/ancestry_ieua_chr1.py \
      --reference /local-scratch/data/opengwas/ancestry_reference/ref_freqs.hg38.chr1.tsv.gz \
      --groups /local-scratch/data/opengwas/ancestry_reference/ancestry_groups.tsv \
      --igd /local-scratch/data/opengwas/igd \
      --out docs/benchmark-output/ancestry_ieua_chr1_catalogue.tsv
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

from opengwasdb.ancestry import (
    Gates,
    load_liftover,
    load_reference,
)
from opengwasdb.ancestry.catalogue import CatalogueRow, write_catalogue
from opengwasdb.ancestry.mixture import assign_from_vcf

log = logging.getLogger("ancestry_ieua_chr1")

# The three studies from the hybrid benchmark (all declared "Mixed"), plus a few
# genuinely non-European studies so the τ/δ scatter spans labelled ancestries.
STUDIES = [
    "ieu-a-2",     # BMI (GIANT) — declared Mixed, European in truth
    "ieu-a-7",     # Coronary heart disease (CARDIoGRAMplusC4D) — declared Mixed
    "ieu-a-300",   # LDL cholesterol (GLGC) — declared Mixed
    "ieu-a-1015",  # CRP (Japanese)
    "ieu-a-11",    # Crohn's disease (East Asian)
    "ieu-a-1099",  # Serum creatinine (African American)
    "ieu-a-1020",  # Percent emphysema (Hispanic)
    "bbj-a-2",     # BMI (Biobank Japan) — East Asian; same trait as ieu-a-2, contrast
]


def _reported(igd: Path, study: str) -> tuple[str, str]:
    """(reported_population, trait) from the study's OpenGWAS metadata JSON."""
    meta = igd / study / f"{study}.json"
    if not meta.exists():
        return "", study
    d = json.loads(meta.read_text())
    return str(d.get("population", "")), str(d.get("trait", study))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--groups", type=Path, required=True)
    ap.add_argument("--igd", type=Path, default=Path("/local-scratch/data/opengwas/igd"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--region", default="1", help="chromosome to assign on")
    ap.add_argument("--tau", type=float, default=0.90)
    ap.add_argument("--delta", type=float, default=0.20)
    ap.add_argument("--n-min", type=int, default=20_000)
    ap.add_argument("--residual-max", type=float, default=0.06)
    args = ap.parse_args()

    log.info("Loading chr%s ancestry reference…", args.region)
    reference = load_reference(args.reference, args.groups, maf_floor=0.01)
    log.info("Reference: %d variants, %d groups, super-pops %s",
             reference.n_variants, len(reference.groups), reference.superpops)
    gates = Gates(tau=args.tau, delta=args.delta, n_min=args.n_min, residual_max=args.residual_max)
    lo = load_liftover("hg19", "hg38")

    rows: list[CatalogueRow] = []
    for study in STUDIES:
        vcf = args.igd / study / f"{study}.vcf.gz"
        if not vcf.exists():
            log.warning("skipping %s (no VCF at %s)", study, vcf)
            continue
        reported, trait = _reported(args.igd, study)
        t0 = time.monotonic()
        a = assign_from_vcf(vcf, reference, gates, region=args.region, liftover=lo)
        log.info(
            "%-12s reported=%-18s assigned=%-10s dom=%s:%.3f margin=%.3f overlap=%d "
            "resid=%.4f gate=%s (%.1fs)",
            study, reported, a.assigned_ancestry, a.dominant_superpop, a.dominant_proportion,
            a.runner_up_margin, a.af_overlap, a.residual, a.gate_reason, time.monotonic() - t0,
        )
        rows.append(CatalogueRow(study, str(vcf), trait, 0, reported, a))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_catalogue(
        args.out, rows, reference.superpops,
        catalogue_version=f"ieua-chr{args.region}-dev",
        ancestry_reference_version="prive2022-hg38-chr" + args.region,
        gates=gates,
    )
    log.info("Wrote Catalogue: %s (%d studies)", args.out, len(rows))


if __name__ == "__main__":
    main()
