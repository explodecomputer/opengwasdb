#!/usr/bin/env python3
"""Genome-wide ancestry assignment across the ieu-a / ieu-b / bbj-a collections.

Study VCFs are GRCh37, so this assigns against the **hg19** reference
(``ref_freqs.hg19.tsv.gz``, built with ``--no-liftover``) — study and reference
share coordinates, so no per-variant liftover is needed. Extraction streams each
study once and filters to reference sites by O(1) membership (no bcftools ``-R``
region reload per study). Studies are processed across a fork pool (the ~1 GB
reference is loaded once in the parent and fork-inherited); results are gathered
in manifest order so the Catalogue is worker-count independent.

Usage:
  uv run python benchmarks/ancestry_ieu_batch.py \
      --reference /local-scratch/data/opengwas/ancestry_reference/ref_freqs.hg19.tsv.gz \
      --groups    /local-scratch/data/opengwas/ancestry_reference/ancestry_groups.tsv \
      --manifest  /local-scratch/data/opengwas/ancestry_reference/ieu_manifest.tsv \
      --out       docs/benchmark-output/ancestry_ieu_catalogue.tsv \
      --workers 24
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from opengwasdb.ancestry import Gates, load_reference
from opengwasdb.ancestry.catalogue import CatalogueRow, write_catalogue
from opengwasdb.ancestry.extract import extract_af_at_sites
from opengwasdb.ancestry.mixture import assign_ancestry
from opengwasdb.ancestry.pipeline import read_source_manifest

log = logging.getLogger("ancestry_ieu_batch")

_REF = None  # fork-inherited AncestryReference
_GATES: Gates | None = None


def _annotate(file_path: str):
    assert _REF is not None and _GATES is not None
    study_af = extract_af_at_sites(file_path, _REF.index)  # stream; no -R, no liftover
    return assign_ancestry(study_af, _REF, _GATES)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%H:%M:%S"
    )
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--groups", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--maf-floor", type=float, default=0.01)
    ap.add_argument("--tau", type=float, default=0.50)
    ap.add_argument("--delta", type=float, default=0.20)
    ap.add_argument("--n-min", type=int, default=5_000)
    ap.add_argument("--residual-max", type=float, default=0.06)
    ap.add_argument("--catalogue-version", default="ieu-abj-genome-v1")
    args = ap.parse_args()

    global _REF, _GATES
    t0 = time.monotonic()
    log.info("Loading hg19 ancestry reference…")
    _REF = load_reference(args.reference, args.groups, maf_floor=args.maf_floor)
    log.info("Reference: %d variants, super-pops %s (%.0fs)",
             _REF.n_variants, _REF.superpops, time.monotonic() - t0)
    _GATES = Gates(tau=args.tau, delta=args.delta, n_min=args.n_min, residual_max=args.residual_max)

    source_rows = read_source_manifest(args.manifest)
    n = len(source_rows)
    log.info("Assigning %d studies with %d workers…", n, args.workers)

    results: dict[int, object] = {}
    t1 = time.monotonic()
    fork_ctx = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=fork_ctx) as pool:
        futures = {pool.submit(_annotate, r.file_path): i for i, r in enumerate(source_rows)}
        done = 0
        from concurrent.futures import as_completed
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                results[i] = fut.result()
            except Exception as exc:  # noqa: BLE001 — one bad VCF must not kill the batch
                log.warning("study %s failed: %s", source_rows[i].trait_id, exc)
                results[i] = None
            done += 1
            if done % 50 == 0 or done == n:
                rate = done / (time.monotonic() - t1)
                eta = (n - done) / rate if rate else 0
                log.info("  %d/%d done (%.1f studies/s, ETA %.0fs)", done, n, rate, eta)

    rows = []
    for i, src in enumerate(source_rows):
        a = results.get(i)
        if a is None:
            continue
        rows.append(CatalogueRow(src.trait_id, src.file_path, src.trait_name, src.n,
                                 src.reported_population, a))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_catalogue(
        args.out, rows, _REF.superpops,
        catalogue_version=args.catalogue_version,
        ancestry_reference_version="prive2022-hg19",
        gates=_GATES,
    )

    assigned = Counter(r.assignment.assigned_ancestry or "Unassigned" for r in rows)
    log.info("Wrote Catalogue: %s (%d studies, %.0fs total)",
             args.out, len(rows), time.monotonic() - t0)
    log.info("Assigned-ancestry tally: %s", dict(assigned.most_common()))


if __name__ == "__main__":
    main()
