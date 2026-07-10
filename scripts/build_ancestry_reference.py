#!/usr/bin/env python3
"""Build the hg38 Ancestry Reference Panel from Privé (2022)'s UK Biobank reference.

Downloads bigsnpr's ``ref_freqs.csv.gz`` (the ``snp_ancestry_summary`` reference:
allele frequencies for ~5.8M variants across 21 UK-Biobank-derived ancestry groups,
on **GRCh37**), lifts each variant to **GRCh38**, re-keys it to the canonical ALID
convention used by the stores (``chr:pos:A1:A2`` with ``A1 = min(a0, a1)``), and
orients every group's frequency to that canonical A1. Emits:

  * ``ref_freqs.hg38.tsv.gz`` — one row per lifted variant: ``alid, chromosome,
    position, effect_allele, other_allele, rsid`` + one canonical-A1 frequency
    column per group.
  * ``ancestry_groups.tsv`` — the fine-group → super-population map used to
    aggregate composition to the routing granularity (ADR 0028; fit-fine,
    label-coarse).

The output feeds the ancestry annotator (issues 062/063). It is a large static
artifact — write it to a data directory, not the git repo.

Usage:
  uv run python scripts/build_ancestry_reference.py \
      --out /local-scratch/data/opengwas/ancestry_reference/ref_freqs.hg38.tsv.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
import time
import urllib.request
from pathlib import Path

log = logging.getLogger("build_ancestry_reference")

# bigsnpr snp_ancestry_summary reference frequencies (figshare file 31620968).
# The ndownloader.figshare.com host 302-redirects straight to a signed S3 URL;
# figshare.com/ndownloader returns a 202 interstitial, so use this host.
DEFAULT_URL = "https://ndownloader.figshare.com/files/31620968"
_UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/124.0 Safari/537.36"

# The 21 fine groups (raw ref_freqs.csv.gz column order) → routing super-population.
# European groups aggregate to EUR (the only LD panel we have today); Middle East
# and North Africa are kept as their own coarse groups (no clean 1000G super-pop
# and no LD panel), so studies dominated by them fall to Unassigned — correct,
# since we cannot yet complete them. Edit this map as panels are added.
GROUP_TO_SUPERPOP: dict[str, str] = {
    "Africa (West)": "AFR",
    "Africa (South)": "AFR",
    "Africa (East)": "AFR",
    "Africa (North)": "NAF",
    "Middle East": "MID",
    "Ashkenazi": "EUR",
    "Italy": "EUR",
    "Europe (South East)": "EUR",
    "Europe (North East)": "EUR",
    "Finland": "EUR",
    "Scandinavia": "EUR",
    "United Kingdom": "EUR",
    "Ireland": "EUR",
    "Europe (South West)": "EUR",
    "South America": "AMR",
    "Sri Lanka": "SAS",
    "Pakistan": "SAS",
    "Bangladesh": "SAS",
    "Asia (East)": "EAS",
    "Japan": "EAS",
    "Philippines": "EAS",
}

_STANDARD_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def _bare_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        log.info("Using cached raw reference: %s (%.0f MB)", dest, dest.stat().st_size / 1e6)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    log.info("Downloading %s → %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req) as resp, open(tmp, "wb") as fh:  # noqa: S310
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)
    log.info("Downloaded %.0f MB", dest.stat().st_size / 1e6)


def _load_liftover(from_build: str, to_build: str, chain_file: str | None):
    try:
        from pyliftover import LiftOver  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pyliftover is required: pip install pyliftover") from exc
    if chain_file:
        return LiftOver(chain_file)
    return LiftOver(from_build, to_build)


def build_reference(
    raw_path: Path,
    out_path: Path,
    groups_out: Path,
    *,
    from_build: str,
    to_build: str,
    chain_file: str | None,
) -> None:
    lo = _load_liftover(from_build, to_build, chain_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_out = n_fail_lift = n_nonstd = n_bad = 0
    t0 = time.monotonic()

    with gzip.open(raw_path, "rt", newline="") as fin, \
            gzip.open(out_path, "wt", newline="") as fout:
        reader = csv.reader(fin)
        header = next(reader)
        # Fixed leading columns, then every remaining column is a group frequency.
        lead = {name: header.index(name) for name in ("chr", "pos", "a0", "a1", "rsid")}
        group_cols = [(i, name) for i, name in enumerate(header) if name not in lead]
        groups = [name for _i, name in group_cols]
        unmapped = [g for g in groups if g not in GROUP_TO_SUPERPOP]
        if unmapped:
            raise ValueError(f"reference groups without a super-pop mapping: {unmapped}")

        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            ["alid", "chromosome", "position", "effect_allele", "other_allele", "rsid", *groups]
        )

        i_chr, i_pos, i_a0, i_a1, i_rsid = (
            lead["chr"], lead["pos"], lead["a0"], lead["a1"], lead["rsid"]
        )
        for row in reader:
            n_in += 1
            try:
                chrom = _bare_chrom(row[i_chr])
                pos = int(row[i_pos])
                a0, a1 = row[i_a0], row[i_a1]
            except (IndexError, ValueError):
                n_bad += 1
                continue

            mapped = lo.convert_coordinate(f"chr{chrom}", pos - 1)  # 1-based → 0-based
            if not mapped:
                n_fail_lift += 1
                continue
            new_chrom = _bare_chrom(mapped[0][0])
            if new_chrom not in _STANDARD_CHROMS:
                n_nonstd += 1
                continue
            new_pos = int(mapped[0][1]) + 1  # 0-based → 1-based

            # Canonical A1 = min(a0, a1); orient every group's freq (freq of a1) to A1.
            allele_1, allele_2 = (a0, a1) if a0 < a1 else (a1, a0)
            flip = a1 != allele_1  # bigsnpr freq is of a1; flip if a1 is the non-canonical allele
            freqs = []
            for i, _name in group_cols:
                try:
                    f = float(row[i])
                except (IndexError, ValueError):
                    f = float("nan")
                freqs.append(f"{(1.0 - f) if flip else f:.6g}")

            alid = f"{new_chrom}:{new_pos}:{allele_1}:{allele_2}"
            writer.writerow(
                [alid, new_chrom, new_pos, allele_1, allele_2, row[i_rsid], *freqs]
            )
            n_out += 1

            if n_in % 500_000 == 0:
                log.info(
                    "  %d read, %d written (%.0fs elapsed)", n_in, n_out, time.monotonic() - t0
                )

    with open(groups_out, "w", newline="") as gf:
        w = csv.writer(gf, delimiter="\t")
        w.writerow(["group", "super_pop"])
        for g in groups:
            w.writerow([g, GROUP_TO_SUPERPOP[g]])

    log.info(
        "Done: %d in → %d out (%d liftover-failed, %d non-standard contig, %d malformed) in %.0fs",
        n_in, n_out, n_fail_lift, n_nonstd, n_bad, time.monotonic() - t0,
    )
    log.info("Reference: %s", out_path)
    log.info("Group map: %s", groups_out)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%H:%M:%S"
    )
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", default=DEFAULT_URL, help="ref_freqs.csv.gz download URL")
    ap.add_argument(
        "--raw-cache", type=Path, default=None,
        help="Where to cache the downloaded ref_freqs.csv.gz (default: alongside --out)",
    )
    ap.add_argument("--out", type=Path, required=True, help="Output hg38 reference .tsv.gz")
    ap.add_argument(
        "--groups-out", type=Path, default=None,
        help="Output fine→super-pop map TSV (default: ancestry_groups.tsv beside --out)",
    )
    ap.add_argument("--from-build", default="hg19")
    ap.add_argument("--to-build", default="hg38")
    ap.add_argument("--chain", default=None, help="Optional liftover chain file (overrides builds)")
    args = ap.parse_args()

    raw = args.raw_cache or args.out.with_name("ref_freqs.csv.gz")
    groups_out = args.groups_out or args.out.with_name("ancestry_groups.tsv")

    _download(args.url, raw)
    build_reference(
        raw, args.out, groups_out,
        from_build=args.from_build, to_build=args.to_build, chain_file=args.chain,
    )


if __name__ == "__main__":
    sys.exit(main())
