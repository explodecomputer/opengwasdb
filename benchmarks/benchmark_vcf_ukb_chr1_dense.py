#!/usr/bin/env python3
"""Build and benchmark a Dense Observed-Only Store from the UKB chr1 VCF dataset.

Uses the 100-trait chr1 GWAS-VCF dataset at
  /home/gh13047/repo/besdq/data/vcf-ukb/
which contains bgzipped, GRCh37/hg19 GWAS-VCFs in GWAS-SSF format.

The store is built via the two-pass ``build_dense_from_vcf_manifest`` pipeline
with inline hg19 → hg38 liftover.  Five query patterns are benchmarked (10
repetitions each) and compared against the besdq zarr baseline in
  /home/gh13047/repo/besdq/data/ukb-chr1_zarr_benchmark.json

Outputs:
  docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.json

Usage:
  python benchmarks/benchmark_vcf_ukb_chr1_dense.py [--rebuild] [--reps N]
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from pathlib import Path

import numpy as np

from opengwasdb.index import connect
from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.query import query_store
from opengwasdb.validation import validate_store

DEFAULT_MANIFEST = Path("/home/gh13047/repo/besdq/data/vcf-ukb/manifest.tsv")
DEFAULT_STORE = Path("/home/gh13047/repo/opengwasdb/docs/benchmark-output/vcf_ukb_chr1_dense.opengwasdb")
DEFAULT_OUTPUT = Path("/home/gh13047/repo/opengwasdb/docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.json")
DEFAULT_QMD = Path("/home/gh13047/repo/opengwasdb/docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.qmd")
BESDQ_BASELINE = Path("/home/gh13047/repo/besdq/data/ukb-chr1_zarr_benchmark.json")
ROW_BASELINE = Path("/home/gh13047/repo/opengwasdb/docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.json")
SLOWDOWN_FLAG_THRESHOLD = 2.0  # flag patterns >2× slower than besdq zstd_bitshuffle

# Match besdq's benchmark region exactly (dense_05_query_benchmark.py REGION constant)
BESDQ_REGION = ("1", 100_000_000, 101_000_000)

RNG = np.random.default_rng(42)


def main() -> None:
    args = _parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    build_seconds = None
    liftover_failure_count = 0

    if args.rebuild and args.store.exists():
        shutil.rmtree(args.store)

    if args.rebuild or not args.store.exists():
        print(f"Building store from {args.manifest} ...")
        t0 = time.perf_counter()
        result = build_dense_from_vcf_manifest(
            args.manifest,
            args.store,
            store_id="ukb-chr1-vcf",
            release_id="dense-observed-vcf-v1",
        )
        build_seconds = time.perf_counter() - t0
        print(
            f"Build complete in {build_seconds:.1f}s: "
            f"{result.n_variants:,} variants × {result.n_analyses} analyses"
        )

    print("Validating store ...")
    validation = validate_store(args.store)
    if not validation.ok:
        raise SystemExit(f"Store is invalid: {validation.errors}")
    print("Validation passed.")

    print(f"Running benchmark ({args.reps} reps per pattern) ...")
    results = _run_benchmark(args.store, args.reps, build_seconds, liftover_failure_count)

    besdq_baseline = _load_besdq_baseline(BESDQ_BASELINE)
    row_baseline = _load_besdq_baseline(args.row_baseline) if args.row_baseline else None
    if besdq_baseline:
        _annotate_slowdowns(results["timings"], besdq_baseline, row_baseline)

    args.output.write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Results written to {args.output}")

    # Copy besdq baseline alongside our output so the comparison QMD is self-contained.
    besdq_copy = args.output.parent / "besdq_ukb_chr1_benchmark.json"
    if BESDQ_BASELINE.exists():
        shutil.copy2(BESDQ_BASELINE, besdq_copy)
        print(f"besdq baseline copied to {besdq_copy}")

    qmd_path = args.qmd
    _write_qmd(results, qmd_path, args.output, besdq_baseline, row_baseline, args.reps)
    print(f"QMD written to {qmd_path}")

    _print_summary(results, besdq_baseline)


def _run_benchmark(
    store_path: Path,
    n_reps: int,
    build_seconds: float | None,
    liftover_failure_count: int,
) -> dict:
    query = query_store(store_path)
    selection = _choose_queries(store_path)

    query_specs = {
        "regional": lambda: query.range(
            selection["region_chrom"], selection["region_start"], selection["region_end"]
        ),
        "phewas": lambda: query.phewas(selection["phewas_alid"]),
        "tophits": lambda: query.top_hits(threshold=selection["top_hit_threshold"]),
        "bulk": lambda: query.analysis(selection["bulk_analysis_id"]),
        "random_lookup": lambda: query.lookup(
            selection["random_alids"], selection["random_analysis_ids"]
        ),
    }

    timings = []
    for name, fn in query_specs.items():
        median_ms, p95_ms, result_count = _bench(fn, n_reps)
        timings.append(
            {
                "query": name,
                "median_ms": round(median_ms, 3),
                "p95_ms": round(p95_ms, 3),
                "result_count": result_count,
            }
        )
        print(f"  {name}: median={median_ms:.2f}ms  p95={p95_ms:.2f}ms  rows={result_count:,}")

    with connect(store_path / "index.sqlite") as conn:
        n_variants = conn.execute(
            "SELECT COUNT(*) FROM variants"
        ).fetchone()[0]
        n_analyses = conn.execute(
            "SELECT COUNT(*) FROM analyses"
        ).fetchone()[0]

    store_bytes = sum(f.stat().st_size for f in store_path.rglob("*") if f.is_file())

    return {
        "dataset": {
            "n_variants": n_variants,
            "n_analyses": n_analyses,
        },
        "build": {
            "store_path": str(store_path),
            "build_seconds": round(build_seconds, 2) if build_seconds is not None else None,
            "liftover_failure_count": liftover_failure_count,
        },
        "storage": {
            "store_bytes": store_bytes,
            "store_mb": round(store_bytes / 1_000_000, 2),
        },
        "selection": selection,
        "timings": timings,
    }


def _choose_queries(store_path: Path) -> dict:
    with connect(store_path / "index.sqlite") as conn:
        variants = conn.execute(
            "SELECT variant_index, alid, chromosome, position FROM variants ORDER BY variant_index"
        ).fetchall()
        analyses = conn.execute(
            "SELECT analysis_index, analysis_id FROM analyses ORDER BY analysis_index"
        ).fetchall()

    import zarr

    root = zarr.open_group(str(store_path / "data.zarr"), mode="r")
    z = root["z"][:].astype("float32")
    finite = np.isfinite(z)

    phewas_row = int(np.argmax(finite.sum(axis=1)))
    bulk_col = int(np.argmax(finite.sum(axis=0)))

    # Use the same window as besdq's benchmark for a fair comparison.
    region_chrom, region_start, region_end = BESDQ_REGION

    # Top-hit threshold: pick the loosest threshold that returns results
    top_hit_threshold = 5e-8
    for threshold in (5e-8, 5e-6, 5e-4):
        key = f"p_{threshold:.0e}".replace("-", "_").replace("+", "")
        if "top_hits" in root and key in root["top_hits"] and len(root["top_hits"][key]["z"]) > 0:
            top_hit_threshold = threshold
            break

    n_random_variants = min(100, len(variants))
    n_random_analyses = min(10, len(analyses))
    random_rows = sorted(RNG.choice(len(variants), n_random_variants, replace=False))
    random_cols = sorted(RNG.choice(len(analyses), n_random_analyses, replace=False))

    return {
        "region_chrom": region_chrom,
        "region_start": region_start,
        "region_end": region_end,
        "phewas_alid": str(variants[phewas_row]["alid"]),
        "bulk_analysis_id": str(analyses[bulk_col]["analysis_id"]),
        "top_hit_threshold": top_hit_threshold,
        "random_alids": [str(variants[r]["alid"]) for r in random_rows],
        "random_analysis_ids": [str(analyses[c]["analysis_id"]) for c in random_cols],
    }


def _bench(fn, n_reps: int) -> tuple[float, float, int]:
    result = fn()  # warm-up
    times = []
    for _ in range(n_reps):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return (
        float(np.median(times) * 1000),
        float(np.percentile(times, 95) * 1000),
        int(len(result["z"])),
    )


def _load_besdq_baseline(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _annotate_slowdowns(
    timings: list[dict],
    baseline: dict,
    row_baseline: dict | None = None,
) -> None:
    zstd = baseline.get("zstd_bitshuffle", {})
    row_by_query = {}
    if row_baseline:
        row_by_query = {t["query"]: t for t in row_baseline.get("timings", [])}
    for row in timings:
        name = row["query"]
        besdq_median = zstd.get(name, {}).get("median_ms")
        if besdq_median is not None and besdq_median > 0:
            ratio = row["median_ms"] / besdq_median
            row["besdq_zstd_median_ms"] = besdq_median
            row["ratio_vs_besdq"] = round(ratio, 2)
            if ratio > SLOWDOWN_FLAG_THRESHOLD:
                row["notes"] = (
                    f"WARNING: {ratio:.1f}× slower than besdq zstd_bitshuffle baseline."
                )
        if name in row_by_query:
            old_median = row_by_query[name].get("median_ms")
            if old_median and old_median > 0:
                speedup = old_median / row["median_ms"]
                row["row_api_median_ms"] = old_median
                row["speedup_vs_row_api"] = round(speedup, 2)


def _print_summary(results: dict, baseline: dict | None) -> None:
    print("\n--- Summary ---")
    print(f"  Variants: {results['dataset']['n_variants']:,}")
    print(f"  Analyses: {results['dataset']['n_analyses']}")
    if results["build"]["build_seconds"] is not None:
        print(f"  Build time: {results['build']['build_seconds']:.1f}s")
    print(f"  Store size: {results['storage']['store_mb']:.1f} MB")
    print("")
    print(f"  {'Query':<15} {'Median ms':>10} {'p95 ms':>10} {'vs besdq':>10}")
    for row in results["timings"]:
        ratio = row.get("ratio_vs_besdq", "")
        ratio_str = f"{ratio:.2f}×" if ratio != "" else "n/a"
        flag = " ⚠" if "notes" in row else ""
        print(f"  {row['query']:<15} {row['median_ms']:>10.2f} {row['p95_ms']:>10.2f} {ratio_str:>10}{flag}")


def _write_qmd(
    results: dict,
    qmd_path: Path,
    json_path: Path,
    besdq_baseline: dict | None,
    row_baseline: dict | None,
    n_reps: int,
) -> None:
    build = results["build"]
    dataset = results["dataset"]
    storage = results["storage"]
    timings = results["timings"]

    build_line = (
        f"{build['build_seconds']:.1f}s" if build["build_seconds"] is not None
        else "not rebuilt — existing store reused"
    )
    liftover_line = str(build.get("liftover_failure_count", "n/a"))

    timing_rows = []
    for t in timings:
        row_str = f"{t.get('row_api_median_ms', 'n/a')}"
        speedup_str = (
            f"{t['speedup_vs_row_api']:.1f}×" if "speedup_vs_row_api" in t else "n/a"
        )
        besdq_str = f"{t.get('besdq_zstd_median_ms', 'n/a')}"
        ratio_str = (
            f"{t['ratio_vs_besdq']:.1f}×" if "ratio_vs_besdq" in t else "n/a"
        )
        flag = " ⚠" if "notes" in t else ""
        timing_rows.append(
            f"| {t['query']} | {t['median_ms']:.1f} | {t['p95_ms']:.1f}"
            f" | {t['result_count']:,} | {row_str} | {speedup_str}"
            f" | {besdq_str} | {ratio_str}{flag} |"
        )
    timing_table = "\n".join(timing_rows)

    interp_lines = []
    for t in timings:
        if "speedup_vs_row_api" in t:
            interp_lines.append(
                f"- **{t['query']}**: {t['speedup_vs_row_api']:.1f}× faster than row API "
                f"({t['row_api_median_ms']:.1f} ms → {t['median_ms']:.1f} ms)"
            )
        if "notes" in t:
            interp_lines.append(f"- {t['notes']}")
    interp_block = "\n".join(interp_lines) if interp_lines else "_No baseline comparisons available._"

    text = f"""\
---
title: "UKB chr1 Dense Store — Array Query API Benchmark"
subtitle: "chr1 × 100 UKB analyses · sparse flat array return contract"
date: today
format:
  html:
    toc: true
    toc-depth: 3
    embed-resources: true
---

## Dataset

| Field | Value |
|---|---:|
| Variants (hg38) | {dataset['n_variants']:,} |
| Analyses | {dataset['n_analyses']} |
| Store size | {storage['store_mb']:.1f} MB |
| Build time | {build_line} |
| Liftover failures | {liftover_line} |
| Store path | `{build['store_path']}` |

## Query API

All five query methods return `dict[str, np.ndarray]` with four parallel sparse
flat arrays (`variant_index`, `analysis_index`, `z`, `se`), matching the format
of the internal top-hit index. `AssociationResult` Python object materialisation
is completely eliminated.

## Timings

Each pattern ran {n_reps} repetitions after one warm-up call.
Row API column shows the previous `AssociationResult` baseline for comparison.
besdq column is the besdq zarr `zstd_bitshuffle` baseline.

| Query | Median ms | p95 ms | Result count | Row API ms | Speedup | besdq ms | vs besdq |
|---|---:|---:|---:|---:|---:|---:|---:|
{timing_table}

## Interpretation

{interp_block}

## Provenance

- JSON results: `{json_path.name}`
- Row-API baseline: `{row_baseline['build']['store_path'] if row_baseline else 'n/a'}`
- besdq baseline: `ukb-chr1_zarr_benchmark.json`
"""
    qmd_path.parent.mkdir(parents=True, exist_ok=True)
    qmd_path.write_text(text, encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--store", type=Path, default=DEFAULT_STORE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qmd", type=Path, default=DEFAULT_QMD)
    parser.add_argument("--row-baseline", type=Path, default=None,
                        help="Path to row-API baseline JSON for speedup comparison")
    parser.add_argument("--rebuild", action="store_true", default=False)
    parser.add_argument("--reps", type=int, default=10, help="Query repetitions (default 10)")
    return parser.parse_args()


if __name__ == "__main__":
    main()
