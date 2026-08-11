#!/usr/bin/env python3
"""Build and benchmark a Dense Observed-Only Store from UKB chr1 VCFs.

Uses the chr1 GWAS-VCF dataset at
  /home/gh13047/repo/besdq/data/vcf-ukb-1000/
which contains bgzipped, GRCh37/hg19 GWAS-VCFs in GWAS-SSF format.

The store is built via the two-pass ``build_dense_from_vcf_manifest`` pipeline
with inline hg19 → hg38 liftover.  Five query patterns are benchmarked (10
repetitions each).

Outputs:
  docs/benchmark-output/opengwasdb_vcf_ukb_chr1_1000_benchmark.json
  docs/benchmark-output/opengwasdb_vcf_ukb_chr1_128_benchmark.json

Usage:
  python benchmarks/benchmark_vcf_ukb_chr1_1000_dense.py [--analysis-count N] [--rebuild] [--reps N]
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

DEFAULT_MANIFEST = Path("/home/gh13047/repo/besdq/data/vcf-ukb-1000/manifest.tsv")
DEFAULT_STORE_DIR = Path("/local-scratch/data/opengwas/opengwasdb")
DEFAULT_OUTPUT_DIR = Path("/home/gh13047/repo/opengwasdb/docs/benchmark-output")
MANIFEST_SCRATCH_DIR = Path("/tmp/opengwasdb-benchmark-manifests")
DEFAULT_ANALYSIS_COUNT = 1000
SLOWDOWN_FLAG_THRESHOLD = 2.0

# Match besdq's benchmark region (dense_05_query_benchmark.py REGION constant)
BESDQ_REGION = ("1", 100_000_000, 101_000_000)

RNG = np.random.default_rng(42)


def main() -> None:
    args = _parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build_manifest = (
        _manifest_for_analysis_count(args.manifest, args.analysis_count, MANIFEST_SCRATCH_DIR)
        if args.manifest.exists()
        else args.manifest
    )

    build_seconds = None

    if args.rebuild and args.store.exists():
        shutil.rmtree(args.store)

    if args.rebuild or not args.store.exists():
        print(f"Building store from {build_manifest} ...")
        t0 = time.perf_counter()
        result = build_dense_from_vcf_manifest(
            build_manifest,
            args.store,
            store_id=f"ukb-chr1-vcf-{args.analysis_count}",
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
    results = _run_benchmark(args.store, args.reps, build_seconds, build_manifest)

    args.output.write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Results written to {args.output}")
    _print_summary(results)


def _run_benchmark(
    store_path: Path,
    n_reps: int,
    build_seconds: float | None,
    manifest_path: Path,
) -> dict:
    query = query_store(store_path)
    selection = _choose_queries(store_path)

    query_specs = {
        "regional": lambda: query.range_phewas(
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
        n_variants = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        n_analyses = conn.execute("SELECT COUNT(*) FROM analyses").fetchone()[0]

    store_bytes = sum(f.stat().st_size for f in store_path.rglob("*") if f.is_file())
    raw_vcf_bytes, n_source_files = _raw_vcf_bytes(manifest_path)

    return {
        "dataset": {
            "n_variants": n_variants,
            "n_analyses": n_analyses,
        },
        "build": {
            "store_path": str(store_path),
            "build_seconds": round(build_seconds, 2) if build_seconds is not None else None,
        },
        "storage": {
            "store_bytes": store_bytes,
            "store_mb": round(store_bytes / 1_000_000, 2),
            "raw_vcf_bytes": raw_vcf_bytes,
            "raw_vcf_mb": round(raw_vcf_bytes / 1_000_000, 2),
            "n_source_files": n_source_files,
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

    region_chrom, region_start, region_end = BESDQ_REGION

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


def _manifest_for_analysis_count(manifest_path: Path, analysis_count: int, output_dir: Path) -> Path:
    """Return a manifest limited to the first ``analysis_count`` rows."""
    if analysis_count <= 0:
        raise ValueError("--analysis-count must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(manifest_path, encoding="utf-8") as src:
        header = src.readline()
        rows = [line for _, line in zip(range(analysis_count), src, strict=False)]

    if len(rows) < analysis_count:
        raise ValueError(
            f"{manifest_path} contains only {len(rows)} analyses, "
            f"cannot select {analysis_count}"
        )

    limited = output_dir / f"ukb_chr1_{analysis_count}_manifest.tsv"
    limited.write_text(header + "".join(rows), encoding="utf-8")
    return limited


def _raw_vcf_bytes(manifest_path: Path) -> tuple[int, int]:
    if not manifest_path.exists():
        return 0, 0
    paths = []
    with open(manifest_path, encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                paths.append(Path(parts[1]))
    total = 0
    for path in paths:
        try:
            total += path.stat().st_size
        except OSError:
            pass
    return total, len(paths)


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


def _print_summary(results: dict) -> None:
    print("\n--- Summary ---")
    print(f"  Variants: {results['dataset']['n_variants']:,}")
    print(f"  Analyses: {results['dataset']['n_analyses']}")
    bs = results["build"]["build_seconds"]
    if bs is not None:
        print(f"  Build time: {bs:.1f}s")
    print(f"  Store size: {results['storage']['store_mb']:.1f} MB")
    print("")
    print(f"  {'Query':<15} {'Median ms':>10} {'p95 ms':>10} {'Result count':>14}")
    for row in results["timings"]:
        print(f"  {row['query']:<15} {row['median_ms']:>10.2f} {row['p95_ms']:>10.2f} {row['result_count']:>14,}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--analysis-count", type=int, default=DEFAULT_ANALYSIS_COUNT)
    parser.add_argument("--store", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--rebuild", action="store_true", default=False)
    parser.add_argument("--reps", type=int, default=10)
    args = parser.parse_args()
    if args.store is None:
        args.store = DEFAULT_STORE_DIR / f"vcf_ukb_chr1_{args.analysis_count}_dense.opengwasdb"
    if args.output is None:
        args.output = (
            DEFAULT_OUTPUT_DIR
            / f"opengwasdb_vcf_ukb_chr1_{args.analysis_count}_benchmark.json"
        )
    return args


if __name__ == "__main__":
    main()
