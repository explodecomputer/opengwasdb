#!/usr/bin/env python3
"""Compare query latency of two genome-wide ukb-b Dense stores that differ only
in zarr chunk shape: 1000x1000 vs 1000x128.

Both stores hold the same 9.85M variants x 2,514 analyses, so an identical query
selection (portable ALIDs / analysis ids / region) is run against each. The
narrower analysis chunk should speed up the per-analysis (bulk) read, which is
chunk-bound, without materially changing the per-variant (phewas) read.

Writes docs/benchmark-output/opengwasdb_chunk_comparison_benchmark.json, rendered
by the companion QMD.

Usage:
  uv run python benchmarks/benchmark_chunk_comparison.py [--reps N]
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.query import query_store

CONFIGS = [
    ("1000x1000", Path("/local-scratch/data/opengwas/opengwasdb/ukb-b.opengwasdb")),
    ("1000x128", Path("/local-scratch/data/opengwas/opengwasdb/ukb-b-c128.opengwasdb")),
]
OUTPUT = Path(
    "/home/gh13047/repo/opengwasdb/docs/benchmark-output/"
    "opengwasdb_chunk_comparison_benchmark.json"
)

# Reuse the ukb-b benchmark's anchors so the two reports line up.
EXPOSURE = "ukb-b-10912"
REGION = ("19", 44_500_000, 45_500_000)


def _median_ms(fn, reps: int) -> tuple[float, float, int]:
    fn()  # warm-up
    times = []
    count = 0
    for _ in range(reps):
        import time
        t0 = time.perf_counter()
        res = fn()
        times.append((time.perf_counter() - t0) * 1000.0)
        count = len(res["z"])
    times.sort()
    return times[len(times) // 2], times[min(len(times) - 1, int(0.95 * len(times)))], count


def _dir_bytes(path: Path) -> int:
    out = subprocess.run(["du", "-sb", str(path)], capture_output=True, text=True, check=True)
    return int(out.stdout.split()[0])


def _choose_selection() -> dict:
    """Pick one portable selection from the baseline store, applied to both."""
    q = query_store(CONFIGS[0][1])
    an = q.analyses_table()
    analyses_by_id = {v["analysis_id"]: k for k, v in an.items()}
    n_variants = int(q._root["z"].shape[0])
    n_analyses = len(an)

    th = q.top_hits(threshold=5e-8)
    m = th["analysis_index"] == analyses_by_id[EXPOSURE]
    strong_vi = int(th["variant_index"][m][np.argmax(np.abs(th["z"][m]))])
    phewas_alid = q._variant_axis.by_index(strong_vi).alid

    rng = np.random.default_rng(0)
    rand_vi = rng.choice(n_variants, size=100, replace=False)
    rand_alids = [
        r.alid for r in (q._variant_axis.by_index(int(v)) for v in rand_vi) if r is not None
    ]
    rand_a = rng.choice(n_analyses, size=10, replace=False)
    rand_analyses = [an[int(a)]["analysis_id"] for a in rand_a]
    return {
        "bulk_analysis_id": EXPOSURE,
        "phewas_alid": phewas_alid,
        "region": {"chrom": REGION[0], "start": REGION[1], "end": REGION[2]},
        "rand_alids": rand_alids,
        "rand_analyses": rand_analyses,
        "n_variants": n_variants,
        "n_analyses": n_analyses,
    }


def _bench_store(store: Path, sel: dict, reps: int) -> dict:
    q = query_store(store)
    patterns = {
        "bulk": lambda: q.analysis(sel["bulk_analysis_id"]),
        "phewas": lambda: q.phewas(sel["phewas_alid"]),
        "regional": lambda: q.range_phewas(
            sel["region"]["chrom"], sel["region"]["start"], sel["region"]["end"]
        ),
        "tophits": lambda: q.top_hits(threshold=5e-8),
        "random_lookup": lambda: q.lookup(sel["rand_alids"], sel["rand_analyses"]),
    }
    timings = []
    for name, fn in patterns.items():
        med, p95, cnt = _median_ms(fn, reps)
        timings.append({"query": name, "median_ms": round(med, 3),
                        "p95_ms": round(p95, 3), "result_count": cnt})
        print(f"    {name:15s} median={med:9.2f} ms  count={cnt:,}")
    chunks = list(zarr.open_group(str(store / "data.zarr"), mode="r")["z"].chunks)
    store_bytes = _dir_bytes(store)
    return {"store": str(store), "chunks": chunks,
            "store_gb": round(store_bytes / 1e9, 2), "timings": timings}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=7)
    args = ap.parse_args()

    sel = _choose_selection()
    print(f"selection: bulk={sel['bulk_analysis_id']} phewas={sel['phewas_alid']} "
          f"region={sel['region']['chrom']}:{sel['region']['start']}-{sel['region']['end']}")

    configs = []
    for name, store in CONFIGS:
        if not store.exists():
            raise SystemExit(f"store missing: {store}")
        print(f"benchmarking {name} ({store.name}) ...")
        entry = {"name": name, **_bench_store(store, sel, args.reps)}
        configs.append(entry)

    result = {
        "dataset": {"n_variants": sel["n_variants"], "n_analyses": sel["n_analyses"],
                    "reference_assembly": "GRCh38"},
        "reps": args.reps,
        "selection": {k: sel[k] for k in ("bulk_analysis_id", "phewas_alid", "region")},
        "configs": configs,
    }
    OUTPUT.write_text(json.dumps(result, indent=2))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
