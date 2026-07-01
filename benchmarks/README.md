# Benchmarks

Each script builds (or reuses) a store, runs timed queries, and writes its results
to `docs/benchmark-output/` as a JSON file.  The comparison Quarto document
(`docs/benchmark-output/opengwasdb_vs_besdq_comparison.qmd`) reads those JSON
files at render time — re-running a script and re-rendering the QMD is all that
is needed to reproduce or update the report.

---

## Scripts

### `benchmark_vcf_ukb_chr1_dense.py`

Builds and benchmarks a dense observed-only store from the 100 UKB chr1 GWAS-VCF
dataset.  Requires the VCFs at `/home/gh13047/repo/besdq/data/vcf-ukb/` and the
besdq repo at `/home/gh13047/repo/besdq/`.

**Output files written to `docs/benchmark-output/`:**

| File | Description |
|---|---|
| `opengwasdb_vcf_ukb_chr1_benchmark.json` | Post-optimisation query timings |
| `besdq_ukb_chr1_benchmark.json` | besdq baseline (copied from besdq repo) |
| `opengwasdb_vcf_ukb_chr1_benchmark.qmd` | Per-run standalone QMD |

**Usage:**

```bash
# First run — build the store and benchmark (takes ~10 min)
conda run -n snakemake python benchmarks/benchmark_vcf_ukb_chr1_dense.py --rebuild --reps 10

# Subsequent runs — reuse existing store, re-benchmark only
conda run -n snakemake python benchmarks/benchmark_vcf_ukb_chr1_dense.py --reps 10
```

The `--row-baseline` flag accepts a path to an earlier JSON to show speedup ratios:

```bash
python benchmarks/benchmark_vcf_ukb_chr1_dense.py \
    --row-baseline docs/benchmark-output/opengwasdb_vcf_ukb_chr1_array_benchmark.json \
    --reps 10
```

---

### `benchmark_ragged_besd.py`

Benchmarks a ragged observed-only store built from BESD files.  Six query
patterns are timed and a storage comparison against the source BESD is included.
Defaults to the pre-built eqtlgen-cis store.

**Output files written to `docs/benchmark-output/`:**

| File | Description |
|---|---|
| `opengwasdb_eqtlgen_ragged_benchmark.json` | Query timings + storage comparison |
| `opengwasdb_eqtlgen_ragged_benchmark.qmd` | Self-contained Quarto report (rendered to HTML) |

**Usage:**

```bash
# Benchmark existing store (no rebuild)
conda run -n snakemake python benchmarks/benchmark_ragged_besd.py --reps 5

# Force a full rebuild then benchmark
conda run -n snakemake python benchmarks/benchmark_ragged_besd.py --rebuild --reps 5

# Use a different BESD source (e.g. hg19 with liftover)
conda run -n snakemake python benchmarks/benchmark_ragged_besd.py \
    --besd /path/to/prefix \
    --store /path/to/out.opengwasdb \
    --source-build hg19 \
    --tissue Whole_Blood
```

**Render the QMD:**

```bash
cd docs/benchmark-output
QUARTO_PYTHON=/home/gh13047/miniforge3/envs/snakemake/bin/python \
  /home/gh13047/miniforge3/bin/quarto render opengwasdb_eqtlgen_ragged_benchmark.qmd
```

**Query patterns:**

| Pattern | Description |
|---|---|
| `analysis` | All cis associations for one probe (O(1) CSR slice) |
| `range_by_probe` | All analyses whose TSS falls in a 2 Mb window |
| `range` | All associations where the variant falls in a 2 Mb window (O(n) scan) |
| `phewas` | One variant across all analyses (O(n) scan) |
| `tophits` | Top-10 hits by \|z\| from precomputed index |
| `random_lookup` | 100 random variants × 10 random analyses |

---

## Comparison document

After all JSONs are present in `docs/benchmark-output/`, render the comparison report:

```bash
QUARTO_PYTHON=/home/gh13047/miniforge3/envs/snakemake/bin/python \
  /home/gh13047/miniforge3/bin/quarto render docs/benchmark-output/opengwasdb_vs_besdq_comparison.qmd
```

The rendered HTML is written to the same directory.

---

## JSON schema

All opengwasdb benchmark JSONs share a common top-level structure:

```json
{
  "dataset":  { "n_variants": int, "n_analyses": int },
  "build":    { "store_path": str, "build_seconds": float|null, "liftover_failure_count": int },
  "storage":  { "store_bytes": int, "store_mb": float },
  "selection": { ... query parameters used ... },
  "timings": [
    {
      "query": str,
      "median_ms": float,
      "p95_ms": float,
      "result_count": int,
      "besdq_zstd_median_ms": float,   // present if besdq baseline loaded
      "ratio_vs_besdq": float,         // present if besdq baseline loaded
      "row_api_median_ms": float,      // present if --row-baseline supplied
      "speedup_vs_row_api": float,     // present if --row-baseline supplied
      "notes": str                     // present if ratio > 2×
    }
  ]
}
```

The besdq baseline JSON (`besdq_ukb_chr1_benchmark.json`) uses a different
structure produced by `besdq/scripts/dense_05_query_benchmark.py`:

```json
{
  "zstd_bitshuffle": {
    "regional":      { "median_ms": float, ... },
    "phewas":        { "median_ms": float, ... },
    ...
  },
  "raw_float16": { ... }
}
```
