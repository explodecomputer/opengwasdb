## Parent PRD

`issues/prd-vcf-dense-build-with-liftover.md`

## What to build

A benchmark script that builds a Dense Observed-Only Store from the 100-trait chr1 UKB GWAS-VCF dataset at `/home/gh13047/repo/besdq/data/vcf-ukb/` and measures build time, storage size, and query latency across the five query patterns used by the besdq prototype benchmark.

The script records results to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.json` in the same structure used by the existing `opengwasdb_38714679_benchmark.json`. The five query patterns to benchmark are the same as `ukb-chr1_zarr_benchmark.json` in besdq:

1. **regional** — all variants in a 1 Mb window on chr1 across all traits
2. **phewas** — single variant across all traits
3. **tophits** — genome-wide-significant variants for one trait
4. **bulk** — all variants for one trait
5. **random_lookup** — 100 random variants × 10 random traits

Each query pattern runs 10 repetitions; record median and p95 latency in milliseconds. Also record: build time in seconds, store size on disk in MB, n_variants, n_analyses, and liftover failure count.

Compare results against the besdq zarr baseline in `besdq/data/ukb-chr1_zarr_benchmark.json`. The comparison is manual inspection, not an automated assertion. If cyvcf2 proves materially slower than bcftools for the build, note this explicitly in the benchmark output so the decision recorded in ADR-0019 can be revisited.

This is a benchmark script, not a test. It lives in `benchmarks/` alongside the existing benchmark scripts, not in `tests/`.

## Acceptance criteria

- [ ] The build completes successfully on the 100-trait chr1 UKB VCF dataset using `build_dense_from_vcf_manifest`.
- [ ] `validate_store` passes on the resulting store.
- [ ] All five query patterns execute without error and return non-empty results.
- [ ] Benchmark results are written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.json`.
- [ ] The JSON includes build time, store size, n_variants, n_analyses, liftover failure count, and per-query-pattern median/p95 latency.
- [ ] Results are compared against `besdq/data/ukb-chr1_zarr_benchmark.json`; any pattern where opengwasdb is >2× slower is flagged in a `notes` field in the output JSON.
- [ ] A Quarto report is written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_benchmark.qmd` following the `write_qmd()` pattern from `benchmark_38714679_dense.py`.

## Blocked by

- `issues/015-two-pass-dense-build-from-vcf-manifest.md`

## User stories addressed

- User story 7 (build time and memory comparable to besdq prototype — verified empirically)
