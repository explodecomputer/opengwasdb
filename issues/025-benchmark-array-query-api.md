## Parent PRD

`issues/prd-query-array-api.md`

## What to build

Re-run `benchmarks/benchmark_vcf_ukb_chr1_dense.py` against the existing chr1 × 100-analysis store (no rebuild needed — the store is already on disk) using the new array query methods. Record results to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_array_benchmark.json`.

The benchmark script needs updating: the five query lambdas currently call the row-returning methods which no longer exist. Replace each lambda with the equivalent sparse-array call. Result count is `len(result["z"])` rather than `len(result)`.

Compare the new timings against both:
- The previous row-materialisation baseline (`opengwasdb_vcf_ukb_chr1_benchmark.json`)
- The besdq zstd_bitshuffle baseline (`besdq/data/ukb-chr1_zarr_benchmark.json`)

Flag any pattern still >2× slower than besdq. If `top_hits` or `phewas` are still >2× slower the gap is likely zarr read overhead or SQLite lookup overhead rather than object materialisation — note this explicitly.

A Quarto report is written alongside the JSON to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_array_benchmark.qmd`, following the same `write_qmd()` pattern used in `benchmark_38714679_dense.py`. The report includes dataset summary, build provenance, storage footprint, query selection, timing table, and interpretation notes — including the comparison against the row-materialisation baseline and besdq.

## Acceptance criteria

- [ ] Benchmark script runs without error using the array query API.
- [ ] Results written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_array_benchmark.json`.
- [ ] Quarto report written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_array_benchmark.qmd`.
- [ ] JSON includes timings for all five patterns plus a comparison against both baselines.
- [ ] QMD renders without error (`quarto render <path>`).
- [ ] Any pattern still >2× slower than besdq has a `notes` field identifying the likely bottleneck.
- [ ] `bulk` and `range` show meaningful improvement over the row-materialisation baseline (expected large gains from removing Python object construction).

## Blocked by

- `issues/023-sparse-array-query-contract.md`

## User stories addressed

- Validates the performance target stated in `issues/prd-query-array-api.md`
