## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Rebuild the chr1 UKB store (to pick up the new sidecar files) and re-run the benchmark.
Update `opengwasdb_vs_besdq_comparison.qmd` with the post-fix timings. Record result
counts alongside timings so regional comparisons note the variant count difference
between benchmark windows.

Expected post-fix results:
- regional: ~30 ms (down from 78 ms) — driven by `range_indices()` removing object allocation
- random_lookup: ~70 ms (down from 114 ms) — driven by mmap index removing sequential fetches
- All other patterns: unchanged or marginally faster

## Acceptance criteria

- [ ] Store rebuilt with sidecar files present
- [ ] Benchmark runs without error; results written to JSON and QMD
- [ ] `opengwasdb_vs_besdq_comparison.qmd` updated with post-fix timings and a
      before/after table
- [ ] regional result count and variant count noted in benchmark output

## Blocked by

- `issues/028-range-indices-variant-axis.md`
- `issues/030-mmap-alid-sidecar-read.md`
- `issues/031-validate-alid-sidecar.md`

## User stories addressed

- User story 4 (benchmark reflects accurate per-pattern timings)
