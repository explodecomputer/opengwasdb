## Parent PRD

`issues/prd.md`

## What to build

Update the dense benchmark to measure the new variant-axis design against the current dense vertical-slice workload. The benchmark should use the three full GWAS-SSF files from `besdq/data/ebi_input` and should not use or compare against the sparse `38714679` dataset.

## Acceptance criteria

- [ ] The benchmark reports source input size, Zarr size, SQLite size, variant TSV size, tabix index size, row-offset sidecar size, top-hit size, and total store size.
- [ ] The benchmark reports build timing split by normalisation, variant-axis writing, SQLite metadata writing, dense array writing, and top-hit indexing where practical.
- [ ] The benchmark reports the existing query families: regional, PheWAS, top hits, bulk-analysis arrays, and random variant-by-analysis lookup.
- [ ] The benchmark uses `besdq/data/ebi_input` as the dense vertical-slice dataset.
- [ ] The benchmark does not run, compare, or document the `38714679` sparse dataset as part of dense vertical-slice evaluation.
- [ ] Benchmark output makes storage regressions attributable to a specific store component.

## Blocked by

- Blocked by `issues/020-end-to-end-cli-store-build.md`

## User stories addressed

- User story 12
- User story 13
- User story 37
