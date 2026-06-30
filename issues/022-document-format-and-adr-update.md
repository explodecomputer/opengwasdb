## Parent PRD

`issues/prd.md`

## What to build

Update the store-format specification, ADRs, and benchmark Quarto documentation after the benchmark confirms the tabix-backed variant axis is the preferred dense variant-axis backend. The active dense vertical-slice documentation should describe the EBI full GWAS-SSF benchmark only.

## Acceptance criteria

- [ ] The store-format specification describes the dense variant-axis files and their relationship to dense Zarr arrays.
- [ ] ADRs record why high-cardinality dense variant metadata moved out of SQLite.
- [ ] Quarto benchmark documentation names `besdq/data/ebi_input` as the dense vertical-slice benchmark dataset.
- [ ] Quarto benchmark documentation does not compare dense vertical-slice results against the sparse `38714679` dataset.
- [ ] The docs index links only the active dense vertical-slice benchmark report.
- [ ] Archived sparse smoke-test artifacts, if retained, are clearly labelled as not part of dense vertical-slice evaluation.

## Blocked by

- Blocked by `issues/021-benchmark-storage-and-build-phases.md`

## User stories addressed

- User story 8
- User story 14
- User story 16
- User story 37
