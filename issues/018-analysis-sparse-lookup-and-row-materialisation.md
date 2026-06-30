## Parent PRD

`issues/prd.md`

## What to build

Use the row-offset sidecar to materialise variant metadata by dense row index and row-index batches. This should preserve full-analysis extraction and targeted variant-by-analysis lookup without reintroducing a SQLite variants table.

## Acceptance criteria

- [ ] A single variant row can be materialised by dense row index using the offset sidecar.
- [ ] A batch of variant rows can be materialised by dense row index using the offset sidecar.
- [ ] `analysis_arrays()` continues to return dense arrays without materialising every variant row.
- [ ] `analysis()` returns finite materialised rows for one analysis using the new variant axis.
- [ ] `lookup()` resolves requested variants and analyses, reads the dense block, and returns finite materialised rows.
- [ ] Missing dense cells remain excluded from `analysis()` and `lookup()` results.
- [ ] Query result fields remain stable for downstream callers.

## Blocked by

- Blocked by `issues/016-canonical-alid-and-range-queries.md`

## User stories addressed

- User story 7
- User story 15
- User story 19
- User story 26
- User story 27
- User story 29
- User story 30
