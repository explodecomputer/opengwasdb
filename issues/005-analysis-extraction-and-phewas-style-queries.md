## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Add analysis-centric query paths: extract all available associations for one Analysis and extract one variant across all Analyses. These should use the same layout-independent query facade and return the same result contract as exact/range queries.

## Acceptance criteria

- [ ] A full-analysis extraction returns all finite associations for one Analysis.
- [ ] A variant-across-analyses query returns all finite associations for one variant across Analyses.
- [ ] Query results use the same field names and derived-statistic behaviour as exact/range queries.
- [ ] Missing Dense cells do not corrupt output ordering or result counts.
- [ ] Tests cover one full-analysis extraction and one PheWAS-style variant extraction from a tiny store.

## Blocked by

- Blocked by `issues/004-layout-independent-exact-variant-and-range-queries.md`

## User stories addressed

- User story 23
- User story 24
- User story 25
- User story 30
- User story 31
- User story 34
- User story 37

