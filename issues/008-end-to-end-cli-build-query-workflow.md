## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Expose a demoable CLI workflow for the first vertical slice: build a tiny Dense Observed-Only store, validate/show info for it, and run basic queries against an explicit store path.

## Acceptance criteria

- [ ] CLI can build a Dense Observed-Only store from fixture-style input into an explicit output path.
- [ ] CLI can validate and print info for the resulting store.
- [ ] CLI can run exact variant, range, full-analysis, variant-across-analyses, and top-hit queries.
- [ ] CLI outputs are deterministic enough for tests.
- [ ] Tests exercise the CLI workflow end to end using a temporary directory.

## Blocked by

- Blocked by `issues/003-build-tiny-dense-observed-only-store.md`
- Blocked by `issues/004-layout-independent-exact-variant-and-range-queries.md`
- Blocked by `issues/005-analysis-extraction-and-phewas-style-queries.md`
- Blocked by `issues/007-dense-top-hit-index-and-query.md`

## User stories addressed

- User story 1
- User story 3
- User story 21
- User story 22
- User story 23
- User story 24
- User story 25
- User story 26
- User story 29
- User story 33

