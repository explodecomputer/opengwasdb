## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Implement a Dense top-hit index and layout-independent top-hit query using the accepted thresholds. Reference the `besdq` dense top-hit index prototype, but adapt it to the OpenGWASDB `z` and `se` contract and future-compatible result stream.

## Acceptance criteria

- [ ] A built Dense store can create top-hit indexes for the accepted thresholds.
- [ ] The index excludes NaN Z values.
- [ ] Top-hit queries return associations ranked by significance.
- [ ] Top-hit results use the same result contract as other query paths, including derived beta and p-value when requested.
- [ ] Tests verify thresholds, offsets/index contents, result ordering, and NaN exclusion on a tiny store.

## Blocked by

- Blocked by `issues/003-build-tiny-dense-observed-only-store.md`
- Blocked by `issues/006-derived-beta-and-pvalue-result-contract.md`

## User stories addressed

- User story 26
- User story 27
- User story 28
- User story 30
- User story 31
- User story 34
- User story 37

