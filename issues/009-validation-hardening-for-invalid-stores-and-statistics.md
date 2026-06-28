## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Harden validation beyond the happy path. The validator should reject invalid Dense Observed-Only stores and invalid statistic arrays with actionable errors.

## Acceptance criteria

- [ ] Validator rejects duplicate canonical variants.
- [ ] Validator rejects missing required arrays or dimension mismatches.
- [ ] Validator rejects negative finite SE values.
- [ ] Validator rejects inconsistent Z/SE missingness.
- [ ] Validator rejects invalid Stored Effect Scale values.
- [ ] Validator rejects top-hit indexes inconsistent with stored Z values when top-hit indexes are present.
- [ ] Tests cover each invalid-store case with clear expected error messages.

## Blocked by

- Blocked by `issues/003-build-tiny-dense-observed-only-store.md`
- Blocked by `issues/006-derived-beta-and-pvalue-result-contract.md`
- Blocked by `issues/007-dense-top-hit-index-and-query.md`

## User stories addressed

- User story 7
- User story 12
- User story 15
- User story 17
- User story 18
- User story 39

