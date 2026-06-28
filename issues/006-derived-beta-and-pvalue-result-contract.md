## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Implement and test the derived statistic contract for query results. Beta is derived as `z * se`; p-value is derived from Z. The behaviour must be shared by all v0.1 query paths and compatible with future layouts.

## Acceptance criteria

- [ ] Query results can include beta derived from `z * se`.
- [ ] Query results can include p-value derived from Z.
- [ ] SE is validated as non-negative for all finite stored values.
- [ ] Derived fields behave consistently across exact variant, range, analysis, and variant-across-analyses queries.
- [ ] Tests cover positive, negative, zero, and missing Z/SE cases.

## Blocked by

- Blocked by `issues/004-layout-independent-exact-variant-and-range-queries.md`

## User stories addressed

- User story 11
- User story 12
- User story 13
- User story 14
- User story 30
- User story 31

