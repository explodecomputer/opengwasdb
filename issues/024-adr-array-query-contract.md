## Parent PRD

`issues/prd-query-array-api.md`

## What to build

Write a new ADR documenting the removal of Python-object row materialisation from the query engine and the adoption of sparse flat numpy arrays as the public query contract.

The ADR should record:
- What was removed: `AssociationResult` frozen dataclass, per-row `beta` and `p_value` computation, `analysis_arrays()` companion method
- Why: chr1 × 100-analysis benchmark showed 28–562× slowdown vs besdq prototype; root cause was constructing O(result_count) Python objects with field access and arithmetic per row
- What replaced it: sparse flat arrays `{variant_index, analysis_index, z, se}` matching the format the top-hit index already uses
- What was explicitly considered and rejected: keeping row methods alongside array methods (two APIs for the same thing, no clear ownership boundary); deprecation wrappers (no existing external callers at this stage)
- Relationship to ADR-0006 (layout-independent query engine): the facade contract is amended; the layout-independence guarantee is preserved

## Acceptance criteria

- [ ] ADR file written to `docs/adr/` with the next available number.
- [ ] ADR records the performance evidence (benchmark ratios) that motivated the change.
- [ ] ADR notes the top-hit index format as the prior art for the sparse array shape.
- [ ] ADR references ADR-0006 as the contract it amends.

## Blocked by

- `issues/023-sparse-array-query-contract.md`

## User stories addressed

- Implicit in all implementation decisions in `issues/prd-query-array-api.md`
