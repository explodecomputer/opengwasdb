## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

Extend `validate_store` to validate a Hybrid store: run the existing dense and
ragged validators on the two components, and add the hybrid-specific invariants
(PRD "Validation"). Reuse the component validators; add only the cross-component
checks.

Hybrid invariants:
- The shared Store Variant Table covers every `variant_index` referenced by both
  components.
- **Disjoint partition**: on-panel variants (Dense Component rows) never appear in
  the overflow; every overflow variant is off-panel (not a dense row).
- Every Ragged Overflow association has `association_status = observed` (the
  overflow is never imputed), including in a Reference-Completed Hybrid store.
- `primary_layout = hybrid` and `association_coverage = full` in the manifest.

## Acceptance criteria

- [ ] A well-formed Hybrid store (observed-only and reference-completed) validates.
- [ ] The dense and ragged component validators are reused, not reimplemented.
- [ ] Corruption is caught: an overflow entry for an on-panel variant; a dense row
      that is also in the overflow; an imputed overflow association; a
      `variant_index` outside the shared table — each fails with an actionable
      message.
- [ ] A store without the hybrid layout is unaffected.

## Blocked by

- Blocked by `issues/055-hybrid-tracer-bullet.md`

## User stories addressed

- User story 9
