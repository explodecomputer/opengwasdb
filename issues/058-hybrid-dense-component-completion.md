## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

Reference-complete a Hybrid store by imputing **only the Dense Component**,
leaving the Ragged Overflow Component untouched. The Dense Component's axis is the
LD reference panel, so completion is exactly the existing dense completion problem
— reuse the dense completion pipeline (ADR 0022/0023) unchanged on that component.

The overflow is off-panel (no LD structure) and stays observed-only. A completed
Hybrid store's queries then return imputed + observed cells from the dense
component and observed cells from the overflow (PRD "Completion").

## Acceptance criteria

- [ ] `complete` on a Hybrid store produces a new Reference-Completed Hybrid
      release whose Dense Component has imputed cells (and `imputed`/`on_panel`
      arrays) and whose Ragged Overflow Component is byte-identical to the source
      (observed-only, untouched).
- [ ] Completion runs the existing dense pipeline on the dense component — no new
      imputation logic and no overflow LD handling.
- [ ] Queries on the completed store report `association_status = imputed` only
      for dense-component cells; every overflow association stays `observed`.
- [ ] The completed store declares `completion_state = reference_completed` and
      validates (see `issues/060-...` once available).

## Blocked by

- Blocked by `issues/055-hybrid-tracer-bullet.md`

## User stories addressed

- User story 8
- User story 14
