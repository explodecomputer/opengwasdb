## Parent PRD

`issues/prd-ancestry-assignment.md`

## What to build

Make reference completion **ancestry-matched per Analysis** (ADR 0028,
"Ancestry-Matched Completion"): given an LD panel of a particular ancestry,
completion imputes only the Analyses whose **Assigned Ancestry** matches that panel
and leaves the others observed-only, recording `completed_against` (the applied
panel ancestry, or null) per Analysis. Store-level Completion State stays a coarse
release flag; per-cell Association Status remains the ground truth for what was
imputed.

For the homogeneous EUR store (065) this is a no-op (every Analysis matches), so the
slice is demonstrated on a **synthetic mixed store** where only the EUR-assigned
Analyses are imputed and the rest stay observed-only.

## Acceptance criteria

- [ ] Completion reads per-Analysis Assigned Ancestry and imputes only Analyses whose
      ancestry matches the applied panel; non-matching Analyses are untouched.
- [ ] Each Analysis records `completed_against` (panel ancestry or null); store-level
      Completion State is unchanged in shape.
- [ ] On a synthetic mixed store, only the matching-ancestry Analyses gain imputed
      cells (Association Status `imputed`); the rest report only `observed`.
- [ ] The EUR homogeneous store (065) completes identically to today (all Analyses
      match), i.e. no regression.

## Blocked by

- Blocked by `issues/065-catalogue-eur-hybrid-build.md`

## User stories addressed

- User story 8
