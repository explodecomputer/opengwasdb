## Parent PRD

`issues/prd-ancestry-assignment.md`

## What to build

Turn the calibrated Catalogue into a EUR **Hybrid** store. Filter the Catalogue to
`Assigned Ancestry == EUR` to produce a build manifest (a row-subset of the same
TSV — no format translation), run the existing `build-hybrid` + `complete-hybrid`
(EUR LD panel) **unchanged**, and record provenance linking the store back to the
Catalogue (PRD user stories 6, 7, 9; ADR 0027).

- The store's analyses table records **per-Analysis Assigned Ancestry** (inherited
  from the subset), plus `catalogue_version` and the subset filter as provenance.
- Non-EUR/Unassigned Analyses are simply absent from the store — they remain parked
  in the Catalogue.

## Acceptance criteria

- [ ] A EUR build manifest is produced purely by row-filtering the Catalogue; the
      unchanged manifest reader consumes it.
- [ ] `build-hybrid` + `complete-hybrid` produce a EUR hybrid store that validates.
- [ ] The store records `catalogue_version`, the subset filter, and each Analysis's
      Assigned Ancestry.
- [ ] Only EUR-assigned Analyses appear in the store; non-EUR/Unassigned are absent
      (still present in the Catalogue).

## Blocked by

- Blocked by `issues/063-assign-ancestry-cli.md`
- Blocked by `issues/064-calibrate-ancestry-thresholds.md`

## User stories addressed

- User story 6
- User story 7
- User story 9
