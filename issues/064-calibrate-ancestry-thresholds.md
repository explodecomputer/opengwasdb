## Parent PRD

`issues/prd-ancestry-assignment.md`

## What to build

Calibrate the admission thresholds (τ, δ) against the **Reported Population**, using
the full Catalogue from `issues/063-assign-ancestry-cli.md`. Fetch each Analysis's
declared `population` (OpenGWAS metadata / the `ieu-*.json` files) into the
Catalogue, cross-tabulate **Assigned vs Reported Population** across the collection,
and choose the operating point where reported-European Analyses clear the EUR gate
and reported-Mixed/other fall to Unassigned. Emit a **disagreement report** (Assigned
≠ Reported) for human review (PRD "Calibration"; ADR 0028 — Reported Population
calibrates/audits, never routes).

This is **HITL**: a human inspects the cross-tab and disagreement report and picks
τ/δ; the chosen values are recorded (e.g. as Catalogue provenance / config) so 065
builds from a calibrated Catalogue.

## Acceptance criteria

- [ ] `reported_population` is populated for each Analysis from the source metadata.
- [ ] A cross-tabulation of Assigned vs Reported Population over the collection is
      produced, plus a disagreement report listing Analyses where they conflict.
- [ ] Chosen τ/δ are recorded as provenance and re-applied to (re)label the Catalogue.
- [ ] The operating point admits (nearly) all reported-European as EUR and routes
      reported-Mixed/other to Unassigned, documented with the counts.

## Blocked by

- Blocked by `issues/063-assign-ancestry-cli.md`

## User stories addressed

- User story 4
