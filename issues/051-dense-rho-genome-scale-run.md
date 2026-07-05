## Parent PRD

`issues/prd-dense-rho-matrix.md`

## What to build

Run the Rho Matrix build on the real genome-wide ukb-b Dense store and confirm it
is correct at scale and scientifically plausible. This is a **HITL** slice: the
sanity checks require human judgement about the resulting Rho distribution.

- Run `build-dense-rho` on the ukb-b store (2,514 Analyses); record wall-clock
  time, peak memory, and `n_variants_used`.
- Validate the resulting store (`opengwasdb validate`, issue 049).
- Sanity-check the Rho distribution: Analyses on overlapping samples / related
  traits show elevated Rho; unrelated, non-overlapping Analyses cluster near 0;
  the per-pair support (`n_null`) is well above `min_nulls` for most pairs.
- Produce a short report (a benchmark JSON/QMD in the existing `docs/benchmark-
  output/` style, or a note) summarising timing, support distribution, and a few
  spot-checked trait pairs.

## Acceptance criteria

- [ ] `build-dense-rho` completes on the ukb-b store within a reasonable
      wall-clock (minutes, not hours) and bounded memory, both recorded.
- [ ] The completed store validates.
- [ ] The Rho distribution passes human review: known overlapping/related trait
      pairs are elevated, unrelated pairs are near 0, and support counts are
      sensible; findings are written up.
- [ ] Any surprises (unexpectedly high Rho for unrelated pairs, many
      low-support/NaN pairs) are surfaced with a hypothesis (e.g. `window_bp` too
      small, residual LD) rather than silently accepted.

## Blocked by

- Blocked by `issues/048-dense-rho-wide-query-api.md`
- Blocked by `issues/049-dense-rho-validation.md`
- Blocked by `issues/050-dense-rho-parallel-and-cli.md`

## User stories addressed

- User story 1 (validated at scale)
- User story 9 (validated at scale)
