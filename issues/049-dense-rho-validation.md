## Parent PRD

`issues/prd-dense-rho-matrix.md`

## What to build

Extend `validate_store` to check the `data.zarr/rho` group when it is present (a
store without a Rho Matrix stays valid). Validation is structural/consistency
only — it does not re-run the estimator. See PRD "Validation".

Checks:
- Packed `rho` and `n_null` lengths both equal `n_analyses (n_analyses − 1) / 2`,
  and `n_analyses` matches the store's analyses count.
- `rho` is finite exactly where `n_null >= min_nulls`, and NaN where
  `n_null < min_nulls` (drawn from the stored `min_nulls` attr).
- Every finite `rho` lies within `[-1, 1]`.
- Required provenance attrs are present (`z_thresh`, `min_nulls`, `grid_n`,
  `window_bp`, `n_variants_used`, `observed_only`, `method`, `n_analyses`), and
  the `variant_index` provenance array length equals `n_variants_used`.

## Acceptance criteria

- [ ] A store with no `rho` group validates unchanged (feature is optional).
- [ ] A well-formed Rho Matrix (from issue 047) passes validation.
- [ ] Corruption is caught: wrong packed length, a finite `rho` with
      `n_null < min_nulls` (or NaN with sufficient support), an out-of-range
      `rho`, a missing provenance attr, and a `variant_index` length mismatch each
      fail validation with an actionable message.
- [ ] The check is added to the dense-validation path and covered by tests in
      `tests/test_validation.py` (or the dense validation test module).

## Blocked by

- Blocked by `issues/047-dense-rho-tracer-bullet.md`

## User stories addressed

- User story 11
