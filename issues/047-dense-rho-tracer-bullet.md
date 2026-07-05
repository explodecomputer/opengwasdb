## Parent PRD

`issues/prd-dense-rho-matrix.md`

## What to build

The thin end-to-end path for the Dense Rho Matrix, serial, verified on a small
store. Compute Rho for every Analysis pair and store + query it:

- **Variant selection**: distance-thin the store's own axis — fixed-bp windows
  (`window_bp`, default sized to ≈200k genome-wide), keep the first variant per
  window (PRD "Independent variant selection").
- **Inputs**: observed Z only (imputed treated as missing); per pair, the
  both-null subset `isfinite(z_j) & isfinite(z_k) & |z_j| < z_thresh &
  |z_k| < z_thresh`.
- **Estimator**: port `estimate_rho_cml` (sufficient stats `A,B,C,n`; truncated
  bivariate-normal NLL minimised over ρ with the `_GRID_N=2000` normalising grid;
  NaN when `n < min_nulls`) from pleiodb `rho.py`, reused verbatim (PRD
  "Estimator"). Defaults `z_thresh=1.0`, `min_nulls=500`.
- **Storage**: write a `data.zarr/rho` group — packed strict lower triangle
  `rho` (float16) + `n_null` (int32), group provenance attrs (`z_thresh`,
  `min_nulls`, `grid_n`, `window_bp`, `n_variants_used`, `observed_only`,
  `method`, `n_analyses`), and the thinned `variant_index` array (PRD "Storage").
- **Query (minimal)**: `rho(*ids)` returning long format
  (`analysis_id_a`, `analysis_id_b`, `rho`, `n_null`), self-pairs excluded,
  diagonal implied = 1 (PRD "Query API").

Serial compute is fine here; parallelism and the CLI come in `issues/050-*`.

## Acceptance criteria

- [ ] A `build_dense_rho(store_path, *, window_bp, z_thresh, min_nulls)` Python
      entry point produces a `data.zarr/rho` group on a small synthetic Dense store.
- [ ] The ported `estimate_rho_cml` reproduces pleiodb's output on shared
      fixtures (same `z_j`, `z_k`, `z_thresh`, `min_nulls`) to float tolerance,
      including the `NaN` return when support `< min_nulls`.
- [ ] Rho is symmetric; a pair's stored value equals a direct CML computation on
      that pair's observed both-null Z-scores at the thinned variants.
- [ ] `n_null` equals the both-null shared-variant count; `rho` is NaN exactly
      when `n_null < min_nulls`.
- [ ] Imputed cells (if any) are excluded from the estimate (observed-only).
- [ ] `rho(*ids)` returns the correct long-format pairwise result for 2 and for
      N ids, excluding self-pairs.
- [ ] Provenance attrs and the thinned `variant_index` array are written and
      round-trip.

## Blocked by

None - can start immediately.

## User stories addressed

- User story 1
- User story 2
- User story 5
- User story 6
- User story 7
- User story 12
- User story 13
- User story 14
