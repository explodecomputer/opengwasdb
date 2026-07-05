# PRD: Dense Store Pairwise Rho Matrix

## Problem Statement

Analyses in a Dense Store are frequently correlated for reasons unrelated to
shared biology: they are computed on **overlapping samples** (the same UK
Biobank participants appear in thousands of GWAS) and their traits are
**phenotypically correlated**. Downstream cross-trait methods — multi-trait
conditional analysis, pleiotropy testing, cross-trait LD-score-style
corrections, and any method that combines Z-scores across Analyses — need to
know this nuisance correlation to avoid false positives. Today OpenGWASDB stores
no such quantity, so every downstream user must recompute it from raw
summary statistics, which is slow, non-reproducible, and inconsistent across
callers.

The needed quantity is **Rho** (see CONTEXT.md): the correlation between two
Analyses' statistics under the null, equal to phenotypic correlation × the
proportion of overlapping samples. It is recoverable from pairs of
**non-significant** Z-scores at approximately independent variants, and the
established estimator is the conditional maximum-likelihood (CML) method in
[`pleiodb/src/pleiodb/rho.py`](https://github.com/explodecomputer/pleiodb/blob/main/src/pleiodb/rho.py).

## Solution

Add an optional, opt-in, post-build step for **Dense stores** that estimates
Rho between every pair of Analyses and stores the result as a compact,
queryable **Rho Matrix** inside the release. A new `build-dense-rho` command
(and Python entry point) selects ~200k approximately-independent variants by
distance-thinning the store's own variant axis, estimates Rho for every
Analysis pair using the pleiodb CML estimator on their shared **observed**
non-significant Z-scores, and writes a packed lower-triangle float16 matrix plus
per-pair support counts and provenance into `data.zarr/rho`. The query facade
exposes long-format, one-vs-all, and wide-matrix accessors. The source store is
mutated in place by adding the `rho` group (as the top-hit index is); no new
release is minted.

Rho is a nuisance/summary quantity, not association data: it does not change any
existing query result and is absent unless the step has been run.

## User Stories

1. As a statistical geneticist, I want the pairwise Rho between two Analyses so
   that I can correct a cross-trait test for sample overlap and phenotypic
   correlation without recomputing it from raw data.
2. As a statistical geneticist, I want to pass a vector of Analysis IDs and get
   the pairwise Rho for all of them in **long format** (id_a, id_b, rho,
   support) so that I can feed it straight into a dataframe.
3. As a statistical geneticist, I want one Analysis's Rho against all others so
   that I can screen a focal trait for overlap-driven correlation.
4. As a statistical geneticist, I want the full (or a submatrix) Rho matrix in
   **wide format** so that I can do matrix algebra (e.g. whitening) or draw a
   heatmap.
5. As a statistical geneticist, I want the **per-pair support count** (how many
   shared null variants backed each estimate) so that I can down-weight or
   discard marginally-supported Rho values.
6. As a statistical geneticist, I want pairs with fewer than `min_nulls` shared
   null variants to be **NaN** rather than a noisy number so that I never
   silently use an unreliable estimate.
7. As a statistical geneticist, I want Rho estimated from **observed** Z-scores
   only (imputed cells excluded) so that the estimate reflects true phenotypic
   overlap, not the imputation/LD model.
8. As a data engineer, I want a `build-dense-rho <store>` CLI command with
   tunable `--window-bp`, `--z-thresh`, `--min-nulls`, and `--n-workers` so that
   I can produce the Rho Matrix without writing Python.
9. As a data engineer, I want the Rho Matrix computed **in parallel over trait
   pairs** so that a 2,514-Analysis store finishes in minutes, not hours.
10. As a data engineer, I want the Rho Matrix written into the existing store
    (like the top-hit index), with provenance attrs recording `z_thresh`,
    `min_nulls`, `window_bp`, `n_variants_used`, `observed_only`, and `method`,
    so that the estimate is reproducible and auditable.
11. As a data engineer, I want `opengwasdb validate` to check the Rho Matrix
    when present (shape, packed length, support/NaN consistency, provenance) so
    that a corrupt or mismatched matrix is caught.
12. As a developer, I want the estimator to reuse the pleiodb CML logic
    (`estimate_rho_cml`, truncated bivariate-normal normalising grid) verbatim
    so that results match the reference implementation and are not re-derived.
13. As a developer, I want Rho stored as a packed **strict lower triangle**
    (diagonal implied = 1, symmetric) in float16, with an int32 support triangle,
    so that a 2,514-Analysis matrix is a few MB.
14. As a developer, I want the ~200k approximately-independent variants selected
    by **first-variant-per-fixed-bp-window** distance thinning of the store's own
    axis so that no external pruned variant set or LD panel is required.

## Implementation Decisions

### Artifact and trigger

Dense-only, opt-in. A new `build_dense_rho(store_path, *, window_bp, z_thresh,
min_nulls, n_workers)` function and a `build-dense-rho` CLI command write a `rho`
group into the existing `data.zarr` of a built Dense store — the same
add-in-place pattern as the top-hit index. Ragged stores are out of scope (no
shared variant axis).

### Independent variant selection

Distance-thin the store's own variant axis: partition each chromosome into
fixed-width windows (`window_bp`, default chosen to yield ≈200k variants
genome-wide, e.g. ~15 kb) and keep the **first** variant in each window.
Deterministic, needs no external file and no extra full-matrix pass. The chosen
`variant_index` array is stored as provenance.

### Estimation inputs

**Observed Z only.** For a Reference-Completed store, imputed cells
(`imputed == 1`) are treated as missing. For an Analysis pair (j, k), the
estimator uses variants where both are observed **and** both are non-significant:
`isfinite(z_j) & isfinite(z_k) & |z_j| < z_thresh & |z_k| < z_thresh`.

### Estimator (ported from pleiodb `rho.py`)

`estimate_rho_cml(z_j, z_k, z_thresh=1.0, min_nulls=500)`:
- Sufficient statistics on the both-null subset: `A = Σ z_j²`, `B = Σ z_j z_k`,
  `C = Σ z_k²`, `n`.
- If `n < min_nulls` → **NaN**.
- Otherwise minimise the truncated-bivariate-normal negative log-likelihood over
  ρ ∈ (−1 + 1e-6, 1 − 1e-6) via `minimize_scalar` (bounded), using a precomputed
  grid (`_GRID_N = 2000`) of the normalising constant `P(|X| < z_thresh,
  |Y| < z_thresh; ρ)` keyed on `z_thresh`.

Defaults match pleiodb (`z_thresh = 1.0`, `min_nulls = 500`); all are configurable.

### Compute strategy and parallelism

Load the observed, nulls-zeroed Z (`Z`, NaN/significant → 0) and the null mask
`M` (`isfinite & observed & |z| < z_thresh`) at the ~200k thinned variants as
`(n_variants × n_analyses)` arrays. The pair sufficient statistics `A, B, C, n`
are computed with batched linear algebra (they reduce to masked dot products);
the **per-pair MLE is parallelised over pairs** with a `ProcessPoolExecutor`
(`n_workers`), mirroring the single-level pool of ADR 0023. Parallel and serial
runs must produce identical matrices.

### Storage

A `data.zarr/rho` group with:
- `rho` — float16, packed **strict lower triangle** (row-major, `i > j`), length
  `n_analyses (n_analyses − 1) / 2`. Diagonal is implied = 1.0.
- `n_null` — int32, same packing: shared null-variant support per pair.
- Group attrs: `z_thresh`, `min_nulls`, `grid_n`, `window_bp`,
  `n_variants_used`, `observed_only`, `method = "pleiodb-cml"`, `n_analyses`.
- `variant_index` — int32 array of the thinned variants used (provenance).

### Query API (facade)

Symmetric access, diagonal Rho = 1.0, each accessor returns Rho and its support:
- `rho(*ids)` — **long format**: all unique pairs among the given Analysis IDs
  (positional or a single iterable), `{analysis_id_a, analysis_id_b, rho,
  n_null}` arrays; self-pairs excluded. Two IDs → one row.
- `rho_row(analysis_id)` — one Analysis vs all others (vectors + support).
- `rho_matrix(ids=None)` — **wide format**: full symmetric matrix, or the dense
  submatrix block for a vector of IDs.

### Validation

`validate_store` checks the `rho` group when present: packed lengths match
`n_analyses`, `rho` finite ⇔ `n_null ≥ min_nulls`, `rho` within [−1, 1] where
finite, and required provenance attrs present.

## Out of Scope

- Ragged stores (no shared variant axis).
- External pruned variant sets or LD-panel-based pruning (distance thinning only).
- Cross-store Rho or Rho between Analyses in different releases.
- Recomputing Rho incrementally when a store gains Analyses.

## Known Limitation

Compute and storage are O(n_analyses²). This is fine at ukb-b scale (2,514
Analyses → ~3.16M pairs, ~6 MB) but does not extend to stores with tens of
thousands of Analyses without a different design; documented, not engineered
for now.
