## Context

`validate_store` on a Dense store is now the worst-behaved step in the pipeline —
**worse than a single matrix load**. Measured on the real ukb-b build+validate
(9,847,701 × 2,514): it ran **80+ minutes at 340+ GB RSS and was still climbing**
when it had to be killed (shared box hit 709 GB used). The build itself peaks at
~43 GB (issue 043), so validation is ~8× the build's footprint and never finished.

Two compounding causes, both in `opengwasdb/validation/validate.py`:

1. **Stacked full-matrix float32 loads.** For an Observed-Only dense store,
   `_validate_dense_store` calls `_validate_dense_arrays` (loads `z` + `se` as
   float32 — validate.py:474-475, ~198 GB) and then `_validate_top_hits` which
   loads the **whole matrix again** (`z_matrix = root["z"][:].astype("float32")`,
   validate.py:491, +99 GB). The first arrays are not freed before the second
   load, so RSS stacks toward ~300 GB+.

2. **Pure-Python loop over every top-hit cell (the real killer).**
   `_validate_top_hits` (validate.py:~500-505) iterates the top-hit index entries
   and does `stored_z = float(z_matrix[int(row), int(col)])` per cell to check it
   against the matrix. At genome scale that is **~58M iterations** (41.9M + 11M +
   5.5M across the three tiers), single-cell Python indexing into a 99 GB array —
   the same per-element-loop antipattern that made the old top-hit *build*
   pathological. This is what makes it take well over an hour.

## What to change

- [x] **Vectorise the top-hit cross-check.** Done in `794e756`. `_validate_top_hits`
      now runs one streamed band pass with a vectorised `z_band[rows-r0, cols]`
      gather + `np.isclose` per band instead of ~58M per-cell Python iterations.
- [x] **Don't reload the matrix per sub-check / free between checks.** Done. The
      observed-array checks and the top-hit cross-check each stream `z`/`se` in
      row-bands; the completion checks (below) now share the `_validate_dense_arrays`
      band loop so z/se are read once for a completed store, not twice.
- [x] **Stream the dense finite/NaN checks in bands** and avoid the float32 upcast.
      Done in `794e756` for observed-only stores (finite/sign checks run on float16
      directly). The Reference-Completed path (`_validate_dense_completion`) was
      **missed by `794e756`** — it still loaded full `imputed` + `z`/`se` float32
      (~250 GB) and is the path a completed ukb-b store actually hits. Now fixed:
      `_validate_completion_metadata` does only the cheap structural/table checks and
      hands the `imputed` (lazy) + `on_panel` (1-D) arrays to `_validate_dense_arrays`,
      whose band loop folds in the imputed-0/1, imputed-cell-finite-z/se, and
      off-panel-never-imputed checks — peak = one band, not the full matrices.

## Acceptance criteria

- [x] `validate_store` on the ukb-b store completes in minutes at ~one band of
      memory, not 80+ min at 300+ GB. (Observed-only path fixed in `794e756`; the
      completed path — the remaining ~250 GB balloon — is fixed here. No genome-scale
      completed store built yet to time end-to-end, but no full matrix is resident in
      any dense sub-check.)
- [x] Validation results unchanged on existing stores (dense-validation tests stay
      green); a completion cross-check test locks in the streamed path — a NaN in an
      imputed cell (`test_corrupt_imputed_nan_z_fails`) is caught by the band loop, in
      addition to the existing off-panel-imputed and valid-store tests.

## Notes

- The ukb-b store built on branch `dense-pass2-memory-streaming` is confirmed
  correct by direct query checks (shapes, per-tier top-hit counts + min|z|, phewas
  / analysis / top_hits queries) — the formal `validate_store` was killed only
  because it malfunctions at scale, not because anything is wrong with the store.
- Sibling of [[043-dense-vcf-build-memory-footprint]] and
  [[044-dense-completion-memory-footprint]] — same whole-matrix-in-RAM class plus a
  per-element loop, on the validate side.
- The ragged validators (`_validate_ragged_store` / `_validate_ragged_completion`)
  still load their full `z`/`se`, but those are 1-D `(n_assoc,)` arrays — the actual
  observed data, not a dense `(n_variants × n_analyses)` matrix — so they are O(data)
  by nature and out of scope for this dense-matrix issue.
