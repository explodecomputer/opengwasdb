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

- [ ] **Vectorise the top-hit cross-check.** Replace the per-cell Python loop with
      a single `z_matrix[rows, cols]` gather (fancy indexing) compared to the
      stored `z` array in one `np.array_equal`/`allclose`. Turns ~58M Python
      iterations into one vectorised op.
- [ ] **Don't reload the matrix per sub-check / free between checks.** Load `z`
      (and `se`) once, share across `_validate_dense_arrays` and
      `_validate_top_hits`, or better, stream in chunk-row bands.
- [ ] **Stream the dense finite/NaN checks in bands** and avoid the float32 upcast
      where `np.isfinite` on float16 suffices — reuse the `_write_dense_bands`
      chunk-band iteration pattern from the build.

## Acceptance criteria

- [ ] `validate_store` on the ukb-b store completes in minutes at ~one band of
      memory, not 80+ min at 300+ GB.
- [ ] Validation results unchanged on existing stores (dense-validation tests in
      `tests/test_validation.py` stay green); add a top-hit cross-check test on a
      small store to lock in the vectorised path.

## Notes

- The ukb-b store built on branch `dense-pass2-memory-streaming` is confirmed
  correct by direct query checks (shapes, per-tier top-hit counts + min|z|, phewas
  / analysis / top_hits queries) — the formal `validate_store` was killed only
  because it malfunctions at scale, not because anything is wrong with the store.
- Sibling of [[043-dense-vcf-build-memory-footprint]] and
  [[044-dense-completion-memory-footprint]] — same whole-matrix-in-RAM class plus a
  per-element loop, on the validate side. Now the highest-impact remaining item.
