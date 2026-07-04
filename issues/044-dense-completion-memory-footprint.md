## Context

Dense reference completion (`opengwasdb/layouts/dense/complete.py`,
`complete_dense_store`) materialises the full `(n_variants × n_analyses)` matrices
in RAM the same way the VCF builder used to (issue 043), but worse: it holds
**several float32 matrices at once**.

Observed in `complete_dense_store` (line refs approximate):

- `src_z = src_root["z"][:].astype(np.float32)` and
  `src_se = ...astype(np.float32)` — the full source matrices upcast float16→float32
  (complete.py:431-432).
- `z = np.full((n_variants, n_analyses), np.nan, np.float32)` and matching `se`
  (complete.py:434-435) — the full output matrices.
- `imputed = ...` uint8 full matrix (~complete.py:527).
- Final write via `_write_dense_zarr` (complete.py:~623-651) uses the same one-shot
  `create_dataset(..., data=...)` as the old dense build.

So peak ≈ `4 × (n_variants × n_analyses × 4 bytes)` for z/se (src + output) plus a
uint8 mask. For a completed ukb-b (9.85M × 2514) that is **~396 GB of float32
matrices** (99 GB each) + ~25 GB mask — before the LD-block worker pool's own
footprint. This will not fit / will thrash at genome scale.

Two compounding costs, both avoidable:

1. **float16 → float32 upcast of the whole source** (`src_z`/`src_se`). Doubles the
   source matrices to 198 GB when they are 99 GB on disk; only needed transiently
   per row/block, not for the whole matrix at once.
2. **Whole-matrix-in-RAM then write-once**, holding source and output
   simultaneously — the same pattern issue 043 removes from the VCF builder via
   band-streaming.

## What to change

- [ ] **Band-stream the completion write** the same way as the VCF builder (see the
      band-streaming design landed for issue 043): create the `z`/`se`/`imputed`/
      `on_panel` zarr datasets empty with `fill_value`, and write them in
      analysis-bands (or LD-block row-bands) so the full output matrix is never
      resident.
- [ ] **Avoid upcasting the entire source to float32.** Read source z/se per
      block/row-range as needed rather than `src_root["z"][:].astype(float32)` over
      the whole array; keep working precision local to each block.
- [ ] Consider whether `src_z`/`src_se` need to be resident at all once the LD-block
      workers read their own row slices directly from the source zarr (the workers
      already open the source store).

## Acceptance criteria

- [ ] A completed store at genome scale (or a scaled proxy) peaks well below the
      current ~4× matrix footprint — target ≈ one band + block-worker memory.
- [ ] Output store byte-identical to the current completion (guard with the
      existing dense-completion parity tests, e.g. `test_resume_matches_fresh_run`
      / `test_two_workers_matches_serial` in `tests/test_dense_completion.py`).

## Notes

- Sibling of [[043-dense-vcf-build-memory-footprint]]; the band-streaming machinery
  built there should be factored so completion can reuse it rather than
  reimplementing.
- Lower urgency than 043 while no genome-scale Dense completion is being run, but it
  is the same class of bug and will block Dense reference completion of ukb-b-scale
  stores.
