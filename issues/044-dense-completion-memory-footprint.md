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

- [x] **Band-stream the completion write.** Done: `_create_completed_zarr` creates
      empty `z`/`se`/`imputed`/`on_panel` (fill_value), and `_write_dense_bands`
      seeds from the source, applies the sorted imputed fills, and writes
      z/se/imputed one row-band at a time — the full matrix is never resident. The
      analyses-table write (needs `n_missing_off_panel`) is reordered to after the
      band pass.
- [x] **Don't hold the source as float32.** Done: the full `src_z`/`src_se` load
      is gone; each band reads only the source rows it needs
      (`src_z.oindex[srows]`) via an `out_to_src` inverse map.

Validated on chr1×100 (869,476 × 100, 34.6M imputed): parent peak ~8.6 GB (fills +
block results), band write never holds the matrix; store validates. Two adjacent
fixes landed with it: a region-based imputation z-cap (QC, per pleiodb) and
canonicalisation of the real LD panel's `chr:pos_ref_alt` ALIDs (the completion
otherwise matched nothing), plus a BLAS thread cap so the LD-block pool isn't
oversubscribed (~50 min → ~2 min).

## Follow-up (not blocking)
- [x] **Resolved.** The accumulated `fills` are no longer held as Python tuples in
      the parent. Workers checkpoint fills to disk and return an empty-fills marker
      (no fills over the pool result queue); Phase 3 reads each checkpoint back as
      packed numpy arrays and resolves fill ALIDs → union rows with one vectorised
      `searchsorted`, so parent fill memory is packed arrays (~12–16 B/cell), not
      ~O(n_imputed) tuples. Verified on chr1×100: parent peak 8.6 GB → 3.03 GB,
      same n_imputed, same time, store validates and is byte-identical across two
      runs. Parity tests now assert z/se **value** equality (not just cell counts),
      so the streaming path is guarded byte-for-byte against the serial path.

## Acceptance criteria

- [x] A completed store at genome scale (or a scaled proxy) peaks well below the
      current ~4× matrix footprint — target ≈ one band + block-worker memory.
      (chr1×100 proxy: parent peak 3.03 GB; band-write never holds the matrix.)
- [x] Output store byte-identical to the current completion — now guarded by the
      dense-completion parity tests, which assert z/se/imputed/on_panel **value**
      equality (not just counts) in `test_resume_matches_fresh_run` /
      `test_two_workers_matches_serial` (`tests/test_dense_completion.py`).

## Notes

- Sibling of [[043-dense-vcf-build-memory-footprint]]; the band-streaming machinery
  built there should be factored so completion can reuse it rather than
  reimplementing.
- Lower urgency than 043 while no genome-scale Dense completion is being run, but it
  is the same class of bug and will block Dense reference completion of ukb-b-scale
  stores.
