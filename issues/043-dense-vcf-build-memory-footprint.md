## Context

The parallel two-pass Dense VCF builder (`opengwasdb/layouts/dense/build_vcf.py`,
see ADR-0024) works correctly and completes the full ukb-b build (9.85M variants
× 2514 analyses) in ~4h. But Pass 2 is far more memory-intensive than the ~99 GB
you would expect from the dense `z`/`se` matrices alone.

Measured on the live ukb-b build (32 workers, mid Pass 2), kernel counters:

- `AnonPages: ~492 GB` — total anonymous RAM across the process group
- `z` + `se` matrices: 9.85M × 2514 × 2 bytes × 2 arrays = **99 GB logical**
- `free`: ~498 GB used, ~515 GB available, **swap untouched** (fits the 1 TB box
  with headroom, so this is not currently blocking — but it is ~5× the logical
  matrix size, and would become the binding constraint at larger analysis counts
  or with more workers)

Two fork-related effects inflate the footprint, neither of which is fundamental:

1. **COW-doubling of the output matrices (~+99 GB).** `z_mat`/`se_mat` are
   allocated in the parent *before* the worker pool forks, so all 32 children
   inherit them copy-on-write. The parent then scatters writes into them during
   Pass 2; every written page COW-faults into a fresh private parent copy while
   the children still pin the untouched originals — so the arrays approach ~2×
   their logical size in physical RAM. The workers never read or write
   `z_mat`/`se_mat` at all, so inheriting them is pure waste.

2. **Refcount-COW duplication of the lookup dicts (~+100 GB).** `hg19_lookup`
   and `variant_index` (~9.85M entries each) are inherited COW and meant to be
   read-only in workers. But CPython reference counting *writes* to an object's
   header on every access, so a worker calling `hg19_lookup.get(...)` across the
   dataset progressively private-copies most of the dict's pages. With 32 workers
   this duplicates the lookups many times over.

## What to change

Both fixes are independent and can land separately.

- [x] **Keep `z_mat`/`se_mat` out of the workers' inherited COW space.** Done, then
      superseded: an interim MAP_SHARED memmap removed the COW-doubling; the final
      design (below) removes the full matrix entirely. The build now **band-streams
      the zarr** — Pass 2 spills each analysis column to disk and holds no matrix;
      a post-Pass-2 band-write phase fills the zarr in chunk-column bands (peak =
      one band ≈ `n_variants × chunk-analysis-width`, ~40 GB at ukb-b scale), so the
      99 GB matrix is never resident alongside the workers and COW-doubling is moot.
      (`_create_dense_zarr` + `_write_dense_bands` in `build_vcf.py`.)
- [x] **Represent the Pass 2 lookups as fork-safe structures.** Done: the two dicts
      are composed into a single sorted numpy byte-key array (`chrom:pos:ref:alt`)
      + int32 `rows`, binary-searched in workers (`_build_variant_key_index` /
      `_resolve_column`). Numpy buffers have no per-element refcounts, so forked
      workers reading them don't refcount-COW — the ~100 GB per-worker dict
      duplication is gone. Output byte-identical
      (`test_two_workers_matches_serial`, `test_short_final_band_matches_single_band`,
      `test_last_wins_dedup_on_collision`, `test_absent_variant_not_mismapped`).

## Acceptance criteria

- [ ] A full-scale build (or the chr1×1000 proxy) shows `AnonPages` peak close to
      the logical matrix size (~99 GB for ukb-b) rather than ~5×
- [ ] Output store is byte-identical to the current builder (validate parity the
      way `test_two_workers_matches_serial` does)
- [ ] Peak-RSS assertion or a documented measurement added to the benchmark so
      regressions are visible

## Notes

- Worker transient memory (each file's `last_by_row` dict, ~2 GB) is bounded and
  fine; the disk-spill design already keeps per-file *results* off the parent's
  heap (spill dir stayed at ~289 MB during the run).
- Not urgent while the target hardware has ~1 TB RAM, but this is the reason a
  bigger analysis count would hit a memory wall before a CPU one. Related:
  [[ADR-0024]] documents the build design these fixes apply to.
