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

- [ ] **Allocate `z_mat`/`se_mat` after the worker pool has forked** (or otherwise
      keep them out of the workers' inherited address space). Workers don't touch
      them, so this removes the COW-doubling with no behavioural change. Low-risk,
      high-value (~99 GB).
- [ ] **Represent the Pass 2 lookups as fork-safe structures** so reading them in a
      worker doesn't trigger refcount-COW. Options: pack `variant_index` as sorted
      numpy arrays + binary search (no per-element Python objects → no per-element
      refcounts), and/or resolve `(chrom,pos,ref,alt) → row` via the existing ALID
      byte index rather than a Python dict. More involved (~100+ GB).

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
