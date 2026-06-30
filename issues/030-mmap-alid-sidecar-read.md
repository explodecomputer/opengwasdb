## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Wire the sidecar files into `VariantAxis` for fast identifier resolution:

1. In `VariantAxis.__init__`, memory-map `variant_alid_bytes.npy` and
   `variant_alid_rows.npy` using `np.load(..., mmap_mode="r")`.

2. Replace the tabix-per-call path in `by_identifier()` for canonical ALID-format
   identifiers with a `np.searchsorted` lookup on the mmap'd byte array. Return the
   matching `VariantRecord` via the existing `by_index()` (offset sidecar path).

3. Update `by_alid()` to use the mmap index directly instead of calling `range()`.

4. rsid / alias identifiers continue to fall back to SQLite as before.

The result: a single `by_identifier()` call for a canonical ALID drops from ~0.15 ms
(pysam tabix fetch) to ~5 µs (`np.searchsorted` on a mmap'd array).

## Acceptance criteria

- [ ] `by_identifier()` with a canonical ALID uses the mmap index (no pysam fetch)
- [ ] `by_identifier()` with an rsid still resolves via SQLite alias table
- [ ] `by_identifier()` returns `None` for an unknown ALID (no error)
- [ ] `phewas()` and `lookup()` query results are identical before and after this change
- [ ] All existing tests pass

## Blocked by

- `issues/029-mmap-alid-sidecar-write.md`

## User stories addressed

- User story 2 (sub-millisecond identifier resolution)
- User story 3 (OS page cache handles warmup without explicit preload)
