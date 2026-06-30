## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Add `range_indices(chromosome, start, end) -> np.ndarray[int32]` to `VariantAxis`.
It parses only the `variant_index` column from tabix output — no `VariantRecord`
allocation. Wire `StoreQuery.range()` to use this method internally instead of
constructing a list of `VariantRecord` objects and extracting indices from them.

The existing `VariantAxis.range()` method (returning `list[VariantRecord]`) is kept
unchanged for callers that need full metadata.

Profiling baseline: `VariantAxis.range()` for 5 463 variants takes 51.6 ms;
raw tabix fetch takes 2.3 ms. The object construction is the entire gap.

## Acceptance criteria

- [ ] `VariantAxis.range_indices(chrom, start, end)` returns `np.ndarray[int32]` of
      variant indices in ascending order
- [ ] Returns an empty `np.ndarray` for regions with no variants or unknown chromosomes
- [ ] `range_indices` returns the same indices as `[v.variant_index for v in range()]`
      (unit test)
- [ ] `StoreQuery.range()` uses `range_indices()` internally and returns the same sparse
      dict as before
- [ ] All existing tests pass

## Blocked by

None — can start immediately.

## User stories addressed

- User story 1 (responsive regional queries)
- User story 7 (facade never allocates Python variant objects during query execution)
