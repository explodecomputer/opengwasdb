## Context

`StoreQuery.top_hits()` is much slower on the reference-completed ukb-b store
than on the observed-only store:

- observed-only: `618.734 ms`, `5,493,678` hits
- reference-completed: `6,026.986 ms`, `6,441,285` hits

The result set only grew by `947,607` cells (`1.17x`), so cardinality does not
explain the `~9.7x` latency regression.

Direct timing shows the dense top-hit index itself is still fast to read. The
slow path is the completed-store status labelling:

```python
imp = self._imputed_pairs(variant_indices, analysis_indices)
```

For completed stores this performs millions of point reads into
`data.zarr/imputed` via `vindex` so `top_hits()` can label each hit as
`observed` or `imputed`. On the completed ukb-b benchmark this lookup took
about `5.0s` for `6.4M` top-hit cells. The sum of those imputed flags was
`947,607`, exactly matching the top-hit count increase over observed-only.

The top-hit index already stores compact contiguous arrays for
`variant_index`, `analysis_index`, `z`, `se`, `abs_z`, and `p_value`. The
imputed flag is the only top-hit result field that still requires random access
back into the dense matrix hierarchy.

## What to change

- [ ] Add an optional `imputed` array to each dense top-hit index group:
      `data.zarr/top_hits/<threshold>/imputed`, `uint8`, same length and chunking
      as `variant_index`.
- [ ] When building top-hit indexes for a reference-completed dense store, write
      the imputed flag for every indexed cell alongside `z` and `se`.
- [ ] Keep observed-only stores compact: either omit the `imputed` index array or
      write all zeros only if that materially simplifies the code.
- [ ] Update `StoreQuery.top_hits()` so it reads `group["imputed"][:]` when
      present and falls back to `_imputed_pairs(...)` for older completed stores
      that do not yet have the denormalised field.
- [ ] Preserve `observed_only=True` semantics by filtering on the indexed
      imputed flags when available.
- [ ] Update dense validation so that, when a top-hit `imputed` array exists, it
      checks length consistency and verifies sampled/indexed values against
      `data.zarr/imputed`.
- [ ] Provide a rebuild path for existing completed stores so top-hit indexes can
      be regenerated without rebuilding the whole reference-completed store.

## Acceptance criteria

- [ ] `top_hits(threshold=5e-8)` on the completed ukb-b store no longer performs
      millions of `data.zarr/imputed.vindex[...]` reads in the hot path.
- [ ] Completed-store `top_hits()` returns byte-identical `variant_index`,
      `analysis_index`, `z`, `se`, and `association_status` arrays before and
      after the index change.
- [ ] `top_hits(observed_only=True)` returns the same observed subset before and
      after the index change.
- [ ] Older completed stores without `top_hits/<threshold>/imputed` remain
      readable via the current fallback path.
- [ ] A small reference-completed dense fixture tests both paths: with indexed
      imputed flags and without them.
- [ ] Re-running the ukb-b completed benchmark shows top-hit latency close to
      contiguous index-read cost rather than the current `~6s` random-mask-read
      cost.

## Notes

- This is a small denormalisation in the query index, not a change to the core
  association matrices.
- Storage cost should be tiny relative to the existing top-hit index: one
  `uint8` per top-hit cell before compression.
- This is distinct from top-hit membership correctness
  (`046-top-hit-index-float16-consistency.md`) and from dense lookup surgical
  reads (`052-dense-lookup-surgical-column-read.md`). Here the problem is
  repeated random status lookup for already-indexed top-hit cells.
