## Context

`StoreQuery.lookup` reads the **full analysis width** for its selected variant
rows and then discards the unwanted columns in numpy
(`opengwasdb/query/facade.py:197-198`):

```python
z_block = self._root["z"].oindex[row_indices, :][:, col_indices].astype("float32")
se_block = self._root["se"].oindex[row_indices, :][:, col_indices].astype("float32")
```

The `[row_indices, :]` slice pulls **every** column-chunk for each selected row,
regardless of how few analyses were requested. The chunk is the atomic
read+decompress unit in zarr, so the floor for a scattered lookup is "read every
chunk that contains a requested cell"; this code reads well above that floor on
the column axis.

Under the old square `(1000, 1000)` chunks this was free — only 3 column-chunks
existed, and a handful of requested analyses hit all of them anyway, so
full-width *was* the surgical set. With a narrow analysis chunk (`(1000, 128)` →
20 column-chunks) it forces reading all 20 column-chunks per row when the
requested analyses actually touch only ~8, which is a large part of why
`random_lookup` is 2.4× slower on the 1000×128 store (825 ms → 1981 ms; see
`docs/benchmark-output/opengwasdb_chunk_comparison_benchmark.*`).

This is a query-layer inefficiency independent of chunk shape: an orthogonal fix
that helps point lookups under any layout, and specifically de-risks choosing a
narrow analysis chunk.

## What to change

- [ ] In `lookup`, replace the full-width read + numpy subset with orthogonal
      indexing on **both** axes so only column-chunks intersecting the requested
      analyses are read — e.g. `oindex[row_indices, col_indices]` (outer-product
      / orthogonal selection, not `vindex` paired coordinates). Results must be
      identical; this is a pure read-path change, no format change.
- [ ] Audit the other full-width gather at `facade.py:240` in `top_hits`
      (`oindex[unique_rows, :]` then column subset) and apply the same surgical
      read if it over-reads for the same reason.
- [ ] Leave genuinely full-width reads alone: `analysis` (`[:, col]`, one whole
      column — bulk), `range_phewas`/`phewas` (one variant row × all analyses —
      the full width *is* the requested set).

## Acceptance criteria

- [ ] `lookup` returns byte-identical `z` / `se` / `association_status` /
      `variant_index` / `analysis_index` arrays to the current implementation,
      including scattered rows × columns and cells that are missing/imputed, on a
      small store and on the real ukb-b store (both chunk shapes).
- [ ] `lookup` on the `(1000, 128)` store reads fewer column-chunks and is
      measurably faster than the current full-width path, with no regression on
      the `(1000, 1000)` store.
- [ ] A test covers `lookup` correctness for scattered `rows × cols` including
      missing cells and `observed_only=True`.
- [ ] The chunk-comparison benchmark's `random_lookup` row is re-measured and the
      report notes the post-fix numbers.

## Notes

- The row axis is not the issue: 100 scattered variants each sit in their own
  1000-tall chunk, so ≥100 row-chunks are read whole regardless — inherent to
  chunked-compressed storage, out of scope here.
- Confirm zarr's orthogonal-indexing path (`oindex` with two index arrays) is not
  itself slower per element than the full-width read on square chunks; if it is,
  gate the surgical path on the column-chunk count (only worthwhile when the
  analysis axis has many chunks).
- Sibling context: the chunk-shape trade-off is a separate lever
  (`benchmarks/benchmark_chunk_comparison.py`); this issue is the orthogonal
  query-side optimisation flagged there.
