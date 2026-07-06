## Context

`StoreQuery.range_phewas()` regressed on the reference-completed ukb-b store:

- observed-only: `620.663 ms`, `5,231,847` associations
- reference-completed: `1,996.685 ms`, `9,669,654` associations

The completed result is larger (`1.85x`) because the same region contains
`4,437,807` additional imputed associations. That explains part, but not all,
of the `~3.2x` latency increase.

The completed variant axis is not fragmented. In the benchmark region
`19:44,500,000-45,500,000`, the completed store returns `4,197` contiguous
variant rows with monotonic positions. New LD-panel variants are physically
interspersed with original variants in genomic order.

The extra avoidable cost is in status labelling:

```python
imp = self._imputed_pairs(rows, cols)
```

For completed stores this performs millions of point reads into
`data.zarr/imputed` via `vindex`, even though regional queries already read
contiguous row blocks from `z` and `se`.

## What to change

- [ ] In dense `range_phewas()`, read the `imputed` matrix as the same row block
      as `z` and `se`, then apply the finite-value mask to derive per-cell
      statuses.
- [ ] When `range_indices()` returns a contiguous interval, use direct row
      slicing (`r0:r1`) rather than orthogonal indexing for `z`, `se`, and
      `imputed`.
- [ ] Keep a fallback for non-contiguous row index arrays so the query remains
      layout-independent.
- [ ] Preserve `observed_only=True` semantics by filtering on the block-read
      imputed flags.
- [ ] Add a regression test proving completed dense `range_phewas()` no longer
      calls `_imputed_pairs()`.

## Acceptance criteria

- [ ] Completed dense regional queries return the same `variant_index`,
      `analysis_index`, `z`, `se`, and `association_status` arrays as the
      current implementation.
- [ ] Completed dense regional queries avoid millions of
      `data.zarr/imputed.vindex[...]` point reads in the hot path.
- [ ] Observed-only stores still return only `"observed"` statuses without
      requiring an `imputed` matrix.
- [ ] `observed_only=True` excludes imputed regional associations.
- [ ] Re-running the ukb-b completed benchmark shows regional latency closer to
      output-size scaling rather than the previous `~3.2x` regression.

## Notes

- This is a query-layer read-path optimisation only. It does not change store
  format or completion output.
- The observed slowdown is not caused by imputed variants being appended out of
  genomic order; they are already interspersed correctly.
