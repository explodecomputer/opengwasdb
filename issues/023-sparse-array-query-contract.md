## Parent PRD

`issues/prd-query-array-api.md`

## What to build

Replace the five public query methods on `StoreQuery` (`analysis`, `phewas`, `range`, `lookup`, `top_hits`) so they all return `dict[str, np.ndarray]` with four parallel sparse flat arrays:

- `variant_index` — int32, store-local row index
- `analysis_index` — int32, store-local column index
- `z` — float32
- `se` — float32

All four arrays have the same length `k` = number of finite cells in the result. Only finite (non-NaN) z and se values are included.

Add two public metadata accessors to `StoreQuery`:

- `variants_table()` → `dict[int, dict]`, keyed by `variant_index`, each dict containing `alid`, `chromosome`, `position`, `effect_allele`, `other_allele`, `rsid`
- `analyses_table()` → `dict[int, dict]`, keyed by `analysis_index`, each dict containing `analysis_id`, `phenotype_id`, `phenotype_label`, `analysis_label`, `stored_effect_scale`

Remove `AssociationResult` (and its `from_rows` factory and `to_dict`) entirely. Remove `analysis_arrays()` — it is superseded by `analysis()`.

The five query methods are implemented as follows:

- **`analysis(analysis_id)`** — reads one zarr column, masks finite cells, builds sparse arrays from non-NaN positions
- **`phewas(identifier)`** — reads one zarr row, same approach
- **`range(chromosome, start, end)`** — fetches variant indices from SQLite, reads the row sub-block from zarr via `oindex`, flattens to sparse arrays
- **`lookup(identifiers, analysis_ids)`** — same sub-block approach as range
- **`top_hits(threshold, limit)`** — reads the pre-indexed zarr `top_hits/{key}` arrays directly; no per-row SQLite lookup; variant and analysis metadata are not loaded inline

All query tests in `test_dense_vertical_slice.py` and `test_dense_vcf_build.py` are rewritten to assert against the sparse array contract. The `conftest.py` fixture `dense_store_path` is unchanged.

## Acceptance criteria

- [ ] All five query methods return `dict` with keys `variant_index`, `analysis_index`, `z`, `se` and no other keys.
- [ ] `variant_index` and `analysis_index` are dtype int32; `z` and `se` are dtype float32.
- [ ] All four arrays have equal length; length equals the number of finite associations in the result.
- [ ] `variants_table()` returns a dict keyed by variant index with correct ALID and coordinate fields.
- [ ] `analyses_table()` returns a dict keyed by analysis index with correct analysis metadata.
- [ ] `top_hits()` does not issue per-row SQLite queries for variant or analysis metadata.
- [ ] `AssociationResult` class is removed from the codebase (no import, no reference).
- [ ] `analysis_arrays()` method is removed.
- [ ] All tests in `test_dense_vertical_slice.py` pass against the new contract.
- [ ] All query-related tests in `test_dense_vcf_build.py` pass against the new contract.
- [ ] The full test suite (`pytest tests/`) passes.

## Blocked by

None — can start immediately.

## User stories addressed

- User stories 1–10 from `issues/prd-query-array-api.md`
