## Parent PRD

`issues/prd-ragged-reference-completion.md`

## What to build

Extend `RaggedStoreQuery` in `opengwasdb/query/facade.py` to support reference-completed stores. All existing query methods should return imputed associations by default; an `observed_only=True` parameter on each method filters them out. Result rows gain an `association_status` field (`"observed"` | `"imputed"` | `"missing"`).

Changes needed:
- Detect `completion_state == "reference_completed"` from the manifest at store-open time.
- Read the `imputed` array from `data.zarr/ragged/imputed` when present.
- Add `observed_only: bool = False` parameter to `analysis()`, `range_by_analysis()`, `phewas()`, `range_phewas()`, `top_hits()`, `lookup()`.
- When `observed_only=True`, filter out rows where `imputed[i] == 1`.
- Add `association_status` column to all returned result dicts (0→`"observed"`, 1→`"imputed"`; NaN z → `"missing"`).
- For observed-only stores, `imputed` array is absent; `association_status` is always `"observed"`.
- `top_hits` fast path: top-hit index was built from completed data; `observed_only=True` falls back to full-scan with filter.

## Acceptance criteria

- [ ] `query_store(completed_store_path).analysis("ENSG...")` returns both observed and imputed rows
- [ ] `query_store(completed_store_path).analysis("ENSG...", observed_only=True)` returns only observed rows
- [ ] All result dicts include `association_status` key
- [ ] `association_status` is `"observed"` for all rows from an observed-only store (no KeyError)
- [ ] `range_by_analysis` with `observed_only=True` excludes imputed rows
- [ ] NaN z rows from missing reference variants have `association_status == "missing"` and are excluded by `observed_only=True`
- [ ] Tests using the integration fixture built in `issues/040`

## Blocked by

- `issues/040-ragged-csr-completion-writer.md`

## User stories addressed

- User story 1 (retrieve all reference-panel variants in cis window)
- User story 2 (distinguish imputed from observed)
- User story 4 (include imputed by default)
- User story 5 (observed-only mode)
