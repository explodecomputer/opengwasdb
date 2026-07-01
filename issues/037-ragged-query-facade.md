## Parent PRD

`issues/prd-ragged-observed-only.md`

## What to build

The query facade for Ragged Observed-Only stores, exposing the same interface as the Dense query facade.

Implement `opengwasdb/layouts/ragged/query.py` with:

```python
def query_analysis(store, analysis_id: str) -> AssociationResult
    """All associations for one analysis (probe_id or analysis_index lookup)."""

def query_region_by_probe(store, chr: str, start: int, end: int) -> AssociationResult
    """All associations where the probe/gene falls in the given region."""

def query_region_by_variant(store, chr: str, start: int, end: int) -> AssociationResult
    """All associations where the variant falls in the given region."""

def query_tophits(store, threshold: float) -> AssociationResult
    """Associations passing a significance threshold, across all analyses."""
```

### Query implementations

- `query_analysis`: resolve `probe_id` → `analysis_index` via SQLite `analyses` table; read CSR slice via `RaggedCSRReader.get_analysis()`; join to variant metadata via `VariantAxisReader`
- `query_region_by_probe`: tabix range on `traits.tsv.gz` → list of `analysis_index` values → CSR slices → flatten; join to variant metadata
- `query_region_by_variant`: tabix range on `variants.tsv.gz` → list of `variant_index` values; for each CSR analysis, filter to matching variant indices (deferred until variant-centric CSR is added; for now: load all analyses in probe region and filter)
- `query_tophits`: scan CSR arrays for |z| above threshold; return top hits with analysis and variant metadata

The facade is opened via the same `open_store()` entry point as Dense, dispatching on `primary_layout` from `manifest.json`.

## Acceptance criteria

- [ ] `query_analysis("ENSG00000000003")` returns correct (variant_index, z, se, chr, pos, alid) rows
- [ ] `query_region_by_probe("1", 1_000_000, 2_000_000)` returns associations for all probes with TSS in that window
- [ ] `query_region_by_variant("1", 1_000_000, 2_000_000)` returns associations where the variant is in that window
- [ ] `query_tophits(threshold=5e-8)` returns associations with p-value below threshold
- [ ] All queries return `AssociationResult` with `analysis_id`, `probe_id`, `alid`, `z`, `se`, `pval` fields
- [ ] `open_store()` dispatches correctly on `primary_layout: "ragged"`
- [ ] Integration test using the eqtlgen smoke test store from `issues/036`

## Blocked by

- `issues/036-ragged-build-besd.md`

## User stories addressed

- User story 1
- User story 2
- User story 3
- User story 4
