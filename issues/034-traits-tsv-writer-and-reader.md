## Parent PRD

`issues/prd-ragged-observed-only.md`

## What to build

`traits.tsv.gz` format: a bgzipped, optionally tabix-indexed TSV that records one row per Analysis, analogous to `variants.tsv.gz` for the variant axis.

Implement a writer that produces the file and a reader that supports:
- lookup by `analysis_index` (direct row access)
- lookup by `probe_id` (linear scan or SQLite)
- regional range query by `probe_chr` / `probe_bp` via tabix (when coordinates are present)

Schema (tab-separated, `#` header):

```
#probe_chr  probe_bp  analysis_index  analysis_id  probe_id  n  gene_id  gene_name  tissue  context
```

`probe_chr` and `probe_bp` are `NA` for GWAS traits; in that case the file is bgzipped but no `.tbi` index is created. When coordinates are present, the file is sorted by `probe_chr` / `probe_bp` and tabix-indexed identically to `variants.tsv.gz`.

The writer lives in `opengwasdb/traits/axis.py` (new module, parallel to `opengwasdb/variants/axis.py`). The module exposes a `TraitRecord` dataclass, `write_traits_axis()`, and `TraitsAxisReader`.

## Acceptance criteria

- [ ] `TraitRecord` dataclass with all nine non-header fields; `probe_chr` and `probe_bp` are `Optional`
- [ ] `write_traits_axis()` writes bgzipped TSV sorted by `probe_chr` / `probe_bp` when coordinates present; falls back to `analysis_index` order when coordinates absent
- [ ] Tabix index (`.tbi`) created when all records have `probe_chr` / `probe_bp`; skipped otherwise
- [ ] `TraitsAxisReader.range(chr, start, end)` returns matching `TraitRecord` list using tabix
- [ ] `TraitsAxisReader.by_probe_id(probe_id)` returns matching records (via linear scan or SQLite-backed index — implementation choice)
- [ ] Round-trip test: write 5 records covering two chromosomes, query by region and by probe_id, assert exact match
- [ ] File is valid bgzip (can be opened with `pysam.TabixFile`)

## Blocked by

None — can start immediately.

## User stories addressed

- User story 6
- User story 8
- User story 12
