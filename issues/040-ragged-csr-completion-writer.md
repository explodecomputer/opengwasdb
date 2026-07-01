## Parent PRD

`issues/prd-ragged-reference-completion.md`

## What to build

Implement `opengwasdb/layouts/ragged/complete.py`: the per-analysis CSR assembly loop that reads an observed-only ragged store and writes a reference-completed one. This is the core of the enhancement pipeline.

The function `complete_ragged_store(source_path, dest_path, ld_dir, ancestry, cis_window_bp, min_cor, thresh, store_id, release_id, overwrite)` should:

1. Open the source observed-only store (manifest, RaggedCSRReader, TraitsAxisReader, VariantAxis).
2. Collect all LD reference panel variants intersecting any analysis's cis window; merge with observed variant set; assign new store-local variant indices; write extended `variants.tsv.gz` and ALID sidecar.
3. For each analysis (indexed by CSR offset):
   - Extract observed z+se from the source CSR.
   - Find LD blocks intersecting cis window; call `match_variants`, `load_ld_eigenvectors`.
   - For each matched block: call `impute_z_block`; if pearson_r ≥ min_cor, fill imputed z and SE via scalar-N; record quality row.
   - Merge observed + reference-panel associations sorted by variant_index.
   - Append to new CSR arrays (offsets, variant_index, z, se, imputed).
4. Write `index.sqlite` with `analyses` table (copied) and `completion_quality` table.
5. Write `traits.tsv.gz` (copied from source).
6. Write `manifest.json` with `completion_state: "reference_completed"`, `ld_panel_id`, `reference_completion_method`, `cis_window_bp`.
7. Build top-hit index and validate.

The new `imputed` array is a uint8 zarr array alongside `z` and `se` in `data.zarr/ragged/`.

## Acceptance criteria

- [ ] `complete_ragged_store` produces a store that passes `validate_store`
- [ ] Observed associations in the completed store have `imputed=0` and match source z+se exactly
- [ ] Imputed associations have `imputed=1`, finite z+se, and p-value ≤ 5e-8 is not required (any value is fine)
- [ ] Reference-panel variants with failed imputation have `imputed=0`, NaN z, NaN se
- [ ] `completion_quality` table has one row per (analysis_index, block_id); `pearson_r` is NULL for non-imputed blocks
- [ ] `manifest.json` `completion_state` is `"reference_completed"`
- [ ] Calling with an existing dest_path raises `FileExistsError` unless `overwrite=True`
- [ ] Integration test: build tiny observed-only store from BESD fixture + synthetic LD panel → complete → validate
- [ ] `opengwasdb/model/enums.py` has `CompletionState.REFERENCE_COMPLETED`

## Blocked by

- `issues/038-ragged-imputation-kernel.md`
- `issues/039-ld-panel-block-reader.md`

## User stories addressed

- User story 3 (NaN rows for missing reference variants)
- User story 7 (new release, source unchanged)
- User story 8 (completion_quality table)
- User story 9 (variant table includes all LD panel variants used)
- User story 11 (imputed mask uint8 array)
- User story 12 (scalar-N SE)
