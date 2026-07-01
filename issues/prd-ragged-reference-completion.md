## Problem Statement

Observed-only ragged stores (e.g. eqtlgen, GTEx) only contain variants that were tested in the source dataset. Downstream analyses such as colocalisation and fine-mapping require statistics at every variant in a reference panel within the cis window — including variants not present in the original summary statistics. Missing reference-panel variants currently require slow, bespoke LD-proxy lookup at query time, which is error-prone and not reproducible across callers.

## Solution

Implement a reference completion enhancement for ragged cis stores. Starting from a built and validated observed-only ragged store, a new Reference-Completed Store Release is produced that extends each analysis's CSR association sequence to cover all LD reference panel variants within the declared cis window (trait_bp ± 1 Mb). Imputed associations are filled using elastic-net regression on LD eigenvectors (the algorithm used in pleiodb), with SE derived from the scalar-N model using LD panel EAF. Reference panel variants that cannot be imputed (block quality below threshold) are retained as NaN rows. The completed store is a new immutable release with the same store identity; the observed-only source remains unchanged.

## User Stories

1. As a bioinformatician, I want to retrieve association statistics at every reference-panel variant in a gene's cis window so that I can run colocalisation without a manual LD-proxy step.
2. As a bioinformatician, I want imputed associations to be distinguishable from observed associations so that I can run observed-only analyses when needed.
3. As a bioinformatician, I want missing reference-panel variants (failed imputation) to appear as NaN rows rather than being silently absent so that query semantics are consistent across completed and observed-only stores.
4. As a bioinformatician, I want the query API to include imputed associations by default so that colocalisation workflows need no changes for reference-completed stores.
5. As a bioinformatician, I want to request observed-only results from a reference-completed store so that I can compare with the original summary statistics.
6. As a data engineer, I want a CLI command `complete-ragged` that takes an observed-only store and an LD panel path and writes a new completed store so that completion can be scripted without writing Python.
7. As a data engineer, I want the enhancement to produce a new release rather than mutating the source store so that the observed-only release remains queryable and immutable.
8. As a data engineer, I want per-block quality metrics (Pearson r, n_imputed, n_missing) stored in SQLite so that I can audit imputation quality without re-running completion.
9. As a data engineer, I want the variant table in the completed store to include all LD reference panel variants used in any cis window so that variant_index references are self-contained within the release.
10. As a data engineer, I want the completion region for each analysis to be trait_bp ± 1 Mb so that it matches the declared cis window used by the source dataset.
11. As a developer, I want the imputed mask stored as a uint8 zarr array aligned to the ragged association sequence so that query code can identify association provenance in O(1).
12. As a developer, I want the SE for imputed associations derived from the scalar-N model using LD panel EAF so that the Z+SE pair is consistent and the SE contract from ADR 0013 is met.
13. As a developer, I want the completion to inherit the LD eigenvector imputation kernel from pleiodb (`_ld_pca`, `_elastic_net_impute`, `_poly_rescale`) so that the algorithm is not re-implemented from scratch.
14. As a developer, I want the completed store to validate cleanly with `opengwasdb validate` so that build correctness is verifiable without re-running imputation.
15. As a developer, I want a quality gate (Pearson r ≥ 0.7) so that poorly imputed blocks leave reference variants as NaN rather than injecting low-quality associations.
16. As a developer, I want per-analysis completion bounded to LD blocks intersecting the cis window so that completion does not silently extend beyond the declared region.

## Implementation Decisions

### Enhancement model

Completion reads an existing observed-only ragged store (via `RaggedCSRReader`) and writes a new store directory with:
- The same `store_id`; a new `release_id` supplied at completion time
- `manifest.json` with `completion_state: "reference_completed"` and the declared LD panel as `ld_panel_id`
- An extended `data.zarr/ragged/` with four existing arrays plus one new array:
  - `offsets`, `variant_index`, `z`, `se` — extended with new entries for reference-panel variants
  - `imputed` — new uint8 array (0 = observed, 1 = imputed), aligned to the ragged sequence
- An extended `variants.tsv.gz` / `variant_alid_bytes.npy` covering all new reference-panel variants
- An extended `index.sqlite` with a `completion_quality` table
- Existing `traits.tsv.gz` and tabix index are copied unchanged

### LD reference panel

The existing panel at `ld_dir/{ancestry}/{chr}/{block_name}.{tsv,unphased.vcor1.gz,ldeig.npz}` is used. Block layout is the flat production layout from pleiodb. Block TSVs provide CHR, SNP (ALID), OA, EA, EAF, BP.

### Completion region

For each analysis, the completion region is `[trait_bp − 1_000_000, trait_bp + 1_000_000]` on `trait_chr`. LD blocks that intersect this window are processed. Only blocks with ≥ 2 observed variants after matching against the store variant table are imputable; the remaining blocks still contribute their reference-panel variants as NaN rows.

### Reference variant set

All LD panel variants falling within any analysis's cis window are added to the completed store's variant table. Store-local variant indices are reassigned in the new release (no cross-release stability). Observed variants that are already in the store retain their observed status; LD panel variants not present in the observed-only store are new entries.

### Imputation kernel (Z)

Ported from `pleiodb/src/pleiodb/impute.py`:
- `_ld_pca`: eigendecompose the LD matrix, retain components explaining ≥ 90% cumulative variance (threshold configurable)
- `_elastic_net_impute`: ElasticNetCV (α = 0.5, cv = min(5, n_obs−1)) on observed z-scores using eigenvectors as features
- `_poly_rescale`: degree-3 polynomial rescaling to calibrate predictions to observed scale; Cook's-distance outlier removal
- Quality gate: Pearson r ≥ 0.7 (configurable) between imputed and observed z-scores at observed positions; blocks below threshold yield NaN rows

### SE imputation (scalar-N model)

For each analysis × LD block, SE for imputed variants is derived as:

```
se_scale = median(se_obs * sqrt(2 * EAF_obs * (1 − EAF_obs)))    # from observed variants in block
se_imputed[i] = se_scale / sqrt(2 * EAF_ref[i] * (1 − EAF_ref[i]))  # EAF from LD panel TSV
```

EAF for observed variants is matched from the LD panel TSV by ALID. If no observed-variant EAF is available, the block's SE imputation is skipped and SE is left as NaN alongside the NaN z-scores.

### Association encoding

Per ADR 0013: finite z + se with `imputed=0` = observed; finite z + se with `imputed=1` = imputed; NaN z + NaN se with `imputed=0` = missing. Inconsistent states (one NaN, imputed=1 with NaN) are rejected by validation.

### CSR assembly per analysis

For each analysis:
1. Read observed associations from the source CSR (sorted by variant_index).
2. Find all LD panel variants in the cis window; determine which are new (not in observed store).
3. Merge observed and new reference-panel variants, sorted by variant_index.
4. Run imputation per intersecting LD block; populate z_imp, se_imp, imputed_mask.
5. Append the merged association sequence to the new CSR.

### Reference Completion Quality

A `completion_quality` table in `index.sqlite`:

```sql
CREATE TABLE completion_quality (
    analysis_index INTEGER NOT NULL,
    block_id       TEXT NOT NULL,
    pearson_r      REAL,
    n_imputed      INTEGER NOT NULL,
    n_missing      INTEGER NOT NULL,
    PRIMARY KEY (analysis_index, block_id)
);
```

`pearson_r` is NULL when the block was not imputable (< 2 observed variants or quality gate failed).

### Manifest additions

`manifest.json` gains:
- `completion_state`: `"reference_completed"` (was `"observed_only"`)
- `ld_panel_id`: string identifying the LD panel used (e.g. `"eur-hg38-gpm"`)
- `reference_completion_method`: `"elastic_net_eigenvectors_v1"`
- `cis_window_bp`: 1000000

### Modules

- `opengwasdb/layouts/ragged/complete.py` — new: `complete_ragged_store(source, dest, ld_dir, ...)`
- `opengwasdb/layouts/ragged/impute.py` — new: ported elastic-net kernel + scalar-N SE; thin wrapper around pleiodb math functions
- `opengwasdb/cli/main.py` — new `complete-ragged` command
- `opengwasdb/validation/validate.py` — extend `_validate_ragged_store` to check `imputed` array alignment and `completion_quality` table when `completion_state == "reference_completed"`
- `opengwasdb/query/facade.py` — extend `RaggedStoreQuery` to honour `observed_only` flag and expose `association_status` in result rows
- `opengwasdb/model/enums.py` — add `CompletionState.REFERENCE_COMPLETED`

## Testing Decisions

- Good tests assert external behaviour through build, validate, and query APIs; they do not check internal intermediate state.
- Unit tests for the imputation kernel: known LD matrix + z-scores → verify imputed values recover known signal; verify SE scalar-N formula; verify quality gate rejects low-r blocks.
- Unit tests for CSR assembly: synthetic 3-analysis BESD + tiny LD panel fixture → verify merged variant order, imputed mask values, NaN rows for missing reference variants.
- Integration test: build observed-only store from BESD fixture, run `complete_ragged_store`, validate, query — assert imputed associations appear by default and are excluded in observed-only mode.
- Validation test: inject an inconsistent imputed mask (imputed=1 with NaN z); assert validator rejects it.
- Query test: `range_by_analysis` on a completed store returns both observed and imputed associations; `observed_only=True` excludes imputed.
- Prior art: `tests/test_ragged_build_besd.py` for the observed-only build pattern; `pleiodb/tests/test_impute.py` for the kernel unit-test pattern.

## Out of Scope

- Dense reference completion (separate PRD)
- Ragged overflow for Dense stores
- Multi-ancestry LD panels or per-analysis ancestry selection
- Significant trans region completion (only cis ±1 Mb)
- Re-imputation or partial update of an existing completed store
- VCF-based dense imputation mode (pleiodb "dense mode") — only sparse mode (LD eigenvectors on observed variants in the store)
- Storing per-association imputation INFO score or R²

## Further Notes

- The imputation kernel is ported from `pleiodb/src/pleiodb/impute.py` (`_ld_pca`, `_elastic_net_impute`, `_poly_rescale`, `_ratio_outlier_mask`). The SE logic (`_se_outliers`) is not ported — SE is handled by the scalar-N model instead.
- LD panel is at `/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR/`. Block layout: `{ancestry}/{chr}/{block_name}.{tsv,unphased.vcor1.gz,ldeig.npz}`.
- ADR 0011–0015 define the conceptual contract; this PRD defines the implementation.
- The first target dataset is `eqtlgen-cis.opengwasdb` at `/local-scratch/data/opengwas/opengwasdb/eqtlgen-cis.opengwasdb`.
