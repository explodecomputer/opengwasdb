## Parent PRD

`issues/prd-ragged-reference-completion.md`

## What to build

Wire up the CLI command and extend the validator for reference-completed ragged stores.

**CLI** (`opengwasdb/cli/main.py`):

```
opengwasdb complete-ragged <source_store> <dest_store> \
    --ld-panel <path> \
    [--ancestry EUR] \
    [--cis-window-bp 1000000] \
    [--min-cor 0.7] \
    [--release-id <string>] \
    [--overwrite]
```

Calls `complete_ragged_store(...)` and prints a summary (n_analyses completed, total imputed associations, wall time).

**Validator** (`opengwasdb/validation/validate.py`), extending `_validate_ragged_store`:

- When `completion_state == "reference_completed"`:
  - `data.zarr/ragged/imputed` must exist and have the same length as `z`
  - All `imputed` values must be 0 or 1
  - Where `imputed == 1`: z and se must both be finite
  - Where z is NaN: se must also be NaN and `imputed` must be 0
  - `completion_quality` table must exist in `index.sqlite` with columns `(analysis_index, block_id, pearson_r, n_imputed, n_missing)`
  - `manifest.json` must have `ld_panel_id`, `reference_completion_method`, `cis_window_bp`

## Acceptance criteria

- [ ] `opengwasdb complete-ragged --help` lists all options
- [ ] Running against the eqtlgen observed-only store produces a validated completed store
- [ ] `validate_store` on a completed store with a corrupted `imputed` array (imputed=1, NaN z) returns an error
- [ ] `validate_store` on a completed store with missing `completion_quality` table returns an error
- [ ] `validate_store` on an observed-only store is unchanged (no regression)
- [ ] Tests for the two new validation failure cases

## Blocked by

- `issues/040-ragged-csr-completion-writer.md`

## User stories addressed

- User story 6 (CLI command)
- User story 14 (validate cleanly)
