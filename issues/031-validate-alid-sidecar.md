## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Add sidecar presence and consistency checks to `validate_store()`:

- `variant_alid_bytes.npy` exists
- `variant_alid_rows.npy` exists
- Both arrays have the same length
- That length matches the variant count in SQLite

Stores built before issue 029 will fail this check with a clear message pointing to the
need for a rebuild.

## Acceptance criteria

- [ ] `validate_store()` reports an error if either sidecar file is absent
- [ ] `validate_store()` reports an error if sidecar length ≠ variant count
- [ ] `validate_store()` passes for a freshly built store
- [ ] Error message names the missing file and suggests a rebuild

## Blocked by

- `issues/029-mmap-alid-sidecar-write.md`

## User stories addressed

- User story 5 (stores built without the sidecar fail loudly)
