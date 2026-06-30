## Parent PRD

`issues/prd.md`

## What to build

Teach store opening and validation to understand the tabix-backed dense variant axis. A valid tiny store should be opened and validated without relying on a SQLite variants table for variant counts or genomic lookup.

## Acceptance criteria

- [ ] Store opening exposes paths for the variant table, tabix index, and row-offset sidecar.
- [ ] Validation fails clearly when `variants.tsv.gz` is missing.
- [ ] Validation fails clearly when `variants.tsv.gz.tbi` is missing.
- [ ] Validation fails clearly when `variant_offsets.npy` is missing or has the wrong length.
- [ ] Validation confirms the dense array row count matches the variant table row count.
- [ ] Validation detects duplicate canonical variants in the variant table.
- [ ] Validation confirms representative tabix range fetches agree with sequential variant-table reads.
- [ ] Validation confirms alias metadata points only to existing variants.

## Blocked by

- Blocked by `issues/014-write-tiny-tabix-variant-axis.md`

## User stories addressed

- User story 14
- User story 16
- User story 31
- User story 32
- User story 33
- User story 34
- User story 35
