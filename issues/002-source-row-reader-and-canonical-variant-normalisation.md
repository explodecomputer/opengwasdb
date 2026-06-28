## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Implement the first source-row model and canonical variant normalisation path. The slice should read tiny fixture summary-stat rows, construct canonical ALIDs, orient signed statistics to canonical A1, and expose normalised records for builders. Reference the `besdq` GWAS-SSF reader and fast reader as prior art, but keep OpenGWASDB terminology and interfaces.

## Acceptance criteria

- [ ] A tiny fixture source file can be streamed into normalised association records.
- [ ] ALID is constructed as `chr:pos:A1:A2` with alphabetically first A1.
- [ ] Signed statistics are negated when source effect allele orientation is reversed.
- [ ] SE remains non-negative after normalisation.
- [ ] Invalid or ambiguous rows are rejected or skipped with test-covered behaviour.
- [ ] Tests cover no-flip, flip, missing rsid, invalid rows, and canonical ALID output.

## Blocked by

None - can start immediately.

## User stories addressed

- User story 6
- User story 7
- User story 11
- User story 12
- User story 38
- User story 39
- User story 40

