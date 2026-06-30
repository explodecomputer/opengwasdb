## Parent PRD

`issues/prd.md`

## What to build

Route canonical variant lookup, PheWAS lookup, and genomic range lookup through the tabix-backed variant axis while preserving the existing public query facade and result contract.

## Acceptance criteria

- [ ] `variant()` resolves a canonical ALID by parsing coordinates and alleles, fetching the single position through tabix, and filtering alleles.
- [ ] `phewas()` continues to return the same results as `variant()` for a canonical ALID.
- [ ] `range()` streams variant rows through tabix and returns finite dense associations in genomic order.
- [ ] Missing dense cells remain excluded from materialised result rows.
- [ ] Query results still include variant metadata, analysis metadata, `z`, `se`, derived beta, and derived p-value.
- [ ] The query facade does not require callers to know the physical variant-axis backend.
- [ ] Tests no longer assert that canonical variant lookup depends on SQLite.

## Blocked by

- Blocked by `issues/015-open-and-validate-tabix-variant-axis.md`

## User stories addressed

- User story 15
- User story 17
- User story 20
- User story 22
- User story 23
- User story 25
- User story 29
- User story 30
