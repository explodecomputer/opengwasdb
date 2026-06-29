## Parent PRD

`issues/prd.md`

## What to build

Build a tiny dense observed-only store that writes the variant axis as bgzipped, tabix-indexed TSV plus row-offset sidecar, while keeping dense statistic arrays in Zarr and analysis metadata in compact SQLite. This should prove the new store envelope can be produced end to end for the existing tiny fixture.

## Acceptance criteria

- [ ] A dense build writes `variants.tsv.gz`, `variants.tsv.gz.tbi`, and `variant_offsets.npy`.
- [ ] The SQLite database no longer contains a high-cardinality canonical variants table for the new dense variant-axis format.
- [ ] Variant rows are sorted by canonical chromosome, position, effect allele, other allele, and variant index.
- [ ] `variant_index` values are assigned from sorted row order, starting at zero.
- [ ] The dense `z` and `se` arrays use the same row indices as the written variant table.
- [ ] The writer records one row offset per variant row.
- [ ] The writer does not insert every canonical ALID as an alias.
- [ ] A partial failed build does not leave a valid-looking release directory.

## Blocked by

- Blocked by `issues/013-variant-store-contract-and-dependency.md`

## User stories addressed

- User story 1
- User story 2
- User story 3
- User story 4
- User story 5
- User story 6
- User story 7
- User story 10
