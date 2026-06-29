## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Upgrade v0.1 fixture-grade variant normalisation into production-grade variant
normalisation before ingesting complex indel-heavy datasets.

The first vertical slice validates simple A/C/G/T alleles and canonicalises ALID
orientation, but it deliberately does not yet implement full trimming,
left-alignment, long-allele hashing, or retention of full normalised alleles for
hashed ALIDs. Those behaviours are required by `docs/spec/store-format.md` and
ADR 0017 before production ingestion.

## Acceptance criteria

- [ ] Alleles are trimmed and left-aligned before ALID assignment.
- [ ] Long alleles can use deterministic compact hashed ALIDs.
- [ ] Full normalised alleles are retained once per variant when compact ALIDs
      are used.
- [ ] Multi-allelic source rows are split, rejected, or otherwise handled with
      explicit test-covered behaviour.
- [ ] Tests cover SNPs, insertions, deletions, left-alignment edge cases,
      multi-allelic rows, and long-allele hashing.

## Blocked by

None.

## User stories addressed

- User story 6
- User story 7
- User story 39
