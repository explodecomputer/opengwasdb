## Parent PRD

`issues/prd.md`

## What to build

Support rsid and other compact alias lookups without duplicating every canonical ALID in SQLite. Alias lookup should use the documented ambiguity policy from the variant-store contract and then materialise variants through the new variant axis.

## Acceptance criteria

- [ ] The build records rsid aliases separately from canonical ALID identity.
- [ ] Alias lookup resolves aliases to existing variant indices or canonical identities only.
- [ ] Ambiguous alias behaviour matches the documented policy.
- [ ] `variant()` supports rsid lookup for aliases present in the store.
- [ ] Absent aliases return an empty query result rather than raising.
- [ ] Validation rejects alias metadata that points to missing variants.

## Blocked by

- Blocked by `issues/016-canonical-alid-and-range-queries.md`

## User stories addressed

- User story 10
- User story 21
- User story 24
- User story 35
