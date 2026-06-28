## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Implement the layout-independent query facade and the first Dense adapter for exact variant and genomic range queries. Callers should query through the public interface without depending on Dense-specific array details.

## Acceptance criteria

- [ ] A Store Release opened from an explicit path can create a layout-independent query object.
- [ ] Exact variant queries return expected associations for a canonical variant identity or alias supported by the index.
- [ ] Range queries return expected associations within chromosome and coordinate bounds.
- [ ] Result rows include variant identity, analysis identity, Z, SE, stored effect scale, and enough metadata to interpret the result.
- [ ] Missing Dense cells are excluded or represented according to the query contract and documented behaviour.
- [ ] Tests build a tiny store and verify exact/range query outputs against expected fixture values.

## Blocked by

- Blocked by `issues/003-build-tiny-dense-observed-only-store.md`

## User stories addressed

- User story 19
- User story 20
- User story 21
- User story 22
- User story 30
- User story 31
- User story 32
- User story 34
- User story 37

