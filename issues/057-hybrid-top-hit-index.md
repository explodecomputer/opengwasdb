## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

A `top_hits` query for a Hybrid store that ranks across **both** components. The
Dense Component and the Ragged Overflow Component each already have a top-hit
index (dense and ragged layouts respectively); build/merge them so the hybrid
`top_hits` returns the combined ranked set at a given threshold.

Reuse the existing dense and ragged top-hit index builders and readers
(PRD "Component reuse"); the new work is emitting the overflow's index during the
hybrid build and merging the two ranked result sets at query time.

## Acceptance criteria

- [ ] The hybrid build writes a top-hit index for the Dense Component and one for
      the Ragged Overflow Component, using the existing builders.
- [ ] `top_hits(threshold)` returns the union of both components' hits, ranked by
      significance, with correct `association_status` (a top hit may be an
      off-panel overflow association).
- [ ] Per-tier counts and min|z| equal the two components' indexes combined.
- [ ] No matrix/CSR rescan is introduced beyond the existing index machinery.

## Blocked by

- Blocked by `issues/055-hybrid-tracer-bullet.md`
- Blocked by `issues/056-hybrid-unified-query-surface.md`

## User stories addressed

- User story 6
- User story 14
