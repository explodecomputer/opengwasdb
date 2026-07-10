## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

Complete the query surface for a Hybrid store on top of the tracer
(`issues/055-hybrid-tracer-bullet.md`): `analysis`, `phewas`, `range_phewas`,
and `lookup` all return the union of the Dense Component and the Ragged Overflow
Component with correct `association_status`.

The hybrid facade **dispatches to and composes the existing dense `StoreQuery`
and `RaggedStoreQuery`** per component and concatenates — it does not reimplement
matrix or CSR reads (PRD "Query", "Component reuse"). Dispatch is by panel
membership: a variant is on-panel (served by the dense component) xor off-panel
(served by the overflow), so results are a plain concatenation with no dedup.
`phewas` for an off-panel variant inherits the ragged O(n) scan (PRD "Off-panel
PheWAS"); on-panel `phewas` uses the fast dense row read.

## Acceptance criteria

- [ ] `analysis(id)` returns the analysis's dense-column finite cells **plus** its
      overflow sequence, unified.
- [ ] `phewas(variant)` dispatches to the dense row (on-panel) or the overflow
      (off-panel) and returns all analyses with that variant.
- [ ] `range_phewas(chrom, start, end)` returns both components' associations in
      the range.
- [ ] `lookup` (from 055) continues to pass; all methods carry correct
      `association_status` and honour `observed_only`.
- [ ] Results equal a direct read of each underlying component (dense
      `StoreQuery` / `RaggedStoreQuery`) unioned — the hybrid facade adds dispatch,
      not new storage reads.

## Blocked by

- Blocked by `issues/055-hybrid-tracer-bullet.md`

## User stories addressed

- User story 5
- User story 6
- User story 7
- User story 14
