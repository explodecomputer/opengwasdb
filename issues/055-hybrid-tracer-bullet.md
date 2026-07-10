## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

The thin end-to-end path for the Hybrid Layout, on a tiny synthetic store.
Introduce `PrimaryStorageLayout = hybrid` and build one integrated store whose
Dense Component axis is a supplied reference panel and whose Ragged Overflow
Component holds each study's off-panel observed variants, over a single shared
Store Variant Table — then read a cell back through a unifying query.

- **Model/manifest**: add the `hybrid` layout; the manifest declares two
  components over one variant table (PRD "Representation").
- **Integrated build**: extend the existing two-pass builder to route each
  study's variants by panel membership — on-panel → dense fill, off-panel →
  ragged CSR — in a single read. Reuse the dense builder's spill/band-write and
  `RaggedCSRWriter`; the only new logic is the routing and the shared variant
  table = panel ∪ off-panel-observed (PRD "Generation", "Component reuse").
- **Unifying query**: a `lookup` that dispatches by panel membership (on-panel →
  dense, off-panel → ragged) and concatenates, with `association_status`.

Serial and small is fine here; parallel/scale and the CLI come later.

## Acceptance criteria

- [ ] A hybrid store builds from a small manifest + a small reference panel, with
      a `data.zarr` dense matrix (panel axis) and a `data.zarr/ragged` overflow,
      over one shared variant table; manifest declares `primary_layout = hybrid`.
- [ ] On-panel source associations land in the dense matrix; off-panel ones land
      in the overflow; the partition is disjoint (no variant in both).
- [ ] The build reads each study once (routing, not two passes over the data).
- [ ] `lookup` returns the correct z/se/status for a requested variant × analysis
      whether the cell lives in the dense matrix or the overflow, matching a
      direct read of the underlying component.
- [ ] The build composes the dense two-pass machinery and `RaggedCSRWriter`
      rather than reimplementing matrix/CSR writes.

## Blocked by

None - can start immediately.

## User stories addressed

- User story 1
- User story 2
- User story 3
- User story 4
- User story 5
- User story 12
- User story 14
