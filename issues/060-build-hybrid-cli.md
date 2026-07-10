## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

A `build-hybrid` CLI command that generates a Hybrid store from a source manifest
and a reference panel, wiring the integrated build (issue 055) to the CLI and
reusing the dense build's **parallel, band-streamed** machinery so it runs at
genome scale within bounded memory (PRD "Generation", "Component reuse";
ADR 0024).

Options mirror the dense builder plus the panel: manifest path, output path,
`--reference-panel` (the dense axis / completion panel), `--store-id`,
`--release-id`, `--n-workers`, `--chunk-variants`/`--chunk-analyses`, and the
liftover controls. Sources may be GWAS-VCF or delimited text (reuse the existing
VCF and normalised-association readers), since heterogeneous collections arrive in
mixed formats.

## Acceptance criteria

- [ ] `opengwasdb build-hybrid <manifest> <out> --reference-panel <panel> ...`
      builds a Hybrid store whose dense axis is the panel and whose overflow holds
      the off-panel tail; the result queries (055/056) and validates (059).
- [ ] The build uses the dense builder's fork-pool + spill + band-write path
      (`--n-workers` is a pure runtime knob; results independent of worker count).
- [ ] Chunk-shape options are honoured on the Dense Component.
- [ ] Both VCF and delimited-text sources are accepted via the existing readers.

## Blocked by

- Blocked by `issues/055-hybrid-tracer-bullet.md`

## User stories addressed

- User story 10
- User story 11
- User story 14
