## Parent PRD

`issues/prd.md`

## What to build

Decide and encode the initial contract for the tabix-backed dense variant axis. This slice should make the format, dependency choice, manifest metadata, and alias ambiguity policy explicit enough that the writer, reader, validator, and benchmarks can be implemented independently against the same target.

## Acceptance criteria

- [ ] The dense store contract names the active variant-axis files: `variants.tsv.gz`, `variants.tsv.gz.tbi`, and `variant_offsets.npy`.
- [ ] The manifest records enough variant-axis metadata for readers to reject unsupported dense variant-axis formats clearly.
- [ ] The runtime dependency for BGZF/tabix operations is added through normal project dependency management.
- [ ] Type-checking configuration is updated if the tabix dependency lacks complete typing metadata.
- [ ] Canonical ALID, rsid alias, and ambiguous alias behaviour are documented before query implementation starts.
- [ ] Existing dense statistic compression, dtype, and chunk defaults are unchanged.

## Blocked by

None - can start immediately

## User stories addressed

- User story 8
- User story 14
- User story 16
- User story 21
- User story 36
- User story 38
