## Parent PRD

`issues/prd-ragged-observed-only.md`

## What to build

A builder that reads filtered GWAS-SSF files (one `.tsv.gz` per analysis, as
produced by the opengwasdb-stores download+filter step) plus an analyses
manifest, and produces a valid Ragged Observed-Only Store — the GWAS-SSF
sibling of `036-ragged-build-besd.md`'s BESD builder, for eQTL/pQTL-style
Source Collections distributed as filtered summary-stat files rather than
BESD triples.

`opengwasdb/layouts/ragged/build_ssf.py` already contains a working draft of
`build_ragged_from_ssf()` (adapter mirroring `build_ragged_from_besd`: reads
the manifest, groups filtered rows per analysis, canonicalises variants,
writes the variant/traits axes + `index.sqlite` + ragged CSR + manifest). It
has never been tested, is not wired into the CLI, and is not imported
anywhere else in the codebase. This ticket is about bringing that draft up to
the same bar as its BESD sibling, not writing the builder from scratch —
review the draft's per-row filtering (finite/positive `se`, finite `beta`,
orientation via `orient_to_canonical`) and its `index.sqlite` schema (kept
compatible with the completion pipeline and query paths, per its own comment)
for correctness before trusting it as-is.

Add a CLI command `build-ragged-ssf` to `opengwasdb/cli/main.py`, following
the existing `build-ragged-besd` command's shape:

```
opengwasdb build-ragged-ssf <manifest_path> <filtered_dir> <output_path> --store-id <id> --release-id <id>
```

## Acceptance criteria

- [ ] `build_ragged_from_ssf(manifest_path, filtered_dir, output_path, store_id, release_id)` produces a complete, valid store directory
- [ ] `manifest.json` contains `primary_layout: "ragged"` and `completion_state: "observed_only"`
- [ ] `variants.tsv.gz` contains the canonical union of variants across all filtered files, tabix-indexed
- [ ] `traits.tsv.gz` contains one row per manifest analysis, tabix-indexed by `trait_chr`/`trait_bp`
- [ ] `data.zarr/ragged/` contains valid CSR arrays with correct offsets
- [ ] `opengwasdb validate <output_path>` exits 0
- [ ] `build-ragged-ssf` CLI command is reachable via `opengwasdb --help`
- [ ] Integration test: build from a synthetic multi-analysis filtered-GWAS-SSF fixture (including a row with non-finite/non-positive `se` and a row needing allele-orientation flip), assert variant count, analysis count, and round-trip z/se values
- [ ] `analysis_index` validation (must be dense `0..n-1`) is covered by a test with a malformed manifest

## Blocked by

None — the storage/axis/top-hit primitives it composes (`RaggedCSRWriter`, `write_variant_axis`, `write_traits_axis`, `build_ragged_top_hit_indexes`) are already shipped and used by `build_ragged_from_besd`.
