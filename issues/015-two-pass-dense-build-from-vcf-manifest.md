## Parent PRD

`issues/prd-vcf-dense-build-with-liftover.md`

## What to build

The two-pass dense build orchestrator in `opengwasdb/layouts/dense/build_vcf.py`, exposing:

```
build_dense_from_vcf_manifest(
    manifest_path,
    output_path,
    *,
    chain_file,
    store_id,
    release_id,
    liftover_failure_threshold=0.01,
)
```

**Pass 1** — streams all VCFs via `stream_vcf_variants` from `013`, collects the union of `(chrom, pos, ref, alt)` tuples, calls `build_liftover_lookup` from `014` once on the deduplicated set. Sorts the resulting hg38 ALIDs by `(chromosome_sort_key, pos, A1, A2)`, writes the SQLite variant table (same schema as `build_dense_observed_store`), and pre-allocates zarr z and se arrays of shape `(n_variants, n_analyses)` filled with NaN.

**Pass 2** — iterates the manifest row-by-row. For each trait, reads the VCF via `stream_vcf_associations`, resolves each `(chrom, pos, ref, alt)` through the hg19→hg38 lookup dict built in pass 1 (no second liftover calls), and writes z and se values into the correct zarr column. Handles allele flip: if the association's A1 (after ALID normalisation) is not the hg38 ALID's A1, z is negated.

After pass 2, calls `build_top_hit_indexes`, writes `manifest.json` with `reference_assembly: GRCh38`, `completion_state: observed_only`, and provenance recording the builder name and chain file basename.

Peak memory is O(n_variants): the liftover lookup dict plus the zarr arrays. The full association list is never materialised.

The output format must be byte-for-byte compatible with the existing store envelope — same zarr layout, same SQLite schema, same manifest structure — so that `validate_store`, `query_store`, and all existing query patterns work without modification.

## Acceptance criteria

- [ ] `build_dense_from_vcf_manifest` completes without error on a synthetic manifest pointing to two or three tiny GWAS-VCF fixtures (five to ten chr1 variants each).
- [ ] `validate_store` passes on the output store.
- [ ] `query_store` supports range, phewas, top-hits, and bulk analysis queries on the output store with correct z/se values.
- [ ] A variant where the hg19→hg38 lookup produces an ALID with flipped allele orientation has its z negated in the correct zarr cell.
- [ ] The output `manifest.json` contains `reference_assembly: GRCh38` and `completion_state: observed_only`.
- [ ] `stored_effect_scale` in the SQLite analyses table matches the VCF header `StudyType` for each trait.
- [ ] Peak memory during the build does not require holding all associations simultaneously (verified by inspection of the implementation, not a runtime assertion).
- [ ] A manifest with a VCF whose liftover failure rate exceeds `liftover_failure_threshold` causes the build to raise before writing any output.
- [ ] End-to-end test uses only synthetic GWAS-VCF fixtures in `tmp_path`; no dependency on besdq data files.

## Blocked by

- `issues/013-gwas-vcf-source-reader.md`
- `issues/014-inline-hg19-hg38-liftover.md`

## User stories addressed

- User story 1 (build from GWAS-VCF without modifying source files)
- User story 7 (build time and memory comparable to besdq prototype)
- User story 8 (output passes validate_store)
- User story 9 (query_store API unchanged)
- User story 10 (range, phewas, top-hits, bulk all work)
- User story 11 (allele orientation normalised end-to-end)
- User story 12 (manifest format compatible with besdq manifest TSV)
- User story 17 (output format identical to build_dense_observed_store)
- User story 19 (unit tests use synthetic VCF fixtures)
- User story 20 (end-to-end test builds and queries a tiny store)
