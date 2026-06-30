# PRD: Dense Store Build from GWAS-VCF with Inline Liftover

## Problem Statement

The existing dense store builder only accepts pre-normalised TSV/CSV source files. Large GWAS collections (e.g. IEU OpenGWAS UKB ~4,000 traits) are distributed as GWAS-VCF files in GRCh37 (hg19). A researcher wanting to build a Dense Observed-Only Store from these files today must either convert each VCF to a custom TSV format (touching source files) or write ad-hoc ingestion scripts. Neither option scales: the TSV conversion duplicates storage, and the in-memory build pipeline becomes infeasible as trait count or variant count grows because it holds every association as a Python object before writing a single byte to zarr.

There is also no mechanism for lifting variant coordinates from hg19 to hg38 during ingestion. Per ADR-0003, each Store Release declares exactly one Reference Assembly, and hg38 is the target assembly for new releases.

## Solution

Add a two-pass streaming build path that ingests GWAS-VCF files directly, lifts hg19 variant coordinates to hg38 inline using pyliftover, and writes zarr arrays one trait-column at a time. Source files are never modified or copied. The output format is identical to `build_dense_observed_store` so all existing queries, validation, and top-hit indexes work without changes.

The build is driven by a manifest TSV (same format as the besdq prototype: `trait_id`, `file_path`, `trait_name`, `n`) and a single Python function. Per-trait metadata (`stored_effect_scale`) is inferred from the `StudyType` field in the VCF `##SAMPLE` header line rather than carried as a manifest column.

## User Stories

1. As a researcher, I want to build a Dense Observed-Only Store from a directory of GWAS-VCF files without converting or copying them, so that I avoid duplicating tens of gigabytes of source data.
2. As a researcher, I want variant coordinates lifted from GRCh37 to GRCh38 automatically during the build, so that the resulting store conforms to the hg38 reference assembly without manual preprocessing.
3. As a researcher, I want the build to fail loudly if liftover failures exceed a configurable threshold, so that a wrong chain file or assembly mismatch produces an error rather than a silently truncated store.
4. As a researcher, I want liftover failures below the threshold to be dropped and logged, so that a small number of unmappable telomeric or patch variants do not abort an otherwise successful build.
5. As a researcher, I want `stored_effect_scale` inferred automatically from the VCF header `StudyType` field, so that I do not need to annotate case-control vs continuous traits manually in the manifest.
6. As a researcher, I want the build to skip multi-allelic variants, so that only biallelic SNPs and indels are included in the store.
7. As a researcher, I want build time and peak memory to be comparable to the besdq prototype for the same input data, so that the new pipeline is practically usable at scale.
8. As a researcher, I want the resulting store to pass existing `validate_store` checks, so that I can be confident the output is a well-formed Store Release.
9. As a researcher, I want to query the VCF-built store with the same `query_store` API as any other store, so that downstream analysis code does not need to know how the store was built.
10. As a researcher, I want range queries, phewas queries, top-hit queries, and bulk analysis extraction to work on the VCF-built store, so that all existing query patterns are available.
11. As a researcher, I want allele orientation to be normalised to the canonical ALID convention (A1 = alphabetically first) during VCF ingestion, so that z-scores are consistent and cross-store variant matching works without allele harmonisation.
12. As a researcher, I want the manifest format to remain compatible with the existing besdq manifest TSV, so that I can reuse existing manifest files without modification.
13. As a researcher, I want the build to handle VCF files with both bare (`1`) and `chr`-prefixed (`chr1`) chromosome names, so that ingestion is robust to CHROM convention differences.
14. As a researcher, I want the `FORMAT/EZ` field used in preference to `FORMAT/ES / FORMAT/SE` when computing z-scores, so that the most accurate z-score representation is used when available.
15. As a developer, I want the VCF reader to be a layout-independent module separate from the dense orchestrator, so that it can be reused by future Ragged or Reference-Completed build paths.
16. As a developer, I want the two-pass orchestrator to reuse the hg19→hg38 position mapping built in pass 1 during pass 2, so that pyliftover is called once per unique variant rather than once per variant per trait.
17. As a developer, I want the VCF build function to produce identical zarr layout, SQLite schema, and manifest format as `build_dense_observed_store`, so that the query engine, validation, and top-hit index require no changes.
18. As a developer, I want the liftover failure threshold to be configurable with a sensible default (1%), so that different datasets with different expected failure rates can be accommodated.
19. As a developer, I want unit tests for the VCF reader that use synthetic in-memory GWAS-VCF fixtures, so that tests run without requiring real VCF files or network access.
20. As a developer, I want an end-to-end test that builds a tiny store from synthetic VCF fixtures and queries it with the existing query API, so that the full pipeline is verified without depending on the besdq data directory.

## Implementation Decisions

### New modules

**Layout-independent VCF reader** (`opengwasdb/build/vcf_source.py`):
- Uses `cyvcf2` to iterate records from a GWAS-VCF file.
- Exposes two functions:
  - One that reads only CHROM/POS/REF/ALT from a VCF, used in pass 1 to collect variants.
  - One that reads CHROM/POS/REF/ALT/ES/SE (with EZ fallback) from a VCF, used in pass 2 to fill stats.
- Parses `StudyType` from the `##SAMPLE` VCF meta-line and maps `CaseControl → LOG_ODDS`, `Continuous → SD_UNITS`. Raises loudly if the field is absent.
- Handles bare and `chr`-prefixed CHROM values transparently.
- Skips multi-allelic records (ALT containing `,`).
- Does not perform liftover itself; callers pass in a resolved position lookup or liftover function.

**Two-pass dense build orchestrator** (`opengwasdb/layouts/dense/build_vcf.py`):
- Exposes `build_dense_from_vcf_manifest(manifest_path, output_path, *, chain_file, store_id, release_id, liftover_failure_threshold=0.01)`.
- **Pass 1**: streams all VCFs via the VCF reader, collects union of (hg19 CHROM, POS, REF, ALT) tuples. For each unique variant, calls pyliftover once to produce a hg38 ALID. Builds a `dict[hg19_key → hg38_alid]` lookup in memory. Counts failures; raises if rate exceeds threshold. Sorts variants by (chromosome, position, A1, A2) and writes the SQLite variant table and pre-allocates zarr z and se arrays filled with NaN.
- **Pass 2**: iterates the manifest trait-by-trait. For each trait, reads the VCF using the VCF reader, resolves each record's hg19 position through the pre-built lookup (no second pyliftover calls), looks up the zarr column index, and fills z and se values. Handles allele orientation flip (negate z when effect allele is not the canonical A1).
- After pass 2, calls `build_top_hit_indexes` and writes the manifest JSON, matching the output format of `build_dense_observed_store` exactly.

### Liftover
- Uses `pyliftover.LiftOver` with a user-supplied chain file path (e.g. `hg19ToHg38.over.chain.gz`). The `LiftOver` object is created once per build.
- VCF positions are 1-based; pyliftover expects 0-based; the reader subtracts 1 before calling `convert_coordinate` and adds 1 to the result.
- Chromosome names are normalised to `chr`-prefixed form before liftover and stripped for the output ALID.

### Manifest format
- Input: TSV with columns `trait_id`, `file_path`, `trait_name`, `n` (same as besdq).
- `trait_id` maps to `analysis_id` in the Store.
- `trait_name` maps to `analysis_label`.
- `n` is stored as provenance in the manifest JSON but is not used for z/SE derivation (SE is read directly from the VCF FORMAT/SE field).

### StoredEffectScale inference
- `CaseControl` in VCF header → `StoredEffectScale.LOG_ODDS`.
- `Continuous` → `StoredEffectScale.SD_UNITS`.
- Any other or missing value → build raises `ValueError` with the trait ID.

### Output format
- Zarr layout: identical to `build_dense_observed_store` (z and se float16 arrays, Blosc zstd bitshuffle compression, same chunk shape defaults).
- SQLite schema: identical (same `variants`, `variant_aliases`, `analyses`, `metadata` tables).
- Manifest JSON: `reference_assembly` set to `GRCh38`, `completion_state` set to `observed_only`, provenance records `builder: opengwasdb.v0.1_dense_vcf`, chain file basename, and liftover failure count.

### Performance
- Peak memory is O(n_variants) — the hg19→hg38 lookup dict and the zarr arrays only. No full association list is ever held in memory.
- The cyvcf2 vs bcftools query performance trade-off is acknowledged; benchmarking against the besdq prototype will determine whether a bcftools-backed path is needed.

## Testing Decisions

Good tests for this feature verify observable outputs (zarr array values, SQLite row counts, manifest fields, query results, error messages) without asserting on internal call counts, class structure, or intermediate state.

### VCF reader unit tests
- Synthetic GWAS-VCF fixtures written to `tmp_path` using `cyvcf2`-compatible plain-text VCF format or bgzipped VCF — no real data files required.
- Tests verify: correct z derivation from ES/SE, EZ field preferred over ES/SE when present, z negation when effect allele is not canonical A1, StudyType parsing to stored_effect_scale, multi-allelic records skipped, missing FORMAT fields produce NaN, both bare and chr-prefixed CHROM handled.
- Prior art: `test_source_normalisation.py` which uses `tmp_path` TSV fixtures to test the TSV reader at the same level of abstraction.

### Two-pass build and liftover tests
- End-to-end test: build a store from synthetic GWAS-VCF fixtures (two or three traits, five to ten variants on chr1), then query with `query_store` and assert on specific z/se values, analysis metadata, and `validate_store` passing.
- Test liftover threshold enforcement: a synthetic set of variants where >1% fail liftover should raise.
- Test correct allele flipping through the full pipeline: a variant where the VCF has REF > ALT alphabetically should store the negated z in the correct ALID cell.
- Test that `reference_assembly: GRCh38` is present in the output manifest.
- Prior art: `test_dense_vertical_slice.py` which builds a store via `conftest.py` fixtures and then exercises the full query API against it.

### What is not tested
- Real VCF files from the besdq data directory — those belong in a benchmark script, not in the test suite.
- Internal pass boundary state (intermediate dicts, partial zarr writes).
- pyliftover internals — the liftover is tested by verifying that output ALIDs contain hg38 coordinates, not by mocking the library.

## Out of Scope

- Full-genome (all chromosomes) ingestion — the initial implementation targets chr1 as a like-for-like benchmark against the besdq prototype.
- Reference-Completed builds from VCF — only Observed-Only completion state is in scope.
- Ragged layout VCF builds — only Dense layout is in scope.
- Streaming ingestion from remote URLs — all VCF paths must be local.
- Support for VCF files without a CSI or TBI index — indexed files are assumed.
- A CLI subcommand exposing `build_dense_from_vcf_manifest` — Python API only for now.
- Parallel trait ingestion (multi-threaded pass 2) — single-threaded initially; parallelism can be added later if benchmarks show it is needed.
- Validation that the liftover chain file matches the declared source assembly — the VCF header is trusted as the source of truth.

## Further Notes

- The benchmark target is `ukb-chr1_zarr_benchmark.json` in the besdq repo, which recorded regional query median 2.7 ms, phewas median 0.9 ms, and bulk median 490 ms on a 100-trait chr1 store.
- The chain file (`hg19ToHg38.over.chain.gz`) must be obtained separately; the standard UCSC file is used by pyliftover automatically when given the build name pair.
- pleiodb (`~/repo/pleiodb`) contains a working reference implementation of pyliftover-based inline liftover (`pleiodb/src/pleiodb/liftover.py`) and a cyvcf2-based GWAS-VCF reader (`pleiodb/src/pleiodb/vcf.py`) that can be consulted during implementation.
- ADR-0019 records the architectural rationale for the two-pass approach, inline liftover, cyvcf2 choice, and liftover failure threshold policy.
