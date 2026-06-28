# PRD: OpenGWASDB v0.1 Dense Observed-Only Vertical Slice

## Problem Statement

OpenGWAS needs a standalone, storage-efficient, fast-querying summary-statistics store that can support local Python use and later service-backed federation. The current `besdq` repository contains useful prototypes and lessons, but its implementation carries legacy BESD concepts, sparse SQLite assumptions, and prototype scripts that do not cleanly express the new OpenGWASDB store contract.

The first OpenGWASDB vertical slice must prove the core store architecture end to end without taking on every future mode. It should build and query a **Dense Observed-Only** Store Release using the new specification:

- `manifest.json` for release identity and contract declaration;
- `index.sqlite` for metadata and variant/analysis lookup;
- `data.zarr/` for compressed `z` and `se` arrays;
- canonical variant identity and allele orientation;
- layout-independent query API;
- top-hit indexing;
- validation sufficient to catch invalid stores early.

This first slice must also preserve the broader design path: Ragged layout, Reference Completion, imputation masks, Dense+Overflow, and service deployment are out of implementation scope for v0.1 but already recorded in the spec and ADRs. The v0.1 implementation must not bake in assumptions that make those later designs difficult.

## Solution

Build the first production OpenGWASDB vertical slice around Dense Observed-Only stores.

The solution should:

1. Define a clear Store Release writer that creates the standard envelope.
2. Build a dense source-faithful variant-by-analysis matrix from a small tabular/GWAS summary-stat input.
3. Store `z` and `se` as the canonical statistics.
4. Derive beta and p-value at query time.
5. Store variant and analysis metadata in SQLite.
6. Store numerical arrays in Zarr with compression/chunking selected from the existing dense benchmarks.
7. Provide a layout-independent query API even though only Dense is implemented initially.
8. Provide CLI commands that operate on explicit input/output store paths.
9. Provide validation tooling that confirms the release conforms to the v0.1 spec.
10. Provide tests using tiny deterministic fixtures so the store contract can evolve safely.

The implementation should treat `besdq` as prior art, not as a codebase to copy wholesale. Specific modules/scripts in `besdq` are useful references for parsing, allele normalisation, dense Zarr writing, top-hit indexing, and tests, but OpenGWASDB should keep its own clean terminology and module boundaries.

## User Stories

1. As an OpenGWAS data engineer, I want to build a Dense Observed-Only Store Release from a small batch of summary-statistic files, so that I can validate the new storage contract end to end.
2. As an OpenGWAS data engineer, I want the build command to write `manifest.json`, `index.sqlite`, and `data.zarr/`, so that every output follows the standard store envelope.
3. As an OpenGWAS data engineer, I want build and query commands to accept explicit store paths, so that store naming and catalogue concerns remain outside v0.1.
4. As an OpenGWAS data engineer, I want the manifest to declare layout, coverage, completion state, reference assembly, and provenance, so that a downloaded store is self-describing.
5. As an OpenGWAS data engineer, I want Dense Observed-Only stores to be source-faithful, so that they do not require an LD reference panel or reference-completed axis.
6. As an OpenGWAS data engineer, I want variants normalised to canonical ALID orientation, so that associations are comparable and signed statistics are consistent.
7. As an OpenGWAS data engineer, I want multi-allelic or ambiguous alleles handled explicitly, so that invalid source records do not silently corrupt the store.
8. As an OpenGWAS data engineer, I want the variant table to include Store-local variant indices, so that Zarr arrays can be compact and fast.
9. As an OpenGWAS data engineer, I want variant indices treated as release-local, so that future releases can rebuild or reorder variants without cross-store identity breakage.
10. As an OpenGWAS data engineer, I want analysis metadata stored separately from trait metadata, so that multiple analyses of the same trait can be represented correctly.
11. As an OpenGWAS data engineer, I want `z` and `se` stored directly, so that beta can be derived without EAF or sample size reconstruction.
12. As an OpenGWAS data engineer, I want SE stored as non-negative, so that effect direction is carried only by signed Z.
13. As an OpenGWAS data engineer, I want beta derived as `z * se`, so that users can still request beta while the store preserves the canonical statistic pair.
14. As an OpenGWAS data engineer, I want p-values derived from Z, so that p-values do not need separate storage in v0.1.
15. As an OpenGWAS data engineer, I want Dense arrays to use canonical NaN for unavailable observed-only cells, so that missingness is represented consistently and compressibly.
16. As an OpenGWAS data engineer, I want Zarr compression and chunking to be configurable but defaulted sensibly, so that the first implementation is performant without hard-coding one benchmark forever.
17. As an OpenGWAS data engineer, I want the builder to emit enough metadata to validate dimensions and array names, so that broken partial stores are detected.
18. As an OpenGWAS data engineer, I want validation to check manifest, SQLite, Zarr dimensions, statistic constraints, and variant/analysis counts, so that invalid stores fail fast.
19. As an OpenGWAS data engineer, I want the query engine to open exactly the path provided, so that default release selection is left to a future catalogue.
20. As an OpenGWAS data engineer, I want a layout-independent query object, so that Dense-specific details do not leak into public APIs.
21. As an OpenGWAS data engineer, I want range queries, so that I can retrieve associations for genomic regions.
22. As an OpenGWAS data engineer, I want exact variant queries, so that I can retrieve associations for specific variants.
23. As an OpenGWAS data engineer, I want analysis queries, so that I can extract all available associations for a given analysis.
24. As an OpenGWAS data engineer, I want full-analysis extraction, so that downstream pipelines can export one GWAS efficiently.
25. As an OpenGWAS data engineer, I want PheWAS-style variant extraction across analyses, so that I can quickly query one variant against many traits.
26. As an OpenGWAS data engineer, I want top-hit queries, so that high-significance associations can be retrieved without scanning full arrays.
27. As an OpenGWAS data engineer, I want top-hit indexes stored inside the store, so that top-hit lookup is fast and portable with the release.
28. As an OpenGWAS data engineer, I want top-hit thresholds to match the accepted design, so that later Dense/Ragged merging can share the same query contract.
29. As a Python package user, I want to install OpenGWASDB with normal Python packaging tools, so that querying does not require external database servers.
30. As a Python package user, I want query results to include standard fields such as variant identity, analysis identity, Z, SE, beta, p-value, and effect scale, so that results are directly usable.
31. As a Python package user, I want query result fields to be consistent with future Ragged and Reference-Completed stores, so that client code does not need to be rewritten later.
32. As a Python package user, I want observed-only stores to omit association status fields unless needed, so that the v0.1 result contract remains simple.
33. As a Python package user, I want CLI `info` or validation output to explain what is inside a store, so that I can diagnose build outputs quickly.
34. As a future API developer, I want the local query engine to be reusable by an API service, so that the service can wrap the same core query semantics.
35. As a future API developer, I want the store to be standalone, so that data can be deployed in multiple locations and federated later.
36. As a future Reference Completion implementer, I want v0.1 to avoid assuming observed-only forever, so that imputed associations and masks can be added without rewriting the public API.
37. As a future Ragged implementer, I want v0.1 to keep layout-specific code behind adapters, so that Ragged can be added without changing caller-facing query methods.
38. As a future builder implementer, I want source parsing and normalisation separated from physical layout writing, so that Dense and Ragged writers can share ingestion logic.
39. As a future builder implementer, I want clear test fixtures and expected outputs, so that adding new input formats does not break the store contract.
40. As a maintainer, I want the v0.1 implementation to reference existing `besdq` prior art explicitly, so that useful parsing and dense benchmark knowledge is not lost.

## Implementation Decisions

- Implement the first vertical slice as Dense Observed-Only only.
- Keep Ragged layout, Reference Completion, Dense+Overflow, imputed masks, and API/catalogue service out of v0.1 implementation scope, while preserving their design constraints in the spec and ADRs.
- Use the standard store envelope: manifest, SQLite index, and Zarr data hierarchy.
- Use explicit store paths for build and query CLI commands.
- Use `z` and `se` as the canonical stored statistic pair.
- Derive beta as `z * se`.
- Derive p-value from Z.
- Treat SE as non-negative and Z as the signed effect-direction statistic.
- Use canonical NaN in `z` and `se` for unavailable Dense Observed-Only cells.
- Do not add association status or imputed mask in v0.1 Observed-Only stores.
- Implement a layout-independent query facade from the start, even though the only backing adapter is Dense.
- Implement a layout-independent build orchestration layer, with Dense writer as the first physical writer.
- Store metadata and lookup structures in SQLite.
- Store Dense numerical arrays in Zarr.
- Include a top-hit index for Dense stores using the accepted significance thresholds.
- Keep top-hit index semantics compatible with future Ragged top-hit indexes.
- Implement validation as a first-class module, not as ad hoc builder checks.
- Use a tiny fixture-driven test suite to prove store creation, validation, and query behaviour.
- Port only selected prior art from `besdq`; do not copy legacy BESD terminology or schemas wholesale.

## Proposed Modules

- Domain model: manifest, enums, analysis metadata, variant metadata, query result types.
- Store opening: path-based store loading, manifest validation, and component path resolution.
- SQLite index: schema creation, migrations/versioning, variant lookup, analysis lookup, range lookup.
- Variant normalisation: ALID construction, allele orientation, signed-statistic flipping, long-allele policy hooks.
- Source readers: minimal summary-stat/GWAS-SSF-compatible reader for v0.1 fixtures, with room for additional input adapters.
- Common build pipeline: source ingestion, normalisation, metadata assembly, validation, then dispatch to Dense writer.
- Dense writer: build variant axis, create Zarr `z`/`se`, fill columns, write metadata/indexes.
- Dense query adapter: exact variant, range, analysis, full-analysis, and top-hit retrieval from Dense arrays.
- Query facade: layout-independent public interface selected from manifest layout.
- Top-hit index builder: layout-specific dense index using shared thresholds/result contract.
- Validation: manifest/index/Zarr consistency, statistic constraints, missingness rules, and top-hit index consistency.
- CLI: build, validate, info, and query entry points operating on explicit paths.

## Existing `besdq` Prior Art to Reference

Use the following as implementation references, not as direct architecture templates:

- GWAS-SSF row model and allele orientation logic: `/Users/gh13047/repo/besdq/besdq/gwas_ssf_reader.py`.
- Faster candidate scanning and partial-split parsing ideas: `/Users/gh13047/repo/besdq/besdq/gwas_ssf_fast_reader.py`.
- Two-pass GWAS-SSF build lessons, trait metadata, retained-row accounting, and intermediate handling: `/Users/gh13047/repo/besdq/besdq/gwas_ssf_builder.py`.
- Current SQLite query patterns and beta/SE reconstruction behaviour to avoid or simplify: `/Users/gh13047/repo/besdq/besdq/sqlite_query.py`.
- Significance filtering logic and tests for cis/significant/suggestive concepts, useful later for Ragged Cis-and-Signals: `/Users/gh13047/repo/besdq/besdq/significance_filter.py`.
- Dense manifest extraction prototype: `/Users/gh13047/repo/besdq/scripts/dense_01_make_manifest.py`.
- Dense union variant table prototype and ALID row-index map idea: `/Users/gh13047/repo/besdq/scripts/dense_02_build_variant_table.py`.
- Dense Zarr array writer prototype, including NaN initialisation and chunk/compression lessons: `/Users/gh13047/repo/besdq/scripts/dense_03_build_stats_arrays.py`.
- Dense top-hit index prototype using flat hit arrays plus offsets: `/Users/gh13047/repo/besdq/scripts/dense_11_build_sig_index.py`.
- Dense Zarr vs TileDB and chunking benchmark material: `/Users/gh13047/repo/besdq/docs/zarr-vs-tiledb-benchmark.qmd`, `/Users/gh13047/repo/besdq/docs/chunk-benchmark.qmd`, and `/Users/gh13047/repo/besdq/docs/sig-index-benchmark.qmd`.
- Existing tests for parser, allele flipping, significance filtering, stage build, and query behaviour: `/Users/gh13047/repo/besdq/tests/test_gwas_ssf_import.py`, `/Users/gh13047/repo/besdq/tests/test_gwas_ssf_fast_reader.py`, `/Users/gh13047/repo/besdq/tests/test_queries.py`, `/Users/gh13047/repo/besdq/tests/test_stage1.py`, and `/Users/gh13047/repo/besdq/tests/test_stage2.py`.

Important differences from `besdq`:

- OpenGWASDB must not use `probe`, `epi`, `esi`, or BESD legacy terminology in public APIs.
- OpenGWASDB stores `z` and `se`, not `beta` and `z`.
- OpenGWASDB uses manifest + SQLite + Zarr as the core envelope.
- OpenGWASDB query API must be layout-independent from day one.
- OpenGWASDB Observed-Only Dense stores are source-faithful and do not require LD reference panels.

## Testing Decisions

Tests should exercise external behaviour and store-contract validity, not internal implementation details. A good test should build or open a tiny store, issue a query or validation command, and assert the externally visible result.

Test modules:

- Manifest tests: load valid manifest, reject missing/invalid required fields, reject unsupported enum values.
- Store opening tests: open exactly the supplied path; do not perform discovery or default-release selection.
- Variant normalisation tests: canonical ALID construction, allele flipping, signed Z flipping, EAF flipping when EAF is present.
- Dense build tests: build a tiny two-analysis store with overlapping and missing variants; verify manifest, SQLite metadata, Zarr arrays, and NaN missingness.
- Dense query tests: exact variant query, range query, full-analysis extraction, and PheWAS-style variant extraction.
- Statistic tests: beta derivation from `z * se`, p-value derivation from Z, non-negative SE validation.
- Top-hit index tests: build index from tiny Z matrix; verify thresholds, offsets, returned variants, and NaN exclusion.
- Validation tests: reject inconsistent dimensions, missing arrays, negative SE, non-canonical missingness, duplicate variants, and invalid effect scale.
- CLI tests: build/validate/info commands operate on explicit paths and return useful errors.

Prior art for tests:

- `besdq` parser and allele-flip tests demonstrate useful fixture style.
- `besdq` significance-filter tests are useful later for Ragged Cis-and-Signals but should not drive v0.1 Dense Observed-Only scope.
- OpenGWASDB should use pytest-style tests consistently rather than mixing unittest patterns unless there is a strong reason.

## Out of Scope

- Reference Completion implementation.
- LD reference panel ingestion.
- Imputation of Z or SE.
- Imputed mask storage.
- Ragged primary layout implementation.
- Ragged Reference-Completed region implementation.
- Dense Ragged Overflow implementation.
- Remote/object-store direct querying.
- Multi-store catalogue, default-release selection, and deployment naming conventions.
- Full OpenGWAS API service.
- Legacy BESD file compatibility.
- Migration of existing `besdq` SQLite stores.
- Performance optimisation beyond basic chunk/compression choices needed for a usable first slice.

## Further Notes

The v0.1 implementation should be built against the current OpenGWASDB specification and ADRs:

- `docs/spec/store-format.md`
- `docs/adr/0001-hybrid-sqlite-zarr-store-envelope.md`
- `docs/adr/0005-store-z-and-se-as-canonical-statistics.md`
- `docs/adr/0006-layout-independent-query-engine.md`
- `docs/adr/0007-layout-independent-build-pipeline.md`
- `docs/adr/0008-implement-dense-observed-only-first.md`

The core delivery criterion is an end-to-end vertical slice:

```text
tiny source inputs
  -> build Dense Observed-Only Store Release
  -> validate release
  -> query variants/ranges/analyses/top hits
  -> return correct z, se, beta, p-value, and metadata
```

This PRD intentionally keeps the first implementation small. The broader imputation and ragged designs are already captured in ADRs and the spec; v0.1 should preserve those extension points without implementing them.
