# Alternative Variant Store for Dense Observed-Only Releases

## Problem Statement

The current Dense Observed-Only vertical slice stores the high-cardinality variant axis in SQLite. That made the first build from three GWAS-SSF files slower than expected and produced a store footprint close to the compressed source inputs rather than the expected multi-fold reduction. The numerical Zarr arrays are already using the intended dense compression strategy, but the SQLite variant table and its lookup indexes dominate storage size at realistic variant counts.

The team needs a replacement variant-store structure that keeps the dense array model intact while moving genomic range and exact variant lookup to a more compact, sequentially writable, genomic-indexed representation inspired by the besdq prototype. The replacement must still support the public query behaviours already proven in the vertical slice: variant/PheWAS lookup, genomic range lookup, analysis extraction, sparse lookup across variants and analyses, and top-hit materialisation.

## Solution

Introduce a dedicated variant-axis store for dense releases. The dense association arrays remain in the compressed array hierarchy, and SQLite remains available for low-cardinality relational metadata such as analyses, release metadata, and small alias maps. The Store Variant Table moves out of SQLite into a bgzipped, tabix-indexed TSV.

The new dense release envelope will contain:

- `manifest.json`: release contract, including the declared variant table format and dependency-free metadata needed to open the store.
- `index.sqlite`: compact metadata database for analyses, key/value metadata, and optional low-cardinality alias records.
- `data.zarr/`: compressed dense statistic arrays and dense numerical sidecar indexes such as top hits.
- `variants.tsv.gz`: canonical variant table, sorted by genomic coordinate and allele, bgzip-compressed.
- `variants.tsv.gz.tbi`: tabix index over chromosome and position.
- `variant_offsets.npy`: compact row-index to BGZF virtual-offset sidecar for direct materialisation of variants by dense row index.

The variant table will use store-local variant indices as the join key to dense arrays. A proposed column contract is:

- chromosome
- position
- variant_index
- effect_allele
- other_allele
- alid
- rsid

The table is sorted by chromosome, position, effect allele, other allele, and variant index. `variant_index` is assigned in exactly this sorted order, starting at zero. `alid` remains the canonical within-store variant identity. `rsid` remains an alias and must not become primary identity.

Writes become a streaming-friendly two-phase process:

1. Collect or stream-normalise association records into a canonical variant set and analysis set.
2. Sort the variant set once using the store's canonical chromosome and allele ordering.
3. Assign variant indices from that sorted variant order.
4. Write `variants.tsv.gz` as BGZF, recording one virtual offset per row into `variant_offsets.npy`.
5. Build the tabix index over chromosome and position.
6. Write compact analysis metadata and any small alias records to SQLite.
7. Populate dense `z` and `se` arrays by mapping canonical ALIDs to assigned variant indices during the array fill step.
8. Build existing dense top-hit indexes using the final row indices and analysis indices.

Queries become layout-independent over a small `VariantAxis` interface rather than depending directly on SQLite rows:

- Exact canonical ALID lookup parses `chr:pos:A1:A2`, fetches `chr:pos-pos` through tabix, and filters by the two alleles.
- RSID or other alias lookup uses SQLite only to resolve aliases to canonical variant identity or variant index, then materialises the variant through tabix or row offset.
- Genomic range lookup streams matching variant rows from tabix and reads the corresponding dense array rows.
- Analysis extraction reads one dense array column and streams variant metadata in row order when row materialisation is requested.
- Sparse lookup resolves requested variants through tabix/aliases and requested analyses through SQLite, then reads the dense block from Zarr.
- Top-hit lookup uses top-hit row indices plus `variant_offsets.npy` to materialise only the required variant rows.

The required new runtime dependency is a Python binding for BGZF and tabix operations, preferably `pysam`, so the implementation does not depend on shelling out to external `bgzip` and `tabix` binaries. The dependency update should include any necessary type-checking configuration if the package does not provide complete typing metadata.

## User Stories

1. As an OpenGWASDB builder, I want the dense variant table stored outside SQLite, so that high-cardinality variant metadata does not dominate store size.
2. As an OpenGWASDB builder, I want variant rows written in canonical sorted order, so that dense row indices are deterministic within a release.
3. As an OpenGWASDB builder, I want to assign store-local variant indices while writing the variant table, so that Zarr arrays and query results share one compact row coordinate.
4. As an OpenGWASDB builder, I want to bgzip the variant table during write, so that the store does not need to hold an uncompressed variant TSV on disk after build.
5. As an OpenGWASDB builder, I want to create a tabix index as part of the store build, so that range queries work without post-processing.
6. As an OpenGWASDB builder, I want failures during variant-table writing to leave no partial valid-looking store, so that broken releases are easy to detect and clean up.
7. As an OpenGWASDB builder, I want row offsets recorded while writing BGZF, so that later query paths can materialise variants by dense row index without a SQLite variant table.
8. As an OpenGWASDB builder, I want the manifest to declare the variant-axis format, so that future readers can reject unsupported store variants clearly.
9. As an OpenGWASDB builder, I want analysis metadata to remain in SQLite, so that analysis lookup stays simple and compact.
10. As an OpenGWASDB builder, I want rsid aliases stored separately from canonical ALID identity, so that alias collisions do not corrupt the canonical variant axis.
11. As an OpenGWASDB builder, I want the dense Zarr arrays to keep their existing compression and chunking defaults, so that this experiment isolates variant-store changes from statistic-array changes.
12. As an OpenGWASDB builder, I want the benchmark suite to report variant table, tabix index, row-offset, SQLite, Zarr, and total store sizes separately, so that storage regressions can be attributed quickly.
13. As an OpenGWASDB builder, I want build timing split between normalisation, variant-axis writing, SQLite metadata writing, dense array writing, and top-hit indexing, so that performance bottlenecks are visible.
14. As an OpenGWASDB maintainer, I want a deep variant-axis module with a small public interface, so that query code does not care whether variants live in SQLite, tabix, or a future backend.
15. As an OpenGWASDB maintainer, I want the existing query facade to keep its public methods, so that callers are not forced to change for this storage experiment.
16. As an OpenGWASDB maintainer, I want validation to verify the variant table, tabix index, offset sidecar, SQLite metadata, and Zarr shapes together, so that invalid stores fail before query time.
17. As an OpenGWASDB maintainer, I want canonical ALID lookup to avoid a full ALID index, so that exact lookup does not recreate the SQLite footprint problem.
18. As an OpenGWASDB maintainer, I want top-hit queries to materialise variants directly from row offsets, so that top-hit latency remains close to the current vertical slice.
19. As an OpenGWASDB maintainer, I want analysis extraction to stream variant metadata in row order, so that full-analysis export stays memory-conscious.
20. As an OpenGWASDB maintainer, I want range queries to use tabix, so that genomic slicing is handled by a proven genomic index rather than a general SQL B-tree.
21. As an OpenGWASDB maintainer, I want rsid lookup semantics documented, so that ambiguous rsids can be handled deliberately rather than accidentally.
22. As an OpenGWASDB maintainer, I want the old SQLite variant table assumptions removed from tests, so that tests protect the new contract rather than the previous implementation.
23. As an OpenGWASDB query user, I want to query a variant by canonical ALID, so that PheWAS-style lookup still works.
24. As an OpenGWASDB query user, I want to query a variant by rsid when available, so that common user-facing identifiers still work.
25. As an OpenGWASDB query user, I want to query a genomic range, so that regional association views still work.
26. As an OpenGWASDB query user, I want to retrieve all finite associations for one analysis, so that full-study extraction remains available.
27. As an OpenGWASDB query user, I want to retrieve a small variant-by-analysis lookup set, so that targeted extraction remains fast.
28. As an OpenGWASDB query user, I want top-hit queries to return full variant metadata, so that top-hit indexes remain useful without a separate join step.
29. As an OpenGWASDB query user, I want missing dense cells to remain excluded from materialised results, so that query semantics do not change with the storage backend.
30. As an OpenGWASDB query user, I want derived beta and p-value fields to remain present in result rows, so that downstream code sees the same result contract.
31. As an OpenGWASDB release validator, I want to count variants from the variant-axis metadata rather than a SQLite variants table, so that validation matches the new physical layout.
32. As an OpenGWASDB release validator, I want to confirm row offsets point to the expected variant indices, so that random row materialisation cannot silently drift.
33. As an OpenGWASDB release validator, I want to confirm tabix range fetches agree with sequential variant-table reads on representative regions, so that genomic indexing is trustworthy.
34. As an OpenGWASDB release validator, I want to confirm dense array row counts match the variant table row count, so that array and metadata axes cannot diverge.
35. As an OpenGWASDB release validator, I want to confirm alias records point to existing variants, so that alias lookup cannot return orphaned rows.
36. As an OpenGWASDB developer, I want a dependency-managed Python tabix implementation, so that builds work in local and CI environments without hidden system binaries.
37. As an OpenGWASDB developer, I want tests for small fixtures and realistic benchmark-sized stores, so that correctness and performance are both visible.
38. As an OpenGWASDB developer, I want the experiment isolated on a development branch, so that the current vertical slice remains recoverable while the alternative is evaluated.

## Implementation Decisions

- The dense statistic arrays remain the source of association values. This PRD changes the variant metadata axis and lookup structures only.
- The Store Variant Table moves from a high-cardinality SQLite table to a bgzipped TSV plus tabix index.
- SQLite remains in the release, but its role is narrowed to low-cardinality metadata, analysis lookup, key/value metadata, and optional alias information.
- Store-local variant indices remain mandatory and continue to connect variant metadata to dense array rows.
- The variant table must include `variant_index` explicitly rather than relying only on physical line number.
- The variant table must include canonical ALID, chromosome, position, effect allele, other allele, and optional rsid.
- The variant table must be sorted by the same canonical order used to assign dense row indices.
- A compact row-offset sidecar is part of the proposed format so arbitrary dense row indices can be materialised without a SQLite variant table.
- The row-offset sidecar should store BGZF virtual offsets using a fixed-width numeric dtype.
- Canonical ALID exact lookup should parse the ALID into genomic coordinates and alleles, tabix-fetch that single position, and filter allele columns.
- Alias lookup should not require duplicating every canonical ALID in SQLite.
- Ambiguous aliases should be represented explicitly. The initial implementation may return all matching variants for an alias internally while the public single-identifier API can preserve its current first-match behaviour only if documented and tested.
- The query facade should depend on a variant-axis abstraction that can return variants by identifier, range, row index, row-index batch, and sequential row order.
- Query result construction should accept a variant-record object or mapping rather than requiring a SQLite row.
- Validation should no longer derive variant count from a SQLite variants table for dense stores using the new format.
- The manifest should identify the variant-axis backend and enough format version information for future migrations.
- The implementation should prefer Python library calls for BGZF and tabix operations over shelling out to command-line binaries.
- The benchmark suite should compare the old and new physical layouts using the same source inputs and query selections.
- Success should be judged by total store size, build time, and query latency together, not by any single metric.
- The first acceptance target is to remove the large SQLite variant table and its ALID/range indexes from the dense store footprint.
- The second acceptance target is to keep regional, PheWAS, top-hit, random lookup, and bulk-analysis timings within an acceptable range of the current vertical slice, with any regressions explained by benchmark output.

## Testing Decisions

- Tests should assert external behaviour through build, validation, and query APIs rather than checking private implementation details.
- Unit tests should cover the variant-axis module using tiny BGZF/tabix fixtures with canonical ALID lookup, range lookup, row-offset lookup, sequential streaming, and missing identifiers.
- Unit tests should cover alias semantics, including rsid aliases, absent aliases, and ambiguous aliases.
- Integration tests should update the current dense vertical-slice fixture so the store validates and all existing query behaviours pass with the new variant backend.
- Validation tests should cover missing variant table, missing tabix index, missing or wrong-length row-offset sidecar, row offsets pointing to the wrong variant index, duplicate canonical variants, and dense array shape mismatch.
- Build tests should verify the written variant indices match dense array row order.
- Build tests should verify the implementation does not insert every canonical ALID into the alias table.
- Query tests should verify canonical ALID lookup works through tabix without a SQLite variants table.
- Query tests should verify range lookup preserves genomic ordering and excludes missing dense cells.
- Query tests should verify top-hit materialisation returns the same public result fields as before.
- Query tests should verify full-analysis extraction can return dense arrays without materialising every variant row.
- Benchmark tests should report storage breakdown for source inputs, Zarr arrays, SQLite metadata, variant TSV, tabix index, row-offset sidecar, top-hit arrays, and total store.
- Benchmark tests should report timing breakdown for build phases and the existing query families.
- Prior art for these tests is the current dense vertical-slice test suite, store validation tests, source normalisation tests, and benchmark scripts.

## Out of Scope

- Changing dense statistic compression, dtype, or chunk shape defaults.
- Adding reference-completed imputation support.
- Reworking ragged layout storage.
- Replacing Zarr for statistic arrays.
- Introducing a remote service or catalogue layer.
- Adding variant-scoped EAF or INFO unless its provenance and scope are explicitly established.
- Guaranteeing variant indices are stable across releases.
- Solving global cross-release variant identity beyond the existing reference-assembly plus ALID contract.
- Optimising every query path before the storage-footprint experiment has been benchmarked.

## Further Notes

The besdq prototype demonstrated that bgzipped, tabix-indexed variant metadata is a plausible better fit for the dense variant axis than a large relational table. OpenGWASDB should adopt that idea while avoiding a large object-array key sidecar as the primary lookup structure. Tabix should handle genomic lookup; compact BGZF row offsets should handle dense row-index materialisation.

The current benchmark shows the compressed statistic arrays are not the main footprint problem. The new PRD therefore focuses on removing duplicated variant identity and range indexes from SQLite while preserving the public query contract.

The implementation should update the store-format documentation and ADRs once benchmark results confirm the new structure is the preferred dense variant-axis backend.
