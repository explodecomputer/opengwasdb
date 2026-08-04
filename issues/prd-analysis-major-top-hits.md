# Analysis-Major Top-Hit Index Experiment

## Problem Statement

OpenGWASDB's precomputed top-hit index makes genome-wide significance queries
much cheaper than scanning the association matrix, but its physical ordering is
not suited to the most common lookup: retrieving significant associations for
one analysis. The current Dense query API has no analysis selector, and the
global index must be materialised and filtered in memory. On the full UK
Biobank Dense stores, obtaining 7,390–8,580 hits for one analysis consequently
takes roughly 300–400 ms because 5.5–6.4 million global index entries are read.

This work is an experiment. It must establish whether an analysis-addressable
top-hit representation produces a substantial latency improvement before that
representation becomes mandatory, migration tooling is designed, or the work
is expanded to the Ragged and Hybrid layouts.

## Solution

Add an experimental analysis-major top-hit index for Dense stores. Each of the
existing significance-threshold tiers will store associations grouped by
analysis, with a compact offset table identifying the contiguous slice for each
analysis. Within an analysis, associations will be ordered by canonical genomic
position. A global read will concatenate analyses in analysis-index order and
retain genomic ordering within each analysis.

Extend the layout-independent query contract with an optional `analysis_id`
selector on top-hit queries. When supplied for a Dense store, the query will
resolve the analysis through the existing analysis metadata and read only that
analysis's indexed slice. An unknown analysis identifier will follow the
existing query-facade convention and return an empty sparse-array result. Calls
without an analysis selector will remain supported; their precise global
ranking is not an optimization target for this experiment.

The index will continue to carry the values needed to satisfy the complete
public association result contract without falling back to random matrix reads.
Reference-completed Dense stores must retain indexed imputation status, and the
analysis-filtered query must support the existing observed-only behavior.

The experiment will benchmark observed-only and reference-completed full UK
Biobank Dense stores. Its warm-cache target is under 10 ms for the representative
single-analysis query at the genome-wide significance tier. Results must also
show the improvement relative to the current global-index filtering path and
report cold or first-read behavior separately where practical.

The work has two gated phases:

1. Implement, validate, and benchmark the analysis-major index for Dense
   observed-only and reference-completed stores.
2. If the Dense benchmark demonstrates a substantial improvement, separately
   plan and implement equivalent query behavior for Ragged and Hybrid stores.

Whether the representation becomes mandatory for new stores, whether old
stores require an in-place rebuild command, and the exact Ragged and Hybrid
physical representations will be decided only after reviewing the Dense
experiment.

## User Stories

1. As an OpenGWASDB query user, I want to request top hits for one analysis, so that I do not materialise millions of unrelated associations.
2. As an OpenGWASDB query user, I want to identify the analysis by its public analysis ID, so that I do not need to know its internal numeric coordinate.
3. As an OpenGWASDB query user, I want an unknown analysis ID to behave consistently with other query endpoints, so that missing identifiers are predictable.
4. As an OpenGWASDB query user, I want single-analysis hits returned in chromosome-position order, so that I can inspect and process them genomically.
5. As an OpenGWASDB query user, I want global top-hit queries to remain available, so that the experiment does not remove existing behavior.
6. As an OpenGWASDB query user, I want global results grouped by analysis-index order and genomic order within each analysis, so that their ordering is deterministic.
7. As an OpenGWASDB query user, I want the existing threshold choices retained, so that significance semantics do not change during the indexing experiment.
8. As an OpenGWASDB query user, I want the existing sparse-array result fields retained, so that downstream consumers need only opt into the new selector.
9. As an OpenGWASDB query user, I want a limit to apply to the index's chromosome-position order, so that its behavior remains deterministic even though limited queries are uncommon.
10. As an OpenGWASDB query user, I want observed-only filtering on a completed store to remain available, so that imputed associations can be excluded explicitly.
11. As an OpenGWASDB builder, I want significant associations grouped contiguously by analysis, so that one analysis can be read as a narrow slice.
12. As an OpenGWASDB builder, I want an offset table for every threshold tier, so that lookup does not scan an analysis-index array.
13. As an OpenGWASDB builder, I want offsets to cover analyses with zero hits, so that every valid analysis coordinate can be addressed directly.
14. As an OpenGWASDB builder, I want deterministic genomic ordering within each analysis, so that serial, parallel, harvested, and rebuilt indexes agree.
15. As an OpenGWASDB builder, I want the existing precomputed significance tiers retained, so that this experiment isolates physical ordering from threshold policy.
16. As an OpenGWASDB builder, I want observed and completed indexes to use the same addressing scheme, so that query behavior does not depend on completion state.
17. As an OpenGWASDB builder, I want completed indexes to include imputation flags, so that sliced reads do not trigger random access to the dense imputation matrix.
18. As an OpenGWASDB maintainer, I want top-hit writing and reading encapsulated behind a small module interface, so that physical index details do not spread through builders and query facades.
19. As an OpenGWASDB maintainer, I want existing global top-hit calls to remain valid, so that the experimental API extension is backward compatible.
20. As an OpenGWASDB maintainer, I want index validation to detect malformed offsets and incorrectly ordered slices, so that fast reads cannot return corrupt results silently.
21. As an OpenGWASDB maintainer, I want harvested and full-scan index construction to produce equivalent query results, so that build strategy does not affect behavior.
22. As an OpenGWASDB maintainer, I want analysis selection tested through the public query facade, so that tests protect user-visible behavior rather than private array access.
23. As an OpenGWASDB maintainer, I want the global and analysis-filtered results to agree for each analysis, so that the new fast path is demonstrably complete.
24. As an OpenGWASDB maintainer, I want chunking selected for narrow analysis slices, so that a small logical lookup does not decompress an unnecessarily large global chunk.
25. As an OpenGWASDB maintainer, I want index storage size reported even though space is not a primary constraint, so that the experiment's trade-off remains visible.
26. As an OpenGWASDB performance investigator, I want the benchmark to query a representative analysis rather than the full global index, so that it reflects usual use.
27. As an OpenGWASDB performance investigator, I want warm-cache latency reported separately, so that the under-10-ms target is assessed consistently.
28. As an OpenGWASDB performance investigator, I want first-read or cold behavior reported where practical, so that decompression and filesystem effects are not hidden.
29. As an OpenGWASDB performance investigator, I want the old global-read-and-filter timing retained as a comparison, so that the achieved speedup is explicit.
30. As an OpenGWASDB performance investigator, I want observed-only and completed Dense stores benchmarked, so that imputation metadata costs are measured.
31. As an OpenGWASDB decision maker, I want a documented success gate before changing the required store format, so that an experiment does not silently become a migration obligation.
32. As an OpenGWASDB decision maker, I want Dense evidence before expanding to Ragged and Hybrid, so that later layout work is justified by measured value.
33. As an OpenGWASDB decision maker, I want the Dense results recorded in the benchmark report, so that the adoption decision has reproducible evidence.
34. As a future Ragged or Hybrid implementer, I want the public analysis-filtered query semantics settled independently of the physical layout, so that each layout can choose an appropriate representation.

## Implementation Decisions

- Phase one applies to observed-only and reference-completed Dense stores.
- Ragged and Hybrid behavior is part of the intended feature but is gated on a successful Dense experiment; their physical implementation is not decided by this PRD.
- The public top-hit query gains an optional analysis-ID selector while preserving calls that omit it.
- Unknown analysis IDs return the standard empty sparse-array result, consistent with existing analysis-addressed query behavior.
- The existing precomputed significance-threshold tiers are retained.
- Each Dense threshold tier stores a packed analysis-major association sequence plus an offsets array of length analysis-count plus one.
- The offsets array is the direct addressing structure for per-analysis lookup; a query must not scan the complete analysis-index array to find a slice.
- Associations are ordered by analysis index globally and by canonical chromosome-position order within an analysis. Canonical variant order or variant index may implement the genomic ordering where the store contract guarantees their equivalence.
- A limit applies after selecting the requested scope and therefore returns the first entries in the defined analysis/genomic order. Strongest-first limited queries are not a requirement.
- Global top-hit ordering may change from significance ranking to deterministic analysis/genomic ordering. Exact global ordering is not a performance goal, but it must be documented and tested.
- The index continues to store variant coordinate, statistic, standard error, and any completion-status information required by the public result contract.
- The reader and writer should be encapsulated as a deep Dense top-hit-index module with a narrow interface for building a tier, reading all hits, and reading one analysis.
- Dense build paths that harvest hits inline and paths that rebuild from the matrix must use the same writer contract.
- Chunk sizing should be tuned for typical per-analysis slices rather than only global sequential reads. The experiment should compare sensible smaller chunks with the current large chunks rather than assuming a value without measurement.
- Phase one may introduce an experimental format marker or optional index capability. It must not make the representation mandatory for all future stores until the benchmark gate is reviewed.
- No in-place migration or rebuild command for existing stores is required during the experiment. That decision follows the benchmark.
- The primary success target is a warm-cache median below 10 ms for a representative single-analysis genome-wide-significant query on each full UK Biobank Dense store.
- Success also requires a substantial improvement over the current 300–400-ms global-index filtering path, with correctness and global-query regressions reported alongside latency.
- If the Dense gate succeeds, follow-up work will extend the public behavior to Ragged and Hybrid. The Dense representation must not be assumed to be optimal for those layouts.

## Testing Decisions

- Tests will assert observable build, validation, and query behavior rather than coupling themselves to private helper functions or exact storage library calls.
- Isolated index round-trip tests will cover analyses with many, one, and zero hits; every retained threshold; genomic ordering; deterministic ties; limits; and completed-store imputation flags.
- Dense build integration tests will cover both inline-harvested and full-matrix-rebuilt indexes and require their public query results to agree.
- Query-facade integration tests will cover a known analysis ID, an unknown analysis ID, omitted analysis selection, limits, observed-only filtering, and reference-completed results.
- Equivalence tests will compare each analysis-filtered result with the corresponding subset of the global result, without requiring the implementation to obtain both through the same physical path.
- Validation tests will cover missing offsets, wrong offset length, non-monotonic offsets, a final offset inconsistent with array lengths, slices assigned to the wrong analysis, and associations out of genomic order.
- Corruption tests should exercise the public validator and report actionable failures rather than asserting private validator structure.
- Performance evidence will come from the full UK Biobank Dense benchmark rather than timing-sensitive unit tests in CI.
- The benchmark will report warm-up policy, repetitions, median and tail latency, result count, chunk configuration, index size, and comparison with the previous global-filter path.
- The benchmark will cover the same representative analysis in observed-only and reference-completed stores.
- The warm-cache acceptance target is under 10 ms for the single-analysis query. First-read behavior will be recorded where practical but has no fixed acceptance threshold in this experiment.
- Prior art includes the existing Dense top-hit harvest/rebuild equivalence tests, Dense completion tests, query-facade vertical-slice tests, validation corruption tests, and the full UK Biobank Dense benchmark report.

## Out of Scope

- Making the experimental representation mandatory for every new store.
- Providing an in-place migration or index-rebuild command for existing releases.
- Removing or changing the existing significance-threshold tiers.
- Supporting arbitrary significance thresholds from a single loose index.
- Optimizing strongest-first limited queries.
- Changing the sparse-array association result contract beyond adding the optional analysis selector.
- Changing the dense association matrix representation, compression, precision, or general query chunking.
- Implementing the Ragged or Hybrid physical index before the Dense success gate is reviewed.
- Choosing the eventual Ragged or Hybrid index representation.
- Treating index storage space as a primary optimization constraint.
- Establishing a hard cold-cache latency target.

## Further Notes

The current top-hit arrays occupy only tens of megabytes on the full Dense UK
Biobank stores, so modest duplication or smaller chunks are acceptable if they
produce a clear latency win. Packed analysis-major arrays plus offsets are
preferred over one dataset or group per analysis because they provide direct
addressing without creating a large number of small Zarr objects.

The current index writer's ordering and comments should be reviewed carefully:
the experimental format must define ordering as part of its observable contract
and use a sort whose primary key is analysis index. Since Dense variant indices
follow canonical genomic order, they are a natural compact secondary key if
validation confirms that invariant.

The adoption decision should consider more than the headline median. A useful
experiment will also reveal the effect of Zarr chunk size, decompression,
completion-status reads, global-query behavior, build time, and validation cost.
If the Dense target is not met, the benchmark and profiling evidence should be
retained rather than proceeding automatically to Ragged and Hybrid.
