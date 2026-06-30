## Problem Statement

The current query facade materialises every association result as a Python `AssociationResult` dataclass with 17 fields (including computed `beta` and `p_value`) regardless of result size. At chr1 × 100-analysis scale, a bulk analysis query materialises 762,080 Python objects and a regional query materialises 324,048. Benchmark results show this is 28–562× slower than the besdq prototype, which returns numpy arrays directly. The row-materialisation path is also the only public query contract, so callers cannot avoid the overhead.

## Solution

Replace all five query methods (`analysis`, `phewas`, `range`, `lookup`, `top_hits`) with a sparse flat-array return contract. Each method returns `dict[str, np.ndarray]` with four parallel arrays — `variant_index`, `analysis_index`, `z`, `se` — covering only finite cells. Callers access variant and analysis metadata separately via new `variants_table()` and `analyses_table()` accessors on `StoreQuery`. The `AssociationResult` dataclass and the existing `analysis_arrays()` companion method are removed.

## User Stories

1. As a researcher, I want to extract all associations for one analysis as numpy arrays, so that I can feed them into downstream statistical tools without per-row Python overhead.
2. As a researcher, I want a PheWAS result as sparse arrays, so that I can filter and sort across analyses without constructing Python objects.
3. As a researcher, I want a genomic range query to return sparse arrays, so that I can overlay results from multiple analyses in a window efficiently.
4. As a researcher, I want top-hit associations as sparse arrays, so that I can post-process significance thresholds without object allocation overhead.
5. As a researcher, I want a random variant × analysis lookup as sparse arrays, so that batch extraction is fast regardless of result size.
6. As a developer, I want variant and analysis metadata available as a table (keyed by index), so that I can join sparse array results to ALIDs and analysis identifiers when needed.
7. As a developer, I want all five query methods to return results in the same sparse format, so that downstream code that combines query types does not need to handle multiple result shapes.
8. As a developer, I want the sparse arrays to use float32 for z/se and int32 for indices, so that memory usage matches the top-hit index format already stored in zarr.
9. As a developer, I want the existing test suite to validate query correctness against the sparse array contract, so that correctness guarantees are not lost in the migration.
10. As a developer, I want `analysis_arrays()` removed, so that there is only one bulk-analysis query path and no confusion about which to use.

## Implementation Decisions

- All five query methods return `dict[str, np.ndarray]` with keys `variant_index` (int32), `analysis_index` (int32), `z` (float32), `se` (float32). Only finite cells are included.
- `top_hits` reads the pre-indexed zarr arrays directly without per-row SQLite lookups; the other four methods compute finite masks from zarr slices and build sparse index arrays.
- Two new public accessors are added to `StoreQuery`: `variants_table()` returns a list of dicts (or equivalent) keyed by `variant_index`, and `analyses_table()` does the same for `analysis_index`. These are light SQLite reads.
- `AssociationResult` and `analysis_arrays()` are removed entirely. No deprecation wrapper.
- Existing tests are rewritten against the new contract. Test assertions check array shapes, dtypes, and specific finite-cell values rather than object attributes.
- The `beta` and `p_value` derivations (`beta_from_z_se`, `p_value_from_z`) are no longer called on the query path; they remain in `opengwasdb.stats` for callers that need them.
- An ADR is written to record why Python-object row materialisation was removed (performance at chr1+ scale is the trigger; the top_hit index format is the precedent).

## Testing Decisions

- Good tests check the array output directly: shape `(k,)`, dtype, and spot-check known finite cells from the synthetic fixture.
- The conftest `dense_store_path` fixture remains unchanged; only the assertions in the test files change.
- Prior art: `test_dense_vertical_slice.py` — rewrite its query assertions to work with sparse arrays.
- `test_vcf_source.py` and `test_liftover.py` are unaffected.
- `test_dense_vcf_build.py` query assertions must be rewritten.

## Out of Scope

- bcftools-based VCF reader (separate PRD).
- Streaming / lazy result objects.
- DataFrame (pandas/polars) return types.
- Any change to the zarr layout, SQLite schema, or manifest format.
- The `analysis_arrays()` method signature; it is simply removed, not migrated.

## Further Notes

- ADR-0006 (layout-independent query engine) documents the facade contract; the new array return type is a change to that contract.
- The top-hit zarr group already stores `variant_index`, `analysis_index`, `z`, `se` arrays — the new query contract mirrors this format exactly.
- Benchmark target: chr1 × 100-analysis query timings should be within 2× of besdq zstd_bitshuffle baseline after this change.
