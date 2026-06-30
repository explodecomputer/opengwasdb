## Problem Statement

After implementing the sparse array query contract (ADR-0020) and fixing zarr chunk shape
(ADR-0021), two query patterns remain substantially slower than the besdq prototype:
regional (15.6×) and random_lookup (3.76×). Profiling has isolated the root causes to the
VariantAxis layer, not zarr reads.

**Regional** (78ms ours vs ~5ms besdq): `VariantAxis.range()` spends 51.6ms constructing
5 463 frozen Python `VariantRecord` objects from tabix lines. The facade's `range()` query
uses only `variant_index` from each object, so the full parse is wasted work. Raw tabix
fetch costs only 2.3ms; zarr oindex for the same rows costs 25ms. The Python object
allocation is the dominant cost.

**random_lookup** (114ms ours vs 30ms besdq): Each call to `by_identifier()` does one
pysam tabix fetch per identifier (100 sequential fetches × ~0.15ms = 15ms). The
remaining ~95ms is zarr scattered-chunk access, which is irreducible given our two-array
(z + se) contract. The identifier lookup overhead is avoidable.

The bulk gap (1.68×) is structural: we read two zarr arrays (z, se) where besdq's
benchmark reads one (beta). This is correct per ADR-0020 and is not addressed here.

## Solution

Two targeted fixes to `VariantAxis`:

1. **`range_indices()` method** — returns `np.ndarray[int32]` of variant indices for a
   genomic range by parsing only the index column from tabix output, with no `VariantRecord`
   allocation. The facade's `range()` query uses this internally.

2. **mmap'd sorted ALID sidecar** — at build time, write a sorted numpy array of
   fixed-width ALID byte strings paired with their variant indices as a sidecar file
   alongside `variants.tsv.gz`. At open time, `VariantAxis` memory-maps this file.
   `by_identifier()` uses `np.searchsorted` for O(log n) lookup with near-zero startup
   cost and no explicit preload. The OS page cache handles warmup naturally.
   rsid/alias lookups continue to fall back to SQLite as before.

Both fixes work correctly for CLI use (no startup penalty) and service deployments
(pages warm up across requests) without an explicit preload flag.

## User Stories

1. As a researcher running regional queries via the CLI, I want them to complete in under
   10ms for typical 1 Mb windows, so that interactive exploration is responsive.
2. As a researcher running phewas queries via the CLI, I want sub-millisecond identifier
   resolution, so that the lookup overhead is negligible.
3. As a service operator, I want identifier lookups to use OS-cached memory-mapped data
   without explicit preloading, so that the VariantAxis scales without additional
   configuration.
4. As a developer, I want the benchmark to reflect accurate per-pattern timings comparable
   to besdq, so that regressions are detectable.
5. As a developer, I want the new sidecar file to be validated at store-open time, so that
   stores built without it fail loudly rather than silently falling back to a slow path.
6. As a developer building a store, I want `write_variant_axis()` to write the mmap sidecar
   automatically, so no build-path code needs updating.
7. As a developer, I want `range_indices()` to return a numpy array directly, so that the
   facade never allocates Python variant objects during query execution.

## Implementation Decisions

- **`range_indices(chromosome, start, end) → np.ndarray[int32]`** added to `VariantAxis`.
  Parses only the `variant_index` column (column 2) from tabix output. The existing
  `range()` method (returning `list[VariantRecord]`) is kept for callers that need
  full metadata; the facade switches to `range_indices()` internally.

- **Sidecar filename**: `variant_alid_index.npy` stored at the store root alongside
  `variants.tsv.gz`. Contains a structured numpy array sorted by ALID bytes.

- **ALID encoding**: fixed-width `|S32` byte strings — sufficient for canonical ALID
  format (`chr:pos:A1:A2`, max length ~30 chars). Sorted lexicographically to enable
  `np.searchsorted`.

- **Two parallel arrays** rather than a structured array: `alid_bytes.npy` (sorted `|S32`)
  and `alid_indices.npy` (`int32`), so `searchsorted` operates on a contiguous byte array
  without strided access.

- **`by_identifier()` flow**: parse ALID → binary search on mmap'd bytes → return
  `VariantRecord` via `by_index()` (uses offset sidecar, already fast). rsid/alias path
  unchanged (SQLite fallback).

- **`by_alid()` flow**: no longer calls `range()` internally; uses the mmap index directly.

- **Validation**: `validate_store()` checks `variant_alid_index.npy` exists and its
  length matches the variant count. Stores built before this change are invalid until
  rebuilt.

- **`write_variant_axis()`** writes both sidecar files in one pass — no API change to
  build callers.

- The benchmark region in `benchmark_vcf_ukb_chr1_dense.py` currently selects chr1:4-5Mb
  (5 463 variants) while besdq uses chr1:100-101Mb (2 751 variants). After these fixes
  the benchmark should use a consistent region selection methodology and note the variant
  count alongside timing, so cross-benchmark comparisons are fair.

## Testing Decisions

- Good tests exercise the public contract: correct variant indices returned, correct
  sparse dict output, correct behaviour on empty regions and unknown identifiers.
- Tests should not inspect internal sidecar files or verify that a particular code path
  was taken — assert on outputs, not mechanism.
- `test_dense_vertical_slice.py` and `test_dense_vcf_build.py` are the primary
  integration test targets; the existing fixtures are small enough to run fast.
- Add one focused unit test for `range_indices()` asserting it returns the same indices
  as `[v.variant_index for v in range()]`.
- Stores built by earlier test fixtures will be missing the new sidecar; tests must use
  the updated build path (which happens automatically via `write_variant_axis()`).

## Out of Scope

- Reducing the bulk gap (1.68×): inherent in reading two arrays per our contract.
- Parallelising zarr chunk reads across threads.
- Reducing the phewas gap (1.57×) beyond what the mmap index already provides.
- Any change to the sparse array return contract (ADR-0020).
- Preloading the full variant table into a Python dict.

## Further Notes

- ADR-0020: sparse array query contract
- ADR-0021: zarr chunk shape clipping
- besdq prototype uses `load_keys()` (full in-memory dict) — we deliberately choose mmap
  to avoid the CLI startup penalty at 12M-variant scale.
- `variant_offsets.npy` (existing) covers random-access by index; the new sidecar covers
  random-access by ALID.
