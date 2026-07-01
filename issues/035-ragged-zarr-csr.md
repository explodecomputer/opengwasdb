## Parent PRD

`issues/prd-ragged-observed-only.md`

## What to build

The zarr CSR infrastructure for the Ragged Layout: a writer that accumulates per-analysis association sequences and a reader that slices them by analysis index.

Physical layout inside `data.zarr/ragged/`:

```
offsets        [n_analyses + 1]  int64   — offsets[i]:offsets[i+1] is the slice for analysis i
variant_index  [n_total_assoc]   int32   — store-local variant row indices
z              [n_total_assoc]   float16
se             [n_total_assoc]   float16
```

All four arrays are written together after all analyses are known. The implementation lives in `opengwasdb/layouts/ragged/zarr_csr.py`.

### Writer

`RaggedCSRWriter` accepts analyses one at a time (streaming), buffers them, and flushes to zarr on `close()` or as a context manager. Internally it accumulates lists of arrays and concatenates at flush time.

### Reader

`RaggedCSRReader` wraps the zarr group and exposes:

```python
def get_analysis(self, analysis_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (variant_index, z, se) for the given analysis."""
```

Access pattern: read `offsets[analysis_index:analysis_index+2]` (two int64 values), then slice the three payload arrays. Both operations are O(1) zarr reads with no Python-level deserialization.

## Acceptance criteria

- [ ] `RaggedCSRWriter` accumulates (variant_index, z, se) per analysis and writes valid zarr arrays on close
- [ ] `RaggedCSRReader.get_analysis(i)` returns correct numpy arrays for any valid index
- [ ] Empty analysis (0 associations) is handled: `offsets[i] == offsets[i+1]`, all three arrays return empty
- [ ] Round-trip test: write 1000 analyses with varying association counts (0 to 500), read back all, assert exact values
- [ ] Zarr arrays use float16 for z/se, int32 for variant_index, int64 for offsets
- [ ] Compression: zstd with bitshuffle filter (consistent with dense layout defaults)

## Blocked by

None — can start immediately.

## User stories addressed

- User story 1
- User story 3
- User story 4
