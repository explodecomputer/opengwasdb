## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Reduce the numpy dtype of `variant_alid_bytes.npy` from `|S64` to `|S32`.

The sidecar is a sorted array of fixed-width byte strings used for O(log n)
`np.searchsorted` ALID lookups.  It is currently written with `_ALID_DTYPE = "|S64"`,
allocating 64 bytes per entry regardless of actual string length.

Measured against the full-genome ukb-b smoke test store (9.85M variants):

| Metric | Value |
|---|---|
| Max ALID length | 31 chars (`14_GL000009v2_random:104542:A:G`) |
| p99.9 ALID length | 16 chars |
| Mean ALID length | 14.6 chars |
| Current dtype | `\|S64` → 64 bytes/entry → **602 MB** |
| Target dtype | `\|S32` → 32 bytes/entry → **~301 MB** |
| Saving | ~300 MB per store |

`|S32` is the minimum safe width: it fits the longest observed ALID (31 chars)
with 1 byte to spare.  The canonical ALID format (`chrom:pos:ref:alt`) will not
produce longer strings for standard chromosomes; the 31-char outlier is an
unplaced contig name.

The change is in `opengwasdb/variants/axis.py`: update the `_ALID_DTYPE` constant
and any tests that assert on the dtype.  Stores built before this change remain
readable — numpy pads shorter dtypes on load, so old `|S64` files loaded with a
`|S32` reader will silently truncate; existing stores should be rebuilt.

## Acceptance criteria

- [ ] `_ALID_DTYPE` constant changed to `"|S32"`
- [ ] `variant_alid_bytes.npy` written at 32 bytes/entry in new stores
- [ ] `np.searchsorted` lookups still return correct results
- [ ] Existing test in `test_dense_vertical_slice.py` passes unchanged
- [ ] `variant_alid_bytes.npy` in a rebuilt ukb-b-test10 store is ≤ 310 MB

## Blocked by

None — can start immediately.

## User stories addressed

- Reduces fixed per-store overhead by ~300 MB at no query-performance cost.
