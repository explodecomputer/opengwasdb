## Parent PRD

`issues/prd-query-variant-axis-performance.md`

## What to build

Extend `write_variant_axis()` to write two additional sidecar files alongside
`variants.tsv.gz`:

- `variant_alid_bytes.npy` — numpy array of shape `(n_variants,)`, dtype `|S32`,
  containing ALID strings as fixed-width bytes, **sorted lexicographically**.
- `variant_alid_rows.npy` — numpy array of shape `(n_variants,)`, dtype `int32`,
  containing the variant index corresponding to each sorted ALID entry.

The sort is by ALID bytes so `np.searchsorted` can be used for O(log n) lookup.
Both arrays are written in a single pass over the variant list — sort once, write both.

No changes to build callers (`build.py`, `build_vcf.py`) are needed; they both go
through `write_variant_axis()`.

## Acceptance criteria

- [ ] After any `write_variant_axis()` call, `variant_alid_bytes.npy` and
      `variant_alid_rows.npy` exist in the store root
- [ ] `variant_alid_bytes.npy` is sorted (assert `np.all(arr[:-1] <= arr[1:])`)
- [ ] `variant_alid_rows.npy[i]` is the variant index of the ALID at `variant_alid_bytes.npy[i]`
- [ ] Length of both arrays equals the number of variants written
- [ ] Existing tests pass (the small fixture stores now get the sidecar too)

## Blocked by

None — can start immediately.

## User stories addressed

- User story 6 (build path writes sidecar automatically)
