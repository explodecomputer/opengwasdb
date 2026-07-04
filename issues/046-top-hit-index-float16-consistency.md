## Context

The Dense top-hit index is harvested from **float32** z-scores but the store keeps
z as **float16**. Near a threshold's `z_critical` cutoff the two disagree, so the
index and the stored matrix are inconsistent — and `validate_store` (correctly)
flags it.

Measured on the ukb-b store (9,847,701 × 2,514, built on branch
`dense-pass2-memory-streaming`):

- **p_5e_04** (`z_crit = 3.4808`): **48,839 of 41,926,521 index cells (0.12%)**
  have a harvested float32 `|z| >= z_crit` that rounds to `< z_crit` once stored
  as float16 → their stored value falls below the cutoff.
- **p_5e_06 / p_5e_08**: the reverse — matrix cells whose float16 value rounds
  *up* to `>= z_crit` but whose float32 value was below, so they were never
  harvested → the index is missing them (count mismatch).

Validation surfaces this as:
```
top-hit index p_5e_04 contains z value inconsistent with z array
top-hit index p_5e_06 does not match stored z values
top-hit index p_5e_08 does not match stored z values
```

Root cause: harvest membership is decided on the pre-storage float32 value
(`_resolve_column` / the `_TOP_HIT_Z_CRIT` mask in `build_vcf.py`), while the
store — and therefore any query reading the `z` array — sees the float16 value.
The index should be defined by the **stored** value so it is self-consistent with
what queries return.

## What to change

- [ ] **Harvest on the float16-rounded z.** In the Pass 2 harvest, threshold on
      `z.astype(float16).astype(float32)` and store that rounded value as the index
      `z`, so index membership and index `z` exactly match the stored matrix. Then
      the strict validation passes with no tolerance fudge.
- [ ] **Rebuild the existing ukb-b top-hit index** consistently — band-stream the
      float16 matrix (reuse the validation/`_write_dense_bands` band pattern),
      threshold on the stored values, and rewrite `data.zarr/top_hits/*`. This
      avoids a full ~3.5h rebuild; only the index is regenerated (~10 min, low
      memory).
- [ ] Add a test: a store whose float32 z sits within a float16 ULP of a cutoff
      round-trips through build → validate cleanly.

## Alternative (rejected)

Loosening validation to tolerate a float16-ULP band around each cutoff would let
the store pass, but it hides real inconsistencies and leaves the index disagreeing
with what `analysis`/`phewas` queries return from the `z` array. Defining the
index by the stored value (above) is the correct, exact fix.

## Notes

- Functionally minor: ~0.12% of the loosest tier near the cutoff; the store's
  actual z/se data is correct. But `validate_store` fails until the index is made
  consistent, so `build.sh` cannot report a clean pass.
- Distinct from [[043-dense-vcf-build-memory-footprint]] /
  [[045-dense-validation-memory-footprint]] (those are memory; this is
  correctness). Surfaced only because the 045 fix made validation actually finish.
