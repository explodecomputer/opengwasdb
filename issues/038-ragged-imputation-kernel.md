## Parent PRD

`issues/prd-ragged-reference-completion.md`

## What to build

Port the elastic-net imputation math from pleiodb into `opengwasdb/layouts/ragged/impute.py` as a set of pure functions with no filesystem I/O. Add the scalar-N SE derivation alongside the ported Z kernel. This slice has no zarr or store dependencies — it is purely a numeric library.

Functions to expose:
- `ld_pca(ld_matrix, thresh)` — eigendecompose LD matrix, return (eigenvalues, eigenvectors) truncated to cumulative variance threshold
- `elastic_net_impute(z, eigenvectors, n_comp)` — fit ElasticNetCV on observed z, predict for all positions; return array or None
- `poly_rescale(truth, predicted, npoly)` — polynomial rescaling with Cook's-distance outlier removal; return (rescaled, pearson_r)
- `impute_z_block(z_obs_dense, eigenvectors)` — orchestrate PCA + elastic net + rescale for one LD block; return (z_imp, pearson_r)
- `scalar_n_se(se_obs, eaf_obs, eaf_ref)` — compute se_scale from observed, return se array for reference-panel positions

## Acceptance criteria

- [ ] `opengwasdb/layouts/ragged/impute.py` exists with all five public functions
- [ ] `ld_pca` returns eigenvalues in descending order, all non-negative, cumulative variance ≥ thresh
- [ ] `elastic_net_impute` returns None when < 2 observed or all z identical; otherwise returns full-length array with no NaN
- [ ] `poly_rescale` returns (array, pearson_r); pearson_r is NaN when < npoly+1 clean points
- [ ] `scalar_n_se` returns NaN for positions where EAF is 0 or 1
- [ ] Unit tests in `tests/test_ragged_impute.py` covering all five functions; prior art in `pleiodb/tests/test_impute.py`
- [ ] No dependency on pleiodb at import time (functions are self-contained)

## Blocked by

None — can start immediately.

## User stories addressed

- User story 13 (inherit kernel from pleiodb)
- User story 12 (scalar-N SE model)
- User story 15 (quality gate pearson_r)
