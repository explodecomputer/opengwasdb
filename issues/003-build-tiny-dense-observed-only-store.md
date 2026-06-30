## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Build a tiny Dense Observed-Only Store Release end to end from deterministic fixture inputs. The output must contain `manifest.json`, `index.sqlite`, and `data.zarr/` with `z` and `se` arrays, using the source-faithful dense variant axis required by the v0.1 specification.

Use `/Users/gh13047/repo/besdq/data/ebi_input/` as the realistic Dense-mode source reference, but the test fixture used by this issue should live in OpenGWASDB so the test suite is self-contained. Do not use or compare the sparse `38714679` fixture for the Dense vertical slice.

## Acceptance criteria

- [ ] A build command or build API writes the standard store envelope to an explicit output path.
- [ ] SQLite contains variant metadata, analysis metadata, and lookup indexes needed for later queries.
- [ ] Zarr contains `z` and `se` arrays with shape `n_variants x n_analyses`.
- [ ] Missing source associations are represented by canonical NaN in both `z` and `se`.
- [ ] Store-local variant indices connect SQLite variant rows to Zarr array rows.
- [ ] A built tiny store passes the minimal validator from `issues/001-minimal-store-release-open-info-validate.md`.
- [ ] Tests verify array values, NaN cells, variant ordering, analysis metadata, and manifest fields.
- [ ] Test fixture provenance documents whether rows were generated or reduced from `besdq/data`.

## Blocked by

- Blocked by `issues/001-minimal-store-release-open-info-validate.md`
- Blocked by `issues/002-source-row-reader-and-canonical-variant-normalisation.md`

## User stories addressed

- User story 1
- User story 2
- User story 5
- User story 8
- User story 9
- User story 10
- User story 11
- User story 12
- User story 15
- User story 17
- User story 38
- User story 39
