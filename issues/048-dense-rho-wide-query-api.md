## Parent PRD

`issues/prd-dense-rho-matrix.md`

## What to build

Complete the Rho query surface on the facade, on top of the `data.zarr/rho`
group written by `issues/047-dense-rho-tracer-bullet.md`:

- **`rho_row(analysis_id)`** — one Analysis vs all others, returning aligned
  arrays of the other Analysis IDs, `rho`, and `n_null` support.
- **`rho_matrix(ids=None)`** — wide format: the full symmetric `n_analyses ×
  n_analyses` matrix, or the dense submatrix block for a given vector of IDs.

Both unpack the stored strict lower triangle into symmetric form with diagonal
Rho = 1.0, and return the matching support (PRD "Query API"). ID resolution
reuses the facade's existing analysis-id → index mapping.

## Acceptance criteria

- [ ] `rho_row(a)` returns Rho and support of `a` against every other Analysis,
      with `a`'s self entry excluded (or reported as 1.0/self — documented).
- [ ] `rho_matrix()` returns the full symmetric matrix with a diagonal of 1.0 and
      `rho[i,j] == rho[j,i]`, matching the packed lower-triangle values.
- [ ] `rho_matrix(ids)` returns the correct pairwise submatrix (and support) for
      an arbitrary vector of IDs, in the given ID order.
- [ ] Unknown / missing Analysis IDs are handled with a clear error or documented
      empty result, consistently with the other facade methods.
- [ ] Results agree cell-for-cell with `rho(*ids)` long format from issue 047.

## Blocked by

- Blocked by `issues/047-dense-rho-tracer-bullet.md`

## User stories addressed

- User story 3
- User story 4
