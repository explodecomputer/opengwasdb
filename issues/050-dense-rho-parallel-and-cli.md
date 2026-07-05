## Parent PRD

`issues/prd-dense-rho-matrix.md`

## What to build

Make the Rho Matrix build production-usable: parallelise the per-pair estimation
and expose it on the CLI. Builds on the serial `build_dense_rho` from
`issues/047-dense-rho-tracer-bullet.md`.

- **Parallel compute**: precompute the observed nulls-zeroed `Z` and null mask
  `M` at the thinned variants once; parallelise the per-pair MLE over pairs with a
  single-level `ProcessPoolExecutor` (`n_workers`), mirroring ADR 0023. The
  cheap sufficient-statistics (`A,B,C,n`) may be batched with linear algebra. The
  produced matrix must be identical to the serial path (PRD "Compute strategy").
- **CLI**: a `build-dense-rho <store_path>` command with options `--window-bp`,
  `--z-thresh`, `--min-nulls`, and `--n-workers`, writing the `rho` group in
  place (same add-in-place pattern as the top-hit index).

## Acceptance criteria

- [ ] `build_dense_rho(..., n_workers=N)` produces a byte-identical `rho`/`n_null`
      result to the serial (`n_workers=1`) path on a small store.
- [ ] `n_workers` has no effect on the computed values (pure runtime knob).
- [ ] BLAS thread usage is capped per worker so the pool is not oversubscribed
      (as in the completion pipeline).
- [ ] `opengwasdb build-dense-rho <store> --window-bp ... --z-thresh ...
      --min-nulls ... --n-workers ...` builds and writes the `rho` group, and the
      result validates (issue 049) and queries (issues 047/048).
- [ ] CLI defaults match the pleiodb defaults (`z_thresh=1.0`, `min_nulls=500`)
      and the PRD's `window_bp` default.

## Blocked by

- Blocked by `issues/047-dense-rho-tracer-bullet.md`

## User stories addressed

- User story 8
- User story 9
- User story 10
