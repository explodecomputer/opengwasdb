## Parent PRD

`issues/prd-ragged-reference-completion.md`

## What to build

Implement `opengwasdb/layouts/ragged/ld_panel.py`: a read-only interface for the LD reference panel on disk. Covers block discovery, variant matching, LD matrix loading, and eigenvector cache reading. No imputation logic here — this is a pure I/O layer.

Functions to expose:
- `LDBlock(dataclass)` — `block_id`, `tsv_path`, `ld_path`, `ldeig_npz_path | None`, `n_ld_snps`, `chrom`, `start_bp`, `end_bp`, `snp_ids: list[str]`, `eaf: np.ndarray`
- `find_blocks(ld_dir, ancestry, chrom, start_bp, end_bp)` → `list[LDBlock]` — find all blocks whose genomic extent intersects the query window
- `match_variants(block, store_alids)` → `(variant_indices, ld_row_indices)` — match store ALIDs to LD panel rows; return parallel index arrays; require ≥ 2 matches to include a block
- `load_ld_eigenvectors(block, thresh)` → `(eigenvalues, eigenvectors)` — load from `.ldeig.npz` cache if present, else load `.unphased.vcor1.gz` and call `ld_pca`
- `load_ld_matrix(block, row_indices)` → `np.ndarray` — load submatrix for matched rows only

The flat panel layout (`ld_dir/{ancestry}/{chr}/{block_name}.{tsv,unphased.vcor1.gz,ldeig.npz}`) is the only supported layout.

## Acceptance criteria

- [ ] `find_blocks` returns all blocks whose `[start_bp, end_bp]` overlaps the query window; returns empty list for unknown chrom
- [ ] `match_variants` handles `chr`-prefixed SNP IDs in the TSV (strip prefix before ALID comparison)
- [ ] `match_variants` returns empty result (not an error) when fewer than 2 variants match
- [ ] `load_ld_eigenvectors` uses `.ldeig.npz` when present; falls back to loading the full LD matrix + `ld_pca` when absent
- [ ] `LDBlock.eaf` is a float64 array aligned to `snp_ids`; NaN for rows where EAF is absent
- [ ] Unit tests using a synthetic LD panel fixture (same fixture format as in `pleiodb/tests/test_impute.py`)
- [ ] No dependency on pleiodb at import time

## Blocked by

- `issues/038-ragged-imputation-kernel.md` (needs `ld_pca` for the fallback path)

## User stories addressed

- User story 13 (reuse pleiodb panel format)
- User story 16 (bound completion to blocks intersecting cis window)
- User story 10 (completion region ±1 Mb)
