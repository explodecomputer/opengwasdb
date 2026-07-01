## Parent PRD

`issues/prd-ragged-observed-only.md`

## What to build

A builder that reads BESD files (`.esi`, `.epi`, `.besd`) and produces a valid Ragged Observed-Only Store.

Port the BESD reader from besdq (`besd_reader.py`) into `opengwasdb/layouts/ragged/besd_reader.py`, then implement `build_ragged_from_besd()` in `opengwasdb/layouts/ragged/build_besd.py`.

Build sequence:

1. Read ESI → construct the Store Variant Table → write `variants.tsv.gz` + tabix index + `variant_alid_bytes.npy`
2. Read EPI → construct `TraitRecord` list → write `traits.tsv.gz` + tabix index
3. For each probe (in EPI order): read BESD associations → map SNP row indices to `variant_index` → accumulate in `RaggedCSRWriter`
4. Flush `RaggedCSRWriter` → write `data.zarr/ragged/`
5. Write `manifest.json` with `primary_layout: "ragged"`, `completion_state: "observed-only"`, format version, store_id, release_id, creation time
6. Write `index.sqlite` with `analyses` table (probe_id, gene_id, gene_name, tissue, analysis_index, probe_chr, probe_bp, n)

Statistics: store z and SE directly. z = beta / se from BESD; SE is the source-reported value. Both stored as float16. No AF reconstruction, no scalar-N approximation.

The BESD datasets in `/local-scratch/data/hg38/` are already GRCh38; no liftover is applied.

Add a CLI command `build-ragged-besd` to `opengwasdb/cli/main.py`:

```
opengwasdb build-ragged-besd <besd_prefix> <output_path> --store-id <id> --release-id <id>
```

where `besd_prefix` is the path without extension (e.g. `/local-scratch/data/hg38/eqtlgen/cis-eQTL`).

## Acceptance criteria

- [ ] `build_ragged_from_besd(besd_prefix, output_path, store_id, release_id)` produces a complete store directory
- [ ] `manifest.json` contains `primary_layout: "ragged"` and `completion_state: "observed-only"`
- [ ] `variants.tsv.gz` contains all SNPs from the ESI file, tabix-indexed
- [ ] `traits.tsv.gz` contains one row per probe from the EPI file, tabix-indexed by probe_chr/probe_bp
- [ ] `data.zarr/ragged/` contains valid CSR arrays with correct offsets
- [ ] `opengwasdb validate <output_path>` exits 0
- [ ] `build-ragged-besd` CLI command is reachable via `opengwasdb --help`
- [ ] Integration test: build from a synthetic 3-probe BESD fixture, assert variant count, analysis count, and round-trip z/se values
- [ ] eqtlgen smoke test: build from eqtlgen BESD, assert n_analyses ≈ expected gene count, store validates

## Blocked by

- `issues/034-traits-tsv-writer-and-reader.md`
- `issues/035-ragged-zarr-csr.md`

## User stories addressed

- User story 5
- User story 7
- User story 9
- User story 10
