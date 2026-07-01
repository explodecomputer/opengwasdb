## Problem Statement

Molecular QTL datasets (eQTL, methylation QTL, pQTL) cannot be stored in the Dense Layout because each analysis (one gene-tissue or CpG-tissue pair) has associations for only a small cis window — typically 1,000–10,000 variants out of ~10M. A Dense store would be 99.9% NaN, wasting storage and query time.

A Ragged Layout is needed: each analysis stores only its retained associations, with the variant axis shared across the store but each analysis maintaining its own sequence of (variant_index, z, se) tuples.

## Solution

Implement a Ragged Observed-Only Store: a new Primary Storage Layout backed by zarr Compressed Sparse Row (CSR) arrays, a tabix-indexed `traits.tsv.gz` for analysis metadata, and the same `variants.tsv.gz` / `index.sqlite` envelope as the Dense Layout. Build the first production store from the eqtlgen BESD dataset.

## User Stories

1. As a bioinformatician, I want to retrieve all associations for a specific gene-tissue pair so that I can run colocalisation against a GWAS locus.
2. As a bioinformatician, I want to query all associations in a genomic region (by probe location) so that I can find which genes have cis evidence overlapping a GWAS signal.
3. As a bioinformatician, I want to query all associations in a genomic region (by variant location) so that I can extract the full cis window around a locus.
4. As a developer, I want the ragged store to expose the same query facade as the dense store so that downstream tools need no layout-specific branching.
5. As a data engineer, I want to build a ragged store from BESD files so that existing eQTL resources (GTEx, eqtlgen, BrainMeta, godmc) can be ingested without format conversion.
6. As a data engineer, I want a `traits.tsv.gz` file I can inspect with standard tools so that I understand what analyses a store contains without opening SQLite.
7. As a data engineer, I want the ragged store to record z and SE directly (not reconstructed from AF and N) so that statistics are lossless and reproducible.
8. As a data engineer, I want the traits.tsv schema to accommodate GWAS traits (no probe position), eQTL (gene TSS), methylation QTL (CpG position), and pQTL (protein/gene) so that one format covers the full OpenGWAS data catalogue.
9. As a data engineer, I want a CLI command `build-ragged-besd` so that ragged stores can be built from scripts without writing Python.
10. As a data engineer, I want the store to validate cleanly with `opengwasdb validate` so that build correctness is verifiable.
11. As a developer, I want the ragged layout to live under `opengwasdb/layouts/ragged/` so that it is structurally parallel to `layouts/dense/`.
12. As a developer, I want `traits.tsv.gz` to be tabix-indexed by probe_chr/probe_bp (when coordinates are available) so that regional-by-probe queries are O(log n).

## Implementation Decisions

### Primary Storage Layout — zarr CSR (analysis-centric)

Associations are stored as four parallel zarr arrays sorted by analysis:

```
data.zarr/ragged/
  offsets        [n_analyses + 1]  int64   — exclusive end index per analysis
  variant_index  [n_total_assoc]   int32   — store-local variant row index
  z              [n_total_assoc]   float16
  se             [n_total_assoc]   float16
```

`offsets[i]` to `offsets[i+1]` gives the slice for analysis `i`. Access is O(1): read two int64 values, then a zarr slice. No deserialization overhead — direct memory-mapped numpy views.

Variant-centric (phewas-style) access is explicitly deferred. The phewas pattern is less critical for molecular QTL data (most probes have no data for a given variant) and can be added later by duplicating the data in a `by_variant/` CSR group.

### Statistics

z and SE are stored directly as float16. No AF or scalar-N reconstruction (unlike besdq ScalarN mode). SE is always the `se_vector` equivalent — lossless per-association values from the source.

### traits.tsv.gz schema

Tabix-indexed on `probe_chr` / `probe_bp` when coordinates are present. Column order matches tabix convention (chr/pos first):

```
#probe_chr  probe_bp  analysis_index  analysis_id  probe_id  n  gene_id  gene_name  tissue  context
```

- `probe_chr` / `probe_bp`: NA for GWAS traits (no tabix index created); TSS for eQTL; CpG position for mQTL
- `analysis_index`: store-local integer (links to zarr CSR offsets)
- `analysis_id`: unique string within store (constructed at build time, e.g. `ENSG00000000003::Whole_Blood`)
- `probe_id`: source-native identifier (Ensembl gene ID, CpG ID, protein ID, GWAS trait ID)
- `n`: sample size
- `gene_id` / `gene_name`: populated for eQTL and pQTL; NA for mQTL and GWAS
- `tissue`: populated for all molecular QTL; NA for GWAS
- `context`: optional (sex-stratified, ancestry label, etc.)

The file is bgzipped and tabix-indexed using the same tooling as `variants.tsv.gz`. For GWAS dense stores, the file exists but `probe_chr`/`probe_bp` are NA and no tabix index is created.

### Store envelope

A Ragged Store Release uses the same top-level envelope as Dense:

```
manifest.json          — primary_layout: "ragged", format_version, store_id, etc.
index.sqlite           — analyses table (probe_id, gene_id, tissue, analysis_index) + variants alias table
variants.tsv.gz        — store-wide variant table (same format as Dense)
variants.tsv.gz.tbi
variant_alid_bytes.npy — mmap'd ALID sidecar for O(log n) lookups
traits.tsv.gz          — analysis metadata (new; tabix-indexed)
traits.tsv.gz.tbi      — (present only when probe coordinates available)
data.zarr/ragged/      — CSR arrays
```

### Build pipeline — BESD input

Port the BESD reader from besdq (`besd_reader.py`): `IndexReader.read_esi()`, `IndexReader.read_epi()`, `BESDReader`. Build proceeds in one pass per probe, writing CSR arrays incrementally via zarr append or pre-allocated chunks. The store variant table is built from ESI before any probe data is written. BESD datasets in `/local-scratch/data/hg38/` are already in GRCh38; no liftover is applied.

### Query facade

Ragged stores expose the same query interface as Dense: `query_region`, `query_analysis`, `query_tophits`. Under the hood, `query_region` by probe uses tabix on `traits.tsv.gz`; `query_region` by variant uses tabix on `variants.tsv.gz` then CSR slice lookup. `query_analysis` uses the SQLite `analyses` table to resolve probe_id → analysis_index, then a single CSR slice.

### Modules

- `opengwasdb/layouts/ragged/` — new package: `build_besd.py`, `query.py`, `zarr_csr.py`
- `opengwasdb/traits/axis.py` — new module: traits.tsv.gz writer and tabix-indexed reader (parallel to `variants/axis.py`)
- `opengwasdb/cli/main.py` — new `build-ragged-besd` command
- `manifest.json` schema — add `primary_layout` field (already in CONTEXT.md; enforce in code)

### First target dataset

eqtlgen cis-eQTL BESD from `/local-scratch/data/hg38/eqtlgen/` (or equivalent prefix). This is the integration test that validates the full pipeline end-to-end.

### Deferred

- Variant-centric (phewas) CSR duplication
- VCF input format for ragged stores
- Reference-Completed Ragged stores (imputed cis regions)
- Ragged Overflow from Dense stores

## Testing Decisions

- Unit tests for the zarr CSR writer and reader: write N analyses with known associations, read back by index, assert exact values.
- Unit tests for `traits.tsv.gz` writer and tabix reader: write a small traits file, query by probe_id and by region, assert correct rows returned.
- Integration test: build a ragged store from a small BESD fixture (synthetic or a slice of eqtlgen), run `validate`, assert store structure is correct.
- Prior art: `tests/test_dense_vertical_slice.py` and `tests/test_dense_vcf_build.py` show the pattern for build + validate + query round-trips.

## Out of Scope

- Variant-centric (phewas) access path
- VCF input format for ragged builds
- Reference completion for ragged stores
- Ragged Overflow from Dense stores
- Multi-tissue store merging
- Trait-level significance filtering or LD clumping at build time

## Further Notes

- besdq `builder.py` and `sqlite_query.py` are the reference implementation for BESD parsing and the probe-centric query pattern.
- besdq stores z-scores only with AF/N reconstruction; OpenGWASDB stores z + SE directly — this is a deliberate departure from the besdq model.
- `issues/prd-v0.1-dense-observed-only.md` is the prior art for the Dense Layout PRD structure.
- The `traits.tsv.gz` format is intentionally designed to accommodate all molecular QTL probe types (eQTL, mQTL, pQTL) and GWAS traits in a single schema.
