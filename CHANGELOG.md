# Changelog

All notable changes to the `opengwasdb` package.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning is described in [`CONTRIBUTING.md`](CONTRIBUTING.md#versioning) —
in particular, **the package version and a Store Release's `format_version`
are different things** and move independently. See the compatibility table at
the end of this file.

## [Unreleased]

Work lands on `dev` and appears here under *Unreleased* until `dev` merges to
`main`, at which point it is cut into a version.

### Changed

- `CONTRIBUTING.md` gains a Documentation section: what goes stale, what
  regenerates it, that benchmark numbers are re-run rather than edited, that
  `opengwasdb-stores` holds docs depending on this package's CLI surface, and
  a pre-merge checklist.

### Fixed

- `benchmarks/README.md` documented seven commands as `uv run`, which stopped
  working when the project moved to Pixi.
- `README.md` described a two-month-old codebase as "newly scaffolded".

## [0.2.0] — 2026-08-22

The first tagged version. `0.1.0` was the placeholder the project carried from
its first commit and was never released, so this entry backfills the changes a
user would notice across the whole history to date.

### Added

- **Effect allele frequency is stored** (ADR 0036, #106). Every source format
  reported it and every reader parsed it; no layout kept it. Now a per
  (variant, Analysis) plane in Dense, Ragged and Hybrid, oriented to the stored
  effect allele, declared per Analysis by a new `eaf_scope` column, and exposed
  through the query facade and `--variant-info`.
- **rsid lookups work** (#109). An rsid search index is written by
  `write_variant_axis`, so every layout gets one. Collision policy is explicit:
  an rsid names every row it appears on.
- **Human-readable query output by default** (#104), with `--variant-info` for
  rsid and eaf, and `--format json` for the raw index-keyed result.
- **`overview.html`** — a generated, store-wide browsable rendering of
  `analyses.tsv` with persisted Top-Hit Counts (ADR 0032).
- **Hybrid layout** — a Dense Component over a reference panel plus a Ragged
  Overflow Component for off-panel observations (ADR 0026).
- **Reference Completion** for Dense, Ragged and Hybrid, with ancestry-matched
  imputation (ADR 0028) and per-block quality records.
- **Ancestry assignment from allele frequencies** (ADR 0029) and the Analysis
  Catalogue ingestion hub (ADR 0027).
- **Rho matrix** for Dense stores (ADR 0025).
- Validation: a closed-envelope rule (#80), rsid index coverage, and EAF range
  and shape checks.

### Changed

- **`analyses.tsv` is the single Analytical Metadata contract for every layout**
  (ADR 0030, ADR 0034). `phenotype_id`/`phenotype_label`/`trait_id` retired in
  favour of `analysis_label`/`trait_ontology_id`/`trait_ontology_label`;
  `gene_id`/`gene_name` retired too (ADR 0035). Ragged's divergent SQLite
  `analyses` table is gone — a leftover one is now a validation failure.
- **Manifest metadata reaches the store.** Dense and Hybrid in #86, Ragged in
  #83 — fifteen shared-core columns that manifests supplied and builders
  silently dropped, now carried through one shared `PassthroughMetadata`.
- The `SourceReader` interface carries variant identifiers (`SourceVariant`)
  and per-association EAF, both of which sources provided and the interface
  previously had nowhere to put.
- Query results are a documented adapter contract across layouts (ADR 0033).

### Fixed

- **Hybrid `analyses_table()` undercounted Top-Hit Counts** by delegating to the
  Dense Component (#107). Hid 4,476 hits on the gwas-catalog pilot; one Analysis
  reported 10,727 instead of 14,706.
- **Ragged and Hybrid Reference Completion silently dropped data** on rebuild —
  rsids (#109) and observed EAF (found in review of #106). The eqtlgen pilot had
  49,967 rsids observed and none in its completed sibling.
- **Dense and Hybrid never captured rsids at all** (#109) — 0 in 50,000 rows of
  both the finngen and gwas-catalog pilots.
- Disjoint-partition violation on Hybrid LD-panel extension (#99).
- Duplicate canonical-variant rows within a Ragged Analysis (#101); the
  `trait_id` requirement dropped from the Ragged SSF builder (#100).
- Reference completion now warns when no Analysis carries `assigned_ancestry`
  and it would otherwise impute everything against one panel (#108).

### Known limitations

- `z` and `se` are `float16`. At the significant tail this costs real accuracy:
  worst-case p error is 26% at |z| = 30, and the 568 cells above |z| = 64 in a
  UK Biobank-scale store have no meaningful precision. Tracked as #114; see
  ADR 0037.
- Reference-Completed releases carry no EAF on imputed cells (#113, #116).
- Existing Store Releases predate every build-time fix above and must be
  rebuilt to gain them (#117).

## Store format compatibility

The package version and the store `format_version` are independent. A store
records the `format_version` it was written against; the package records which
it can read.

| package | writes `format_version` | reads |
|---|---|---|
| 0.2.0 | 0.1 | 0.1 |

[Unreleased]: https://github.com/opengwas/opengwasdb/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/opengwas/opengwasdb/releases/tag/v0.2.0
