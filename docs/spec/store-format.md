# OpenGWASDB Store Format Specification

Status: draft  
Format version described: `0.1`

This document defines the contract for valid OpenGWASDB Store Releases. The v0.1 implementation target is **Dense Observed-Only**, but the specification also records accepted design semantics for Ragged layout and Reference-Completed releases so the first implementation does not block later extensions.

Normative language:

- **MUST** indicates a requirement for a valid store.
- **SHOULD** indicates a recommended default that implementations may override with documented reason.
- **MAY** indicates an optional feature.

## 1. Store release envelope

A Store Release is a self-contained directory.

```text
manifest.json
index.sqlite
data.zarr/
variants.tsv.gz
variants.tsv.gz.tbi
variant_offsets.npy
```

`manifest.json` identifies the release and declares how to interpret it.
`index.sqlite` stores compact relational metadata such as analyses, key/value
metadata, and small alias maps. `data.zarr/` stores compressed numerical
association arrays and layout-specific numerical indexes. Dense Observed-Only
releases store their high-cardinality variant axis in `variants.tsv.gz`, indexed
by `variants.tsv.gz.tbi`, with `variant_offsets.npy` mapping Store-local Variant
Indices to BGZF row offsets.

Build and query commands operate on an explicit Store Release path. Directory naming, multi-store catalogues, default release selection, and remote API deployment are outside the store-format contract.

## 2. Manifest contract

Every Store Release MUST contain `manifest.json`.

Required fields:

| Field | Meaning |
|---|---|
| `store_id` | Stable identity of the logical Store |
| `release_id` | Identity of this immutable release |
| `format_version` | Store format version |
| `primary_layout` | `dense` or `ragged` |
| `association_coverage` | `full` or `cis_and_signals` |
| `completion_state` | `observed_only` or `reference_completed` |
| `reference_assembly` | One genome assembly for all coordinates in the release |
| `created_at` | Release creation timestamp |
| `provenance` | Source and build provenance object |

Reference-Completed releases MUST additionally declare:

| Field | Meaning |
|---|---|
| `ld_reference_panel` | The panel defining the Reference Variant Set and LD resources |
| `reference_completion_method` | Algorithm, software, version, and parameters used for imputation |

Observed-Only and Reference-Completed releases for the same source collection SHOULD share `store_id` and use different `release_id` values.

Published releases are immutable. Enhancing an Observed-Only release to Reference-Completed produces a new release.

## 3. Identity and terminology

An **Analysis** is one statistical analysis of one Trait. An Analysis produces associations between that Trait and variants.

A **Trait** is a measured or derived outcome. `phenotype` is treated as a synonym in user-facing descriptions, but the canonical term is Trait.

A Store Release contains one or more Analyses. Analysis metadata MUST be sufficient to interpret stored effect scales, sample-size semantics, and source provenance.

## 4. Variant identity and variant table

Each Store Release MUST contain a Store Variant Table.

The canonical within-Store variant key is ALID:

```text
chr:pos:A1:A2
```

where:

- coordinates are on the release's declared Reference Assembly;
- alleles are trimmed and left-aligned before identity assignment;
- `A1` is alphabetically first;
- `A2` is the other allele.

Cross-Store identity is:

```text
reference_assembly + ALID
```

`rsid` is an alias, not primary identity.

Every Store Release assigns compact Store-local Variant Indices. Variant Indices MUST NOT be assumed stable across releases or stores.

Dense Observed-Only releases use a tabix-backed Store Variant Table:

```text
variants.tsv.gz
variants.tsv.gz.tbi
variant_offsets.npy
```

The `variants.tsv.gz` table MUST contain one row per Store-local Variant Index,
sorted by chromosome, position, effect allele, other allele, and variant index.
The v0.1 dense column contract is:

```text
chromosome
position
variant_index
effect_allele
other_allele
alid
rsid
```

`variants.tsv.gz.tbi` MUST index chromosome and position for genomic range and
single-position lookup. `variant_offsets.npy` MUST contain one fixed-width
integer offset per variant row so row-index materialisation does not require a
large SQL variant table.

Long alleles MAY use deterministic hashed ALIDs as compact identifiers, but the complete normalised alleles MUST be retained once per variant so exact export and validation do not depend on an irreversible hash.

Reference-Completed releases MUST expose Reference Panel Membership as variant metadata, distinguishing Reference Variant Set variants from observed off-panel variants.

## 5. Allele orientation and signed statistics

Every association MUST be normalised to the Store's canonical allele orientation.

If the source effect allele is not canonical `A1`, the builder MUST swap orientation and negate signed statistics. Z is signed and carries effect direction. SE is non-negative.

Derived beta is:

```text
beta = z * se
```

## 6. Stored statistics

OpenGWASDB stores the canonical statistic pair:

```text
z
se
```

Requirements:

- `z` is signed.
- `se` is non-negative.
- `se` is on the same Stored Effect Scale as beta.
- beta is queryable but derived from `z * se`.
- p-value is queryable but derived from Z.
- EAF, INFO, and sample size MUST NOT be required to reconstruct beta, SE, Z, or p-value.

The v0.1 target dtype for dense statistics is `float16`, subject to validation benchmarks.

## 7. Effect scale

Each Analysis MUST declare Stored Effect Scale from the controlled vocabulary:

```text
sd_units
log_or
log_hazard
```

There is no `other`, `unknown`, or original-units Stored Effect Scale in v0.1. Unsupported stored scales MUST fail ingestion until the vocabulary is deliberately extended.

Original Effect Scale MAY be recorded as free-text provenance. For continuous traits, builders SHOULD store effects in SD Units when phenotype standard deviation is available or can be derived with acceptable provenance.

## 8. Sample size metadata

Sample size semantics are represented by kind and scope.

Sample Size Kind:

```text
participants
case_control
effective
unknown
```

Sample Size Scope:

```text
analysis
variant
none
```

OpenGWASDB MUST NOT present sample size inferred from effect statistics as an observed participant count.

Physical sample-size encoding is an implementation detail. Builders SHOULD choose the smallest lossless representation compatible with source data, including scalar counts, total N plus constant case fraction, sparse residuals, or full per-variant counts.

## 9. EAF and INFO metadata

EAF Scope and INFO Scope use the same values:

```text
absent
variant
association
```

EAF and INFO are optional. They are not required for statistical reconstruction.

Variant-scoped EAF or INFO is valid only when the builder can establish that one value is genuinely shared. Builders MUST NOT average differing association values into variant-scoped values.

For imputed associations, EAF comes from the LD Reference Panel when stored. In v0.1, EAF provenance is inferred from Association Status:

- observed association: source EAF;
- imputed association: reference-panel EAF.

## 10. Dense Observed-Only layout

Dense Observed-Only is the first implementation target.

Dense layout stores arrays with shape:

```text
n_variants x n_analyses
```

Required statistic arrays:

```text
data.zarr/z
data.zarr/se
```

In Observed-Only Dense stores:

- every finite cell represents an observed association;
- unavailable associations are represented by canonical NaN in both `z` and `se`;
- no imputed mask is required;
- the dense variant axis is source-faithful and does not require an LD Reference Panel.

Recommended compression for initial implementation is Zarr with Zstandard and bitshuffle, using benchmarked chunking appropriate for mixed range, variant, PheWAS, and full-analysis extraction workloads.

## 11. Ragged layout

Ragged layout stores Analysis-specific association sequences referencing the Store Variant Table.

Each retained association row contains at least:

```text
analysis index or analysis offset
variant_idx
z
se
```

Ragged layout is used when Analyses do not share one dense source variant axis or when Association Coverage is Cis-and-Signals.

For Observed-Only Ragged stores, absence from an Analysis sequence means the association is not retained by that Store Release.

## 12. Reference completion model

Reference Completion is an optional build phase that produces a new Reference-Completed Store Release from an Observed-Only source release.

Reference Completion exists to avoid query-time LD proxy lookup and reduce missingness by filling gaps once at build time.

Reference-Completed releases MUST:

- declare an LD Reference Panel;
- declare a Reference Completion Method;
- expose Association Status;
- produce both Z and SE for imputed associations;
- record Reference Completion Quality at LD-block-by-Analysis granularity;
- preserve observed associations not present in the Reference Variant Set.

There is no Z-only completion mode.

## 13. LD Reference Panel requirements

An LD Reference Panel used for Reference Completion MUST define:

| Requirement | Meaning |
|---|---|
| `reference_panel_id` | Stable panel identity |
| `reference_panel_version` | Panel version |
| `reference_assembly` | Genome assembly, matching the Store Release |
| `ancestry` | Population or ancestry label |
| variant list | Canonical ALIDs defining the Reference Variant Set |
| allele orientation | Orientation compatible with Store canonical ALID convention |
| LD blocks | Block definitions used by the completion method |
| LD representation | Files/data structures required by the completion method |
| checksums | Integrity checks for reference files |
| provenance | Source, build, and filtering provenance for the panel |

If EAF is emitted for imputed associations, the LD Reference Panel MUST provide EAF for panel variants or declare why EAF is absent.

The Store Release reference assembly MUST match the LD Reference Panel reference assembly.

## 14. Reference completion method requirements

The Reference Completion Method MUST record:

| Field | Meaning |
|---|---|
| `method_name` | Algorithm name |
| `method_version` | Method version |
| `software` | Software/package and version |
| `parameters` | Parameters affecting imputed Z or SE |
| `required_reference_inputs` | LD reference files/data required |
| `deterministic` | Whether same inputs produce same output |
| `quality_metric` | Definition of recorded quality values |
| `quality_thresholds` | Pass/fail or reporting thresholds, if used |

The initial intended method family is LD-eigenvector based imputation of Z and SE. The store format is not locked to one algorithm: different methods are valid if they record provenance, satisfy the `z + se` output contract, and pass validation.

The initial Reference Completion implementation should use the existing `pleiodb` imputation work as prior art:

- method source: `https://github.com/explodecomputer/pleiodb/blob/main/src/pleiodb/impute.py`;
- benchmark and reference-file location notes: `https://github.com/explodecomputer/pleiodb/blob/main/scratch/imputation_benchmark.qmd`;
- benchmarked LD panel root recorded there: `/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR`.

That prior art uses LD block TSV files, `.unphased.vcor1.gz` LD matrices, optional precomputed `.ldeig.rds` eigenfactor files, and cached `.npz` extracts. These file types are implementation inputs to the Reference Completion Method and MUST be captured through LD Reference Panel and Reference Completion Method provenance when used.

## 15. Association Status encoding

Reference-Completed Stores encode Association Status using statistic NaNs plus an imputed mask.

State derivation:

| Z | SE | imputed mask | Association Status |
|---|---|---|---|
| finite | finite | false | observed |
| finite | finite | true | imputed |
| NaN | NaN | false | missing |

Invalid states:

- only one of Z or SE is NaN;
- imputed mask is true while Z or SE is NaN;
- NaN payloads are non-canonical.

Builders and validators MUST reject invalid states.

The imputed mask is:

- dense boolean or uint8 Zarr for Dense Reference-Completed grids;
- association-aligned boolean or uint8 Zarr for Ragged Reference-Completed sequences;
- chunk-aligned with Z and SE arrays.

The imputed mask is not a sparse offsets index like the top-hit significance index.

## 16. Dense Reference-Completed layout and overflow

For Dense Reference-Completed releases:

- the dense matrix axis MUST contain only Reference Variant Set variants;
- observed off-panel associations MUST be stored in Ragged Overflow rather than discarded;
- the dense axis SHOULD be identical across Stores completed with the same LD Reference Panel.

Observed-Only Dense releases remain source-faithful and do not require an LD Reference Panel or panel-defined axis. A later Reference-Completed release MAY use a different dense axis defined by the LD Reference Panel.

Queries against Dense Reference-Completed releases with Ragged Overflow include overflow results by default for exact and range queries. Returned rows MUST expose Query Component when multiple components are involved.

## 17. Ragged Reference-Completed regions

For Ragged Cis-and-Signals molecular Stores, Reference Completion is bounded to retained regions:

- complete cis regions within existing cis boundaries;
- complete significant trans regions within existing trans-region boundaries;
- do not expand singleton suggestive associations.

Each completed region contains:

- the full slice of Reference Variant Set variants within the region boundary;
- observed off-panel variants inside the same boundary;
- NaN statistic rows for reference-panel variants that were neither observed nor imputed.

Observed, imputed, and missing rows belong to the same ragged association sequence. Ragged Reference-Completed stores do not create a separate imputed-ragged component.

Observed off-panel variants inside a Ragged Reference-Completed region have ordinary observed Association Status and are not labelled as a separate Query Component.

## 18. Top-hit query contract

Top-Hit Queries return associations ranked by statistical significance.

Dense, Ragged, Ragged Overflow, observed, and imputed associations have equal priority. Ranking is by significance, not storage component or Association Status.

Both Dense and Ragged components SHOULD provide Top-Hit Indexes using the same thresholds and result contract. The physical index encoding MAY differ by layout.

Default thresholds:

```text
5e-8
5e-6
5e-4
```

Top-hit results from Reference-Completed releases include imputed associations by default. Observed-Only Query mode applies to Top-Hit Queries as well as exact and range queries.

## 19. Query result contract

Query results SHOULD expose:

| Field | Required when |
|---|---|
| variant identity | always |
| analysis identity | always |
| z | always when available |
| se | always when available |
| beta | when requested, derived as `z * se` |
| p-value | when requested, derived from Z |
| stored_effect_scale | always or through Analysis metadata |
| association_status | Reference-Completed releases |
| query_component | multi-component releases |
| eaf | when stored/requested |
| info | when stored/requested |
| sample size | when stored/requested |

Reference-Completed queries include imputed associations by default. Query APIs MUST provide an Observed-Only mode that excludes imputed associations.

Detailed Reference Completion Quality is not joined onto every imputed association by default. Query APIs MAY expose it when requested.

## 20. Validation rules

Validators MUST check at least:

- `manifest.json` contains required fields;
- `reference_assembly` is single-valued across the release;
- canonical variant identity is valid and unique within the Store Variant Table;
- signed statistics have been normalised to canonical allele orientation;
- `se >= 0` for all finite SE values;
- Z and SE missingness is consistent;
- Stored Effect Scale values are in the controlled vocabulary;
- Dense arrays match declared dimensions;
- Reference-Completed releases declare LD Reference Panel and Reference Completion Method;
- Reference-Completed Dense axes match the Reference Variant Set;
- imputed mask is consistent with Z and SE;
- Ragged Reference-Completed regions include all Reference Variant Set variants within completed boundaries;
- top-hit indexes, when present, are consistent with stored Z values.

## 21. Compatibility

`format_version` describes compatibility of the store representation. It does not describe biological data release version or source publication version.

Readers MUST reject Store Releases with unsupported major format versions. Readers MAY support older minor versions when validation can establish compatibility.

Future format versions may add fields, arrays, or indexes, but MUST preserve explicit manifest-based feature discovery.
