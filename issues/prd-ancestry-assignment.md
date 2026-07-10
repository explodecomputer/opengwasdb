# PRD: Ancestry Assignment & the Analysis Catalogue (ieu-a/b)

Design basis: `docs/adr/0027-analysis-catalogue-ingestion-hub.md` and
`docs/adr/0028-ancestry-assignment-from-allele-frequencies.md`.

## Problem Statement

Reference completion imputes against a single **ancestry-specific LD Reference
Panel**, so an Analysis is only correctly completed against a panel of its own
ancestry. The `ieu-a` / `ieu-b` consortium collections mix ancestries, and their
source-declared `population` metadata is coarse, often `"Mixed"`, and sometimes
wrong. To build a usable (imputable) store from these collections we must decide
each Analysis's ancestry ourselves — from its summary statistics — record it, and
build a store from a single-ancestry subset. Only a **EUR** LD panel exists today,
so the concrete near-term output is one EUR hybrid store; everything else must be
annotated and parked, not lost.

## Solution

An **annotate-then-subset** ingestion pipeline (ADR 0027):

1. **Ancestry annotator** (ADR 0028): for every `ieu-a`/`ieu-b` Analysis, extract
   allele frequencies at reference sites, estimate **Ancestry Composition** by NNLS
   mixture against a fine **Ancestry Reference Panel** (Privé 2022's UK Biobank
   21-group reference), aggregate to super-populations, and apply a multi-gate rule
   (τ, δ, overlap, residual) to
   produce an **Assigned Ancestry** or leave the Analysis **Unassigned**.
2. **Analysis Catalogue**: a versioned TSV — a superset of the build manifest —
   with one row per Analysis carrying the build columns plus the annotations
   (Assigned Ancestry, Ancestry Composition, gate results, Reported Population,
   AF-overlap).
3. **Subset & build**: filter the Catalogue to `Assigned Ancestry == EUR` → a build
   manifest → the existing `build-hybrid` + `complete-hybrid` (EUR panel). Non-EUR
   and Unassigned Analyses stay annotated in the Catalogue, parked.

Calibration uses the **Reported Population** across the collection to choose τ/δ and
to audit disagreements — never to route.

## User Stories

1. As a data engineer, I want the ancestry of each `ieu-a`/`ieu-b` Analysis inferred
   from its summary-stat allele frequencies, so I do not have to trust the coarse,
   often-`"Mixed"` declared population.
2. As a data engineer, I want a fine-reference (Privé 2022 UK Biobank 21-group) NNLS
   mixture aggregated to super-populations, so admixed/edge Analyses are modelled
   accurately but routed at the granularity our LD panels use.
3. As a data engineer, I want a multi-gate admission rule (dominant proportion,
   runner-up margin, SNP overlap, fit residual) so that only confidently
   single-ancestry Analyses are assigned and the rest are explicitly Unassigned.
4. As a developer, I want τ/δ calibrated against the declared population with a
   disagreement report, so the thresholds are evidence-based and mislabels surface.
5. As a data engineer, I want every Analysis's ancestry annotation written once into
   a versioned **Analysis Catalogue** TSV, so store membership is decided from one
   recorded source of truth.
6. As a data engineer, I want the Catalogue to be a superset of the build manifest,
   so a build is a row-filter with no format translation.
7. As a data engineer, I want to filter the Catalogue to EUR and build a homogeneous
   EUR hybrid store with the existing `build-hybrid`/`complete-hybrid`, unchanged.
8. As a developer, I want completion to be **ancestry-matched per Analysis** — impute
   only Analyses whose Assigned Ancestry matches the applied panel, recording
   `completed_against` — so a store need not be ancestry-homogeneous even though we
   build one that is.
9. As a data engineer, I want non-EUR and Unassigned Analyses parked in the Catalogue
   (annotated, not dropped), so they are re-routable when their panel exists without
   re-extracting allele frequencies.
10. As a developer, I want an `assign-ancestry` CLI that turns a raw source manifest
    into the annotated Catalogue, so ingestion is scriptable and reproducible.

## Implementation Decisions

### Ancestry Reference Panel
Privé (2022)'s UK Biobank "global reference": allele frequencies for ~5.8M variants
across 21 fine groups (bigsnpr `ref_freqs.csv.gz`, figshare file 31620968), on
GRCh37 — lifted to hg38 and re-keyed to canonical ALID — with a fine→super-population
grouping map. A static artifact (prerequisite; sourced/lifted separately). Built by
`scripts/build_ancestry_reference.py` → `ref_freqs.hg38.tsv.gz` (columns `alid`,
`chromosome`, `position`, `effect_allele`, `other_allele`, `rsid`, then one
A1-oriented frequency column per group) plus `ancestry_groups.tsv` (fine→super-pop
map). **Built artifact (this machine):**
`/local-scratch/data/opengwas/ancestry_reference/ref_freqs.hg38.tsv.gz` +
`ancestry_groups.tsv` (5,810,529 variants, 21 groups). The IEU 1000G v3 plink
reference (`fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz`) is a separate resource, useful
later for building non-EUR LD panels — not for ancestry assignment.

### Allele-frequency extraction
Targeted read of GWAS-VCF `FORMAT/AF` at the reference sites (not a full scan),
lifted hg19→hg38 and oriented to the canonical A1 (reuse `orient_to_canonical`'s
EAF flip). Exclude palindromic (A/T, C/G) variants; restrict to common
(reference MAF ≥ ~1%); require overlap ≥ N_min (~20k); no heavy LD-pruning.

### Inference & gates
NNLS mixture (α ≥ 0, Σα = 1) over shared variants → fine proportions → aggregate to
super-populations. Assigned Ancestry = dominant super-pop if proportion ≥ τ AND
margin ≥ δ AND overlap ≥ N_min AND residual ≤ gate; else Unassigned.

### Calibration
Cross-tabulate Assigned vs Reported Population over the collection; choose τ/δ at the
operating point that admits reported-European as EUR and sends reported-Mixed/other
to Unassigned; emit a disagreement report.

### Catalogue
Versioned TSV, one row per Analysis: build columns (`trait_id`, `file_path`,
`trait_name`, `n`) + annotations (`assigned_ancestry`, composition vector,
`reported_population`, gate results, `af_overlap`, `catalogue_version`,
`ancestry_reference_version`). Subset = row-filter.

### Ancestry-matched completion
Store records Assigned Ancestry per Analysis; completion imputes only matching
Analyses and records `completed_against`. Store-level Completion State stays a coarse
release flag; per-cell Association Status is ground truth.

### Code shape
`opengwasdb/ancestry/` (AF extraction, NNLS mixture, gating, calibration) + an
`assign-ancestry` CLI writing the Catalogue.

## Out of Scope

- Sourcing/lifting the Ancestry Reference Panel (Privé 2022 UKB reference,
  GRCh37→hg38) — a prerequisite artifact done separately from this pipeline
  (built by `scripts/build_ancestry_reference.py`; artifact on this machine at
  `/local-scratch/data/opengwas/ancestry_reference/`).
- Building non-EUR stores (no non-EUR LD panel yet) — non-EUR Analyses are parked.
- What ultimately happens to parked non-EUR/Unassigned Analyses (future work).
- A SQLite-backed Catalogue and all-OpenGWAS scope (TSV, ieu-a/b only for now).
- PC-projection / individual-level ancestry (we have summary statistics).
- Changing the Hybrid layout, dense/ragged builders, or completion internals beyond
  adding the per-Analysis ancestry-match filter.
