## Parent PRD

`issues/prd-ancestry-assignment.md`

## What to build

The thin end-to-end path for ancestry assignment, on tiny synthetic data. Given a
small **Ancestry Reference Panel** (a fixed-format table of per-population allele
frequencies over a handful of variants) and a few synthetic GWAS-VCFs carrying
`FORMAT/AF`, run the whole chain and emit one **Analysis Catalogue** TSV row per
Analysis.

The chain (PRD "Solution", ADR 0028):
- **Reference loader**: read the Ancestry Reference Panel (fine populations) and the
  fine→super-population grouping map; define the on-disk format the assigner expects.
- **AF extraction**: pull `FORMAT/AF` at the reference sites, orient to the canonical
  A1 (reuse `orient_to_canonical`'s EAF flip), drop palindromic (A/T, C/G) and
  below-MAF-floor variants.
- **Inference**: NNLS mixture (α ≥ 0, Σα = 1) over shared variants → fine
  proportions → aggregate to super-populations (**Ancestry Composition**).
- **Gates**: Assigned Ancestry = dominant super-pop if proportion ≥ τ AND margin ≥ δ
  AND overlap ≥ N_min AND residual ≤ gate; else **Unassigned**.
- **Catalogue writer**: emit a TSV that is a superset of the build manifest — build
  columns (`trait_id`, `file_path`, `trait_name`, `n`) plus `assigned_ancestry`, the
  composition vector, `reported_population`, gate results, and `af_overlap`.

Serial and small is fine; the CLI, scale, and the real reference come later.

## Acceptance criteria

- [ ] A clearly-single-ancestry synthetic study (allele frequencies matching one
      reference population) is assigned that super-population.
- [ ] An admixed synthetic study (a deliberate blend of two populations) is left
      **Unassigned** by the margin/proportion gate, not force-labelled.
- [ ] A study with too few overlapping reference sites is **Unassigned** by the
      overlap gate; a study with corrupt/mis-oriented AF is Unassigned by the
      residual gate.
- [ ] Palindromic variants are excluded from the fit; AF is oriented to the
      canonical allele so it is comparable to the reference.
- [ ] The output is a valid Catalogue TSV: a superset of the build manifest whose
      extra columns are ignored by the existing manifest reader.

## Blocked by

None - can start immediately.

## User stories addressed

- User story 1
- User story 2
- User story 3
- User story 5
- User story 6
