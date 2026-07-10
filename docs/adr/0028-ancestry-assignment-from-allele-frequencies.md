# Ancestry assignment from allele frequencies; ancestry-matched completion

Reference completion imputes against a single **ancestry-specific LD Reference
Panel**, so an Analysis can only be completed against a panel of its own ancestry
— imputing a European GWAS against an East-Asian LD structure is wrong. The
`ieu-a` / `ieu-b` consortium collections are a mixture of ancestries, and their
source-declared `population` metadata is coarse, often literally `"Mixed"`, and
sometimes wrong. We need to decide each Analysis's ancestry *ourselves*, from its
summary statistics, and route it accordingly.

## Decision

**Recover each Analysis's ancestry from its summary-statistic allele frequencies**
(GWAS-VCF `FORMAT/AF`, present throughout this harmonised collection), and let that
govern routing. The method follows Privé (2022): model an Analysis's allele
frequencies as a **non-negative, sum-to-one mixture** of reference-population
frequencies over a common variant set, solved by NNLS, yielding **Ancestry
Composition**. We **fit against a fine reference** (1000G+HGDP, ~21 groups —
the **Ancestry Reference Panel**) for accuracy on admixed and edge populations,
then **aggregate the proportions up to 1000G super-populations** for the routing
label, because our LD panels are keyed to super-populations. Fit-fine,
label-coarse.

An Analysis's **Assigned Ancestry** is the dominant super-population, admitted only
if it clears a **multi-gate rule**: proportion ≥ τ, margin over the runner-up ≥ δ,
reference-SNP overlap ≥ N_min, and NNLS fit residual below a gate. Failing any gate
leaves the Analysis **Unassigned** rather than force-labelled. τ/δ are **calibrated**
against the **Reported Population** across the collection — which is used only to
choose the operating point and to audit disagreements, **never to route**.

**Ancestry-Matched Completion:** a store records **Assigned Ancestry per Analysis**,
and reference completion imputes an Analysis only when its Assigned Ancestry matches
the panel being applied; other Analyses are left observed-only (`completed_against =
null`). Store-level Completion State stays a coarse release flag ("a completion pass
ran"); per-cell Association Status remains the ground truth for what was imputed.
This makes ancestry-homogeneity **optional at the store contract level** — but for
the first `ieu-a`/`ieu-b` build we deliberately build a **homogeneous EUR store**
from a EUR subset of the Catalogue, because mixing ancestries in one store inflates
missingness (non-overlapping off-panel tails) and only EUR has an LD panel today.

## Considered options

- **Trust the source-declared `population`.** Rejected: coarse, frequently `"Mixed"`,
  sometimes wrong — it is the thing AF-matching exists to improve on. Kept only as a
  calibration/audit signal.
- **PC-projection ancestry (individual-level).** Not applicable: we have summary
  statistics, not genotypes. The AF-mixture is the summary-statistic analogue and is
  what Privé's `snp_ancestry_summary` implements.
- **Fit directly against 5 super-population means (coarse-only).** Rejected:
  super-population mean frequencies model admixed samples poorly (AMR is itself
  admixed), so the EUR-vs-rest gate is miscalibrated. Fitting fine and aggregating
  fixes this at the cost of a larger static reference table.
- **A single dominant-proportion threshold.** Rejected in favour of the multi-gate
  rule: overlap and residual gates catch corrupt / mis-oriented AF that a bare
  proportion threshold would force-label; the margin gate catches near-ties.
- **Store-level ancestry as a homogeneity requirement.** Rejected: a per-Analysis
  attribute + Ancestry-Matched Completion is strictly more general (supports future
  multi-panel completion) and costs nothing here, where we still choose a homogeneous
  build.

## Consequences

- Palindromic (A/T, C/G) variants are excluded from the fit (strand-ambiguous,
  unalignable); the fit uses common (reference MAF ≥ ~1%) variants in the
  reference ∩ study-AF intersection, without heavy LD-pruning (NNLS tolerates
  correlated SNPs).
- The **Ancestry Reference Panel** (1000G+HGDP frequencies, on hg38 with canonical
  ALIDs, plus a fine→super-population map) is a new static artifact to source and
  normalise.
- Ancestry assignment is an **annotator that writes into the Analysis Catalogue**
  (ADR 0027): Assigned Ancestry, Ancestry Composition, gate results, and the
  Reported-Population comparison. Stores inherit Assigned Ancestry via the subset
  they are built from.
- Unassigned and non-target-ancestry Analyses are annotated and **parked in the
  Catalogue**, not dropped — re-routable when their panel exists, without
  re-extracting allele frequencies.
- Completion becomes ancestry-aware: it reads per-Analysis Assigned Ancestry and
  imputes only matching Analyses, recording `completed_against` per Analysis.
