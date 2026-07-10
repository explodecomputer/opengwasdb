## Parent PRD

`issues/prd-ancestry-assignment.md`

## What to build

Wire the tracer (`issues/062-ancestry-assignment-tracer.md`) into an
`assign-ancestry` CLI that produces the **full versioned Analysis Catalogue** for
the `ieu-a`/`ieu-b` collection. It reads a raw source manifest (all candidate
Analyses), extracts allele frequencies at the reference sites **in parallel**
(targeted `bcftools`, not a full scan), runs the mixture + gates, and writes the
Catalogue TSV with every Analysis annotated (PRD "Catalogue", "Code shape").

- Version stamps: `catalogue_version` and `ancestry_reference_version` columns/rows.
- **Parking**: non-EUR and Unassigned Analyses are present and annotated in the
  Catalogue, not dropped (PRD user story 9).
- Uses the real **Ancestry Reference Panel** (1000G+HGDP, hg38, canonical ALID) —
  an external prerequisite artifact, not built here.

## Acceptance criteria

- [ ] `opengwasdb assign-ancestry <raw-manifest> <catalogue.tsv> --ancestry-reference
      <panel> ...` writes a Catalogue with one annotated row per Analysis.
- [ ] AF extraction is parallel and targeted (reference sites only); a full genome
      scan is not performed per study.
- [ ] The Catalogue carries the version stamps and every annotation column from 062;
      non-EUR/Unassigned Analyses are retained and labelled.
- [ ] Results are independent of worker count; re-running with the same inputs and
      versions reproduces the Catalogue.

## Blocked by

- Blocked by `issues/062-ancestry-assignment-tracer.md`

## User stories addressed

- User story 1
- User story 5
- User story 9
- User story 10
