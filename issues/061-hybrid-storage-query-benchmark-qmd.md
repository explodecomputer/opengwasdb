## Parent PRD

`issues/prd-hybrid-layout.md`

## What to build

The capstone: a storage + query benchmark comparing a **Hybrid** store against a
**Dense-of-union** store built from the *same* real, heterogeneous collection,
rendered as a QMD and added to `docs/index.html`. Placed last so it exercises a
**complete, CLI-built, reference-completed, validated, top-hit-indexed** hybrid
store — an accurate representation, not a stub (PRD user story 13).

**Collection** (small, real, heterogeneous — different consortia/panels, hg19):
- `ieu-a-2` — BMI (GIANT)
- `ieu-a-7` — CAD (CARDIoGRAMplusC4D)
- `ieu-a-300` — LDL cholesterol

Under `/local-scratch/data/opengwas/igd/`. These arrive as raw consortium text on
GRCh37, so the build exercises per-study normalisation + hg19→hg38 liftover across
heterogeneous formats (a realistic hybrid input); a small per-study source
manifest / column mapping is expected data-prep.

**What to measure** (Hybrid vs Dense-of-union, both built and completed against the
same reference panel via the CLIs):
- Storage size of each store (and the dense/overflow split for the hybrid), plus
  the on-panel vs off-panel variant counts.
- Query latency across `analysis`, `phewas`, `range_phewas`, `lookup`,
  `top_hits`.
- **Positive-control MR, imputed-only** — following the ukb-b benchmark's MR
  section (`docs/benchmark-output/opengwasdb_ukbb_dense_benchmark.qmd`): run
  **LDL → CAD** and **BMI → CAD** IVW MR restricted to **imputed instruments
  only** (both exposure and outcome cells imputed by completion of the Dense
  Component), and compare each against the observed-only estimate. Both are
  expected positive; agreement between the imputed-only and observed estimates
  validates that the completed hybrid store's imputed variants reproduce the
  causal signal.
- **Regional imputation plot** — mirror that same qmd's "Regional imputation
  sanity check": a ~1 Mb window centred on a strong imputed instrument, plotting
  every association returned for the window with **imputed associations
  highlighted vs observed** (via `association_status`), so the imputed fill is
  visible against the observed backbone.

Reuse the existing benchmark harness + QMD/`index.html` pattern
(`benchmarks/`, `docs/benchmark-output/`), including the imputed-only MR and
regional-imputation-plot code already in the ukb-b benchmark.

## Acceptance criteria

- [ ] A Hybrid store and a Dense-of-union store are built from `ieu-a-2/7/300`
      against the same reference panel (via `build-hybrid` and the dense CLI),
      reference-completed, and both validate.
- [ ] The QMD reports storage size (with the hybrid dense/overflow split and
      on-/off-panel counts) and per-query latency for both stores, with a table +
      chart.
- [ ] The QMD runs **LDL→CAD** and **BMI→CAD** IVW MR using **imputed-only
      instruments** and compares each to the observed-only estimate; all are in the
      expected (positive) direction and the imputed-only vs observed estimates
      agree.
- [ ] The QMD includes a **regional plot** of a ~1 Mb window with imputed
      associations highlighted against observed (matching the ukb-b benchmark's
      regional imputation sanity check).
- [ ] The rendered HTML is added as a card + source link in `docs/index.html`.
- [ ] The benchmark is reproducible from a committed script under `benchmarks/`.

## Blocked by

- Blocked by `issues/056-hybrid-unified-query-surface.md`
- Blocked by `issues/057-hybrid-top-hit-index.md`
- Blocked by `issues/058-hybrid-dense-component-completion.md`
- Blocked by `issues/059-hybrid-store-validation.md`
- Blocked by `issues/060-build-hybrid-cli.md`

## User stories addressed

- User story 13
- (demonstrates end-to-end: 5, 6, 7, 8)
