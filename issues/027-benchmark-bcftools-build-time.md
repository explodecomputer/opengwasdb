## Parent PRD

`issues/prd-bcftools-vcf-reader.md`

## What to build

Rebuild the chr1 × 100-analysis store from scratch using the bcftools-based reader and record build time. Run the same five query benchmarks to confirm query performance is unchanged (queries do not use the VCF reader, so no regression expected). Write results to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_bcftools_benchmark.json`.

Update ADR-0019 to record the outcome: if bcftools is materially faster, amend the decision to "bcftools is the chosen reader; cyvcf2 is no longer used". If the improvement is small, document the measured times and note that the bottleneck lies elsewhere (likely liftover or numpy allocation, not VCF parsing).

A Quarto report is written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_bcftools_benchmark.qmd` using the same `write_qmd()` pattern as the existing benchmark scripts. The report includes build time, store size, query timings, and an interpretation section comparing bcftools vs cyvcf2 build times and opengwasdb vs besdq query timings.

## Acceptance criteria

- [ ] Store rebuilds successfully using `build_dense_from_vcf_manifest` with the bcftools reader.
- [ ] `validate_store` passes on the rebuilt store.
- [ ] Build time recorded in `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_bcftools_benchmark.json` alongside the cyvcf2 baseline (876.6s).
- [ ] Quarto report written to `docs/benchmark-output/opengwasdb_vcf_ukb_chr1_bcftools_benchmark.qmd`.
- [ ] QMD renders without error (`quarto render <path>`).
- [ ] Five query timings included in the output JSON (confirming no query regression).
- [ ] ADR-0019 updated with the measured build time comparison and a clear statement of the new decision.

## Blocked by

- `issues/026-bcftools-vcf-reader.md`

## User stories addressed

- User story 9 from `issues/prd-bcftools-vcf-reader.md`
