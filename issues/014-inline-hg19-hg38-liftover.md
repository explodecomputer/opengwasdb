## Parent PRD

`issues/prd-vcf-dense-build-with-liftover.md`

## What to build

A liftover utility in `opengwasdb/build/liftover.py` that converts a batch of hg19 variant positions to hg38 ALIDs in a single pass, suitable for reuse across build paths.

The core function takes an iterable of `(bare_chrom, pos, ref, alt)` tuples in hg19 coordinates and returns a `dict` mapping each input tuple to its hg38 ALID (`chr:pos:A1:A2` where A1 is alphabetically first). A single `pyliftover.LiftOver` object is created once for the batch. Variants whose position cannot be mapped are omitted from the output dict (NaN in any trait that references them). After processing, if the fraction of failures exceeds `failure_threshold` (default `0.01`), the function raises `LiftoverFailureError` reporting the count, total, and rate.

VCF positions are 1-based; pyliftover expects 0-based input; the function converts before calling `convert_coordinate` and converts back. CHROM is normalised to `chr`-prefixed form for pyliftover and stripped for the output ALID.

See the pleiodb reference implementation at `~/repo/pleiodb/src/pleiodb/liftover.py` for prior art.

## Acceptance criteria

- [ ] `build_liftover_lookup(variants, from_build, to_build, failure_threshold)` returns a `dict[(bare_chrom, pos, ref, alt) → hg38_alid]` for all successfully lifted variants.
- [ ] Variants that fail liftover are omitted from the dict; a warning is logged with the count and rate.
- [ ] When the failure rate exceeds `failure_threshold`, `LiftoverFailureError` is raised with count, total, and rate in the message.
- [ ] A single `LiftOver` object is created per call (not per variant).
- [ ] Positions are correctly converted 1-based ↔ 0-based at the pyliftover boundary.
- [ ] Both bare (`1`) and `chr`-prefixed (`chr1`) input chrom forms are accepted.
- [ ] Output ALIDs use bare chromosome and have A1 as the alphabetically first allele (allele ordering from the input tuples is preserved; no reordering happens here — the caller normalised alleles when building the input).
- [ ] Unit tests cover: successful lift of a handful of real chr1 positions, failure rate below threshold (no raise), failure rate above threshold (raises), zero-variant input (returns empty dict).
- [ ] The module imports cleanly without pyliftover installed, with the `ImportError` deferred to first call.

## Blocked by

None — can start immediately.

## User stories addressed

- User story 2 (coordinates lifted from GRCh37 to GRCh38 automatically)
- User story 3 (build fails if failure rate exceeds configurable threshold)
- User story 4 (failures below threshold are dropped and logged)
- User story 16 (liftover called once per unique variant, not once per variant per trait)
- User story 18 (threshold configurable with sensible default)
