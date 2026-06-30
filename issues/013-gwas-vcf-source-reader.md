## Parent PRD

`issues/prd-vcf-dense-build-with-liftover.md`

## What to build

A layout-independent GWAS-VCF reader in `opengwasdb/build/vcf_source.py`. The module streams records from a bgzipped GWAS-VCF file (GWAS-SSF VCF spec: Lyon et al. 2021) using cyvcf2 and exposes two functions:

- One that yields `(chrom, pos, ref, alt)` tuples for every biallelic record — used by pass 1 of the dense build to collect the union variant set.
- One that yields per-record association data `(chrom, pos, ref, alt, z, se, stored_effect_scale)` — used by pass 2 to fill zarr columns.

Both functions handle bare (`1`) and `chr`-prefixed (`chr1`) CHROM values by normalising to bare form. Multi-allelic records (ALT containing `,`) are skipped silently. `stored_effect_scale` is inferred once per file from the `StudyType` field in the `##SAMPLE` VCF meta-line (`CaseControl → LOG_ODDS`, `Continuous → SD_UNITS`); the function raises `ValueError` with the file path if the field is absent or unrecognised.

Z-score is taken from `FORMAT/EZ` when present and finite; otherwise computed as `FORMAT/ES / FORMAT/SE`. When the effect allele is not the canonical A1 (alphabetically first), z is negated. SE is always taken from `FORMAT/SE` and is never sign-flipped. Records with SE ≤ 0 or non-finite z are skipped.

This module has no knowledge of liftover, zarr, or the dense layout. It is a pure streaming reader. See the pleiodb reference implementation at `~/repo/pleiodb/src/pleiodb/vcf.py` for prior art on EZ/ES/SE extraction and cyvcf2 usage.

## Acceptance criteria

- [ ] `stream_vcf_variants(path)` yields `(bare_chrom: str, pos: int, ref: str, alt: str)` for every biallelic record in a GWAS-VCF; multi-allelic records are silently skipped.
- [ ] `stream_vcf_associations(path)` yields `(bare_chrom, pos, ref, alt, z, se, stored_effect_scale)` with z derived preferentially from EZ, falling back to ES/SE.
- [ ] A record where the effect allele is alphabetically second (canonical A1 = REF > ALT) has its z negated; SE is unchanged.
- [ ] Both functions handle bare and `chr`-prefixed CHROM forms, normalising to bare.
- [ ] `StudyType=CaseControl` → `StoredEffectScale.LOG_ODDS`; `StudyType=Continuous` → `StoredEffectScale.SD_UNITS`.
- [ ] Missing or unrecognised `StudyType` raises `ValueError` naming the file.
- [ ] Records with SE ≤ 0 or non-finite z are skipped without raising.
- [ ] Unit tests pass using synthetic GWAS-VCF fixtures written to `tmp_path` (no real data files required).
- [ ] The module imports cleanly without cyvcf2 installed, with the `ImportError` deferred to first call.

## Blocked by

None — can start immediately.

## User stories addressed

- User story 5 (StudyType → stored_effect_scale inference)
- User story 6 (multi-allelic records skipped)
- User story 11 (allele orientation normalised to ALID convention)
- User story 13 (bare and chr-prefixed CHROM handled)
- User story 14 (EZ preferred over ES/SE)
- User story 15 (layout-independent reader reusable by future build paths)
