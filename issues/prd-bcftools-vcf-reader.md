## Problem Statement

The current GWAS-VCF reader (`vcf_source.py`) uses cyvcf2 to iterate records in Python. The chr1 × 100-analysis build took 876 seconds — roughly 14 minutes — because cyvcf2 processes each VCF record as a Python object with per-field extraction overhead. The besdq prototype used `bcftools query` subprocess piping, which runs in C with multi-threaded decompression. The build time gap is the build-pipeline counterpart to the query-side object-materialisation problem.

## Solution

Replace the cyvcf2 record-iteration internals of `stream_vcf_variants` and `stream_vcf_associations` with `bcftools query` subprocesses piped through stdout. The public function signatures are unchanged so `build_dense_from_vcf_manifest` requires no modification. `StudyType` is still parsed from the VCF header via `bcftools view -h`. cyvcf2 is removed as a runtime dependency.

## User Stories

1. As a developer building a Dense Store from GWAS-VCFs, I want pass 1 (variant collection) to stream records via bcftools, so that the union variant set is collected faster than with cyvcf2.
2. As a developer building a Dense Store from GWAS-VCFs, I want pass 2 (association fill) to stream records via bcftools, so that the zarr matrix fill completes in substantially less than 876 seconds.
3. As a developer, I want `StudyType` inference to work with the bcftools-based reader using the same VCF header parsing logic as before, so that `stored_effect_scale` is still correctly assigned.
4. As a developer, I want the build to raise a clear error if `bcftools` is not found in PATH, so that misconfigured environments fail fast with an actionable message rather than a confusing subprocess error.
5. As a developer, I want multi-allelic records to still be silently skipped in the bcftools path, so that the ALID contract is preserved.
6. As a developer, I want EZ to still be preferred over ES/SE when present in the bcftools output, so that z-score fidelity is maintained.
7. As a developer, I want the allele-flip logic (negate z when ALT > REF) to still be applied in the bcftools reader, so that ALID orientation is correct.
8. As a developer, I want the existing unit tests for `stream_vcf_variants` and `stream_vcf_associations` to pass unchanged against the bcftools implementation (same synthetic VCF fixtures), so that correctness is not regressed.
9. As a developer, I want the chr1 × 100-analysis build time to be measured and recorded in `docs/benchmark-output/` after the bcftools change, so that the improvement (or lack thereof) is documented.

## Implementation Decisions

- `stream_vcf_variants` uses `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n"` via `subprocess.Popen`, iterating stdout lines.
- `stream_vcf_associations` uses `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%EZ\t%ES\t%SE]\n"`, parsing six fields per line; missing values (`.`) are treated as NaN.
- `read_vcf_study_type` uses `bcftools view -h` piped into the existing `_STUDY_TYPE_RE` regex — no change to the regex or the `StoredEffectScale` mapping.
- If `bcftools` is not on PATH, all three functions raise `RuntimeError` with a clear install message. The check is done at first call, not at import time.
- cyvcf2 and its lazy-import guard (`_get_cyvcf2`) are removed from the module.
- ADR-0019 is updated to record that cyvcf2 was replaced by bcftools after benchmarking confirmed the performance gap.
- No change to `build_dense_from_vcf_manifest`, `build_liftover_lookup`, or any query code.

## Testing Decisions

- The existing test suite in `test_vcf_source.py` uses synthetic plain-text VCF fixtures written to `tmp_path`. bcftools can read uncompressed VCFs so the fixtures require no change.
- The only module under test is `vcf_source.py`; all 14 existing tests should pass without modification other than replacing the cyvcf2 dependency.
- One new test: `bcftools` unavailable → `RuntimeError` with helpful message.
- Prior art: the existing `test_vcf_source.py` fixture pattern.

## Out of Scope

- Query array API (separate PRD).
- bcftools regional filtering (tabix `-R` pre-filter for pass 2). This is a further optimisation that could be added later.
- Parallelising pass 1 or pass 2 across multiple VCFs.
- Any change to `build_liftover_lookup` or `build_dense_from_vcf_manifest`.

## Further Notes

- Reference implementation: `besdq/scripts/dense_03_build_stats_arrays.py` uses `subprocess.Popen(['bcftools', 'query', ...])` and iterates stdout lines — the direct prior art.
- ADR-0019 currently says "VCF parsing uses cyvcf2; performance against bcftools query is to be benchmarked and this choice may be revisited." This PRD is that revisit.
