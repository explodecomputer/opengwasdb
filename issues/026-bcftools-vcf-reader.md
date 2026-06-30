## Parent PRD

`issues/prd-bcftools-vcf-reader.md`

## What to build

Rewrite the internals of `opengwasdb/build/vcf_source.py` to use `bcftools` subprocesses instead of cyvcf2 for all three public functions: `stream_vcf_variants`, `stream_vcf_associations`, and `read_vcf_study_type`. Public signatures are unchanged.

**`stream_vcf_variants(path)`**
Uses `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" <path>` via `subprocess.Popen`, iterates stdout lines, skips multi-allelic records (ALT containing `,`), and normalises CHROM to bare form.

**`stream_vcf_associations(path)`**
Uses `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%EZ\t%ES\t%SE]\n" <path>`. Parses six tab-separated fields per line. Missing values (`.`) become `float('nan')`. Applies the same EZ-preferred-over-ES/SE logic, SE ≤ 0 skip, and allele-flip (negate z when ALT > REF). Infers `stored_effect_scale` by calling `read_vcf_study_type` once per file before the loop.

**`read_vcf_study_type(path)`**
Uses `subprocess.run(['bcftools', 'view', '-h', path], capture_output=True, text=True)` and applies the existing `_STUDY_TYPE_RE` regex to the stdout. No change to the regex or the `StoredEffectScale` mapping.

**bcftools availability check**
All three functions check for `bcftools` via `shutil.which('bcftools')` on first call and raise `RuntimeError("bcftools not found in PATH — install via conda: conda install -c bioconda bcftools")` if absent. The check is done at call time, not at import time.

Remove the `_get_cyvcf2()` lazy-import helper and all cyvcf2 references from the module.

The existing 14 tests in `test_vcf_source.py` use synthetic plain-text VCF fixtures written to `tmp_path`. bcftools reads plain VCF files so no fixture changes are needed. Add one new test: when bcftools is not on PATH, all three functions raise `RuntimeError` with a message containing "bcftools".

## Acceptance criteria

- [ ] All 14 existing tests in `test_vcf_source.py` pass against the bcftools implementation.
- [ ] New test: mock `shutil.which` to return `None`; all three functions raise `RuntimeError` mentioning "bcftools".
- [ ] No import of `cyvcf2` anywhere in `vcf_source.py`.
- [ ] `_get_cyvcf2()` helper is removed.
- [ ] `stream_vcf_associations` still negates z when ALT > REF (allele flip test passes).
- [ ] `stream_vcf_associations` still prefers EZ over ES/SE (EZ test passes).
- [ ] `stream_vcf_associations` still skips records with SE ≤ 0.
- [ ] Multi-allelic records (comma in ALT) are skipped in both streaming functions.
- [ ] The full test suite (`pytest tests/`) passes.

## Blocked by

None — can start immediately (independent of issues 023–025).

## User stories addressed

- User stories 1–8 from `issues/prd-bcftools-vcf-reader.md`
