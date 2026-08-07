## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

Extend the observed-only Dense tracer to reference-completed Dense stores. The
same analysis-major addressing and public analysis selector must work while
preserving indexed imputation status and the existing `observed_only=True`
behavior. Completed-store queries must remain self-contained within the top-hit
index and must not fall back to millions of random reads from the dense
imputation matrix.

## Acceptance criteria

- [ ] Reference-completed Dense builds write analysis-major top-hit tiers using the same offsets and ordering contract as observed-only stores.
- [ ] Completed tiers carry imputation status aligned with every indexed association.
- [ ] A completed `top_hits(analysis_id=...)` call returns the same analysis subset as the corresponding global top-hit result.
- [ ] `observed_only=True` excludes imputed associations for both selected-analysis and global calls without reading completion status from the dense matrix.
- [ ] Known, unknown, zero-hit, and limited analysis queries retain the result and ordering semantics established by the observed-only tracer.
- [ ] Validation detects missing, misaligned, or invalid completion-status data in an analysis-major tier.
- [ ] Integration tests cover completed building, selected-analysis reads, global equivalence, observed-only filtering, limits, and corrupt completion metadata.

## Blocked by

- Blocked by `issues/068-dense-analysis-major-top-hits.md`

## User stories addressed

- User story 10
- User stories 16–17
- User stories 21–23
- User story 30
