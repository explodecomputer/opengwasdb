## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

Deliver the observed-only Dense tracer for analysis-addressable top hits. A
Dense build must produce the analysis-major representation described in the
parent PRD, and the public top-hit query must accept an optional analysis ID and
read only that analysis's contiguous indexed slice. Preserve calls without an
analysis selector, the sparse-array result contract, existing threshold tiers,
and deterministic analysis/genomic ordering.

This slice includes the narrow reader/writer abstraction, all Dense build paths,
query-facade behavior, format validation, and external-behavior tests needed to
make the observed-only path independently usable and verifiable.

## Acceptance criteria

- [ ] Every existing Dense top-hit threshold tier is written in analysis-index order, with canonical genomic ordering within each analysis.
- [ ] Each tier contains a valid offsets structure of length analysis-count plus one, including valid empty slices for analyses with no hits.
- [ ] Inline-harvested and full-matrix-rebuilt indexes produce equivalent public results and deterministic ordering.
- [ ] `top_hits(analysis_id=...)` resolves a known public analysis ID and reads only its indexed slice rather than materialising and filtering the global index.
- [ ] An unknown analysis ID returns the standard empty sparse-array result, consistent with other analysis-addressed queries.
- [ ] `top_hits()` without an analysis selector remains supported and returns results ordered by analysis index and then canonical genomic position.
- [ ] `limit` is applied after scope selection and returns the first entries in the documented analysis/genomic order.
- [ ] All existing threshold and sparse-array result fields retain their current semantics.
- [ ] Dense validation rejects missing, wrong-length, non-monotonic, and terminally inconsistent offsets, incorrect analysis slices, and non-genomic slice ordering.
- [ ] Round-trip and Dense integration tests cover analyses with many, one, and zero hits; known and unknown analysis IDs; global reads; limits; and all threshold tiers.
- [ ] The top-hit reader and writer expose a small testable interface that keeps physical addressing details out of builders and the public query facade.

## Blocked by

None - can start immediately.

## User stories addressed

- User stories 1–9
- User stories 11–16
- User stories 18–24
