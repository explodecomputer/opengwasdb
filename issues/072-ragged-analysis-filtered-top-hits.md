## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

After the Dense adoption gate approves expansion, deliver equivalent
analysis-filtered top-hit behavior for Ragged stores. Choose a
Ragged-appropriate physical representation while preserving the public
analysis-ID selector, threshold tiers, empty-result behavior, sparse-array
contract, analysis/genomic ordering, limit semantics, and validation guarantees
proven by the Dense experiment.

This slice is complete only when a Ragged store can be built, validated, and
queried through the public facade without scanning the global top-hit index for
a selected analysis.

## Acceptance criteria

- [ ] The adoption review explicitly approves expansion before implementation starts.
- [ ] Ragged `top_hits(analysis_id=...)` reads an analysis-addressable index path rather than materialising and filtering every global hit.
- [ ] Known, unknown, zero-hit, global, limited, and every retained threshold query follows the public semantics established by the parent PRD.
- [ ] Results are ordered by analysis index globally and canonical genomic position within an analysis.
- [ ] The chosen Ragged representation is documented and does not assume that the Dense physical layout is automatically appropriate.
- [ ] Ragged build and rebuild paths produce equivalent selected-analysis and global results.
- [ ] Validation detects malformed analysis addressing and incorrect analysis/genomic ordering.
- [ ] Public integration tests cover selected-analysis/global equivalence and preserve all existing Ragged top-hit behavior.
- [ ] A representative Ragged benchmark demonstrates the selected-analysis improvement and reports any global-query or build regressions.

## Blocked by

- Blocked by `issues/071-review-top-hit-index-adoption.md` approving expansion

## User stories addressed

- User stories 1–10
- User story 32
- User story 34
