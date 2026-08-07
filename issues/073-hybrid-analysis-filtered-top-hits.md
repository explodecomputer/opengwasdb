## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

After the Dense adoption gate and Ragged implementation, deliver
analysis-filtered top hits across Hybrid stores. A selected-analysis query must
read only the relevant Dense component and overflow/Ragged component slices,
remap coordinates correctly, merge them into the shared result contract, and
return deterministic genomic ordering without first materialising either
component's global hits.

## Acceptance criteria

- [ ] Hybrid `top_hits(analysis_id=...)` delegates the selector to both component paths and does not materialise either global top-hit result before filtering.
- [ ] Dense-component variant coordinates are remapped to shared Hybrid coordinates correctly.
- [ ] Overflow associations are merged with Dense associations and returned in canonical chromosome-position order for the selected analysis.
- [ ] Global calls remain supported and are ordered by analysis index and genomic position according to the parent PRD.
- [ ] Known, unknown, zero-hit, limited, and every retained threshold query preserves the common sparse-array contract.
- [ ] Observed-only filtering retains correct behavior for the completed Dense component while overflow associations remain correctly classified.
- [ ] Validation detects component/index incompatibility and ordering or coordinate-remapping corruption.
- [ ] Integration fixtures cover hits present only in Dense, only in overflow, and in both components, including selected-analysis/global equivalence.
- [ ] A representative Hybrid benchmark demonstrates that selected-analysis lookup does not regress to global materialisation.

## Blocked by

- Blocked by `issues/071-review-top-hit-index-adoption.md` approving expansion
- Blocked by `issues/072-ragged-analysis-filtered-top-hits.md`

## User stories addressed

- User stories 1–10
- User story 32
- User story 34
