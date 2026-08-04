## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

Review the Dense experiment with a human decision maker and record the adoption
decision. Use the observed-only and completed UK Biobank evidence to decide
whether the improvement is substantial enough to expand, whether the
analysis-major representation should become part of the required store format,
and whether existing releases need migration or rebuild tooling.

This is a HITL decision slice. It does not presume that meeting one latency
number automatically justifies adoption; correctness, build cost, global-query
effects, chunking, validation complexity, and operational consequences must be
considered together.

## Acceptance criteria

- [ ] A human reviewer confirms whether the Dense experiment demonstrates a substantial and reproducible speedup.
- [ ] The review records whether the under-10-ms warm-cache target was met for both observed-only and completed stores.
- [ ] The decision explicitly approves or declines expansion to Ragged and Hybrid.
- [ ] The decision records whether analysis-major indexing becomes mandatory, remains optional/experimental, or is abandoned.
- [ ] The decision records whether existing stores require an in-place rebuild path, release rebuilds only, or no migration work.
- [ ] Any accepted global-ordering, chunking, storage, build-time, or compatibility trade-offs are documented in the project's architectural documentation or a linked follow-up decision artifact.
- [ ] If expansion is declined, the benchmark evidence and reason are retained and the dependent layout issues are closed or marked not approved without implementation.

## Blocked by

- Blocked by `issues/070-benchmark-analysis-major-top-hits.md`

## User stories addressed

- User stories 31–34

## Decision

Approved on 2026-08-04 after review of the full UK Biobank Dense experiment.
Both observed-only and reference-completed stores met the under-10-ms target
(1.17 ms and 1.77 ms warm-cache medians), with 261x and 225x improvements over
global materialisation and filtering. Expansion to Ragged and Hybrid is
approved.

Analysis-major top-hit tiers and their validation contract are adopted for new
Dense, Ragged, and Hybrid builds. Existing stores can be migrated in place with
the layout-specific top-hit rebuild commands; full release rebuilds are not
required. The accepted trade-off is a 16,384-entry association chunk, genomic
ordering within analyses, analysis-index ordering globally, and increased
index storage in exchange for substantially faster selected-analysis access.
