## Parent PRD

`issues/prd-analysis-major-top-hits.md`

## What to build

Tune and benchmark the completed Dense experiment on the full UK Biobank
observed-only and reference-completed stores. Measure the representative
single-analysis query against the previous global-index materialisation and
filtering path, evaluate chunk sizes appropriate to narrow analysis slices, and
render the evidence into the existing Dense benchmark report.

The benchmark is the experiment's deliverable: it must make the adoption gate
decidable rather than merely report the fastest favorable number.

## Acceptance criteria

- [ ] The benchmark queries the same representative analysis and significance tier on observed-only and reference-completed full UK Biobank Dense stores.
- [ ] Warm-cache repetition policy, median, tail latency, result count, and selected chunk configuration are recorded.
- [ ] The previous global-index-read-and-filter latency is retained as a direct comparison and the achieved speedup is calculated explicitly.
- [ ] First-read or cold behavior is reported separately where practical and is not conflated with the warm-cache acceptance target.
- [ ] Multiple sensible chunk sizes are evaluated for narrow analysis slices, with the selected trade-off justified by measured selected-analysis and global-read behavior.
- [ ] Index storage size and build or rebuild cost are reported even though space is not the primary constraint.
- [ ] Correctness checks demonstrate that benchmarked selected-analysis results equal the corresponding global subsets for observed and completed stores.
- [ ] The report states whether each store meets the under-10-ms warm-cache target and whether the improvement over the 300–400-ms baseline is substantial.
- [ ] Benchmark source data, JSON results, explanatory source, and rendered HTML are updated together and do not make unsupported low-millisecond claims.

## Blocked by

- Blocked by `issues/068-dense-analysis-major-top-hits.md`
- Blocked by `issues/069-completed-analysis-major-top-hits.md`

## User stories addressed

- User stories 24–30
- User story 33
