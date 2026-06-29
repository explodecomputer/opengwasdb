## Parent PRD

`issues/prd.md`

## What to build

Update the public local build/query workflow so the dense observed-only vertical slice produces and queries the new tabix-backed variant-axis envelope from the CLI. The existing tiny fixture should remain self-contained, and the dense vertical-slice benchmark should use the EBI full GWAS-SSF inputs only.

## Acceptance criteria

- [ ] The CLI build command writes the new dense variant-axis files.
- [ ] The CLI query commands continue to support canonical ALID, rsid, range, analysis, lookup, and top-hit workflows.
- [ ] The standard validation command accepts the new dense store envelope.
- [ ] The tiny vertical-slice tests exercise the CLI path with the new envelope.
- [ ] Dense vertical-slice docs and tests do not use the sparse `38714679` dataset as a benchmark or comparison target.
- [ ] Dense Zarr compression, dtype, and chunking defaults are unchanged.

## Blocked by

- Blocked by `issues/017-alias-lookup-semantics.md`
- Blocked by `issues/018-analysis-sparse-lookup-and-row-materialisation.md`
- Blocked by `issues/019-top-hits-with-variant-offsets.md`

## User stories addressed

- User story 9
- User story 11
- User story 15
- User story 22
- User story 23
- User story 24
- User story 25
- User story 26
- User story 27
- User story 28
