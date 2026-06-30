## Parent PRD

`issues/prd.md`

## What to build

Preserve dense top-hit query behaviour when top-hit indexes contain row indices and variant metadata lives in the tabix-backed variant axis. The query should materialise only the variants needed for the requested top-hit result set.

## Acceptance criteria

- [ ] Existing top-hit index arrays continue to store variant indices and analysis indices.
- [ ] `top_hits()` materialises required variant rows through the row-offset sidecar.
- [ ] `top_hits()` returns full public association results with variant metadata and analysis metadata.
- [ ] Missing or inconsistent top-hit rows are caught by validation.
- [ ] Top-hit tests pass without a SQLite variants table.

## Blocked by

- Blocked by `issues/018-analysis-sparse-lookup-and-row-materialisation.md`

## User stories addressed

- User story 18
- User story 28
- User story 29
- User story 30
