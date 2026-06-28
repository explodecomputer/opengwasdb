## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Create the first end-to-end path for opening an explicit Store Release directory, reading its manifest, reporting basic information, and validating the minimal release envelope. This implements the path described in the PRD Solution items 1, 8, and 9.

## Acceptance criteria

- [ ] A valid minimal Store Release directory with `manifest.json` can be opened from an explicit path.
- [ ] `info` output reports store ID, release ID, format version, primary layout, association coverage, completion state, and reference assembly.
- [ ] `validate` succeeds for a valid minimal release and returns useful errors for missing or malformed manifest fields.
- [ ] The implementation does not perform discovery, default release selection, or catalogue lookup.
- [ ] Tests cover successful open/info/validate and common manifest failures.

## Blocked by

None - can start immediately.

## User stories addressed

- User story 2
- User story 3
- User story 4
- User story 18
- User story 19
- User story 29
- User story 33
- User story 35

