## Parent PRD

`issues/prd-v0.1-dense-observed-only.md`

## What to build

Choose and document v0.1 Dense Zarr compression and chunking defaults using the existing `besdq` benchmark prior art. The implementation should expose defaults without preventing later benchmarking or configuration.

## Acceptance criteria

- [ ] Dense writer has documented default compressor and chunk shape.
- [ ] Defaults are justified by references to existing `besdq` dense benchmark material.
- [ ] Builder metadata records compressor and chunk shape used.
- [ ] Tests verify defaults are applied and metadata records them.
- [ ] Any configurable override is documented and tested if included.

## Blocked by

- Blocked by `issues/003-build-tiny-dense-observed-only-store.md`

## User stories addressed

- User story 16
- User story 35
- User story 40

