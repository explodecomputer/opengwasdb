# opengwasdb

Standalone storage and query engine for OpenGWAS-scale summary statistic stores.

The project is starting from a clean store contract:

- self-contained store releases;
- embedded SQLite for metadata and lookup indexes;
- Zarr for compressed association arrays;
- layout-independent build and query APIs;
- Dense and Ragged primary layouts;
- optional reference completion using LD reference panels.

The first implementation slice is intentionally narrow: **Dense Observed-Only** stores with `z` and `se` arrays, metadata, validation, and layout-independent queries. Ragged layout, reference completion, and service/catalogue deployment are recorded in the ADRs but are not part of v0.1.

For the broader OpenGWAS platform direction, see [docs/opengwas-roadmap.md](./docs/opengwas-roadmap.md).

## Repository status

This repository is newly scaffolded. The design baseline lives in:

- [CONTEXT.md](./CONTEXT.md)
- [docs/spec/store-format.md](./docs/spec/store-format.md)
- [docs/adr/](./docs/adr/)

## Development

All Python dependencies and native tooling (bcftools) are managed by
[Pixi](https://pixi.sh) from `pyproject.toml`'s `[tool.pixi.*]` tables. A
fresh checkout needs only Pixi installed:

```bash
pixi run -e dev test        # pytest
pixi run -e dev lint        # ruff check .
pixi run -e dev typecheck   # mypy opengwasdb
```

For a one-off script or REPL, use `pixi run -e dev python <script>.py`
rather than invoking a bare interpreter.

The package is still a normal `pip install`-able library for downstream
consumers (`pyproject.toml` + hatchling): `pip install -e ".[dev]"` continues
to work if you'd rather manage the environment yourself, but it won't provide
`bcftools` — install that separately in that case.

