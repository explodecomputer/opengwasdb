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

## Repository status

This repository is newly scaffolded. The design baseline lives in:

- [CONTEXT.md](./CONTEXT.md)
- [docs/spec/store-format.md](./docs/spec/store-format.md)
- [docs/adr/](./docs/adr/)

## Development

```bash
pip install -e ".[dev]"
pytest
```
