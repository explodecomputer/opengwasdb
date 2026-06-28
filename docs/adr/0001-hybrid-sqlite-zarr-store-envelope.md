# Use a hybrid SQLite and Zarr store envelope

Each Store Release uses a small root manifest, an embedded SQLite database for relational metadata and lookup structures, and a Zarr hierarchy for numerical association arrays. This keeps OpenGWASDB stores standalone while using SQLite for metadata/index strengths and Zarr for compressed chunked numerical payloads.

The standard release layout is:

```text
manifest.json
index.sqlite
data.zarr/
```

