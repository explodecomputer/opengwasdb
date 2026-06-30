"""Open local Store Releases."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from opengwasdb.model.manifest import StoreManifest

if TYPE_CHECKING:
    from opengwasdb.query.facade import StoreQuery


@dataclass(frozen=True)
class OpenGWASDBStore:
    """A local Store Release opened from an explicit path."""

    path: Path
    manifest: StoreManifest

    @property
    def index_path(self) -> Path:
        return self.path / "index.sqlite"

    @property
    def data_path(self) -> Path:
        return self.path / "data.zarr"

    @property
    def variant_table_path(self) -> Path:
        return self.path / "variants.tsv.gz"

    @property
    def variant_tabix_path(self) -> Path:
        return self.path / "variants.tsv.gz.tbi"

    @property
    def variant_offsets_path(self) -> Path:
        return self.path / "variant_offsets.npy"

    def query(self) -> StoreQuery:
        from opengwasdb.query.facade import StoreQuery

        return StoreQuery(self)


def open_store(path: str | Path) -> OpenGWASDBStore:
    """Open a local Store Release directory.

    The function intentionally opens exactly the path supplied by the caller.
    Release discovery and default selection belong to a higher-level catalogue.
    """

    store_path = Path(path)
    manifest = StoreManifest.load(store_path)
    return OpenGWASDBStore(path=store_path, manifest=manifest)
