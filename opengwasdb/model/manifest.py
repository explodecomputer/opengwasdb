"""Store manifest model."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from opengwasdb.model.enums import (
    AssociationCoverage,
    CompletionState,
    PrimaryStorageLayout,
)


@dataclass(frozen=True)
class StoreManifest:
    """Minimal manifest required to identify and open a Store Release."""

    store_id: str
    release_id: str
    format_version: str
    primary_layout: PrimaryStorageLayout
    association_coverage: AssociationCoverage
    completion_state: CompletionState
    reference_assembly: str
    created_at: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "StoreManifest":
        return cls(
            store_id=str(data["store_id"]),
            release_id=str(data["release_id"]),
            format_version=str(data["format_version"]),
            primary_layout=PrimaryStorageLayout(data["primary_layout"]),
            association_coverage=AssociationCoverage(data["association_coverage"]),
            completion_state=CompletionState(data["completion_state"]),
            reference_assembly=str(data["reference_assembly"]),
            created_at=data.get("created_at"),
            provenance=dict(data.get("provenance", {})),
        )

    @classmethod
    def load(cls, path: str | Path) -> "StoreManifest":
        manifest_path = Path(path)
        if manifest_path.is_dir():
            manifest_path = manifest_path / "manifest.json"
        with manifest_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        return cls.from_dict(data)

    def to_dict(self) -> dict[str, Any]:
        return {
            "store_id": self.store_id,
            "release_id": self.release_id,
            "format_version": self.format_version,
            "primary_layout": self.primary_layout.value,
            "association_coverage": self.association_coverage.value,
            "completion_state": self.completion_state.value,
            "reference_assembly": self.reference_assembly,
            "created_at": self.created_at,
            "provenance": self.provenance,
        }

