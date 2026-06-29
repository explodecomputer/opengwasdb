"""Build pipeline orchestration."""

from typing import Any

from opengwasdb.build.source import (
    NormalisedAssociation,
    SourceRowError,
    read_normalised_associations,
)


def build_dense_observed_from_sources(*args: Any, **kwargs: Any) -> Any:
    from opengwasdb.build.observed import build_dense_observed_from_sources as _build

    return _build(*args, **kwargs)

__all__ = [
    "NormalisedAssociation",
    "SourceRowError",
    "build_dense_observed_from_sources",
    "read_normalised_associations",
]
