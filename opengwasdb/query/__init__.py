"""Layout-independent query API."""

from opengwasdb.query.facade import StoreQuery, query_store
from opengwasdb.query.results import AssociationResult

__all__ = ["AssociationResult", "StoreQuery", "query_store"]
