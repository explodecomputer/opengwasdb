"""Public query result contract."""

from __future__ import annotations

import sqlite3
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from typing import Any

from opengwasdb.stats import beta_from_z_se, p_value_from_z


@dataclass(frozen=True)
class AssociationResult:
    """One queryable association result row."""

    variant_index: int
    alid: str
    chromosome: str
    position: int
    effect_allele: str
    other_allele: str
    rsid: str | None
    analysis_index: int
    analysis_id: str
    phenotype_id: str | None
    phenotype_label: str | None
    analysis_label: str | None
    stored_effect_scale: str
    z: float
    se: float
    beta: float
    p_value: float

    @classmethod
    def from_rows(
        cls,
        variant: Mapping[str, Any],
        analysis: sqlite3.Row,
        z: float,
        se: float,
    ) -> AssociationResult:
        return cls(
            variant_index=int(variant["variant_index"]),
            alid=str(variant["alid"]),
            chromosome=str(variant["chromosome"]),
            position=int(variant["position"]),
            effect_allele=str(variant["effect_allele"]),
            other_allele=str(variant["other_allele"]),
            rsid=variant["rsid"],
            analysis_index=int(analysis["analysis_index"]),
            analysis_id=str(analysis["analysis_id"]),
            phenotype_id=analysis["phenotype_id"],
            phenotype_label=analysis["phenotype_label"],
            analysis_label=analysis["analysis_label"],
            stored_effect_scale=str(analysis["stored_effect_scale"]),
            z=float(z),
            se=float(se),
            beta=beta_from_z_se(float(z), float(se)),
            p_value=p_value_from_z(float(z)),
        )

    def to_dict(self) -> dict[str, object]:
        return asdict(self)
