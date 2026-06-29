"""Statistic helpers shared by builders, validators, and query adapters."""

from __future__ import annotations

import math


def beta_from_z_se(z: float, se: float) -> float:
    """Derive beta from the canonical stored statistic pair."""

    return z * se


def p_value_from_z(z: float) -> float:
    """Return the two-sided normal p-value implied by a Z score."""

    return math.erfc(abs(z) / math.sqrt(2.0))


def finite_pair(z: float, se: float) -> bool:
    """True when both canonical statistics are finite and therefore queryable."""

    return math.isfinite(z) and math.isfinite(se)
