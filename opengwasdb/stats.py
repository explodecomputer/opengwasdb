"""Statistic helpers shared by builders, validators, and query adapters."""

from __future__ import annotations

import math

import numpy as np
from scipy import stats as _scipy_stats

_LN10 = math.log(10.0)


def beta_from_z_se(z: float, se: float) -> float:
    """Derive beta from the canonical stored statistic pair."""

    return z * se


def p_value_from_z(z: float) -> float:
    """Return the two-sided normal p-value implied by a Z score."""

    return math.erfc(abs(z) / math.sqrt(2.0))


def log10_p_two_sided(z: np.ndarray) -> np.ndarray:
    """Vectorised log10 of `p_value_from_z`, stable past the point the plain
    erfc-based p itself underflows to 0.0 in float64 (|z| >~ 38) -- e.g. the
    FADS1/FADS2 window's |z| = 47.8 (issue #104). Computed in log-space via
    `scipy.stats.norm.logsf` rather than `log(p_value_from_z(z))`, which
    would already have lost the value to underflow before the log is taken.
    """

    abs_z = np.abs(np.asarray(z, dtype="float64"))
    result: np.ndarray = (math.log(2.0) + _scipy_stats.norm.logsf(abs_z)) / _LN10
    return result


def finite_pair(z: float, se: float) -> bool:
    """True when both canonical statistics are finite and therefore queryable."""

    return math.isfinite(z) and math.isfinite(se)
