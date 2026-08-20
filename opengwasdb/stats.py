"""Statistic helpers shared by builders, validators, and query adapters."""

from __future__ import annotations

import math

import numpy as np
from scipy import stats as _scipy_stats

_LN10 = math.log(10.0)


_MISSING_AF = {"", ".", "NA", "NaN", "nan", "None"}


def parse_af(value: str | float | None) -> float | None:
    """A usable allele frequency, or None (ADR 0036).

    One rule for every source format: unparseable, non-finite, and
    out-of-range values are all None rather than clamped or substituted -- a
    frequency outside [0, 1] is a broken row, not a nearly-right one, and a
    fabricated 0.5 would be indistinguishable from a real one downstream.

    Lives here rather than beside any one reader because the GWAS-VCF, tabular
    (GWAS-SSF/FinnGen) and Ragged-SSF paths each need it and each sits in a
    different corner of the import graph; three copies of the same rule is
    three chances for one to drift.
    """
    if value is None:
        return None
    if isinstance(value, str):
        if value.strip() in _MISSING_AF:
            return None
        try:
            af = float(value)
        except ValueError:
            return None
    else:
        af = float(value)
    return af if math.isfinite(af) and 0.0 <= af <= 1.0 else None


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
