"""NNLS ancestry mixture and the multi-gate admission rule (ADR 0028).

Models a study's A1-oriented allele frequencies as a non-negative, sum-to-one
mixture of the fine reference groups (solved by NNLS with a soft sum-to-one
penalty row), aggregates the fitted proportions to super-populations, and admits
a single **Assigned Ancestry** only when the dominant super-population clears the
proportion (τ), margin (δ), overlap (N_min), and NNLS-residual gates. Any failure
leaves the Analysis **Unassigned** rather than force-labelled.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import nnls  # type: ignore[import-untyped]

from opengwasdb.ancestry.extract import extract_af_at_sites
from opengwasdb.ancestry.reference import AncestryReference


@dataclass(frozen=True)
class Gates:
    """Multi-gate admission thresholds (ADR 0028)."""

    tau: float = 0.90  # min dominant super-population proportion
    delta: float = 0.20  # min margin of dominant over runner-up
    n_min: int = 20_000  # min overlapping reference sites
    residual_max: float = 0.06  # max RMS NNLS residual (mis-oriented/corrupt AF)
    sum_to_one_penalty: float = 10.0  # weight of the soft Σα = 1 constraint row


@dataclass(frozen=True)
class AncestryAssignment:
    """Result of fitting one Analysis against the Ancestry Reference Panel."""

    assigned_ancestry: str | None  # super-population, or None if Unassigned
    dominant_superpop: str | None
    dominant_proportion: float
    runner_up_margin: float
    af_overlap: int
    residual: float
    gate_reason: str  # "ok" | "overlap" | "residual" | "proportion" | "margin"
    superpop_composition: dict[str, float] = field(default_factory=dict)
    fine_composition: dict[str, float] = field(default_factory=dict)


def assign_ancestry(
    study_af: dict[str, float],
    reference: AncestryReference,
    gates: Gates | None = None,
) -> AncestryAssignment:
    """Fit ``study_af`` against the reference and apply the multi-gate rule."""
    gates = gates or Gates()

    rows = [reference.index[a] for a in study_af if a in reference.index]
    overlap = len(rows)
    if overlap == 0:
        return AncestryAssignment(
            assigned_ancestry=None,
            dominant_superpop=None,
            dominant_proportion=0.0,
            runner_up_margin=0.0,
            af_overlap=0,
            residual=float("nan"),
            gate_reason="overlap",
        )

    row_idx = np.asarray(rows, dtype=np.int64)
    ref_alids = reference.alids[row_idx]
    b = np.asarray([study_af[a] for a in ref_alids.tolist()], dtype=np.float64)
    A = reference.freqs[row_idx]  # (overlap, G)

    # Drop any reference NaNs (uninformative for this fit).
    valid = ~np.isnan(A).any(axis=1) & ~np.isnan(b)
    A, b = A[valid], b[valid]
    overlap = int(A.shape[0])
    if overlap == 0:
        return AncestryAssignment(
            assigned_ancestry=None,
            dominant_superpop=None,
            dominant_proportion=0.0,
            runner_up_margin=0.0,
            af_overlap=0,
            residual=float("nan"),
            gate_reason="overlap",
        )

    alpha = _fit_mixture(A, b, gates.sum_to_one_penalty)
    residual = float(np.sqrt(np.mean((A @ alpha - b) ** 2)))

    fine_composition = dict(zip(reference.groups, alpha.tolist(), strict=True))
    sp_props = reference.aggregate(alpha)
    superpop_composition = dict(zip(reference.superpops, sp_props.tolist(), strict=True))

    order = np.argsort(sp_props)[::-1]
    dominant_superpop = reference.superpops[int(order[0])]
    dominant_proportion = float(sp_props[order[0]])
    runner_up = float(sp_props[order[1]]) if sp_props.size > 1 else 0.0
    margin = dominant_proportion - runner_up

    gate_reason = apply_gates(
        gates,
        overlap=overlap,
        residual=residual,
        dominant_proportion=dominant_proportion,
        margin=margin,
    )
    assigned = dominant_superpop if gate_reason == "ok" else None

    return AncestryAssignment(
        assigned_ancestry=assigned,
        dominant_superpop=dominant_superpop,
        dominant_proportion=dominant_proportion,
        runner_up_margin=margin,
        af_overlap=overlap,
        residual=residual,
        gate_reason=gate_reason,
        superpop_composition=superpop_composition,
        fine_composition=fine_composition,
    )


def assign_from_vcf(
    vcf_path: str | Path,
    reference: AncestryReference,
    gates: Gates | None = None,
    *,
    regions_file: str | Path | None = None,
) -> AncestryAssignment:
    """Extract AF at reference sites from a GWAS-VCF, then assign ancestry."""
    study_af = extract_af_at_sites(vcf_path, reference.index.keys(), regions_file=regions_file)
    return assign_ancestry(study_af, reference, gates)


def _fit_mixture(A: np.ndarray, b: np.ndarray, penalty: float) -> np.ndarray:
    """NNLS with a soft Σα = 1 row, then renormalise to sum exactly to one."""
    n_groups = A.shape[1]
    if penalty > 0:
        A_aug = np.vstack([A, np.full((1, n_groups), penalty)])
        b_aug = np.concatenate([b, [penalty]])
    else:
        A_aug, b_aug = A, b
    alpha_raw, _ = nnls(A_aug, b_aug)
    alpha = np.asarray(alpha_raw, dtype=np.float64)
    total = alpha.sum()
    if total > 0:
        alpha = alpha / total
    return alpha


def apply_gates(
    gates: Gates,
    *,
    overlap: int,
    residual: float,
    dominant_proportion: float,
    margin: float,
) -> str:
    """Return the first failing gate, or ``"ok"``. Order matters (ADR 0028).

    Pure function of the summary statistics, so a calibrated τ/δ pick can relabel
    a Catalogue from its stored numbers without re-extracting allele frequencies.
    """
    if overlap < gates.n_min:
        return "overlap"
    if not np.isfinite(residual) or residual > gates.residual_max:
        return "residual"
    if dominant_proportion < gates.tau:
        return "proportion"
    if margin < gates.delta:
        return "margin"
    return "ok"
