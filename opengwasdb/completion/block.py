"""The per-(LD block, Analysis) Reference Completion kernel.

This is the algorithm dense and ragged completion each used to implement
separately: given one LD block's eigendecomposition and one Analysis's
observed z/se at that block's positions, impute the missing positions and
apply the region-based z-cap QC. Dense calls this once per Analysis for
every block in a vectorised sweep over its z/se matrix columns; ragged calls
it once per (block, Analysis) pair for the analyses whose cis window touches
that block. Both build the same dense NaN-filled input arrays first, so the
kernel itself does not need to know which layout is calling it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from opengwasdb.completion.impute import impute_z_block, scalar_n_se

# The completion method recorded in a completed release's manifest
# (provenance["completion"]["method"]). Both layouts now run the identical
# elastic-net-on-eigenvectors kernel below, including the region z-cap QC, so
# they share one method string rather than recording two that imply a
# difference which no longer exists.
COMPLETION_METHOD = "elastic_net_eigenvectors_v2_regioncap"

# Region-based imputation z-cap (QC): an imputed |z| may not exceed the largest
# observed |z| within +/- this many bp of it (same LD block). Caps spurious
# large imputed z-scores from LD extrapolation. See pleiodb.
REGION_CAP_BP = 1_000_000


@dataclass(frozen=True)
class BlockAnalysisResult:
    """Result of completing one LD block for one Analysis.

    ``pearson_r`` is ``None`` when imputation was never attempted (fewer than
    two observed positions in the block) and a float (possibly NaN) once it
    was -- the same three-way distinction the ``completion_quality`` table
    records. ``fills`` holds only positions actually filled, as
    ``(local_block_index, z, se)`` -- the position's index into the block's
    own SNP list, not a store-wide variant index, since that mapping is the
    caller's responsibility (dense resolves it against the union axis by
    ALID; ragged resolves it the same way). ``n_missing`` is the count of
    block positions still unfilled after this call (the full block-missing
    count when nothing could be imputed).
    """

    pearson_r: float | None
    fills: list[tuple[int, float, float]]
    n_missing: int


def _region_capped_z(
    zv: float,
    pos: int,
    obs_pos_sorted: np.ndarray,
    obs_absz_sorted: np.ndarray,
    window_bp: int = REGION_CAP_BP,
) -> float:
    """Clamp an imputed z-score so its magnitude does not exceed the largest
    observed |z| within +/- ``window_bp`` of ``pos`` (QC, per pleiodb). Returns
    ``zv`` unchanged when there is no observed variant in the window.

    ``obs_pos_sorted`` must be ascending, with ``obs_absz_sorted`` the matching
    |z| of the observed variants.
    """
    lo = int(np.searchsorted(obs_pos_sorted, pos - window_bp, side="left"))
    hi = int(np.searchsorted(obs_pos_sorted, pos + window_bp, side="right"))
    if hi <= lo:
        return zv
    z_region_max = float(obs_absz_sorted[lo:hi].max())
    if abs(zv) <= z_region_max:
        return zv
    return z_region_max if zv >= 0 else -z_region_max


def complete_block_for_analysis(
    z_dense: np.ndarray,
    se_dense: np.ndarray,
    eaf: np.ndarray,
    positions: np.ndarray,
    eigenvectors: np.ndarray,
    eigenvalues: np.ndarray,
    *,
    min_cor: float,
    region_cap_bp: int | None = REGION_CAP_BP,
) -> BlockAnalysisResult:
    """Impute one Analysis's missing positions within one LD block.

    Parameters
    ----------
    z_dense, se_dense : float64, length = n_ld_snps_in_block.
        NaN at positions not observed for this Analysis. Both must share the
        same missingness pattern (a validated store invariant).
    eaf : float64, length = n_ld_snps_in_block. Reference-panel EAF.
    positions : int64, length = n_ld_snps_in_block. Base-pair position of
        each block SNP, for the region z-cap.
    eigenvectors, eigenvalues : the block's LD eigendecomposition.
    min_cor : Pearson r quality gate (see ``impute_z_block``).
    region_cap_bp : window for the region z-cap QC, or ``None`` to disable it.
    """
    obs_mask = np.isfinite(z_dense)
    n_obs = int(obs_mask.sum())
    n_miss_block = int((~obs_mask).sum())
    if n_obs < 2:
        return BlockAnalysisResult(pearson_r=None, fills=[], n_missing=n_miss_block)

    z_imp_arr, corr = impute_z_block(z_dense, eigenvectors, eigenvalues, min_cor=min_cor)
    if z_imp_arr is None:
        pearson_r = float(corr) if np.isfinite(corr) else None
        return BlockAnalysisResult(pearson_r=pearson_r, fills=[], n_missing=n_miss_block)

    se_all = scalar_n_se(se_dense[obs_mask], eaf[obs_mask], eaf)

    obs_pos_s = obs_absz_s = None
    if region_cap_bp is not None:
        obs_idx = np.where(obs_mask)[0]
        o_order = np.argsort(positions[obs_idx])
        obs_pos_s = positions[obs_idx][o_order]
        obs_absz_s = np.abs(z_dense[obs_idx])[o_order]

    fills: list[tuple[int, float, float]] = []
    for i in np.where(~obs_mask)[0]:
        zv, sev = float(z_imp_arr[i]), float(se_all[i])
        if not (np.isfinite(zv) and np.isfinite(sev)):
            continue
        if region_cap_bp is not None:
            zv = _region_capped_z(zv, int(positions[i]), obs_pos_s, obs_absz_s, region_cap_bp)
        fills.append((int(i), zv, sev))

    return BlockAnalysisResult(
        pearson_r=float(corr), fills=fills, n_missing=n_miss_block - len(fills)
    )
