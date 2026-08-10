"""Tests for opengwasdb.build.phenotype_sd (issue #15)."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from opengwasdb.build.phenotype_sd import estimate_phenotype_sd
from opengwasdb.model.enums import OriginalSdMethod


def _synthetic_se(af: np.ndarray, sd_true: float, sample_size: float) -> np.ndarray:
    """se implied by a per-variant asymptotic SE model for a near-null,
    common-variant continuous trait: se = sd / sqrt(2*af*(1-af)*N) -- the
    inverse of the ADR-0029 estimator formula, so recovering sd_true from
    these se/af pairs exercises the real inversion, not a tautology."""
    return sd_true / np.sqrt(2.0 * af * (1.0 - af) * sample_size)


def test_recovers_known_sd_from_source_maf():
    rng = np.random.default_rng(1)
    af = rng.uniform(0.05, 0.5, 500)
    sample_size = 20_000.0
    sd_true = 2.5
    se = _synthetic_se(af, sd_true, sample_size)

    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, sample_size, se=se, af=af
    )

    assert result.method is OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF
    assert result.sd == pytest.approx(sd_true, rel=1e-6)
    assert result.dispersion == pytest.approx(0.0, abs=1e-6)


def test_reports_reference_maf_method_tier():
    rng = np.random.default_rng(2)
    af = rng.uniform(0.05, 0.5, 200)
    sample_size = 10_000.0
    se = _synthetic_se(af, 1.8, sample_size)

    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_REFERENCE_MAF, sample_size, se=se, af=af
    )

    assert result.method is OriginalSdMethod.ESTIMATED_FROM_REFERENCE_MAF
    assert result.sd == pytest.approx(1.8, rel=1e-6)


def test_reports_beta_distribution_method_tier():
    rng = np.random.default_rng(3)
    beta = rng.normal(loc=0.0, scale=1.4, size=1000)

    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_BETA_DISTRIBUTION, 10_000.0, beta=beta
    )

    assert result.method is OriginalSdMethod.ESTIMATED_FROM_BETA_DISTRIBUTION
    assert result.sd == pytest.approx(1.4, rel=0.1)


def test_dispersion_rises_with_corrupted_allele_frequency():
    """A clean f -> 1-f flip cannot move this diagnostic: 2*f*(1-f) is
    symmetric under allele exchange, so that is not a useful corruption to
    test with. A genuinely wrong-ancestry AF substitutes an unrelated
    frequency for the same variant -- that does perturb the per-variant
    implied SD, which is what the dispersion diagnostic exists to catch."""
    rng = np.random.default_rng(4)
    n = 500
    sample_size = 20_000.0
    sd_true = 2.0
    af_true = rng.uniform(0.05, 0.5, n)
    se = _synthetic_se(af_true, sd_true, sample_size)

    clean = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, sample_size, se=se, af=af_true
    )

    # A narrow, systematically-low af range (rather than another draw from the
    # same distribution the true af came from) makes the corruption's effect
    # unambiguous: real "wrong ancestry" AF is systematically biased, not
    # merely re-randomised from the same population.
    af_corrupted = af_true.copy()
    corrupt_idx = rng.choice(n, size=n // 2, replace=False)
    af_corrupted[corrupt_idx] = rng.uniform(0.01, 0.03, len(corrupt_idx))
    corrupted = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, sample_size, se=se, af=af_corrupted
    )

    assert clean.dispersion == pytest.approx(0.0, abs=1e-6)
    # Clearly non-zero and "measurable", not just above the noise floor.
    assert corrupted.dispersion > 0.05


@pytest.mark.parametrize("sample_size", [None, 0.0, -5.0, float("nan")])
def test_refuses_missing_or_unusable_sample_size(sample_size):
    se = np.array([0.1, 0.2, 0.3])
    af = np.array([0.2, 0.3, 0.4])

    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, sample_size, se=se, af=af
    )

    assert result.method is OriginalSdMethod.UNAVAILABLE
    assert np.isnan(result.sd)
    assert result.notes


def test_refuses_when_no_variant_has_usable_se_or_af():
    se = np.array([np.nan, np.nan])
    af = np.array([0.2, 0.3])

    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, 10_000.0, se=se, af=af
    )

    assert result.method is OriginalSdMethod.UNAVAILABLE


def test_refuses_beta_distribution_with_fewer_than_two_finite_values():
    result = estimate_phenotype_sd(
        OriginalSdMethod.ESTIMATED_FROM_BETA_DISTRIBUTION,
        10_000.0,
        beta=np.array([1.0, np.nan]),
    )

    assert result.method is OriginalSdMethod.UNAVAILABLE


def test_rejects_a_non_estimation_method():
    with pytest.raises(ValueError, match="estimated_from_source_maf"):
        estimate_phenotype_sd(
            OriginalSdMethod.DECLARED_STANDARDISED,
            10_000.0,
            se=np.array([0.1]),
            af=np.array([0.2]),
        )


def test_rejects_af_method_missing_se_or_af():
    with pytest.raises(ValueError, match="requires se and af"):
        estimate_phenotype_sd(OriginalSdMethod.ESTIMATED_FROM_SOURCE_MAF, 10_000.0)


def test_rejects_beta_method_missing_beta():
    with pytest.raises(ValueError, match="requires beta"):
        estimate_phenotype_sd(OriginalSdMethod.ESTIMATED_FROM_BETA_DISTRIBUTION, 10_000.0)


def test_scalar_n_se_shares_implementation_with_phenotype_sd():
    """opengwasdb.completion.impute.scalar_n_se must call this module's
    se_scale_constant rather than recompute the same median independently --
    the issue #15 consolidation requirement. Patching the name as imported
    into the impute module (not the original) proves scalar_n_se actually
    calls through, rather than merely producing a coincidentally identical
    result via its own copy of the formula."""
    from opengwasdb.build.phenotype_sd import se_scale_constant
    from opengwasdb.completion.impute import scalar_n_se

    se_obs = np.array([0.1, 0.2, 0.15])
    eaf_obs = np.array([0.2, 0.3, 0.4])
    eaf_ref = np.array([0.25, 0.35])

    with patch(
        "opengwasdb.completion.impute.se_scale_constant", wraps=se_scale_constant
    ) as spy:
        scalar_n_se(se_obs, eaf_obs, eaf_ref)

    spy.assert_called_once()
    called_se, called_af = spy.call_args[0]
    np.testing.assert_array_equal(called_se, se_obs)
    np.testing.assert_array_equal(called_af, eaf_obs)
