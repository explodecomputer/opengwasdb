"""Unit tests for the ragged imputation kernel (issue 038)."""

from __future__ import annotations

import numpy as np
import pytest

from opengwasdb.layouts.ragged.impute import (
    elastic_net_impute,
    impute_z_block,
    ld_pca,
    poly_rescale,
    scalar_n_se,
)


def _random_pd_ld(n: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((n, n))
    return A @ A.T + np.eye(n) * n * 0.1


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def small_ld(rng):
    return _random_pd_ld(20, rng)


class TestLdPca:
    def test_eigenvalues_descending(self, small_ld):
        vals, _ = ld_pca(small_ld, thresh=0.9)
        assert np.all(vals[:-1] >= vals[1:])

    def test_cumvar_reaches_thresh(self, small_ld):
        thresh = 0.85
        vals, _ = ld_pca(small_ld, thresh=thresh)
        all_vals, _ = ld_pca(small_ld, thresh=1.0)
        cumvar = np.cumsum(vals) / all_vals.sum()
        assert cumvar[-1] >= thresh

    def test_fewer_components_than_variants(self, small_ld):
        vals, _ = ld_pca(small_ld, thresh=0.5)
        assert len(vals) < small_ld.shape[0]

    def test_eigenvector_shape(self, small_ld):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        assert vecs.shape == (small_ld.shape[0], len(vals))

    def test_non_negative_eigenvalues(self, small_ld):
        vals, _ = ld_pca(small_ld, thresh=0.9)
        assert np.all(vals >= 0)


class TestElasticNetImpute:
    def test_returns_full_length(self, small_ld, rng):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64)
        z[::4] = np.nan
        result = elastic_net_impute(z, vecs, len(vals))
        assert result is not None
        assert len(result) == small_ld.shape[0]

    def test_no_nans_in_result(self, small_ld, rng):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64)
        z[:5] = np.nan
        result = elastic_net_impute(z, vecs, len(vals))
        assert result is not None
        assert np.all(np.isfinite(result))

    def test_returns_none_all_missing(self, small_ld):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = np.full(small_ld.shape[0], np.nan)
        assert elastic_net_impute(z, vecs, len(vals)) is None

    def test_returns_none_constant_z(self, small_ld):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = np.ones(small_ld.shape[0])
        z[:3] = np.nan
        assert elastic_net_impute(z, vecs, len(vals)) is None

    def test_returns_none_single_observed(self, small_ld):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = np.full(small_ld.shape[0], np.nan)
        z[0] = 1.5
        assert elastic_net_impute(z, vecs, len(vals)) is None


class TestPolyRescale:
    def test_high_correlation(self, rng):
        n = 100
        pred = rng.standard_normal(n)
        truth = 2.5 * pred + rng.standard_normal(n) * 0.1
        truth[10:20] = np.nan
        _, corr = poly_rescale(truth, pred, npoly=1)
        assert np.isfinite(corr) and corr > 0.9

    def test_output_length_matches_input(self, rng):
        n = 50
        pred = rng.standard_normal(n)
        truth = pred * 1.5
        truth[:5] = np.nan
        adj, _ = poly_rescale(truth, pred)
        assert len(adj) == n

    def test_degenerate_too_few_points(self):
        truth = np.array([1.0, np.nan, np.nan, np.nan])
        pred = np.array([0.5, 1.0, 1.5, 2.0])
        adj, corr = poly_rescale(truth, pred, npoly=3)
        assert len(adj) == 4

    def test_all_missing_truth_returns_nan_corr(self, rng):
        truth = np.full(20, np.nan)
        pred = rng.standard_normal(20)
        _, corr = poly_rescale(truth, pred)
        assert not np.isfinite(corr)


class TestImputeZBlock:
    def test_returns_imputed_array_and_corr(self, small_ld, rng):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64)
        z[:5] = np.nan
        z_imp, corr = impute_z_block(z, vecs, vals, min_cor=0.0)
        assert z_imp is not None
        assert np.isfinite(corr)
        assert len(z_imp) == small_ld.shape[0]

    def test_quality_gate_rejects_low_corr(self, small_ld, rng):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64)
        z[:5] = np.nan
        z_imp, _ = impute_z_block(z, vecs, vals, min_cor=1.0)
        assert z_imp is None

    def test_clamps_to_observed_range(self, small_ld, rng):
        vals, vecs = ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64) * 3
        z[:3] = np.nan
        obs_max = float(np.nanmax(np.abs(z)))
        z_imp, _ = impute_z_block(z, vecs, vals, min_cor=0.0)
        if z_imp is not None:
            assert float(np.max(np.abs(z_imp))) <= obs_max + 1e-4


class TestScalarNSe:
    def test_basic_output(self):
        se_obs = np.array([0.05, 0.06, 0.04])
        eaf_obs = np.array([0.3, 0.4, 0.5])
        eaf_ref = np.array([0.2, 0.3, 0.5])
        se_ref = scalar_n_se(se_obs, eaf_obs, eaf_ref)
        assert len(se_ref) == 3
        assert np.all(np.isfinite(se_ref))
        assert np.all(se_ref > 0)

    def test_nan_for_eaf_zero_or_one(self):
        se_obs = np.array([0.05])
        eaf_obs = np.array([0.3])
        eaf_ref = np.array([0.0, 0.5, 1.0])
        se_ref = scalar_n_se(se_obs, eaf_obs, eaf_ref)
        assert not np.isfinite(se_ref[0])
        assert np.isfinite(se_ref[1])
        assert not np.isfinite(se_ref[2])

    def test_nan_when_no_valid_observed(self):
        se_obs = np.array([0.05])
        eaf_obs = np.array([np.nan])
        eaf_ref = np.array([0.3])
        se_ref = scalar_n_se(se_obs, eaf_obs, eaf_ref)
        assert not np.isfinite(se_ref[0])

    def test_scale_consistent_with_eaf_formula(self):
        # If se_obs * sqrt(2*eaf*(1-eaf)) = C for all obs,
        # then se_ref = C / sqrt(2*eaf_ref*(1-eaf_ref))
        eaf_obs = np.array([0.3, 0.4])
        C = 0.5
        se_obs = C / np.sqrt(2 * eaf_obs * (1 - eaf_obs))
        eaf_ref = np.array([0.2, 0.6])
        se_ref = scalar_n_se(se_obs, eaf_obs, eaf_ref)
        expected = C / np.sqrt(2 * eaf_ref * (1 - eaf_ref))
        np.testing.assert_allclose(se_ref, expected, rtol=1e-3)
