import pytest
from pytest import raises as assert_raises, warns as assert_warns

import numpy as np
from numpy.testing import assert_approx_equal, assert_allclose, assert_equal

from scipy.spatial.distance import cdist
from scipy import stats

class TestMGCErrorWarnings:
    """ Tests errors and warnings derived from MGC.
    """
    def test_error_notndarray(self):
        # raises error if x or y is not a ndarray
        x = np.arange(20)
        y = [5] * 20
        assert_raises(ValueError, stats.multiscale_graphcorr, x, y)
        assert_raises(ValueError, stats.multiscale_graphcorr, y, x)

    def test_error_shape(self):
        # raises error if number of samples different (n)
        x = np.arange(100).reshape(25, 4)
        y = x.reshape(10, 10)
        assert_raises(ValueError, stats.multiscale_graphcorr, x, y)

    def test_error_lowsamples(self):
        # raises error if samples are low (< 3)
        x = np.arange(3)
        y = np.arange(3)
        assert_raises(ValueError, stats.multiscale_graphcorr, x, y)

    def test_error_nans(self):
        # raises error if inputs contain NaNs
        x = np.arange(20, dtype=float)
        x[0] = np.nan
        assert_raises(ValueError, stats.multiscale_graphcorr, x, x)

        y = np.arange(20)
        assert_raises(ValueError, stats.multiscale_graphcorr, x, y)

    def test_error_wrongdisttype(self):
        # raises error if metric is not a function
        x = np.arange(20)
        compute_distance = 0
        assert_raises(ValueError, stats.multiscale_graphcorr, x, x,
                      compute_distance=compute_distance)

    @pytest.mark.parametrize("reps", [
        -1,    # reps is negative
        '1',   # reps is not integer
    ])
    def test_error_reps(self, reps):
        # raises error if reps is negative
        x = np.arange(20)
        assert_raises(ValueError, stats.multiscale_graphcorr, x, x, reps=reps)

    def test_warns_reps(self):
        # raises warning when reps is less than 1000
        x = np.arange(20)
        reps = 100
        assert_warns(RuntimeWarning, stats.multiscale_graphcorr, x, x, reps=reps)

    def test_error_infty(self):
        # raises error if input contains infinities
        x = np.arange(20)
        y = np.ones(20) * np.inf
        assert_raises(ValueError, stats.multiscale_graphcorr, x, y)


class TestMGCStat:
    """ Test validity of MGC test statistic
    """
    def setup_method(self):
        self.rng = np.random.default_rng(1266219746)

    def _simulations(self, samps=100, dims=1, sim_type="", rng=None):
        rng = rng or self.rng
        # linear simulation
        if sim_type == "linear":
            x = rng.uniform(-1, 1, size=(samps, 1))
            y = x + 0.3 * rng.random(size=(x.size, 1))

        # spiral simulation
        elif sim_type == "nonlinear":
            unif = np.array(self.rng.uniform(0, 5, size=(samps, 1)))
            x = unif * np.cos(np.pi * unif)
            y = (unif * np.sin(np.pi * unif) +
                 0.4*self.rng.random(size=(x.size, 1)))

        # independence (tests type I simulation)
        elif sim_type == "independence":
            u = self.rng.standard_normal(size=(samps, 1))
            v = self.rng.standard_normal(size=(samps, 1))
            u_2 = self.rng.binomial(1, p=0.5, size=(samps, 1))
            v_2 = self.rng.binomial(1, p=0.5, size=(samps, 1))
            x = u/3 + 2*u_2 - 1
            y = v/3 + 2*v_2 - 1

        # raises error if not approved sim_type
        else:
            raise ValueError("sim_type must be linear, nonlinear, or "
                             "independence")

        # add dimensions of noise for higher dimensions
        if dims > 1:
            dims_noise = self.rng.standard_normal(size=(samps, dims-1))
            x = np.concatenate((x, dims_noise), axis=1)

        return x, y

    @pytest.mark.xslow
    @pytest.mark.parametrize("sim_type, obs_stat, obs_pvalue", [
        # Reference values produced by `multiscale_graphcorr` ->
        # this only guards against unintended changes
        ("linear", 0.970137188683, 0.000999000999),      # test linear simulation
        ("nonlinear", 0.149843048137, 0.000999000999),   # test spiral simulation
        ("independence", -0.003965353120, 0.4675324675)  # test independence simulation
    ])
    def test_oned(self, sim_type, obs_stat, obs_pvalue):
        # generate x and y
        rng = np.random.default_rng(8157117705)
        x, y = self._simulations(samps=100, dims=1, sim_type=sim_type, rng=rng)

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, random_state=rng)
        assert_allclose(stat, obs_stat)
        assert_allclose(pvalue, obs_pvalue)

    @pytest.mark.xslow
    @pytest.mark.parametrize("sim_type, obs_stat, obs_pvalue", [
        # Reference values produced by `multiscale_graphcorr` ->
        # this only guards against unintended changes
        ("linear", 0.193063672972, 0.000999000999),     # test linear simulation
        ("nonlinear", 0.010042844939, 0.216783216783),  # test spiral simulation
    ])
    def test_fived(self, sim_type, obs_stat, obs_pvalue):
        # generate x and y
        rng = np.random.default_rng(8157117705)
        x, y = self._simulations(samps=100, dims=5, sim_type=sim_type, rng=rng)

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, random_state=rng)
        assert_allclose(stat, obs_stat)
        assert_allclose(pvalue, obs_pvalue)

    @pytest.mark.xslow
    def test_twosamp(self):
        # generate x and y
        x = self.rng.binomial(100, 0.5, size=(100, 5))
        y = self.rng.standard_normal(size=(80, 5))

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, random_state=self.rng)
        assert_approx_equal(stat, 1.0, significant=1)
        assert_approx_equal(pvalue, 0.001, significant=1)

        # generate x and y
        y = self.rng.standard_normal(size=(100, 5))

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, is_twosamp=True,
                                                     random_state=self.rng)
        assert_approx_equal(stat, 1.0, significant=1)
        assert_approx_equal(pvalue, 0.001, significant=1)

    @pytest.mark.xslow
    def test_workers(self):
        # generate x and y
        x, y = self._simulations(samps=100, dims=1, sim_type="linear")

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, workers=2,
                                                     random_state=self.rng)
        assert_approx_equal(stat, 0.97, significant=1)
        assert_approx_equal(pvalue, 0.001, significant=1)

    @pytest.mark.xslow
    def test_random_state(self):
        # generate x and y
        x, y = self._simulations(samps=100, dims=1, sim_type="linear")

        # test stat and pvalue
        stat, pvalue, _ = stats.multiscale_graphcorr(x, y, random_state=1)
        assert_approx_equal(stat, 0.97, significant=1)
        assert_approx_equal(pvalue, 0.001, significant=1)

    @pytest.mark.xslow
    def test_dist_perm(self):
        # generate x and y
        x, y = self._simulations(samps=100, dims=1, sim_type="nonlinear")
        distx = cdist(x, x, metric="euclidean")
        disty = cdist(y, y, metric="euclidean")

        stat_dist, pvalue_dist, _ = stats.multiscale_graphcorr(distx, disty,
                                                               compute_distance=None,
                                                               random_state=1)
        assert_approx_equal(stat_dist, 0.163, significant=1)
        assert_approx_equal(pvalue_dist, 0.001, significant=1)

    @pytest.mark.fail_slow(20)  # all other tests are XSLOW; we need at least one to run
    @pytest.mark.slow
    def test_pvalue_literature(self):
        # generate x and y
        x, y = self._simulations(samps=100, dims=1, sim_type="linear")

        # test stat and pvalue
        _, pvalue, _ = stats.multiscale_graphcorr(x, y, random_state=1)
        assert_allclose(pvalue, 1/1001)

    @pytest.mark.xslow
    def test_alias(self):
        # generate x and y
        x, y = self._simulations(samps=100, dims=1, sim_type="linear")

        res = stats.multiscale_graphcorr(x, y, random_state=1)
        assert_equal(res.stat, res.statistic)
