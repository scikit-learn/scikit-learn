import warnings

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from scipy._lib._util import rng_integers
from scipy._lib._array_api import (is_numpy, make_xp_test_case, xp_default_dtype,
                                   xp_size, array_namespace, _xp_copy_to_numpy)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal
from scipy._lib import array_api_extra as xpx
from scipy import stats, special
from scipy.fft.tests.test_fftlog import skip_xp_backends
from scipy.optimize import root

from scipy.stats import bootstrap, monte_carlo_test, permutation_test, power
import scipy.stats._resampling as _resampling


@make_xp_test_case(bootstrap)
class TestBootstrap:
    @skip_xp_backends('numpy', reason='NumPy does not raise')
    def test_bootstrap_iv_other(self, xp):
        message = f"When using array library {xp.__name__}"
        with pytest.raises(TypeError, match=message):
            bootstrap((xp.asarray([1, 2, 3]),), lambda x: xp.mean(x))

    def test_bootstrap_iv(self, xp):
        sample = xp.asarray([1, 2, 3])

        message = "`data` must contain at least one sample."
        with pytest.raises(ValueError, match=message):
            bootstrap(tuple(), xp.mean)

        message = "each sample in `data` must contain two or more observations..."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample, xp.asarray([1])), xp.mean)

        message = ("When `paired is True`, all samples must have the same length ")
        with pytest.raises(ValueError, match=message):
            bootstrap((sample, xp.asarray([1, 2, 3, 4])), xp.mean, paired=True)

        message = "`vectorized` must be `True`, `False`, or `None`."
        with pytest.raises(ValueError, match=message):
            bootstrap(sample, xp.mean, vectorized='ekki')

        message = "`axis` must be an integer."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, axis=1.5)

        message = "could not convert string to float"
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, confidence_level='ni')

        message = "`n_resamples` must be a non-negative integer."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, n_resamples=-1000)

        message = "`n_resamples` must be a non-negative integer."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, n_resamples=1000.5)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, batch=-1000)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, batch=1000.5)

        message = "`method` must be in"
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, method='ekki')

        message = "`bootstrap_result` must have attribute `bootstrap_distribution'"
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, bootstrap_result=10)

        message = "Either `bootstrap_result.bootstrap_distribution.size`"
        with pytest.raises(ValueError, match=message):
            bootstrap((sample,), xp.mean, n_resamples=0)

        message = "SeedSequence expects int or sequence of ints"
        with pytest.raises(TypeError, match=message):
            bootstrap((sample,), xp.mean, rng='herring')

    @pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_bootstrap_batch(self, method, axis, xp):
        # for one-sample statistics, batch size shouldn't affect the result
        rng = np.random.RandomState(0)

        x = rng.rand(10, 11, 12)
        # SPEC-007 leave one call with random_state to ensure it still works
        res1 = bootstrap((xp.asarray(x),), xp.mean, batch=None, method=method,
                         random_state=0, axis=axis, n_resamples=100)
        rng = np.random.RandomState(0)
        res2 = bootstrap((xp.asarray(x),), xp.mean, batch=10, method=method,
                         axis=axis, n_resamples=100, random_state=rng)

        xp_assert_equal(res2.confidence_interval.low, res1.confidence_interval.low)
        xp_assert_equal(res2.confidence_interval.high, res1.confidence_interval.high)
        xp_assert_equal(res2.standard_error, res1.standard_error)

    @pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
    def test_bootstrap_paired(self, method, xp):
        # test that `paired` works as expected
        rng = np.random.RandomState(0)
        n = 100
        x = xp.asarray(rng.rand(n))
        y = xp.asarray(rng.rand(n))

        def my_statistic(x, y, axis=-1):
            return xp.mean((x-y)**2, axis=axis)

        def my_paired_statistic(i, axis=-1):
            a = x[i]
            b = y[i]
            res = my_statistic(a, b)
            return res

        i = xp.arange(x.shape[0])

        res1 = bootstrap((i,), my_paired_statistic, rng=0)
        res2 = bootstrap((x, y), my_statistic, paired=True, rng=0)

        xp_assert_close(res1.confidence_interval.low, res2.confidence_interval.low)
        xp_assert_close(res1.confidence_interval.high, res2.confidence_interval.high)
        xp_assert_close(res1.standard_error, res2.standard_error)

    @pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
    @pytest.mark.parametrize("axis", [0, 1, 2])
    @pytest.mark.parametrize("paired", [True, False])
    def test_bootstrap_vectorized(self, method, axis, paired, xp):
        # test that paired is vectorized as expected: when samples are tiled,
        # CI and standard_error of each axis-slice is the same as those of the
        # original 1d sample

        rng = np.random.RandomState(0)

        def my_statistic(x, y, z, axis=-1):
            return xp.mean(x, axis=axis) + xp.mean(y, axis=axis) + xp.mean(z, axis=axis)

        shape = 10, 11, 12
        n_samples = shape[axis]

        x = rng.rand(n_samples)
        y = rng.rand(n_samples)
        z = rng.rand(n_samples)
        x, y, z = xp.asarray(x), xp.asarray(y), xp.asarray(z)

        res1 = bootstrap((x, y, z), my_statistic, paired=paired, method=method,
                         rng=0, axis=0, n_resamples=100)
        assert (res1.bootstrap_distribution.shape
                == res1.standard_error.shape + (100,))

        reshape = [1, 1, 1]
        reshape[axis] = n_samples
        reshape = tuple(reshape)
        x = xp.broadcast_to(xp.reshape(x, reshape), shape)
        y = xp.broadcast_to(xp.reshape(y, reshape), shape)
        z = xp.broadcast_to(xp.reshape(z, reshape), shape)
        res2 = bootstrap((x, y, z), my_statistic, paired=paired, method=method,
                         rng=0, axis=axis, n_resamples=100)

        result_shape = list(shape)
        result_shape.pop(axis)

        ref_ci_low = xp.broadcast_to(res2.confidence_interval.low, result_shape)
        ref_ci_high = xp.broadcast_to(res2.confidence_interval.high, result_shape)
        ref_ci_standard_error = xp.broadcast_to(res2.standard_error, result_shape)

        xp_assert_close(res2.confidence_interval.low, ref_ci_low)
        xp_assert_close(res2.confidence_interval.high, ref_ci_high)
        xp_assert_close(res2.standard_error, ref_ci_standard_error)

    @pytest.mark.slow
    @pytest.mark.xfail_on_32bit("MemoryError with BCa observed in CI")
    @pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
    def test_bootstrap_against_theory(self, method, xp):
        # based on https://www.statology.org/confidence-intervals-python/
        rng = np.random.default_rng(2442101192988600726)
        data = stats.norm.rvs(loc=5, scale=2, size=5000, random_state=rng)
        alpha = 0.95
        dist = stats.t(df=len(data)-1, loc=np.mean(data), scale=stats.sem(data))
        ref_low, ref_high = dist.interval(confidence=alpha)
        ref_se = dist.std()

        config = dict(data=(xp.asarray(data),), statistic=xp.mean, n_resamples=5000,
                      method=method, rng=rng)
        res = bootstrap(**config, confidence_level=alpha)

        xp_assert_close(res.confidence_interval.low, xp.asarray(ref_low), rtol=5e-4)
        xp_assert_close(res.confidence_interval.high, xp.asarray(ref_high), rtol=5e-4)
        xp_assert_close(res.standard_error, xp.asarray(ref_se), atol=3e-4)

        config.update(dict(n_resamples=0, bootstrap_result=res))
        res = bootstrap(**config, confidence_level=alpha, alternative='less')
        xp_assert_close(res.confidence_interval.high,
                        xp.asarray(dist.ppf(alpha)), rtol=5e-4)

        config.update(dict(n_resamples=0, bootstrap_result=res))
        res = bootstrap(**config, confidence_level=alpha, alternative='greater')
        xp_assert_close(res.confidence_interval.low,
                        xp.asarray(dist.ppf(1-alpha)), rtol=5e-4)

    tests_R = [("basic", 23.77, 79.12),
               ("percentile", 28.86, 84.21),
               ("BCa", 32.31, 91.43)]

    @pytest.mark.parametrize("method, ref_low, ref_high", tests_R)
    def test_bootstrap_against_R(self, method, ref_low, ref_high, xp):
        # Compare against R's "boot" library
        # library(boot)

        # stat <- function (x, a) {
        #     mean(x[a])
        # }

        # x <- c(10, 12, 12.5, 12.5, 13.9, 15, 21, 22,
        #        23, 34, 50, 81, 89, 121, 134, 213)

        # # Use a large value so we get a few significant digits for the CI.
        # n = 1000000
        # bootresult = boot(x, stat, n)
        # result <- boot.ci(bootresult)
        # print(result)
        x = xp.asarray([10, 12, 12.5, 12.5, 13.9, 15, 21, 22,
                        23, 34, 50, 81, 89, 121, 134, 213])
        res = bootstrap((x,), xp.mean, n_resamples=1000000, method=method, rng=0)
        xp_assert_close(res.confidence_interval.low, xp.asarray(ref_low), rtol=0.005)
        xp_assert_close(res.confidence_interval.high, xp.asarray(ref_high), rtol=0.005)

    def test_multisample_BCa_against_R(self, xp):
        # Because bootstrap is stochastic, it's tricky to test against reference
        # behavior. Here, we show that SciPy's BCa CI matches R wboot's BCa CI
        # much more closely than the other SciPy CIs do.

        # arbitrary skewed data
        x = xp.asarray([0.75859206, 0.5910282, -0.4419409, -0.36654601,
                        0.34955357, -1.38835871, 0.76735821])
        y = xp.asarray([1.41186073, 0.49775975, 0.08275588, 0.24086388,
                        0.03567057, 0.52024419, 0.31966611, 1.32067634])

        # a multi-sample statistic for which the BCa CI tends to be different
        # from the other CIs
        def statistic(x, y, axis):
            s1 = stats.skew(x, axis=axis)
            s2 = stats.skew(y, axis=axis)
            return s1 - s2

        # compute confidence intervals using each method
        rng = np.random.default_rng(468865032284792692)

        res_basic = stats.bootstrap((x, y), statistic, method='basic',
                                    batch=100, rng=rng)
        res_percent = stats.bootstrap((x, y), statistic, method='percentile',
                                      batch=100, rng=rng)
        res_bca = stats.bootstrap((x, y), statistic, method='bca',
                                  batch=100, rng=rng)

        # compute midpoints so we can compare just one number for each
        mid_basic = xp.mean(xp.stack(res_basic.confidence_interval))
        mid_percent = xp.mean(xp.stack(res_percent.confidence_interval))
        mid_bca = xp.mean(xp.stack(res_bca.confidence_interval))

        # reference for BCA CI computed using R wboot package:
        # library(wBoot)
        # library(moments)

        # x = c(0.75859206, 0.5910282, -0.4419409, -0.36654601,
        #       0.34955357, -1.38835871,  0.76735821)
        # y = c(1.41186073, 0.49775975, 0.08275588, 0.24086388,
        #       0.03567057, 0.52024419, 0.31966611, 1.32067634)

        # twoskew <- function(x1, y1) {skewness(x1) - skewness(y1)}
        # boot.two.bca(x, y, skewness, conf.level = 0.95,
        #              R = 9999, stacked = FALSE)
        mid_wboot = -1.5519

        # compute percent difference relative to wboot BCA method
        diff_basic = (mid_basic - mid_wboot)/abs(mid_wboot)
        diff_percent = (mid_percent - mid_wboot)/abs(mid_wboot)
        diff_bca = (mid_bca - mid_wboot)/abs(mid_wboot)

        # SciPy's BCa CI midpoint is much closer than that of the other methods
        assert diff_basic < -0.15
        assert diff_percent > 0.15
        assert abs(diff_bca) < 0.03

    def test_BCa_acceleration_against_reference(self, xp):
        # Compare the (deterministic) acceleration parameter for a multi-sample
        # problem against a reference value. The example is from [1], but Efron's
        # value seems inaccurate. Straightforward code for computing the
        # reference acceleration (0.011008228344026734) is available at:
        # https://github.com/scipy/scipy/pull/16455#issuecomment-1193400981

        y = xp.asarray([10., 27., 31., 40., 46., 50., 52., 104., 146.])
        z = xp.asarray([16., 23., 38., 94., 99., 141., 197.])

        def statistic(z, y, axis=0):
            return xp.mean(z, axis=axis) - xp.mean(y, axis=axis)

        data = [z, y]
        res = stats.bootstrap(data, statistic)

        axis = -1
        alpha = 0.95
        theta_hat_b = res.bootstrap_distribution
        batch = 100
        _, _, a_hat = _resampling._bca_interval(data, statistic, axis, alpha,
                                                theta_hat_b, batch, xp)
        xp_assert_close(a_hat, xp.asarray(0.011008228344026734))

    tests_against_itself_1samp = {"basic": 1789,
                                  "percentile": 1790,
                                  "BCa": 1789}

    @pytest.mark.slow
    @pytest.mark.parametrize("method, expected",
                             tests_against_itself_1samp.items())
    def test_bootstrap_against_itself_1samp(self, method, expected, xp):
        # The expected values in this test were generated using bootstrap
        # to check for unintended changes in behavior. The test also makes sure
        # that bootstrap works with multi-sample statistics and that the
        # `axis` argument works as expected / function is vectorized.
        rng = np.random.default_rng(9123847)

        n = 100  # size of sample
        n_resamples = 999  # number of bootstrap resamples used to form each CI
        confidence_level = 0.9

        # The true mean is 5
        dist = stats.norm(loc=5, scale=1)
        stat_true = float(dist.mean())

        # Do the same thing 2000 times. (The code is fully vectorized.)
        n_replications = 2000
        data = dist.rvs(size=(n_replications, n), random_state=rng)
        res = bootstrap((xp.asarray(data),),
                        statistic=xp.mean,
                        confidence_level=confidence_level,
                        n_resamples=n_resamples,
                        batch=50,
                        method=method,
                        axis=-1,
                        rng=rng)
        ci = res.confidence_interval

        # ci contains vectors of lower and upper confidence interval bounds
        ci_contains_true = xp.count_nonzero((ci[0] < stat_true) & (stat_true < ci[1]))
        assert ci_contains_true == expected

        # ci_contains_true is not inconsistent with confidence_level
        pvalue = stats.binomtest(int(ci_contains_true), n_replications,
                                 confidence_level).pvalue
        assert pvalue > 0.1

    tests_against_itself_2samp = {"basic": 892,
                                  "percentile": 890}

    @pytest.mark.slow
    @pytest.mark.parametrize("method, expected",
                             tests_against_itself_2samp.items())
    def test_bootstrap_against_itself_2samp(self, method, expected, xp):
        # The expected values in this test were generated using bootstrap
        # to check for unintended changes in behavior. The test also makes sure
        # that bootstrap works with multi-sample statistics and that the
        # `axis` argument works as expected / function is vectorized.
        rng = np.random.RandomState(0)

        n1 = 100  # size of sample 1
        n2 = 120  # size of sample 2
        n_resamples = 999  # number of bootstrap resamples used to form each CI
        confidence_level = 0.9

        # The statistic we're interested in is the difference in means
        def my_stat(data1, data2, axis=-1):
            mean1 = xp.mean(data1, axis=axis)
            mean2 = xp.mean(data2, axis=axis)
            return mean1 - mean2

        # The true difference in the means is -0.1
        dist1 = stats.norm(loc=0, scale=1)
        dist2 = stats.norm(loc=0.1, scale=1)
        stat_true = float(dist1.mean() - dist2.mean())

        # Do the same thing 1000 times. (The code is fully vectorized.)
        n_replications = 1000
        data1 = dist1.rvs(size=(n_replications, n1), random_state=rng)
        data2 = dist2.rvs(size=(n_replications, n2), random_state=rng)
        res = bootstrap((xp.asarray(data1), xp.asarray(data2)),
                        statistic=my_stat,
                        confidence_level=confidence_level,
                        n_resamples=n_resamples,
                        batch=50,
                        method=method,
                        axis=-1,
                        random_state=rng)
        ci = res.confidence_interval

        # ci contains vectors of lower and upper confidence interval bounds
        ci_contains_true = xp.count_nonzero((ci[0] < stat_true) & (stat_true < ci[1]))
        assert ci_contains_true == expected

        # ci_contains_true is not inconsistent with confidence_level
        pvalue = stats.binomtest(int(ci_contains_true), n_replications,
                                 confidence_level).pvalue
        assert pvalue > 0.1

    @pytest.mark.parametrize("method", ["basic", "percentile", "BCa"])
    @pytest.mark.parametrize("axis", [0, 1])
    @pytest.mark.parametrize("dtype", ['float32', 'float64'])
    def test_bootstrap_vectorized_3samp(self, method, axis, dtype, xp):
        def statistic(*data, axis=0):
            # an arbitrary, vectorized statistic
            return sum(xp.mean(sample, axis=axis) for sample in data)

        def statistic_1d(*data):
            # the same statistic, not vectorized
            for sample in data:
                assert sample.ndim == 1
            return np.asarray(sum(sample.mean() for sample in data), dtype=dtype)

        rng = np.random.RandomState(0)
        x = rng.rand(4, 5).astype(dtype)
        y = rng.rand(4, 5).astype(dtype)
        z = rng.rand(4, 5).astype(dtype)
        res1 = bootstrap((xp.asarray(x), xp.asarray(y), xp.asarray(z)),
                         statistic, vectorized=True, axis=axis, n_resamples=100,
                         method=method, rng=0)
        res2 = bootstrap((x, y, z), statistic_1d, vectorized=False,
                         axis=axis, n_resamples=100, method=method, rng=0)

        rtol = 1e-6 if dtype == 'float32' else 1e-14
        xp_assert_close(res1.confidence_interval.low,
                        xp.asarray(res2.confidence_interval.low), rtol=rtol)
        xp_assert_close(res1.confidence_interval.high,
                        xp.asarray(res2.confidence_interval.high), rtol=rtol)
        xp_assert_close(res1.standard_error, xp.asarray(res2.standard_error), rtol=rtol)

    @pytest.mark.xfail_on_32bit("Failure is not concerning; see gh-14107")
    @pytest.mark.parametrize("method", ["basic", "percentile", "BCa"])
    @pytest.mark.parametrize("axis", [0, 1])
    def test_bootstrap_vectorized_1samp(self, method, axis, xp):
        def statistic(x, axis=0):
            # an arbitrary, vectorized statistic
            return xp.mean(x, axis=axis)

        def statistic_1d(x):
            # the same statistic, not vectorized
            assert x.ndim == 1
            return x.mean(axis=0)

        rng = np.random.default_rng(7939952824)
        x = rng.random((2, 3))
        res1 = bootstrap((xp.asarray(x),), statistic, vectorized=True, axis=axis,
                         n_resamples=100, batch=None, method=method, rng=0)
        res2 = bootstrap((x,), statistic_1d, vectorized=False, axis=axis,
                         n_resamples=100, batch=10, method=method, rng=0)
        xp_assert_close(res1.confidence_interval.low,
                        xp.asarray(res2.confidence_interval.low))
        xp_assert_close(res1.confidence_interval.high,
                        xp.asarray(res2.confidence_interval.high))
        xp_assert_close(res1.standard_error, xp.asarray(res2.standard_error))

    @pytest.mark.parametrize("method", ["basic", "percentile", "BCa"])
    def test_bootstrap_degenerate(self, method, xp):
        data = xp.full((35,), 10000.)
        if method == "BCa":
            with np.errstate(invalid='ignore'):
                msg = "The BCa confidence interval cannot be calculated"
                with pytest.warns(stats.DegenerateDataWarning, match=msg):
                    res = bootstrap([data, ], xp.mean, method=method)
                    xp_assert_equal(res.confidence_interval.low, xp.asarray(xp.nan))
                    xp_assert_equal(res.confidence_interval.high, xp.asarray(xp.nan))
        else:
            res = bootstrap([data, ], xp.mean, method=method)
            xp_assert_equal(res.confidence_interval.low, xp.asarray(10000.))
            xp_assert_equal(res.confidence_interval.high, xp.asarray(10000.))
        xp_assert_equal(res.standard_error, xp.asarray(0.))

    @pytest.mark.parametrize("method", ["BCa", "basic", "percentile"])
    def test_bootstrap_gh15678(self, method, xp):
        # Check that gh-15678 is fixed: when statistic function returned a Python
        # float, method="BCa" failed when trying to add a dimension to the float
        rng = np.random.default_rng(354645618886684)
        dist = stats.norm(loc=2, scale=4)
        data = dist.rvs(size=100, random_state=rng)
        res = bootstrap((xp.asarray(data),), stats.skew, method=method, n_resamples=100,
                        rng=np.random.default_rng(9563))
        # this always worked because np.apply_along_axis returns NumPy data type
        ref = bootstrap((data,), stats.skew, method=method, n_resamples=100,
                        rng=np.random.default_rng(9563), vectorized=False)
        xp_assert_close(res.confidence_interval.low,
                        xp.asarray(ref.confidence_interval.low))
        xp_assert_close(res.confidence_interval.high,
                        xp.asarray(ref.confidence_interval.high))
        xp_assert_close(res.standard_error, xp.asarray(ref.standard_error))

    def test_bootstrap_min(self, xp):
        # Check that gh-15883 is fixed: percentileofscore should
        # behave according to the 'mean' behavior and not trigger nan for BCa
        rng = np.random.default_rng(1891289180021102)
        dist = stats.norm(loc=2, scale=4)
        data = dist.rvs(size=100, random_state=rng)
        true_min = np.min(data)
        data = xp.asarray(data)
        res = bootstrap((data,), xp.min, method="BCa", n_resamples=100,
                        rng=np.random.default_rng(3942))
        xp_assert_equal(res.confidence_interval.low, xp.asarray(true_min))
        res2 = bootstrap((-data,), xp.max, method="BCa", n_resamples=100,
                         rng=np.random.default_rng(3942))
        xp_assert_close(-res.confidence_interval.low, res2.confidence_interval.high)
        xp_assert_close(-res.confidence_interval.high, res2.confidence_interval.low)

    @pytest.mark.parametrize("additional_resamples", [0, 1000])
    def test_re_bootstrap(self, additional_resamples, xp):
        # Test behavior of parameter `bootstrap_result`
        rng = np.random.default_rng(8958153316228384)
        x = rng.random(size=100)

        n1 = 1000
        n2 = additional_resamples
        n3 = n1 + additional_resamples

        rng = np.random.default_rng(296689032789913033)
        res = stats.bootstrap((xp.asarray(x),), xp.mean, n_resamples=n1, rng=rng,
                              confidence_level=0.95, method='percentile')
        res = stats.bootstrap((xp.asarray(x),), xp.mean, n_resamples=n2, rng=rng,
                              confidence_level=0.90, method='BCa',
                              bootstrap_result=res)

        rng = np.random.default_rng(296689032789913033)
        ref = stats.bootstrap((xp.asarray(x),), xp.mean, n_resamples=n3, rng=rng,
                              confidence_level=0.90, method='BCa')

        xp_assert_close(res.confidence_interval.low,
                        xp.asarray(ref.confidence_interval.low),
                        rtol=1e-14)
        xp_assert_close(res.confidence_interval.high,
                        xp.asarray(ref.confidence_interval.high),
                        rtol=1e-14)
        xp_assert_close(res.standard_error, xp.asarray(ref.standard_error), rtol=1e-14)

    @pytest.mark.xfail_on_32bit("Sensitive to machine precision")
    @pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
    def test_bootstrap_alternative(self, method, xp):
        rng = np.random.default_rng(5894822712842015040)
        dist = stats.norm(loc=2, scale=4)
        data = (xp.asarray(dist.rvs(size=(100), random_state=rng)),)

        config = dict(data=data, statistic=xp.std, rng=rng, axis=-1)
        t = stats.bootstrap(**config, confidence_level=0.9)

        config.update(dict(n_resamples=0, bootstrap_result=t))
        l = stats.bootstrap(**config, confidence_level=0.95, alternative='less')
        g = stats.bootstrap(**config, confidence_level=0.95, alternative='greater')

        xp_assert_close(l.confidence_interval.high, t.confidence_interval.high,
                        rtol=1e-14)
        xp_assert_close(g.confidence_interval.low, t.confidence_interval.low,
                        rtol=1e-14)
        assert xp.isinf(l.confidence_interval.low) and l.confidence_interval.low < 0
        assert xp.isinf(g.confidence_interval.high) and g.confidence_interval.high > 0

        with pytest.raises(ValueError, match='`alternative` must be one of'):
            stats.bootstrap(**config, alternative='ekki-ekki')

    def test_jackknife_resample(self, xp):
        shape = 3, 4, 5, 6
        rng = np.random.default_rng(5274950392)
        x = rng.random(size=shape)
        y = next(_resampling._jackknife_resample(xp.asarray(x), xp=xp))

        for i in range(shape[-1]):
            # each resample is indexed along second to last axis
            # (last axis is the one the statistic will be taken over / consumed)
            slc = y[..., i, :]
            expected = np.delete(x, i, axis=-1)

            xp_assert_equal(slc, xp.asarray(expected))

        y2 = list(_resampling._jackknife_resample(xp.asarray(x), batch=2, xp=xp))
        xp_assert_equal(xp.concat(y2, axis=-2), y)

    @pytest.mark.skip_xp_backends("array_api_strict",
                                  reason="Test uses ... + fancy indexing")
    @pytest.mark.parametrize("rng_name", ["RandomState", "default_rng"])
    def test_bootstrap_resample(self, rng_name, xp):
        rng_gen = getattr(np.random, rng_name)
        rng1 = rng_gen(3949441460)
        rng2 = rng_gen(3949441460)

        n_resamples = 10
        shape = 3, 4, 5, 6

        rng = np.random.default_rng(5274950392)
        x = xp.asarray(rng.random(shape))
        y = _resampling._bootstrap_resample(x, n_resamples, rng=rng1, xp=xp)

        for i in range(n_resamples):
            # each resample is indexed along second to last axis
            # (last axis is the one the statistic will be taken over / consumed)
            slc = y[..., i, :]

            js = xp.asarray(rng_integers(rng2, 0, shape[-1], shape[-1]))
            expected = x[..., js]

            xp_assert_equal(slc, expected)

    @pytest.mark.parametrize("score", [0, 0.5, 1])
    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_percentile_of_score(self, score, axis, xp):
        shape = 10, 20, 30
        rng = np.random.default_rng(5903363153)
        x = rng.random(shape)
        dtype = xp_default_dtype(xp)
        p = _resampling._percentile_of_score(xp.asarray(x, dtype=dtype),
                                             xp.asarray(score, dtype=dtype),
                                             axis=-1, xp=xp)

        def vectorized_pos(a, score, axis):
            return np.apply_along_axis(stats.percentileofscore, axis, a, score)

        p2 = vectorized_pos(x, score, axis=-1)/100

        xp_assert_close(p, xp.asarray(p2, dtype=dtype), rtol=1e-15)

    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_vectorize_statistic(self, axis):
        # test that _vectorize_statistic vectorizes a statistic along `axis`
        # only used with NumPy arrays

        def statistic(*data, axis):
            # an arbitrary, vectorized statistic
            return sum(sample.mean(axis) for sample in data)

        def statistic_1d(*data):
            # the same statistic, not vectorized
            for sample in data:
                assert sample.ndim == 1
            return statistic(*data, axis=0)

        # vectorize the non-vectorized statistic
        statistic2 = _resampling._vectorize_statistic(statistic_1d)

        rng = np.random.RandomState(0)
        x = rng.rand(4, 5, 6)
        y = rng.rand(4, 1, 6)
        z = rng.rand(1, 5, 6)

        res1 = statistic(x, y, z, axis=axis)
        res2 = statistic2(x, y, z, axis=axis)
        assert_allclose(res1, res2)

    @pytest.mark.slow
    @pytest.mark.parametrize("method", ["basic", "percentile", "BCa"])
    def test_vector_valued_statistic(self, method, xp):
        # Generate 95% confidence interval around MLE of normal distribution
        # parameters. Repeat 100 times, each time on sample of size 100.
        # Check that confidence interval contains true parameters ~95 times.
        # Confidence intervals are estimated and stochastic; a test failure
        # does not necessarily indicate that something is wrong. More important
        # than values of `counts` below is that the shapes of the outputs are
        # correct.

        rng = np.random.default_rng(2196847219)
        params = 1, 0.5
        sample = xp.asarray(rng.normal(*params, size=(100, 100)))

        def statistic(data, axis):
            return xp.stack([xp.mean(data, axis=axis),
                             xp.std(data, axis=axis, correction=1)])

        res = bootstrap((sample,), statistic, method=method, axis=-1,
                        n_resamples=9999, batch=200, random_state=rng)

        params = xp.asarray([1, 0.5])
        counts = xp.count_nonzero((res.confidence_interval.low.T < params)
                                  & (res.confidence_interval.high.T > params),
                                  axis=0)
        assert xp.all(counts >= 90)
        assert xp.all(counts <= 100)
        assert res.confidence_interval.low.shape == (2, 100)
        assert res.confidence_interval.high.shape == (2, 100)
        assert res.standard_error.shape == (2, 100)
        assert res.bootstrap_distribution.shape == (2, 100, 9999)

    @pytest.mark.slow
    @pytest.mark.filterwarnings('ignore::RuntimeWarning')
    def test_vector_valued_statistic_gh17715(self):
        # gh-17715 reported a mistake introduced in the extension of BCa to
        # multi-sample statistics; a `len` should have been `.shape[-1]`. Check
        # that this is resolved.
        # If this is resolved for NumPy, it's resolved for the rest,
        # let's only run this for NumPy.

        rng = np.random.default_rng(141921000979291141)

        def concordance(x, y, axis):
            xm = x.mean(axis)
            ym = y.mean(axis)
            cov = ((x - xm[..., None]) * (y - ym[..., None])).mean(axis)
            return (2 * cov) / (x.var(axis) + y.var(axis) + (xm - ym) ** 2)

        def statistic(tp, tn, fp, fn, axis):
            actual = tp + fp
            expected = tp + fn
            return np.nan_to_num(concordance(actual, expected, axis))

        def statistic_extradim(*args, axis):
            return statistic(*args, axis)[np.newaxis, ...]

        data = [[4, 0, 0, 2],  # (tp, tn, fp, fn)
                [2, 1, 2, 1],
                [0, 6, 0, 0],
                [0, 6, 3, 0],
                [0, 8, 1, 0]]
        data = np.array(data).T

        res = bootstrap(data, statistic_extradim, rng=rng, paired=True)
        ref = bootstrap(data, statistic, rng=rng, paired=True)
        assert_allclose(res.confidence_interval.low[0],
                        ref.confidence_interval.low, atol=1e-15)
        assert_allclose(res.confidence_interval.high[0],
                        ref.confidence_interval.high, atol=1e-15)

    # torch doesn't have T distribution CDF and can't fall back to NumPy on GPU
    @pytest.mark.skip_xp_backends(np_only=True, exceptions=['cupy'])
    def test_gh_20850(self, xp):
        rng = np.random.default_rng(2085020850)
        x = xp.asarray(rng.random((10, 2)))
        y = xp.asarray(rng.random((11, 2)))
        def statistic(x, y, axis):
            return stats.ttest_ind(x, y, axis=axis).statistic

        # The shapes do *not* need to be the same along axis
        stats.bootstrap((x, y), statistic)
        stats.bootstrap((x.T, y.T), statistic, axis=1)
        # But even when the shapes *are* the same along axis, the lengths
        # along other dimensions have to be the same (or `bootstrap` warns).
        message = "Array shapes are incompatible for broadcasting."
        with pytest.raises(ValueError, match=message):
            stats.bootstrap((x, y[:10, 0]), statistic)  # this won't work after 1.16
        stats.bootstrap((x, y[:10, 0:1]), statistic)  # this will
        stats.bootstrap((x.T, y.T[0:1, :10]), statistic, axis=1)  # this will

# --- Test Monte Carlo Hypothesis Test --- #

@make_xp_test_case(monte_carlo_test)
class TestMonteCarloHypothesisTest:
    atol = 2.5e-2  # for comparing p-value

    def get_rvs(self, rvs_in, rs, dtype=None, xp=np):
        return lambda *args, **kwds: xp.asarray(rvs_in(*args, random_state=rs, **kwds),
                                                dtype=dtype)

    def get_statistic(self, xp):
        def statistic(x, axis):
            m = xp.mean(x, axis=axis)
            v = xp.var(x, axis=axis, correction=1)
            n = x.shape[axis]
            return m / (v/n)**0.5
            # return stats.ttest_1samp(x, popmean=0., axis=axis).statistic)
        return statistic

    def test_input_validation(self, xp):
        # test that the appropriate error messages are raised for invalid input

        data = xp.asarray([1., 2., 3.])
        def stat(x, axis=None):
            return xp.mean(x, axis=axis)

        message = "Array shapes are incompatible for broadcasting."
        temp = (xp.zeros((2, 5)), xp.zeros((3, 5)))
        rvs = (stats.norm.rvs, stats.norm.rvs)
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(temp, rvs, lambda x, y, axis: 1, axis=-1)

        message = "`axis` must be an integer."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, axis=1.5)

        message = "`vectorized` must be `True`, `False`, or `None`."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, vectorized=1.5)

        message = "`rvs` must be callable or sequence of callables."
        with pytest.raises(TypeError, match=message):
            monte_carlo_test(data, None, stat)
        with pytest.raises(TypeError, match=message):
            temp = xp.asarray([[1., 2.], [3., 4.]])
            monte_carlo_test(temp, [lambda x: x, None], stat)

        message = "If `rvs` is a sequence..."
        with pytest.raises(ValueError, match=message):
            temp = xp.asarray([[1., 2., 3.]])
            monte_carlo_test(temp, [lambda x: x, lambda x: x], stat)

        message = "`statistic` must be callable."
        with pytest.raises(TypeError, match=message):
            monte_carlo_test(data, stats.norm.rvs, None)

        message = "`n_resamples` must be a positive integer."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, n_resamples=-1000)

        message = "`n_resamples` must be a positive integer."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, n_resamples=1000.5)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, batch=-1000)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, batch=1000.5)

        message = "`alternative` must be in..."
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(data, stats.norm.rvs, stat, alternative='ekki')

        # *If* this raises a value error, make sure it has the intended message
        message = "Signature inspection of statistic"
        def rvs(size):
            return xp.asarray(stats.norm.rvs(size=size))
        try:
            monte_carlo_test(data, rvs, xp.mean)
        except ValueError as e:
            assert str(e).startswith(message)

    def test_input_validation_xp(self, xp):
        def non_vectorized_statistic(x):
            return xp.mean(x)

        message = "`statistic` must be vectorized..."
        sample = xp.asarray([1., 2., 3.])
        if is_numpy(xp):
            monte_carlo_test(sample, stats.norm.rvs, non_vectorized_statistic)
            return

        with pytest.raises(ValueError, match=message):
            monte_carlo_test(sample, stats.norm.rvs, non_vectorized_statistic)
        with pytest.raises(ValueError, match=message):
            monte_carlo_test(sample, stats.norm.rvs, xp.mean, vectorized=False)

    @pytest.mark.xslow
    def test_batch(self, xp):
        # make sure that the `batch` parameter is respected by checking the
        # maximum batch size provided in calls to `statistic`
        rng = np.random.default_rng(23492340193)
        x = xp.asarray(rng.standard_normal(size=10))

        def statistic(x, axis):
            batch_size = 1 if x.ndim == 1 else x.shape[0]
            statistic.batch_size = max(batch_size, statistic.batch_size)
            statistic.counter += 1
            return self.get_statistic(xp)(x, axis=axis)
        statistic.counter = 0
        statistic.batch_size = 0

        kwds = {'sample': x, 'statistic': statistic,
                'n_resamples': 1000, 'vectorized': True}

        kwds['rvs'] = self.get_rvs(stats.norm.rvs, np.random.default_rng(328423), xp=xp)
        res1 = monte_carlo_test(batch=1, **kwds)
        assert_equal(statistic.counter, 1001)
        assert_equal(statistic.batch_size, 1)

        kwds['rvs'] = self.get_rvs(stats.norm.rvs, np.random.default_rng(328423), xp=xp)
        statistic.counter = 0
        res2 = monte_carlo_test(batch=50, **kwds)
        assert_equal(statistic.counter, 21)
        assert_equal(statistic.batch_size, 50)

        kwds['rvs'] = self.get_rvs(stats.norm.rvs, np.random.default_rng(328423), xp=xp)
        statistic.counter = 0
        res3 = monte_carlo_test(**kwds)
        assert_equal(statistic.counter, 2)
        assert_equal(statistic.batch_size, 1000)

        xp_assert_equal(res1.pvalue, res3.pvalue)
        xp_assert_equal(res2.pvalue, res3.pvalue)

    @pytest.mark.parametrize('axis', range(-3, 3))
    def test_axis_dtype(self, axis, xp):
        # test that Nd-array samples are handled correctly for valid values
        # of the `axis` parameter; also make sure non-default dtype is maintained
        rng = np.random.default_rng(2389234)
        size = [2, 3, 4]
        size[axis] = 100

        # Determine non-default dtype
        dtype_default = xp.asarray(1.).dtype
        dtype_str = 'float32'if ("64" in str(dtype_default)) else 'float64'
        dtype_np = getattr(np, dtype_str)
        dtype = getattr(xp, dtype_str)

        # ttest_1samp is CPU array-API compatible, but it would be good to
        # include CuPy in this test. We'll perform ttest_1samp with a
        # NumPy array, but all the rest with be done with fully array-API
        # compatible code.
        x = rng.standard_normal(size=size, dtype=dtype_np)
        expected = stats.ttest_1samp(x, popmean=0., axis=axis)

        x = xp.asarray(x, dtype=dtype)
        statistic = self.get_statistic(xp)
        rvs = self.get_rvs(stats.norm.rvs, rng, dtype=dtype, xp=xp)

        res = monte_carlo_test(x, rvs, statistic, vectorized=True,
                               n_resamples=20000, axis=axis)

        ref_statistic = xp.asarray(expected.statistic, dtype=dtype)
        ref_pvalue = xp.asarray(expected.pvalue, dtype=dtype)
        xp_assert_close(res.statistic, ref_statistic)
        xp_assert_close(res.pvalue, ref_pvalue, atol=self.atol)

    @pytest.mark.parametrize('alternative', ("two-sided", "less", "greater"))
    def test_alternative(self, alternative, xp):
        # test that `alternative` is working as expected
        rng = np.random.default_rng(65723433)

        x = rng.standard_normal(size=30)
        ref = stats.ttest_1samp(x, 0., alternative=alternative)

        x = xp.asarray(x)
        statistic = self.get_statistic(xp)
        rvs = self.get_rvs(stats.norm.rvs, rng, xp=xp)

        res = monte_carlo_test(x, rvs, statistic, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(ref.statistic))
        xp_assert_close(res.pvalue, xp.asarray(ref.pvalue), atol=self.atol)


    # Tests below involve statistics that are not yet array-API compatible.
    # They can be converted when the statistics are converted.
    @pytest.mark.slow
    @pytest.mark.parametrize('alternative', ("less", "greater"))
    @pytest.mark.parametrize('a', np.linspace(-0.5, 0.5, 5))  # skewness
    def test_against_ks_1samp(self, alternative, a):
        # test that monte_carlo_test can reproduce pvalue of ks_1samp
        rng = np.random.default_rng(65723433)

        x = stats.skewnorm.rvs(a=a, size=30, random_state=rng)
        expected = stats.ks_1samp(x, stats.norm.cdf, alternative=alternative)

        def statistic1d(x):
            return stats.ks_1samp(x, stats.norm.cdf, mode='asymp',
                                  alternative=alternative).statistic

        norm_rvs = self.get_rvs(stats.norm.rvs, rng)
        res = monte_carlo_test(x, norm_rvs, statistic1d,
                               n_resamples=1000, vectorized=False,
                               alternative=alternative)

        assert_allclose(res.statistic, expected.statistic)
        if alternative == 'greater':
            assert_allclose(res.pvalue, expected.pvalue, atol=self.atol)
        elif alternative == 'less':
            assert_allclose(1-res.pvalue, expected.pvalue, atol=self.atol)

    @pytest.mark.parametrize('hypotest', (stats.skewtest, stats.kurtosistest))
    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    @pytest.mark.parametrize('a', np.linspace(-2, 2, 5))  # skewness
    def test_against_normality_tests(self, hypotest, alternative, a):
        # test that monte_carlo_test can reproduce pvalue of normality tests
        rng = np.random.default_rng(85723405)

        x = stats.skewnorm.rvs(a=a, size=150, random_state=rng)
        expected = hypotest(x, alternative=alternative)

        def statistic(x, axis):
            return hypotest(x, axis=axis).statistic

        norm_rvs = self.get_rvs(stats.norm.rvs, rng)
        res = monte_carlo_test(x, norm_rvs, statistic, vectorized=True,
                               alternative=alternative)

        assert_allclose(res.statistic, expected.statistic)
        assert_allclose(res.pvalue, expected.pvalue, atol=self.atol)

    @pytest.mark.parametrize('a', np.arange(-2, 3))  # skewness parameter
    def test_against_normaltest(self, a):
        # test that monte_carlo_test can reproduce pvalue of normaltest
        rng = np.random.default_rng(12340513)

        x = stats.skewnorm.rvs(a=a, size=150, random_state=rng)
        expected = stats.normaltest(x)

        def statistic(x, axis):
            return stats.normaltest(x, axis=axis).statistic

        norm_rvs = self.get_rvs(stats.norm.rvs, rng)
        res = monte_carlo_test(x, norm_rvs, statistic, vectorized=True,
                               alternative='greater')

        assert_allclose(res.statistic, expected.statistic)
        assert_allclose(res.pvalue, expected.pvalue, atol=self.atol)

    @pytest.mark.xslow
    @pytest.mark.parametrize('a', np.linspace(-0.5, 0.5, 5))  # skewness
    def test_against_cramervonmises(self, a):
        # test that monte_carlo_test can reproduce pvalue of cramervonmises
        rng = np.random.default_rng(234874135)

        x = stats.skewnorm.rvs(a=a, size=30, random_state=rng)
        expected = stats.cramervonmises(x, stats.norm.cdf)

        def statistic1d(x):
            return stats.cramervonmises(x, stats.norm.cdf).statistic

        norm_rvs = self.get_rvs(stats.norm.rvs, rng)
        res = monte_carlo_test(x, norm_rvs, statistic1d,
                               n_resamples=1000, vectorized=False,
                               alternative='greater')

        assert_allclose(res.statistic, expected.statistic)
        assert_allclose(res.pvalue, expected.pvalue, atol=self.atol)

    @pytest.mark.slow
    @pytest.mark.parametrize('dist_name', ['norm', 'logistic'])
    @pytest.mark.parametrize('target_statistic', [0.6, 0.7, 0.8])
    def test_against_anderson(self, dist_name, target_statistic):
        # test that monte_carlo_test can reproduce results of `anderson`.

        # find the skewness for which the sample statistic is within the range of
        # values tubulated by `anderson`

        def fun(a):
            rng = np.random.default_rng(394295467)
            x = stats.tukeylambda.rvs(a, size=100, random_state=rng)
            expected = stats.anderson(x, dist_name, method='interpolate')
            return expected.statistic - target_statistic
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            sol = root(fun, x0=0)
        assert sol.success

        # get the significance level (p-value) associated with that critical
        # value
        a = sol.x[0]
        rng = np.random.default_rng(394295467)
        x = stats.tukeylambda.rvs(a, size=100, random_state=rng)
        expected = stats.anderson(x, dist_name, method='interpolate')
        expected_stat = expected.statistic
        expected_p = expected.pvalue

        # perform equivalent Monte Carlo test and compare results
        def statistic1d(x):
            return stats.anderson(x, dist_name, method='interpolate').statistic

        dist_rvs = self.get_rvs(getattr(stats, dist_name).rvs, rng)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            res = monte_carlo_test(x, dist_rvs,
                                   statistic1d, n_resamples=1000,
                                   vectorized=False, alternative='greater')

        assert_allclose(res.statistic, expected_stat)
        assert_allclose(res.pvalue, expected_p, atol=2*self.atol)

    def test_p_never_zero(self):
        # Use biased estimate of p-value to ensure that p-value is never zero
        # per monte_carlo_test reference [1]
        rng = np.random.default_rng(2190176673029737545)
        x = np.zeros(100)
        res = monte_carlo_test(x, rng.random, np.mean,
                               vectorized=True, alternative='less')
        assert res.pvalue == 0.0001

    def test_against_ttest_ind(self):
        # test that `monte_carlo_test` can reproduce results of `ttest_ind`.
        rng = np.random.default_rng(219017667302737545)
        data = rng.random(size=(2, 5)), rng.random(size=7)  # broadcastable
        rvs = rng.normal, rng.normal
        def statistic(x, y, axis):
            return stats.ttest_ind(x, y, axis=axis).statistic

        res = stats.monte_carlo_test(data, rvs, statistic, axis=-1)
        ref = stats.ttest_ind(data[0], [data[1]], axis=-1)
        assert_allclose(res.statistic, ref.statistic)
        assert_allclose(res.pvalue, ref.pvalue, rtol=2e-2)

    def test_against_f_oneway(self):
        # test that `monte_carlo_test` can reproduce results of `f_oneway`.
        rng = np.random.default_rng(219017667302737545)
        data = (rng.random(size=(2, 100)), rng.random(size=(2, 101)),
                rng.random(size=(2, 102)), rng.random(size=(2, 103)))
        rvs = rng.normal, rng.normal, rng.normal, rng.normal

        def statistic(*args, axis):
            return stats.f_oneway(*args, axis=axis).statistic

        res = stats.monte_carlo_test(data, rvs, statistic, axis=-1,
                                     alternative='greater')
        ref = stats.f_oneway(*data, axis=-1)

        assert_allclose(res.statistic, ref.statistic)
        assert_allclose(res.pvalue, ref.pvalue, atol=1e-2)

    @pytest.mark.fail_slow(2)
    @pytest.mark.xfail_on_32bit("Statistic may not depend on sample order on 32-bit")
    def test_finite_precision_statistic(self):
        # Some statistics return numerically distinct values when the values
        # should be equal in theory. Test that `monte_carlo_test` accounts
        # for this in some way.
        rng = np.random.default_rng(2549824598234528)
        n_resamples = 9999
        def rvs(size):
            return 1. * stats.bernoulli(p=0.333).rvs(size=size, random_state=rng)

        x = rvs(100)
        res = stats.monte_carlo_test(x, rvs, np.var, alternative='less',
                                     n_resamples=n_resamples)
        # show that having a tolerance matters
        c0 = np.sum(res.null_distribution <= res.statistic)
        c1 = np.sum(res.null_distribution <= res.statistic*(1+1e-15))
        assert c0 != c1
        assert res.pvalue == (c1 + 1)/(n_resamples + 1)


@make_xp_test_case(power)
class TestPower:
    def xp_normal(self, rng, *, xp, dtype=None):
        dtype = xp_default_dtype(xp) if dtype is None else dtype
        return lambda size: xp.asarray(rng.normal(size=size), dtype=dtype)

    def test_input_validation(self, xp):
        # test that the appropriate error messages are raised for invalid input
        rng = np.random.default_rng(8519895914314711673)

        test = stats.ttest_ind
        rvs = (self.xp_normal(rng, xp=xp), self.xp_normal(rng, xp=xp))
        n_observations = (xp.asarray(10), xp.asarray(12))

        message = "`vectorized` must be `True`, `False`, or `None`."
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, vectorized=1.5)

        message = "`rvs` must be callable or sequence of callables."
        with pytest.raises(TypeError, match=message):
            power(test, None, n_observations)
        with pytest.raises(TypeError, match=message):
            power(test, (self.xp_normal(rng, xp=xp), 'ekki'), n_observations)

        message = "If `rvs` is a sequence..."
        with pytest.raises(ValueError, match=message):
            power(test, (self.xp_normal(rng, xp=xp),), n_observations)
        with pytest.raises(ValueError, match=message):
            power(test, rvs, (10,))

        message = "`significance` must contain floats between 0 and 1."
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, significance=2)
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, significance=xp.linspace(-1, 1, 50))

        message = "`kwargs` must be a dictionary"
        with pytest.raises(TypeError, match=message):
            power(test, rvs, n_observations, kwargs=(1, 2, 3))

        message = "not be broadcast|Chunks do not add|Incompatible shapes"
        with pytest.raises((ValueError, RuntimeError), match=message):
            power(test, rvs, (xp.asarray([10, 11]), xp.asarray([12, 13, 14])))
        with pytest.raises((ValueError, RuntimeError), match=message):
            power(test, rvs, (xp.asarray([10, 11]), xp.asarray([12, 13])),
                  kwargs={'x': xp.asarray([1, 2, 3])})

        message = "`test` must be callable"
        with pytest.raises(TypeError, match=message):
            power(None, rvs, n_observations)

        message = "`n_resamples` must be a positive integer"
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, n_resamples=-10)
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, n_resamples=10.5)

        message = "`batch` must be a positive integer"
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, batch=-10)
        with pytest.raises(ValueError, match=message):
            power(test, rvs, n_observations, batch=10.5)

    @pytest.mark.slow
    def test_batch(self, xp):
        # make sure that the `batch` parameter is respected by checking the
        # maximum batch size provided in calls to `test`
        rng = np.random.default_rng(23492340193)

        def test(x, axis):
            batch_size = 1 if x.ndim == 1 else x.shape[0]
            test.batch_size = max(batch_size, test.batch_size)
            test.counter += 1
            return stats.ttest_1samp(x, xp.asarray(0.), axis=axis).pvalue
        test.counter = 0
        test.batch_size = 0

        kwds = dict(test=test,
                    n_observations=xp.asarray(10),
                    n_resamples=1000)

        rng = np.random.default_rng(23492340193)
        res1 = power(**kwds, rvs=self.xp_normal(rng, xp=xp), batch=1)
        assert test.counter == 1000
        assert test.batch_size == 1

        rng = np.random.default_rng(23492340193)
        test.counter = 0
        res2 = power(**kwds, rvs=self.xp_normal(rng, xp=xp), batch=50)
        assert test.counter == 20
        assert test.batch_size == 50

        rng = np.random.default_rng(23492340193)
        test.counter = 0
        res3 = power(**kwds, rvs=self.xp_normal(rng, xp=xp), batch=1000)
        assert test.counter == 1
        assert test.batch_size == 1000

        xp_assert_equal(res1.power, res3.power)
        xp_assert_equal(res2.power, res3.power)

    @pytest.mark.slow
    def test_vectorization(self, xp):
        # Test that `power` is vectorized as expected
        rng = np.random.default_rng(25495254834552)
        alternatives = {-1: 'less', 0:'two-sided', 1: 'greater'}

        # Single vectorized call
        popmeans = xp.asarray([0, 0.2])
        def test(x, alternative, axis=-1):
            # ensure that popmeans axis is zeroth and orthogonal to the rest
            popmeans_expanded = xpx.expand_dims(popmeans,
                                                axis=tuple(range(1, x.ndim + 1)))
            alternative = alternatives[int(alternative)]
            return stats.ttest_1samp(x, popmeans_expanded, alternative=alternative,
                                     axis=axis)

        # nx and kwargs broadcast against one another
        nx = xp.asarray([10, 15, 20, 50, 100])[:, xp.newaxis]
        kwargs = {'alternative': xp.asarray([-1, 0, 1])}

        # This dimension is added to the beginning
        significance = xp.asarray([0.01, 0.025, 0.05, 0.1])
        res = stats.power(test, self.xp_normal(rng, xp=xp), nx,
                          significance=significance, kwargs=kwargs)

        # Looping over all combinations
        ref = []
        for significance_i in significance:
            for i in range(nx.shape[0]):
                nx_i = nx[i, ...]
                for alternative_i in kwargs['alternative']:
                    for j in range(popmeans.shape[0]):
                        popmean_j = popmeans[j]
                        def test2(x, axis=-1):
                            return stats.ttest_1samp(x, popmean_j, axis=axis,
                                alternative=alternatives[int(alternative_i)])

                        tmp = stats.power(test2, self.xp_normal(rng, xp=xp), nx_i,
                                          significance=significance_i)
                        ref.append(tmp.power)
        ref = xp.reshape(xp.stack(ref), res.power.shape)

        # Show that results are similar
        xp_assert_close(res.power, ref, rtol=2e-2, atol=1e-2)

    @pytest.mark.skip_xp_backends(cpu_only=True,
                                  exceptions=['cupy', 'jax.numpy', 'dask.array'])
    def test_ttest_ind_null(self, xp):
        # Check that the p-values of `ttest_ind` are uniformly distributed under
        # the null hypothesis
        rng = np.random.default_rng(254952548345528)

        test = stats.ttest_ind
        n_observations = (xp.asarray(rng.integers(10, 100, size=(10))),
                          xp.asarray(rng.integers(10, 100, size=(10))))
        rvs = self.xp_normal(rng, xp=xp), self.xp_normal(rng, xp=xp)
        significance = xp.asarray([0.01, 0.05, 0.1])
        res = stats.power(test, rvs, n_observations, significance=significance)
        significance = xp.broadcast_to(significance[:, xp.newaxis], res.power.shape)
        xp_assert_close(res.power, significance, atol=1e-2)

    @pytest.mark.skip_xp_backends('array_api_strict',
                                  reason='currently combines integer and float arrays')
    @pytest.mark.skip_xp_backends(cpu_only=True,
                                  exceptions=['cupy', 'jax.numpy', 'dask.array'])
    def test_ttest_1samp_power(self, xp):
        # Check simulated ttest_1samp power against reference
        rng = np.random.default_rng(254952548345528)

        # Reference values computed with statmodels
        # import numpy as np
        # from statsmodels.stats.power import tt_solve_power
        # tt_solve_power = np.vectorize(tt_solve_power)
        # tt_solve_power([0.1, 0.5, 0.9], [[10], [20]], [[[0.01]], [[0.05]]])
        ref = [[[0.0126515 , 0.10269751, 0.40415802],
                [0.01657775, 0.29734608, 0.86228288]],
               [[0.0592903 , 0.29317561, 0.71718121],
                [0.07094116, 0.56450441, 0.96815163]]]

        kwargs = {'popmean': xp.asarray([0.1, 0.5, 0.9])}
        n_observations = xp.asarray([[10], [20]])
        significance = xp.asarray([0.01, 0.05])
        res = stats.power(stats.ttest_1samp, self.xp_normal(rng, xp=xp), n_observations,
                          significance=significance, kwargs=kwargs)
        xp_assert_close(res.power, xp.asarray(ref), atol=1e-2)


@make_xp_test_case(permutation_test)
class TestPermutationTest:

    rtol = 1e-14

    def setup_method(self):
        self.rng = np.random.default_rng(7170559330470561044)

    # -- Input validation -- #
    def test_permutation_test_iv(self, xp):

        def stat(x, y, axis):
            return stats.ttest_ind((x, y), axis).statistic

        data = (xp.asarray([1, 2, 3]), xp.asarray([1, 2, 3]))

        message = "each sample in `data` must contain two or more ..."
        with pytest.raises(ValueError, match=message):
            permutation_test((data[0], xp.asarray([0])), stat)

        message = "`data` must be a tuple containing at least two samples"
        with pytest.raises(ValueError, match=message):
            permutation_test((1,), stat)
        with pytest.raises(TypeError, match=message):
            permutation_test(1, stat)

        message = "`axis` must be an integer."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, axis=1.5)

        message = "`permutation_type` must be in..."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat,
                             permutation_type="ekki")

        message = "`vectorized` must be `True`, `False`, or `None`."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, vectorized=1.5)

        message = "`n_resamples` must be a positive integer."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, n_resamples=-1000)

        message = "`n_resamples` must be a positive integer."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, n_resamples=1000.5)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, batch=-1000)

        message = "`batch` must be a positive integer or None."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, batch=1000.5)

        message = "`alternative` must be in..."
        with pytest.raises(ValueError, match=message):
            permutation_test(data, stat, alternative='ekki')

        message = "SeedSequence expects int or sequence of ints"
        with pytest.raises(TypeError, match=message):
            permutation_test(data, stat, rng='herring')

    # -- Test Parameters -- #
    # SPEC-007 leave one call with seed to check it still works
    @pytest.mark.parametrize('random_state', [np.random.RandomState,
                                              np.random.default_rng])
    @pytest.mark.parametrize('permutation_type',
                             ['pairings', 'samples', 'independent'])
    def test_batch(self, permutation_type, random_state, xp):
        # make sure that the `batch` parameter is respected by checking the
        # maximum batch size provided in calls to `statistic`
        x = xp.asarray(self.rng.random(10))
        y = xp.asarray(self.rng.random(10))

        def statistic(x, y, axis):
            batch_size = 1 if x.ndim == 1 else x.shape[0]
            statistic.batch_size = max(batch_size, statistic.batch_size)
            statistic.counter += 1
            return xp.mean(x, axis=axis) - xp.mean(y, axis=axis)
        statistic.counter = 0
        statistic.batch_size = 0

        kwds = {'n_resamples': 1000, 'permutation_type': permutation_type,
                'vectorized': True}
        res1 = stats.permutation_test((x, y), statistic, batch=1,
                                      random_state=random_state(0), **kwds)
        assert statistic.counter == 1001
        assert statistic.batch_size == 1

        statistic.counter = 0
        res2 = stats.permutation_test((x, y), statistic, batch=50,
                                      random_state=random_state(0), **kwds)
        assert statistic.counter == 21
        assert statistic.batch_size == 50

        statistic.counter = 0
        res3 = stats.permutation_test((x, y), statistic, batch=1000,
                                      random_state=random_state(0), **kwds)
        assert statistic.counter == 2
        assert statistic.batch_size == 1000

        xp_assert_equal(res1.pvalue, res3.pvalue)
        xp_assert_equal(res2.pvalue, res3.pvalue)

    # SPEC-007 leave at least one call with seed to check it still works
    @pytest.mark.parametrize('random_state', [np.random.RandomState,
                                              np.random.default_rng])
    @pytest.mark.parametrize('permutation_type, exact_size',
                             [('pairings', special.factorial(3)**2),
                              ('samples', 2**3),
                              ('independent', special.binom(6, 3))])
    def test_permutations(self, permutation_type, exact_size, random_state, xp):
        # make sure that the `permutations` parameter is respected by checking
        # the size of the null distribution
        x = xp.asarray(self.rng.random(3))
        y = xp.asarray(self.rng.random(3))

        def statistic(x, y, axis):
            return xp.mean(x, axis=axis) - xp.mean(y, axis=axis)

        kwds = {'permutation_type': permutation_type,
                'vectorized': True}
        res = stats.permutation_test((x, y), statistic, n_resamples=3,
                                     random_state=random_state(0), **kwds)
        assert xp_size(res.null_distribution) == 3

        res = stats.permutation_test((x, y), statistic, **kwds)
        assert xp_size(res.null_distribution) == exact_size

    # -- Randomized Permutation Tests -- #

    # To get reasonable accuracy, these next three tests are somewhat slow.
    # Originally, I had them passing for all combinations of permutation type,
    # alternative, and RNG, but that takes too long for CI. Instead, split
    # into three tests, each testing a particular combination of the three
    # parameters.

    def test_randomized_test_against_exact_both(self, xp):
        # check that the randomized and exact tests agree to reasonable
        # precision for permutation_type='both

        alternative, rng = 'less', 0

        nx, ny, permutations = 8, 9, 24000
        assert special.binom(nx + ny, nx) > permutations

        rng = np.random.default_rng(8235259808)
        x = xp.asarray(rng.standard_normal(size=nx))
        y = xp.asarray(rng.standard_normal(size=ny))
        data = x, y

        def statistic(x, y, axis):
            return xp.mean(x, axis=axis) - xp.mean(y, axis=axis)

        kwds = {'vectorized': True, 'permutation_type': 'independent',
                'batch': 100, 'alternative': alternative, 'rng': rng}
        res = permutation_test(data, statistic, n_resamples=permutations, **kwds)
        res2 = permutation_test(data, statistic, n_resamples=xp.inf, **kwds)

        assert res.statistic == res2.statistic
        xp_assert_close(res.pvalue, res2.pvalue, atol=1e-2)

    @pytest.mark.slow()
    def test_randomized_test_against_exact_samples(self, xp):
        # check that the randomized and exact tests agree to reasonable
        # precision for permutation_type='samples'

        alternative, rng = 'greater', None

        nx, ny, permutations = 15, 15, 32000
        assert 2**nx > permutations

        rng = np.random.default_rng(8235259808)
        x = xp.asarray(rng.standard_normal(size=nx))
        y = xp.asarray(rng.standard_normal(size=ny))
        data = x, y

        def statistic(x, y, axis):
            return xp.mean(x - y, axis=axis)

        kwds = {'vectorized': True, 'permutation_type': 'samples',
                'batch': 100, 'alternative': alternative, 'rng': rng}
        res = permutation_test(data, statistic, n_resamples=permutations, **kwds)
        res2 = permutation_test(data, statistic, n_resamples=xp.inf, **kwds)

        assert res.statistic == res2.statistic
        xp_assert_close(res.pvalue, res2.pvalue, atol=1e-2)

    # I only need to skip torch on GPU because it doesn't have betaincc for pearsonr
    @pytest.mark.skip_xp_backends(cpu_only=True, exceptions=['cupy', 'jax.numpy'])
    def test_randomized_test_against_exact_pairings(self, xp):
        # check that the randomized and exact tests agree to reasonable
        # precision for permutation_type='pairings'

        alternative, rng = 'two-sided', self.rng

        nx, ny, permutations = 8, 8, 40000
        assert special.factorial(nx) > permutations

        rng = np.random.default_rng(8235259808)
        x = xp.asarray(rng.standard_normal(size=nx))
        y = xp.asarray(rng.standard_normal(size=ny))
        data = [x]

        def statistic(x, axis):
            return stats.pearsonr(x, y, axis=axis).statistic

        kwds = {'vectorized': True, 'permutation_type': 'samples',
                'batch': 100, 'alternative': alternative, 'rng': rng}
        res = permutation_test(data, statistic, n_resamples=permutations, **kwds)
        res2 = permutation_test(data, statistic, n_resamples=xp.inf, **kwds)

        assert res.statistic == res2.statistic
        xp_assert_close(res.pvalue, res2.pvalue, atol=1e-2)

    # -- Independent (Unpaired) Sample Tests -- #

    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    def test_against_ks_2samp(self, alternative, xp):
        x = self.rng.normal(size=4, scale=1)
        y = self.rng.normal(size=5, loc=3, scale=3)

        expected = stats.ks_2samp(x, y, alternative=alternative, mode='exact')

        def statistic(x, y, axis):
            # todo: use `xp` as backend when `ks_2samp` is translated to array API
            x, y = _xp_copy_to_numpy(x), _xp_copy_to_numpy(y)
            res = stats.ks_2samp(x, y, axis=axis, mode='asymp', alternative=alternative)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        # ks_2samp is always a one-tailed 'greater' test
        # it's the statistic that changes (D+ vs D- vs max(D+, D-))
        x, y = xp.asarray(x), xp.asarray(y)
        res = permutation_test((x, y), statistic, n_resamples=np.inf,
                               alternative='greater', rng=self.rng)

        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    def test_against_ansari(self, alternative, xp):
        x = self.rng.normal(size=4, scale=1)
        y = self.rng.normal(size=5, scale=3)

        # ansari has a different convention for 'alternative'
        alternative_correspondence = {"less": "greater",
                                      "greater": "less",
                                      "two-sided": "two-sided"}
        alternative_scipy = alternative_correspondence[alternative]
        expected = stats.ansari(x, y, alternative=alternative_scipy)

        def statistic(x, y, axis):
            # todo: use `xp` as backend when `ansari` is translated to array API
            x, y = _xp_copy_to_numpy(x), _xp_copy_to_numpy(y)
            res = stats.ansari(x, y, axis=axis)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        x, y = xp.asarray(x), xp.asarray(y)
        res = permutation_test((x, y), statistic, n_resamples=np.inf,
                               alternative=alternative, rng=self.rng)

        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    @skip_xp_backends('cupy', reason='needs mannwhitneyu')
    @skip_xp_backends(eager_only=True)  # mannwhitneyu does input validation
    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    def test_against_mannwhitneyu(self, alternative, xp):
        x = stats.uniform.rvs(size=(3, 5, 2), loc=0, random_state=self.rng)
        y = stats.uniform.rvs(size=(3, 5, 2), loc=0.05, random_state=self.rng)

        expected = stats.mannwhitneyu(x, y, axis=1, alternative=alternative)

        def statistic(x, y, axis):
            return stats.mannwhitneyu(x, y, axis=axis, method='asymptotic').statistic

        x, y = xp.asarray(x), xp.asarray(y)
        res = permutation_test((x, y), statistic, vectorized=True,
                               n_resamples=xp.inf, alternative=alternative,
                               axis=1, rng=self.rng)

        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    @skip_xp_backends('cupy', reason='needs cramervonmises_2samp')
    @skip_xp_backends(eager_only=True)  # cramervonmises_2samp does input validation
    @skip_xp_backends(cpu_only=True)  # torch doesn't have `kv`
    def test_against_cvm(self, xp):
        x = stats.norm.rvs(size=4, scale=1, random_state=self.rng)
        y = stats.norm.rvs(size=5, loc=3, scale=3, random_state=self.rng)

        expected = stats.cramervonmises_2samp(x, y, method='exact')

        def statistic(x, y, axis):
            res = stats.cramervonmises_2samp(x, y, axis=axis, method='asymptotic')
            return res.statistic

        # cramervonmises_2samp has only one alternative, greater
        x, y = xp.asarray(x), xp.asarray(y)
        res = permutation_test((x, y), statistic, n_resamples=np.inf,
                               alternative='greater', rng=self.rng)

        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    @skip_xp_backends('cupy', reason='needs kruskal')
    @skip_xp_backends(eager_only=True)  # kruskal does input validation
    @pytest.mark.parametrize('axis', (-1, 2))
    def test_vectorized_nsamp_ptype_both(self, axis, xp):
        # statistic only available for NumPy

        # Test that permutation_test with permutation_type='independent' works
        # properly for a 3-sample statistic with nd array samples of different
        # (but compatible) shapes and ndims. Show that exact permutation test
        # and random permutation tests approximate SciPy's asymptotic pvalues
        # and that exact and random permutation test results are even closer
        # to one another (than they are to the asymptotic results).

        # Three samples, different (but compatible) shapes with different ndims
        rng = np.random.default_rng(6709265303529651545)
        x = rng.random(size=(3))
        y = rng.random(size=(1, 3, 2))
        z = rng.random(size=(2, 1, 4))
        data = (x, y, z)

        expected = stats.kruskal(*data, axis=axis)

        # Define the statistic (and pvalue for comparison)
        def statistic(*data, axis):
            return stats.kruskal(*data, axis=axis).statistic

        # Calculate exact and randomized permutation results
        kwds = {'axis': axis, 'alternative': 'greater',
                'permutation_type': 'independent', 'rng': rng}
        data = [xp.asarray(data_) for data_ in data]
        res = permutation_test(data, statistic, n_resamples=xp.inf, **kwds)
        res2 = permutation_test(data, statistic, n_resamples=1000, **kwds)

        # Check results
        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol*5)
        xp_assert_close(res.statistic, res2.statistic, rtol=self.rtol*5)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), atol=6e-2)
        xp_assert_close(res.pvalue, res2.pvalue, atol=3e-2)

    # -- Paired-Sample Tests -- #

    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    def test_against_wilcoxon(self, alternative, xp):
        x = stats.uniform.rvs(size=(3, 6, 2), loc=0, random_state=self.rng)
        y = stats.uniform.rvs(size=(3, 6, 2), loc=0.05, random_state=self.rng)
        expected = stats.wilcoxon(x, y, alternative=alternative, axis=1)

        # We'll check both 1- and 2-sample versions of the same test;
        # we expect identical results to wilcoxon in all cases.
        def statistic_1samp_1d(z, axis):
            # todo: use `xp` as backend when `wilcoxon` is translated to array API
            # 'less' ensures we get the same of two statistics every time
            z = _xp_copy_to_numpy(z)
            res = stats.wilcoxon(z, alternative='less', axis=axis)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        def statistic_2samp_1d(x, y, axis):
            # todo: use `xp` as backend when `wilcoxon` is translated to array API
            x, y = _xp_copy_to_numpy(x), _xp_copy_to_numpy(y)
            res = stats.wilcoxon(x, y, alternative='less', axis=axis)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        x, y = xp.asarray(x), xp.asarray(y)
        kwds = {'axis': 1, 'alternative': alternative, 'permutation_type': 'samples',
                'rng': self.rng, 'n_resamples': np.inf}
        res1 = permutation_test((x-y,), statistic_1samp_1d, **kwds)
        res2 = permutation_test((x, y), statistic_2samp_1d, **kwds)

        # `wilcoxon` returns a different statistic with 'two-sided'
        xp_assert_close(res1.statistic, res2.statistic, rtol=self.rtol)
        if alternative != 'two-sided':
            xp_assert_close(res2.statistic, xp.asarray(expected.statistic),
                            rtol=self.rtol)

        xp_assert_close(res2.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)
        xp_assert_close(res1.pvalue, res2.pvalue, rtol=self.rtol)

    @pytest.mark.parametrize('alternative', ("less", "greater", "two-sided"))
    def test_against_binomtest(self, alternative, xp):

        x = self.rng.integers(0, 2, size=10)
        x[x == 0] = -1
        # More naturally, the test would flip elements between 0 and one.
        # However, permutation_test will flip the _signs_ of the elements.
        # So we have to work with +1/-1 instead of 1/0.

        def statistic(x, axis=0):
            xp_ = array_namespace(x)
            return xp_.count_nonzero(x > 0, axis=axis)

        k, n, p = statistic(x), 10, 0.5
        expected = stats.binomtest(k, n, p, alternative=alternative)

        res = stats.permutation_test((xp.asarray(x, dtype=xp.float64),), statistic,
                                     vectorized=True,
                                     permutation_type='samples', n_resamples=xp.inf,
                                     rng=self.rng, alternative=alternative)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    # -- Exact Association Tests -- #

    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    def test_against_kendalltau(self, xp):
        x = self.rng.normal(size=6)
        y = x + self.rng.normal(size=6)

        expected = stats.kendalltau(x, y, method='exact')

        def statistic(x, axis):
            # todo: use `xp` as backend when `kendalltau` is translated to array API
            x = _xp_copy_to_numpy(x)
            res = stats.kendalltau(x, y, method='asymptotic', axis=axis)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        # kendalltau currently has only one alternative, two-sided
        x = xp.asarray(x)
        res = permutation_test((x,), statistic, permutation_type='pairings',
                               n_resamples=np.inf, rng=self.rng)

        xp_assert_close(res.statistic, xp.asarray(expected.statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected.pvalue), rtol=self.rtol)

    @pytest.mark.parametrize('alternative', ('less', 'greater', 'two-sided'))
    def test_against_fisher_exact(self, alternative, xp):
        # x and y are binary random variables with some dependence
        rng = np.random.default_rng(6235696159000529929)
        x = (rng.random(7) > 0.6).astype(float)
        y = (rng.random(7) + 0.25*x > 0.6).astype(float)
        tab = stats.contingency.crosstab(x, y)[1]

        x, y = xp.asarray(x), xp.asarray(y)
        def statistic(x, axis):
            return xp.count_nonzero((x == 1) & (y == 1), axis=axis)

        res = permutation_test((x,), statistic, permutation_type='pairings',
                               n_resamples=xp.inf, alternative=alternative,
                               rng=rng)
        res2 = stats.fisher_exact(tab, alternative=alternative)

        xp_assert_close(res.pvalue, xp.asarray(res2.pvalue, dtype=x.dtype))

    @pytest.mark.xslow()
    @pytest.mark.parametrize('axis', (-2, 1))
    def test_vectorized_nsamp_ptype_samples(self, axis):
        # statistic only available for NumPy, and it's a pain to vectorize

        # Test that permutation_test with permutation_type='samples' works
        # properly for a 3-sample statistic with nd array samples of different
        # (but compatible) shapes and ndims. Show that exact permutation test
        # reproduces SciPy's exact pvalue and that random permutation test
        # approximates it.

        x = self.rng.random(size=(2, 4, 3))
        y = self.rng.random(size=(1, 4, 3))
        z = self.rng.random(size=(2, 4, 1))
        x = stats.rankdata(x, axis=axis)
        y = stats.rankdata(y, axis=axis)
        z = stats.rankdata(z, axis=axis)
        y = y[0]  # to check broadcast with different ndim
        data = (x, y, z)

        def statistic1d(*data):
            return stats.page_trend_test(data, ranked=True,
                                         method='asymptotic').statistic

        def pvalue1d(*data):
            return stats.page_trend_test(data, ranked=True,
                                         method='exact').pvalue

        statistic = _resampling._vectorize_statistic(statistic1d)
        pvalue = _resampling._vectorize_statistic(pvalue1d)

        expected_statistic = statistic(*np.broadcast_arrays(*data), axis=axis)
        expected_pvalue = pvalue(*np.broadcast_arrays(*data), axis=axis)

        # Let's forgive this use of an integer seed, please.
        kwds = {'vectorized': False, 'axis': axis, 'alternative': 'greater',
                'permutation_type': 'pairings', 'rng': 0}
        res = permutation_test(data, statistic1d, n_resamples=np.inf, **kwds)
        res2 = permutation_test(data, statistic1d, n_resamples=5000, **kwds)

        assert_allclose(res.statistic, expected_statistic, rtol=self.rtol)
        assert_allclose(res.statistic, res2.statistic, rtol=self.rtol)
        assert_allclose(res.pvalue, expected_pvalue, rtol=self.rtol)
        assert_allclose(res.pvalue, res2.pvalue, atol=3e-2)

    # -- Test Against External References -- #

    tie_case_1 = {'x': [1, 2, 3, 4], 'y': [1.5, 2, 2.5],
                  'expected_less': 0.2000000000,
                  'expected_2sided': 0.4,  # 2*expected_less
                  'expected_Pr_gte_S_mean': 0.3428571429,  # see note below
                  'expected_statistic': 7.5,
                  'expected_avg': 9.142857, 'expected_std': 1.40698}
    tie_case_2 = {'x': [111, 107, 100, 99, 102, 106, 109, 108],
                  'y': [107, 108, 106, 98, 105, 103, 110, 105, 104],
                  'expected_less': 0.1555738379,
                  'expected_2sided': 0.3111476758,
                  'expected_Pr_gte_S_mean': 0.2969971205,  # see note below
                  'expected_statistic': 32.5,
                  'expected_avg': 38.117647, 'expected_std': 5.172124}

    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    @pytest.mark.xslow()  # only the second case is slow, really
    @pytest.mark.parametrize('case', (tie_case_1, tie_case_2))
    def test_with_ties(self, case, xp):
        """
        Results above from SAS PROC NPAR1WAY, e.g.

        DATA myData;
        INPUT X Y;
        CARDS;
        1 1
        1 2
        1 3
        1 4
        2 1.5
        2 2
        2 2.5
        ods graphics on;
        proc npar1way AB data=myData;
            class X;
            EXACT;
        run;
        ods graphics off;

        Note: SAS provides Pr >= |S-Mean|, which is different from our
        definition of a two-sided p-value.

        """

        x = case['x']
        y = case['y']

        expected_statistic = xp.asarray(case['expected_statistic'])
        expected_less = xp.asarray(case['expected_less'])
        expected_2sided = xp.asarray(case['expected_2sided'])
        expected_Pr_gte_S_mean = xp.asarray(case['expected_Pr_gte_S_mean'])
        expected_avg = xp.asarray(case['expected_avg'])
        expected_std = xp.asarray(case['expected_std'])

        def statistic(x, y, axis):
            # todo: use `xp` as backend when `ansari` is translated to array API
            x, y = _xp_copy_to_numpy(x), _xp_copy_to_numpy(y)
            res = stats.ansari(x, y, axis=axis)
            res = xp.asarray(res.statistic)
            return res[()] if res.ndim == 0 else res

        dtype = xp_default_dtype(xp)
        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "Ties preclude use of exact statistic", UserWarning)
            res = permutation_test((x, y), statistic, n_resamples=np.inf,
                                   alternative='less')
            res2 = permutation_test((x, y), statistic, n_resamples=np.inf,
                                    alternative='two-sided')

        xp_assert_close(res.statistic, expected_statistic, rtol=self.rtol)
        xp_assert_close(res.pvalue, expected_less, atol=1e-10)
        xp_assert_close(res2.pvalue, expected_2sided, atol=1e-10)
        xp_assert_close(xp.mean(res2.null_distribution), expected_avg, rtol=1e-6)
        xp_assert_close(xp.std(res2.null_distribution), expected_std, rtol=1e-6)

        # SAS provides Pr >= |S-Mean|; might as well check against that, too
        S = res.statistic
        mean = xp.mean(res.null_distribution)
        n = res.null_distribution.shape[0]
        Pr_gte_S_mean = xp.astype(xp.count_nonzero(
            xp.abs(res.null_distribution-mean) >= xp.abs(S-mean)), S.dtype) / n
        xp_assert_close(Pr_gte_S_mean, expected_Pr_gte_S_mean)

    @pytest.mark.slow
    @pytest.mark.parametrize('alternative, expected_pvalue',
                             (('less', 0.9708333333333),
                              ('greater', 0.05138888888889),
                              ('two-sided', 0.1027777777778)))
    # I only need to skip torch on GPU because it doesn't have betaincc for pearsonr
    @pytest.mark.skip_xp_backends(cpu_only=True, exceptions=['cupy', 'jax.numpy'])
    @pytest.mark.skip_xp_backends(eager_only=True)  # TODO: change to jax_jit=False
    def test_against_spearmanr_in_R(self, alternative, expected_pvalue, xp):
        """
        Results above from R cor.test, e.g.

        options(digits=16)
        x <- c(1.76405235, 0.40015721, 0.97873798,
               2.2408932, 1.86755799, -0.97727788)
        y <- c(2.71414076, 0.2488, 0.87551913,
               2.6514917, 2.01160156, 0.47699563)
        cor.test(x, y, method = "spearm", alternative = "t")
        """
        # data comes from
        # np.random.seed(0)
        # x = stats.norm.rvs(size=6)
        # y = x + stats.norm.rvs(size=6)
        x = xp.asarray([1.76405235, 0.40015721, 0.97873798,
                        2.2408932, 1.86755799, -0.97727788])
        y = xp.asarray([2.71414076, 0.2488, 0.87551913,
                        2.6514917, 2.01160156, 0.47699563])
        expected_statistic = 0.7714285714285715

        y = xp.asarray(stats.rankdata(_xp_copy_to_numpy(y)))
        def statistic(x, axis):
            # `spearmanr` is not array api compatible, but `pearsonr` is. So for now
            # use _xp_copy_to_numpy just for ranking so we can run this test w/ CuPy.
            # TODO: use `xp` as backend when cupy works with `rankdata`
            x = xp.asarray(stats.rankdata(_xp_copy_to_numpy(x), axis=axis))
            return stats.pearsonr(x, y, axis=axis).statistic

        res = permutation_test((x,), statistic, permutation_type='pairings',
                               n_resamples=xp.inf, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(expected_statistic), rtol=self.rtol)
        xp_assert_close(res.pvalue, xp.asarray(expected_pvalue), atol=1e-13)

    @pytest.mark.parametrize("batch", (-1, 0))
    def test_batch_generator_iv(self, batch):
        with pytest.raises(ValueError, match="`batch` must be positive."):
            list(_resampling._batch_generator([1, 2, 3], batch))

    batch_generator_cases = [(range(0), 3, []),
                             (range(6), 3, [[0, 1, 2], [3, 4, 5]]),
                             (range(8), 3, [[0, 1, 2], [3, 4, 5], [6, 7]])]

    @pytest.mark.parametrize("iterable, batch, expected",
                             batch_generator_cases)
    def test_batch_generator(self, iterable, batch, expected):
        got = list(_resampling._batch_generator(iterable, batch))
        assert got == expected

    @pytest.mark.fail_slow(2)
    # I only need to skip torch on GPU because it doesn't have betaincc for pearsonr
    @pytest.mark.skip_xp_backends(cpu_only=True, exceptions=['cupy', 'jax.numpy'])
    def test_finite_precision_statistic(self, xp):
        # Some statistics return numerically distinct values when the values
        # should be equal in theory. Test that `permutation_test` accounts
        # for this in some way.
        x = xp.asarray([1., 2., 4., 3.], dtype=xp.float64)
        y = xp.asarray([2., 4., 6., 8.], dtype=xp.float64)

        def statistic(x, y, axis):
            return stats.pearsonr(x, y, axis=axis)[0]

        res = stats.permutation_test((x, y), statistic,
                                     permutation_type='pairings')
        r, pvalue, null = res.statistic, res.pvalue, res.null_distribution

        correct_p = 2 * float(xp.count_nonzero(null >= r - 1e-14)) / null.shape[0]
        assert pvalue == correct_p == 1/3
        # Compare against other exact correlation tests using R corr.test
        # options(digits=16)
        # x = c(1, 2, 4, 3)
        # y = c(2, 4, 6, 8)
        # cor.test(x, y, alternative = "t", method = "spearman")  # 0.333333333
        # cor.test(x, y, alternative = "t", method = "kendall")  # 0.333333333


def test_all_partitions_concatenated():
    # make sure that _all_paritions_concatenated produces the correct number
    # of partitions of the data into samples of the given sizes and that
    # all are unique
    n = np.array([3, 2, 4], dtype=int)
    nc = np.cumsum(n)

    all_partitions = set()
    counter = 0
    for partition_concatenated in _resampling._all_partitions_concatenated(n):
        counter += 1
        partitioning = np.split(partition_concatenated, nc[:-1])
        all_partitions.add(tuple([frozenset(i) for i in partitioning]))

    expected = np.prod([special.binom(sum(n[i:]), sum(n[i+1:]))
                        for i in range(len(n)-1)])

    assert_equal(counter, expected)
    assert_equal(len(all_partitions), expected)


@pytest.mark.parametrize('fun_name',
                         ['bootstrap', 'permutation_test', 'monte_carlo_test'])
def test_parameter_vectorized(fun_name):
    # Check that parameter `vectorized` is working as desired for all
    # resampling functions. Results don't matter; just don't fail asserts.
    rng = np.random.default_rng(75245098234592)
    sample = rng.random(size=10)

    def rvs(size):  # needed by `monte_carlo_test`
        return stats.norm.rvs(size=size, random_state=rng)

    fun_options = {'bootstrap': {'data': (sample,), 'rng': rng,
                                 'method': 'percentile'},
                   'permutation_test': {'data': (sample,), 'rng': rng,
                                        'permutation_type': 'samples'},
                   'monte_carlo_test': {'sample': sample, 'rvs': rvs}}
    common_options = {'n_resamples': 100}

    fun = getattr(stats, fun_name)
    options = fun_options[fun_name]
    options.update(common_options)

    def statistic(x, axis):
        assert x.ndim > 1 or np.array_equal(x, sample)
        return np.mean(x, axis=axis)
    fun(statistic=statistic, vectorized=None, **options)
    fun(statistic=statistic, vectorized=True, **options)

    def statistic(x):
        assert x.ndim == 1
        return np.mean(x)
    fun(statistic=statistic, vectorized=None, **options)
    fun(statistic=statistic, vectorized=False, **options)


class TestMonteCarloMethod:
    def test_rvs_and_random_state(self):
        message = "Use of `rvs` and `rng` are mutually exclusive."
        rng = np.random.default_rng(34982345)
        with pytest.raises(ValueError, match=message):
            stats.MonteCarloMethod(rvs=rng.random, rng=rng)
