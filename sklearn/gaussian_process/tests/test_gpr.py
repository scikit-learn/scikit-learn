"""Testing for Gaussian process regression """

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

import numpy as np

from scipy.optimize import approx_fprime

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

from sklearn.utils.testing \
    import (assert_true, assert_greater, assert_array_less,
            assert_almost_equal, assert_equal)


def f(x):
    return x * np.sin(x)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = f(X).ravel()


kernels = [RBF(l=1.0), RBF(l=1.0, l_bounds=(1e-3, 1e3)),
           (1e-2, 1.0, 1e2)*RBF(l=1.0, l_bounds=(1e-3, 1e3)),
           (1e-2, 1.0, 1e2)*RBF(l=1.0, l_bounds=(1e-3, 1e3)) + (0.0, 0.0, 1e2),
           (1e-2, 0.1, 1e2)*RBF(l=1.0, l_bounds=(1e-3, 1e3)) + (0.0, 0.0, 1e2)]


def test_gpr_interpolation():
    """Test the interpolating property for different kernels."""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)
        y_pred, y_cov = gpr.predict(X, return_cov=True)

        assert_true(np.allclose(y_pred, y))
        assert_true(np.allclose(np.diag(y_cov), 0.))


def test_lml_improving():
    """ Test that hyperparameter-tuning improves log-marginal likelihood. """
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)
        assert_greater(gpr.log_marginal_likelihood(gpr.kernel_.theta),
                       gpr.log_marginal_likelihood(kernel.theta))


def test_converged_to_local_maximum():
    """ Test that we are in local maximum after hyperparameter-optimization."""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)

        lml, lml_gradient = \
            gpr.log_marginal_likelihood(gpr.kernel_.theta, True)

        assert_true(np.all((np.abs(lml_gradient) < 1e-5)
                           | (gpr.kernel_.theta == gpr.kernel_.bounds[:, 0])
                           | (gpr.kernel_.theta == gpr.kernel_.bounds[:, 1])))


def test_solution_inside_bounds():
    """ Test that hyperparameter-optimization remains in bounds"""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)

        bounds = gpr.kernel_.bounds
        max_ = np.finfo(gpr.kernel_.theta.dtype).max
        tiny = 1e-10
        bounds[~np.isfinite(bounds[:, 1]), 1] = max_

        assert_array_less(bounds[:, 0], gpr.kernel_.theta + tiny)
        assert_array_less(gpr.kernel_.theta, bounds[:, 1] + tiny)


def test_lml_gradient():
    """ Compare analytic and numeric gradient of log marginal likelihood. """
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)

        lml, lml_gradient = gpr.log_marginal_likelihood(kernel.theta, True)
        lml_gradient_approx = \
            approx_fprime(kernel.theta,
                          lambda theta: gpr.log_marginal_likelihood(theta,
                                                                    False),
                          1e-10)

        assert_almost_equal(lml_gradient, lml_gradient_approx, 3)


def test_prior():
    """ Test that GP prior has mean 0 and identical variances."""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel)

        y_mean, y_cov = gpr.predict(X, return_cov=True)

        assert_almost_equal(y_mean, 0, 5)
        if len(gpr.kernel.theta) > 1:
            # XXX: quite hacky, works only for current kernels
            assert_almost_equal(np.diag(y_cov), kernel.theta[0] , 5)
        else:
            assert_almost_equal(np.diag(y_cov), 1, 5)


def test_sample_statistics():
    """ Test that statistics of samples drawn from GP are correct."""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)

        y_mean, y_cov = gpr.predict(X2, return_cov=True)

        samples = gpr.sample_y(X2, 1000000)

        # More digits accuracy would require many more samples
        assert_almost_equal(y_mean, np.mean(samples, 1), 2)
        assert_almost_equal(np.diag(y_cov) / np.diag(y_cov).max(),
                            np.var(samples, 1) / np.diag(y_cov).max(), 1)


def test_no_optimizer():
    """ Test that kernel parameters are unmodified when optimizer is None."""
    kernel = RBF(1.0)
    gpr = GaussianProcessRegressor(kernel=kernel, optimizer=None).fit(X, y)
    assert_equal(gpr.kernel_.theta, 1.0)
    assert_equal(gpr.theta_, 1.0)


def test_predict_cov_vs_std():
    """ Test that predicted std.-dev. is consistent with cov's diagonal."""
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)
        print gpr.kernel_
        y_mean, y_cov = gpr.predict(X2, return_cov=True)
        y_mean, y_std = gpr.predict(X2, return_std=True)
        assert_almost_equal(np.sqrt(np.diag(y_cov)), y_std)
