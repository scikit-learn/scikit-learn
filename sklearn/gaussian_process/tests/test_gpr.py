"""Testing for Gaussian process regression """

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

from copy import deepcopy

import numpy as np

from scipy.optimize import approx_fprime

from sklearn.gaussian_process import GaussianProcessRegression
from sklearn.gaussian_process.kernels import RBF

from sklearn.utils.testing import (assert_true, assert_greater,
                                   assert_almost_equal)


f = lambda x: x * np.sin(x)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = f(X).ravel()


kernels = [RBF(0.1), RBF((1e-3, 1.0, 1e3)),
           (1e-2, 1.0, 1e2) * RBF((1e-3, 0.1, 1e3)),
           (1e-2, 1.0, 1e2) * RBF((1e-3, 0.1, 1e3)) + (None, 0.0, None),
           (1e-2, 0.1, 1e2) * RBF((1e-3, 0.1, 1e3)) + (None, 0.0, None)]

def test_gpr_interpolation():
    """Test the interpolating property for different kernels."""
    for kernel in kernels:
        kernel = deepcopy(kernel)
        gpr = GaussianProcessRegression(kernel=kernel).fit(X, y)
        y_pred, y_cov = gpr.predict(X, return_cov=True)

        assert_true(np.allclose(y_pred, y))
        assert_true(np.allclose(np.diag(y_cov), 0.))


def test_lml_improving():
    """ Test that hyperparameter-tuning improves log-marginal likelihood. """
    for kernel in kernels:
        if not kernel.has_bounds:
            continue
        kernel = deepcopy(kernel)
        params_initial = kernel.params
        gpr = GaussianProcessRegression(kernel=kernel).fit(X, y)
        assert_greater(gpr.log_marginal_likelihood(kernel.params),
                       gpr.log_marginal_likelihood(params_initial))


def test_converged_to_local_maximum():
    """ Test that we are in local maximum after hyperparameter-optimization. """
    for kernel in kernels:
        if not kernel.has_bounds:
            continue
        kernel = deepcopy(kernel)
        gpr = GaussianProcessRegression(kernel=kernel).fit(X, y)

        lml, lml_gradient = gpr.log_marginal_likelihood(kernel.params, True)

        assert_almost_equal(lml_gradient, 0, 5)


def test_lml_gradient():
    """ Compare analytic and numeric gradient of log marginal likelihood. """
    for kernel in kernels:
        kernel = deepcopy(kernel)
        params = kernel.params
        gpr = GaussianProcessRegression(kernel=kernel).fit(X, y)

        lml, lml_gradient = gpr.log_marginal_likelihood(params, True)
        lml_gradient_approx = \
            approx_fprime(params,
                          lambda theta: gpr.log_marginal_likelihood(theta,
                                                                    False),
                          1e-10)

        assert_almost_equal(lml_gradient, lml_gradient_approx, 3)


def test_prior():
    """ Test that GP prior has mean 0 and identical variances."""
    for kernel in kernels:
        kernel = deepcopy(kernel)
        gpr = GaussianProcessRegression(kernel=kernel)

        y_mean, y_cov = gpr.predict(X, return_cov=True)

        assert_almost_equal(y_mean, 0, 5)
        if len(kernel.params) > 1:
            # XXX: quite hacky, works only for current kernels
            assert_almost_equal(np.diag(y_cov), kernel.params[0] , 5)
        else:
            assert_almost_equal(np.diag(y_cov), 1, 5)


def test_sample_statistics():
    """ Test that statistics of samples drawn from GP are correct."""
    for kernel in kernels:
        kernel = deepcopy(kernel)
        gpr = GaussianProcessRegression(kernel=kernel).fit(
            X, y)

        y_mean, y_cov = gpr.predict(X2, return_cov=True)

        samples = gpr.sample(X2, 1000000)

        # More digits accuracy would require many more samples
        assert_almost_equal(y_mean, np.mean(samples, 1), 2)
        assert_almost_equal(np.diag(y_cov) / np.diag(y_cov).max(),
                            np.var(samples, 1) / np.diag(y_cov).max(), 1)
