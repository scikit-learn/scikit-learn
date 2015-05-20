"""Testing for Gaussian process regression """

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

import numpy as np

from scipy.optimize import approx_fprime

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels \
    import RBF, ConstantKernel as C, WhiteKernel

from sklearn.utils.testing \
    import (assert_true, assert_greater, assert_array_less,
            assert_almost_equal, assert_equal)


def f(x):
    return x * np.sin(x)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = f(X).ravel()


kernels = [RBF(l=1.0), RBF(l=1.0, l_bounds=(1e-3, 1e3)),
           C(1.0, (1e-2, 1e2)) * RBF(l=1.0, l_bounds=(1e-3, 1e3)),
           C(1.0, (1e-2, 1e2)) * RBF(l=1.0, l_bounds=(1e-3, 1e3))
               + C(0.0, (0.0, 1e2)),
           C(0.1, (1e-2, 1e2)) * RBF(l=1.0, l_bounds=(1e-3, 1e3))
               + C(0.0, (0.0, 1e2))]


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
        y_mean, y_cov = gpr.predict(X2, return_cov=True)
        y_mean, y_std = gpr.predict(X2, return_std=True)
        assert_almost_equal(np.sqrt(np.diag(y_cov)), y_std)


def test_anisotropic_kernel():
    """ Test that GPR can identify meaningful anisotropic length-scales. """
    # We learn a function which varies in one dimension ten-times slower
    # than in the other. The corresponding length-scales should differ by at
    # least a factor 5
    rng = np.random.RandomState(0)
    X = rng.uniform(-1, 1, (50, 2))
    y = X[:, 0] +  0.1 * X[:, 1]

    kernel = RBF([1.0, 1.0])
    gpr = GaussianProcessRegressor(kernel=kernel).fit(X, y)
    assert_greater(gpr.kernel_.theta[1], gpr.kernel_.theta[0] * 5)

def test_random_starts():
    """
    Test that an increasing number of random-starts of GP fitting only
    increases the log marginal likelihood of the chosen theta.
    """
    n_samples, n_features = 25, 3
    np.random.seed(0)
    rng = np.random.RandomState(0)
    X = rng.randn(n_samples, n_features) * 2 - 1
    y = np.sin(X).sum(axis=1) + np.sin(3 * X).sum(axis=1) \
        + rng.normal(scale=0.1, size=n_samples)

    kernel = C(1.0, (1e-2, 1e2)) \
        * RBF(l=[1.0] * n_features, l_bounds=[(1e-4, 1e+2)] * n_features) \
        + WhiteKernel(c=1e-5, c_bounds=(1e-5, 1e1))
    last_lml = -np.inf
    for n_restarts_optimizer in range(1, 10):
        gp = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=n_restarts_optimizer,
            random_state=0,).fit(X, y)
        lml = gp.log_marginal_likelihood(gp.theta_)
        assert_greater(lml, last_lml - np.finfo(np.float32).eps)
        last_lml = lml

def test_y_normalization():
    """ Test normalization of the target values in GP

    Fitting non-normalizing GP on normalized y and fitting normalizing GP
    on unnormalized y should yield identical results
    """
    y_mean, y_std = y.mean(0), y.std(0)
    y_norm = (y - y_mean) / y_std
    for kernel in kernels:
        # Fit non-normalizing GP on normalized y
        gpr = GaussianProcessRegressor(kernel=kernel,
                                       sigma_squared_n=1e-10 / y_std**2)
        gpr.fit(X, y_norm)
        # Fit normalizing GP on unnormalized y
        gpr_norm = GaussianProcessRegressor(kernel=kernel,
                                            sigma_squared_n=1e-10,
                                            normalize_y=True)
        gpr_norm.fit(X, y)

        # Compare predicted mean, std-devs and covariances
        y_pred, y_pred_std = gpr.predict(X2, return_std=True)
        y_pred = y_mean + y_pred * y_std
        y_pred_std *= y_std
        y_pred_norm, y_pred_std_norm = gpr_norm.predict(X2, return_std=True)

        assert_almost_equal(y_pred, y_pred_norm)
        assert_almost_equal(y_pred_std, y_pred_std_norm)

        _, y_cov = gpr.predict(X2, return_cov=True)
        y_cov *= y_std ** 2
        _, y_cov_norm = gpr_norm.predict(X2, return_cov=True)
        assert_almost_equal(y_cov, y_cov_norm)


def test_y_multioutput():
    """
    """
    y_2d = np.vstack((y, y*2)).T

    # Test for fixed kernel that first dimension of 2d GP equals the output
    # of 1d GP and that second dimension is twice as large
    kernel = RBF(l=1.0)

    gpr = GaussianProcessRegressor(kernel=kernel, optimizer=None,
                                   normalize_y=False)
    gpr.fit(X, y)

    gpr_2d = GaussianProcessRegressor(kernel=kernel, optimizer=None,
                                      normalize_y=False)
    gpr_2d.fit(X, y_2d)

    y_pred_1d, y_std_1d = gpr.predict(X2, return_std=True)
    y_pred_2d, y_std_2d = gpr_2d.predict(X2, return_std=True)
    _, y_cov_1d = gpr.predict(X2, return_cov=True)
    _, y_cov_2d = gpr_2d.predict(X2, return_cov=True)

    assert_almost_equal(y_pred_1d, y_pred_2d[:, 0])
    assert_almost_equal(y_pred_1d, y_pred_2d[:, 1] / 2)

    # Standard deviation and covariance do not depend on output
    assert_almost_equal(y_std_1d, y_std_2d)
    assert_almost_equal(y_cov_1d, y_cov_2d)

    y_sample_1d = gpr.sample_y(X2, n_samples=10)
    y_sample_2d = gpr_2d.sample_y(X2, n_samples=10)
    assert_almost_equal(y_sample_1d, y_sample_2d[:, 0])

    # Test hyperparameter optimization
    for kernel in kernels:
        gpr = GaussianProcessRegressor(kernel=kernel, normalize_y=True)
        gpr.fit(X, y)

        gpr_2d = GaussianProcessRegressor(kernel=kernel, normalize_y=True)
        gpr_2d.fit(X, np.vstack((y, y)).T)

        assert_almost_equal(gpr.kernel_.theta, gpr_2d.kernel_.theta, 4)
