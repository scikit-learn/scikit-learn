"""Testing for Gaussian process classification """

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

from copy import deepcopy

import numpy as np

from scipy.optimize import approx_fprime

from sklearn.gaussian_process import GaussianProcessClassification
from sklearn.gaussian_process.kernels import RBF

from sklearn.utils.testing import (assert_true, assert_greater,
                                   assert_almost_equal, assert_array_equal)


f = lambda x: x * np.sin(x)
X = np.atleast_2d(np.linspace(0, 10, 30)).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = np.array(f(X).ravel() > 0, dtype=int)


kernels = [RBF(0.1), RBF((1e-3, 1.0, 1e3)),
           (1e-2, 1.0, 1e2) * RBF((1e-3, 0.1, 1e3))]


def test_predict_consistent():
    """ Check binary predict decision has also predicted probability above 0.5.
    """
    for kernel in kernels:
        kernel = deepcopy(kernel)
        gpc = GaussianProcessClassification(kernel=kernel).fit(X, y)
        assert_array_equal(gpc.predict(X),
                           gpc.predict_proba(X) >=0.5)


def test_lml_improving():
    """ Test that hyperparameter-tuning improves log-marginal likelihood. """
    for kernel in kernels:
        if not kernel.has_bounds:
            continue
        kernel = deepcopy(kernel)
        params_initial = kernel.params
        gpc = GaussianProcessClassification(kernel=kernel).fit(X, y)
        assert_greater(gpc.log_marginal_likelihood(kernel.params),
                       gpc.log_marginal_likelihood(params_initial))


def test_converged_to_local_maximum():
    """ Test that we are in local maximum after hyperparameter-optimization. """
    for kernel in kernels:
        if not kernel.has_bounds:
            continue
        kernel = deepcopy(kernel)
        gpc = GaussianProcessClassification(kernel=kernel).fit(X, y)

        lml, lml_gradient = gpc.log_marginal_likelihood(kernel.params, True)

        assert_almost_equal(lml_gradient, 0, 2)


def test_lml_gradient():
    """ Compare analytic and numeric gradient of log marginal likelihood. """
    for kernel in kernels:
        kernel = deepcopy(kernel)
        params = kernel.params
        gpc = GaussianProcessClassification(kernel=kernel).fit(X, y)

        lml, lml_gradient = gpc.log_marginal_likelihood(params, True)
        lml_gradient_approx = \
            approx_fprime(params,
                          lambda theta: gpc.log_marginal_likelihood(theta,
                                                                    False),
                          1e-10)

        assert_almost_equal(lml_gradient, lml_gradient_approx, 3)
