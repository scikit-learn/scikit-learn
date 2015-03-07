"""Testing for kernels for Gaussian processes."""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

from copy import deepcopy

import numpy as np

from scipy.optimize import approx_fprime

from sklearn.metrics.pairwise import PAIRWISE_KERNEL_FUNCTIONS
from sklearn.gaussian_process.kernels \
    import (RBF, RationalQuadratic, ExpSineSquared, DotProduct,
            ConstantKernel, WhiteKernel, PairwiseKernel)

from sklearn.utils.testing import assert_equal, assert_almost_equal


X = np.random.normal(0, 1, (10, 2))

kernels = [RBF(2.0), RBF([0.5, 2.0]),
           ConstantKernel(10.0),
           2.0 * RBF(0.5), RBF(2.0) + WhiteKernel(1.0),
           RationalQuadratic([1.0, 1.0]),
           ExpSineSquared([1.0, 1.0]),
           DotProduct(1.0), DotProduct(1.0, degree=2)]
for metric in PAIRWISE_KERNEL_FUNCTIONS:
    if metric in ["additive_chi2", "chi2"]:
        continue
    kernels.append(PairwiseKernel(1.0, metric=metric))


def test_kernel_gradient():
    """ Compare analytic and numeric gradient of kernels. """
    for kernel in kernels:
        K, K_gradient = kernel(X, eval_gradient=True)

        assert_equal(K_gradient.shape[0], X.shape[0])
        assert_equal(K_gradient.shape[1], X.shape[0])
        assert_equal(K_gradient.shape[2], kernel.params.shape[0])

        K_gradient_approx = np.empty_like(K_gradient)
        for i in range(K.shape[0]):
            for j in range(K.shape[1]):
                def eval_kernel_ij_for_theta(theta):
                    kernel_copy = deepcopy(kernel)
                    kernel_copy.params = theta
                    K = kernel_copy(X, eval_gradient=False)
                    return K[i, j]
                K_gradient_approx[i, j] = \
                    approx_fprime(kernel.params, eval_kernel_ij_for_theta,
                                  1e-10)

        assert_almost_equal(K_gradient, K_gradient_approx, 4)


def test_auto_vs_cross():
    """ Auto-correlation and cross-correlation should be consistent. """
    for kernel in kernels:
        K_auto = kernel(X)
        K_cross = kernel(X, X)
        assert_almost_equal(K_auto, K_cross, 5)

def test_kernel_operator_commutative():
    """ Adding kernels and multiplying kernels should be commutative. """
    # Check addition
    assert_almost_equal((RBF(2.0) + 1.0)(X),
                        (1.0 + RBF(2.0))(X))

    # Check multiplication
    assert_almost_equal((3.0 * RBF(2.0))(X),
                        (RBF(2.0) * 3.0)(X))


def test_kernel_anisotropic():
    """ Anisotropic kernel should be consistent with isotropic kernels."""
    K = RBF([0.5, 2.0])(X)
    X1 = np.array(X)
    X1[:, 0] *= 4
    K1 = RBF(2.0)(X1)
    assert_almost_equal(K, K1)

    X2 = np.array(X)
    X2[:, 1] /= 4
    K2 = RBF(0.5)(X2)
    assert_almost_equal(K, K2)


def test_kernel_stationary():
    """ Test stationarity of kernels."""
    for kernel in kernels:
        if not kernel.is_stationary():
            continue
        K = kernel(X, X + 1)
        assert_almost_equal(K[0, 0], np.diag(K))
