import numpy as np

from sklearn.datasets import make_classification
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics.pairwise import pairwise_kernels

from sklearn.utils.testing import assert_array_almost_equal


X, y = make_classification(n_classes=2, random_state=0)
Y = np.array([y, y]).T


def test_kernel_ridge():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, y).predict(X)
    pred2 = KernelRidge(kernel="linear", alpha=1).fit(X, y).predict(X)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_singula():
    # alpha=0 causes a LinAlgError in computing the dual coeeficients,
    # which causes a fall back to a lstsq solver. This is tested here.
    # XXX: Other solvers than lsqr give different results for Ridge
    pred = Ridge(alpha=0, fit_intercept=False, solver="lsqr").fit(X, y).predict(X)
    pred2 = KernelRidge(kernel="linear", alpha=0).fit(X, y).predict(X)
    # XXX: we are only equal within the first two decimals
    assert_array_almost_equal(pred, pred2, decimal=2)


def test_kernel_ridge_precomputed():
    for kernel in ["linear", "rbf", "poly", "cosine"]:
        K = pairwise_kernels(X, X, metric=kernel)
        pred = KernelRidge(kernel=kernel).fit(X, y).predict(X)
        pred2 = KernelRidge(kernel="precomputed").fit(K, y).predict(K)
        assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed_kernel_unchanged():
    K = np.dot(X, X.T)
    K2 = K.copy()
    KernelRidge(kernel="precomputed").fit(K, y)
    assert_array_almost_equal(K, K2)


def test_kernel_ridge_sample_weights():
    K = np.dot(X, X.T)  # precomputed kernel
    sw = np.random.RandomState(0).rand(X.shape[0])

    pred = Ridge(alpha=1,
                 fit_intercept=False).fit(X, y, sample_weight=sw).predict(X)
    pred2 = KernelRidge(kernel="linear",
                        alpha=1).fit(X, y, sample_weight=sw).predict(X)
    pred3 = KernelRidge(kernel="precomputed",
                        alpha=1).fit(K, y, sample_weight=sw).predict(K)
    assert_array_almost_equal(pred, pred2)
    assert_array_almost_equal(pred, pred3)


def test_kernel_ridge_multi_output():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, Y).predict(X)
    pred2 = KernelRidge(kernel="linear", alpha=1).fit(X, Y).predict(X)
    assert_array_almost_equal(pred, pred2)

    pred3 = KernelRidge(kernel="linear", alpha=1).fit(X, y).predict(X)
    pred3 = np.array([pred3, pred3]).T
    assert_array_almost_equal(pred2, pred3)
