import numpy as np
import scipy.sparse as sp

from sklearn.datasets import make_regression
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KRR, KernelRidge
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.testing import assert_array_almost_equal


X, y = make_regression(n_features=10, random_state=0)
Xcsr = sp.csr_matrix(X)
Xcsc = sp.csc_matrix(X)
Y = np.array([y, y]).T


def test_check_estimator():
    check_estimator(KernelRidge)
    check_estimator(KRR)


def test_kernel_ridge():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, y).predict(X)
    pred2 = KRR(kernel="linear", alpha=1).fit(X, y).predict(X)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csr():
    pred = Ridge(alpha=1, fit_intercept=False,
                 solver="cholesky").fit(Xcsr, y).predict(Xcsr)
    pred2 = KRR(kernel="linear", alpha=1).fit(Xcsr, y).predict(Xcsr)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csc():
    pred = Ridge(alpha=1, fit_intercept=False,
                 solver="cholesky").fit(Xcsc, y).predict(Xcsc)
    pred2 = KRR(kernel="linear", alpha=1).fit(Xcsc, y).predict(Xcsc)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_singular_kernel():
    # alpha=0 causes a LinAlgError in computing the dual coefficients,
    # which causes a fallback to a lstsq solver. This is tested here.
    pred = Ridge(alpha=0, fit_intercept=False).fit(X, y).predict(X)
    kr = KRR(kernel="linear", alpha=0)
    ignore_warnings(kr.fit)(X, y)
    pred2 = kr.predict(X)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed():
    for kernel in ["linear", "rbf", "poly", "cosine"]:
        K = pairwise_kernels(X, X, metric=kernel)
        pred = KRR(kernel=kernel).fit(X, y).predict(X)
        pred2 = KRR(kernel="precomputed").fit(K, y).predict(K)
        assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed_kernel_unchanged():
    K = np.dot(X, X.T)
    K2 = K.copy()
    KRR(kernel="precomputed").fit(K, y)
    assert_array_almost_equal(K, K2)


def test_kernel_ridge_sample_weights():
    K = np.dot(X, X.T)  # precomputed kernel
    sw = np.random.RandomState(0).rand(X.shape[0])

    pred = Ridge(alpha=1,
                 fit_intercept=False).fit(X, y, sample_weight=sw).predict(X)
    pred2 = KRR(kernel="linear",
                        alpha=1).fit(X, y, sample_weight=sw).predict(X)
    pred3 = KRR(kernel="precomputed",
                        alpha=1).fit(K, y, sample_weight=sw).predict(K)
    assert_array_almost_equal(pred, pred2)
    assert_array_almost_equal(pred, pred3)


def test_kernel_ridge_multi_output():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, Y).predict(X)
    pred2 = KRR(kernel="linear", alpha=1).fit(X, Y).predict(X)
    assert_array_almost_equal(pred, pred2)

    pred3 = KRR(kernel="linear", alpha=1).fit(X, y).predict(X)
    pred3 = np.array([pred3, pred3]).T
    assert_array_almost_equal(pred2, pred3)
