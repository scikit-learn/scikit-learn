import numpy as np

from sklearn.datasets import make_classification
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge

from sklearn.utils.testing import assert_array_almost_equal


X, y = make_classification(n_classes=2, random_state=0)
Y = np.array([y, y]).T


def test_kernel_ridge():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, y).predict(X)
    pred2 = KernelRidge(kernel="linear", alpha=1).fit(X, y).predict(X)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed():
    K = np.dot(X, X.T)
    pred = KernelRidge(kernel="linear").fit(X, y).predict(X)
    pred2 = KernelRidge(kernel="precomputed").fit(K, y).predict(K)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed_sample_weight():
    K = np.dot(X, X.T)
    K2 = K.copy()
    sw = np.ones(X.shape[0]) / float(X.shape[0])
    KernelRidge(kernel="precomputed").fit(K, y, sample_weight=sw)
    assert_array_almost_equal(K, K2)


def test_kernel_ridge_multi_output():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X, Y).predict(X)
    pred2 = KernelRidge(kernel="linear", alpha=1).fit(X, Y).predict(X)
    assert_array_almost_equal(pred, pred2)

    pred3 = KernelRidge(kernel="linear", alpha=1).fit(X, y).predict(X)
    pred3 = np.array([pred3, pred3]).T
    assert_array_almost_equal(pred2, pred3)
