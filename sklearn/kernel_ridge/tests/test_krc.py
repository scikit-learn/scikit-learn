"""
Testing for Kernel Ridge Classifier
"""
import numpy as np
import scipy.sparse as sp

from sklearn import datasets
from sklearn.linear_model import RidgeClassifier
from sklearn.kernel_ridge.krc import KernelRidgeClassifier
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.utils import check_random_state

from sklearn.utils._testing import assert_array_almost_equal


# also load the iris dataset
iris = datasets.load_iris()
rng = check_random_state(42)
perm = rng.permutation(iris.target.size)
X = iris.data[perm]
y = iris.target[perm]
Xcsr = sp.csr_matrix(X)
Xcsc = sp.csc_matrix(X)


def test_kernel_RidgeClassifier():
    pred = RidgeClassifier(alpha=1, fit_intercept=False).fit(X, y).predict(X)
    pred2 = KernelRidgeClassifier(kernel="linear", alpha=1).fit(X, y).predict(X)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csr():
    pred = (
        RidgeClassifier(alpha=1, fit_intercept=False, solver="cholesky")
        .fit(Xcsr, y)
        .predict(Xcsr)
    )
    pred2 = KernelRidgeClassifier(kernel="linear", alpha=1).fit(Xcsr, y).predict(Xcsr)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csc():
    pred = (
        RidgeClassifier(alpha=1, fit_intercept=False, solver="cholesky")
        .fit(Xcsc, y)
        .predict(Xcsc)
    )
    pred2 = KernelRidgeClassifier(kernel="linear", alpha=1).fit(Xcsc, y).predict(Xcsc)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed():
    for kernel in ["linear", "rbf", "poly", "cosine"]:
        K = pairwise_kernels(X, X, metric=kernel)
        pred = KernelRidgeClassifier(kernel=kernel).fit(X, y).predict(X)
        pred2 = KernelRidgeClassifier(kernel="precomputed").fit(K, y).predict(K)
        assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_precomputed_kernel_unchanged():
    K = np.dot(X, X.T)
    K2 = K.copy()
    KernelRidgeClassifier(kernel="precomputed").fit(K, y)
    assert_array_almost_equal(K, K2)


def test_kernel_ridge_sample_weights():
    K = np.dot(X, X.T)  # precomputed kernel
    sw = np.random.RandomState(0).rand(X.shape[0])

    pred = (
        RidgeClassifier(alpha=1, fit_intercept=False)
        .fit(X, y, sample_weight=sw)
        .predict(X)
    )
    pred2 = (
        KernelRidgeClassifier(kernel="linear", alpha=1)
        .fit(X, y, sample_weight=sw)
        .predict(X)
    )
    pred3 = (
        KernelRidgeClassifier(kernel="precomputed", alpha=1)
        .fit(K, y, sample_weight=sw)
        .predict(K)
    )
    assert_array_almost_equal(pred, pred2)
    assert_array_almost_equal(pred, pred3)
