"""
Testing for Kernel Ridge Classifier
"""
import pytest

import numpy as np
import scipy.sparse as sp

from sklearn import datasets
from sklearn.linear_model import RidgeClassifier
from sklearn.kernel_ridge.krc import KernelRidgeClassifier
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.utils._testing import ignore_warnings
from sklearn.utils import check_random_state

from sklearn.utils._testing import assert_array_almost_equal

# mypy error: Module 'sklearn.svm' has no attribute '_libsvm'
from sklearn.svm import _libsvm  # type: ignore

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

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

    pred = RidgeClassifier(alpha=1, fit_intercept=False).fit(X, y, sample_weight=sw).predict(X)
    pred2 = KernelRidgeClassifier(kernel="linear", alpha=1).fit(X, y, sample_weight=sw).predict(X)
    pred3 = (
        KernelRidgeClassifier(kernel="precomputed", alpha=1)
        .fit(K, y, sample_weight=sw)
        .predict(K)
    )
    assert_array_almost_equal(pred, pred2)
    assert_array_almost_equal(pred, pred3)


# TODO: Remove in 1.1
def test_kernel_ridge_pairwise_is_deprecated():
    k_ridge = KernelRidgeClassifier(kernel="precomputed")
    msg = r"Attribute `_pairwise` was deprecated in version 0\.24"
    with pytest.warns(FutureWarning, match=msg):
        k_ridge._pairwise
