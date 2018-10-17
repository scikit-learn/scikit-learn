import numpy as np

from sklearn.datasets import make_classification
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.fast_kernel import FKC_EigenPro

np.random.seed(1)
X, y = make_classification()
X2, y2 = make_classification(n_samples=300, n_features=500, n_informative=170)
X3, y3 = make_classification(n_redundant=18)
X4, y4 = make_classification(shift=1, hypercube=False)

X5, y5 = make_classification(n_features=200, n_repeated=50)
X6, y6 = np.concatenate([X, X]), np.concatenate([y, 1-y])


def test_fast_kernel_classification_gaussian():
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_laplace():
    FKC_prediction = FKC_EigenPro(
        kernel="laplace", n_epoch=50, bandwidth=13).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_regression_cauchy():
    FKR_prediction = FKC_EigenPro(
        kernel="cauchy", n_epoch=50, bandwidth=10, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(FKR_prediction, y)


def test_fast_kernel_classification_complex():
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X2, y2).predict(X2)
    assert_array_almost_equal(FKC_prediction, y2)


def test_fast_kernel_classification_redundant():
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=.01, n_epoch=50).fit(X3, y3).predict(X3)
    assert_array_almost_equal(FKC_prediction, y3)


def test_fast_kernel_classification_shift():
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X4, y4).predict(X4)
    assert_array_almost_equal(FKC_prediction, y4)

def test_fast_kernel_classification_duplicate_data():
    FKR_prediction = FKC_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=.1, random_state=0).fit(X5, y5).predict(X5)
    assert_array_almost_equal(FKR_prediction, y5)


def test_fast_kernel_classification_conflict_data():
    try:
        FKR_prediction = FKC_EigenPro(
            kernel="gaussian", n_epoch=5, bandwidth=5, random_state=0).fit(X6, y6).predict(X6)
        assert_array_almost_equal(FKR_prediction, y6)
        raise AssertionError("Predicted impossibly well")
    except AssertionError:
        return
