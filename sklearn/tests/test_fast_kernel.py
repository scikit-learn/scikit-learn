import numpy as np

from sklearn.datasets import make_regression, make_classification
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.fast_kernel import FKR_EigenPro, FKC_EigenPro

np.random.seed(1)
# Tests for Fast Kernel Regression and Classification


def test_fast_kernel_regression_gaussian():
    X, y = make_regression(n_features=100, random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=100, bandwidth=10,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_laplace():
    X, y = make_regression(random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="laplace", n_epoch=100, bandwidth=8,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_cauchy():
    X, y = make_regression(random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="cauchy", n_epoch=100, bandwidth=10, subsample_size=1000,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_2d():
    X, y = make_regression(n_features=200, n_targets=30, random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=100, bandwidth=14,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_many_features():
    X, y = make_regression(n_features=10000, random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=100, bandwidth=1,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_simple():
    X, y = make_regression(n_features=100, n_informative=1,
                           random_state=1)
    FKR_prediction = FKR_EigenPro(
        batch_size=500, kernel="gaussian", n_epoch=100, bandwidth=10,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_complex():
    X, y = make_regression(n_samples=500, n_informative=100,
                           random_state=1)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=60, bandwidth=10,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_duplicate_data():
    X, y = make_regression(random_state=1)
    X, y = np.concatenate([X, X]), np.concatenate([y, y])
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=100, bandwidth=1,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y)/2.0, .5, decimal=2)


def test_fast_kernel_regression_conflict_data():
    X, y = make_regression(random_state=1)
    y = np.reshape(y, (-1, 1))
    X, y = X, np.hstack([y, y+2])
    # Make sure we don't throw an error when fitting or predicting
    FKR_EigenPro(kernel="linear", n_epoch=5, bandwidth=1,
                 random_state=1).fit(X, y).predict(X)


# Tests for FastKernelClassification


def test_fast_kernel_classification_gaussian():
    X, y = make_classification(n_samples=10, hypercube=False,
                               random_state=1)
    FKC_prediction = FKC_EigenPro(
        batch_size=9, kernel="gaussian", bandwidth=2.5,
        n_epoch=100, random_state=1)\
        .fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_laplace():
    X, y = make_classification(random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="laplace", n_epoch=100,
        bandwidth=13, random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_cauchy():
    X, y = make_classification(random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="cauchy", n_epoch=100, bandwidth=10,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_complex():
    X, y = make_classification(n_samples=500, n_features=500,
                               n_informative=170, scale=30, shift=6,
                               random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="gaussian", bandwidth=5,
        n_epoch=100, random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_redundant():
    X, y = make_classification(n_redundant=18, random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="laplace", bandwidth=1,
        n_epoch=100, random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_shift():
    X, y = make_classification(shift=1, hypercube=False,
                               random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="gaussian", bandwidth=5,
        n_epoch=100, random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_duplicate_data():
    X, y = make_classification(n_features=200, n_repeated=50,
                               random_state=1)
    FKC_prediction = FKC_EigenPro(
        kernel="gaussian", n_epoch=60, bandwidth=1,
        random_state=1).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_conflict_data():
    X, y = make_classification(random_state=1)
    X, y = np.concatenate([X, X]), np.concatenate([y, 1-y])
    # Make sure we don't throw an error when fitting or predicting
    FKC_EigenPro(kernel="linear", n_epoch=5, bandwidth=5,
                 random_state=1).fit(X, y).predict(X)
