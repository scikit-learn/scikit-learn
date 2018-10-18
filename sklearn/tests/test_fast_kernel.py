import numpy as np

from sklearn.datasets import make_regression, make_classification
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.fast_kernel import FKR_EigenPro, FKC_EigenPro


from sklearn.utils.estimator_checks import check_estimator
# check_estimator(FKR_EigenPro)
# check_estimator(FKC_EigenPro)

np.random.seed(1)
#           Tests for Fast Kernel Classification


def test_fast_kernel_regression_gaussian():
    X, y = make_regression()
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_laplace():
    X, y = make_regression()
    FKR_prediction = FKR_EigenPro(
        kernel="laplace", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_cauchy():
    X, y = make_regression()
    FKR_prediction = FKR_EigenPro(
        kernel="cauchy", n_epoch=50, bandwidth=10, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_2d():
    X, Y = make_regression(n_targets=30)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X,Y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / Y), 1, decimal=2)


def test_fast_kernel_regression_many_features():
    X, y = make_regression(n_features=10000)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X,y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_simple():
    X,y = make_regression(n_informative=1)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=200, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_complex():
    X, y = make_regression(n_samples=500, n_informative=500)
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(y)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_duplicate_data():
    X, y = make_regression()
    X, y = np.concatenate([X,X]), np.concatenate([y,y])
    FKR_prediction = FKR_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_conflict_data():
    X, y = make_regression()
    X, y = X, np.concatenate([y, y+2])
    try:
        FKR_prediction = FKR_EigenPro(
            kernel="linear", n_epoch=5, bandwidth=1, random_state=0).fit(X, y).predict(X)
        assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=1)
        raise AssertionError("Predicted impossibly well")
    except AssertionError:
        return


#       Tests for FastKernelClassification


def test_fast_kernel_classification_gaussian():
    X, y = make_classification()
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_laplace():
    X, y = make_classification()
    FKC_prediction = FKC_EigenPro(
        kernel="laplace", n_epoch=50, bandwidth=13).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_cauchy():
    X, y = make_classification()
    FKC_prediction = FKC_EigenPro(
        kernel="cauchy", n_epoch=50, bandwidth=10, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_complex():
    X, y = make_classification(n_samples=300, n_features=500, n_informative=170)
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_redundant():
    X, y = make_classification(n_redundant=18)
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=.01, n_epoch=50).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_shift():
    X, y = make_classification(shift=1, hypercube=False)
    FKC_prediction = FKC_EigenPro(kernel="gaussian", bandwidth=5, n_epoch=50).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_duplicate_data():
    X, y = make_classification(n_features=200, n_repeated=50)
    FKC_prediction = FKC_EigenPro(
        kernel="gaussian", n_epoch=50, bandwidth=.1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(FKC_prediction, y)


def test_fast_kernel_classification_conflict_data():
    X, y = make_classification()
    X, y = np.concatenate([X, X]), np.concatenate([y, 1-y])
    try:
        FKR_prediction = FKC_EigenPro(
            kernel="linear", n_epoch=5, bandwidth=5, random_state=0).fit(X, y).predict(X)
        assert_array_almost_equal(FKR_prediction, y)
        raise AssertionError("Predicted impossibly well")
    except AssertionError:
        return
