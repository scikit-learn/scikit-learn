import numpy as np

from sklearn.datasets import make_regression
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.fast_kernel_regression import FastKernelRegression

np.random.seed(1)
X, y = make_regression()
X2, Y2 = make_regression(n_targets=30)
X3, y3 = make_regression(n_features=10000)
X4, y4 = make_regression(n_informative=1)
X5, y5 = make_regression(n_samples=500, n_informative=500)
X6, y6 = np.concatenate([X,X]), np.concatenate([y,y])
X7, y7 = X6, np.concatenate([y, y+2])

def test_fast_kernel_regression_gaussian():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_laplace():
    FKR_prediction = FastKernelRegression(
        kernel="laplace", n_epoch=50, bandwidth=1, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_cauchy():
    FKR_prediction = FastKernelRegression(
        kernel="cauchy", n_epoch=50, bandwidth=10, random_state=0).fit(X, y).predict(X)
    assert_array_almost_equal(abs(FKR_prediction / y), 1, decimal=2)


def test_fast_kernel_regression_2d():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X2,Y2).predict(X2)
    assert_array_almost_equal(abs(FKR_prediction / Y2), 1, decimal=2)

def test_fast_kernel_regression_many_features():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X3,y3).predict(X3)
    assert_array_almost_equal(abs(FKR_prediction / y3), 1, decimal=2)

def test_fast_kernel_regression_simple():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=200, bandwidth=1, random_state=0).fit(X4, y4).predict(X4)
    assert_array_almost_equal(abs(FKR_prediction / y4), 1, decimal=2)

def test_fast_kernel_regression_complex():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X5, y5).predict(X5)
    assert_array_almost_equal(abs(FKR_prediction / y5), 1, decimal=2)

def test_fast_kernel_regression_duplicate_data():
    FKR_prediction = FastKernelRegression(
        kernel="gaussian", n_epoch=50, bandwidth=1, random_state=0).fit(X6, y6).predict(X6)
    assert_array_almost_equal(abs(FKR_prediction / y6), 1, decimal=2)

def test_fast_kernel_regression_conflict_data():
    try:
        FKR_prediction = FastKernelRegression(
            kernel="gaussian", n_epoch=5, bandwidth=1, random_state=0).fit(X7, y7).predict(X7)
        assert_array_almost_equal(abs(FKR_prediction / y7), 1, decimal=1)
        raise AssertionError("Predicted impossibly well")
    except AssertionError:
        return
