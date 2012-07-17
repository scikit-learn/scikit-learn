# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_array_almost_equal, assert_equal
import numpy as np
from scipy import sparse

from ..base import LinearRegression
from ...utils import check_random_state
from ...datasets.samples_generator import make_sparse_uncorrelated
from ...datasets.samples_generator import make_regression


def test_linear_regression():
    """
    Test LinearRegression on a simple dataset.
    """
    # a simple dataset
    X = [[1], [2]]
    Y = [1, 2]

    clf = LinearRegression()
    clf.fit(X, Y)

    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [1, 2])

    # test it also for degenerate input
    X = [[1]]
    Y = [0]

    clf = LinearRegression()
    clf.fit(X, Y)
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [0])


def test_fit_intercept():
    """
    Test assertions on betas shape.
    """
    X2 = np.array([[0.38349978, 0.61650022],
                   [0.58853682, 0.41146318]])
    X3 = np.array([[0.27677969, 0.70693172, 0.01628859],
                   [0.08385139, 0.20692515, 0.70922346]])
    y = np.array([1, 1])

    lr2_without_intercept = LinearRegression(fit_intercept=False).fit(X2, y)
    lr2_with_intercept = LinearRegression(fit_intercept=True).fit(X2, y)

    lr3_without_intercept = LinearRegression(fit_intercept=False).fit(X3, y)
    lr3_with_intercept = LinearRegression(fit_intercept=True).fit(X3, y)

    assert_equal(lr2_with_intercept.coef_.shape,
                 lr2_without_intercept.coef_.shape)
    assert_equal(lr3_with_intercept.coef_.shape,
                 lr3_without_intercept.coef_.shape)
    assert_equal(lr2_without_intercept.coef_.ndim,
                 lr3_without_intercept.coef_.ndim)


def test_linear_regression_sparse(random_state=0):
    "Test that linear regression also works with sparse data"
    random_state = check_random_state(random_state)
    n = 100
    X = sparse.eye(n, n)
    beta = random_state.rand(n)
    y = X * beta[:, np.newaxis]

    ols = LinearRegression()
    ols.fit(X, y.ravel())
    assert_array_almost_equal(beta, ols.coef_ + ols.intercept_)
    assert_array_almost_equal(ols.residues_, 0)


def test_linear_regression_multiple_outcome(random_state=0):
    "Test multiple-outcome linear regressions"
    X, y = make_regression(random_state=random_state)

    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    clf = LinearRegression(fit_intercept=True)
    clf.fit((X), Y)
    assert_equal(clf.coef_.shape, (2, n_features))
    Y_pred = clf.predict(X)
    clf.fit(X, y)
    y_pred = clf.predict(X)
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T, Y_pred, decimal=3)


def test_linear_regression_sparse_multiple_outcome(random_state=0):
    "Test multiple-outcome linear regressions with sparse data"
    random_state = check_random_state(random_state)
    X, y = make_sparse_uncorrelated(random_state=random_state)
    X = sparse.coo_matrix(X)
    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    ols = LinearRegression()
    ols.fit(X, Y)
    assert_equal(ols.coef_.shape, (2, n_features))
    Y_pred = ols.predict(X)
    ols.fit(X, y.ravel())
    y_pred = ols.predict(X)
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T, Y_pred, decimal=3)
