# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_array_almost_equal
import numpy as np
from scipy import sparse

from ..base import LinearRegression
from ...utils import check_random_state


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
