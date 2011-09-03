# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_array_almost_equal

from ..base import LinearRegression

def test_LinearRegression():
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

