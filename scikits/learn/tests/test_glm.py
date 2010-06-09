from numpy.testing import *

from scikits.learn import glm

X = [[1, 0, -.1],
     [0, 0, 0],
     [0, 1, .1]]
Y = [1, 0, -1]


def test_1():
    """
    Very simple test
    """
    clf = glm.LassoLARS().fit(X, Y)
    assert_array_almost_equal(clf.coef_, [-1.4142, -1.4142, 0], decimal=4)



