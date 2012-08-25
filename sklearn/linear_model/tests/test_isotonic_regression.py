import numpy as np
from numpy.testing import assert_array_equal

from sklearn.linear_model.isotonic_regression_ import isotonic_regression
from sklearn.linear_model import IsotonicRegression

from nose.tools import assert_raises


def test_isotonic_regression():
    X = np.array([3, 7, 5, 9, 8, 7, 10])
    Y = np.array([3, 6, 6, 8, 8, 8, 10])
    assert_array_equal(Y, isotonic_regression(X))

def assert_raises_exceptions():
    assert_raises(ValueError, IsotonicRegression.fit([5, 7, 3],
                                                     w=[0.1, 0.6]))
