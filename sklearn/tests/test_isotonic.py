import numpy as np

from sklearn.isotonic import isotonic_regression, IsotonicRegression

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_array_equal


def test_isotonic_regression():
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    y_ = np.array([3, 6, 6, 8, 8, 8, 10])
    assert_array_equal(y_, isotonic_regression(y))

    x = np.arange(len(y))
    ir = IsotonicRegression(y_min=0., y_max=1.)
    ir.fit(x, y)
    assert_array_equal(ir.fit(x, y).transform(x), ir.fit_transform(x, y))
    assert_array_equal(ir.transform(x), ir.predict(x))

    # check that it is immune to permutation
    perm = np.random.permutation(len(y))
    ir = IsotonicRegression(y_min=0., y_max=1.)
    assert_array_equal(ir.fit_transform(x[perm], y[perm]),
                       ir.fit_transform(x, y)[perm])
    assert_array_equal(ir.transform(x[perm]), ir.transform(x)[perm])

    # check it doesn't change y when all x are equal:
    ir = IsotonicRegression()
    assert_array_equal(ir.fit_transform(np.ones(len(x)), y), y)


def test_assert_raises_exceptions():
    ir = IsotonicRegression()
    rng = np.random.RandomState(42)
    assert_raises(ValueError, ir.fit, [0, 1, 2], [5, 7, 3], [0.1, 0.6])
    assert_raises(ValueError, ir.fit, [0, 1, 2], [5, 7])
    assert_raises(ValueError, ir.fit, rng.randn(3, 10), [0, 1, 2])
    assert_raises(ValueError, ir.transform, rng.randn(3, 10))
