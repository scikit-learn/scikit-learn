import numpy as np
from sklearn.utils.weighted_quantile import weighted_quantile

from numpy.testing import assert_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_raises


def test_quantile_equal_weights():
    rng = np.random.RandomState(0)
    x = rng.randn(10)
    weights = 0.1 * np.ones(10)

    # since weights are equal, quantiles lie in the midpoint.
    sorted_x = np.sort(x)
    expected = 0.5 * (sorted_x[1:] + sorted_x[:-1])
    actual = np.asarray([weighted_quantile(x, q, weights) for q in np.arange(0.1, 1.0, 0.1)])

    assert_array_almost_equal(expected, actual)

    # check quantiles at (0.05, 0.95) at intervals of 0.1
    actual = np.asarray([weighted_quantile(x, q, weights) for q in np.arange(0.05, 1.05, 0.1)])
    assert_array_almost_equal(sorted_x, actual)

    # it should be the same the calculated all quantiles at the same time instead of looping over them
    assert_array_almost_equal(actual, weighted_quantile(x, weights=weights, q=np.arange(0.05, 1.05, 0.1)))


def test_quantile_toy_data():
    x = [1, 2, 3]
    weights = [1, 4, 5]

    assert_equal(weighted_quantile(x, 0.0, weights), 1)
    assert_equal(weighted_quantile(x, 1.0, weights), 3)

    assert_equal(weighted_quantile(x, 0.05, weights), 1)
    assert_almost_equal(weighted_quantile(x, 0.30, weights), 2)
    assert_equal(weighted_quantile(x, 0.75, weights), 3)
    assert_almost_equal(weighted_quantile(x, 0.50, weights), 2.44, 2)


def test_zero_weights():
    x = [1, 2, 3, 4, 5]
    w = [0, 0, 0, 0.1, 0.1]

    for q in np.arange(0.0, 1.10, 0.1):
        assert_equal(
            weighted_quantile(x, q, w),
            weighted_quantile([4, 5], q, [0.1, 0.1])
        )


def test_xd_shapes():
    rng = np.random.RandomState(0)
    x = rng.randn(100, 10, 20)
    weights = 0.01 * np.ones_like(x)

    # shape should be the same as the output of np.quantile
    assert weighted_quantile(x, 0.5, weights, axis=0).shape == np.quantile(x, 0.5, axis=0).shape
    assert weighted_quantile(x, 0.5, weights, axis=1).shape == np.quantile(x, 0.5, axis=1).shape
    assert weighted_quantile(x, 0.5, weights, axis=2).shape == np.quantile(x, 0.5, axis=2).shape
    assert isinstance(weighted_quantile(x, 0.5, weights, axis=None), float)
    assert weighted_quantile(x, (0.5, 0.8), weights, axis=0).shape == np.quantile(x, (0.5, 0.8), axis=0).shape

    # axis should be integer
    assert_raises(NotImplementedError, weighted_quantile, x, 0.5, weights, axis=(1, 2))

    # weighted_quantile should yield very similar results to np.quantile
    assert np.allclose(weighted_quantile(x, 0.5, weights, axis=2), np.quantile(x, q=0.5, axis=2))


if __name__ == "sklearn.utils.tests.test_utils":
    print("Test utils")
    test_quantile_equal_weights()
    test_quantile_toy_data()
    test_zero_weights()
    test_xd_shapes()
