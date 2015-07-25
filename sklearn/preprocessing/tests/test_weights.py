from sklearn.preprocessing._weights import _balance_weights
from sklearn.utils.testing import assert_array_equal


def test_balance_weights():
    weights = _balance_weights([0, 0, 1, 1])
    assert_array_equal(weights, [1., 1., 1., 1.])

    weights = _balance_weights([0, 1, 1, 1, 1])
    assert_array_equal(weights, [1., 0.25, 0.25, 0.25, 0.25])

    weights = _balance_weights([0, 0])
    assert_array_equal(weights, [1., 1.])
