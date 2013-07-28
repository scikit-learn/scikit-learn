from sklearn.preprocessing import balance_weights
from sklearn.utils.testing import assert_array_equal


def test_balance_weights():
    weights = balance_weights([0, 0, 1, 1])
    assert_array_equal(weights, [1., 1., 1., 1.])

    weights = balance_weights([0, 1, 1, 1, 1])
    assert_array_equal(weights, [1., 0.25, 0.25, 0.25, 0.25])

    weights = balance_weights([0, 0])
    assert_array_equal(weights, [1., 1.])
