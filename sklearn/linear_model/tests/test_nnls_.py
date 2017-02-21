from sklearn.linear_model.nnls_ import NNLS
from sklearn.utils.testing import assert_array_almost_equal


def test_nnls():
    # Test a Non-negative least squares on a simple dataset.

    X = [[1], [2]]
    y = [1, 2]

    reg = NNLS()
    reg.fit(X, y)
    assert_array_almost_equal(reg.coef_, [1])
    assert_array_almost_equal(reg.intercept_, [0])
    assert_array_almost_equal(reg.predict(X), [1, 2])
