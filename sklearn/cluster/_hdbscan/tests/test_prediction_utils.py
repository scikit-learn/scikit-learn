import pytest

from sklearn.cluster._hdbscan._prediction_utils import safe_always_positive_division


@pytest.mark.parametrize("denominator", [-1, 0, 1])
def test_safe_always_positive_division(denominator):
    numerator = 1
    # Given negative, zero and positive denominator and positive numerator
    value = safe_always_positive_division(numerator, denominator)
    # Make sure safe division is always positive and doesn't raise ZeroDivision error
    assert value >= 0
