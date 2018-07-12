import pytest
from sklearn.utils.testing import assert_array_equal, ignore_warnings

from sklearn.utils.stats import rankdata


_cases = (
    # values, method, expected
    ([100], 'max', [1.0]),
    ([100, 100, 100], 'max', [3.0, 3.0, 3.0]),
    ([100, 300, 200], 'max', [1.0, 3.0, 2.0]),
    ([100, 200, 300, 200], 'max', [1.0, 3.0, 4.0, 3.0]),
    ([100, 200, 300, 200, 100], 'max', [2.0, 4.0, 5.0, 4.0, 2.0]),
)


@pytest.mark.parametrize("values, method, expected", _cases)
def test_cases_rankdata(values, method, expected):

    # Test deprecated backport to be removed in 0.21
    with ignore_warnings():
        r = rankdata(values, method=method)
        assert_array_equal(r, expected)
