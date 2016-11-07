from sklearn.utils.testing import assert_array_equal

from sklearn.utils.stats import rankdata


_cases = (
    # values, method, expected
    ([100], 'max', [1.0]),
    ([100, 100, 100], 'max', [3.0, 3.0, 3.0]),
    ([100, 300, 200], 'max', [1.0, 3.0, 2.0]),
    ([100, 200, 300, 200], 'max', [1.0, 3.0, 4.0, 3.0]),
    ([100, 200, 300, 200, 100], 'max', [2.0, 4.0, 5.0, 4.0, 2.0]),
)


def test_cases():

    def check_case(values, method, expected):
        r = rankdata(values, method=method)
        assert_array_equal(r, expected)

    for values, method, expected in _cases:
        yield check_case, values, method, expected
