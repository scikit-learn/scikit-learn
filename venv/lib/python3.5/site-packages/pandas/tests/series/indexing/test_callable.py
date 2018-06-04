import pandas as pd
import pandas.util.testing as tm


def test_getitem_callable():
    # GH 12533
    s = pd.Series(4, index=list('ABCD'))
    result = s[lambda x: 'A']
    assert result == s.loc['A']

    result = s[lambda x: ['A', 'B']]
    tm.assert_series_equal(result, s.loc[['A', 'B']])

    result = s[lambda x: [True, False, True, True]]
    tm.assert_series_equal(result, s.iloc[[0, 2, 3]])


def test_setitem_callable():
    # GH 12533
    s = pd.Series([1, 2, 3, 4], index=list('ABCD'))
    s[lambda x: 'A'] = -1
    tm.assert_series_equal(s, pd.Series([-1, 2, 3, 4], index=list('ABCD')))


def test_setitem_other_callable():
    # GH 13299
    inc = lambda x: x + 1

    s = pd.Series([1, 2, -1, 4])
    s[s < 0] = inc

    expected = pd.Series([1, 2, inc, 4])
    tm.assert_series_equal(s, expected)
