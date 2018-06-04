# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

import numpy as np
import pandas as pd

from pandas import (Series, Timestamp)

from pandas.compat import lrange
from pandas.util.testing import (assert_series_equal)


def test_loc_getitem(test_data):
    inds = test_data.series.index[[3, 4, 7]]
    assert_series_equal(
        test_data.series.loc[inds],
        test_data.series.reindex(inds))
    assert_series_equal(test_data.series.iloc[5::2], test_data.series[5::2])

    # slice with indices
    d1, d2 = test_data.ts.index[[5, 15]]
    result = test_data.ts.loc[d1:d2]
    expected = test_data.ts.truncate(d1, d2)
    assert_series_equal(result, expected)

    # boolean
    mask = test_data.series > test_data.series.median()
    assert_series_equal(test_data.series.loc[mask], test_data.series[mask])

    # ask for index value
    assert test_data.ts.loc[d1] == test_data.ts[d1]
    assert test_data.ts.loc[d2] == test_data.ts[d2]


def test_loc_getitem_not_monotonic(test_data):
    d1, d2 = test_data.ts.index[[5, 15]]

    ts2 = test_data.ts[::2][[1, 2, 0]]

    pytest.raises(KeyError, ts2.loc.__getitem__, slice(d1, d2))
    pytest.raises(KeyError, ts2.loc.__setitem__, slice(d1, d2), 0)


def test_loc_getitem_setitem_integer_slice_keyerrors():
    s = Series(np.random.randn(10), index=lrange(0, 20, 2))

    # this is OK
    cp = s.copy()
    cp.iloc[4:10] = 0
    assert (cp.iloc[4:10] == 0).all()

    # so is this
    cp = s.copy()
    cp.iloc[3:11] = 0
    assert (cp.iloc[3:11] == 0).values.all()

    result = s.iloc[2:6]
    result2 = s.loc[3:11]
    expected = s.reindex([4, 6, 8, 10])

    assert_series_equal(result, expected)
    assert_series_equal(result2, expected)

    # non-monotonic, raise KeyError
    s2 = s.iloc[lrange(5) + lrange(5, 10)[::-1]]
    pytest.raises(KeyError, s2.loc.__getitem__, slice(3, 11))
    pytest.raises(KeyError, s2.loc.__setitem__, slice(3, 11), 0)


def test_loc_getitem_iterator(test_data):
    idx = iter(test_data.series.index[:10])
    result = test_data.series.loc[idx]
    assert_series_equal(result, test_data.series[:10])


def test_loc_setitem_boolean(test_data):
    mask = test_data.series > test_data.series.median()

    result = test_data.series.copy()
    result.loc[mask] = 0
    expected = test_data.series
    expected[mask] = 0
    assert_series_equal(result, expected)


def test_loc_setitem_corner(test_data):
    inds = list(test_data.series.index[[5, 8, 12]])
    test_data.series.loc[inds] = 5
    pytest.raises(Exception, test_data.series.loc.__setitem__,
                  inds + ['foo'], 5)


def test_basic_setitem_with_labels(test_data):
    indices = test_data.ts.index[[5, 10, 15]]

    cp = test_data.ts.copy()
    exp = test_data.ts.copy()
    cp[indices] = 0
    exp.loc[indices] = 0
    assert_series_equal(cp, exp)

    cp = test_data.ts.copy()
    exp = test_data.ts.copy()
    cp[indices[0]:indices[2]] = 0
    exp.loc[indices[0]:indices[2]] = 0
    assert_series_equal(cp, exp)

    # integer indexes, be careful
    s = Series(np.random.randn(10), index=lrange(0, 20, 2))
    inds = [0, 4, 6]
    arr_inds = np.array([0, 4, 6])

    cp = s.copy()
    exp = s.copy()
    s[inds] = 0
    s.loc[inds] = 0
    assert_series_equal(cp, exp)

    cp = s.copy()
    exp = s.copy()
    s[arr_inds] = 0
    s.loc[arr_inds] = 0
    assert_series_equal(cp, exp)

    inds_notfound = [0, 4, 5, 6]
    arr_inds_notfound = np.array([0, 4, 5, 6])
    pytest.raises(Exception, s.__setitem__, inds_notfound, 0)
    pytest.raises(Exception, s.__setitem__, arr_inds_notfound, 0)

    # GH12089
    # with tz for values
    s = Series(pd.date_range("2011-01-01", periods=3, tz="US/Eastern"),
               index=['a', 'b', 'c'])
    s2 = s.copy()
    expected = Timestamp('2011-01-03', tz='US/Eastern')
    s2.loc['a'] = expected
    result = s2.loc['a']
    assert result == expected

    s2 = s.copy()
    s2.iloc[0] = expected
    result = s2.iloc[0]
    assert result == expected

    s2 = s.copy()
    s2['a'] = expected
    result = s2['a']
    assert result == expected
