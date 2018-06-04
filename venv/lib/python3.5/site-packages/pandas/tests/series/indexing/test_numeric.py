# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame)

from pandas.compat import lrange, range
from pandas.util.testing import (assert_series_equal)

import pandas.util.testing as tm


def test_get():
    # GH 6383
    s = Series(np.array([43, 48, 60, 48, 50, 51, 50, 45, 57, 48, 56, 45,
                         51, 39, 55, 43, 54, 52, 51, 54]))

    result = s.get(25, 0)
    expected = 0
    assert result == expected

    s = Series(np.array([43, 48, 60, 48, 50, 51, 50, 45, 57, 48, 56,
                         45, 51, 39, 55, 43, 54, 52, 51, 54]),
               index=pd.Float64Index(
                   [25.0, 36.0, 49.0, 64.0, 81.0, 100.0,
                    121.0, 144.0, 169.0, 196.0, 1225.0,
                    1296.0, 1369.0, 1444.0, 1521.0, 1600.0,
                    1681.0, 1764.0, 1849.0, 1936.0],
                   dtype='object'))

    result = s.get(25, 0)
    expected = 43
    assert result == expected

    # GH 7407
    # with a boolean accessor
    df = pd.DataFrame({'i': [0] * 3, 'b': [False] * 3})
    vc = df.i.value_counts()
    result = vc.get(99, default='Missing')
    assert result == 'Missing'

    vc = df.b.value_counts()
    result = vc.get(False, default='Missing')
    assert result == 3

    result = vc.get(True, default='Missing')
    assert result == 'Missing'


def test_get_nan():
    # GH 8569
    s = pd.Float64Index(range(10)).to_series()
    assert s.get(np.nan) is None
    assert s.get(np.nan, default='Missing') == 'Missing'


def test_get_nan_multiple():
    # GH 8569
    # ensure that fixing "test_get_nan" above hasn't broken get
    # with multiple elements
    s = pd.Float64Index(range(10)).to_series()

    idx = [2, 30]
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        assert_series_equal(s.get(idx),
                            Series([2, np.nan], index=idx))

    idx = [2, np.nan]
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        assert_series_equal(s.get(idx),
                            Series([2, np.nan], index=idx))

    # GH 17295 - all missing keys
    idx = [20, 30]
    assert(s.get(idx) is None)

    idx = [np.nan, np.nan]
    assert(s.get(idx) is None)


def test_delitem():
    # GH 5542
    # should delete the item inplace
    s = Series(lrange(5))
    del s[0]

    expected = Series(lrange(1, 5), index=lrange(1, 5))
    assert_series_equal(s, expected)

    del s[1]
    expected = Series(lrange(2, 5), index=lrange(2, 5))
    assert_series_equal(s, expected)

    # empty
    s = Series()

    def f():
        del s[0]

    pytest.raises(KeyError, f)

    # only 1 left, del, add, del
    s = Series(1)
    del s[0]
    assert_series_equal(s, Series(dtype='int64', index=Index(
        [], dtype='int64')))
    s[0] = 1
    assert_series_equal(s, Series(1))
    del s[0]
    assert_series_equal(s, Series(dtype='int64', index=Index(
        [], dtype='int64')))

    # Index(dtype=object)
    s = Series(1, index=['a'])
    del s['a']
    assert_series_equal(s, Series(dtype='int64', index=Index(
        [], dtype='object')))
    s['a'] = 1
    assert_series_equal(s, Series(1, index=['a']))
    del s['a']
    assert_series_equal(s, Series(dtype='int64', index=Index(
        [], dtype='object')))


def test_slice_float64():
    values = np.arange(10., 50., 2)
    index = Index(values)

    start, end = values[[5, 15]]

    s = Series(np.random.randn(20), index=index)

    result = s[start:end]
    expected = s.iloc[5:16]
    assert_series_equal(result, expected)

    result = s.loc[start:end]
    assert_series_equal(result, expected)

    df = DataFrame(np.random.randn(20, 3), index=index)

    result = df[start:end]
    expected = df.iloc[5:16]
    tm.assert_frame_equal(result, expected)

    result = df.loc[start:end]
    tm.assert_frame_equal(result, expected)


def test_getitem_negative_out_of_bounds():
    s = Series(tm.rands_array(5, 10), index=tm.rands_array(10, 10))

    pytest.raises(IndexError, s.__getitem__, -11)
    pytest.raises(IndexError, s.__setitem__, -11, 'foo')


def test_getitem_regression():
    s = Series(lrange(5), index=lrange(5))
    result = s[lrange(5)]
    assert_series_equal(result, s)


def test_getitem_setitem_slice_bug():
    s = Series(lrange(10), lrange(10))
    result = s[-12:]
    assert_series_equal(result, s)

    result = s[-7:]
    assert_series_equal(result, s[3:])

    result = s[:-12]
    assert_series_equal(result, s[:0])

    s = Series(lrange(10), lrange(10))
    s[-12:] = 0
    assert (s == 0).all()

    s[:-12] = 5
    assert (s == 0).all()


def test_getitem_setitem_slice_integers():
    s = Series(np.random.randn(8), index=[2, 4, 6, 8, 10, 12, 14, 16])

    result = s[:4]
    expected = s.reindex([2, 4, 6, 8])
    assert_series_equal(result, expected)

    s[:4] = 0
    assert (s[:4] == 0).all()
    assert not (s[4:] == 0).any()


def test_setitem_float_labels():
    # note labels are floats
    s = Series(['a', 'b', 'c'], index=[0, 0.5, 1])
    tmp = s.copy()

    s.loc[1] = 'zoo'
    tmp.iloc[2] = 'zoo'

    assert_series_equal(s, tmp)


def test_slice_float_get_set(test_data):
    pytest.raises(TypeError, lambda: test_data.ts[4.0:10.0])

    def f():
        test_data.ts[4.0:10.0] = 0

    pytest.raises(TypeError, f)

    pytest.raises(TypeError, test_data.ts.__getitem__, slice(4.5, 10.0))
    pytest.raises(TypeError, test_data.ts.__setitem__, slice(4.5, 10.0), 0)


def test_slice_floats2():
    s = Series(np.random.rand(10), index=np.arange(10, 20, dtype=float))

    assert len(s.loc[12.0:]) == 8
    assert len(s.loc[12.5:]) == 7

    i = np.arange(10, 20, dtype=float)
    i[2] = 12.2
    s.index = i
    assert len(s.loc[12.0:]) == 8
    assert len(s.loc[12.5:]) == 7


def test_int_indexing():
    s = Series(np.random.randn(6), index=[0, 0, 1, 1, 2, 2])

    pytest.raises(KeyError, s.__getitem__, 5)

    pytest.raises(KeyError, s.__getitem__, 'c')

    # not monotonic
    s = Series(np.random.randn(6), index=[2, 2, 0, 0, 1, 1])

    pytest.raises(KeyError, s.__getitem__, 5)

    pytest.raises(KeyError, s.__getitem__, 'c')


def test_getitem_int64(test_data):
    idx = np.int64(5)
    assert test_data.ts[idx] == test_data.ts[5]
