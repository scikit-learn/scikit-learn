# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

from datetime import datetime

import pandas as pd
import numpy as np

from numpy import nan

from pandas import compat

from pandas import (Series, date_range, isna, Categorical)
from pandas.compat import lrange, range

from pandas.util.testing import (assert_series_equal)
import pandas.util.testing as tm


@pytest.mark.parametrize(
    'first_slice,second_slice', [
        [[2, None], [None, -5]],
        [[None, 0], [None, -5]],
        [[None, -5], [None, 0]],
        [[None, 0], [None, 0]]
    ])
@pytest.mark.parametrize('fill', [None, -1])
def test_align(test_data, first_slice, second_slice, join_type, fill):
    a = test_data.ts[slice(*first_slice)]
    b = test_data.ts[slice(*second_slice)]

    aa, ab = a.align(b, join=join_type, fill_value=fill)

    join_index = a.index.join(b.index, how=join_type)
    if fill is not None:
        diff_a = aa.index.difference(join_index)
        diff_b = ab.index.difference(join_index)
        if len(diff_a) > 0:
            assert (aa.reindex(diff_a) == fill).all()
        if len(diff_b) > 0:
            assert (ab.reindex(diff_b) == fill).all()

    ea = a.reindex(join_index)
    eb = b.reindex(join_index)

    if fill is not None:
        ea = ea.fillna(fill)
        eb = eb.fillna(fill)

    assert_series_equal(aa, ea)
    assert_series_equal(ab, eb)
    assert aa.name == 'ts'
    assert ea.name == 'ts'
    assert ab.name == 'ts'
    assert eb.name == 'ts'


@pytest.mark.parametrize(
    'first_slice,second_slice', [
        [[2, None], [None, -5]],
        [[None, 0], [None, -5]],
        [[None, -5], [None, 0]],
        [[None, 0], [None, 0]]
    ])
@pytest.mark.parametrize('method', ['pad', 'bfill'])
@pytest.mark.parametrize('limit', [None, 1])
def test_align_fill_method(test_data,
                           first_slice, second_slice,
                           join_type, method, limit):
    a = test_data.ts[slice(*first_slice)]
    b = test_data.ts[slice(*second_slice)]

    aa, ab = a.align(b, join=join_type, method=method, limit=limit)

    join_index = a.index.join(b.index, how=join_type)
    ea = a.reindex(join_index)
    eb = b.reindex(join_index)

    ea = ea.fillna(method=method, limit=limit)
    eb = eb.fillna(method=method, limit=limit)

    assert_series_equal(aa, ea)
    assert_series_equal(ab, eb)


def test_align_nocopy(test_data):
    b = test_data.ts[:5].copy()

    # do copy
    a = test_data.ts.copy()
    ra, _ = a.align(b, join='left')
    ra[:5] = 5
    assert not (a[:5] == 5).any()

    # do not copy
    a = test_data.ts.copy()
    ra, _ = a.align(b, join='left', copy=False)
    ra[:5] = 5
    assert (a[:5] == 5).all()

    # do copy
    a = test_data.ts.copy()
    b = test_data.ts[:5].copy()
    _, rb = a.align(b, join='right')
    rb[:3] = 5
    assert not (b[:3] == 5).any()

    # do not copy
    a = test_data.ts.copy()
    b = test_data.ts[:5].copy()
    _, rb = a.align(b, join='right', copy=False)
    rb[:2] = 5
    assert (b[:2] == 5).all()


def test_align_same_index(test_data):
    a, b = test_data.ts.align(test_data.ts, copy=False)
    assert a.index is test_data.ts.index
    assert b.index is test_data.ts.index

    a, b = test_data.ts.align(test_data.ts, copy=True)
    assert a.index is not test_data.ts.index
    assert b.index is not test_data.ts.index


def test_align_multiindex():
    # GH 10665

    midx = pd.MultiIndex.from_product([range(2), range(3), range(2)],
                                      names=('a', 'b', 'c'))
    idx = pd.Index(range(2), name='b')
    s1 = pd.Series(np.arange(12, dtype='int64'), index=midx)
    s2 = pd.Series(np.arange(2, dtype='int64'), index=idx)

    # these must be the same results (but flipped)
    res1l, res1r = s1.align(s2, join='left')
    res2l, res2r = s2.align(s1, join='right')

    expl = s1
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = pd.Series([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)

    res1l, res1r = s1.align(s2, join='right')
    res2l, res2r = s2.align(s1, join='left')

    exp_idx = pd.MultiIndex.from_product([range(2), range(2), range(2)],
                                         names=('a', 'b', 'c'))
    expl = pd.Series([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = pd.Series([0, 0, 1, 1] * 2, index=exp_idx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)


def test_reindex(test_data):
    identity = test_data.series.reindex(test_data.series.index)

    # __array_interface__ is not defined for older numpies
    # and on some pythons
    try:
        assert np.may_share_memory(test_data.series.index, identity.index)
    except AttributeError:
        pass

    assert identity.index.is_(test_data.series.index)
    assert identity.index.identical(test_data.series.index)

    subIndex = test_data.series.index[10:20]
    subSeries = test_data.series.reindex(subIndex)

    for idx, val in compat.iteritems(subSeries):
        assert val == test_data.series[idx]

    subIndex2 = test_data.ts.index[10:20]
    subTS = test_data.ts.reindex(subIndex2)

    for idx, val in compat.iteritems(subTS):
        assert val == test_data.ts[idx]
    stuffSeries = test_data.ts.reindex(subIndex)

    assert np.isnan(stuffSeries).all()

    # This is extremely important for the Cython code to not screw up
    nonContigIndex = test_data.ts.index[::2]
    subNonContig = test_data.ts.reindex(nonContigIndex)
    for idx, val in compat.iteritems(subNonContig):
        assert val == test_data.ts[idx]

    # return a copy the same index here
    result = test_data.ts.reindex()
    assert not (result is test_data.ts)


def test_reindex_nan():
    ts = Series([2, 3, 5, 7], index=[1, 4, nan, 8])

    i, j = [nan, 1, nan, 8, 4, nan], [2, 0, 2, 3, 1, 2]
    assert_series_equal(ts.reindex(i), ts.iloc[j])

    ts.index = ts.index.astype('object')

    # reindex coerces index.dtype to float, loc/iloc doesn't
    assert_series_equal(ts.reindex(i), ts.iloc[j], check_index_type=False)


def test_reindex_series_add_nat():
    rng = date_range('1/1/2000 00:00:00', periods=10, freq='10s')
    series = Series(rng)

    result = series.reindex(lrange(15))
    assert np.issubdtype(result.dtype, np.dtype('M8[ns]'))

    mask = result.isna()
    assert mask[-5:].all()
    assert not mask[:-5].any()


def test_reindex_with_datetimes():
    rng = date_range('1/1/2000', periods=20)
    ts = Series(np.random.randn(20), index=rng)

    result = ts.reindex(list(ts.index[5:10]))
    expected = ts[5:10]
    tm.assert_series_equal(result, expected)

    result = ts[list(ts.index[5:10])]
    tm.assert_series_equal(result, expected)


def test_reindex_corner(test_data):
    # (don't forget to fix this) I think it's fixed
    test_data.empty.reindex(test_data.ts.index, method='pad')  # it works

    # corner case: pad empty series
    reindexed = test_data.empty.reindex(test_data.ts.index, method='pad')

    # pass non-Index
    reindexed = test_data.ts.reindex(list(test_data.ts.index))
    assert_series_equal(test_data.ts, reindexed)

    # bad fill method
    ts = test_data.ts[::2]
    pytest.raises(Exception, ts.reindex, test_data.ts.index, method='foo')


def test_reindex_pad():
    s = Series(np.arange(10), dtype='int64')
    s2 = s[::2]

    reindexed = s2.reindex(s.index, method='pad')
    reindexed2 = s2.reindex(s.index, method='ffill')
    assert_series_equal(reindexed, reindexed2)

    expected = Series([0, 0, 2, 2, 4, 4, 6, 6, 8, 8], index=np.arange(10))
    assert_series_equal(reindexed, expected)

    # GH4604
    s = Series([1, 2, 3, 4, 5], index=['a', 'b', 'c', 'd', 'e'])
    new_index = ['a', 'g', 'c', 'f']
    expected = Series([1, 1, 3, 3], index=new_index)

    # this changes dtype because the ffill happens after
    result = s.reindex(new_index).ffill()
    assert_series_equal(result, expected.astype('float64'))

    result = s.reindex(new_index).ffill(downcast='infer')
    assert_series_equal(result, expected)

    expected = Series([1, 5, 3, 5], index=new_index)
    result = s.reindex(new_index, method='ffill')
    assert_series_equal(result, expected)

    # inference of new dtype
    s = Series([True, False, False, True], index=list('abcd'))
    new_index = 'agc'
    result = s.reindex(list(new_index)).ffill()
    expected = Series([True, True, False], index=list(new_index))
    assert_series_equal(result, expected)

    # GH4618 shifted series downcasting
    s = Series(False, index=lrange(0, 5))
    result = s.shift(1).fillna(method='bfill')
    expected = Series(False, index=lrange(0, 5))
    assert_series_equal(result, expected)


def test_reindex_nearest():
    s = Series(np.arange(10, dtype='int64'))
    target = [0.1, 0.9, 1.5, 2.0]
    actual = s.reindex(target, method='nearest')
    expected = Series(np.around(target).astype('int64'), target)
    assert_series_equal(expected, actual)

    actual = s.reindex_like(actual, method='nearest')
    assert_series_equal(expected, actual)

    actual = s.reindex_like(actual, method='nearest', tolerance=1)
    assert_series_equal(expected, actual)
    actual = s.reindex_like(actual, method='nearest',
                            tolerance=[1, 2, 3, 4])
    assert_series_equal(expected, actual)

    actual = s.reindex(target, method='nearest', tolerance=0.2)
    expected = Series([0, 1, np.nan, 2], target)
    assert_series_equal(expected, actual)

    actual = s.reindex(target, method='nearest',
                       tolerance=[0.3, 0.01, 0.4, 3])
    expected = Series([0, np.nan, np.nan, 2], target)
    assert_series_equal(expected, actual)


def test_reindex_backfill():
    pass


def test_reindex_int(test_data):
    ts = test_data.ts[::2]
    int_ts = Series(np.zeros(len(ts), dtype=int), index=ts.index)

    # this should work fine
    reindexed_int = int_ts.reindex(test_data.ts.index)

    # if NaNs introduced
    assert reindexed_int.dtype == np.float_

    # NO NaNs introduced
    reindexed_int = int_ts.reindex(int_ts.index[::2])
    assert reindexed_int.dtype == np.int_


def test_reindex_bool(test_data):
    # A series other than float, int, string, or object
    ts = test_data.ts[::2]
    bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)

    # this should work fine
    reindexed_bool = bool_ts.reindex(test_data.ts.index)

    # if NaNs introduced
    assert reindexed_bool.dtype == np.object_

    # NO NaNs introduced
    reindexed_bool = bool_ts.reindex(bool_ts.index[::2])
    assert reindexed_bool.dtype == np.bool_


def test_reindex_bool_pad(test_data):
    # fail
    ts = test_data.ts[5:]
    bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)
    filled_bool = bool_ts.reindex(test_data.ts.index, method='pad')
    assert isna(filled_bool[:5]).all()


def test_reindex_categorical():
    index = date_range('20000101', periods=3)

    # reindexing to an invalid Categorical
    s = Series(['a', 'b', 'c'], dtype='category')
    result = s.reindex(index)
    expected = Series(Categorical(values=[np.nan, np.nan, np.nan],
                                  categories=['a', 'b', 'c']))
    expected.index = index
    tm.assert_series_equal(result, expected)

    # partial reindexing
    expected = Series(Categorical(values=['b', 'c'], categories=['a', 'b',
                                                                 'c']))
    expected.index = [1, 2]
    result = s.reindex([1, 2])
    tm.assert_series_equal(result, expected)

    expected = Series(Categorical(
        values=['c', np.nan], categories=['a', 'b', 'c']))
    expected.index = [2, 3]
    result = s.reindex([2, 3])
    tm.assert_series_equal(result, expected)


def test_reindex_like(test_data):
    other = test_data.ts[::2]
    assert_series_equal(test_data.ts.reindex(other.index),
                        test_data.ts.reindex_like(other))

    # GH 7179
    day1 = datetime(2013, 3, 5)
    day2 = datetime(2013, 5, 5)
    day3 = datetime(2014, 3, 5)

    series1 = Series([5, None, None], [day1, day2, day3])
    series2 = Series([None, None], [day1, day3])

    result = series1.reindex_like(series2, method='pad')
    expected = Series([5, np.nan], index=[day1, day3])
    assert_series_equal(result, expected)


def test_reindex_fill_value():
    # -----------------------------------------------------------
    # floats
    floats = Series([1., 2., 3.])
    result = floats.reindex([1, 2, 3])
    expected = Series([2., 3., np.nan], index=[1, 2, 3])
    assert_series_equal(result, expected)

    result = floats.reindex([1, 2, 3], fill_value=0)
    expected = Series([2., 3., 0], index=[1, 2, 3])
    assert_series_equal(result, expected)

    # -----------------------------------------------------------
    # ints
    ints = Series([1, 2, 3])

    result = ints.reindex([1, 2, 3])
    expected = Series([2., 3., np.nan], index=[1, 2, 3])
    assert_series_equal(result, expected)

    # don't upcast
    result = ints.reindex([1, 2, 3], fill_value=0)
    expected = Series([2, 3, 0], index=[1, 2, 3])
    assert issubclass(result.dtype.type, np.integer)
    assert_series_equal(result, expected)

    # -----------------------------------------------------------
    # objects
    objects = Series([1, 2, 3], dtype=object)

    result = objects.reindex([1, 2, 3])
    expected = Series([2, 3, np.nan], index=[1, 2, 3], dtype=object)
    assert_series_equal(result, expected)

    result = objects.reindex([1, 2, 3], fill_value='foo')
    expected = Series([2, 3, 'foo'], index=[1, 2, 3], dtype=object)
    assert_series_equal(result, expected)

    # ------------------------------------------------------------
    # bools
    bools = Series([True, False, True])

    result = bools.reindex([1, 2, 3])
    expected = Series([False, True, np.nan], index=[1, 2, 3], dtype=object)
    assert_series_equal(result, expected)

    result = bools.reindex([1, 2, 3], fill_value=False)
    expected = Series([False, True, False], index=[1, 2, 3])
    assert_series_equal(result, expected)


def test_rename():
    # GH 17407
    s = Series(range(1, 6), index=pd.Index(range(2, 7), name='IntIndex'))
    result = s.rename(str)
    expected = s.rename(lambda i: str(i))
    assert_series_equal(result, expected)

    assert result.name == expected.name


def test_drop():
    # unique
    s = Series([1, 2], index=['one', 'two'])
    expected = Series([1], index=['one'])
    result = s.drop(['two'])
    assert_series_equal(result, expected)
    result = s.drop('two', axis='rows')
    assert_series_equal(result, expected)

    # non-unique
    # GH 5248
    s = Series([1, 1, 2], index=['one', 'two', 'one'])
    expected = Series([1, 2], index=['one', 'one'])
    result = s.drop(['two'], axis=0)
    assert_series_equal(result, expected)
    result = s.drop('two')
    assert_series_equal(result, expected)

    expected = Series([1], index=['two'])
    result = s.drop(['one'])
    assert_series_equal(result, expected)
    result = s.drop('one')
    assert_series_equal(result, expected)

    # single string/tuple-like
    s = Series(range(3), index=list('abc'))
    pytest.raises(KeyError, s.drop, 'bc')
    pytest.raises(KeyError, s.drop, ('a',))

    # errors='ignore'
    s = Series(range(3), index=list('abc'))
    result = s.drop('bc', errors='ignore')
    assert_series_equal(result, s)
    result = s.drop(['a', 'd'], errors='ignore')
    expected = s.iloc[1:]
    assert_series_equal(result, expected)

    # bad axis
    pytest.raises(ValueError, s.drop, 'one', axis='columns')

    # GH 8522
    s = Series([2, 3], index=[True, False])
    assert s.index.is_object()
    result = s.drop(True)
    expected = Series([3], index=[False])
    assert_series_equal(result, expected)

    # GH 16877
    s = Series([2, 3], index=[0, 1])
    with tm.assert_raises_regex(KeyError, 'not contained in axis'):
        s.drop([False, True])
