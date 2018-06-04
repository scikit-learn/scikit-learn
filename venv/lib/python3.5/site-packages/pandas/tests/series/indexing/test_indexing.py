# coding=utf-8
# pylint: disable-msg=E1101,W0612

""" test get/set & misc """

import pytest

from datetime import timedelta

import numpy as np
import pandas as pd

from pandas.core.dtypes.common import is_scalar
from pandas import (Series, DataFrame, MultiIndex,
                    Timestamp, Timedelta, Categorical)
from pandas.tseries.offsets import BDay

from pandas.compat import lrange, range

from pandas.util.testing import (assert_series_equal)
import pandas.util.testing as tm


def test_basic_indexing():
    s = Series(np.random.randn(5), index=['a', 'b', 'a', 'a', 'b'])

    pytest.raises(IndexError, s.__getitem__, 5)
    pytest.raises(IndexError, s.__setitem__, 5, 0)

    pytest.raises(KeyError, s.__getitem__, 'c')

    s = s.sort_index()

    pytest.raises(IndexError, s.__getitem__, 5)
    pytest.raises(IndexError, s.__setitem__, 5, 0)


def test_basic_getitem_with_labels(test_data):
    indices = test_data.ts.index[[5, 10, 15]]

    result = test_data.ts[indices]
    expected = test_data.ts.reindex(indices)
    assert_series_equal(result, expected)

    result = test_data.ts[indices[0]:indices[2]]
    expected = test_data.ts.loc[indices[0]:indices[2]]
    assert_series_equal(result, expected)

    # integer indexes, be careful
    s = Series(np.random.randn(10), index=lrange(0, 20, 2))
    inds = [0, 2, 5, 7, 8]
    arr_inds = np.array([0, 2, 5, 7, 8])
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = s[inds]
    expected = s.reindex(inds)
    assert_series_equal(result, expected)

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = s[arr_inds]
    expected = s.reindex(arr_inds)
    assert_series_equal(result, expected)

    # GH12089
    # with tz for values
    s = Series(pd.date_range("2011-01-01", periods=3, tz="US/Eastern"),
               index=['a', 'b', 'c'])
    expected = Timestamp('2011-01-01', tz='US/Eastern')
    result = s.loc['a']
    assert result == expected
    result = s.iloc[0]
    assert result == expected
    result = s['a']
    assert result == expected


def test_getitem_setitem_ellipsis():
    s = Series(np.random.randn(10))

    np.fix(s)

    result = s[...]
    assert_series_equal(result, s)

    s[...] = 5
    assert (result == 5).all()


def test_getitem_get(test_data):
    test_series = test_data.series
    test_obj_series = test_data.objSeries

    idx1 = test_series.index[5]
    idx2 = test_obj_series.index[5]

    assert test_series[idx1] == test_series.get(idx1)
    assert test_obj_series[idx2] == test_obj_series.get(idx2)

    assert test_series[idx1] == test_series[5]
    assert test_obj_series[idx2] == test_obj_series[5]

    assert test_series.get(-1) == test_series.get(test_series.index[-1])
    assert test_series[5] == test_series.get(test_series.index[5])

    # missing
    d = test_data.ts.index[0] - BDay()
    pytest.raises(KeyError, test_data.ts.__getitem__, d)

    # None
    # GH 5652
    for s in [Series(), Series(index=list('abc'))]:
        result = s.get(None)
        assert result is None


def test_getitem_fancy(test_data):
    slice1 = test_data.series[[1, 2, 3]]
    slice2 = test_data.objSeries[[1, 2, 3]]
    assert test_data.series.index[2] == slice1.index[1]
    assert test_data.objSeries.index[2] == slice2.index[1]
    assert test_data.series[2] == slice1[1]
    assert test_data.objSeries[2] == slice2[1]


def test_getitem_generator(test_data):
    gen = (x > 0 for x in test_data.series)
    result = test_data.series[gen]
    result2 = test_data.series[iter(test_data.series > 0)]
    expected = test_data.series[test_data.series > 0]
    assert_series_equal(result, expected)
    assert_series_equal(result2, expected)


def test_type_promotion():
    # GH12599
    s = pd.Series()
    s["a"] = pd.Timestamp("2016-01-01")
    s["b"] = 3.0
    s["c"] = "foo"
    expected = Series([pd.Timestamp("2016-01-01"), 3.0, "foo"],
                      index=["a", "b", "c"])
    assert_series_equal(s, expected)


@pytest.mark.parametrize(
    'result_1, duplicate_item, expected_1',
    [
        [
            pd.Series({1: 12, 2: [1, 2, 2, 3]}), pd.Series({1: 313}),
            pd.Series({1: 12, }, dtype=object),
        ],
        [
            pd.Series({1: [1, 2, 3], 2: [1, 2, 2, 3]}),
            pd.Series({1: [1, 2, 3]}), pd.Series({1: [1, 2, 3], }),
        ],
    ])
def test_getitem_with_duplicates_indices(
        result_1, duplicate_item, expected_1):
    # GH 17610
    result = result_1.append(duplicate_item)
    expected = expected_1.append(duplicate_item)
    assert_series_equal(result[1], expected)
    assert result[2] == result_1[2]


def test_getitem_out_of_bounds(test_data):
    # don't segfault, GH #495
    pytest.raises(IndexError, test_data.ts.__getitem__, len(test_data.ts))

    # GH #917
    s = Series([])
    pytest.raises(IndexError, s.__getitem__, -1)


def test_getitem_setitem_integers():
    # caused bug without test
    s = Series([1, 2, 3], ['a', 'b', 'c'])

    assert s.iloc[0] == s['a']
    s.iloc[0] = 5
    tm.assert_almost_equal(s['a'], 5)


def test_getitem_box_float64(test_data):
    value = test_data.ts[5]
    assert isinstance(value, np.float64)


def test_series_box_timestamp():
    rng = pd.date_range('20090415', '20090519', freq='B')
    ser = Series(rng)

    assert isinstance(ser[5], pd.Timestamp)

    rng = pd.date_range('20090415', '20090519', freq='B')
    ser = Series(rng, index=rng)
    assert isinstance(ser[5], pd.Timestamp)

    assert isinstance(ser.iat[5], pd.Timestamp)


def test_getitem_ambiguous_keyerror():
    s = Series(lrange(10), index=lrange(0, 20, 2))
    pytest.raises(KeyError, s.__getitem__, 1)
    pytest.raises(KeyError, s.loc.__getitem__, 1)


def test_getitem_unordered_dup():
    obj = Series(lrange(5), index=['c', 'a', 'a', 'b', 'b'])
    assert is_scalar(obj['c'])
    assert obj['c'] == 0


def test_getitem_dups_with_missing():
    # breaks reindex, so need to use .loc internally
    # GH 4246
    s = Series([1, 2, 3, 4], ['foo', 'bar', 'foo', 'bah'])
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        expected = s.loc[['foo', 'bar', 'bah', 'bam']]

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = s[['foo', 'bar', 'bah', 'bam']]
    assert_series_equal(result, expected)


def test_getitem_dups():
    s = Series(range(5), index=['A', 'A', 'B', 'C', 'C'], dtype=np.int64)
    expected = Series([3, 4], index=['C', 'C'], dtype=np.int64)
    result = s['C']
    assert_series_equal(result, expected)


def test_setitem_ambiguous_keyerror():
    s = Series(lrange(10), index=lrange(0, 20, 2))

    # equivalent of an append
    s2 = s.copy()
    s2[1] = 5
    expected = s.append(Series([5], index=[1]))
    assert_series_equal(s2, expected)

    s2 = s.copy()
    s2.loc[1] = 5
    expected = s.append(Series([5], index=[1]))
    assert_series_equal(s2, expected)


def test_getitem_dataframe():
    rng = list(range(10))
    s = pd.Series(10, index=rng)
    df = pd.DataFrame(rng, index=rng)
    pytest.raises(TypeError, s.__getitem__, df > 5)


def test_setitem(test_data):
    test_data.ts[test_data.ts.index[5]] = np.NaN
    test_data.ts[[1, 2, 17]] = np.NaN
    test_data.ts[6] = np.NaN
    assert np.isnan(test_data.ts[6])
    assert np.isnan(test_data.ts[2])
    test_data.ts[np.isnan(test_data.ts)] = 5
    assert not np.isnan(test_data.ts[2])

    # caught this bug when writing tests
    series = Series(tm.makeIntIndex(20).astype(float),
                    index=tm.makeIntIndex(20))

    series[::2] = 0
    assert (series[::2] == 0).all()

    # set item that's not contained
    s = test_data.series.copy()
    s['foobar'] = 1

    app = Series([1], index=['foobar'], name='series')
    expected = test_data.series.append(app)
    assert_series_equal(s, expected)

    # Test for issue #10193
    key = pd.Timestamp('2012-01-01')
    series = pd.Series()
    series[key] = 47
    expected = pd.Series(47, [key])
    assert_series_equal(series, expected)

    series = pd.Series([], pd.DatetimeIndex([], freq='D'))
    series[key] = 47
    expected = pd.Series(47, pd.DatetimeIndex([key], freq='D'))
    assert_series_equal(series, expected)


def test_setitem_dtypes():
    # change dtypes
    # GH 4463
    expected = Series([np.nan, 2, 3])

    s = Series([1, 2, 3])
    s.iloc[0] = np.nan
    assert_series_equal(s, expected)

    s = Series([1, 2, 3])
    s.loc[0] = np.nan
    assert_series_equal(s, expected)

    s = Series([1, 2, 3])
    s[0] = np.nan
    assert_series_equal(s, expected)

    s = Series([False])
    s.loc[0] = np.nan
    assert_series_equal(s, Series([np.nan]))

    s = Series([False, True])
    s.loc[0] = np.nan
    assert_series_equal(s, Series([np.nan, 1.0]))


def test_set_value(test_data):
    idx = test_data.ts.index[10]
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        res = test_data.ts.set_value(idx, 0)
    assert res is test_data.ts
    assert test_data.ts[idx] == 0

    # equiv
    s = test_data.series.copy()
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        res = s.set_value('foobar', 0)
    assert res is s
    assert res.index[-1] == 'foobar'
    assert res['foobar'] == 0

    s = test_data.series.copy()
    s.loc['foobar'] = 0
    assert s.index[-1] == 'foobar'
    assert s['foobar'] == 0


def test_setslice(test_data):
    sl = test_data.ts[5:20]
    assert len(sl) == len(sl.index)
    assert sl.index.is_unique


def test_basic_getitem_setitem_corner(test_data):
    # invalid tuples, e.g. td.ts[:, None] vs. td.ts[:, 2]
    with tm.assert_raises_regex(ValueError, 'tuple-index'):
        test_data.ts[:, 2]
    with tm.assert_raises_regex(ValueError, 'tuple-index'):
        test_data.ts[:, 2] = 2

    # weird lists. [slice(0, 5)] will work but not two slices
    result = test_data.ts[[slice(None, 5)]]
    expected = test_data.ts[:5]
    assert_series_equal(result, expected)

    # OK
    pytest.raises(Exception, test_data.ts.__getitem__,
                  [5, slice(None, None)])
    pytest.raises(Exception, test_data.ts.__setitem__,
                  [5, slice(None, None)], 2)


@pytest.mark.parametrize('tz', ['US/Eastern', 'UTC', 'Asia/Tokyo'])
def test_setitem_with_tz(tz):
    orig = pd.Series(pd.date_range('2016-01-01', freq='H', periods=3,
                                   tz=tz))
    assert orig.dtype == 'datetime64[ns, {0}]'.format(tz)

    # scalar
    s = orig.copy()
    s[1] = pd.Timestamp('2011-01-01', tz=tz)
    exp = pd.Series([pd.Timestamp('2016-01-01 00:00', tz=tz),
                     pd.Timestamp('2011-01-01 00:00', tz=tz),
                     pd.Timestamp('2016-01-01 02:00', tz=tz)])
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[1] = pd.Timestamp('2011-01-01', tz=tz)
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[1] = pd.Timestamp('2011-01-01', tz=tz)
    tm.assert_series_equal(s, exp)

    # vector
    vals = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                      pd.Timestamp('2012-01-01', tz=tz)], index=[1, 2])
    assert vals.dtype == 'datetime64[ns, {0}]'.format(tz)

    s[[1, 2]] = vals
    exp = pd.Series([pd.Timestamp('2016-01-01 00:00', tz=tz),
                     pd.Timestamp('2011-01-01 00:00', tz=tz),
                     pd.Timestamp('2012-01-01 00:00', tz=tz)])
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)


def test_setitem_with_tz_dst():
    # GH XXX
    tz = 'US/Eastern'
    orig = pd.Series(pd.date_range('2016-11-06', freq='H', periods=3,
                                   tz=tz))
    assert orig.dtype == 'datetime64[ns, {0}]'.format(tz)

    # scalar
    s = orig.copy()
    s[1] = pd.Timestamp('2011-01-01', tz=tz)
    exp = pd.Series([pd.Timestamp('2016-11-06 00:00-04:00', tz=tz),
                     pd.Timestamp('2011-01-01 00:00-05:00', tz=tz),
                     pd.Timestamp('2016-11-06 01:00-05:00', tz=tz)])
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[1] = pd.Timestamp('2011-01-01', tz=tz)
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[1] = pd.Timestamp('2011-01-01', tz=tz)
    tm.assert_series_equal(s, exp)

    # vector
    vals = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                      pd.Timestamp('2012-01-01', tz=tz)], index=[1, 2])
    assert vals.dtype == 'datetime64[ns, {0}]'.format(tz)

    s[[1, 2]] = vals
    exp = pd.Series([pd.Timestamp('2016-11-06 00:00', tz=tz),
                     pd.Timestamp('2011-01-01 00:00', tz=tz),
                     pd.Timestamp('2012-01-01 00:00', tz=tz)])
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.loc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.iloc[[1, 2]] = vals
    tm.assert_series_equal(s, exp)


def test_categorial_assigning_ops():
    orig = Series(Categorical(["b", "b"], categories=["a", "b"]))
    s = orig.copy()
    s[:] = "a"
    exp = Series(Categorical(["a", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[1] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[s.index > 0] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s[[False, True]] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
    tm.assert_series_equal(s, exp)

    s = orig.copy()
    s.index = ["x", "y"]
    s["y"] = "a"
    exp = Series(Categorical(["b", "a"], categories=["a", "b"]),
                 index=["x", "y"])
    tm.assert_series_equal(s, exp)

    # ensure that one can set something to np.nan
    s = Series(Categorical([1, 2, 3]))
    exp = Series(Categorical([1, np.nan, 3], categories=[1, 2, 3]))
    s[1] = np.nan
    tm.assert_series_equal(s, exp)


def test_slice(test_data):
    numSlice = test_data.series[10:20]
    numSliceEnd = test_data.series[-10:]
    objSlice = test_data.objSeries[10:20]

    assert test_data.series.index[9] not in numSlice.index
    assert test_data.objSeries.index[9] not in objSlice.index

    assert len(numSlice) == len(numSlice.index)
    assert test_data.series[numSlice.index[0]] == numSlice[numSlice.index[0]]

    assert numSlice.index[1] == test_data.series.index[11]
    assert tm.equalContents(numSliceEnd, np.array(test_data.series)[-10:])

    # Test return view.
    sl = test_data.series[10:20]
    sl[:] = 0

    assert (test_data.series[10:20] == 0).all()


def test_slice_can_reorder_not_uniquely_indexed():
    s = Series(1, index=['a', 'a', 'b', 'b', 'c'])
    s[::-1]  # it works!


def test_ix_setitem(test_data):
    inds = test_data.series.index[[3, 4, 7]]

    result = test_data.series.copy()
    result.loc[inds] = 5

    expected = test_data.series.copy()
    expected[[3, 4, 7]] = 5
    assert_series_equal(result, expected)

    result.iloc[5:10] = 10
    expected[5:10] = 10
    assert_series_equal(result, expected)

    # set slice with indices
    d1, d2 = test_data.series.index[[5, 15]]
    result.loc[d1:d2] = 6
    expected[5:16] = 6  # because it's inclusive
    assert_series_equal(result, expected)

    # set index value
    test_data.series.loc[d1] = 4
    test_data.series.loc[d2] = 6
    assert test_data.series[d1] == 4
    assert test_data.series[d2] == 6


def test_setitem_na():
    # these induce dtype changes
    expected = Series([np.nan, 3, np.nan, 5, np.nan, 7, np.nan, 9, np.nan])
    s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])
    s[::2] = np.nan
    assert_series_equal(s, expected)

    # gets coerced to float, right?
    expected = Series([np.nan, 1, np.nan, 0])
    s = Series([True, True, False, False])
    s[::2] = np.nan
    assert_series_equal(s, expected)

    expected = Series([np.nan, np.nan, np.nan, np.nan, np.nan, 5, 6, 7, 8,
                       9])
    s = Series(np.arange(10))
    s[:5] = np.nan
    assert_series_equal(s, expected)


def test_timedelta_assignment():
    # GH 8209
    s = Series([])
    s.loc['B'] = timedelta(1)
    tm.assert_series_equal(s, Series(Timedelta('1 days'), index=['B']))

    s = s.reindex(s.index.insert(0, 'A'))
    tm.assert_series_equal(s, Series(
        [np.nan, Timedelta('1 days')], index=['A', 'B']))

    result = s.fillna(timedelta(1))
    expected = Series(Timedelta('1 days'), index=['A', 'B'])
    tm.assert_series_equal(result, expected)

    s.loc['A'] = timedelta(1)
    tm.assert_series_equal(s, expected)

    # GH 14155
    s = Series(10 * [np.timedelta64(10, 'm')])
    s.loc[[1, 2, 3]] = np.timedelta64(20, 'm')
    expected = pd.Series(10 * [np.timedelta64(10, 'm')])
    expected.loc[[1, 2, 3]] = pd.Timedelta(np.timedelta64(20, 'm'))
    tm.assert_series_equal(s, expected)


def test_underlying_data_conversion():
    # GH 4080
    df = DataFrame({c: [1, 2, 3] for c in ['a', 'b', 'c']})
    df.set_index(['a', 'b', 'c'], inplace=True)
    s = Series([1], index=[(2, 2, 2)])
    df['val'] = 0
    df
    df['val'].update(s)

    expected = DataFrame(
        dict(a=[1, 2, 3], b=[1, 2, 3], c=[1, 2, 3], val=[0, 1, 0]))
    expected.set_index(['a', 'b', 'c'], inplace=True)
    tm.assert_frame_equal(df, expected)

    # GH 3970
    # these are chained assignments as well
    pd.set_option('chained_assignment', None)
    df = DataFrame({"aa": range(5), "bb": [2.2] * 5})
    df["cc"] = 0.0

    ck = [True] * len(df)

    df["bb"].iloc[0] = .13

    # TODO: unused
    df_tmp = df.iloc[ck]  # noqa

    df["bb"].iloc[0] = .15
    assert df['bb'].iloc[0] == 0.15
    pd.set_option('chained_assignment', 'raise')

    # GH 3217
    df = DataFrame(dict(a=[1, 3], b=[np.nan, 2]))
    df['c'] = np.nan
    df['c'].update(pd.Series(['foo'], index=[0]))

    expected = DataFrame(dict(a=[1, 3], b=[np.nan, 2], c=['foo', np.nan]))
    tm.assert_frame_equal(df, expected)


def test_preserve_refs(test_data):
    seq = test_data.ts[[5, 10, 15]]
    seq[1] = np.NaN
    assert not np.isnan(test_data.ts[10])


def test_cast_on_putmask():
    # GH 2746

    # need to upcast
    s = Series([1, 2], index=[1, 2], dtype='int64')
    s[[True, False]] = Series([0], index=[1], dtype='int64')
    expected = Series([0, 2], index=[1, 2], dtype='int64')

    assert_series_equal(s, expected)


def test_type_promote_putmask():
    # GH8387: test that changing types does not break alignment
    ts = Series(np.random.randn(100), index=np.arange(100, 0, -1)).round(5)
    left, mask = ts.copy(), ts > 0
    right = ts[mask].copy().map(str)
    left[mask] = right
    assert_series_equal(left, ts.map(lambda t: str(t) if t > 0 else t))

    s = Series([0, 1, 2, 0])
    mask = s > 0
    s2 = s[mask].map(str)
    s[mask] = s2
    assert_series_equal(s, Series([0, '1', '2', 0]))

    s = Series([0, 'foo', 'bar', 0])
    mask = Series([False, True, True, False])
    s2 = s[mask]
    s[mask] = s2
    assert_series_equal(s, Series([0, 'foo', 'bar', 0]))


def test_multilevel_preserve_name():
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                              'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['first', 'second'])
    s = Series(np.random.randn(len(index)), index=index, name='sth')

    result = s['foo']
    result2 = s.loc['foo']
    assert result.name == s.name
    assert result2.name == s.name


def test_setitem_scalar_into_readonly_backing_data():
    # GH14359: test that you cannot mutate a read only buffer

    array = np.zeros(5)
    array.flags.writeable = False  # make the array immutable
    series = Series(array)

    for n in range(len(series)):
        with pytest.raises(ValueError):
            series[n] = 1

        assert array[n] == 0


def test_setitem_slice_into_readonly_backing_data():
    # GH14359: test that you cannot mutate a read only buffer

    array = np.zeros(5)
    array.flags.writeable = False  # make the array immutable
    series = Series(array)

    with pytest.raises(ValueError):
        series[1:3] = 1

    assert not array.any()


"""
miscellaneous methods
"""


def test_select(test_data):
    # deprecated: gh-12410
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        n = len(test_data.ts)
        result = test_data.ts.select(lambda x: x >= test_data.ts.index[n // 2])
        expected = test_data.ts.reindex(test_data.ts.index[n // 2:])
        assert_series_equal(result, expected)

        result = test_data.ts.select(lambda x: x.weekday() == 2)
        expected = test_data.ts[test_data.ts.index.weekday == 2]
        assert_series_equal(result, expected)


def test_pop():
    # GH 6600
    df = DataFrame({'A': 0, 'B': np.arange(5, dtype='int64'), 'C': 0, })
    k = df.iloc[4]

    result = k.pop('B')
    assert result == 4

    expected = Series([0, 0], index=['A', 'C'], name=4)
    assert_series_equal(k, expected)


def test_take():
    s = Series([-1, 5, 6, 2, 4])

    actual = s.take([1, 3, 4])
    expected = Series([5, 2, 4], index=[1, 3, 4])
    tm.assert_series_equal(actual, expected)

    actual = s.take([-1, 3, 4])
    expected = Series([4, 2, 4], index=[4, 3, 4])
    tm.assert_series_equal(actual, expected)

    pytest.raises(IndexError, s.take, [1, 10])
    pytest.raises(IndexError, s.take, [2, 5])

    with tm.assert_produces_warning(FutureWarning):
        s.take([-1, 3, 4], convert=False)


def test_take_categorical():
    # https://github.com/pandas-dev/pandas/issues/20664
    s = Series(pd.Categorical(['a', 'b', 'c']))
    result = s.take([-2, -2, 0])
    expected = Series(pd.Categorical(['b', 'b', 'a'],
                      categories=['a', 'b', 'c']),
                      index=[1, 1, 0])
    assert_series_equal(result, expected)


def test_head_tail(test_data):
    assert_series_equal(test_data.series.head(), test_data.series[:5])
    assert_series_equal(test_data.series.head(0), test_data.series[0:0])
    assert_series_equal(test_data.series.tail(), test_data.series[-5:])
    assert_series_equal(test_data.series.tail(0), test_data.series[0:0])
