# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

import pandas as pd
import numpy as np

from pandas import (Series, date_range, isna, Index, Timestamp)
from pandas.compat import lrange, range
from pandas.core.dtypes.common import is_integer

from pandas.core.indexing import IndexingError
from pandas.tseries.offsets import BDay

from pandas.util.testing import (assert_series_equal)
import pandas.util.testing as tm


def test_getitem_boolean(test_data):
    s = test_data.series
    mask = s > s.median()

    # passing list is OK
    result = s[list(mask)]
    expected = s[mask]
    assert_series_equal(result, expected)
    tm.assert_index_equal(result.index, s.index[mask])


def test_getitem_boolean_empty():
    s = Series([], dtype=np.int64)
    s.index.name = 'index_name'
    s = s[s.isna()]
    assert s.index.name == 'index_name'
    assert s.dtype == np.int64

    # GH5877
    # indexing with empty series
    s = Series(['A', 'B'])
    expected = Series(np.nan, index=['C'], dtype=object)
    result = s[Series(['C'], dtype=object)]
    assert_series_equal(result, expected)

    s = Series(['A', 'B'])
    expected = Series(dtype=object, index=Index([], dtype='int64'))
    result = s[Series([], dtype=object)]
    assert_series_equal(result, expected)

    # invalid because of the boolean indexer
    # that's empty or not-aligned
    def f():
        s[Series([], dtype=bool)]

    pytest.raises(IndexingError, f)

    def f():
        s[Series([True], dtype=bool)]

    pytest.raises(IndexingError, f)


def test_getitem_boolean_object(test_data):
    # using column from DataFrame

    s = test_data.series
    mask = s > s.median()
    omask = mask.astype(object)

    # getitem
    result = s[omask]
    expected = s[mask]
    assert_series_equal(result, expected)

    # setitem
    s2 = s.copy()
    cop = s.copy()
    cop[omask] = 5
    s2[mask] = 5
    assert_series_equal(cop, s2)

    # nans raise exception
    omask[5:10] = np.nan
    pytest.raises(Exception, s.__getitem__, omask)
    pytest.raises(Exception, s.__setitem__, omask, 5)


def test_getitem_setitem_boolean_corner(test_data):
    ts = test_data.ts
    mask_shifted = ts.shift(1, freq=BDay()) > ts.median()

    # these used to raise...??

    pytest.raises(Exception, ts.__getitem__, mask_shifted)
    pytest.raises(Exception, ts.__setitem__, mask_shifted, 1)
    # ts[mask_shifted]
    # ts[mask_shifted] = 1

    pytest.raises(Exception, ts.loc.__getitem__, mask_shifted)
    pytest.raises(Exception, ts.loc.__setitem__, mask_shifted, 1)
    # ts.loc[mask_shifted]
    # ts.loc[mask_shifted] = 2


def test_setitem_boolean(test_data):
    mask = test_data.series > test_data.series.median()

    # similar indexed series
    result = test_data.series.copy()
    result[mask] = test_data.series * 2
    expected = test_data.series * 2
    assert_series_equal(result[mask], expected[mask])

    # needs alignment
    result = test_data.series.copy()
    result[mask] = (test_data.series * 2)[0:5]
    expected = (test_data.series * 2)[0:5].reindex_like(test_data.series)
    expected[-mask] = test_data.series[mask]
    assert_series_equal(result[mask], expected[mask])


def test_get_set_boolean_different_order(test_data):
    ordered = test_data.series.sort_values()

    # setting
    copy = test_data.series.copy()
    copy[ordered > 0] = 0

    expected = test_data.series.copy()
    expected[expected > 0] = 0

    assert_series_equal(copy, expected)

    # getting
    sel = test_data.series[ordered > 0]
    exp = test_data.series[test_data.series > 0]
    assert_series_equal(sel, exp)


def test_where_unsafe():
    # unsafe dtype changes
    for dtype in [np.int8, np.int16, np.int32, np.int64, np.float16,
                  np.float32, np.float64]:
        s = Series(np.arange(10), dtype=dtype)
        mask = s < 5
        s[mask] = lrange(2, 7)
        expected = Series(lrange(2, 7) + lrange(5, 10), dtype=dtype)
        assert_series_equal(s, expected)
        assert s.dtype == expected.dtype

    # these are allowed operations, but are upcasted
    for dtype in [np.int64, np.float64]:
        s = Series(np.arange(10), dtype=dtype)
        mask = s < 5
        values = [2.5, 3.5, 4.5, 5.5, 6.5]
        s[mask] = values
        expected = Series(values + lrange(5, 10), dtype='float64')
        assert_series_equal(s, expected)
        assert s.dtype == expected.dtype

    # GH 9731
    s = Series(np.arange(10), dtype='int64')
    mask = s > 5
    values = [2.5, 3.5, 4.5, 5.5]
    s[mask] = values
    expected = Series(lrange(6) + values, dtype='float64')
    assert_series_equal(s, expected)

    # can't do these as we are forced to change the itemsize of the input
    # to something we cannot
    for dtype in [np.int8, np.int16, np.int32, np.float16, np.float32]:
        s = Series(np.arange(10), dtype=dtype)
        mask = s < 5
        values = [2.5, 3.5, 4.5, 5.5, 6.5]
        pytest.raises(Exception, s.__setitem__, tuple(mask), values)

    # GH3235
    s = Series(np.arange(10), dtype='int64')
    mask = s < 5
    s[mask] = lrange(2, 7)
    expected = Series(lrange(2, 7) + lrange(5, 10), dtype='int64')
    assert_series_equal(s, expected)
    assert s.dtype == expected.dtype

    s = Series(np.arange(10), dtype='int64')
    mask = s > 5
    s[mask] = [0] * 4
    expected = Series([0, 1, 2, 3, 4, 5] + [0] * 4, dtype='int64')
    assert_series_equal(s, expected)

    s = Series(np.arange(10))
    mask = s > 5

    def f():
        s[mask] = [5, 4, 3, 2, 1]

    pytest.raises(ValueError, f)

    def f():
        s[mask] = [0] * 5

    pytest.raises(ValueError, f)

    # dtype changes
    s = Series([1, 2, 3, 4])
    result = s.where(s > 2, np.nan)
    expected = Series([np.nan, np.nan, 3, 4])
    assert_series_equal(result, expected)

    # GH 4667
    # setting with None changes dtype
    s = Series(range(10)).astype(float)
    s[8] = None
    result = s[8]
    assert isna(result)

    s = Series(range(10)).astype(float)
    s[s > 8] = None
    result = s[isna(s)]
    expected = Series(np.nan, index=[9])
    assert_series_equal(result, expected)


def test_where_raise_on_error_deprecation():
    # gh-14968
    # deprecation of raise_on_error
    s = Series(np.random.randn(5))
    cond = s > 0
    with tm.assert_produces_warning(FutureWarning):
        s.where(cond, raise_on_error=True)
    with tm.assert_produces_warning(FutureWarning):
        s.mask(cond, raise_on_error=True)


def test_where():
    s = Series(np.random.randn(5))
    cond = s > 0

    rs = s.where(cond).dropna()
    rs2 = s[cond]
    assert_series_equal(rs, rs2)

    rs = s.where(cond, -s)
    assert_series_equal(rs, s.abs())

    rs = s.where(cond)
    assert (s.shape == rs.shape)
    assert (rs is not s)

    # test alignment
    cond = Series([True, False, False, True, False], index=s.index)
    s2 = -(s.abs())

    expected = s2[cond].reindex(s2.index[:3]).reindex(s2.index)
    rs = s2.where(cond[:3])
    assert_series_equal(rs, expected)

    expected = s2.abs()
    expected.iloc[0] = s2[0]
    rs = s2.where(cond[:3], -s2)
    assert_series_equal(rs, expected)


def test_where_error():
    s = Series(np.random.randn(5))
    cond = s > 0

    pytest.raises(ValueError, s.where, 1)
    pytest.raises(ValueError, s.where, cond[:3].values, -s)

    # GH 2745
    s = Series([1, 2])
    s[[True, False]] = [0, 1]
    expected = Series([0, 2])
    assert_series_equal(s, expected)

    # failures
    pytest.raises(ValueError, s.__setitem__, tuple([[[True, False]]]),
                  [0, 2, 3])
    pytest.raises(ValueError, s.__setitem__, tuple([[[True, False]]]),
                  [])


@pytest.mark.parametrize('klass', [list, tuple, np.array, Series])
def test_where_array_like(klass):
    # see gh-15414
    s = Series([1, 2, 3])
    cond = [False, True, True]
    expected = Series([np.nan, 2, 3])

    result = s.where(klass(cond))
    assert_series_equal(result, expected)


@pytest.mark.parametrize('cond', [
    [1, 0, 1],
    Series([2, 5, 7]),
    ["True", "False", "True"],
    [Timestamp("2017-01-01"), pd.NaT, Timestamp("2017-01-02")]
])
def test_where_invalid_input(cond):
    # see gh-15414: only boolean arrays accepted
    s = Series([1, 2, 3])
    msg = "Boolean array expected for the condition"

    with tm.assert_raises_regex(ValueError, msg):
        s.where(cond)

    msg = "Array conditional must be same shape as self"
    with tm.assert_raises_regex(ValueError, msg):
        s.where([True])


def test_where_ndframe_align():
    msg = "Array conditional must be same shape as self"
    s = Series([1, 2, 3])

    cond = [True]
    with tm.assert_raises_regex(ValueError, msg):
        s.where(cond)

    expected = Series([1, np.nan, np.nan])

    out = s.where(Series(cond))
    tm.assert_series_equal(out, expected)

    cond = np.array([False, True, False, True])
    with tm.assert_raises_regex(ValueError, msg):
        s.where(cond)

    expected = Series([np.nan, 2, np.nan])

    out = s.where(Series(cond))
    tm.assert_series_equal(out, expected)


def test_where_setitem_invalid():
    # GH 2702
    # make sure correct exceptions are raised on invalid list assignment

    # slice
    s = Series(list('abc'))

    def f():
        s[0:3] = list(range(27))

    pytest.raises(ValueError, f)

    s[0:3] = list(range(3))
    expected = Series([0, 1, 2])
    assert_series_equal(s.astype(np.int64), expected, )

    # slice with step
    s = Series(list('abcdef'))

    def f():
        s[0:4:2] = list(range(27))

    pytest.raises(ValueError, f)

    s = Series(list('abcdef'))
    s[0:4:2] = list(range(2))
    expected = Series([0, 'b', 1, 'd', 'e', 'f'])
    assert_series_equal(s, expected)

    # neg slices
    s = Series(list('abcdef'))

    def f():
        s[:-1] = list(range(27))

    pytest.raises(ValueError, f)

    s[-3:-1] = list(range(2))
    expected = Series(['a', 'b', 'c', 0, 1, 'f'])
    assert_series_equal(s, expected)

    # list
    s = Series(list('abc'))

    def f():
        s[[0, 1, 2]] = list(range(27))

    pytest.raises(ValueError, f)

    s = Series(list('abc'))

    def f():
        s[[0, 1, 2]] = list(range(2))

    pytest.raises(ValueError, f)

    # scalar
    s = Series(list('abc'))
    s[0] = list(range(10))
    expected = Series([list(range(10)), 'b', 'c'])
    assert_series_equal(s, expected)


@pytest.mark.parametrize('size', range(2, 6))
@pytest.mark.parametrize('mask', [
    [True, False, False, False, False],
    [True, False],
    [False]
])
@pytest.mark.parametrize('item', [
    2.0, np.nan, np.finfo(np.float).max, np.finfo(np.float).min
])
# Test numpy arrays, lists and tuples as the input to be
# broadcast
@pytest.mark.parametrize('box', [
    lambda x: np.array([x]),
    lambda x: [x],
    lambda x: (x,)
])
def test_broadcast(size, mask, item, box):
    selection = np.resize(mask, size)

    data = np.arange(size, dtype=float)

    # Construct the expected series by taking the source
    # data or item based on the selection
    expected = Series([item if use_item else data[
        i] for i, use_item in enumerate(selection)])

    s = Series(data)
    s[selection] = box(item)
    assert_series_equal(s, expected)

    s = Series(data)
    result = s.where(~selection, box(item))
    assert_series_equal(result, expected)

    s = Series(data)
    result = s.mask(selection, box(item))
    assert_series_equal(result, expected)


def test_where_inplace():
    s = Series(np.random.randn(5))
    cond = s > 0

    rs = s.copy()

    rs.where(cond, inplace=True)
    assert_series_equal(rs.dropna(), s[cond])
    assert_series_equal(rs, s.where(cond))

    rs = s.copy()
    rs.where(cond, -s, inplace=True)
    assert_series_equal(rs, s.where(cond, -s))


def test_where_dups():
    # GH 4550
    # where crashes with dups in index
    s1 = Series(list(range(3)))
    s2 = Series(list(range(3)))
    comb = pd.concat([s1, s2])
    result = comb.where(comb < 2)
    expected = Series([0, 1, np.nan, 0, 1, np.nan],
                      index=[0, 1, 2, 0, 1, 2])
    assert_series_equal(result, expected)

    # GH 4548
    # inplace updating not working with dups
    comb[comb < 1] = 5
    expected = Series([5, 1, 2, 5, 1, 2], index=[0, 1, 2, 0, 1, 2])
    assert_series_equal(comb, expected)

    comb[comb < 2] += 10
    expected = Series([5, 11, 2, 5, 11, 2], index=[0, 1, 2, 0, 1, 2])
    assert_series_equal(comb, expected)


def test_where_numeric_with_string():
    # GH 9280
    s = pd.Series([1, 2, 3])
    w = s.where(s > 1, 'X')

    assert not is_integer(w[0])
    assert is_integer(w[1])
    assert is_integer(w[2])
    assert isinstance(w[0], str)
    assert w.dtype == 'object'

    w = s.where(s > 1, ['X', 'Y', 'Z'])
    assert not is_integer(w[0])
    assert is_integer(w[1])
    assert is_integer(w[2])
    assert isinstance(w[0], str)
    assert w.dtype == 'object'

    w = s.where(s > 1, np.array(['X', 'Y', 'Z']))
    assert not is_integer(w[0])
    assert is_integer(w[1])
    assert is_integer(w[2])
    assert isinstance(w[0], str)
    assert w.dtype == 'object'


def test_where_timedelta_coerce():
    s = Series([1, 2], dtype='timedelta64[ns]')
    expected = Series([10, 10])
    mask = np.array([False, False])

    rs = s.where(mask, [10, 10])
    assert_series_equal(rs, expected)

    rs = s.where(mask, 10)
    assert_series_equal(rs, expected)

    rs = s.where(mask, 10.0)
    assert_series_equal(rs, expected)

    rs = s.where(mask, [10.0, 10.0])
    assert_series_equal(rs, expected)

    rs = s.where(mask, [10.0, np.nan])
    expected = Series([10, None], dtype='object')
    assert_series_equal(rs, expected)


def test_where_datetime_conversion():
    s = Series(date_range('20130102', periods=2))
    expected = Series([10, 10])
    mask = np.array([False, False])

    rs = s.where(mask, [10, 10])
    assert_series_equal(rs, expected)

    rs = s.where(mask, 10)
    assert_series_equal(rs, expected)

    rs = s.where(mask, 10.0)
    assert_series_equal(rs, expected)

    rs = s.where(mask, [10.0, 10.0])
    assert_series_equal(rs, expected)

    rs = s.where(mask, [10.0, np.nan])
    expected = Series([10, None], dtype='object')
    assert_series_equal(rs, expected)

    # GH 15701
    timestamps = ['2016-12-31 12:00:04+00:00',
                  '2016-12-31 12:00:04.010000+00:00']
    s = Series([pd.Timestamp(t) for t in timestamps])
    rs = s.where(Series([False, True]))
    expected = Series([pd.NaT, s[1]])
    assert_series_equal(rs, expected)


def test_mask():
    # compare with tested results in test_where
    s = Series(np.random.randn(5))
    cond = s > 0

    rs = s.where(~cond, np.nan)
    assert_series_equal(rs, s.mask(cond))

    rs = s.where(~cond)
    rs2 = s.mask(cond)
    assert_series_equal(rs, rs2)

    rs = s.where(~cond, -s)
    rs2 = s.mask(cond, -s)
    assert_series_equal(rs, rs2)

    cond = Series([True, False, False, True, False], index=s.index)
    s2 = -(s.abs())
    rs = s2.where(~cond[:3])
    rs2 = s2.mask(cond[:3])
    assert_series_equal(rs, rs2)

    rs = s2.where(~cond[:3], -s2)
    rs2 = s2.mask(cond[:3], -s2)
    assert_series_equal(rs, rs2)

    pytest.raises(ValueError, s.mask, 1)
    pytest.raises(ValueError, s.mask, cond[:3].values, -s)

    # dtype changes
    s = Series([1, 2, 3, 4])
    result = s.mask(s > 2, np.nan)
    expected = Series([1, 2, np.nan, np.nan])
    assert_series_equal(result, expected)


def test_mask_inplace():
    s = Series(np.random.randn(5))
    cond = s > 0

    rs = s.copy()
    rs.mask(cond, inplace=True)
    assert_series_equal(rs.dropna(), s[~cond])
    assert_series_equal(rs, s.mask(cond))

    rs = s.copy()
    rs.mask(cond, -s, inplace=True)
    assert_series_equal(rs, s.mask(cond, -s))
