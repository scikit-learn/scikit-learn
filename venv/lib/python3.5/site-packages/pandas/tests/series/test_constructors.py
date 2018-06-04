# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

from datetime import datetime, timedelta
from collections import OrderedDict

from numpy import nan
import numpy as np
import numpy.ma as ma
import pandas as pd

from pandas.api.types import CategoricalDtype
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_datetime64tz_dtype)
from pandas import (Index, Series, isna, date_range, Timestamp,
                    NaT, period_range, timedelta_range, MultiIndex,
                    IntervalIndex, Categorical, DataFrame)

from pandas._libs import lib
from pandas._libs.tslib import iNaT

from pandas.compat import lrange, range, zip, long, PY36
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesConstructors(TestData):

    def test_invalid_dtype(self):
        # GH15520
        msg = 'not understood'
        invalid_list = [pd.Timestamp, 'pd.Timestamp', list]
        for dtype in invalid_list:
            with tm.assert_raises_regex(TypeError, msg):
                Series([], name='time', dtype=dtype)

    def test_scalar_conversion(self):

        # Pass in scalar is disabled
        scalar = Series(0.5)
        assert not isinstance(scalar, float)

        # Coercion
        assert float(Series([1.])) == 1.0
        assert int(Series([1.])) == 1
        assert long(Series([1.])) == 1

    def test_constructor(self):
        assert self.ts.index.is_all_dates

        # Pass in Series
        derived = Series(self.ts)
        assert derived.index.is_all_dates

        assert tm.equalContents(derived.index, self.ts.index)
        # Ensure new index is not created
        assert id(self.ts.index) == id(derived.index)

        # Mixed type Series
        mixed = Series(['hello', np.NaN], index=[0, 1])
        assert mixed.dtype == np.object_
        assert mixed[1] is np.NaN

        assert not self.empty.index.is_all_dates
        assert not Series({}).index.is_all_dates
        pytest.raises(Exception, Series, np.random.randn(3, 3),
                      index=np.arange(3))

        mixed.name = 'Series'
        rs = Series(mixed).name
        xp = 'Series'
        assert rs == xp

        # raise on MultiIndex GH4187
        m = MultiIndex.from_arrays([[1, 2], [3, 4]])
        pytest.raises(NotImplementedError, Series, m)

    @pytest.mark.parametrize('input_class', [list, dict, OrderedDict])
    def test_constructor_empty(self, input_class):
        empty = Series()
        empty2 = Series(input_class())

        # these are Index() and RangeIndex() which don't compare type equal
        # but are just .equals
        assert_series_equal(empty, empty2, check_index_type=False)

        # With explicit dtype:
        empty = Series(dtype='float64')
        empty2 = Series(input_class(), dtype='float64')
        assert_series_equal(empty, empty2, check_index_type=False)

        # GH 18515 : with dtype=category:
        empty = Series(dtype='category')
        empty2 = Series(input_class(), dtype='category')
        assert_series_equal(empty, empty2, check_index_type=False)

        if input_class is not list:
            # With index:
            empty = Series(index=lrange(10))
            empty2 = Series(input_class(), index=lrange(10))
            assert_series_equal(empty, empty2)

            # With index and dtype float64:
            empty = Series(np.nan, index=lrange(10))
            empty2 = Series(input_class(), index=lrange(10), dtype='float64')
            assert_series_equal(empty, empty2)

            # GH 19853 : with empty string, index and dtype str
            empty = Series('', dtype=str, index=range(3))
            empty2 = Series('', index=range(3))
            assert_series_equal(empty, empty2)

    @pytest.mark.parametrize('input_arg', [np.nan, float('nan')])
    def test_constructor_nan(self, input_arg):
        empty = Series(dtype='float64', index=lrange(10))
        empty2 = Series(input_arg, index=lrange(10))

        assert_series_equal(empty, empty2, check_index_type=False)

    @pytest.mark.parametrize('dtype', [
        'f8', 'i8', 'M8[ns]', 'm8[ns]', 'category', 'object',
        'datetime64[ns, UTC]',
    ])
    @pytest.mark.parametrize('index', [None, pd.Index([])])
    def test_constructor_dtype_only(self, dtype, index):
        # GH-20865
        result = pd.Series(dtype=dtype, index=index)
        assert result.dtype == dtype
        assert len(result) == 0

    def test_constructor_no_data_index_order(self):
        result = pd.Series(index=['b', 'a', 'c'])
        assert result.index.tolist() == ['b', 'a', 'c']

    def test_constructor_series(self):
        index1 = ['d', 'b', 'a', 'c']
        index2 = sorted(index1)
        s1 = Series([4, 7, -5, 3], index=index1)
        s2 = Series(s1, index=index2)

        assert_series_equal(s2, s1.sort_index())

    def test_constructor_iterator(self):

        expected = Series(list(range(10)), dtype='int64')
        result = Series(range(10), dtype='int64')
        assert_series_equal(result, expected)

    def test_constructor_list_like(self):

        # make sure that we are coercing different
        # list-likes to standard dtypes and not
        # platform specific
        expected = Series([1, 2, 3], dtype='int64')
        for obj in [[1, 2, 3], (1, 2, 3),
                    np.array([1, 2, 3], dtype='int64')]:
            result = Series(obj, index=[0, 1, 2])
            assert_series_equal(result, expected)

    @pytest.mark.parametrize('input_vals', [
        ([1, 2]),
        ([1.0, 2.0, np.nan]),
        (['1', '2']),
        (list(pd.date_range('1/1/2011', periods=2, freq='H'))),
        (list(pd.date_range('1/1/2011', periods=2, freq='H',
                            tz='US/Eastern'))),
        ([pd.Interval(left=0, right=5)]),
    ])
    def test_constructor_list_str(self, input_vals):
        # GH 16605
        # Ensure that data elements from a list are converted to strings
        # when dtype is str, 'str', or 'U'

        for dtype in ['str', str, 'U']:
            result = Series(input_vals, dtype=dtype)
            expected = Series(input_vals).astype(dtype)
            assert_series_equal(result, expected)

    def test_constructor_generator(self):
        gen = (i for i in range(10))

        result = Series(gen)
        exp = Series(lrange(10))
        assert_series_equal(result, exp)

        gen = (i for i in range(10))
        result = Series(gen, index=lrange(10, 20))
        exp.index = lrange(10, 20)
        assert_series_equal(result, exp)

    def test_constructor_map(self):
        # GH8909
        m = map(lambda x: x, range(10))

        result = Series(m)
        exp = Series(lrange(10))
        assert_series_equal(result, exp)

        m = map(lambda x: x, range(10))
        result = Series(m, index=lrange(10, 20))
        exp.index = lrange(10, 20)
        assert_series_equal(result, exp)

    def test_constructor_categorical(self):
        cat = pd.Categorical([0, 1, 2, 0, 1, 2], ['a', 'b', 'c'],
                             fastpath=True)
        res = Series(cat)
        tm.assert_categorical_equal(res.values, cat)

        # GH12574
        pytest.raises(
            ValueError, lambda: Series(pd.Categorical([1, 2, 3]),
                                       dtype='int64'))
        cat = Series(pd.Categorical([1, 2, 3]), dtype='category')
        assert is_categorical_dtype(cat)
        assert is_categorical_dtype(cat.dtype)
        s = Series([1, 2, 3], dtype='category')
        assert is_categorical_dtype(s)
        assert is_categorical_dtype(s.dtype)

    def test_constructor_categorical_with_coercion(self):
        factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        # test basic creation / coercion of categoricals
        s = Series(factor, name='A')
        assert s.dtype == 'category'
        assert len(s) == len(factor)
        str(s.values)
        str(s)

        # in a frame
        df = DataFrame({'A': factor})
        result = df['A']
        tm.assert_series_equal(result, s)
        result = df.iloc[:, 0]
        tm.assert_series_equal(result, s)
        assert len(df) == len(factor)
        str(df.values)
        str(df)

        df = DataFrame({'A': s})
        result = df['A']
        tm.assert_series_equal(result, s)
        assert len(df) == len(factor)
        str(df.values)
        str(df)

        # multiples
        df = DataFrame({'A': s, 'B': s, 'C': 1})
        result1 = df['A']
        result2 = df['B']
        tm.assert_series_equal(result1, s)
        tm.assert_series_equal(result2, s, check_names=False)
        assert result2.name == 'B'
        assert len(df) == len(factor)
        str(df.values)
        str(df)

        # GH8623
        x = DataFrame([[1, 'John P. Doe'], [2, 'Jane Dove'],
                       [1, 'John P. Doe']],
                      columns=['person_id', 'person_name'])
        x['person_name'] = Categorical(x.person_name
                                       )  # doing this breaks transform

        expected = x.iloc[0].person_name
        result = x.person_name.iloc[0]
        assert result == expected

        result = x.person_name[0]
        assert result == expected

        result = x.person_name.loc[0]
        assert result == expected

    def test_constructor_categorical_dtype(self):
        result = pd.Series(['a', 'b'],
                           dtype=CategoricalDtype(['a', 'b', 'c'],
                                                  ordered=True))
        assert is_categorical_dtype(result) is True
        tm.assert_index_equal(result.cat.categories, pd.Index(['a', 'b', 'c']))
        assert result.cat.ordered

        result = pd.Series(['a', 'b'], dtype=CategoricalDtype(['b', 'a']))
        assert is_categorical_dtype(result)
        tm.assert_index_equal(result.cat.categories, pd.Index(['b', 'a']))
        assert result.cat.ordered is False

        # GH 19565 - Check broadcasting of scalar with Categorical dtype
        result = Series('a', index=[0, 1],
                        dtype=CategoricalDtype(['a', 'b'], ordered=True))
        expected = Series(['a', 'a'], index=[0, 1],
                          dtype=CategoricalDtype(['a', 'b'], ordered=True))
        tm.assert_series_equal(result, expected, check_categorical=True)

    def test_categorical_sideeffects_free(self):
        # Passing a categorical to a Series and then changing values in either
        # the series or the categorical should not change the values in the
        # other one, IF you specify copy!
        cat = Categorical(["a", "b", "c", "a"])
        s = Series(cat, copy=True)
        assert s.cat is not cat
        s.cat.categories = [1, 2, 3]
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        exp_cat = np.array(["a", "b", "c", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # setting
        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # however, copy is False by default
        # so this WILL change values
        cat = Categorical(["a", "b", "c", "a"])
        s = Series(cat)
        assert s.values is cat
        s.cat.categories = [1, 2, 3]
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)
        tm.assert_numpy_array_equal(cat.__array__(), exp_s)

        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)
        tm.assert_numpy_array_equal(cat.__array__(), exp_s2)

    def test_unordered_compare_equal(self):
        left = pd.Series(['a', 'b', 'c'],
                         dtype=CategoricalDtype(['a', 'b']))
        right = pd.Series(pd.Categorical(['a', 'b', np.nan],
                                         categories=['a', 'b']))
        tm.assert_series_equal(left, right)

    def test_constructor_maskedarray(self):
        data = ma.masked_all((3, ), dtype=float)
        result = Series(data)
        expected = Series([nan, nan, nan])
        assert_series_equal(result, expected)

        data[0] = 0.0
        data[2] = 2.0
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([0.0, nan, 2.0], index=index)
        assert_series_equal(result, expected)

        data[1] = 1.0
        result = Series(data, index=index)
        expected = Series([0.0, 1.0, 2.0], index=index)
        assert_series_equal(result, expected)

        data = ma.masked_all((3, ), dtype=int)
        result = Series(data)
        expected = Series([nan, nan, nan], dtype=float)
        assert_series_equal(result, expected)

        data[0] = 0
        data[2] = 2
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([0, nan, 2], index=index, dtype=float)
        assert_series_equal(result, expected)

        data[1] = 1
        result = Series(data, index=index)
        expected = Series([0, 1, 2], index=index, dtype=int)
        assert_series_equal(result, expected)

        data = ma.masked_all((3, ), dtype=bool)
        result = Series(data)
        expected = Series([nan, nan, nan], dtype=object)
        assert_series_equal(result, expected)

        data[0] = True
        data[2] = False
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([True, nan, False], index=index, dtype=object)
        assert_series_equal(result, expected)

        data[1] = True
        result = Series(data, index=index)
        expected = Series([True, True, False], index=index, dtype=bool)
        assert_series_equal(result, expected)

        data = ma.masked_all((3, ), dtype='M8[ns]')
        result = Series(data)
        expected = Series([iNaT, iNaT, iNaT], dtype='M8[ns]')
        assert_series_equal(result, expected)

        data[0] = datetime(2001, 1, 1)
        data[2] = datetime(2001, 1, 3)
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([datetime(2001, 1, 1), iNaT,
                           datetime(2001, 1, 3)], index=index, dtype='M8[ns]')
        assert_series_equal(result, expected)

        data[1] = datetime(2001, 1, 2)
        result = Series(data, index=index)
        expected = Series([datetime(2001, 1, 1), datetime(2001, 1, 2),
                           datetime(2001, 1, 3)], index=index, dtype='M8[ns]')
        assert_series_equal(result, expected)

    def test_series_ctor_plus_datetimeindex(self):
        rng = date_range('20090415', '20090519', freq='B')
        data = {k: 1 for k in rng}

        result = Series(data, index=rng)
        assert result.index is rng

    def test_constructor_default_index(self):
        s = Series([0, 1, 2])
        tm.assert_index_equal(s.index, pd.Index(np.arange(3)))

    @pytest.mark.parametrize('input', [[1, 2, 3],
                                       (1, 2, 3),
                                       list(range(3)),
                                       pd.Categorical(['a', 'b', 'a']),
                                       (i for i in range(3)),
                                       map(lambda x: x, range(3))])
    def test_constructor_index_mismatch(self, input):
        # GH 19342
        # test that construction of a Series with an index of different length
        # raises an error
        msg = 'Length of passed values is 3, index implies 4'
        with pytest.raises(ValueError, message=msg):
            Series(input, index=np.arange(4))

    def test_constructor_numpy_scalar(self):
        # GH 19342
        # construction with a numpy scalar
        # should not raise
        result = Series(np.array(100), index=np.arange(4), dtype='int64')
        expected = Series(100, index=np.arange(4), dtype='int64')
        tm.assert_series_equal(result, expected)

    def test_constructor_broadcast_list(self):
        # GH 19342
        # construction with single-element container and index
        # should raise
        pytest.raises(ValueError, Series, ['foo'], index=['a', 'b', 'c'])

    def test_constructor_corner(self):
        df = tm.makeTimeDataFrame()
        objs = [df, df]
        s = Series(objs, index=[0, 1])
        assert isinstance(s, Series)

    def test_constructor_sanitize(self):
        s = Series(np.array([1., 1., 8.]), dtype='i8')
        assert s.dtype == np.dtype('i8')

        s = Series(np.array([1., 1., np.nan]), copy=True, dtype='i8')
        assert s.dtype == np.dtype('f8')

    def test_constructor_copy(self):
        # GH15125
        # test dtype parameter has no side effects on copy=True
        for data in [[1.], np.array([1.])]:
            x = Series(data)
            y = pd.Series(x, copy=True, dtype=float)

            # copy=True maintains original data in Series
            tm.assert_series_equal(x, y)

            # changes to origin of copy does not affect the copy
            x[0] = 2.
            assert not x.equals(y)
            assert x[0] == 2.
            assert y[0] == 1.

    @pytest.mark.parametrize(
        "index",
        [
            pd.date_range('20170101', periods=3, tz='US/Eastern'),
            pd.date_range('20170101', periods=3),
            pd.timedelta_range('1 day', periods=3),
            pd.period_range('2012Q1', periods=3, freq='Q'),
            pd.Index(list('abc')),
            pd.Int64Index([1, 2, 3]),
            pd.RangeIndex(0, 3)],
        ids=lambda x: type(x).__name__)
    def test_constructor_limit_copies(self, index):
        # GH 17449
        # limit copies of input
        s = pd.Series(index)

        # we make 1 copy; this is just a smoke test here
        assert s._data.blocks[0].values is not index

    def test_constructor_pass_none(self):
        s = Series(None, index=lrange(5))
        assert s.dtype == np.float64

        s = Series(None, index=lrange(5), dtype=object)
        assert s.dtype == np.object_

        # GH 7431
        # inference on the index
        s = Series(index=np.array([None]))
        expected = Series(index=Index([None]))
        assert_series_equal(s, expected)

    def test_constructor_pass_nan_nat(self):
        # GH 13467
        exp = Series([np.nan, np.nan], dtype=np.float64)
        assert exp.dtype == np.float64
        tm.assert_series_equal(Series([np.nan, np.nan]), exp)
        tm.assert_series_equal(Series(np.array([np.nan, np.nan])), exp)

        exp = Series([pd.NaT, pd.NaT])
        assert exp.dtype == 'datetime64[ns]'
        tm.assert_series_equal(Series([pd.NaT, pd.NaT]), exp)
        tm.assert_series_equal(Series(np.array([pd.NaT, pd.NaT])), exp)

        tm.assert_series_equal(Series([pd.NaT, np.nan]), exp)
        tm.assert_series_equal(Series(np.array([pd.NaT, np.nan])), exp)

        tm.assert_series_equal(Series([np.nan, pd.NaT]), exp)
        tm.assert_series_equal(Series(np.array([np.nan, pd.NaT])), exp)

    def test_constructor_cast(self):
        pytest.raises(ValueError, Series, ['a', 'b', 'c'], dtype=float)

    def test_constructor_dtype_nocast(self):
        # 1572
        s = Series([1, 2, 3])

        s2 = Series(s, dtype=np.int64)

        s2[1] = 5
        assert s[1] == 5

    def test_constructor_datelike_coercion(self):

        # GH 9477
        # incorrectly inferring on dateimelike looking when object dtype is
        # specified
        s = Series([Timestamp('20130101'), 'NOV'], dtype=object)
        assert s.iloc[0] == Timestamp('20130101')
        assert s.iloc[1] == 'NOV'
        assert s.dtype == object

        # the dtype was being reset on the slicing and re-inferred to datetime
        # even thought the blocks are mixed
        belly = '216 3T19'.split()
        wing1 = '2T15 4H19'.split()
        wing2 = '416 4T20'.split()
        mat = pd.to_datetime('2016-01-22 2019-09-07'.split())
        df = pd.DataFrame(
            {'wing1': wing1,
             'wing2': wing2,
             'mat': mat}, index=belly)

        result = df.loc['3T19']
        assert result.dtype == object
        result = df.loc['216']
        assert result.dtype == object

    def test_constructor_datetimes_with_nulls(self):
        # gh-15869
        for arr in [np.array([None, None, None, None,
                              datetime.now(), None]),
                    np.array([None, None, datetime.now(), None])]:
            result = Series(arr)
            assert result.dtype == 'M8[ns]'

    def test_constructor_dtype_datetime64(self):

        s = Series(iNaT, dtype='M8[ns]', index=lrange(5))
        assert isna(s).all()

        # in theory this should be all nulls, but since
        # we are not specifying a dtype is ambiguous
        s = Series(iNaT, index=lrange(5))
        assert not isna(s).all()

        s = Series(nan, dtype='M8[ns]', index=lrange(5))
        assert isna(s).all()

        s = Series([datetime(2001, 1, 2, 0, 0), iNaT], dtype='M8[ns]')
        assert isna(s[1])
        assert s.dtype == 'M8[ns]'

        s = Series([datetime(2001, 1, 2, 0, 0), nan], dtype='M8[ns]')
        assert isna(s[1])
        assert s.dtype == 'M8[ns]'

        # GH3416
        dates = [
            np.datetime64(datetime(2013, 1, 1)),
            np.datetime64(datetime(2013, 1, 2)),
            np.datetime64(datetime(2013, 1, 3)),
        ]

        s = Series(dates)
        assert s.dtype == 'M8[ns]'

        s.iloc[0] = np.nan
        assert s.dtype == 'M8[ns]'

        # GH3414 related
        pytest.raises(TypeError, lambda x: Series(
            Series(dates).astype('int') / 1000000, dtype='M8[ms]'))
        pytest.raises(TypeError,
                      lambda x: Series(dates, dtype='datetime64'))

        # invalid dates can be help as object
        result = Series([datetime(2, 1, 1)])
        assert result[0] == datetime(2, 1, 1, 0, 0)

        result = Series([datetime(3000, 1, 1)])
        assert result[0] == datetime(3000, 1, 1, 0, 0)

        # don't mix types
        result = Series([Timestamp('20130101'), 1], index=['a', 'b'])
        assert result['a'] == Timestamp('20130101')
        assert result['b'] == 1

        # GH6529
        # coerce datetime64 non-ns properly
        dates = date_range('01-Jan-2015', '01-Dec-2015', freq='M')
        values2 = dates.view(np.ndarray).astype('datetime64[ns]')
        expected = Series(values2, index=dates)

        for dtype in ['s', 'D', 'ms', 'us', 'ns']:
            values1 = dates.view(np.ndarray).astype('M8[{0}]'.format(dtype))
            result = Series(values1, dates)
            assert_series_equal(result, expected)

        # GH 13876
        # coerce to non-ns to object properly
        expected = Series(values2, index=dates, dtype=object)
        for dtype in ['s', 'D', 'ms', 'us', 'ns']:
            values1 = dates.view(np.ndarray).astype('M8[{0}]'.format(dtype))
            result = Series(values1, index=dates, dtype=object)
            assert_series_equal(result, expected)

        # leave datetime.date alone
        dates2 = np.array([d.date() for d in dates.to_pydatetime()],
                          dtype=object)
        series1 = Series(dates2, dates)
        tm.assert_numpy_array_equal(series1.values, dates2)
        assert series1.dtype == object

        # these will correctly infer a datetime
        s = Series([None, pd.NaT, '2013-08-05 15:30:00.000001'])
        assert s.dtype == 'datetime64[ns]'
        s = Series([np.nan, pd.NaT, '2013-08-05 15:30:00.000001'])
        assert s.dtype == 'datetime64[ns]'
        s = Series([pd.NaT, None, '2013-08-05 15:30:00.000001'])
        assert s.dtype == 'datetime64[ns]'
        s = Series([pd.NaT, np.nan, '2013-08-05 15:30:00.000001'])
        assert s.dtype == 'datetime64[ns]'

        # tz-aware (UTC and other tz's)
        # GH 8411
        dr = date_range('20130101', periods=3)
        assert Series(dr).iloc[0].tz is None
        dr = date_range('20130101', periods=3, tz='UTC')
        assert str(Series(dr).iloc[0].tz) == 'UTC'
        dr = date_range('20130101', periods=3, tz='US/Eastern')
        assert str(Series(dr).iloc[0].tz) == 'US/Eastern'

        # non-convertible
        s = Series([1479596223000, -1479590, pd.NaT])
        assert s.dtype == 'object'
        assert s[2] is pd.NaT
        assert 'NaT' in str(s)

        # if we passed a NaT it remains
        s = Series([datetime(2010, 1, 1), datetime(2, 1, 1), pd.NaT])
        assert s.dtype == 'object'
        assert s[2] is pd.NaT
        assert 'NaT' in str(s)

        # if we passed a nan it remains
        s = Series([datetime(2010, 1, 1), datetime(2, 1, 1), np.nan])
        assert s.dtype == 'object'
        assert s[2] is np.nan
        assert 'NaN' in str(s)

    def test_constructor_with_datetime_tz(self):

        # 8260
        # support datetime64 with tz

        dr = date_range('20130101', periods=3, tz='US/Eastern')
        s = Series(dr)
        assert s.dtype.name == 'datetime64[ns, US/Eastern]'
        assert s.dtype == 'datetime64[ns, US/Eastern]'
        assert is_datetime64tz_dtype(s.dtype)
        assert 'datetime64[ns, US/Eastern]' in str(s)

        # export
        result = s.values
        assert isinstance(result, np.ndarray)
        assert result.dtype == 'datetime64[ns]'

        exp = pd.DatetimeIndex(result)
        exp = exp.tz_localize('UTC').tz_convert(tz=s.dt.tz)
        tm.assert_index_equal(dr, exp)

        # indexing
        result = s.iloc[0]
        assert result == Timestamp('2013-01-01 00:00:00-0500',
                                   tz='US/Eastern', freq='D')
        result = s[0]
        assert result == Timestamp('2013-01-01 00:00:00-0500',
                                   tz='US/Eastern', freq='D')

        result = s[Series([True, True, False], index=s.index)]
        assert_series_equal(result, s[0:2])

        result = s.iloc[0:1]
        assert_series_equal(result, Series(dr[0:1]))

        # concat
        result = pd.concat([s.iloc[0:1], s.iloc[1:]])
        assert_series_equal(result, s)

        # short str
        assert 'datetime64[ns, US/Eastern]' in str(s)

        # formatting with NaT
        result = s.shift()
        assert 'datetime64[ns, US/Eastern]' in str(result)
        assert 'NaT' in str(result)

        # long str
        t = Series(date_range('20130101', periods=1000, tz='US/Eastern'))
        assert 'datetime64[ns, US/Eastern]' in str(t)

        result = pd.DatetimeIndex(s, freq='infer')
        tm.assert_index_equal(result, dr)

        # inference
        s = Series([pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),
                    pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Pacific')])
        assert s.dtype == 'datetime64[ns, US/Pacific]'
        assert lib.infer_dtype(s) == 'datetime64'

        s = Series([pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),
                    pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Eastern')])
        assert s.dtype == 'object'
        assert lib.infer_dtype(s) == 'datetime'

        # with all NaT
        s = Series(pd.NaT, index=[0, 1], dtype='datetime64[ns, US/Eastern]')
        expected = Series(pd.DatetimeIndex(['NaT', 'NaT'], tz='US/Eastern'))
        assert_series_equal(s, expected)

    @pytest.mark.parametrize("arr_dtype", [np.int64, np.float64])
    @pytest.mark.parametrize("dtype", ["M8", "m8"])
    @pytest.mark.parametrize("unit", ['ns', 'us', 'ms', 's', 'h', 'm', 'D'])
    def test_construction_to_datetimelike_unit(self, arr_dtype, dtype, unit):
        # tests all units
        # gh-19223
        dtype = "{}[{}]".format(dtype, unit)
        arr = np.array([1, 2, 3], dtype=arr_dtype)
        s = Series(arr)
        result = s.astype(dtype)
        expected = Series(arr.astype(dtype))

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('arg',
                             ['2013-01-01 00:00:00', pd.NaT, np.nan, None])
    def test_constructor_with_naive_string_and_datetimetz_dtype(self, arg):
        # GH 17415: With naive string
        result = Series([arg], dtype='datetime64[ns, CET]')
        expected = Series(pd.Timestamp(arg)).dt.tz_localize('CET')
        assert_series_equal(result, expected)

    def test_construction_interval(self):
        # construction from interval & array of intervals
        index = IntervalIndex.from_breaks(np.arange(3), closed='right')
        result = Series(index)
        repr(result)
        str(result)
        tm.assert_index_equal(Index(result.values), index)

        result = Series(index.values)
        tm.assert_index_equal(Index(result.values), index)

    def test_construction_consistency(self):

        # make sure that we are not re-localizing upon construction
        # GH 14928
        s = Series(pd.date_range('20130101', periods=3, tz='US/Eastern'))

        result = Series(s, dtype=s.dtype)
        tm.assert_series_equal(result, s)

        result = Series(s.dt.tz_convert('UTC'), dtype=s.dtype)
        tm.assert_series_equal(result, s)

        result = Series(s.values, dtype=s.dtype)
        tm.assert_series_equal(result, s)

    def test_constructor_periodindex(self):
        # GH7932
        # converting a PeriodIndex when put in a Series

        pi = period_range('20130101', periods=5, freq='D')
        s = Series(pi)
        expected = Series(pi.astype(object))
        assert_series_equal(s, expected)

        assert s.dtype == 'object'

    def test_constructor_dict(self):
        d = {'a': 0., 'b': 1., 'c': 2.}
        result = Series(d, index=['b', 'c', 'd', 'a'])
        expected = Series([1, 2, nan, 0], index=['b', 'c', 'd', 'a'])
        assert_series_equal(result, expected)

        pidx = tm.makePeriodIndex(100)
        d = {pidx[0]: 0, pidx[1]: 1}
        result = Series(d, index=pidx)
        expected = Series(np.nan, pidx)
        expected.iloc[0] = 0
        expected.iloc[1] = 1
        assert_series_equal(result, expected)

    def test_constructor_dict_order(self):
        # GH19018
        # initialization ordering: by insertion order if python>= 3.6, else
        # order by value
        d = {'b': 1, 'a': 0, 'c': 2}
        result = Series(d)
        if PY36:
            expected = Series([1, 0, 2], index=list('bac'))
        else:
            expected = Series([0, 1, 2], index=list('abc'))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("value", [2, np.nan, None, float('nan')])
    def test_constructor_dict_nan_key(self, value):
        # GH 18480
        d = {1: 'a', value: 'b', float('nan'): 'c', 4: 'd'}
        result = Series(d).sort_values()
        expected = Series(['a', 'b', 'c', 'd'], index=[1, value, np.nan, 4])
        assert_series_equal(result, expected)

        # MultiIndex:
        d = {(1, 1): 'a', (2, np.nan): 'b', (3, value): 'c'}
        result = Series(d).sort_values()
        expected = Series(['a', 'b', 'c'],
                          index=Index([(1, 1), (2, np.nan), (3, value)]))
        assert_series_equal(result, expected)

    def test_constructor_dict_datetime64_index(self):
        # GH 9456

        dates_as_str = ['1984-02-19', '1988-11-06', '1989-12-03', '1990-03-15']
        values = [42544017.198965244, 1234565, 40512335.181958228, -1]

        def create_data(constructor):
            return dict(zip((constructor(x) for x in dates_as_str), values))

        data_datetime64 = create_data(np.datetime64)
        data_datetime = create_data(lambda x: datetime.strptime(x, '%Y-%m-%d'))
        data_Timestamp = create_data(Timestamp)

        expected = Series(values, (Timestamp(x) for x in dates_as_str))

        result_datetime64 = Series(data_datetime64)
        result_datetime = Series(data_datetime)
        result_Timestamp = Series(data_Timestamp)

        assert_series_equal(result_datetime64, expected)
        assert_series_equal(result_datetime, expected)
        assert_series_equal(result_Timestamp, expected)

    def test_constructor_list_of_tuples(self):
        data = [(1, 1), (2, 2), (2, 3)]
        s = Series(data)
        assert list(s) == data

    def test_constructor_tuple_of_tuples(self):
        data = ((1, 1), (2, 2), (2, 3))
        s = Series(data)
        assert tuple(s) == data

    def test_constructor_dict_of_tuples(self):
        data = {(1, 2): 3,
                (None, 5): 6}
        result = Series(data).sort_values()
        expected = Series([3, 6],
                          index=MultiIndex.from_tuples([(1, 2), (None, 5)]))
        tm.assert_series_equal(result, expected)

    def test_constructor_set(self):
        values = set([1, 2, 3, 4, 5])
        pytest.raises(TypeError, Series, values)
        values = frozenset(values)
        pytest.raises(TypeError, Series, values)

    def test_fromDict(self):
        data = {'a': 0, 'b': 1, 'c': 2, 'd': 3}

        series = Series(data)
        assert tm.is_sorted(series.index)

        data = {'a': 0, 'b': '1', 'c': '2', 'd': datetime.now()}
        series = Series(data)
        assert series.dtype == np.object_

        data = {'a': 0, 'b': '1', 'c': '2', 'd': '3'}
        series = Series(data)
        assert series.dtype == np.object_

        data = {'a': '0', 'b': '1'}
        series = Series(data, dtype=float)
        assert series.dtype == np.float64

    def test_fromValue(self):

        nans = Series(np.NaN, index=self.ts.index)
        assert nans.dtype == np.float_
        assert len(nans) == len(self.ts)

        strings = Series('foo', index=self.ts.index)
        assert strings.dtype == np.object_
        assert len(strings) == len(self.ts)

        d = datetime.now()
        dates = Series(d, index=self.ts.index)
        assert dates.dtype == 'M8[ns]'
        assert len(dates) == len(self.ts)

        # GH12336
        # Test construction of categorical series from value
        categorical = Series(0, index=self.ts.index, dtype="category")
        expected = Series(0, index=self.ts.index).astype("category")
        assert categorical.dtype == 'category'
        assert len(categorical) == len(self.ts)
        tm.assert_series_equal(categorical, expected)

    def test_constructor_dtype_timedelta64(self):

        # basic
        td = Series([timedelta(days=i) for i in range(3)])
        assert td.dtype == 'timedelta64[ns]'

        td = Series([timedelta(days=1)])
        assert td.dtype == 'timedelta64[ns]'

        td = Series([timedelta(days=1), timedelta(days=2), np.timedelta64(
            1, 's')])

        assert td.dtype == 'timedelta64[ns]'

        # mixed with NaT
        td = Series([timedelta(days=1), NaT], dtype='m8[ns]')
        assert td.dtype == 'timedelta64[ns]'

        td = Series([timedelta(days=1), np.nan], dtype='m8[ns]')
        assert td.dtype == 'timedelta64[ns]'

        td = Series([np.timedelta64(300000000), pd.NaT], dtype='m8[ns]')
        assert td.dtype == 'timedelta64[ns]'

        # improved inference
        # GH5689
        td = Series([np.timedelta64(300000000), NaT])
        assert td.dtype == 'timedelta64[ns]'

        # because iNaT is int, not coerced to timedelta
        td = Series([np.timedelta64(300000000), iNaT])
        assert td.dtype == 'object'

        td = Series([np.timedelta64(300000000), np.nan])
        assert td.dtype == 'timedelta64[ns]'

        td = Series([pd.NaT, np.timedelta64(300000000)])
        assert td.dtype == 'timedelta64[ns]'

        td = Series([np.timedelta64(1, 's')])
        assert td.dtype == 'timedelta64[ns]'

        # these are frequency conversion astypes
        # for t in ['s', 'D', 'us', 'ms']:
        #    pytest.raises(TypeError, td.astype, 'm8[%s]' % t)

        # valid astype
        td.astype('int64')

        # invalid casting
        pytest.raises(TypeError, td.astype, 'int32')

        # this is an invalid casting
        def f():
            Series([timedelta(days=1), 'foo'], dtype='m8[ns]')

        pytest.raises(Exception, f)

        # leave as object here
        td = Series([timedelta(days=i) for i in range(3)] + ['foo'])
        assert td.dtype == 'object'

        # these will correctly infer a timedelta
        s = Series([None, pd.NaT, '1 Day'])
        assert s.dtype == 'timedelta64[ns]'
        s = Series([np.nan, pd.NaT, '1 Day'])
        assert s.dtype == 'timedelta64[ns]'
        s = Series([pd.NaT, None, '1 Day'])
        assert s.dtype == 'timedelta64[ns]'
        s = Series([pd.NaT, np.nan, '1 Day'])
        assert s.dtype == 'timedelta64[ns]'

    # GH 16406
    def test_constructor_mixed_tz(self):
        s = Series([Timestamp('20130101'),
                    Timestamp('20130101', tz='US/Eastern')])
        expected = Series([Timestamp('20130101'),
                           Timestamp('20130101', tz='US/Eastern')],
                          dtype='object')
        assert_series_equal(s, expected)

    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, iNaT], dtype='M8[ns]')

        val = series[3]
        assert isna(val)

        series[2] = val
        assert isna(series[2])

    def test_NaT_cast(self):
        # GH10747
        result = Series([np.nan]).astype('M8[ns]')
        expected = Series([NaT])
        assert_series_equal(result, expected)

    def test_constructor_name_hashable(self):
        for n in [777, 777., 'name', datetime(2001, 11, 11), (1, ), u"\u05D0"]:
            for data in [[1, 2, 3], np.ones(3), {'a': 0, 'b': 1}]:
                s = Series(data, name=n)
                assert s.name == n

    def test_constructor_name_unhashable(self):
        for n in [['name_list'], np.ones(2), {1: 2}]:
            for data in [['name_list'], np.ones(2), {1: 2}]:
                pytest.raises(TypeError, Series, data, name=n)

    def test_auto_conversion(self):
        series = Series(list(date_range('1/1/2000', periods=10)))
        assert series.dtype == 'M8[ns]'

    def test_convert_non_ns(self):
        # convert from a numpy array of non-ns timedelta64
        arr = np.array([1, 2, 3], dtype='timedelta64[s]')
        s = Series(arr)
        expected = Series(pd.timedelta_range('00:00:01', periods=3, freq='s'))
        assert_series_equal(s, expected)

        # convert from a numpy array of non-ns datetime64
        # note that creating a numpy datetime64 is in LOCAL time!!!!
        # seems to work for M8[D], but not for M8[s]

        s = Series(np.array(['2013-01-01', '2013-01-02',
                             '2013-01-03'], dtype='datetime64[D]'))
        assert_series_equal(s, Series(date_range('20130101', periods=3,
                                                 freq='D')))

        # s = Series(np.array(['2013-01-01 00:00:01','2013-01-01
        # 00:00:02','2013-01-01 00:00:03'],dtype='datetime64[s]'))

        # assert_series_equal(s,date_range('20130101
        # 00:00:01',period=3,freq='s'))

    @pytest.mark.parametrize(
        "index",
        [
            date_range('1/1/2000', periods=10),
            timedelta_range('1 day', periods=10),
            period_range('2000-Q1', periods=10, freq='Q')],
        ids=lambda x: type(x).__name__)
    def test_constructor_cant_cast_datetimelike(self, index):

        # floats are not ok
        msg = "Cannot cast {} to ".format(type(index).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            Series(index, dtype=float)

        # ints are ok
        # we test with np.int64 to get similar results on
        # windows / 32-bit platforms
        result = Series(index, dtype=np.int64)
        expected = Series(index.astype(np.int64))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "index",
        [
            date_range('1/1/2000', periods=10),
            timedelta_range('1 day', periods=10),
            period_range('2000-Q1', periods=10, freq='Q')],
        ids=lambda x: type(x).__name__)
    def test_constructor_cast_object(self, index):
        s = Series(index, dtype=object)
        exp = Series(index).astype(object)
        tm.assert_series_equal(s, exp)

        s = Series(pd.Index(index, dtype=object), dtype=object)
        exp = Series(index).astype(object)
        tm.assert_series_equal(s, exp)

        s = Series(index.astype(object), dtype=object)
        exp = Series(index).astype(object)
        tm.assert_series_equal(s, exp)

    def test_constructor_generic_timestamp_deprecated(self):
        # see gh-15524

        with tm.assert_produces_warning(FutureWarning):
            dtype = np.timedelta64
            s = Series([], dtype=dtype)

            assert s.empty
            assert s.dtype == 'm8[ns]'

        with tm.assert_produces_warning(FutureWarning):
            dtype = np.datetime64
            s = Series([], dtype=dtype)

            assert s.empty
            assert s.dtype == 'M8[ns]'

        # These timestamps have the wrong frequencies,
        # so an Exception should be raised now.
        msg = "cannot convert timedeltalike"
        with tm.assert_raises_regex(TypeError, msg):
            Series([], dtype='m8[ps]')

        msg = "cannot convert datetimelike"
        with tm.assert_raises_regex(TypeError, msg):
            Series([], dtype='M8[ps]')

    @pytest.mark.parametrize('dtype', [None, 'uint8', 'category'])
    def test_constructor_range_dtype(self, dtype):
        # GH 16804
        expected = Series([0, 1, 2, 3, 4], dtype=dtype or 'int64')
        result = Series(range(5), dtype=dtype)
        tm.assert_series_equal(result, expected)
