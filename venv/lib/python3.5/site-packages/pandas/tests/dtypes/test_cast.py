# -*- coding: utf-8 -*-

"""
These test the private routines in types/cast.py

"""

import pytest
from datetime import datetime, timedelta, date
import numpy as np

import pandas as pd
from pandas import (Timedelta, Timestamp, DatetimeIndex,
                    DataFrame, NaT, Period, Series)

from pandas.core.dtypes.cast import (
    maybe_downcast_to_dtype,
    maybe_convert_objects,
    cast_scalar_to_array,
    infer_dtype_from_scalar,
    infer_dtype_from_array,
    maybe_convert_string_to_object,
    maybe_convert_scalar,
    find_common_type,
    construct_1d_object_array_from_listlike,
    construct_1d_arraylike_from_scalar)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    DatetimeTZDtype,
    PeriodDtype)
from pandas.core.dtypes.common import (
    is_dtype_equal)
from pandas.util import testing as tm


class TestMaybeDowncast(object):

    def test_downcast_conv(self):
        # test downcasting

        arr = np.array([8.5, 8.6, 8.7, 8.8, 8.9999999999995])
        result = maybe_downcast_to_dtype(arr, 'infer')
        tm.assert_numpy_array_equal(result, arr)

        arr = np.array([8., 8., 8., 8., 8.9999999999995])
        result = maybe_downcast_to_dtype(arr, 'infer')
        expected = np.array([8, 8, 8, 8, 9], dtype=np.int64)
        tm.assert_numpy_array_equal(result, expected)

        arr = np.array([8., 8., 8., 8., 9.0000000000005])
        result = maybe_downcast_to_dtype(arr, 'infer')
        expected = np.array([8, 8, 8, 8, 9], dtype=np.int64)
        tm.assert_numpy_array_equal(result, expected)

        # GH16875 coercing of bools
        ser = Series([True, True, False])
        result = maybe_downcast_to_dtype(ser, np.dtype(np.float64))
        expected = ser
        tm.assert_series_equal(result, expected)

        # conversions

        expected = np.array([1, 2])
        for dtype in [np.float64, object, np.int64]:
            arr = np.array([1.0, 2.0], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'infer')
            tm.assert_almost_equal(result, expected, check_dtype=False)

        for dtype in [np.float64, object]:
            expected = np.array([1.0, 2.0, np.nan], dtype=dtype)
            arr = np.array([1.0, 2.0, np.nan], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'infer')
            tm.assert_almost_equal(result, expected)

        # empties
        for dtype in [np.int32, np.float64, np.float32, np.bool_,
                      np.int64, object]:
            arr = np.array([], dtype=dtype)
            result = maybe_downcast_to_dtype(arr, 'int64')
            tm.assert_almost_equal(result, np.array([], dtype=np.int64))
            assert result.dtype == np.int64

    def test_datetimelikes_nan(self):
        arr = np.array([1, 2, np.nan])
        exp = np.array([1, 2, np.datetime64('NaT')], dtype='datetime64[ns]')
        res = maybe_downcast_to_dtype(arr, 'datetime64[ns]')
        tm.assert_numpy_array_equal(res, exp)

        exp = np.array([1, 2, np.timedelta64('NaT')], dtype='timedelta64[ns]')
        res = maybe_downcast_to_dtype(arr, 'timedelta64[ns]')
        tm.assert_numpy_array_equal(res, exp)

    def test_datetime_with_timezone(self):
        # GH 15426
        ts = Timestamp("2016-01-01 12:00:00", tz='US/Pacific')
        exp = DatetimeIndex([ts, ts])
        res = maybe_downcast_to_dtype(exp, exp.dtype)
        tm.assert_index_equal(res, exp)

        res = maybe_downcast_to_dtype(exp.asi8, exp.dtype)
        tm.assert_index_equal(res, exp)


class TestInferDtype(object):

    def testinfer_dtype_from_scalar(self):
        # Test that infer_dtype_from_scalar is returning correct dtype for int
        # and float.

        for dtypec in [np.uint8, np.int8, np.uint16, np.int16, np.uint32,
                       np.int32, np.uint64, np.int64]:
            data = dtypec(12)
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == type(data)

        data = 12
        dtype, val = infer_dtype_from_scalar(data)
        assert dtype == np.int64

        for dtypec in [np.float16, np.float32, np.float64]:
            data = dtypec(12)
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == dtypec

        data = np.float(12)
        dtype, val = infer_dtype_from_scalar(data)
        assert dtype == np.float64

        for data in [True, False]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.bool_

        for data in [np.complex64(1), np.complex128(1)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.complex_

        for data in [np.datetime64(1, 'ns'), Timestamp(1),
                     datetime(2000, 1, 1, 0, 0)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == 'M8[ns]'

        for data in [np.timedelta64(1, 'ns'), Timedelta(1),
                     timedelta(1)]:
            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == 'm8[ns]'

        for freq in ['M', 'D']:
            p = Period('2011-01-01', freq=freq)
            dtype, val = infer_dtype_from_scalar(p, pandas_dtype=True)
            assert dtype == 'period[{0}]'.format(freq)
            assert val == p.ordinal

            dtype, val = infer_dtype_from_scalar(p)
            dtype == np.object_
            assert val == p

        # misc
        for data in [date(2000, 1, 1),
                     Timestamp(1, tz='US/Eastern'), 'foo']:

            dtype, val = infer_dtype_from_scalar(data)
            assert dtype == np.object_

    @pytest.mark.parametrize('tz', ['UTC', 'US/Eastern', 'Asia/Tokyo'])
    def testinfer_from_scalar_tz(self, tz):
        dt = Timestamp(1, tz=tz)
        dtype, val = infer_dtype_from_scalar(dt, pandas_dtype=True)
        assert dtype == 'datetime64[ns, {0}]'.format(tz)
        assert val == dt.value

        dtype, val = infer_dtype_from_scalar(dt)
        assert dtype == np.object_
        assert val == dt

    def testinfer_dtype_from_scalar_errors(self):
        with pytest.raises(ValueError):
            infer_dtype_from_scalar(np.array([1]))

    @pytest.mark.parametrize(
        "arr, expected, pandas_dtype",
        [('foo', np.object_, False),
         (b'foo', np.object_, False),
         (1, np.int_, False),
         (1.5, np.float_, False),
         ([1], np.int_, False),
         (np.array([1], dtype=np.int64), np.int64, False),
         ([np.nan, 1, ''], np.object_, False),
         (np.array([[1.0, 2.0]]), np.float_, False),
         (pd.Categorical(list('aabc')), np.object_, False),
         (pd.Categorical([1, 2, 3]), np.int64, False),
         (pd.Categorical(list('aabc')), 'category', True),
         (pd.Categorical([1, 2, 3]), 'category', True),
         (Timestamp('20160101'), np.object_, False),
         (np.datetime64('2016-01-01'), np.dtype('=M8[D]'), False),
         (pd.date_range('20160101', periods=3),
          np.dtype('=M8[ns]'), False),
         (pd.date_range('20160101', periods=3, tz='US/Eastern'),
          'datetime64[ns, US/Eastern]', True),
         (pd.Series([1., 2, 3]), np.float64, False),
         (pd.Series(list('abc')), np.object_, False),
         (pd.Series(pd.date_range('20160101', periods=3, tz='US/Eastern')),
          'datetime64[ns, US/Eastern]', True)])
    def test_infer_dtype_from_array(self, arr, expected, pandas_dtype):

        dtype, _ = infer_dtype_from_array(arr, pandas_dtype=pandas_dtype)
        assert is_dtype_equal(dtype, expected)

    def test_cast_scalar_to_array(self):
        arr = cast_scalar_to_array((3, 2), 1, dtype=np.int64)
        exp = np.ones((3, 2), dtype=np.int64)
        tm.assert_numpy_array_equal(arr, exp)

        arr = cast_scalar_to_array((3, 2), 1.1)
        exp = np.empty((3, 2), dtype=np.float64)
        exp.fill(1.1)
        tm.assert_numpy_array_equal(arr, exp)

        arr = cast_scalar_to_array((2, 3), Timestamp('2011-01-01'))
        exp = np.empty((2, 3), dtype='datetime64[ns]')
        exp.fill(np.datetime64('2011-01-01'))
        tm.assert_numpy_array_equal(arr, exp)

        # pandas dtype is stored as object dtype
        obj = Timestamp('2011-01-01', tz='US/Eastern')
        arr = cast_scalar_to_array((2, 3), obj)
        exp = np.empty((2, 3), dtype=np.object)
        exp.fill(obj)
        tm.assert_numpy_array_equal(arr, exp)

        obj = Period('2011-01-01', freq='D')
        arr = cast_scalar_to_array((2, 3), obj)
        exp = np.empty((2, 3), dtype=np.object)
        exp.fill(obj)
        tm.assert_numpy_array_equal(arr, exp)


class TestMaybe(object):

    def test_maybe_convert_string_to_array(self):
        result = maybe_convert_string_to_object('x')
        tm.assert_numpy_array_equal(result, np.array(['x'], dtype=object))
        assert result.dtype == object

        result = maybe_convert_string_to_object(1)
        assert result == 1

        arr = np.array(['x', 'y'], dtype=str)
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        assert result.dtype == object

        # unicode
        arr = np.array(['x', 'y']).astype('U')
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        assert result.dtype == object

        # object
        arr = np.array(['x', 2], dtype=object)
        result = maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 2], dtype=object))
        assert result.dtype == object

    def test_maybe_convert_scalar(self):

        # pass thru
        result = maybe_convert_scalar('x')
        assert result == 'x'
        result = maybe_convert_scalar(np.array([1]))
        assert result == np.array([1])

        # leave scalar dtype
        result = maybe_convert_scalar(np.int64(1))
        assert result == np.int64(1)
        result = maybe_convert_scalar(np.int32(1))
        assert result == np.int32(1)
        result = maybe_convert_scalar(np.float32(1))
        assert result == np.float32(1)
        result = maybe_convert_scalar(np.int64(1))
        assert result == np.float64(1)

        # coerce
        result = maybe_convert_scalar(1)
        assert result == np.int64(1)
        result = maybe_convert_scalar(1.0)
        assert result == np.float64(1)
        result = maybe_convert_scalar(Timestamp('20130101'))
        assert result == Timestamp('20130101').value
        result = maybe_convert_scalar(datetime(2013, 1, 1))
        assert result == Timestamp('20130101').value
        result = maybe_convert_scalar(Timedelta('1 day 1 min'))
        assert result == Timedelta('1 day 1 min').value

    def test_maybe_infer_to_datetimelike(self):
        # GH16362
        # pandas=0.20.1 raises IndexError: tuple index out of range
        result = DataFrame(np.array([[NaT, 'a', 'b', 0],
                                     [NaT, 'b', 'c', 1]]))
        assert result.size == 8
        # this construction was fine
        result = DataFrame(np.array([[NaT, 'a', 0],
                                     [NaT, 'b', 1]]))
        assert result.size == 6

        # GH19671
        result = Series(['M1701', Timestamp('20130101')])
        assert result.dtype.kind == 'O'


class TestConvert(object):

    def test_maybe_convert_objects_copy(self):
        values = np.array([1, 2])

        out = maybe_convert_objects(values, copy=False)
        assert values is out

        out = maybe_convert_objects(values, copy=True)
        assert values is not out

        values = np.array(['apply', 'banana'])
        out = maybe_convert_objects(values, copy=False)
        assert values is out

        out = maybe_convert_objects(values, copy=True)
        assert values is not out


class TestCommonTypes(object):

    def test_numpy_dtypes(self):
        # (source_types, destination_type)
        testcases = (
            # identity
            ((np.int64,), np.int64),
            ((np.uint64,), np.uint64),
            ((np.float32,), np.float32),
            ((np.object,), np.object),

            # into ints
            ((np.int16, np.int64), np.int64),
            ((np.int32, np.uint32), np.int64),
            ((np.uint16, np.uint64), np.uint64),

            # into floats
            ((np.float16, np.float32), np.float32),
            ((np.float16, np.int16), np.float32),
            ((np.float32, np.int16), np.float32),
            ((np.uint64, np.int64), np.float64),
            ((np.int16, np.float64), np.float64),
            ((np.float16, np.int64), np.float64),

            # into others
            ((np.complex128, np.int32), np.complex128),
            ((np.object, np.float32), np.object),
            ((np.object, np.int16), np.object),

            # bool with int
            ((np.dtype('bool'), np.int64), np.object),
            ((np.dtype('bool'), np.int32), np.object),
            ((np.dtype('bool'), np.int16), np.object),
            ((np.dtype('bool'), np.int8), np.object),
            ((np.dtype('bool'), np.uint64), np.object),
            ((np.dtype('bool'), np.uint32), np.object),
            ((np.dtype('bool'), np.uint16), np.object),
            ((np.dtype('bool'), np.uint8), np.object),

            # bool with float
            ((np.dtype('bool'), np.float64), np.object),
            ((np.dtype('bool'), np.float32), np.object),

            ((np.dtype('datetime64[ns]'), np.dtype('datetime64[ns]')),
             np.dtype('datetime64[ns]')),
            ((np.dtype('timedelta64[ns]'), np.dtype('timedelta64[ns]')),
             np.dtype('timedelta64[ns]')),

            ((np.dtype('datetime64[ns]'), np.dtype('datetime64[ms]')),
             np.dtype('datetime64[ns]')),
            ((np.dtype('timedelta64[ms]'), np.dtype('timedelta64[ns]')),
             np.dtype('timedelta64[ns]')),

            ((np.dtype('datetime64[ns]'), np.dtype('timedelta64[ns]')),
             np.object),
            ((np.dtype('datetime64[ns]'), np.int64), np.object)
        )
        for src, common in testcases:
            assert find_common_type(src) == common

        with pytest.raises(ValueError):
            # empty
            find_common_type([])

    def test_categorical_dtype(self):
        dtype = CategoricalDtype()
        assert find_common_type([dtype]) == 'category'
        assert find_common_type([dtype, dtype]) == 'category'
        assert find_common_type([np.object, dtype]) == np.object

    def test_datetimetz_dtype(self):
        dtype = DatetimeTZDtype(unit='ns', tz='US/Eastern')
        assert find_common_type([dtype, dtype]) == 'datetime64[ns, US/Eastern]'

        for dtype2 in [DatetimeTZDtype(unit='ns', tz='Asia/Tokyo'),
                       np.dtype('datetime64[ns]'), np.object, np.int64]:
            assert find_common_type([dtype, dtype2]) == np.object
            assert find_common_type([dtype2, dtype]) == np.object

    def test_period_dtype(self):
        dtype = PeriodDtype(freq='D')
        assert find_common_type([dtype, dtype]) == 'period[D]'

        for dtype2 in [DatetimeTZDtype(unit='ns', tz='Asia/Tokyo'),
                       PeriodDtype(freq='2D'), PeriodDtype(freq='H'),
                       np.dtype('datetime64[ns]'), np.object, np.int64]:
            assert find_common_type([dtype, dtype2]) == np.object
            assert find_common_type([dtype2, dtype]) == np.object

    @pytest.mark.parametrize('datum1', [1, 2., "3", (4, 5), [6, 7], None])
    @pytest.mark.parametrize('datum2', [8, 9., "10", (11, 12), [13, 14], None])
    def test_cast_1d_array(self, datum1, datum2):
        data = [datum1, datum2]
        result = construct_1d_object_array_from_listlike(data)

        # Direct comparison fails: https://github.com/numpy/numpy/issues/10218
        assert result.dtype == 'object'
        assert list(result) == data

    @pytest.mark.parametrize('val', [1, 2., None])
    def test_cast_1d_array_invalid_scalar(self, val):
        pytest.raises(TypeError, construct_1d_object_array_from_listlike, val)

    def test_cast_1d_arraylike_from_scalar_categorical(self):
        # GH 19565 - Categorical result from scalar did not maintain categories
        # and ordering of the passed dtype
        cats = ['a', 'b', 'c']
        cat_type = CategoricalDtype(categories=cats, ordered=False)
        expected = pd.Categorical(['a', 'a'], categories=cats)
        result = construct_1d_arraylike_from_scalar('a', len(expected),
                                                    cat_type)
        tm.assert_categorical_equal(result, expected,
                                    check_category_order=True,
                                    check_dtype=True)
