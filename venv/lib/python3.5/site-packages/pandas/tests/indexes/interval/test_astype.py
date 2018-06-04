from __future__ import division

import pytest
import numpy as np
from pandas import (
    Index,
    IntervalIndex,
    interval_range,
    CategoricalIndex,
    Timestamp,
    Timedelta,
    NaT)
from pandas.core.dtypes.dtypes import CategoricalDtype, IntervalDtype
import pandas.util.testing as tm


class Base(object):
    """Tests common to IntervalIndex with any subtype"""

    def test_astype_idempotent(self, index):
        result = index.astype('interval')
        tm.assert_index_equal(result, index)

        result = index.astype(index.dtype)
        tm.assert_index_equal(result, index)

    def test_astype_object(self, index):
        result = index.astype(object)
        expected = Index(index.values, dtype='object')
        tm.assert_index_equal(result, expected)
        assert not result.equals(index)

    def test_astype_category(self, index):
        result = index.astype('category')
        expected = CategoricalIndex(index.values)
        tm.assert_index_equal(result, expected)

        result = index.astype(CategoricalDtype())
        tm.assert_index_equal(result, expected)

        # non-default params
        categories = index.dropna().unique().values[:-1]
        dtype = CategoricalDtype(categories=categories, ordered=True)
        result = index.astype(dtype)
        expected = CategoricalIndex(
            index.values, categories=categories, ordered=True)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('dtype', [
        'int64', 'uint64', 'float64', 'complex128', 'period[M]',
        'timedelta64', 'timedelta64[ns]', 'datetime64', 'datetime64[ns]',
        'datetime64[ns, US/Eastern]'])
    def test_astype_cannot_cast(self, index, dtype):
        msg = 'Cannot cast IntervalIndex to dtype'
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)

    def test_astype_invalid_dtype(self, index):
        msg = 'data type "fake_dtype" not understood'
        with tm.assert_raises_regex(TypeError, msg):
            index.astype('fake_dtype')


class TestIntSubtype(Base):
    """Tests specific to IntervalIndex with integer-like subtype"""

    indexes = [
        IntervalIndex.from_breaks(np.arange(-10, 11, dtype='int64')),
        IntervalIndex.from_breaks(
            np.arange(100, dtype='uint64'), closed='left'),
    ]

    @pytest.fixture(params=indexes)
    def index(self, request):
        return request.param

    @pytest.mark.parametrize('subtype', [
        'float64', 'datetime64[ns]', 'timedelta64[ns]'])
    def test_subtype_conversion(self, index, subtype):
        dtype = IntervalDtype(subtype)
        result = index.astype(dtype)
        expected = IntervalIndex.from_arrays(index.left.astype(subtype),
                                             index.right.astype(subtype),
                                             closed=index.closed)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('subtype_start, subtype_end', [
        ('int64', 'uint64'), ('uint64', 'int64')])
    def test_subtype_integer(self, subtype_start, subtype_end):
        index = IntervalIndex.from_breaks(np.arange(100, dtype=subtype_start))
        dtype = IntervalDtype(subtype_end)
        result = index.astype(dtype)
        expected = IntervalIndex.from_arrays(index.left.astype(subtype_end),
                                             index.right.astype(subtype_end),
                                             closed=index.closed)
        tm.assert_index_equal(result, expected)

    @pytest.mark.xfail(reason='GH 15832')
    def test_subtype_integer_errors(self):
        # int64 -> uint64 fails with negative values
        index = interval_range(-10, 10)
        dtype = IntervalDtype('uint64')
        with pytest.raises(ValueError):
            index.astype(dtype)


class TestFloatSubtype(Base):
    """Tests specific to IntervalIndex with float subtype"""

    indexes = [
        interval_range(-10.0, 10.0, closed='neither'),
        IntervalIndex.from_arrays([-1.5, np.nan, 0., 0., 1.5],
                                  [-0.5, np.nan, 1., 1., 3.],
                                  closed='both'),
    ]

    @pytest.fixture(params=indexes)
    def index(self, request):
        return request.param

    @pytest.mark.parametrize('subtype', ['int64', 'uint64'])
    def test_subtype_integer(self, subtype):
        index = interval_range(0.0, 10.0)
        dtype = IntervalDtype(subtype)
        result = index.astype(dtype)
        expected = IntervalIndex.from_arrays(index.left.astype(subtype),
                                             index.right.astype(subtype),
                                             closed=index.closed)
        tm.assert_index_equal(result, expected)

        # raises with NA
        msg = 'Cannot convert NA to integer'
        with tm.assert_raises_regex(ValueError, msg):
            index.insert(0, np.nan).astype(dtype)

    @pytest.mark.xfail(reason='GH 15832')
    def test_subtype_integer_errors(self):
        # float64 -> uint64 fails with negative values
        index = interval_range(-10.0, 10.0)
        dtype = IntervalDtype('uint64')
        with pytest.raises(ValueError):
            index.astype(dtype)

        # float64 -> integer-like fails with non-integer valued floats
        index = interval_range(0.0, 10.0, freq=0.25)
        dtype = IntervalDtype('int64')
        with pytest.raises(ValueError):
            index.astype(dtype)

        dtype = IntervalDtype('uint64')
        with pytest.raises(ValueError):
            index.astype(dtype)

    @pytest.mark.parametrize('subtype', ['datetime64[ns]', 'timedelta64[ns]'])
    def test_subtype_datetimelike(self, index, subtype):
        dtype = IntervalDtype(subtype)
        msg = 'Cannot convert .* to .*; subtypes are incompatible'
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)


class TestDatetimelikeSubtype(Base):
    """Tests specific to IntervalIndex with datetime-like subtype"""

    indexes = [
        interval_range(Timestamp('2018-01-01'), periods=10, closed='neither'),
        interval_range(Timestamp('2018-01-01'), periods=10).insert(2, NaT),
        interval_range(Timestamp('2018-01-01', tz='US/Eastern'), periods=10),
        interval_range(Timedelta('0 days'), periods=10, closed='both'),
        interval_range(Timedelta('0 days'), periods=10).insert(2, NaT),
    ]

    @pytest.fixture(params=indexes)
    def index(self, request):
        return request.param

    @pytest.mark.parametrize('subtype', ['int64', 'uint64'])
    def test_subtype_integer(self, index, subtype):
        dtype = IntervalDtype(subtype)
        result = index.astype(dtype)
        expected = IntervalIndex.from_arrays(index.left.astype(subtype),
                                             index.right.astype(subtype),
                                             closed=index.closed)
        tm.assert_index_equal(result, expected)

    def test_subtype_float(self, index):
        dtype = IntervalDtype('float64')
        msg = 'Cannot convert .* to .*; subtypes are incompatible'
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)

    def test_subtype_datetimelike(self):
        # datetime -> timedelta raises
        dtype = IntervalDtype('timedelta64[ns]')
        msg = 'Cannot convert .* to .*; subtypes are incompatible'

        index = interval_range(Timestamp('2018-01-01'), periods=10)
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)

        index = interval_range(Timestamp('2018-01-01', tz='CET'), periods=10)
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)

        # timedelta -> datetime raises
        dtype = IntervalDtype('datetime64[ns]')
        index = interval_range(Timedelta('0 days'), periods=10)
        with tm.assert_raises_regex(TypeError, msg):
            index.astype(dtype)
