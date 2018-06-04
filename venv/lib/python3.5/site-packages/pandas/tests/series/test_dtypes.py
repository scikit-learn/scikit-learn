# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

from datetime import datetime, timedelta

import sys
import string
import warnings

from numpy import nan
import pandas as pd
import numpy as np

from pandas import (
    Series, Timestamp, Timedelta, DataFrame, date_range,
    Categorical, Index
)
from pandas.api.types import CategoricalDtype
import pandas._libs.tslib as tslib

from pandas.compat import lrange, range, u
from pandas import compat
import pandas.util.testing as tm

from .common import TestData


class TestSeriesDtypes(TestData):

    def test_dt64_series_astype_object(self):
        dt64ser = Series(date_range('20130101', periods=3))
        result = dt64ser.astype(object)
        assert isinstance(result.iloc[0], datetime)
        assert result.dtype == np.object_

    def test_td64_series_astype_object(self):
        tdser = Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')
        result = tdser.astype(object)
        assert isinstance(result.iloc[0], timedelta)
        assert result.dtype == np.object_

    @pytest.mark.parametrize("dtype", ["float32", "float64",
                                       "int64", "int32"])
    def test_astype(self, dtype):
        s = Series(np.random.randn(5), name='foo')
        as_typed = s.astype(dtype)

        assert as_typed.dtype == dtype
        assert as_typed.name == s.name

    def test_asobject_deprecated(self):
        s = Series(np.random.randn(5), name='foo')
        with tm.assert_produces_warning(FutureWarning):
            o = s.asobject
        assert isinstance(o, np.ndarray)

    def test_dtype(self):

        assert self.ts.dtype == np.dtype('float64')
        assert self.ts.dtypes == np.dtype('float64')
        assert self.ts.ftype == 'float64:dense'
        assert self.ts.ftypes == 'float64:dense'
        tm.assert_series_equal(self.ts.get_dtype_counts(),
                               Series(1, ['float64']))
        # GH18243 - Assert .get_ftype_counts is deprecated
        with tm.assert_produces_warning(FutureWarning):
            tm.assert_series_equal(self.ts.get_ftype_counts(),
                                   Series(1, ['float64:dense']))

    @pytest.mark.parametrize("value", [np.nan, np.inf])
    @pytest.mark.parametrize("dtype", [np.int32, np.int64])
    def test_astype_cast_nan_inf_int(self, dtype, value):
        # gh-14265: check NaN and inf raise error when converting to int
        msg = 'Cannot convert non-finite values \\(NA or inf\\) to integer'
        s = Series([value])

        with tm.assert_raises_regex(ValueError, msg):
            s.astype(dtype)

    @pytest.mark.parametrize("dtype", [int, np.int8, np.int64])
    def test_astype_cast_object_int_fail(self, dtype):
        arr = Series(["car", "house", "tree", "1"])
        with pytest.raises(ValueError):
            arr.astype(dtype)

    def test_astype_cast_object_int(self):
        arr = Series(['1', '2', '3', '4'], dtype=object)
        result = arr.astype(int)

        tm.assert_series_equal(result, Series(np.arange(1, 5)))

    def test_astype_datetime(self):
        s = Series(tslib.iNaT, dtype='M8[ns]', index=lrange(5))

        s = s.astype('O')
        assert s.dtype == np.object_

        s = Series([datetime(2001, 1, 2, 0, 0)])

        s = s.astype('O')
        assert s.dtype == np.object_

        s = Series([datetime(2001, 1, 2, 0, 0) for i in range(3)])

        s[1] = np.nan
        assert s.dtype == 'M8[ns]'

        s = s.astype('O')
        assert s.dtype == np.object_

    def test_astype_datetime64tz(self):
        s = Series(date_range('20130101', periods=3, tz='US/Eastern'))

        # astype
        result = s.astype(object)
        expected = Series(s.astype(object), dtype=object)
        tm.assert_series_equal(result, expected)

        result = Series(s.values).dt.tz_localize('UTC').dt.tz_convert(s.dt.tz)
        tm.assert_series_equal(result, s)

        # astype - object, preserves on construction
        result = Series(s.astype(object))
        expected = s.astype(object)
        tm.assert_series_equal(result, expected)

        # astype - datetime64[ns, tz]
        result = Series(s.values).astype('datetime64[ns, US/Eastern]')
        tm.assert_series_equal(result, s)

        result = Series(s.values).astype(s.dtype)
        tm.assert_series_equal(result, s)

        result = s.astype('datetime64[ns, CET]')
        expected = Series(date_range('20130101 06:00:00', periods=3, tz='CET'))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [compat.text_type, np.str_])
    @pytest.mark.parametrize("series", [Series([string.digits * 10,
                                                tm.rands(63),
                                                tm.rands(64),
                                                tm.rands(1000)]),
                                        Series([string.digits * 10,
                                                tm.rands(63),
                                                tm.rands(64), nan, 1.0])])
    def test_astype_str_map(self, dtype, series):
        # see gh-4405
        result = series.astype(dtype)
        expected = series.map(compat.text_type)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [str, compat.text_type])
    def test_astype_str_cast(self, dtype):
        # see gh-9757: test str and unicode on python 2.x
        # and just str on python 3.x
        ts = Series([Timestamp('2010-01-04 00:00:00')])
        s = ts.astype(dtype)

        expected = Series([dtype('2010-01-04')])
        tm.assert_series_equal(s, expected)

        ts = Series([Timestamp('2010-01-04 00:00:00', tz='US/Eastern')])
        s = ts.astype(dtype)

        expected = Series([dtype('2010-01-04 00:00:00-05:00')])
        tm.assert_series_equal(s, expected)

        td = Series([Timedelta(1, unit='d')])
        s = td.astype(dtype)

        expected = Series([dtype('1 days 00:00:00.000000000')])
        tm.assert_series_equal(s, expected)

    def test_astype_unicode(self):
        # see gh-7758: A bit of magic is required to set
        # default encoding to utf-8
        digits = string.digits
        test_series = [
            Series([digits * 10, tm.rands(63), tm.rands(64), tm.rands(1000)]),
            Series([u('データーサイエンス、お前はもう死んでいる')]),
        ]

        former_encoding = None

        if not compat.PY3:
            # In Python, we can force the default encoding for this test
            former_encoding = sys.getdefaultencoding()
            reload(sys)  # noqa

            sys.setdefaultencoding("utf-8")
        if sys.getdefaultencoding() == "utf-8":
            test_series.append(Series([u('野菜食べないとやばい')
                                       .encode("utf-8")]))

        for s in test_series:
            res = s.astype("unicode")
            expec = s.map(compat.text_type)
            tm.assert_series_equal(res, expec)

        # Restore the former encoding
        if former_encoding is not None and former_encoding != "utf-8":
            reload(sys)  # noqa
            sys.setdefaultencoding(former_encoding)

    @pytest.mark.parametrize("dtype_class", [dict, Series])
    def test_astype_dict_like(self, dtype_class):
        # see gh-7271
        s = Series(range(0, 10, 2), name='abc')

        dt1 = dtype_class({'abc': str})
        result = s.astype(dt1)
        expected = Series(['0', '2', '4', '6', '8'], name='abc')
        tm.assert_series_equal(result, expected)

        dt2 = dtype_class({'abc': 'float64'})
        result = s.astype(dt2)
        expected = Series([0.0, 2.0, 4.0, 6.0, 8.0], dtype='float64',
                          name='abc')
        tm.assert_series_equal(result, expected)

        dt3 = dtype_class({'abc': str, 'def': str})
        with pytest.raises(KeyError):
            s.astype(dt3)

        dt4 = dtype_class({0: str})
        with pytest.raises(KeyError):
            s.astype(dt4)

        # GH16717
        # if dtypes provided is empty, it should error
        dt5 = dtype_class({})
        with pytest.raises(KeyError):
            s.astype(dt5)

    def test_astype_categories_deprecation(self):

        # deprecated 17636
        s = Series(['a', 'b', 'a'])
        expected = s.astype(CategoricalDtype(['a', 'b'], ordered=True))
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = s.astype('category', categories=['a', 'b'], ordered=True)
        tm.assert_series_equal(result, expected)

    def test_astype_from_categorical(self):
        l = ["a", "b", "c", "a"]
        s = Series(l)
        exp = Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        l = [1, 2, 3, 1]
        s = Series(l)
        exp = Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        df = DataFrame({"cats": [1, 2, 3, 4, 5, 6],
                        "vals": [1, 2, 3, 4, 5, 6]})
        cats = Categorical([1, 2, 3, 4, 5, 6])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        df = DataFrame({"cats": ['a', 'b', 'b', 'a', 'a', 'd'],
                        "vals": [1, 2, 3, 4, 5, 6]})
        cats = Categorical(['a', 'b', 'b', 'a', 'a', 'd'])
        exp_df = DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        # with keywords
        l = ["a", "b", "c", "a"]
        s = Series(l)
        exp = Series(Categorical(l, ordered=True))
        res = s.astype(CategoricalDtype(None, ordered=True))
        tm.assert_series_equal(res, exp)

        exp = Series(Categorical(l, categories=list('abcdef'), ordered=True))
        res = s.astype(CategoricalDtype(list('abcdef'), ordered=True))
        tm.assert_series_equal(res, exp)

    def test_astype_categorical_to_other(self):

        df = DataFrame({'value': np.random.randint(0, 10000, 100)})
        labels = ["{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500)]
        cat_labels = Categorical(labels, labels)

        df = df.sort_values(by=['value'], ascending=True)
        df['value_group'] = pd.cut(df.value, range(0, 10500, 500),
                                   right=False, labels=cat_labels)

        s = df['value_group']
        expected = s
        tm.assert_series_equal(s.astype('category'), expected)
        tm.assert_series_equal(s.astype(CategoricalDtype()), expected)
        pytest.raises(ValueError, lambda: s.astype('float64'))

        cat = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c']))
        exp = Series(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        tm.assert_series_equal(cat.astype('str'), exp)
        s2 = Series(Categorical(['1', '2', '3', '4']))
        exp2 = Series([1, 2, 3, 4]).astype(int)
        tm.assert_series_equal(s2.astype('int'), exp2)

        # object don't sort correctly, so just compare that we have the same
        # values
        def cmp(a, b):
            tm.assert_almost_equal(
                np.sort(np.unique(a)), np.sort(np.unique(b)))

        expected = Series(np.array(s.values), name='value_group')
        cmp(s.astype('object'), expected)
        cmp(s.astype(np.object_), expected)

        # array conversion
        tm.assert_almost_equal(np.array(s), np.array(s.values))

        # valid conversion
        for valid in [lambda x: x.astype('category'),
                      lambda x: x.astype(CategoricalDtype()),
                      lambda x: x.astype('object').astype('category'),
                      lambda x: x.astype('object').astype(
                          CategoricalDtype())
                      ]:

            result = valid(s)
            # compare series values
            # internal .categories can't be compared because it is sorted
            tm.assert_series_equal(result, s, check_categorical=False)

        # invalid conversion (these are NOT a dtype)
        for invalid in [lambda x: x.astype(Categorical),
                        lambda x: x.astype('object').astype(Categorical)]:
            pytest.raises(TypeError, lambda: invalid(s))

    @pytest.mark.parametrize('name', [None, 'foo'])
    @pytest.mark.parametrize('dtype_ordered', [True, False])
    @pytest.mark.parametrize('series_ordered', [True, False])
    def test_astype_categorical_to_categorical(self, name, dtype_ordered,
                                               series_ordered):
        # GH 10696/18593
        s_data = list('abcaacbab')
        s_dtype = CategoricalDtype(list('bac'), ordered=series_ordered)
        s = Series(s_data, dtype=s_dtype, name=name)

        # unspecified categories
        dtype = CategoricalDtype(ordered=dtype_ordered)
        result = s.astype(dtype)
        exp_dtype = CategoricalDtype(s_dtype.categories, dtype_ordered)
        expected = Series(s_data, name=name, dtype=exp_dtype)
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = s.astype('category', ordered=dtype_ordered)
        tm.assert_series_equal(result, expected)

        # different categories
        dtype = CategoricalDtype(list('adc'), dtype_ordered)
        result = s.astype(dtype)
        expected = Series(s_data, name=name, dtype=dtype)
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = s.astype(
                'category', categories=list('adc'), ordered=dtype_ordered)
        tm.assert_series_equal(result, expected)

        if dtype_ordered is False:
            # not specifying ordered, so only test once
            expected = s
            result = s.astype('category')
            tm.assert_series_equal(result, expected)

    def test_astype_categoricaldtype(self):
        s = Series(['a', 'b', 'a'])
        result = s.astype(CategoricalDtype(['a', 'b'], ordered=True))
        expected = Series(Categorical(['a', 'b', 'a'], ordered=True))
        tm.assert_series_equal(result, expected)

        result = s.astype(CategoricalDtype(['a', 'b'], ordered=False))
        expected = Series(Categorical(['a', 'b', 'a'], ordered=False))
        tm.assert_series_equal(result, expected)

        result = s.astype(CategoricalDtype(['a', 'b', 'c'], ordered=False))
        expected = Series(Categorical(['a', 'b', 'a'],
                                      categories=['a', 'b', 'c'],
                                      ordered=False))
        tm.assert_series_equal(result, expected)
        tm.assert_index_equal(result.cat.categories, Index(['a', 'b', 'c']))

    def test_astype_categoricaldtype_with_args(self):
        s = Series(['a', 'b'])
        type_ = CategoricalDtype(['a', 'b'])

        with pytest.raises(TypeError):
            s.astype(type_, ordered=True)
        with pytest.raises(TypeError):
            s.astype(type_, categories=['a', 'b'])
        with pytest.raises(TypeError):
            s.astype(type_, categories=['a', 'b'], ordered=False)

    def test_astype_generic_timestamp_deprecated(self):
        # see gh-15524
        data = [1]

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            s = Series(data)
            dtype = np.datetime64
            result = s.astype(dtype)
            expected = Series(data, dtype=dtype)
            tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            s = Series(data)
            dtype = np.timedelta64
            result = s.astype(dtype)
            expected = Series(data, dtype=dtype)
            tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", np.typecodes['All'])
    def test_astype_empty_constructor_equality(self, dtype):
        # see gh-15524

        if dtype not in ('S', 'V'):  # poor support (if any) currently
            with warnings.catch_warnings(record=True):
                # Generic timestamp dtypes ('M' and 'm') are deprecated,
                # but we test that already in series/test_constructors.py

                init_empty = Series([], dtype=dtype)
                as_type_empty = Series([]).astype(dtype)
                tm.assert_series_equal(init_empty, as_type_empty)

    def test_complex(self):
        # see gh-4819: complex access for ndarray compat
        a = np.arange(5, dtype=np.float64)
        b = Series(a + 4j * a)

        tm.assert_numpy_array_equal(a, b.real)
        tm.assert_numpy_array_equal(4 * a, b.imag)

        b.real = np.arange(5) + 5
        tm.assert_numpy_array_equal(a + 5, b.real)
        tm.assert_numpy_array_equal(4 * a, b.imag)

    def test_arg_for_errors_in_astype(self):
        # see gh-14878
        s = Series([1, 2, 3])

        with pytest.raises(ValueError):
            s.astype(np.float64, errors=False)

        with tm.assert_produces_warning(FutureWarning):
            s.astype(np.int8, raise_on_error=True)

        s.astype(np.int8, errors='raise')

    def test_intercept_astype_object(self):
        series = Series(date_range('1/1/2000', periods=10))

        # This test no longer makes sense, as
        # Series is by default already M8[ns].
        expected = series.astype('object')

        df = DataFrame({'a': series,
                        'b': np.random.randn(len(series))})
        exp_dtypes = Series([np.dtype('datetime64[ns]'),
                             np.dtype('float64')], index=['a', 'b'])
        tm.assert_series_equal(df.dtypes, exp_dtypes)

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({'a': series, 'b': ['foo'] * len(series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

    def test_series_to_categorical(self):
        # see gh-16524: test conversion of Series to Categorical
        series = Series(['a', 'b', 'c'])

        result = Series(series, dtype='category')
        expected = Series(['a', 'b', 'c'], dtype='category')

        tm.assert_series_equal(result, expected)

    def test_infer_objects_series(self):
        # GH 11221
        actual = Series(np.array([1, 2, 3], dtype='O')).infer_objects()
        expected = Series([1, 2, 3])
        tm.assert_series_equal(actual, expected)

        actual = Series(np.array([1, 2, 3, None], dtype='O')).infer_objects()
        expected = Series([1., 2., 3., np.nan])
        tm.assert_series_equal(actual, expected)

        # only soft conversions, unconvertable pass thru unchanged
        actual = (Series(np.array([1, 2, 3, None, 'a'], dtype='O'))
                  .infer_objects())
        expected = Series([1, 2, 3, None, 'a'])

        assert actual.dtype == 'object'
        tm.assert_series_equal(actual, expected)
