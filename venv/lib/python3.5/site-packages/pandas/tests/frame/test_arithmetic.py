# -*- coding: utf-8 -*-
import pytest
import numpy as np

from pandas.compat import range

import pandas as pd
import pandas.util.testing as tm


# -------------------------------------------------------------------
# Comparisons

class TestFrameComparisons(object):
    def test_df_boolean_comparison_error(self):
        # GH#4576
        # boolean comparisons with a tuple/list give unexpected results
        df = pd.DataFrame(np.arange(6).reshape((3, 2)))

        # not shape compatible
        with pytest.raises(ValueError):
            df == (2, 2)
        with pytest.raises(ValueError):
            df == [2, 2]

    def test_df_float_none_comparison(self):
        df = pd.DataFrame(np.random.randn(8, 3), index=range(8),
                          columns=['A', 'B', 'C'])

        with pytest.raises(TypeError):
            df.__eq__(None)

    def test_df_string_comparison(self):
        df = pd.DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        tm.assert_frame_equal(df[mask_a], df.loc[1:1, :])
        tm.assert_frame_equal(df[-mask_a], df.loc[0:0, :])

        mask_b = df.b == "foo"
        tm.assert_frame_equal(df[mask_b], df.loc[0:0, :])
        tm.assert_frame_equal(df[-mask_b], df.loc[1:1, :])

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types(self, opname):
        # GH#15077, non-empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        result = getattr(df, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types_empty(self, opname):
        # GH#15077 empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        empty = df.iloc[:0]
        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))

    @pytest.mark.parametrize('timestamps', [
        [pd.Timestamp('2012-01-01 13:00:00+00:00')] * 2,
        [pd.Timestamp('2012-01-01 13:00:00')] * 2])
    def test_tz_aware_scalar_comparison(self, timestamps):
        # Test for issue #15966
        df = pd.DataFrame({'test': timestamps})
        expected = pd.DataFrame({'test': [False, False]})
        tm.assert_frame_equal(df == -1, expected)


# -------------------------------------------------------------------
# Arithmetic

class TestFrameFlexArithmetic(object):
    def test_df_add_flex_filled_mixed_dtypes(self):
        # GH#19611
        dti = pd.date_range('2016-01-01', periods=3)
        ser = pd.Series(['1 Day', 'NaT', '2 Days'], dtype='timedelta64[ns]')
        df = pd.DataFrame({'A': dti, 'B': ser})
        other = pd.DataFrame({'A': ser, 'B': ser})
        fill = pd.Timedelta(days=1).to_timedelta64()
        result = df.add(other, fill_value=fill)

        expected = pd.DataFrame(
            {'A': pd.Series(['2016-01-02', '2016-01-03', '2016-01-05'],
                            dtype='datetime64[ns]'),
             'B': ser * 2})
        tm.assert_frame_equal(result, expected)


class TestFrameMulDiv(object):
    """Tests for DataFrame multiplication and division"""
    # ------------------------------------------------------------------
    # Mod By Zero

    def test_df_mod_zero_df(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        # this is technically wrong, as the integer portion is coerced to float
        # ###
        first = pd.Series([0, 0, 0, 0], dtype='float64')
        second = pd.Series([np.nan, np.nan, np.nan, 0])
        expected = pd.DataFrame({'first': first, 'second': second})
        result = df % df
        tm.assert_frame_equal(result, expected)

    def test_df_mod_zero_array(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        # this is technically wrong, as the integer portion is coerced to float
        # ###
        first = pd.Series([0, 0, 0, 0], dtype='float64')
        second = pd.Series([np.nan, np.nan, np.nan, 0])
        expected = pd.DataFrame({'first': first, 'second': second})

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values % df.values
        result2 = pd.DataFrame(arr, index=df.index,
                               columns=df.columns, dtype='float64')
        result2.iloc[0:3, 1] = np.nan
        tm.assert_frame_equal(result2, expected)

    def test_df_mod_zero_int(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        result = df % 0
        expected = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        tm.assert_frame_equal(result, expected)

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values.astype('float64') % 0
        result2 = pd.DataFrame(arr, index=df.index, columns=df.columns)
        tm.assert_frame_equal(result2, expected)

    def test_df_mod_zero_series_does_not_commute(self):
        # GH#3590, modulo as ints
        # not commutative with series
        df = pd.DataFrame(np.random.randn(10, 5))
        ser = df[0]
        res = ser % df
        res2 = df % ser
        assert not res.fillna(0).equals(res2.fillna(0))

    # ------------------------------------------------------------------
    # Division By Zero

    def test_df_div_zero_df(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
        result = df / df

        first = pd.Series([1.0, 1.0, 1.0, 1.0])
        second = pd.Series([np.nan, np.nan, np.nan, 1])
        expected = pd.DataFrame({'first': first, 'second': second})
        tm.assert_frame_equal(result, expected)

    def test_df_div_zero_array(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        first = pd.Series([1.0, 1.0, 1.0, 1.0])
        second = pd.Series([np.nan, np.nan, np.nan, 1])
        expected = pd.DataFrame({'first': first, 'second': second})

        with np.errstate(all='ignore'):
            arr = df.values.astype('float') / df.values
        result = pd.DataFrame(arr, index=df.index,
                              columns=df.columns)
        tm.assert_frame_equal(result, expected)

    def test_df_div_zero_int(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        result = df / 0
        expected = pd.DataFrame(np.inf, index=df.index, columns=df.columns)
        expected.iloc[0:3, 1] = np.nan
        tm.assert_frame_equal(result, expected)

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values.astype('float64') / 0
        result2 = pd.DataFrame(arr, index=df.index,
                               columns=df.columns)
        tm.assert_frame_equal(result2, expected)

    def test_df_div_zero_series_does_not_commute(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame(np.random.randn(10, 5))
        ser = df[0]
        res = ser / df
        res2 = df / ser
        assert not res.fillna(0).equals(res2.fillna(0))


class TestFrameArithmetic(object):

    @pytest.mark.xfail(reason='GH#7996 datetime64 units not converted to nano')
    def test_df_sub_datetime64_not_ns(self):
        df = pd.DataFrame(pd.date_range('20130101', periods=3))
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'
        res = df - dt64
        expected = pd.DataFrame([pd.Timedelta(days=0), pd.Timedelta(days=1),
                                 pd.Timedelta(days=2)])
        tm.assert_frame_equal(res, expected)

    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_radd_str_invalid(self, dtype, data):
        df = pd.DataFrame(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + df

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_int(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([2, 3, 4], dtype=dtype)
        result = 1 + df
        tm.assert_frame_equal(result, expected)
        result = df + 1
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_nan(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([np.nan, np.nan, np.nan], dtype=dtype)
        result = np.nan + df
        tm.assert_frame_equal(result, expected)
        result = df + np.nan
        tm.assert_frame_equal(result, expected)

    def test_df_radd_str(self):
        df = pd.DataFrame(['x', np.nan, 'x'])
        tm.assert_frame_equal('a' + df, pd.DataFrame(['ax', np.nan, 'ax']))
        tm.assert_frame_equal(df + 'a', pd.DataFrame(['xa', np.nan, 'xa']))


class TestPeriodFrameArithmetic(object):

    def test_ops_frame_period(self):
        # GH 13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        assert df['A'].dtype == object
        assert df['B'].dtype == object

        p = pd.Period('2015-03', freq='M')
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2, 1], dtype=object),
                            'B': np.array([14, 13], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -1 * exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        assert df2['A'].dtype == object
        assert df2['B'].dtype == object

        exp = pd.DataFrame({'A': np.array([4, 4], dtype=object),
                            'B': np.array([16, 16], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -1 * exp)
