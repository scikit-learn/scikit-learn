# -*- coding: utf-8 -*-
from datetime import datetime, timedelta
import operator
from decimal import Decimal

import numpy as np
import pytest

from pandas import Series, Timestamp, Timedelta, Period, NaT
from pandas._libs.tslibs.period import IncompatibleFrequency

import pandas as pd
import pandas.util.testing as tm


@pytest.fixture
def tdser():
    """
    Return a Series with dtype='timedelta64[ns]', including a NaT.
    """
    return Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')


# ------------------------------------------------------------------
# Comparisons

class TestSeriesComparison(object):
    def test_compare_invalid(self):
        # GH#8058
        # ops testing
        a = pd.Series(np.random.randn(5), name=0)
        b = pd.Series(np.random.randn(5))
        b.name = pd.Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes(self, opname):
        # GH#15115
        ser = Series([1, 3, 2], index=range(3))
        const = 2

        result = getattr(ser, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes_empty(self, opname):
        # GH#15115 empty Series case
        ser = Series([1, 3, 2], index=range(3))
        empty = ser.iloc[:0]
        const = 2

        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.le, operator.lt,
                                    operator.ge, operator.gt])
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('baz', 'baz', 'baz')])
    def test_ser_cmp_result_names(self, names, op):
        # datetime64 dtype
        dti = pd.date_range('1949-06-07 03:00:00',
                            freq='H', periods=5, name=names[0])
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # datetime64tz dtype
        dti = dti.tz_localize('US/Central')
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # timedelta64 dtype
        tdi = dti - dti.shift(1)
        ser = Series(tdi).rename(names[1])
        result = op(ser, tdi)
        assert result.name == names[2]

        # categorical
        if op in [operator.eq, operator.ne]:
            # categorical dtype comparisons raise for inequalities
            cidx = tdi.astype('category')
            ser = Series(cidx).rename(names[1])
            result = op(ser, cidx)
            assert result.name == names[2]


class TestTimestampSeriesComparison(object):
    def test_dt64ser_cmp_date_invalid(self):
        # GH#19800 datetime.date comparison raises to
        # match DatetimeIndex/Timestamp.  This also matches the behavior
        # of stdlib datetime.datetime
        ser = pd.Series(pd.date_range('20010101', periods=10), name='dates')
        date = ser.iloc[0].to_pydatetime().date()
        assert not (ser == date).any()
        assert (ser != date).all()
        with pytest.raises(TypeError):
            ser > date
        with pytest.raises(TypeError):
            ser < date
        with pytest.raises(TypeError):
            ser >= date
        with pytest.raises(TypeError):
            ser <= date

    def test_dt64ser_cmp_period_scalar(self):
        ser = Series(pd.period_range('2000-01-01', periods=10, freq='D'))
        val = Period('2000-01-04', freq='D')
        result = ser > val
        expected = Series([x > val for x in ser])
        tm.assert_series_equal(result, expected)

        val = ser[5]
        result = ser > val
        expected = Series([x > val for x in ser])
        tm.assert_series_equal(result, expected)

    def test_timestamp_compare_series(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH#4982
        ser = pd.Series(pd.date_range('20010101', periods=10), name='dates')
        s_nat = ser.copy(deep=True)

        ser[0] = pd.Timestamp('nat')
        ser[3] = pd.Timestamp('nat')

        ops = {'lt': 'gt', 'le': 'ge', 'eq': 'eq', 'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            expected = left_f(ser, pd.Timestamp('20010109'))
            result = right_f(pd.Timestamp('20010109'), ser)
            tm.assert_series_equal(result, expected)

            # nats
            expected = left_f(ser, pd.Timestamp('nat'))
            result = right_f(pd.Timestamp('nat'), ser)
            tm.assert_series_equal(result, expected)

            # compare to timestamp with series containing nats
            expected = left_f(s_nat, pd.Timestamp('20010109'))
            result = right_f(pd.Timestamp('20010109'), s_nat)
            tm.assert_series_equal(result, expected)

            # compare to nat with series containing nats
            expected = left_f(s_nat, pd.Timestamp('nat'))
            result = right_f(pd.Timestamp('nat'), s_nat)
            tm.assert_series_equal(result, expected)

    def test_timestamp_equality(self):
        # GH#11034
        ser = pd.Series([pd.Timestamp('2000-01-29 01:59:00'), 'NaT'])
        result = ser != ser
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[0]
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[1]
        tm.assert_series_equal(result, pd.Series([True, True]))

        result = ser == ser
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[0]
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[1]
        tm.assert_series_equal(result, pd.Series([False, False]))


class TestTimedeltaSeriesComparisons(object):
    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)


class TestPeriodSeriesComparisons(object):
    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_cmp_series_period_scalar(self, freq):
        # GH 13200
        base = Series([Period(x, freq=freq) for x in
                       ['2011-01', '2011-02', '2011-03', '2011-04']])
        p = Period('2011-02', freq=freq)

        exp = Series([False, True, False, False])
        tm.assert_series_equal(base == p, exp)
        tm.assert_series_equal(p == base, exp)

        exp = Series([True, False, True, True])
        tm.assert_series_equal(base != p, exp)
        tm.assert_series_equal(p != base, exp)

        exp = Series([False, False, True, True])
        tm.assert_series_equal(base > p, exp)
        tm.assert_series_equal(p < base, exp)

        exp = Series([True, False, False, False])
        tm.assert_series_equal(base < p, exp)
        tm.assert_series_equal(p > base, exp)

        exp = Series([False, True, True, True])
        tm.assert_series_equal(base >= p, exp)
        tm.assert_series_equal(p <= base, exp)

        exp = Series([True, True, False, False])
        tm.assert_series_equal(base <= p, exp)
        tm.assert_series_equal(p >= base, exp)

        # different base freq
        msg = "Input has different freq=A-DEC from Period"
        with tm.assert_raises_regex(IncompatibleFrequency, msg):
            base <= Period('2011', freq='A')

        with tm.assert_raises_regex(IncompatibleFrequency, msg):
            Period('2011', freq='A') >= base

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_cmp_series_period_series(self, freq):
        # GH#13200
        base = Series([Period(x, freq=freq) for x in
                       ['2011-01', '2011-02', '2011-03', '2011-04']])

        ser = Series([Period(x, freq=freq) for x in
                      ['2011-02', '2011-01', '2011-03', '2011-05']])

        exp = Series([False, False, True, False])
        tm.assert_series_equal(base == ser, exp)

        exp = Series([True, True, False, True])
        tm.assert_series_equal(base != ser, exp)

        exp = Series([False, True, False, False])
        tm.assert_series_equal(base > ser, exp)

        exp = Series([True, False, False, True])
        tm.assert_series_equal(base < ser, exp)

        exp = Series([False, True, True, False])
        tm.assert_series_equal(base >= ser, exp)

        exp = Series([True, False, True, True])
        tm.assert_series_equal(base <= ser, exp)

        ser2 = Series([Period(x, freq='A') for x in
                       ['2011', '2011', '2011', '2011']])

        # different base freq
        msg = "Input has different freq=A-DEC from Period"
        with tm.assert_raises_regex(IncompatibleFrequency, msg):
            base <= ser2

    def test_cmp_series_period_series_mixed_freq(self):
        # GH#13200
        base = Series([Period('2011', freq='A'),
                       Period('2011-02', freq='M'),
                       Period('2013', freq='A'),
                       Period('2011-04', freq='M')])

        ser = Series([Period('2012', freq='A'),
                      Period('2011-01', freq='M'),
                      Period('2013', freq='A'),
                      Period('2011-05', freq='M')])

        exp = Series([False, False, True, False])
        tm.assert_series_equal(base == ser, exp)

        exp = Series([True, True, False, True])
        tm.assert_series_equal(base != ser, exp)

        exp = Series([False, True, False, False])
        tm.assert_series_equal(base > ser, exp)

        exp = Series([True, False, False, True])
        tm.assert_series_equal(base < ser, exp)

        exp = Series([False, True, True, False])
        tm.assert_series_equal(base >= ser, exp)

        exp = Series([True, False, True, True])
        tm.assert_series_equal(base <= ser, exp)


# ------------------------------------------------------------------
# Arithmetic

class TestSeriesDivision(object):
    # __div__, __rdiv__, __floordiv__, __rfloordiv__
    # for non-timestamp/timedelta/period dtypes

    def test_divide_decimal(self):
        # resolves issue GH#9787
        expected = Series([Decimal(5)])

        ser = Series([Decimal(10)])
        result = ser / Decimal(2)

        tm.assert_series_equal(result, expected)

        ser = Series([Decimal(10)])
        result = ser // Decimal(2)

        tm.assert_series_equal(result, expected)

    def test_div_equiv_binop(self):
        # Test Series.div as well as Series.__div__
        # float/integer issue
        # GH#7785
        first = Series([1, 0], name='first')
        second = Series([-0.01, -0.02], name='second')
        expected = Series([-0.01, -np.inf])

        result = second.div(first)
        tm.assert_series_equal(result, expected, check_names=False)

        result = second / first
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype2', [
        np.int64, np.int32, np.int16, np.int8,
        np.float64, np.float32, np.float16,
        np.uint64, np.uint32, np.uint16, np.uint8])
    @pytest.mark.parametrize('dtype1', [np.int64, np.float64, np.uint64])
    def test_ser_div_ser(self, dtype1, dtype2):
        # no longer do integer div for any ops, but deal with the 0's
        first = Series([3, 4, 5, 8], name='first').astype(dtype1)
        second = Series([0, 0, 0, 3], name='second').astype(dtype2)

        with np.errstate(all='ignore'):
            expected = Series(first.values.astype(np.float64) / second.values,
                              dtype='float64', name=None)
        expected.iloc[0:3] = np.inf

        result = first / second
        tm.assert_series_equal(result, expected)
        assert not result.equals(second / first)

    def test_rdiv_zero_compat(self):
        # GH#8674
        zero_array = np.array([0] * 5)
        data = np.random.randn(5)
        expected = Series([0.] * 5)

        result = zero_array / Series(data)
        tm.assert_series_equal(result, expected)

        result = Series(zero_array) / data
        tm.assert_series_equal(result, expected)

        result = Series(zero_array) / Series(data)
        tm.assert_series_equal(result, expected)

    def test_div_zero_inf_signs(self):
        # GH#9144, inf signing
        ser = Series([-1, 0, 1], name='first')
        expected = Series([-np.inf, np.nan, np.inf], name='first')

        result = ser / 0
        tm.assert_series_equal(result, expected)

    def test_rdiv_zero(self):
        # GH#9144
        ser = Series([-1, 0, 1], name='first')
        expected = Series([0.0, np.nan, 0.0], name='first')

        result = 0 / ser
        tm.assert_series_equal(result, expected)

    def test_floordiv_div(self):
        # GH#9144
        ser = Series([-1, 0, 1], name='first')

        result = ser // 0
        expected = Series([-np.inf, np.nan, np.inf], name='first')
        tm.assert_series_equal(result, expected)


class TestSeriesArithmetic(object):
    # Standard, numeric, or otherwise not-Timestamp/Timedelta/Period dtypes
    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [Timestamp('2011-01-01'), Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_radd_str_invalid(self, dtype, data):
        ser = Series(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + ser

    # TODO: parametrize, better name
    def test_object_ser_add_invalid(self):
        # invalid ops
        obj_ser = tm.makeObjectSeries()
        obj_ser.name = 'objects'
        with pytest.raises(Exception):
            obj_ser + 1
        with pytest.raises(Exception):
            obj_ser + np.array(1, dtype=np.int64)
        with pytest.raises(Exception):
            obj_ser - 1
        with pytest.raises(Exception):
            obj_ser - np.array(1, dtype=np.int64)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_nan(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=dtype)

        result = np.nan + ser
        tm.assert_series_equal(result, expected)

        result = ser + np.nan
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_int(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([2, 3, 4], dtype=dtype)

        result = 1 + ser
        tm.assert_series_equal(result, expected)

        result = ser + 1
        tm.assert_series_equal(result, expected)

    def test_series_radd_str(self):
        ser = pd.Series(['x', np.nan, 'x'])
        tm.assert_series_equal('a' + ser, pd.Series(['ax', np.nan, 'ax']))
        tm.assert_series_equal(ser + 'a', pd.Series(['xa', np.nan, 'xa']))

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_timedelta(self, dtype):
        # note this test is _not_ aimed at timedelta64-dtyped Series
        ser = pd.Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                         pd.Timedelta('3 days')], dtype=dtype)
        expected = pd.Series([pd.Timedelta('4 days'), pd.Timedelta('5 days'),
                              pd.Timedelta('6 days')])

        result = pd.Timedelta('3 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.Timedelta('3 days')
        tm.assert_series_equal(result, expected)


class TestPeriodSeriesArithmetic(object):
    def test_ops_series_timedelta(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        expected = pd.Series([pd.Period('2015-01-02', freq='D'),
                              pd.Period('2015-01-03', freq='D')], name='xxx')

        result = ser + pd.Timedelta('1 days')
        tm.assert_series_equal(result, expected)

        result = pd.Timedelta('1 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.tseries.offsets.Day()
        tm.assert_series_equal(result, expected)

        result = pd.tseries.offsets.Day() + ser
        tm.assert_series_equal(result, expected)

    def test_ops_series_period(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        per = pd.Period('2015-01-10', freq='D')
        # dtype will be object because of original dtype
        expected = pd.Series([9, 8], name='xxx', dtype=object)
        tm.assert_series_equal(per - ser, expected)
        tm.assert_series_equal(ser - per, -1 * expected)

        s2 = pd.Series([pd.Period('2015-01-05', freq='D'),
                        pd.Period('2015-01-04', freq='D')], name='xxx')
        assert s2.dtype == object

        expected = pd.Series([4, 2], name='xxx', dtype=object)
        tm.assert_series_equal(s2 - ser, expected)
        tm.assert_series_equal(ser - s2, -1 * expected)


class TestTimestampSeriesArithmetic(object):
    def test_timestamp_sub_series(self):
        ser = pd.Series(pd.date_range('2014-03-17', periods=2, freq='D',
                                      tz='US/Eastern'))
        ts = ser[0]

        delta_series = pd.Series([np.timedelta64(0, 'D'),
                                  np.timedelta64(1, 'D')])
        tm.assert_series_equal(ser - ts, delta_series)
        tm.assert_series_equal(ts - ser, -delta_series)

    def test_dt64ser_sub_datetime_dtype(self):
        ts = Timestamp(datetime(1993, 1, 7, 13, 30, 00))
        dt = datetime(1993, 6, 22, 13, 30)
        ser = Series([ts])
        result = pd.to_timedelta(np.abs(ser - dt))
        assert result.dtype == 'timedelta64[ns]'


class TestTimedeltaSeriesAdditionSubtraction(object):
    # Tests for Series[timedelta64[ns]] __add__, __sub__, __radd__, __rsub__

    # ------------------------------------------------------------------
    # Operations with int-like others

    def test_td64series_add_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            tdser + Series([2, 3, 4])

    @pytest.mark.xfail(reason='GH#19123 integer interpreted as nanoseconds')
    def test_td64series_radd_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            Series([2, 3, 4]) + tdser

    def test_td64series_sub_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            tdser - Series([2, 3, 4])

    @pytest.mark.xfail(reason='GH#19123 integer interpreted as nanoseconds')
    def test_td64series_rsub_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            Series([2, 3, 4]) - tdser

    def test_td64_series_add_intlike(self):
        # GH#19123
        tdi = pd.TimedeltaIndex(['59 days', '59 days', 'NaT'])
        ser = Series(tdi)

        other = Series([20, 30, 40], dtype='uint8')

        pytest.raises(TypeError, ser.__add__, 1)
        pytest.raises(TypeError, ser.__sub__, 1)

        pytest.raises(TypeError, ser.__add__, other)
        pytest.raises(TypeError, ser.__sub__, other)

        pytest.raises(TypeError, ser.__add__, other.values)
        pytest.raises(TypeError, ser.__sub__, other.values)

        pytest.raises(TypeError, ser.__add__, pd.Index(other))
        pytest.raises(TypeError, ser.__sub__, pd.Index(other))

    @pytest.mark.parametrize('scalar', [1, 1.5, np.array(2)])
    def test_td64series_add_sub_numeric_scalar_invalid(self, scalar, tdser):
        with pytest.raises(TypeError):
            tdser + scalar
        with pytest.raises(TypeError):
            scalar + tdser
        with pytest.raises(TypeError):
            tdser - scalar
        with pytest.raises(TypeError):
            scalar - tdser

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [
        np.array([1, 2, 3]),
        pd.Index([1, 2, 3]),
        pytest.param(Series([1, 2, 3]),
                     marks=pytest.mark.xfail(reason='GH#19123 integer '
                                                    'interpreted as nanos'))
    ])
    def test_td64series_add_sub_numeric_array_invalid(self, vector,
                                                      dtype, tdser):
        vector = vector.astype(dtype)
        with pytest.raises(TypeError):
            tdser + vector
        with pytest.raises(TypeError):
            vector + tdser
        with pytest.raises(TypeError):
            tdser - vector
        with pytest.raises(TypeError):
            vector - tdser

    # ------------------------------------------------------------------
    # Operations with datetime-like others

    def test_td64series_add_sub_timestamp(self):
        # GH#11925
        tdser = Series(pd.timedelta_range('1 day', periods=3))
        ts = Timestamp('2012-01-01')
        expected = Series(pd.date_range('2012-01-02', periods=3))
        tm.assert_series_equal(ts + tdser, expected)
        tm.assert_series_equal(tdser + ts, expected)

        expected2 = Series(pd.date_range('2011-12-31', periods=3, freq='-1D'))
        tm.assert_series_equal(ts - tdser, expected2)
        tm.assert_series_equal(ts + (-tdser), expected2)

        with pytest.raises(TypeError):
            tdser - ts

    # ------------------------------------------------------------------
    # Operations with timedelta-like others (including DateOffsets)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64_series_with_tdi(self, names):
        # GH#17250 make sure result dtype is correct
        # GH#19043 make sure names are propagated correctly
        tdi = pd.TimedeltaIndex(['0 days', '1 day'], name=names[0])
        ser = Series([Timedelta(hours=3), Timedelta(hours=4)], name=names[1])
        expected = Series([Timedelta(hours=3), Timedelta(days=1, hours=4)],
                          name=names[2])

        result = tdi + ser
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        result = ser + tdi
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        expected = Series([Timedelta(hours=-3), Timedelta(days=1, hours=-4)],
                          name=names[2])

        result = tdi - ser
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        result = ser - tdi
        tm.assert_series_equal(result, -expected)
        assert result.dtype == 'timedelta64[ns]'

    def test_td64_sub_NaT(self):
        # GH#18808
        ser = Series([NaT, Timedelta('1s')])
        res = ser - NaT
        expected = Series([NaT, NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)


class TestTimedeltaSeriesMultiplicationDivision(object):
    # Tests for Series[timedelta64[ns]]
    # __mul__, __rmul__, __div__, __rdiv__, __floordiv__, __rfloordiv__

    # ------------------------------------------------------------------
    # __floordiv__, __rfloordiv__

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_floordiv(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        result = td1 // scalar_td
        expected = Series([0, 0, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_rfloordiv(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan
        result = scalar_td // td1
        expected = Series([1, 1, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_rfloordiv_explicit(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # We can test __rfloordiv__ using this syntax,
        # see `test_timedelta_rfloordiv`
        result = td1.__rfloordiv__(scalar_td)
        expected = Series([1, 1, np.nan])
        tm.assert_series_equal(result, expected)

    # ------------------------------------------------------------------
    # Operations with int-like others

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64series_div_numeric_array(self, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)
        expected = Series(['2.95D', '1D 23H 12m', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser / vector
        tm.assert_series_equal(result, expected)

        with pytest.raises(TypeError):
            vector / tdser

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64series_mul_numeric_array(self, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * vector
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [
        np.array([20, 30, 40]),
        pytest.param(pd.Index([20, 30, 40]),
                     marks=pytest.mark.xfail(reason='__mul__ raises '
                                                    'instead of returning '
                                                    'NotImplemented')),
        Series([20, 30, 40])
    ])
    def test_td64series_rmul_numeric_array(self, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = vector * tdser
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('one', [1, np.array(1), 1.0, np.array(1.0)])
    def test_td64series_mul_numeric_scalar(self, one, tdser):
        # GH#4521
        # divide/multiply by integers
        expected = Series(['-59 Days', '-59 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * (-one)
        tm.assert_series_equal(result, expected)
        result = (-one) * tdser
        tm.assert_series_equal(result, expected)

        expected = Series(['118 Days', '118 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * (2 * one)
        tm.assert_series_equal(result, expected)
        result = (2 * one) * tdser
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('two', [
        2, 2.0,
        pytest.param(np.array(2),
                     marks=pytest.mark.xfail(reason='GH#19011 is_list_like '
                                                    'incorrectly True.')),
        pytest.param(np.array(2.0),
                     marks=pytest.mark.xfail(reason='GH#19011 is_list_like '
                                                    'incorrectly True.')),
    ])
    def test_td64series_div_numeric_scalar(self, two, tdser):
        # GH#4521
        # divide/multiply by integers
        expected = Series(['29.5D', '29.5D', 'NaT'], dtype='timedelta64[ns]')

        result = tdser / two
        tm.assert_series_equal(result, expected)

    # ------------------------------------------------------------------
    # Operations with timedelta-like others

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_tdi_mul_int_series(self, names):
        # GH#19042
        tdi = pd.TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                                name=names[0])
        ser = Series([0, 1, 2, 3, 4], dtype=np.int64, name=names[1])

        expected = Series(['0days', '1day', '4days', '9days', '16days'],
                          dtype='timedelta64[ns]',
                          name=names[2])

        result = ser * tdi
        tm.assert_series_equal(result, expected)

        # The direct operation tdi * ser still needs to be fixed.
        result = ser.__rmul__(tdi)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_float_series_rdiv_tdi(self, names):
        # GH#19042
        # TODO: the direct operation TimedeltaIndex / Series still
        # needs to be fixed.
        tdi = pd.TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                                name=names[0])
        ser = Series([1.5, 3, 4.5, 6, 7.5], dtype=np.float64, name=names[1])

        expected = Series([tdi[n] / ser[n] for n in range(len(ser))],
                          dtype='timedelta64[ns]',
                          name=names[2])

        result = ser.__rdiv__(tdi)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_td64series_mul_timedeltalike_invalid(self, scalar_td):
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot|not supported'
        with tm.assert_raises_regex(TypeError, pattern):
            td1 * scalar_td
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td * td1


class TestTimedeltaSeriesInvalidArithmeticOps(object):
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_td64series_pow_invalid(self, scalar_td):
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot|not supported'
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td ** td1
        with tm.assert_raises_regex(TypeError, pattern):
            td1 ** scalar_td
