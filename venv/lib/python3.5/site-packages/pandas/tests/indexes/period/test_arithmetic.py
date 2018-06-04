# -*- coding: utf-8 -*-
from datetime import timedelta
import operator

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import (Timedelta,
                    period_range, Period, PeriodIndex,
                    _np_version_under1p10)
import pandas.core.indexes.period as period
from pandas.core import ops
from pandas.errors import PerformanceWarning


_common_mismatch = [pd.offsets.YearBegin(2),
                    pd.offsets.MonthBegin(1),
                    pd.offsets.Minute()]


@pytest.fixture(params=[timedelta(minutes=30),
                        np.timedelta64(30, 's'),
                        Timedelta(seconds=30)] + _common_mismatch)
def not_hourly(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Hourly frequencies.
    """
    return request.param


@pytest.fixture(params=[np.timedelta64(4, 'h'),
                        timedelta(hours=23),
                        Timedelta('23:00:00')] + _common_mismatch)
def not_daily(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Daily frequencies.
    """
    return request.param


@pytest.fixture(params=[np.timedelta64(365, 'D'),
                        timedelta(365),
                        Timedelta(days=365)] + _common_mismatch)
def mismatched(request):
    """
    Several timedelta-like and DateOffset instances that are _not_
    compatible with Monthly or Annual frequencies.
    """
    return request.param


@pytest.fixture(params=[pd.offsets.Day(3),
                        timedelta(days=3),
                        np.timedelta64(3, 'D'),
                        pd.offsets.Hour(72),
                        timedelta(minutes=60 * 24 * 3),
                        np.timedelta64(72, 'h'),
                        Timedelta('72:00:00')])
def three_days(request):
    """
    Several timedelta-like and DateOffset objects that each represent
    a 3-day timedelta
    """
    return request.param


@pytest.fixture(params=[pd.offsets.Hour(2),
                        timedelta(hours=2),
                        np.timedelta64(2, 'h'),
                        pd.offsets.Minute(120),
                        timedelta(minutes=120),
                        np.timedelta64(120, 'm')])
def two_hours(request):
    """
    Several timedelta-like and DateOffset objects that each represent
    a 2-hour timedelta
    """
    return request.param


class TestPeriodIndexComparisons(object):
    def test_pi_cmp_period(self):
        idx = period_range('2007-01', periods=20, freq='M')

        result = idx < idx[10]
        exp = idx.values < idx.values[10]
        tm.assert_numpy_array_equal(result, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_pi(self, freq):
        base = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                           freq=freq)
        per = Period('2011-02', freq=freq)

        exp = np.array([False, True, False, False])
        tm.assert_numpy_array_equal(base == per, exp)
        tm.assert_numpy_array_equal(per == base, exp)

        exp = np.array([True, False, True, True])
        tm.assert_numpy_array_equal(base != per, exp)
        tm.assert_numpy_array_equal(per != base, exp)

        exp = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(base > per, exp)
        tm.assert_numpy_array_equal(per < base, exp)

        exp = np.array([True, False, False, False])
        tm.assert_numpy_array_equal(base < per, exp)
        tm.assert_numpy_array_equal(per > base, exp)

        exp = np.array([False, True, True, True])
        tm.assert_numpy_array_equal(base >= per, exp)
        tm.assert_numpy_array_equal(per <= base, exp)

        exp = np.array([True, True, False, False])
        tm.assert_numpy_array_equal(base <= per, exp)
        tm.assert_numpy_array_equal(per >= base, exp)

        idx = PeriodIndex(['2011-02', '2011-01', '2011-03', '2011-05'],
                          freq=freq)

        exp = np.array([False, False, True, False])
        tm.assert_numpy_array_equal(base == idx, exp)

        exp = np.array([True, True, False, True])
        tm.assert_numpy_array_equal(base != idx, exp)

        exp = np.array([False, True, False, False])
        tm.assert_numpy_array_equal(base > idx, exp)

        exp = np.array([True, False, False, True])
        tm.assert_numpy_array_equal(base < idx, exp)

        exp = np.array([False, True, True, False])
        tm.assert_numpy_array_equal(base >= idx, exp)

        exp = np.array([True, False, True, True])
        tm.assert_numpy_array_equal(base <= idx, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_pi_mismatched_freq_raises(self, freq):
        # different base freq
        base = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                           freq=freq)

        msg = "Input has different freq=A-DEC from PeriodIndex"
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= Period('2011', freq='A')

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            Period('2011', freq='A') >= base

        idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='A')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= idx

        # Different frequency
        msg = "Input has different freq=4M from PeriodIndex"
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= Period('2011', freq='4M')

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            Period('2011', freq='4M') >= base

        idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='4M')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= idx

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_nat(self, freq):
        idx1 = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

        result = idx1 > Period('2011-02', freq=freq)
        exp = np.array([False, False, False, True])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('2011-02', freq=freq) < idx1
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == Period('NaT', freq=freq)
        exp = np.array([False, False, False, False])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('NaT', freq=freq) == idx1
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != Period('NaT', freq=freq)
        exp = np.array([True, True, True, True])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('NaT', freq=freq) != idx1
        tm.assert_numpy_array_equal(result, exp)

        idx2 = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq=freq)
        result = idx1 < idx2
        exp = np.array([True, False, False, False])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == idx2
        exp = np.array([False, False, False, False])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != idx2
        exp = np.array([True, True, True, True])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == idx1
        exp = np.array([True, True, False, True])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != idx1
        exp = np.array([False, False, True, False])
        tm.assert_numpy_array_equal(result, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_nat_mismatched_freq_raises(self, freq):
        idx1 = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

        diff = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq='4M')
        msg = "Input has different freq=4M from PeriodIndex"
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            idx1 > diff

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            idx1 == diff

    # TODO: De-duplicate with test_pi_cmp_nat
    @pytest.mark.parametrize('dtype', [object, None])
    def test_comp_nat(self, dtype):
        left = pd.PeriodIndex([pd.Period('2011-01-01'), pd.NaT,
                               pd.Period('2011-01-03')])
        right = pd.PeriodIndex([pd.NaT, pd.NaT, pd.Period('2011-01-03')])

        if dtype is not None:
            left = left.astype(dtype)
            right = right.astype(dtype)

        result = left == right
        expected = np.array([False, False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = left != right
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

        expected = np.array([False, False, False])
        tm.assert_numpy_array_equal(left == pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT == right, expected)

        expected = np.array([True, True, True])
        tm.assert_numpy_array_equal(left != pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT != left, expected)

        expected = np.array([False, False, False])
        tm.assert_numpy_array_equal(left < pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT > left, expected)


class TestPeriodIndexArithmetic(object):

    # -------------------------------------------------------------
    # Invalid Operations

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_pi_add_sub_float(self, op, other):
        dti = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        pi = dti.to_period('D')
        with pytest.raises(TypeError):
            op(pi, other)

    # -----------------------------------------------------------------
    # __add__/__sub__ with ndarray[datetime64] and ndarray[timedelta64]

    def test_pi_add_sub_dt64_array_raises(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=3)
        dti = pd.date_range('2016-01-01', periods=3)
        dtarr = dti.values

        with pytest.raises(TypeError):
            rng + dtarr
        with pytest.raises(TypeError):
            dtarr + rng

        with pytest.raises(TypeError):
            rng - dtarr
        with pytest.raises(TypeError):
            dtarr - rng

    def test_pi_add_sub_td64_array_non_tick_raises(self):
        rng = pd.period_range('1/1/2000', freq='Q', periods=3)
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        with pytest.raises(period.IncompatibleFrequency):
            rng + tdarr
        with pytest.raises(period.IncompatibleFrequency):
            tdarr + rng

        with pytest.raises(period.IncompatibleFrequency):
            rng - tdarr
        with pytest.raises(period.IncompatibleFrequency):
            tdarr - rng

    @pytest.mark.xfail(reason='op with TimedeltaIndex raises, with ndarray OK')
    def test_pi_add_sub_td64_array_tick(self):
        rng = pd.period_range('1/1/2000', freq='Q', periods=3)
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        expected = rng + tdi
        result = rng + tdarr
        tm.assert_index_equal(result, expected)
        result = tdarr + rng
        tm.assert_index_equal(result, expected)

        expected = rng - tdi
        result = rng - tdarr
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            tdarr - rng

    # -----------------------------------------------------------------
    # operations with array/Index of DateOffset objects

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_pi_add_offset_array(self, box):
        # GH#18849
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        offs = box([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                    pd.offsets.QuarterEnd(n=-2, startingMonth=12)])
        expected = pd.PeriodIndex([pd.Period('2015Q2'), pd.Period('2015Q4')])

        with tm.assert_produces_warning(PerformanceWarning):
            res = pi + offs
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = offs + pi
        tm.assert_index_equal(res2, expected)

        unanchored = np.array([pd.offsets.Hour(n=1),
                               pd.offsets.Minute(n=-2)])
        # addition/subtraction ops with incompatible offsets should issue
        # a PerformanceWarning and _then_ raise a TypeError.
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                pi + unanchored
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                unanchored + pi

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_pi_sub_offset_array(self, box):
        # GH#18824
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        other = box([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                     pd.offsets.QuarterEnd(n=-2, startingMonth=12)])

        expected = PeriodIndex([pi[n] - other[n] for n in range(len(pi))])

        with tm.assert_produces_warning(PerformanceWarning):
            res = pi - other
        tm.assert_index_equal(res, expected)

        anchored = box([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        # addition/subtraction ops with anchored offsets should issue
        # a PerformanceWarning and _then_ raise a TypeError.
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                pi - anchored
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                anchored - pi

    def test_pi_add_iadd_pi_raises(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        # previously performed setop union, now raises TypeError (GH14164)
        with pytest.raises(TypeError):
            rng + other

        with pytest.raises(TypeError):
            rng += other

    def test_pi_add_iadd_int(self, one):
        # Variants of `one` for #19012
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng + one
        expected = pd.period_range('2000-01-01 10:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_pi_sub_isub_int(self, one):
        """
        PeriodIndex.__sub__ and __isub__ with several representations of
        the integer 1, e.g. int, long, np.int64, np.uint8, ...
        """
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng - one
        expected = pd.period_range('2000-01-01 08:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= one
        tm.assert_index_equal(rng, expected)

    @pytest.mark.parametrize('five', [5, np.array(5, dtype=np.int64)])
    def test_pi_sub_intlike(self, five):
        rng = period_range('2007-01', periods=50)

        result = rng - five
        exp = rng + (-five)
        tm.assert_index_equal(result, exp)

    def test_pi_sub_isub_pi_raises(self):
        # previously performed setop, now raises TypeError (GH14164)
        # TODO needs to wait on #13077 for decision on result type
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        with pytest.raises(TypeError):
            rng - other

        with pytest.raises(TypeError):
            rng -= other

    def test_pi_sub_isub_offset(self):
        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng - pd.offsets.YearEnd(5)
        expected = pd.period_range('2009', '2019', freq='A')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng - pd.offsets.MonthEnd(5)
        expected = pd.period_range('2013-08', '2016-07', freq='M')
        tm.assert_index_equal(result, expected)

        rng -= pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

    # ---------------------------------------------------------------
    # Timedelta-like (timedelta, timedelta64, Timedelta, Tick)
    # TODO: Some of these are misnomers because of non-Tick DateOffsets

    def test_pi_add_iadd_timedeltalike_daily(self, three_days):
        # Tick
        other = three_days
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        expected = pd.period_range('2014-05-04', '2014-05-18', freq='D')

        result = rng + other
        tm.assert_index_equal(result, expected)

        rng += other
        tm.assert_index_equal(rng, expected)

    def test_pi_sub_isub_timedeltalike_daily(self, three_days):
        # Tick-like 3 Days
        other = three_days
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        expected = pd.period_range('2014-04-28', '2014-05-12', freq='D')

        result = rng - other
        tm.assert_index_equal(result, expected)

        rng -= other
        tm.assert_index_equal(rng, expected)

    def test_pi_add_iadd_timedeltalike_freq_mismatch_daily(self, not_daily):
        other = not_daily
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=D\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other

    def test_pi_sub_timedeltalike_freq_mismatch_daily(self, not_daily):
        other = not_daily
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=D\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other

    def test_pi_add_iadd_timedeltalike_hourly(self, two_hours):
        other = two_hours
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        expected = pd.period_range('2014-01-01 12:00', '2014-01-05 12:00',
                                   freq='H')

        result = rng + other
        tm.assert_index_equal(result, expected)

        rng += other
        tm.assert_index_equal(rng, expected)

    def test_pi_add_timedeltalike_mismatched_freq_hourly(self, not_hourly):
        other = not_hourly
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=H\\)'

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other

    def test_pi_sub_isub_timedeltalike_hourly(self, two_hours):
        other = two_hours
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        expected = pd.period_range('2014-01-01 08:00', '2014-01-05 08:00',
                                   freq='H')

        result = rng - other
        tm.assert_index_equal(result, expected)

        rng -= other
        tm.assert_index_equal(rng, expected)

    def test_add_iadd_timedeltalike_annual(self):
        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng + pd.offsets.YearEnd(5)
        expected = pd.period_range('2019', '2029', freq='A')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

    def test_pi_add_iadd_timedeltalike_freq_mismatch_annual(self, mismatched):
        other = mismatched
        rng = pd.period_range('2014', '2024', freq='A')
        msg = ('Input has different freq(=.+)? '
               'from PeriodIndex\\(freq=A-DEC\\)')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other

    def test_pi_sub_isub_timedeltalike_freq_mismatch_annual(self, mismatched):
        other = mismatched
        rng = pd.period_range('2014', '2024', freq='A')
        msg = ('Input has different freq(=.+)? '
               'from PeriodIndex\\(freq=A-DEC\\)')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng -= other

    def test_pi_add_iadd_timedeltalike_M(self):
        rng = pd.period_range('2014-01', '2016-12', freq='M')
        expected = pd.period_range('2014-06', '2017-05', freq='M')

        result = rng + pd.offsets.MonthEnd(5)
        tm.assert_index_equal(result, expected)

        rng += pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

    def test_pi_add_iadd_timedeltalike_freq_mismatch_monthly(self, mismatched):
        other = mismatched
        rng = pd.period_range('2014-01', '2016-12', freq='M')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=M\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other

    def test_pi_sub_isub_timedeltalike_freq_mismatch_monthly(self, mismatched):
        other = mismatched
        rng = pd.period_range('2014-01', '2016-12', freq='M')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=M\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng -= other

    # ---------------------------------------------------------------
    # PeriodIndex.shift is used by __add__ and __sub__

    def test_pi_shift_ndarray(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        result = idx.shift(np.array([1, 2, 3, 4]))
        expected = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)

        result = idx.shift(np.array([1, -2, 3, -4]))
        expected = PeriodIndex(['2011-02', '2010-12', 'NaT', '2010-12'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)

    def test_shift(self):
        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2002', end='12/1/2010')

        tm.assert_index_equal(pi1.shift(0), pi1)

        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2000', end='12/1/2008')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='2/1/2001', end='1/1/2010')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='12/1/2000', end='11/1/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='1/2/2001', end='12/2/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='12/31/2000', end='11/30/2009')
        assert len(pi1) == len(pi2)
        tm.assert_index_equal(pi1.shift(-1), pi2)

    def test_shift_corner_cases(self):
        # GH#9903
        idx = pd.PeriodIndex([], name='xxx', freq='H')

        with pytest.raises(TypeError):
            # period shift doesn't accept freq
            idx.shift(1, freq='H')

        tm.assert_index_equal(idx.shift(0), idx)
        tm.assert_index_equal(idx.shift(3), idx)

        idx = pd.PeriodIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                              '2011-01-01 12:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(0), idx)
        exp = pd.PeriodIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                              '2011-01-01 15:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(3), exp)
        exp = pd.PeriodIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                              '2011-01-01 09:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(-3), exp)

    def test_shift_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        result = idx.shift(1)
        expected = PeriodIndex(['2011-02', '2011-03', 'NaT', '2011-05'],
                               freq='M', name='idx')
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

    def test_shift_gh8083(self):
        # test shift for PeriodIndex
        # GH#8083
        drange = pd.period_range('20130101', periods=5, freq='D')
        result = drange.shift(1)
        expected = PeriodIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                '2013-01-05', '2013-01-06'], freq='D')
        tm.assert_index_equal(result, expected)


class TestPeriodIndexSeriesMethods(object):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        if isinstance(expected, pd.Index):
            tm.assert_index_equal(result, expected)
        else:
            # comp op results in bool
            tm.assert_numpy_array_equal(result, expected)

        ser = pd.Series(values)
        result = func(ser)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_ops(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        expected = PeriodIndex(['2011-03', '2011-04', '2011-05', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        result = idx - Period('2011-01', freq='M')
        exp = pd.Index([0, 1, 2, 3], name='idx')
        tm.assert_index_equal(result, exp)

        result = Period('2011-01', freq='M') - idx
        exp = pd.Index([0, -1, -2, -3], name='idx')
        tm.assert_index_equal(result, exp)

    @pytest.mark.parametrize('ng', ["str", 1.5])
    def test_pi_ops_errors(self, ng):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')
        ser = pd.Series(idx)

        msg = r"unsupported operand type\(s\)"

        for obj in [idx, ser]:
            with tm.assert_raises_regex(TypeError, msg):
                obj + ng

            with pytest.raises(TypeError):
                # error message differs between PY2 and 3
                ng + obj

            with tm.assert_raises_regex(TypeError, msg):
                obj - ng

            with pytest.raises(TypeError):
                np.add(obj, ng)

            if _np_version_under1p10:
                assert np.add(ng, obj) is NotImplemented
            else:
                with pytest.raises(TypeError):
                    np.add(ng, obj)

            with pytest.raises(TypeError):
                np.subtract(obj, ng)

            if _np_version_under1p10:
                assert np.subtract(ng, obj) is NotImplemented
            else:
                with pytest.raises(TypeError):
                    np.subtract(ng, obj)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        expected = PeriodIndex(['2011-03', '2011-04', 'NaT', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)
        self._check(idx, lambda x: np.add(x, 2), expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        self._check(idx + 2, lambda x: np.subtract(x, 2), idx)

        # freq with mult
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='2M', name='idx')
        expected = PeriodIndex(['2011-07', '2011-08', 'NaT', '2011-10'],
                               freq='2M', name='idx')
        self._check(idx, lambda x: x + 3, expected)
        self._check(idx, lambda x: 3 + x, expected)
        self._check(idx, lambda x: np.add(x, 3), expected)

        self._check(idx + 3, lambda x: x - 3, idx)
        self._check(idx + 3, lambda x: np.subtract(x, 3), idx)

    def test_pi_ops_array_int(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        f = lambda x: x + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.add(x, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT', '2011-06'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT', '2010-12'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.subtract(x, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT', '2011-06'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

    def test_pi_ops_offset(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        f = lambda x: x + pd.offsets.Day()
        exp = PeriodIndex(['2011-01-02', '2011-02-02', '2011-03-02',
                           '2011-04-02'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x + pd.offsets.Day(2)
        exp = PeriodIndex(['2011-01-03', '2011-02-03', '2011-03-03',
                           '2011-04-03'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - pd.offsets.Day(2)
        exp = PeriodIndex(['2010-12-30', '2011-01-30', '2011-02-27',
                           '2011-03-30'], freq='D', name='idx')
        self._check(idx, f, exp)

    def test_pi_offset_errors(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        ser = pd.Series(idx)

        # Series op is applied per Period instance, thus error is raised
        # from Period
        msg_idx = r"Input has different freq from PeriodIndex\(freq=D\)"
        msg_s = r"Input cannot be converted to Period\(freq=D\)"
        for obj, msg in [(idx, msg_idx), (ser, msg_s)]:
            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj + pd.offsets.Hour(2)

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                pd.offsets.Hour(2) + obj

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj - pd.offsets.Hour(2)

    def test_pi_sub_period(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, -11, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(idx, pd.Period('2012-01', freq='M'))
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, 11, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(pd.Period('2012-01', freq='M'), idx)
        if _np_version_under1p10:
            assert result is NotImplemented
        else:
            tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_sub_pdnat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        exp = pd.TimedeltaIndex([pd.NaT] * 4, name='idx')
        tm.assert_index_equal(pd.NaT - idx, exp)
        tm.assert_index_equal(idx - pd.NaT, exp)

    def test_pi_sub_period_nat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, np.nan, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, np.nan, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)
