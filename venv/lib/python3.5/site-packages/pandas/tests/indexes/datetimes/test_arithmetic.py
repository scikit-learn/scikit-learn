# -*- coding: utf-8 -*-
import warnings
from datetime import datetime, timedelta
import operator

import pytest

import numpy as np

import pandas as pd
from pandas.compat.numpy import np_datetime64_compat
import pandas.util.testing as tm
from pandas.errors import PerformanceWarning, NullFrequencyError
from pandas import (Timestamp, Timedelta, Series,
                    DatetimeIndex, TimedeltaIndex,
                    date_range)
from pandas.core import ops
from pandas._libs import tslib
from pandas._libs.tslibs.offsets import shift_months


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo',
                        'US/Eastern', 'dateutil/Asia/Singapore',
                        'dateutil/US/Pacific'])
def tz(request):
    return request.param


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=str)
def delta(request):
    # Several ways of representing two hours
    return request.param


@pytest.fixture(
    params=[
        datetime(2011, 1, 1),
        DatetimeIndex(['2011-01-01', '2011-01-02']),
        DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize('US/Eastern'),
        np.datetime64('2011-01-01'),
        Timestamp('2011-01-01')],
    ids=lambda x: type(x).__name__)
def addend(request):
    return request.param


class TestDatetimeIndexComparisons(object):
    @pytest.mark.parametrize('other', [datetime(2016, 1, 1),
                                       Timestamp('2016-01-01'),
                                       np.datetime64('2016-01-01')])
    def test_dti_cmp_datetimelike(self, other, tz):
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        if tz is not None:
            if isinstance(other, np.datetime64):
                # no tzaware version available
                return
            elif isinstance(other, Timestamp):
                other = other.tz_localize(dti.tzinfo)
            else:
                other = tslib._localize_pydatetime(other, dti.tzinfo)

        result = dti == other
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = dti > other
        expected = np.array([False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = dti >= other
        expected = np.array([True, True])
        tm.assert_numpy_array_equal(result, expected)

        result = dti < other
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)

        result = dti <= other
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

    def dti_cmp_non_datetime(self, tz):
        # GH#19301 by convention datetime.date is not considered comparable
        # to Timestamp or DatetimeIndex.  This may change in the future.
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)

        other = datetime(2016, 1, 1).date()
        assert not (dti == other).any()
        assert (dti != other).all()
        with pytest.raises(TypeError):
            dti < other
        with pytest.raises(TypeError):
            dti <= other
        with pytest.raises(TypeError):
            dti > other
        with pytest.raises(TypeError):
            dti >= other

    @pytest.mark.parametrize('other', [None, np.nan, pd.NaT])
    def test_dti_eq_null_scalar(self, other, tz):
        # GH#19301
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        assert not (dti == other).any()

    @pytest.mark.parametrize('other', [None, np.nan, pd.NaT])
    def test_dti_ne_null_scalar(self, other, tz):
        # GH#19301
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        assert (dti != other).all()

    @pytest.mark.parametrize('other', [None, np.nan])
    def test_dti_cmp_null_scalar_inequality(self, tz, other):
        # GH#19301
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)

        with pytest.raises(TypeError):
            dti < other
        with pytest.raises(TypeError):
            dti <= other
        with pytest.raises(TypeError):
            dti > other
        with pytest.raises(TypeError):
            dti >= other

    def test_dti_cmp_nat(self):
        left = pd.DatetimeIndex([pd.Timestamp('2011-01-01'), pd.NaT,
                                 pd.Timestamp('2011-01-03')])
        right = pd.DatetimeIndex([pd.NaT, pd.NaT, pd.Timestamp('2011-01-03')])

        for lhs, rhs in [(left, right),
                         (left.astype(object), right.astype(object))]:
            result = rhs == lhs
            expected = np.array([False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = lhs != rhs
            expected = np.array([True, True, False])
            tm.assert_numpy_array_equal(result, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs == pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT == rhs, expected)

            expected = np.array([True, True, True])
            tm.assert_numpy_array_equal(lhs != pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT != lhs, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs < pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT > lhs, expected)

    def test_dti_cmp_nat_behaves_like_float_cmp_nan(self):
        fidx1 = pd.Index([1.0, np.nan, 3.0, np.nan, 5.0, 7.0])
        fidx2 = pd.Index([2.0, 3.0, np.nan, np.nan, 6.0, 7.0])

        didx1 = pd.DatetimeIndex(['2014-01-01', pd.NaT, '2014-03-01', pd.NaT,
                                  '2014-05-01', '2014-07-01'])
        didx2 = pd.DatetimeIndex(['2014-02-01', '2014-03-01', pd.NaT, pd.NaT,
                                  '2014-06-01', '2014-07-01'])
        darr = np.array([np_datetime64_compat('2014-02-01 00:00Z'),
                         np_datetime64_compat('2014-03-01 00:00Z'),
                         np_datetime64_compat('nat'), np.datetime64('nat'),
                         np_datetime64_compat('2014-06-01 00:00Z'),
                         np_datetime64_compat('2014-07-01 00:00Z')])

        cases = [(fidx1, fidx2), (didx1, didx2), (didx1, darr)]

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, idx2 in cases:

                result = idx1 < idx2
                expected = np.array([True, False, False, False, True, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx2 > idx1
                expected = np.array([True, False, False, False, True, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= idx2
                expected = np.array([True, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx2 >= idx1
                expected = np.array([True, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == idx2
                expected = np.array([False, False, False, False, False, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != idx2
                expected = np.array([True, True, True, True, True, False])
                tm.assert_numpy_array_equal(result, expected)

        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, np.nan), (didx1, pd.NaT)]:
                result = idx1 < val
                expected = np.array([False, False, False, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, True, True, True, True])
                tm.assert_numpy_array_equal(result, expected)

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, 3), (didx1, datetime(2014, 3, 1))]:
                result = idx1 < val
                expected = np.array([True, False, False, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                expected = np.array([False, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                expected = np.array([True, False, True, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                expected = np.array([False, False, True, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                expected = np.array([False, False, True, False, False, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, False, True, True, True])
                tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    def test_comparison_tzawareness_compat(self, op):
        # GH#18162
        dr = pd.date_range('2016-01-01', periods=6)
        dz = dr.tz_localize('US/Pacific')

        with pytest.raises(TypeError):
            op(dr, dz)
        with pytest.raises(TypeError):
            op(dr, list(dz))
        with pytest.raises(TypeError):
            op(dz, dr)
        with pytest.raises(TypeError):
            op(dz, list(dr))

        # Check that there isn't a problem aware-aware and naive-naive do not
        # raise
        assert (dr == dr).all()
        assert (dr == list(dr)).all()
        assert (dz == dz).all()
        assert (dz == list(dz)).all()

        # Check comparisons against scalar Timestamps
        ts = pd.Timestamp('2000-03-14 01:59')
        ts_tz = pd.Timestamp('2000-03-14 01:59', tz='Europe/Amsterdam')

        assert (dr > ts).all()
        with pytest.raises(TypeError):
            op(dr, ts_tz)

        assert (dz > ts_tz).all()
        with pytest.raises(TypeError):
            op(dz, ts)

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    def test_nat_comparison_tzawareness(self, op):
        # GH#19276
        # tzaware DatetimeIndex should not raise when compared to NaT
        dti = pd.DatetimeIndex(['2014-01-01', pd.NaT, '2014-03-01', pd.NaT,
                                '2014-05-01', '2014-07-01'])
        expected = np.array([op == operator.ne] * len(dti))
        result = op(dti, pd.NaT)
        tm.assert_numpy_array_equal(result, expected)

        result = op(dti.tz_localize('US/Pacific'), pd.NaT)
        tm.assert_numpy_array_equal(result, expected)

    def test_dti_cmp_int_raises(self):
        rng = date_range('1/1/2000', periods=10)

        # raise TypeError for now
        with pytest.raises(TypeError):
            rng < rng[3].value

    def test_dti_cmp_list(self):
        rng = date_range('1/1/2000', periods=10)

        result = rng == list(rng)
        expected = rng == rng
        tm.assert_numpy_array_equal(result, expected)


class TestDatetimeIndexArithmetic(object):

    # -------------------------------------------------------------
    # Invalid Operations

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_dti_add_sub_float(self, op, other):
        dti = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        with pytest.raises(TypeError):
            op(dti, other)

    def test_dti_add_timestamp_raises(self):
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            idx + Timestamp('2011-01-01')

    def test_dti_radd_timestamp_raises(self):
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            Timestamp('2011-01-01') + idx

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and int

    def test_dti_add_int(self, tz, one):
        # Variants of `one` for #19012
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng + one
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng - one
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_isub_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        rng -= one
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # DatetimeIndex.shift is used in integer addition

    def test_dti_shift_tzaware(self, tz):
        # GH#9903
        idx = pd.DatetimeIndex([], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        tm.assert_index_equal(idx.shift(3, freq='H'), idx)

        idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                                '2011-01-01 12:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        exp = pd.DatetimeIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                                '2011-01-01 15:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(3, freq='H'), exp)
        exp = pd.DatetimeIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                                '2011-01-01 09:00'], name='xxx', tz=tz)
        tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

    def test_dti_shift_freqs(self):
        # test shift for DatetimeIndex and non DatetimeIndex
        # GH#8083
        drange = pd.date_range('20130101', periods=5)
        result = drange.shift(1)
        expected = pd.DatetimeIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                     '2013-01-05',
                                     '2013-01-06'], freq='D')
        tm.assert_index_equal(result, expected)

        result = drange.shift(-1)
        expected = pd.DatetimeIndex(['2012-12-31', '2013-01-01', '2013-01-02',
                                     '2013-01-03', '2013-01-04'],
                                    freq='D')
        tm.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D')
        expected = pd.DatetimeIndex(['2013-01-07', '2013-01-08', '2013-01-09',
                                     '2013-01-10',
                                     '2013-01-11'], freq='D')
        tm.assert_index_equal(result, expected)

    def test_dti_shift_int(self):
        rng = date_range('1/1/2000', periods=20)

        result = rng + 5
        expected = rng.shift(5)
        tm.assert_index_equal(result, expected)

        result = rng - 5
        expected = rng.shift(-5)
        tm.assert_index_equal(result, expected)

    def test_dti_shift_no_freq(self):
        # GH#19147
        dti = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01'], freq=None)
        with pytest.raises(NullFrequencyError):
            dti.shift(2)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_shift_localized(self, tzstr):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize(tzstr)

        result = dr_tz.shift(1, '10T')
        assert result.tz == dr_tz.tz

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and timedelta-like

    def test_dti_add_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        result = rng + delta
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        rng += delta
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        result = rng - delta
        tm.assert_index_equal(result, expected)

    def test_dti_isub_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        rng -= delta
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and TimedeltaIndex/array
    def test_dti_add_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # add with TimdeltaIndex
        result = dti + tdi
        tm.assert_index_equal(result, expected)

        result = tdi + dti
        tm.assert_index_equal(result, expected)

        # add with timedelta64 array
        result = dti + tdi.values
        tm.assert_index_equal(result, expected)

        result = tdi.values + dti
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # iadd with TimdeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

        # iadd with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi.values
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # sub with TimedeltaIndex
        result = dti - tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract .*TimedeltaIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi - dti

        # sub with timedelta64 array
        result = dti - tdi.values
        tm.assert_index_equal(result, expected)

        msg = 'cannot perform __neg__ with this index type:'
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values - dti

    def test_dti_isub_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # isub with TimedeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract .*TimedeltaIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi -= dti

        # isub with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi.values
        tm.assert_index_equal(result, expected)

        msg = '|'.join(['cannot perform __neg__ with this index type:',
                        'ufunc subtract cannot use operands with types'])
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values -= dti

    # -------------------------------------------------------------
    # Binary Operations DatetimeIndex and datetime-like
    # TODO: A couple other tests belong in this section.  Move them in
    # A PR where there isn't already a giant diff.

    def test_add_datetimelike_and_dti(self, addend):
        # GH#9631
        dti = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = 'cannot add DatetimeIndex and {0}'.format(
            type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti

    def test_add_datetimelike_and_dti_tz(self, addend):
        # GH#9631
        dti_tz = DatetimeIndex(['2011-01-01',
                                '2011-01-02']).tz_localize('US/Eastern')
        msg = 'cannot add DatetimeIndex and {0}'.format(
            type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti_tz + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti_tz

    # -------------------------------------------------------------
    # __add__/__sub__ with ndarray[datetime64] and ndarray[timedelta64]

    def test_dti_add_dt64_array_raises(self, tz):
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        dtarr = dti.values

        with pytest.raises(TypeError):
            dti + dtarr
        with pytest.raises(TypeError):
            dtarr + dti

    def test_dti_sub_dt64_array_naive(self):
        dti = pd.date_range('2016-01-01', periods=3, tz=None)
        dtarr = dti.values

        expected = dti - dti
        result = dti - dtarr
        tm.assert_index_equal(result, expected)
        result = dtarr - dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_dt64_array_aware_raises(self, tz):
        if tz is None:
            return
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        dtarr = dti.values

        with pytest.raises(TypeError):
            dti - dtarr
        with pytest.raises(TypeError):
            dtarr - dti

    def test_dti_add_td64_array(self, tz):
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        expected = dti + tdi
        result = dti + tdarr
        tm.assert_index_equal(result, expected)
        result = tdarr + dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_td64_array(self, tz):
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        expected = dti - tdi
        result = dti - tdarr
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            tdarr - dti

    # -------------------------------------------------------------

    def test_sub_dti_dti(self):
        # previously performed setop (deprecated in 0.16.0), now changed to
        # return subtraction -> TimeDeltaIndex (GH ...)

        dti = date_range('20130101', periods=3)
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')
        dti_tz2 = date_range('20130101', periods=3).tz_localize('UTC')
        expected = TimedeltaIndex([0, 0, 0])

        result = dti - dti
        tm.assert_index_equal(result, expected)

        result = dti_tz - dti_tz
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            dti_tz - dti

        with pytest.raises(TypeError):
            dti - dti_tz

        with pytest.raises(TypeError):
            dti_tz - dti_tz2

        # isub
        dti -= dti
        tm.assert_index_equal(dti, expected)

        # different length raises ValueError
        dti1 = date_range('20130101', periods=3)
        dti2 = date_range('20130101', periods=4)
        with pytest.raises(ValueError):
            dti1 - dti2

        # NaN propagation
        dti1 = DatetimeIndex(['2012-01-01', np.nan, '2012-01-03'])
        dti2 = DatetimeIndex(['2012-01-02', '2012-01-03', np.nan])
        expected = TimedeltaIndex(['1 days', np.nan, np.nan])
        result = dti2 - dti1
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('freq', [None, 'D'])
    def test_sub_period(self, freq):
        # GH#13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq=freq)

        with pytest.raises(TypeError):
            idx - p

        with pytest.raises(TypeError):
            p - idx

    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '3D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == 'D'

    def test_datetimeindex_sub_timestamp_overflow(self):
        dtimax = pd.to_datetime(['now', pd.Timestamp.max])
        dtimin = pd.to_datetime(['now', pd.Timestamp.min])

        tsneg = Timestamp('1950-01-01')
        ts_neg_variants = [tsneg,
                           tsneg.to_pydatetime(),
                           tsneg.to_datetime64().astype('datetime64[ns]'),
                           tsneg.to_datetime64().astype('datetime64[D]')]

        tspos = Timestamp('1980-01-01')
        ts_pos_variants = [tspos,
                           tspos.to_pydatetime(),
                           tspos.to_datetime64().astype('datetime64[ns]'),
                           tspos.to_datetime64().astype('datetime64[D]')]

        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        expected = pd.Timestamp.max.value - tspos.value
        for variant in ts_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        expected = pd.Timestamp.min.value - tsneg.value
        for variant in ts_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant

    @pytest.mark.parametrize('names', [('foo', None, None),
                                       ('baz', 'bar', None),
                                       ('bar', 'bar', 'bar')])
    @pytest.mark.parametrize('tz', [None, 'America/Chicago'])
    def test_dti_add_series(self, tz, names):
        # GH#13905
        index = DatetimeIndex(['2016-06-28 05:30', '2016-06-28 05:31'],
                              tz=tz, name=names[0])
        ser = Series([Timedelta(seconds=5)] * 2,
                     index=index, name=names[1])
        expected = Series(index + Timedelta(seconds=5),
                          index=index, name=names[2])

        # passing name arg isn't enough when names[2] is None
        expected.name = names[2]
        assert expected.dtype == index.dtype
        result = ser + index
        tm.assert_series_equal(result, expected)
        result2 = index + ser
        tm.assert_series_equal(result2, expected)

        expected = index + Timedelta(seconds=5)
        result3 = ser.values + index
        tm.assert_index_equal(result3, expected)
        result4 = index + ser.values
        tm.assert_index_equal(result4, expected)

    def test_dti_add_offset_array(self, tz):
        # GH#18849
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_add_offset_index(self, tz, names):
        # GH#18849, GH#19744
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    def test_dti_sub_offset_array(self, tz):
        # GH#18824
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_sub_offset_index(self, tz, names):
        # GH#18824, GH#19744
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_with_offset_series(self, tz, names):
        # GH#18849
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = Series([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                       name=names[1])

        expected_add = Series([dti[n] + other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        tm.assert_series_equal(res, expected_add)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_series_equal(res2, expected_add)

        expected_sub = Series([dti[n] - other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res3 = dti - other
        tm.assert_series_equal(res3, expected_sub)

    def test_dti_add_offset_tzaware(self, tz_aware_fixture):
        timezone = tz_aware_fixture
        if timezone == 'US/Pacific':
            dates = date_range('2012-11-01', periods=3, tz=timezone)
            offset = dates + pd.offsets.Hour(5)
            assert dates[0] + pd.offsets.Hour(5) == offset[0]

        dates = date_range('2010-11-01 00:00',
                           periods=3, tz=timezone, freq='H')
        expected = DatetimeIndex(['2010-11-01 05:00', '2010-11-01 06:00',
                                  '2010-11-01 07:00'], freq='H', tz=timezone)

        offset = dates + pd.offsets.Hour(5)
        tm.assert_index_equal(offset, expected)
        offset = dates + np.timedelta64(5, 'h')
        tm.assert_index_equal(offset, expected)
        offset = dates + timedelta(hours=5)
        tm.assert_index_equal(offset, expected)


@pytest.mark.parametrize('klass,assert_func', [
    (Series, tm.assert_series_equal),
    (DatetimeIndex, tm.assert_index_equal)])
def test_dt64_with_offset_array(klass, assert_func):
    # GH#10699
    # array of offsets
    box = Series if klass is Series else pd.Index
    with tm.assert_produces_warning(PerformanceWarning):
        s = klass([Timestamp('2000-1-1'), Timestamp('2000-2-1')])
        result = s + box([pd.offsets.DateOffset(years=1),
                          pd.offsets.MonthEnd()])
        exp = klass([Timestamp('2001-1-1'), Timestamp('2000-2-29')])
        assert_func(result, exp)

        # same offset
        result = s + box([pd.offsets.DateOffset(years=1),
                          pd.offsets.DateOffset(years=1)])
        exp = klass([Timestamp('2001-1-1'), Timestamp('2001-2-1')])
        assert_func(result, exp)


@pytest.mark.parametrize('klass,assert_func', [
    (Series, tm.assert_series_equal),
    (DatetimeIndex, tm.assert_index_equal)])
def test_dt64_with_DateOffsets_relativedelta(klass, assert_func):
    # GH#10699
    vec = klass([Timestamp('2000-01-05 00:15:00'),
                 Timestamp('2000-01-31 00:23:00'),
                 Timestamp('2000-01-01'),
                 Timestamp('2000-03-31'),
                 Timestamp('2000-02-29'),
                 Timestamp('2000-12-31'),
                 Timestamp('2000-05-15'),
                 Timestamp('2001-06-15')])

    # DateOffset relativedelta fastpath
    relative_kwargs = [('years', 2), ('months', 5), ('days', 3),
                       ('hours', 5), ('minutes', 10), ('seconds', 2),
                       ('microseconds', 5)]
    for i, kwd in enumerate(relative_kwargs):
        op = pd.DateOffset(**dict([kwd]))
        assert_func(klass([x + op for x in vec]), vec + op)
        assert_func(klass([x - op for x in vec]), vec - op)
        op = pd.DateOffset(**dict(relative_kwargs[:i + 1]))
        assert_func(klass([x + op for x in vec]), vec + op)
        assert_func(klass([x - op for x in vec]), vec - op)


@pytest.mark.parametrize('cls_and_kwargs', [
    'YearBegin', ('YearBegin', {'month': 5}),
    'YearEnd', ('YearEnd', {'month': 5}),
    'MonthBegin', 'MonthEnd',
    'SemiMonthEnd', 'SemiMonthBegin',
    'Week', ('Week', {'weekday': 3}),
    'BusinessDay', 'BDay', 'QuarterEnd', 'QuarterBegin',
    'CustomBusinessDay', 'CDay', 'CBMonthEnd',
    'CBMonthBegin', 'BMonthBegin', 'BMonthEnd',
    'BusinessHour', 'BYearBegin', 'BYearEnd',
    'BQuarterBegin', ('LastWeekOfMonth', {'weekday': 2}),
    ('FY5253Quarter', {'qtr_with_extra_week': 1,
                       'startingMonth': 1,
                       'weekday': 2,
                       'variation': 'nearest'}),
    ('FY5253', {'weekday': 0, 'startingMonth': 2, 'variation': 'nearest'}),
    ('WeekOfMonth', {'weekday': 2, 'week': 2}),
    'Easter', ('DateOffset', {'day': 4}),
    ('DateOffset', {'month': 5})])
@pytest.mark.parametrize('normalize', [True, False])
@pytest.mark.parametrize('klass,assert_func', [
    (Series, tm.assert_series_equal),
    (DatetimeIndex, tm.assert_index_equal)])
def test_dt64_with_DateOffsets(klass, assert_func, normalize, cls_and_kwargs):
    # GH#10699
    # assert these are equal on a piecewise basis
    vec = klass([Timestamp('2000-01-05 00:15:00'),
                 Timestamp('2000-01-31 00:23:00'),
                 Timestamp('2000-01-01'),
                 Timestamp('2000-03-31'),
                 Timestamp('2000-02-29'),
                 Timestamp('2000-12-31'),
                 Timestamp('2000-05-15'),
                 Timestamp('2001-06-15')])

    if isinstance(cls_and_kwargs, tuple):
        # If cls_name param is a tuple, then 2nd entry is kwargs for
        # the offset constructor
        cls_name, kwargs = cls_and_kwargs
    else:
        cls_name = cls_and_kwargs
        kwargs = {}

    offset_cls = getattr(pd.offsets, cls_name)

    with warnings.catch_warnings(record=True):
        for n in [0, 5]:
            if (cls_name in ['WeekOfMonth', 'LastWeekOfMonth',
                             'FY5253Quarter', 'FY5253'] and n == 0):
                # passing n = 0 is invalid for these offset classes
                continue

            offset = offset_cls(n, normalize=normalize, **kwargs)
            assert_func(klass([x + offset for x in vec]), vec + offset)
            assert_func(klass([x - offset for x in vec]), vec - offset)
            assert_func(klass([offset + x for x in vec]), offset + vec)


# GH 10699
@pytest.mark.parametrize('klass,assert_func', zip([Series, DatetimeIndex],
                                                  [tm.assert_series_equal,
                                                   tm.assert_index_equal]))
def test_datetime64_with_DateOffset(klass, assert_func):
    s = klass(date_range('2000-01-01', '2000-01-31'), name='a')
    result = s + pd.DateOffset(years=1)
    result2 = pd.DateOffset(years=1) + s
    exp = klass(date_range('2001-01-01', '2001-01-31'), name='a')
    assert_func(result, exp)
    assert_func(result2, exp)

    result = s - pd.DateOffset(years=1)
    exp = klass(date_range('1999-01-01', '1999-01-31'), name='a')
    assert_func(result, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.Day()
    result2 = pd.offsets.Day() + s
    exp = klass([Timestamp('2000-01-16 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-16', tz='US/Central')], name='a')
    assert_func(result, exp)
    assert_func(result2, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.MonthEnd()
    result2 = pd.offsets.MonthEnd() + s
    exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-29', tz='US/Central')], name='a')
    assert_func(result, exp)
    assert_func(result2, exp)


@pytest.mark.parametrize('years', [-1, 0, 1])
@pytest.mark.parametrize('months', [-2, 0, 2])
def test_shift_months(years, months):
    s = DatetimeIndex([Timestamp('2000-01-05 00:15:00'),
                       Timestamp('2000-01-31 00:23:00'),
                       Timestamp('2000-01-01'),
                       Timestamp('2000-02-29'),
                       Timestamp('2000-12-31')])
    actual = DatetimeIndex(shift_months(s.asi8, years * 12 + months))

    raw = [x + pd.offsets.DateOffset(years=years, months=months)
           for x in s]
    expected = DatetimeIndex(raw)
    tm.assert_index_equal(actual, expected)
