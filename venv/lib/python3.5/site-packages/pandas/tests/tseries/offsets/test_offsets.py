import os
from distutils.version import LooseVersion
from datetime import date, datetime, timedelta

import pytest
from pandas.compat import range
from pandas import compat

import numpy as np

from pandas.compat.numpy import np_datetime64_compat

from pandas.core.series import Series
from pandas._libs.tslibs.frequencies import (get_freq_code, get_freq_str,
                                             _INVALID_FREQ_ERROR)
from pandas.tseries.frequencies import _offset_map, get_offset
from pandas.core.indexes.datetimes import (
    _to_m8, DatetimeIndex, _daterange_cache)
import pandas._libs.tslibs.offsets as liboffsets
from pandas._libs.tslibs.offsets import WeekDay, CacheableOffset
from pandas.tseries.offsets import (BDay, CDay, BQuarterEnd, BMonthEnd,
                                    BusinessHour, WeekOfMonth, CBMonthEnd,
                                    CustomBusinessHour,
                                    CBMonthBegin, BYearEnd, MonthEnd,
                                    MonthBegin, SemiMonthBegin, SemiMonthEnd,
                                    BYearBegin, QuarterBegin, BQuarterBegin,
                                    BMonthBegin, DateOffset, Week, YearBegin,
                                    YearEnd, Day,
                                    QuarterEnd, BusinessMonthEnd, FY5253,
                                    Nano, Easter, FY5253Quarter,
                                    LastWeekOfMonth)
from pandas.core.tools.datetimes import format, ole2datetime
import pandas.tseries.offsets as offsets
from pandas.io.pickle import read_pickle
from pandas._libs.tslibs import timezones
from pandas._libs.tslib import normalize_date, NaT, Timestamp
import pandas._libs.tslib as tslib
import pandas.util.testing as tm
from pandas.tseries.holiday import USFederalHolidayCalendar

from .common import assert_offset_equal, assert_onOffset


def test_monthrange():
    import calendar
    for y in range(2000, 2013):
        for m in range(1, 13):
            assert tslib.monthrange(y, m) == calendar.monthrange(y, m)

####
# Misc function tests
####


def test_format():
    actual = format(datetime(2008, 1, 15))
    assert actual == '20080115'


def test_ole2datetime():
    actual = ole2datetime(60000)
    assert actual == datetime(2064, 4, 8)

    with pytest.raises(ValueError):
        ole2datetime(60)


def test_normalize_date():
    actual = normalize_date(datetime(2007, 10, 1, 1, 12, 5, 10))
    assert actual == datetime(2007, 10, 1)


def test_to_m8():
    valb = datetime(2007, 10, 1)
    valu = _to_m8(valb)
    assert isinstance(valu, np.datetime64)
    # assert valu == np.datetime64(datetime(2007,10,1))

    # def test_datetime64_box():
    #    valu = np.datetime64(datetime(2007,10,1))
    #    valb = _dt_box(valu)
    #    assert type(valb) == datetime
    #    assert valb == datetime(2007,10,1)

    #####
    # DateOffset Tests
    #####


class Base(object):
    _offset = None

    timezones = [None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
                 'dateutil/Asia/Tokyo', 'dateutil/US/Pacific']

    def _get_offset(self, klass, value=1, normalize=False):
        # create instance from offset class
        if klass is FY5253:
            klass = klass(n=value, startingMonth=1, weekday=1,
                          variation='last', normalize=normalize)
        elif klass is FY5253Quarter:
            klass = klass(n=value, startingMonth=1, weekday=1,
                          qtr_with_extra_week=1, variation='last',
                          normalize=normalize)
        elif klass is LastWeekOfMonth:
            klass = klass(n=value, weekday=5, normalize=normalize)
        elif klass is WeekOfMonth:
            klass = klass(n=value, week=1, weekday=5, normalize=normalize)
        elif klass is Week:
            klass = klass(n=value, weekday=5, normalize=normalize)
        elif klass is DateOffset:
            klass = klass(days=value, normalize=normalize)
        else:
            try:
                klass = klass(value, normalize=normalize)
            except Exception:
                klass = klass(normalize=normalize)
        return klass

    def test_apply_out_of_range(self, tz):
        if self._offset is None:
            return

        # try to create an out-of-bounds result timestamp; if we can't create
        # the offset skip
        try:
            if self._offset in (BusinessHour, CustomBusinessHour):
                # Using 10000 in BusinessHour fails in tz check because of DST
                # difference
                offset = self._get_offset(self._offset, value=100000)
            else:
                offset = self._get_offset(self._offset, value=10000)

            result = Timestamp('20080101') + offset
            assert isinstance(result, datetime)
            assert result.tzinfo is None

            # Check tz is preserved
            t = Timestamp('20080101', tz=tz)
            result = t + offset
            assert isinstance(result, datetime)
            assert t.tzinfo == result.tzinfo

        except tslib.OutOfBoundsDatetime:
            raise
        except (ValueError, KeyError):
            # we are creating an invalid offset
            # so ignore
            pass


class TestCommon(Base):
    # exected value created by Base._get_offset
    # are applied to 2011/01/01 09:00 (Saturday)
    # used for .apply and .rollforward
    expecteds = {'Day': Timestamp('2011-01-02 09:00:00'),
                 'DateOffset': Timestamp('2011-01-02 09:00:00'),
                 'BusinessDay': Timestamp('2011-01-03 09:00:00'),
                 'CustomBusinessDay': Timestamp('2011-01-03 09:00:00'),
                 'CustomBusinessMonthEnd': Timestamp('2011-01-31 09:00:00'),
                 'CustomBusinessMonthBegin': Timestamp('2011-01-03 09:00:00'),
                 'MonthBegin': Timestamp('2011-02-01 09:00:00'),
                 'BusinessMonthBegin': Timestamp('2011-01-03 09:00:00'),
                 'MonthEnd': Timestamp('2011-01-31 09:00:00'),
                 'SemiMonthEnd': Timestamp('2011-01-15 09:00:00'),
                 'SemiMonthBegin': Timestamp('2011-01-15 09:00:00'),
                 'BusinessMonthEnd': Timestamp('2011-01-31 09:00:00'),
                 'YearBegin': Timestamp('2012-01-01 09:00:00'),
                 'BYearBegin': Timestamp('2011-01-03 09:00:00'),
                 'YearEnd': Timestamp('2011-12-31 09:00:00'),
                 'BYearEnd': Timestamp('2011-12-30 09:00:00'),
                 'QuarterBegin': Timestamp('2011-03-01 09:00:00'),
                 'BQuarterBegin': Timestamp('2011-03-01 09:00:00'),
                 'QuarterEnd': Timestamp('2011-03-31 09:00:00'),
                 'BQuarterEnd': Timestamp('2011-03-31 09:00:00'),
                 'BusinessHour': Timestamp('2011-01-03 10:00:00'),
                 'CustomBusinessHour': Timestamp('2011-01-03 10:00:00'),
                 'WeekOfMonth': Timestamp('2011-01-08 09:00:00'),
                 'LastWeekOfMonth': Timestamp('2011-01-29 09:00:00'),
                 'FY5253Quarter': Timestamp('2011-01-25 09:00:00'),
                 'FY5253': Timestamp('2011-01-25 09:00:00'),
                 'Week': Timestamp('2011-01-08 09:00:00'),
                 'Easter': Timestamp('2011-04-24 09:00:00'),
                 'Hour': Timestamp('2011-01-01 10:00:00'),
                 'Minute': Timestamp('2011-01-01 09:01:00'),
                 'Second': Timestamp('2011-01-01 09:00:01'),
                 'Milli': Timestamp('2011-01-01 09:00:00.001000'),
                 'Micro': Timestamp('2011-01-01 09:00:00.000001'),
                 'Nano': Timestamp(np_datetime64_compat(
                                   '2011-01-01T09:00:00.000000001Z'))}

    def test_return_type(self, offset_types):
        offset = self._get_offset(offset_types)

        # make sure that we are returning a Timestamp
        result = Timestamp('20080101') + offset
        assert isinstance(result, Timestamp)

        # make sure that we are returning NaT
        assert NaT + offset is NaT
        assert offset + NaT is NaT

        assert NaT - offset is NaT
        assert (-offset).apply(NaT) is NaT

    def test_offset_n(self, offset_types):
        offset = self._get_offset(offset_types)
        assert offset.n == 1

        neg_offset = offset * -1
        assert neg_offset.n == -1

        mul_offset = offset * 3
        assert mul_offset.n == 3

    def test_offset_freqstr(self, offset_types):
        offset = self._get_offset(offset_types)

        freqstr = offset.freqstr
        if freqstr not in ('<Easter>',
                           "<DateOffset: days=1>",
                           'LWOM-SAT', ):
            code = get_offset(freqstr)
            assert offset.rule_code == code

    def _check_offsetfunc_works(self, offset, funcname, dt, expected,
                                normalize=False):
        offset_s = self._get_offset(offset, normalize=normalize)
        func = getattr(offset_s, funcname)

        result = func(dt)
        assert isinstance(result, Timestamp)
        assert result == expected

        result = func(Timestamp(dt))
        assert isinstance(result, Timestamp)
        assert result == expected

        # see gh-14101
        exp_warning = None
        ts = Timestamp(dt) + Nano(5)

        if (offset_s.__class__.__name__ == 'DateOffset' and
                (funcname == 'apply' or normalize) and
                ts.nanosecond > 0):
            exp_warning = UserWarning

        # test nanosecond is preserved
        with tm.assert_produces_warning(exp_warning,
                                        check_stacklevel=False):
            result = func(ts)
        assert isinstance(result, Timestamp)
        if normalize is False:
            assert result == expected + Nano(5)
        else:
            assert result == expected

        if isinstance(dt, np.datetime64):
            # test tz when input is datetime or Timestamp
            return

        for tz in self.timezones:
            expected_localize = expected.tz_localize(tz)
            tz_obj = timezones.maybe_get_tz(tz)
            dt_tz = tslib._localize_pydatetime(dt, tz_obj)

            result = func(dt_tz)
            assert isinstance(result, Timestamp)
            assert result == expected_localize

            result = func(Timestamp(dt, tz=tz))
            assert isinstance(result, Timestamp)
            assert result == expected_localize

            # see gh-14101
            exp_warning = None
            ts = Timestamp(dt, tz=tz) + Nano(5)

            if (offset_s.__class__.__name__ == 'DateOffset' and
                    (funcname == 'apply' or normalize) and
                    ts.nanosecond > 0):
                exp_warning = UserWarning

            # test nanosecond is preserved
            with tm.assert_produces_warning(exp_warning,
                                            check_stacklevel=False):
                result = func(ts)
            assert isinstance(result, Timestamp)
            if normalize is False:
                assert result == expected_localize + Nano(5)
            else:
                assert result == expected_localize

    def test_apply(self, offset_types):
        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat('2011-01-01 09:00Z')

        for dt in [sdt, ndt]:
            expected = self.expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, 'apply', dt, expected)

            expected = Timestamp(expected.date())
            self._check_offsetfunc_works(offset_types, 'apply', dt, expected,
                                         normalize=True)

    def test_rollforward(self, offset_types):
        expecteds = self.expecteds.copy()

        # result will not be changed if the target is on the offset
        no_changes = ['Day', 'MonthBegin', 'SemiMonthBegin', 'YearBegin',
                      'Week', 'Hour', 'Minute', 'Second', 'Milli', 'Micro',
                      'Nano', 'DateOffset']
        for n in no_changes:
            expecteds[n] = Timestamp('2011/01/01 09:00')

        expecteds['BusinessHour'] = Timestamp('2011-01-03 09:00:00')
        expecteds['CustomBusinessHour'] = Timestamp('2011-01-03 09:00:00')

        # but be changed when normalize=True
        norm_expected = expecteds.copy()
        for k in norm_expected:
            norm_expected[k] = Timestamp(norm_expected[k].date())

        normalized = {'Day': Timestamp('2011-01-02 00:00:00'),
                      'DateOffset': Timestamp('2011-01-02 00:00:00'),
                      'MonthBegin': Timestamp('2011-02-01 00:00:00'),
                      'SemiMonthBegin': Timestamp('2011-01-15 00:00:00'),
                      'YearBegin': Timestamp('2012-01-01 00:00:00'),
                      'Week': Timestamp('2011-01-08 00:00:00'),
                      'Hour': Timestamp('2011-01-01 00:00:00'),
                      'Minute': Timestamp('2011-01-01 00:00:00'),
                      'Second': Timestamp('2011-01-01 00:00:00'),
                      'Milli': Timestamp('2011-01-01 00:00:00'),
                      'Micro': Timestamp('2011-01-01 00:00:00')}
        norm_expected.update(normalized)

        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat('2011-01-01 09:00Z')

        for dt in [sdt, ndt]:
            expected = expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, 'rollforward', dt,
                                         expected)
            expected = norm_expected[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, 'rollforward', dt,
                                         expected, normalize=True)

    def test_rollback(self, offset_types):
        expecteds = {'BusinessDay': Timestamp('2010-12-31 09:00:00'),
                     'CustomBusinessDay': Timestamp('2010-12-31 09:00:00'),
                     'CustomBusinessMonthEnd':
                     Timestamp('2010-12-31 09:00:00'),
                     'CustomBusinessMonthBegin':
                     Timestamp('2010-12-01 09:00:00'),
                     'BusinessMonthBegin': Timestamp('2010-12-01 09:00:00'),
                     'MonthEnd': Timestamp('2010-12-31 09:00:00'),
                     'SemiMonthEnd': Timestamp('2010-12-31 09:00:00'),
                     'BusinessMonthEnd': Timestamp('2010-12-31 09:00:00'),
                     'BYearBegin': Timestamp('2010-01-01 09:00:00'),
                     'YearEnd': Timestamp('2010-12-31 09:00:00'),
                     'BYearEnd': Timestamp('2010-12-31 09:00:00'),
                     'QuarterBegin': Timestamp('2010-12-01 09:00:00'),
                     'BQuarterBegin': Timestamp('2010-12-01 09:00:00'),
                     'QuarterEnd': Timestamp('2010-12-31 09:00:00'),
                     'BQuarterEnd': Timestamp('2010-12-31 09:00:00'),
                     'BusinessHour': Timestamp('2010-12-31 17:00:00'),
                     'CustomBusinessHour': Timestamp('2010-12-31 17:00:00'),
                     'WeekOfMonth': Timestamp('2010-12-11 09:00:00'),
                     'LastWeekOfMonth': Timestamp('2010-12-25 09:00:00'),
                     'FY5253Quarter': Timestamp('2010-10-26 09:00:00'),
                     'FY5253': Timestamp('2010-01-26 09:00:00'),
                     'Easter': Timestamp('2010-04-04 09:00:00')}

        # result will not be changed if the target is on the offset
        for n in ['Day', 'MonthBegin', 'SemiMonthBegin', 'YearBegin', 'Week',
                  'Hour', 'Minute', 'Second', 'Milli', 'Micro', 'Nano',
                  'DateOffset']:
            expecteds[n] = Timestamp('2011/01/01 09:00')

        # but be changed when normalize=True
        norm_expected = expecteds.copy()
        for k in norm_expected:
            norm_expected[k] = Timestamp(norm_expected[k].date())

        normalized = {'Day': Timestamp('2010-12-31 00:00:00'),
                      'DateOffset': Timestamp('2010-12-31 00:00:00'),
                      'MonthBegin': Timestamp('2010-12-01 00:00:00'),
                      'SemiMonthBegin': Timestamp('2010-12-15 00:00:00'),
                      'YearBegin': Timestamp('2010-01-01 00:00:00'),
                      'Week': Timestamp('2010-12-25 00:00:00'),
                      'Hour': Timestamp('2011-01-01 00:00:00'),
                      'Minute': Timestamp('2011-01-01 00:00:00'),
                      'Second': Timestamp('2011-01-01 00:00:00'),
                      'Milli': Timestamp('2011-01-01 00:00:00'),
                      'Micro': Timestamp('2011-01-01 00:00:00')}
        norm_expected.update(normalized)

        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat('2011-01-01 09:00Z')

        for dt in [sdt, ndt]:
            expected = expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, 'rollback', dt,
                                         expected)

            expected = norm_expected[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, 'rollback', dt,
                                         expected, normalize=True)

    def test_onOffset(self, offset_types):
        dt = self.expecteds[offset_types.__name__]
        offset_s = self._get_offset(offset_types)
        assert offset_s.onOffset(dt)

        # when normalize=True, onOffset checks time is 00:00:00
        offset_n = self._get_offset(offset_types, normalize=True)
        assert not offset_n.onOffset(dt)

        if offset_types in (BusinessHour, CustomBusinessHour):
            # In default BusinessHour (9:00-17:00), normalized time
            # cannot be in business hour range
            return
        date = datetime(dt.year, dt.month, dt.day)
        assert offset_n.onOffset(date)

    def test_add(self, offset_types, tz):
        dt = datetime(2011, 1, 1, 9, 0)

        offset_s = self._get_offset(offset_types)
        expected = self.expecteds[offset_types.__name__]

        result_dt = dt + offset_s
        result_ts = Timestamp(dt) + offset_s
        for result in [result_dt, result_ts]:
            assert isinstance(result, Timestamp)
            assert result == expected

        expected_localize = expected.tz_localize(tz)
        result = Timestamp(dt, tz=tz) + offset_s
        assert isinstance(result, Timestamp)
        assert result == expected_localize

        # normalize=True
        offset_s = self._get_offset(offset_types, normalize=True)
        expected = Timestamp(expected.date())

        result_dt = dt + offset_s
        result_ts = Timestamp(dt) + offset_s
        for result in [result_dt, result_ts]:
            assert isinstance(result, Timestamp)
            assert result == expected

        expected_localize = expected.tz_localize(tz)
        result = Timestamp(dt, tz=tz) + offset_s
        assert isinstance(result, Timestamp)
        assert result == expected_localize

    def test_pickle_v0_15_2(self):
        offsets = {'DateOffset': DateOffset(years=1),
                   'MonthBegin': MonthBegin(1),
                   'Day': Day(1),
                   'YearBegin': YearBegin(1),
                   'Week': Week(1)}
        pickle_path = os.path.join(tm.get_data_path(),
                                   'dateoffset_0_15_2.pickle')
        # This code was executed once on v0.15.2 to generate the pickle:
        # with open(pickle_path, 'wb') as f: pickle.dump(offsets, f)
        #
        tm.assert_dict_equal(offsets, read_pickle(pickle_path))


class TestDateOffset(Base):

    def setup_method(self, method):
        self.d = Timestamp(datetime(2008, 1, 2))
        _offset_map.clear()

    def test_repr(self):
        repr(DateOffset())
        repr(DateOffset(2))
        repr(2 * DateOffset())
        repr(2 * DateOffset(months=2))

    def test_mul(self):
        assert DateOffset(2) == 2 * DateOffset(1)
        assert DateOffset(2) == DateOffset(1) * 2

    def test_constructor(self):

        assert ((self.d + DateOffset(months=2)) == datetime(2008, 3, 2))
        assert ((self.d - DateOffset(months=2)) == datetime(2007, 11, 2))

        assert ((self.d + DateOffset(2)) == datetime(2008, 1, 4))

        assert not DateOffset(2).isAnchored()
        assert DateOffset(1).isAnchored()

        d = datetime(2008, 1, 31)
        assert ((d + DateOffset(months=1)) == datetime(2008, 2, 29))

    def test_copy(self):
        assert (DateOffset(months=2).copy() == DateOffset(months=2))

    def test_eq(self):
        offset1 = DateOffset(days=1)
        offset2 = DateOffset(days=365)

        assert offset1 != offset2


class TestBusinessDay(Base):
    _offset = BDay

    def setup_method(self, method):
        self.d = datetime(2008, 1, 1)

        self.offset = BDay()
        self.offset2 = BDay(2)

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = BDay()
        offset2 = BDay()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset) == '<BusinessDay>'
        assert repr(self.offset2) == '<2 * BusinessDays>'

        expected = '<BusinessDay: offset=datetime.timedelta(1)>'
        assert repr(self.offset + timedelta(1)) == expected

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    def testEQ(self):
        assert self.offset2 == self.offset2

    def test_mul(self):
        pass

    def test_hash(self):
        assert hash(self.offset2) == hash(self.offset2)

    def testCall(self):
        assert self.offset2(self.d) == datetime(2008, 1, 3)

    def testRAdd(self):
        assert self.d + self.offset2 == self.offset2 + self.d

    def testSub(self):
        off = self.offset2
        pytest.raises(Exception, off.__sub__, self.d)
        assert 2 * off - off == off

        assert self.d - self.offset2 == self.d + BDay(-2)

    def testRSub(self):
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def testMult1(self):
        assert self.d + 10 * self.offset == self.d + BDay(10)

    def testMult2(self):
        assert self.d + (-5 * BDay(-10)) == self.d + BDay(50)

    def testRollback1(self):
        assert BDay(10).rollback(self.d) == self.d

    def testRollback2(self):
        assert (BDay(10).rollback(datetime(2008, 1, 5)) ==
                datetime(2008, 1, 4))

    def testRollforward1(self):
        assert BDay(10).rollforward(self.d) == self.d

    def testRollforward2(self):
        assert (BDay(10).rollforward(datetime(2008, 1, 5)) ==
                datetime(2008, 1, 7))

    def test_roll_date_object(self):
        offset = BDay()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 14)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 17)

        offset = offsets.Day()
        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 15)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 15)

    def test_onOffset(self):
        tests = [(BDay(), datetime(2008, 1, 1), True),
                 (BDay(), datetime(2008, 1, 5), False)]

        for offset, d, expected in tests:
            assert_onOffset(offset, d, expected)

    apply_cases = []
    apply_cases.append((BDay(), {
        datetime(2008, 1, 1): datetime(2008, 1, 2),
        datetime(2008, 1, 4): datetime(2008, 1, 7),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 8)}))

    apply_cases.append((2 * BDay(), {
        datetime(2008, 1, 1): datetime(2008, 1, 3),
        datetime(2008, 1, 4): datetime(2008, 1, 8),
        datetime(2008, 1, 5): datetime(2008, 1, 8),
        datetime(2008, 1, 6): datetime(2008, 1, 8),
        datetime(2008, 1, 7): datetime(2008, 1, 9)}))

    apply_cases.append((-BDay(), {
        datetime(2008, 1, 1): datetime(2007, 12, 31),
        datetime(2008, 1, 4): datetime(2008, 1, 3),
        datetime(2008, 1, 5): datetime(2008, 1, 4),
        datetime(2008, 1, 6): datetime(2008, 1, 4),
        datetime(2008, 1, 7): datetime(2008, 1, 4),
        datetime(2008, 1, 8): datetime(2008, 1, 7)}))

    apply_cases.append((-2 * BDay(), {
        datetime(2008, 1, 1): datetime(2007, 12, 28),
        datetime(2008, 1, 4): datetime(2008, 1, 2),
        datetime(2008, 1, 5): datetime(2008, 1, 3),
        datetime(2008, 1, 6): datetime(2008, 1, 3),
        datetime(2008, 1, 7): datetime(2008, 1, 3),
        datetime(2008, 1, 8): datetime(2008, 1, 4),
        datetime(2008, 1, 9): datetime(2008, 1, 7)}))

    apply_cases.append((BDay(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2008, 1, 4): datetime(2008, 1, 4),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 7)}))

    @pytest.mark.parametrize('case', apply_cases)
    def test_apply(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + BDay(10)
        assert result == datetime(2012, 11, 6)

        result = dt + BDay(100) - BDay(100)
        assert result == dt

        off = BDay() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 12, 23)
        assert rs == xp

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2011, 12, 26)
        assert rs == xp

        off = BDay() * 10
        rs = datetime(2014, 1, 5) + off  # see #5890
        xp = datetime(2014, 1, 17)
        assert rs == xp

    def test_apply_corner(self):
        pytest.raises(TypeError, BDay().apply, BMonthEnd())

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BDay()
        offset2 = BDay()
        assert not offset1 != offset2


class TestBusinessHour(Base):
    _offset = BusinessHour

    def setup_method(self, method):
        self.d = datetime(2014, 7, 1, 10, 00)

        self.offset1 = BusinessHour()
        self.offset2 = BusinessHour(n=3)

        self.offset3 = BusinessHour(n=-1)
        self.offset4 = BusinessHour(n=-4)

        from datetime import time as dt_time
        self.offset5 = BusinessHour(start=dt_time(11, 0), end=dt_time(14, 30))
        self.offset6 = BusinessHour(start='20:00', end='05:00')
        self.offset7 = BusinessHour(n=-2, start=dt_time(21, 30),
                                    end=dt_time(6, 30))

    def test_constructor_errors(self):
        from datetime import time as dt_time
        with pytest.raises(ValueError):
            BusinessHour(start=dt_time(11, 0, 5))
        with pytest.raises(ValueError):
            BusinessHour(start='AAA')
        with pytest.raises(ValueError):
            BusinessHour(start='14:00:05')

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = self._offset()
        offset2 = self._offset()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset1) == '<BusinessHour: BH=09:00-17:00>'
        assert repr(self.offset2) == '<3 * BusinessHours: BH=09:00-17:00>'
        assert repr(self.offset3) == '<-1 * BusinessHour: BH=09:00-17:00>'
        assert repr(self.offset4) == '<-4 * BusinessHours: BH=09:00-17:00>'

        assert repr(self.offset5) == '<BusinessHour: BH=11:00-14:30>'
        assert repr(self.offset6) == '<BusinessHour: BH=20:00-05:00>'
        assert repr(self.offset7) == '<-2 * BusinessHours: BH=21:30-06:30>'

    def test_with_offset(self):
        expected = Timestamp('2014-07-01 13:00')

        assert self.d + BusinessHour() * 3 == expected
        assert self.d + BusinessHour(n=3) == expected

    def testEQ(self):
        for offset in [self.offset1, self.offset2, self.offset3, self.offset4]:
            assert offset == offset

        assert BusinessHour() != BusinessHour(-1)
        assert BusinessHour(start='09:00') == BusinessHour()
        assert BusinessHour(start='09:00') != BusinessHour(start='09:01')
        assert (BusinessHour(start='09:00', end='17:00') !=
                BusinessHour(start='17:00', end='09:01'))

    def test_hash(self):
        for offset in [self.offset1, self.offset2, self.offset3, self.offset4]:
            assert hash(offset) == hash(offset)

    def testCall(self):
        assert self.offset1(self.d) == datetime(2014, 7, 1, 11)
        assert self.offset2(self.d) == datetime(2014, 7, 1, 13)
        assert self.offset3(self.d) == datetime(2014, 6, 30, 17)
        assert self.offset4(self.d) == datetime(2014, 6, 30, 14)

    def testRAdd(self):
        assert self.d + self.offset2 == self.offset2 + self.d

    def testSub(self):
        off = self.offset2
        pytest.raises(Exception, off.__sub__, self.d)
        assert 2 * off - off == off

        assert self.d - self.offset2 == self.d + self._offset(-3)

    def testRSub(self):
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def testMult1(self):
        assert self.d + 5 * self.offset1 == self.d + self._offset(5)

    def testMult2(self):
        assert self.d + (-3 * self._offset(-2)) == self.d + self._offset(6)

    def testRollback1(self):
        assert self.offset1.rollback(self.d) == self.d
        assert self.offset2.rollback(self.d) == self.d
        assert self.offset3.rollback(self.d) == self.d
        assert self.offset4.rollback(self.d) == self.d
        assert self.offset5.rollback(self.d) == datetime(2014, 6, 30, 14, 30)
        assert self.offset6.rollback(self.d) == datetime(2014, 7, 1, 5, 0)
        assert self.offset7.rollback(self.d) == datetime(2014, 7, 1, 6, 30)

        d = datetime(2014, 7, 1, 0)
        assert self.offset1.rollback(d) == datetime(2014, 6, 30, 17)
        assert self.offset2.rollback(d) == datetime(2014, 6, 30, 17)
        assert self.offset3.rollback(d) == datetime(2014, 6, 30, 17)
        assert self.offset4.rollback(d) == datetime(2014, 6, 30, 17)
        assert self.offset5.rollback(d) == datetime(2014, 6, 30, 14, 30)
        assert self.offset6.rollback(d) == d
        assert self.offset7.rollback(d) == d

        assert self._offset(5).rollback(self.d) == self.d

    def testRollback2(self):
        assert (self._offset(-3).rollback(datetime(2014, 7, 5, 15, 0)) ==
                datetime(2014, 7, 4, 17, 0))

    def testRollforward1(self):
        assert self.offset1.rollforward(self.d) == self.d
        assert self.offset2.rollforward(self.d) == self.d
        assert self.offset3.rollforward(self.d) == self.d
        assert self.offset4.rollforward(self.d) == self.d
        assert (self.offset5.rollforward(self.d) ==
                datetime(2014, 7, 1, 11, 0))
        assert (self.offset6.rollforward(self.d) ==
                datetime(2014, 7, 1, 20, 0))
        assert (self.offset7.rollforward(self.d) ==
                datetime(2014, 7, 1, 21, 30))

        d = datetime(2014, 7, 1, 0)
        assert self.offset1.rollforward(d) == datetime(2014, 7, 1, 9)
        assert self.offset2.rollforward(d) == datetime(2014, 7, 1, 9)
        assert self.offset3.rollforward(d) == datetime(2014, 7, 1, 9)
        assert self.offset4.rollforward(d) == datetime(2014, 7, 1, 9)
        assert self.offset5.rollforward(d) == datetime(2014, 7, 1, 11)
        assert self.offset6.rollforward(d) == d
        assert self.offset7.rollforward(d) == d

        assert self._offset(5).rollforward(self.d) == self.d

    def testRollforward2(self):
        assert (self._offset(-3).rollforward(datetime(2014, 7, 5, 16, 0)) ==
                datetime(2014, 7, 7, 9))

    def test_roll_date_object(self):
        offset = BusinessHour()

        dt = datetime(2014, 7, 6, 15, 0)

        result = offset.rollback(dt)
        assert result == datetime(2014, 7, 4, 17)

        result = offset.rollforward(dt)
        assert result == datetime(2014, 7, 7, 9)

    normalize_cases = []
    normalize_cases.append((BusinessHour(normalize=True), {
        datetime(2014, 7, 1, 8): datetime(2014, 7, 1),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 2),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 2),
        datetime(2014, 7, 1, 23): datetime(2014, 7, 2),
        datetime(2014, 7, 1, 0): datetime(2014, 7, 1),
        datetime(2014, 7, 4, 15): datetime(2014, 7, 4),
        datetime(2014, 7, 4, 15, 59): datetime(2014, 7, 4),
        datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7),
        datetime(2014, 7, 5, 23): datetime(2014, 7, 7),
        datetime(2014, 7, 6, 10): datetime(2014, 7, 7)}))

    normalize_cases.append((BusinessHour(-1, normalize=True), {
        datetime(2014, 7, 1, 8): datetime(2014, 6, 30),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 1),
        datetime(2014, 7, 1, 10): datetime(2014, 6, 30),
        datetime(2014, 7, 1, 0): datetime(2014, 6, 30),
        datetime(2014, 7, 7, 10): datetime(2014, 7, 4),
        datetime(2014, 7, 7, 10, 1): datetime(2014, 7, 7),
        datetime(2014, 7, 5, 23): datetime(2014, 7, 4),
        datetime(2014, 7, 6, 10): datetime(2014, 7, 4)}))

    normalize_cases.append((BusinessHour(1, normalize=True, start='17:00',
                                         end='04:00'), {
        datetime(2014, 7, 1, 8): datetime(2014, 7, 1),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1),
        datetime(2014, 7, 1, 23): datetime(2014, 7, 2),
        datetime(2014, 7, 2, 2): datetime(2014, 7, 2),
        datetime(2014, 7, 2, 3): datetime(2014, 7, 2),
        datetime(2014, 7, 4, 23): datetime(2014, 7, 5),
        datetime(2014, 7, 5, 2): datetime(2014, 7, 5),
        datetime(2014, 7, 7, 2): datetime(2014, 7, 7),
        datetime(2014, 7, 7, 17): datetime(2014, 7, 7)}))

    @pytest.mark.parametrize('case', normalize_cases)
    def test_normalize(self, case):
        offset, cases = case
        for dt, expected in compat.iteritems(cases):
            assert offset.apply(dt) == expected

    on_offset_cases = []
    on_offset_cases.append((BusinessHour(), {
        datetime(2014, 7, 1, 9): True,
        datetime(2014, 7, 1, 8, 59): False,
        datetime(2014, 7, 1, 8): False,
        datetime(2014, 7, 1, 17): True,
        datetime(2014, 7, 1, 17, 1): False,
        datetime(2014, 7, 1, 18): False,
        datetime(2014, 7, 5, 9): False,
        datetime(2014, 7, 6, 12): False}))

    on_offset_cases.append((BusinessHour(start='10:00', end='15:00'), {
        datetime(2014, 7, 1, 9): False,
        datetime(2014, 7, 1, 10): True,
        datetime(2014, 7, 1, 15): True,
        datetime(2014, 7, 1, 15, 1): False,
        datetime(2014, 7, 5, 12): False,
        datetime(2014, 7, 6, 12): False}))

    on_offset_cases.append((BusinessHour(start='19:00', end='05:00'), {
        datetime(2014, 7, 1, 9, 0): False,
        datetime(2014, 7, 1, 10, 0): False,
        datetime(2014, 7, 1, 15): False,
        datetime(2014, 7, 1, 15, 1): False,
        datetime(2014, 7, 5, 12, 0): False,
        datetime(2014, 7, 6, 12, 0): False,
        datetime(2014, 7, 1, 19, 0): True,
        datetime(2014, 7, 2, 0, 0): True,
        datetime(2014, 7, 4, 23): True,
        datetime(2014, 7, 5, 1): True,
        datetime(2014, 7, 5, 5, 0): True,
        datetime(2014, 7, 6, 23, 0): False,
        datetime(2014, 7, 7, 3, 0): False}))

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, cases = case
        for dt, expected in compat.iteritems(cases):
            assert offset.onOffset(dt) == expected

    opening_time_cases = []
    # opening time should be affected by sign of n, not by n's value and
    # end
    opening_time_cases.append(([BusinessHour(), BusinessHour(n=2),
                                BusinessHour(n=4), BusinessHour(end='10:00'),
                                BusinessHour(n=2, end='4:00'),
                                BusinessHour(n=4, end='15:00')], {
        datetime(2014, 7, 1, 11): (datetime(2014, 7, 2, 9),
                                   datetime(2014, 7, 1, 9)),
        datetime(2014, 7, 1, 18): (datetime(2014, 7, 2, 9),
                                   datetime(2014, 7, 1, 9)),
        datetime(2014, 7, 1, 23): (datetime(2014, 7, 2, 9),
                                   datetime(2014, 7, 1, 9)),
        datetime(2014, 7, 2, 8): (datetime(2014, 7, 2, 9),
                                  datetime(2014, 7, 1, 9)),
        # if timestamp is on opening time, next opening time is
        # as it is
        datetime(2014, 7, 2, 9): (datetime(2014, 7, 2, 9),
                                  datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 2, 10): (datetime(2014, 7, 3, 9),
                                   datetime(2014, 7, 2, 9)),
        # 2014-07-05 is saturday
        datetime(2014, 7, 5, 10): (datetime(2014, 7, 7, 9),
                                   datetime(2014, 7, 4, 9)),
        datetime(2014, 7, 4, 10): (datetime(2014, 7, 7, 9),
                                   datetime(2014, 7, 4, 9)),
        datetime(2014, 7, 4, 23): (datetime(2014, 7, 7, 9),
                                   datetime(2014, 7, 4, 9)),
        datetime(2014, 7, 6, 10): (datetime(2014, 7, 7, 9),
                                   datetime(2014, 7, 4, 9)),
        datetime(2014, 7, 7, 5): (datetime(2014, 7, 7, 9),
                                  datetime(2014, 7, 4, 9)),
        datetime(2014, 7, 7, 9, 1): (datetime(2014, 7, 8, 9),
                                     datetime(2014, 7, 7, 9))}))

    opening_time_cases.append(([BusinessHour(start='11:15'),
                                BusinessHour(n=2, start='11:15'),
                                BusinessHour(n=3, start='11:15'),
                                BusinessHour(start='11:15', end='10:00'),
                                BusinessHour(n=2, start='11:15', end='4:00'),
                                BusinessHour(n=3, start='11:15',
                                             end='15:00')], {
        datetime(2014, 7, 1, 11): (datetime(2014, 7, 1, 11, 15),
                                   datetime(2014, 6, 30, 11, 15)),
        datetime(2014, 7, 1, 18): (datetime(2014, 7, 2, 11, 15),
                                   datetime(2014, 7, 1, 11, 15)),
        datetime(2014, 7, 1, 23): (datetime(2014, 7, 2, 11, 15),
                                   datetime(2014, 7, 1, 11, 15)),
        datetime(2014, 7, 2, 8): (datetime(2014, 7, 2, 11, 15),
                                  datetime(2014, 7, 1, 11, 15)),
        datetime(2014, 7, 2, 9): (datetime(2014, 7, 2, 11, 15),
                                  datetime(2014, 7, 1, 11, 15)),
        datetime(2014, 7, 2, 10): (datetime(2014, 7, 2, 11, 15),
                                   datetime(2014, 7, 1, 11, 15)),
        datetime(2014, 7, 2, 11, 15): (datetime(2014, 7, 2, 11, 15),
                                       datetime(2014, 7, 2, 11, 15)),
        datetime(2014, 7, 2, 11, 15, 1): (datetime(2014, 7, 3, 11, 15),
                                          datetime(2014, 7, 2, 11, 15)),
        datetime(2014, 7, 5, 10): (datetime(2014, 7, 7, 11, 15),
                                   datetime(2014, 7, 4, 11, 15)),
        datetime(2014, 7, 4, 10): (datetime(2014, 7, 4, 11, 15),
                                   datetime(2014, 7, 3, 11, 15)),
        datetime(2014, 7, 4, 23): (datetime(2014, 7, 7, 11, 15),
                                   datetime(2014, 7, 4, 11, 15)),
        datetime(2014, 7, 6, 10): (datetime(2014, 7, 7, 11, 15),
                                   datetime(2014, 7, 4, 11, 15)),
        datetime(2014, 7, 7, 5): (datetime(2014, 7, 7, 11, 15),
                                  datetime(2014, 7, 4, 11, 15)),
        datetime(2014, 7, 7, 9, 1): (datetime(2014, 7, 7, 11, 15),
                                     datetime(2014, 7, 4, 11, 15))}))

    opening_time_cases.append(([BusinessHour(-1), BusinessHour(n=-2),
                                BusinessHour(n=-4),
                                BusinessHour(n=-1, end='10:00'),
                                BusinessHour(n=-2, end='4:00'),
                                BusinessHour(n=-4, end='15:00')], {
        datetime(2014, 7, 1, 11): (datetime(2014, 7, 1, 9),
                                   datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 1, 18): (datetime(2014, 7, 1, 9),
                                   datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 1, 23): (datetime(2014, 7, 1, 9),
                                   datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 2, 8): (datetime(2014, 7, 1, 9),
                                  datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 2, 9): (datetime(2014, 7, 2, 9),
                                  datetime(2014, 7, 2, 9)),
        datetime(2014, 7, 2, 10): (datetime(2014, 7, 2, 9),
                                   datetime(2014, 7, 3, 9)),
        datetime(2014, 7, 5, 10): (datetime(2014, 7, 4, 9),
                                   datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 4, 10): (datetime(2014, 7, 4, 9),
                                   datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 4, 23): (datetime(2014, 7, 4, 9),
                                   datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 6, 10): (datetime(2014, 7, 4, 9),
                                   datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 7, 5): (datetime(2014, 7, 4, 9),
                                  datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 7, 9): (datetime(2014, 7, 7, 9),
                                  datetime(2014, 7, 7, 9)),
        datetime(2014, 7, 7, 9, 1): (datetime(2014, 7, 7, 9),
                                     datetime(2014, 7, 8, 9))}))

    opening_time_cases.append(([BusinessHour(start='17:00', end='05:00'),
                                BusinessHour(n=3, start='17:00',
                                             end='03:00')], {
        datetime(2014, 7, 1, 11): (datetime(2014, 7, 1, 17),
                                   datetime(2014, 6, 30, 17)),
        datetime(2014, 7, 1, 18): (datetime(2014, 7, 2, 17),
                                   datetime(2014, 7, 1, 17)),
        datetime(2014, 7, 1, 23): (datetime(2014, 7, 2, 17),
                                   datetime(2014, 7, 1, 17)),
        datetime(2014, 7, 2, 8): (datetime(2014, 7, 2, 17),
                                  datetime(2014, 7, 1, 17)),
        datetime(2014, 7, 2, 9): (datetime(2014, 7, 2, 17),
                                  datetime(2014, 7, 1, 17)),
        datetime(2014, 7, 4, 17): (datetime(2014, 7, 4, 17),
                                   datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 5, 10): (datetime(2014, 7, 7, 17),
                                   datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 4, 10): (datetime(2014, 7, 4, 17),
                                   datetime(2014, 7, 3, 17)),
        datetime(2014, 7, 4, 23): (datetime(2014, 7, 7, 17),
                                   datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 6, 10): (datetime(2014, 7, 7, 17),
                                   datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 7, 5): (datetime(2014, 7, 7, 17),
                                  datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 7, 17, 1): (datetime(2014, 7, 8, 17),
                                      datetime(2014, 7, 7, 17)), }))

    opening_time_cases.append(([BusinessHour(-1, start='17:00', end='05:00'),
                                BusinessHour(n=-2, start='17:00',
                                             end='03:00')], {
        datetime(2014, 7, 1, 11): (datetime(2014, 6, 30, 17),
                                   datetime(2014, 7, 1, 17)),
        datetime(2014, 7, 1, 18): (datetime(2014, 7, 1, 17),
                                   datetime(2014, 7, 2, 17)),
        datetime(2014, 7, 1, 23): (datetime(2014, 7, 1, 17),
                                   datetime(2014, 7, 2, 17)),
        datetime(2014, 7, 2, 8): (datetime(2014, 7, 1, 17),
                                  datetime(2014, 7, 2, 17)),
        datetime(2014, 7, 2, 9): (datetime(2014, 7, 1, 17),
                                  datetime(2014, 7, 2, 17)),
        datetime(2014, 7, 2, 16, 59): (datetime(2014, 7, 1, 17),
                                       datetime(2014, 7, 2, 17)),
        datetime(2014, 7, 5, 10): (datetime(2014, 7, 4, 17),
                                   datetime(2014, 7, 7, 17)),
        datetime(2014, 7, 4, 10): (datetime(2014, 7, 3, 17),
                                   datetime(2014, 7, 4, 17)),
        datetime(2014, 7, 4, 23): (datetime(2014, 7, 4, 17),
                                   datetime(2014, 7, 7, 17)),
        datetime(2014, 7, 6, 10): (datetime(2014, 7, 4, 17),
                                   datetime(2014, 7, 7, 17)),
        datetime(2014, 7, 7, 5): (datetime(2014, 7, 4, 17),
                                  datetime(2014, 7, 7, 17)),
        datetime(2014, 7, 7, 18): (datetime(2014, 7, 7, 17),
                                   datetime(2014, 7, 8, 17))}))

    @pytest.mark.parametrize('case', opening_time_cases)
    def test_opening_time(self, case):
        _offsets, cases = case
        for offset in _offsets:
            for dt, (exp_next, exp_prev) in compat.iteritems(cases):
                assert offset._next_opening_time(dt) == exp_next
                assert offset._prev_opening_time(dt) == exp_prev

    apply_cases = []
    apply_cases.append((BusinessHour(), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 12),
        datetime(2014, 7, 1, 13): datetime(2014, 7, 1, 14),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 1, 19): datetime(2014, 7, 2, 10),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 2, 9),
        datetime(2014, 7, 1, 16, 30, 15): datetime(2014, 7, 2, 9, 30, 15),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 2, 10),
        datetime(2014, 7, 2, 11): datetime(2014, 7, 2, 12),
        # out of business hours
        datetime(2014, 7, 2, 8): datetime(2014, 7, 2, 10),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 3, 10),
        datetime(2014, 7, 2, 23): datetime(2014, 7, 3, 10),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 10),
        # saturday
        datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 10),
        datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 10),
        datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7, 9, 30),
        datetime(2014, 7, 4, 16, 30, 30): datetime(2014, 7, 7, 9, 30, 30)}))

    apply_cases.append((BusinessHour(4), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 15),
        datetime(2014, 7, 1, 13): datetime(2014, 7, 2, 9),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 2, 11),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 2, 12),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 2, 11): datetime(2014, 7, 2, 15),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 3, 13),
        datetime(2014, 7, 2, 23): datetime(2014, 7, 3, 13),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 13),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 13),
        datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 13),
        datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7, 12, 30),
        datetime(2014, 7, 4, 16, 30, 30): datetime(2014, 7, 7, 12, 30, 30)}))

    apply_cases.append((BusinessHour(-1), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 10),
        datetime(2014, 7, 1, 13): datetime(2014, 7, 1, 12),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 1, 14),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 1, 15),
        datetime(2014, 7, 1, 10): datetime(2014, 6, 30, 17),
        datetime(2014, 7, 1, 16, 30, 15): datetime(2014, 7, 1, 15, 30, 15),
        datetime(2014, 7, 1, 9, 30, 15): datetime(2014, 6, 30, 16, 30, 15),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 1, 5): datetime(2014, 6, 30, 16),
        datetime(2014, 7, 2, 11): datetime(2014, 7, 2, 10),
        # out of business hours
        datetime(2014, 7, 2, 8): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 16),
        datetime(2014, 7, 2, 23): datetime(2014, 7, 2, 16),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 2, 16),
        # saturday
        datetime(2014, 7, 5, 15): datetime(2014, 7, 4, 16),
        datetime(2014, 7, 7, 9): datetime(2014, 7, 4, 16),
        datetime(2014, 7, 7, 9, 30): datetime(2014, 7, 4, 16, 30),
        datetime(2014, 7, 7, 9, 30, 30): datetime(2014, 7, 4, 16, 30, 30)}))

    apply_cases.append((BusinessHour(-4), {
        datetime(2014, 7, 1, 11): datetime(2014, 6, 30, 15),
        datetime(2014, 7, 1, 13): datetime(2014, 6, 30, 17),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 1, 11),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 1, 12),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1, 13),
        datetime(2014, 7, 2, 11): datetime(2014, 7, 1, 15),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 1, 13),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 2, 23): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 4, 13),
        datetime(2014, 7, 4, 18): datetime(2014, 7, 4, 13),
        datetime(2014, 7, 7, 9, 30): datetime(2014, 7, 4, 13, 30),
        datetime(2014, 7, 7, 9, 30, 30): datetime(2014, 7, 4, 13, 30, 30)}))

    apply_cases.append((BusinessHour(start='13:00', end='16:00'), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 14),
        datetime(2014, 7, 1, 13): datetime(2014, 7, 1, 14),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 1, 19): datetime(2014, 7, 2, 14),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 2, 14),
        datetime(2014, 7, 1, 15, 30, 15): datetime(2014, 7, 2, 13, 30, 15),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 14),
        datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 14)}))

    apply_cases.append((BusinessHour(n=2, start='13:00', end='16:00'), {
        datetime(2014, 7, 1, 17): datetime(2014, 7, 2, 15),
        datetime(2014, 7, 2, 14): datetime(2014, 7, 3, 13),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 2, 15),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 3, 15),
        datetime(2014, 7, 2, 14, 30): datetime(2014, 7, 3, 13, 30),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 15),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 15),
        datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 15),
        datetime(2014, 7, 4, 14, 30): datetime(2014, 7, 7, 13, 30),
        datetime(2014, 7, 4, 14, 30, 30): datetime(2014, 7, 7, 13, 30, 30)}))

    apply_cases.append((BusinessHour(n=-1, start='13:00', end='16:00'), {
        datetime(2014, 7, 2, 11): datetime(2014, 7, 1, 15),
        datetime(2014, 7, 2, 13): datetime(2014, 7, 1, 15),
        datetime(2014, 7, 2, 14): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 2, 15): datetime(2014, 7, 2, 14),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 15),
        datetime(2014, 7, 2, 16): datetime(2014, 7, 2, 15),
        datetime(2014, 7, 2, 13, 30, 15): datetime(2014, 7, 1, 15, 30, 15),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 4, 15),
        datetime(2014, 7, 7, 11): datetime(2014, 7, 4, 15)}))

    apply_cases.append((BusinessHour(n=-3, start='10:00', end='16:00'), {
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1, 13),
        datetime(2014, 7, 2, 14): datetime(2014, 7, 2, 11),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 1, 13),
        datetime(2014, 7, 2, 13): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 2, 11, 30): datetime(2014, 7, 1, 14, 30),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 2, 13),
        datetime(2014, 7, 4, 10): datetime(2014, 7, 3, 13),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 4, 13),
        datetime(2014, 7, 4, 16): datetime(2014, 7, 4, 13),
        datetime(2014, 7, 4, 12, 30): datetime(2014, 7, 3, 15, 30),
        datetime(2014, 7, 4, 12, 30, 30): datetime(2014, 7, 3, 15, 30, 30)}))

    apply_cases.append((BusinessHour(start='19:00', end='05:00'), {
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1, 20),
        datetime(2014, 7, 2, 14): datetime(2014, 7, 2, 20),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 2, 20),
        datetime(2014, 7, 2, 13): datetime(2014, 7, 2, 20),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 20),
        datetime(2014, 7, 2, 4, 30): datetime(2014, 7, 2, 19, 30),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 1),
        datetime(2014, 7, 4, 10): datetime(2014, 7, 4, 20),
        datetime(2014, 7, 4, 23): datetime(2014, 7, 5, 0),
        datetime(2014, 7, 5, 0): datetime(2014, 7, 5, 1),
        datetime(2014, 7, 5, 4): datetime(2014, 7, 7, 19),
        datetime(2014, 7, 5, 4, 30): datetime(2014, 7, 7, 19, 30),
        datetime(2014, 7, 5, 4, 30, 30): datetime(2014, 7, 7, 19, 30, 30)}))

    apply_cases.append((BusinessHour(n=-1, start='19:00', end='05:00'), {
        datetime(2014, 7, 1, 17): datetime(2014, 7, 1, 4),
        datetime(2014, 7, 2, 14): datetime(2014, 7, 2, 4),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 2, 4),
        datetime(2014, 7, 2, 13): datetime(2014, 7, 2, 4),
        datetime(2014, 7, 2, 20): datetime(2014, 7, 2, 5),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 2, 4),
        datetime(2014, 7, 2, 19, 30): datetime(2014, 7, 2, 4, 30),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 2, 23),
        datetime(2014, 7, 3, 6): datetime(2014, 7, 3, 4),
        datetime(2014, 7, 4, 23): datetime(2014, 7, 4, 22),
        datetime(2014, 7, 5, 0): datetime(2014, 7, 4, 23),
        datetime(2014, 7, 5, 4): datetime(2014, 7, 5, 3),
        datetime(2014, 7, 7, 19, 30): datetime(2014, 7, 5, 4, 30),
        datetime(2014, 7, 7, 19, 30, 30): datetime(2014, 7, 5, 4, 30, 30)}))

    @pytest.mark.parametrize('case', apply_cases)
    def test_apply(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    apply_large_n_cases = []
    # A week later
    apply_large_n_cases.append((BusinessHour(40), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 8, 11),
        datetime(2014, 7, 1, 13): datetime(2014, 7, 8, 13),
        datetime(2014, 7, 1, 15): datetime(2014, 7, 8, 15),
        datetime(2014, 7, 1, 16): datetime(2014, 7, 8, 16),
        datetime(2014, 7, 1, 17): datetime(2014, 7, 9, 9),
        datetime(2014, 7, 2, 11): datetime(2014, 7, 9, 11),
        datetime(2014, 7, 2, 8): datetime(2014, 7, 9, 9),
        datetime(2014, 7, 2, 19): datetime(2014, 7, 10, 9),
        datetime(2014, 7, 2, 23): datetime(2014, 7, 10, 9),
        datetime(2014, 7, 3, 0): datetime(2014, 7, 10, 9),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 14, 9),
        datetime(2014, 7, 4, 18): datetime(2014, 7, 14, 9),
        datetime(2014, 7, 7, 9, 30): datetime(2014, 7, 14, 9, 30),
        datetime(2014, 7, 7, 9, 30, 30): datetime(2014, 7, 14, 9, 30, 30)}))

    # 3 days and 1 hour before
    apply_large_n_cases.append((BusinessHour(-25), {
        datetime(2014, 7, 1, 11): datetime(2014, 6, 26, 10),
        datetime(2014, 7, 1, 13): datetime(2014, 6, 26, 12),
        datetime(2014, 7, 1, 9): datetime(2014, 6, 25, 16),
        datetime(2014, 7, 1, 10): datetime(2014, 6, 25, 17),
        datetime(2014, 7, 3, 11): datetime(2014, 6, 30, 10),
        datetime(2014, 7, 3, 8): datetime(2014, 6, 27, 16),
        datetime(2014, 7, 3, 19): datetime(2014, 6, 30, 16),
        datetime(2014, 7, 3, 23): datetime(2014, 6, 30, 16),
        datetime(2014, 7, 4, 9): datetime(2014, 6, 30, 16),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 6, 18): datetime(2014, 7, 1, 16),
        datetime(2014, 7, 7, 9, 30): datetime(2014, 7, 1, 16, 30),
        datetime(2014, 7, 7, 10, 30, 30): datetime(2014, 7, 2, 9, 30, 30)}))

    # 5 days and 3 hours later
    apply_large_n_cases.append((BusinessHour(28, start='21:00', end='02:00'), {
        datetime(2014, 7, 1, 11): datetime(2014, 7, 9, 0),
        datetime(2014, 7, 1, 22): datetime(2014, 7, 9, 1),
        datetime(2014, 7, 1, 23): datetime(2014, 7, 9, 21),
        datetime(2014, 7, 2, 2): datetime(2014, 7, 10, 0),
        datetime(2014, 7, 3, 21): datetime(2014, 7, 11, 0),
        datetime(2014, 7, 4, 1): datetime(2014, 7, 11, 23),
        datetime(2014, 7, 4, 2): datetime(2014, 7, 12, 0),
        datetime(2014, 7, 4, 3): datetime(2014, 7, 12, 0),
        datetime(2014, 7, 5, 1): datetime(2014, 7, 14, 23),
        datetime(2014, 7, 5, 15): datetime(2014, 7, 15, 0),
        datetime(2014, 7, 6, 18): datetime(2014, 7, 15, 0),
        datetime(2014, 7, 7, 1): datetime(2014, 7, 15, 0),
        datetime(2014, 7, 7, 23, 30): datetime(2014, 7, 15, 21, 30)}))

    @pytest.mark.parametrize('case', apply_large_n_cases)
    def test_apply_large_n(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_apply_nanoseconds(self):
        tests = []

        tests.append((BusinessHour(),
                      {Timestamp('2014-07-04 15:00') + Nano(5): Timestamp(
                          '2014-07-04 16:00') + Nano(5),
                       Timestamp('2014-07-04 16:00') + Nano(5): Timestamp(
                           '2014-07-07 09:00') + Nano(5),
                       Timestamp('2014-07-04 16:00') - Nano(5): Timestamp(
                           '2014-07-04 17:00') - Nano(5)}))

        tests.append((BusinessHour(-1),
                      {Timestamp('2014-07-04 15:00') + Nano(5): Timestamp(
                          '2014-07-04 14:00') + Nano(5),
                       Timestamp('2014-07-04 10:00') + Nano(5): Timestamp(
                           '2014-07-04 09:00') + Nano(5),
                       Timestamp('2014-07-04 10:00') - Nano(5): Timestamp(
                           '2014-07-03 17:00') - Nano(5), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assert_offset_equal(offset, base, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = self._offset()
        offset2 = self._offset()
        assert not offset1 != offset2

    def test_datetimeindex(self):
        idx1 = DatetimeIndex(start='2014-07-04 15:00', end='2014-07-08 10:00',
                             freq='BH')
        idx2 = DatetimeIndex(start='2014-07-04 15:00', periods=12, freq='BH')
        idx3 = DatetimeIndex(end='2014-07-08 10:00', periods=12, freq='BH')
        expected = DatetimeIndex(['2014-07-04 15:00', '2014-07-04 16:00',
                                  '2014-07-07 09:00',
                                  '2014-07-07 10:00', '2014-07-07 11:00',
                                  '2014-07-07 12:00',
                                  '2014-07-07 13:00', '2014-07-07 14:00',
                                  '2014-07-07 15:00',
                                  '2014-07-07 16:00', '2014-07-08 09:00',
                                  '2014-07-08 10:00'],
                                 freq='BH')
        for idx in [idx1, idx2, idx3]:
            tm.assert_index_equal(idx, expected)

        idx1 = DatetimeIndex(start='2014-07-04 15:45', end='2014-07-08 10:45',
                             freq='BH')
        idx2 = DatetimeIndex(start='2014-07-04 15:45', periods=12, freq='BH')
        idx3 = DatetimeIndex(end='2014-07-08 10:45', periods=12, freq='BH')

        expected = DatetimeIndex(['2014-07-04 15:45', '2014-07-04 16:45',
                                  '2014-07-07 09:45',
                                  '2014-07-07 10:45', '2014-07-07 11:45',
                                  '2014-07-07 12:45',
                                  '2014-07-07 13:45', '2014-07-07 14:45',
                                  '2014-07-07 15:45',
                                  '2014-07-07 16:45', '2014-07-08 09:45',
                                  '2014-07-08 10:45'],
                                 freq='BH')
        expected = idx1
        for idx in [idx1, idx2, idx3]:
            tm.assert_index_equal(idx, expected)


class TestCustomBusinessHour(Base):
    _offset = CustomBusinessHour

    def setup_method(self, method):
        # 2014 Calendar to check custom holidays
        #   Sun Mon Tue Wed Thu Fri Sat
        #  6/22  23  24  25  26  27  28
        #    29  30 7/1   2   3   4   5
        #     6   7   8   9  10  11  12
        self.d = datetime(2014, 7, 1, 10, 00)
        self.offset1 = CustomBusinessHour(weekmask='Tue Wed Thu Fri')

        self.holidays = ['2014-06-27', datetime(2014, 6, 30),
                         np.datetime64('2014-07-02')]
        self.offset2 = CustomBusinessHour(holidays=self.holidays)

    def test_constructor_errors(self):
        from datetime import time as dt_time
        with pytest.raises(ValueError):
            CustomBusinessHour(start=dt_time(11, 0, 5))
        with pytest.raises(ValueError):
            CustomBusinessHour(start='AAA')
        with pytest.raises(ValueError):
            CustomBusinessHour(start='14:00:05')

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = self._offset()
        offset2 = self._offset()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset1) == '<CustomBusinessHour: CBH=09:00-17:00>'
        assert repr(self.offset2) == '<CustomBusinessHour: CBH=09:00-17:00>'

    def test_with_offset(self):
        expected = Timestamp('2014-07-01 13:00')

        assert self.d + CustomBusinessHour() * 3 == expected
        assert self.d + CustomBusinessHour(n=3) == expected

    def testEQ(self):
        for offset in [self.offset1, self.offset2]:
            assert offset == offset

        assert CustomBusinessHour() != CustomBusinessHour(-1)
        assert (CustomBusinessHour(start='09:00') ==
                CustomBusinessHour())
        assert (CustomBusinessHour(start='09:00') !=
                CustomBusinessHour(start='09:01'))
        assert (CustomBusinessHour(start='09:00', end='17:00') !=
                CustomBusinessHour(start='17:00', end='09:01'))

        assert (CustomBusinessHour(weekmask='Tue Wed Thu Fri') !=
                CustomBusinessHour(weekmask='Mon Tue Wed Thu Fri'))
        assert (CustomBusinessHour(holidays=['2014-06-27']) !=
                CustomBusinessHour(holidays=['2014-06-28']))

    def test_hash(self):
        assert hash(self.offset1) == hash(self.offset1)
        assert hash(self.offset2) == hash(self.offset2)

    def testCall(self):
        assert self.offset1(self.d) == datetime(2014, 7, 1, 11)
        assert self.offset2(self.d) == datetime(2014, 7, 1, 11)

    def testRAdd(self):
        assert self.d + self.offset2 == self.offset2 + self.d

    def testSub(self):
        off = self.offset2
        pytest.raises(Exception, off.__sub__, self.d)
        assert 2 * off - off == off

        assert self.d - self.offset2 == self.d - (2 * off - off)

    def testRSub(self):
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def testMult1(self):
        assert self.d + 5 * self.offset1 == self.d + self._offset(5)

    def testMult2(self):
        assert self.d + (-3 * self._offset(-2)) == self.d + self._offset(6)

    def testRollback1(self):
        assert self.offset1.rollback(self.d) == self.d
        assert self.offset2.rollback(self.d) == self.d

        d = datetime(2014, 7, 1, 0)

        # 2014/07/01 is Tuesday, 06/30 is Monday(holiday)
        assert self.offset1.rollback(d) == datetime(2014, 6, 27, 17)

        # 2014/6/30 and 2014/6/27 are holidays
        assert self.offset2.rollback(d) == datetime(2014, 6, 26, 17)

    def testRollback2(self):
        assert (self._offset(-3).rollback(datetime(2014, 7, 5, 15, 0)) ==
                datetime(2014, 7, 4, 17, 0))

    def testRollforward1(self):
        assert self.offset1.rollforward(self.d) == self.d
        assert self.offset2.rollforward(self.d) == self.d

        d = datetime(2014, 7, 1, 0)
        assert self.offset1.rollforward(d) == datetime(2014, 7, 1, 9)
        assert self.offset2.rollforward(d) == datetime(2014, 7, 1, 9)

    def testRollforward2(self):
        assert (self._offset(-3).rollforward(datetime(2014, 7, 5, 16, 0)) ==
                datetime(2014, 7, 7, 9))

    def test_roll_date_object(self):
        offset = BusinessHour()

        dt = datetime(2014, 7, 6, 15, 0)

        result = offset.rollback(dt)
        assert result == datetime(2014, 7, 4, 17)

        result = offset.rollforward(dt)
        assert result == datetime(2014, 7, 7, 9)

    def test_normalize(self):
        tests = []

        tests.append((CustomBusinessHour(normalize=True,
                                         holidays=self.holidays),
                      {datetime(2014, 7, 1, 8): datetime(2014, 7, 1),
                       datetime(2014, 7, 1, 17): datetime(2014, 7, 3),
                       datetime(2014, 7, 1, 16): datetime(2014, 7, 3),
                       datetime(2014, 7, 1, 23): datetime(2014, 7, 3),
                       datetime(2014, 7, 1, 0): datetime(2014, 7, 1),
                       datetime(2014, 7, 4, 15): datetime(2014, 7, 4),
                       datetime(2014, 7, 4, 15, 59): datetime(2014, 7, 4),
                       datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7),
                       datetime(2014, 7, 5, 23): datetime(2014, 7, 7),
                       datetime(2014, 7, 6, 10): datetime(2014, 7, 7)}))

        tests.append((CustomBusinessHour(-1, normalize=True,
                                         holidays=self.holidays),
                      {datetime(2014, 7, 1, 8): datetime(2014, 6, 26),
                       datetime(2014, 7, 1, 17): datetime(2014, 7, 1),
                       datetime(2014, 7, 1, 16): datetime(2014, 7, 1),
                       datetime(2014, 7, 1, 10): datetime(2014, 6, 26),
                       datetime(2014, 7, 1, 0): datetime(2014, 6, 26),
                       datetime(2014, 7, 7, 10): datetime(2014, 7, 4),
                       datetime(2014, 7, 7, 10, 1): datetime(2014, 7, 7),
                       datetime(2014, 7, 5, 23): datetime(2014, 7, 4),
                       datetime(2014, 7, 6, 10): datetime(2014, 7, 4)}))

        tests.append((CustomBusinessHour(1, normalize=True, start='17:00',
                                         end='04:00', holidays=self.holidays),
                      {datetime(2014, 7, 1, 8): datetime(2014, 7, 1),
                       datetime(2014, 7, 1, 17): datetime(2014, 7, 1),
                       datetime(2014, 7, 1, 23): datetime(2014, 7, 2),
                       datetime(2014, 7, 2, 2): datetime(2014, 7, 2),
                       datetime(2014, 7, 2, 3): datetime(2014, 7, 3),
                       datetime(2014, 7, 4, 23): datetime(2014, 7, 5),
                       datetime(2014, 7, 5, 2): datetime(2014, 7, 5),
                       datetime(2014, 7, 7, 2): datetime(2014, 7, 7),
                       datetime(2014, 7, 7, 17): datetime(2014, 7, 7)}))

        for offset, cases in tests:
            for dt, expected in compat.iteritems(cases):
                assert offset.apply(dt) == expected

    def test_onOffset(self):
        tests = []

        tests.append((CustomBusinessHour(start='10:00', end='15:00',
                                         holidays=self.holidays),
                      {datetime(2014, 7, 1, 9): False,
                       datetime(2014, 7, 1, 10): True,
                       datetime(2014, 7, 1, 15): True,
                       datetime(2014, 7, 1, 15, 1): False,
                       datetime(2014, 7, 5, 12): False,
                       datetime(2014, 7, 6, 12): False}))

        for offset, cases in tests:
            for dt, expected in compat.iteritems(cases):
                assert offset.onOffset(dt) == expected

    def test_apply(self):
        tests = []

        tests.append((
            CustomBusinessHour(holidays=self.holidays),
            {datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 12),
             datetime(2014, 7, 1, 13): datetime(2014, 7, 1, 14),
             datetime(2014, 7, 1, 15): datetime(2014, 7, 1, 16),
             datetime(2014, 7, 1, 19): datetime(2014, 7, 3, 10),
             datetime(2014, 7, 1, 16): datetime(2014, 7, 3, 9),
             datetime(2014, 7, 1, 16, 30, 15): datetime(2014, 7, 3, 9, 30, 15),
             datetime(2014, 7, 1, 17): datetime(2014, 7, 3, 10),
             datetime(2014, 7, 2, 11): datetime(2014, 7, 3, 10),
             # out of business hours
             datetime(2014, 7, 2, 8): datetime(2014, 7, 3, 10),
             datetime(2014, 7, 2, 19): datetime(2014, 7, 3, 10),
             datetime(2014, 7, 2, 23): datetime(2014, 7, 3, 10),
             datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 10),
             # saturday
             datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 10),
             datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 10),
             datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7, 9, 30),
             datetime(2014, 7, 4, 16, 30, 30): datetime(2014, 7, 7, 9, 30,
                                                        30)}))

        tests.append((
            CustomBusinessHour(4, holidays=self.holidays),
            {datetime(2014, 7, 1, 11): datetime(2014, 7, 1, 15),
             datetime(2014, 7, 1, 13): datetime(2014, 7, 3, 9),
             datetime(2014, 7, 1, 15): datetime(2014, 7, 3, 11),
             datetime(2014, 7, 1, 16): datetime(2014, 7, 3, 12),
             datetime(2014, 7, 1, 17): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 2, 11): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 2, 8): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 2, 19): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 2, 23): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 3, 0): datetime(2014, 7, 3, 13),
             datetime(2014, 7, 5, 15): datetime(2014, 7, 7, 13),
             datetime(2014, 7, 4, 17): datetime(2014, 7, 7, 13),
             datetime(2014, 7, 4, 16, 30): datetime(2014, 7, 7, 12, 30),
             datetime(2014, 7, 4, 16, 30, 30): datetime(2014, 7, 7, 12, 30,
                                                        30)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assert_offset_equal(offset, base, expected)

    def test_apply_nanoseconds(self):
        tests = []

        tests.append((CustomBusinessHour(holidays=self.holidays),
                      {Timestamp('2014-07-01 15:00') + Nano(5): Timestamp(
                          '2014-07-01 16:00') + Nano(5),
                       Timestamp('2014-07-01 16:00') + Nano(5): Timestamp(
                           '2014-07-03 09:00') + Nano(5),
                       Timestamp('2014-07-01 16:00') - Nano(5): Timestamp(
                           '2014-07-01 17:00') - Nano(5)}))

        tests.append((CustomBusinessHour(-1, holidays=self.holidays),
                      {Timestamp('2014-07-01 15:00') + Nano(5): Timestamp(
                          '2014-07-01 14:00') + Nano(5),
                       Timestamp('2014-07-01 10:00') + Nano(5): Timestamp(
                           '2014-07-01 09:00') + Nano(5),
                       Timestamp('2014-07-01 10:00') - Nano(5): Timestamp(
                           '2014-06-26 17:00') - Nano(5), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assert_offset_equal(offset, base, expected)


class TestCustomBusinessDay(Base):
    _offset = CDay

    def setup_method(self, method):
        self.d = datetime(2008, 1, 1)
        self.nd = np_datetime64_compat('2008-01-01 00:00:00Z')

        self.offset = CDay()
        self.offset2 = CDay(2)

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = CDay()
        offset2 = CDay()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset) == '<CustomBusinessDay>'
        assert repr(self.offset2) == '<2 * CustomBusinessDays>'

        expected = '<BusinessDay: offset=datetime.timedelta(1)>'
        assert repr(self.offset + timedelta(1)) == expected

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    def testEQ(self):
        assert self.offset2 == self.offset2

    def test_mul(self):
        pass

    def test_hash(self):
        assert hash(self.offset2) == hash(self.offset2)

    def testCall(self):
        assert self.offset2(self.d) == datetime(2008, 1, 3)
        assert self.offset2(self.nd) == datetime(2008, 1, 3)

    def testRAdd(self):
        assert self.d + self.offset2 == self.offset2 + self.d

    def testSub(self):
        off = self.offset2
        pytest.raises(Exception, off.__sub__, self.d)
        assert 2 * off - off == off

        assert self.d - self.offset2 == self.d + CDay(-2)

    def testRSub(self):
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def testMult1(self):
        assert self.d + 10 * self.offset == self.d + CDay(10)

    def testMult2(self):
        assert self.d + (-5 * CDay(-10)) == self.d + CDay(50)

    def testRollback1(self):
        assert CDay(10).rollback(self.d) == self.d

    def testRollback2(self):
        assert (CDay(10).rollback(datetime(2008, 1, 5)) ==
                datetime(2008, 1, 4))

    def testRollforward1(self):
        assert CDay(10).rollforward(self.d) == self.d

    def testRollforward2(self):
        assert (CDay(10).rollforward(datetime(2008, 1, 5)) ==
                datetime(2008, 1, 7))

    def test_roll_date_object(self):
        offset = CDay()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 14)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 17)

        offset = offsets.Day()
        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 15)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 15)

    on_offset_cases = [(CDay(), datetime(2008, 1, 1), True),
                       (CDay(), datetime(2008, 1, 5), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, d, expected = case
        assert_onOffset(offset, d, expected)

    apply_cases = []
    apply_cases.append((CDay(), {
        datetime(2008, 1, 1): datetime(2008, 1, 2),
        datetime(2008, 1, 4): datetime(2008, 1, 7),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 8)}))

    apply_cases.append((2 * CDay(), {
        datetime(2008, 1, 1): datetime(2008, 1, 3),
        datetime(2008, 1, 4): datetime(2008, 1, 8),
        datetime(2008, 1, 5): datetime(2008, 1, 8),
        datetime(2008, 1, 6): datetime(2008, 1, 8),
        datetime(2008, 1, 7): datetime(2008, 1, 9)}))

    apply_cases.append((-CDay(), {
        datetime(2008, 1, 1): datetime(2007, 12, 31),
        datetime(2008, 1, 4): datetime(2008, 1, 3),
        datetime(2008, 1, 5): datetime(2008, 1, 4),
        datetime(2008, 1, 6): datetime(2008, 1, 4),
        datetime(2008, 1, 7): datetime(2008, 1, 4),
        datetime(2008, 1, 8): datetime(2008, 1, 7)}))

    apply_cases.append((-2 * CDay(), {
        datetime(2008, 1, 1): datetime(2007, 12, 28),
        datetime(2008, 1, 4): datetime(2008, 1, 2),
        datetime(2008, 1, 5): datetime(2008, 1, 3),
        datetime(2008, 1, 6): datetime(2008, 1, 3),
        datetime(2008, 1, 7): datetime(2008, 1, 3),
        datetime(2008, 1, 8): datetime(2008, 1, 4),
        datetime(2008, 1, 9): datetime(2008, 1, 7)}))

    apply_cases.append((CDay(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2008, 1, 4): datetime(2008, 1, 4),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 7)}))

    @pytest.mark.parametrize('case', apply_cases)
    def test_apply(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + CDay(10)
        assert result == datetime(2012, 11, 6)

        result = dt + CDay(100) - CDay(100)
        assert result == dt

        off = CDay() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 12, 23)
        assert rs == xp

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2011, 12, 26)
        assert rs == xp

    def test_apply_corner(self):
        pytest.raises(Exception, CDay().apply, BMonthEnd())

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = CDay()
        offset2 = CDay()
        assert not offset1 != offset2

    def test_holidays(self):
        # Define a TradingDay offset
        holidays = ['2012-05-01', datetime(2013, 5, 1),
                    np.datetime64('2014-05-01')]
        tday = CDay(holidays=holidays)
        for year in range(2012, 2015):
            dt = datetime(year, 4, 30)
            xp = datetime(year, 5, 2)
            rs = dt + tday
            assert rs == xp

    def test_weekmask(self):
        weekmask_saudi = 'Sat Sun Mon Tue Wed'  # Thu-Fri Weekend
        weekmask_uae = '1111001'  # Fri-Sat Weekend
        weekmask_egypt = [1, 1, 1, 1, 0, 0, 1]  # Fri-Sat Weekend
        bday_saudi = CDay(weekmask=weekmask_saudi)
        bday_uae = CDay(weekmask=weekmask_uae)
        bday_egypt = CDay(weekmask=weekmask_egypt)
        dt = datetime(2013, 5, 1)
        xp_saudi = datetime(2013, 5, 4)
        xp_uae = datetime(2013, 5, 2)
        xp_egypt = datetime(2013, 5, 2)
        assert xp_saudi == dt + bday_saudi
        assert xp_uae == dt + bday_uae
        assert xp_egypt == dt + bday_egypt
        xp2 = datetime(2013, 5, 5)
        assert xp2 == dt + 2 * bday_saudi
        assert xp2 == dt + 2 * bday_uae
        assert xp2 == dt + 2 * bday_egypt

    def test_weekmask_and_holidays(self):
        weekmask_egypt = 'Sun Mon Tue Wed Thu'  # Fri-Sat Weekend
        holidays = ['2012-05-01', datetime(2013, 5, 1),
                    np.datetime64('2014-05-01')]
        bday_egypt = CDay(holidays=holidays, weekmask=weekmask_egypt)
        dt = datetime(2013, 4, 30)
        xp_egypt = datetime(2013, 5, 5)
        assert xp_egypt == dt + 2 * bday_egypt

    def test_calendar(self):
        calendar = USFederalHolidayCalendar()
        dt = datetime(2014, 1, 17)
        assert_offset_equal(CDay(calendar=calendar), dt, datetime(2014, 1, 21))

    def test_roundtrip_pickle(self):
        def _check_roundtrip(obj):
            unpickled = tm.round_trip_pickle(obj)
            assert unpickled == obj

        _check_roundtrip(self.offset)
        _check_roundtrip(self.offset2)
        _check_roundtrip(self.offset * 2)

    def test_pickle_compat_0_14_1(self):
        hdays = [datetime(2013, 1, 1) for ele in range(4)]

        pth = tm.get_data_path()

        cday0_14_1 = read_pickle(os.path.join(pth, 'cday-0.14.1.pickle'))
        cday = CDay(holidays=hdays)
        assert cday == cday0_14_1


class CustomBusinessMonthBase(object):

    def setup_method(self, method):
        self.d = datetime(2008, 1, 1)

        self.offset = self._object()
        self.offset2 = self._object(2)

    def testEQ(self):
        assert self.offset2 == self.offset2

    def test_mul(self):
        pass

    def test_hash(self):
        assert hash(self.offset2) == hash(self.offset2)

    def testRAdd(self):
        assert self.d + self.offset2 == self.offset2 + self.d

    def testSub(self):
        off = self.offset2
        pytest.raises(Exception, off.__sub__, self.d)
        assert 2 * off - off == off

        assert self.d - self.offset2 == self.d + self._object(-2)

    def testRSub(self):
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def testMult1(self):
        assert self.d + 10 * self.offset == self.d + self._object(10)

    def testMult2(self):
        assert self.d + (-5 * self._object(-10)) == self.d + self._object(50)

    def test_offsets_compare_equal(self):
        offset1 = self._object()
        offset2 = self._object()
        assert not offset1 != offset2

    def test_roundtrip_pickle(self):
        def _check_roundtrip(obj):
            unpickled = tm.round_trip_pickle(obj)
            assert unpickled == obj

        _check_roundtrip(self._object())
        _check_roundtrip(self._object(2))
        _check_roundtrip(self._object() * 2)

    def test_copy(self):
        # GH 17452
        off = self._object(weekmask='Mon Wed Fri')
        assert off == off.copy()


class TestCustomBusinessMonthEnd(CustomBusinessMonthBase, Base):
    _object = CBMonthEnd

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = CBMonthEnd()
        offset2 = CBMonthEnd()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset) == '<CustomBusinessMonthEnd>'
        assert repr(self.offset2) == '<2 * CustomBusinessMonthEnds>'

    def testCall(self):
        assert self.offset2(self.d) == datetime(2008, 2, 29)

    def testRollback1(self):
        assert (CDay(10).rollback(datetime(2007, 12, 31)) ==
                datetime(2007, 12, 31))

    def testRollback2(self):
        assert CBMonthEnd(10).rollback(self.d) == datetime(2007, 12, 31)

    def testRollforward1(self):
        assert CBMonthEnd(10).rollforward(self.d) == datetime(2008, 1, 31)

    def test_roll_date_object(self):
        offset = CBMonthEnd()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        assert result == datetime(2012, 8, 31)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 28)

        offset = offsets.Day()
        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 15)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 15)

    on_offset_cases = [(CBMonthEnd(), datetime(2008, 1, 31), True),
                       (CBMonthEnd(), datetime(2008, 1, 1), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, d, expected = case
        assert_onOffset(offset, d, expected)

    apply_cases = []
    apply_cases.append((CBMonthEnd(), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 2, 7): datetime(2008, 2, 29)}))

    apply_cases.append((2 * CBMonthEnd(), {
        datetime(2008, 1, 1): datetime(2008, 2, 29),
        datetime(2008, 2, 7): datetime(2008, 3, 31)}))

    apply_cases.append((-CBMonthEnd(), {
        datetime(2008, 1, 1): datetime(2007, 12, 31),
        datetime(2008, 2, 8): datetime(2008, 1, 31)}))

    apply_cases.append((-2 * CBMonthEnd(), {
        datetime(2008, 1, 1): datetime(2007, 11, 30),
        datetime(2008, 2, 9): datetime(2007, 12, 31)}))

    apply_cases.append((CBMonthEnd(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 2, 7): datetime(2008, 2, 29)}))

    @pytest.mark.parametrize('case', apply_cases)
    def test_apply(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + CBMonthEnd(10)
        assert result == datetime(2013, 7, 31)

        result = dt + CDay(100) - CDay(100)
        assert result == dt

        off = CBMonthEnd() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 7, 29)
        assert rs == xp

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2012, 5, 31)
        assert rs == xp

    def test_holidays(self):
        # Define a TradingDay offset
        holidays = ['2012-01-31', datetime(2012, 2, 28),
                    np.datetime64('2012-02-29')]
        bm_offset = CBMonthEnd(holidays=holidays)
        dt = datetime(2012, 1, 1)
        assert dt + bm_offset == datetime(2012, 1, 30)
        assert dt + 2 * bm_offset == datetime(2012, 2, 27)

    def test_datetimeindex(self):
        from pandas.tseries.holiday import USFederalHolidayCalendar
        hcal = USFederalHolidayCalendar()
        freq = CBMonthEnd(calendar=hcal)

        assert (DatetimeIndex(start='20120101', end='20130101',
                              freq=freq).tolist()[0] == datetime(2012, 1, 31))


class TestCustomBusinessMonthBegin(CustomBusinessMonthBase, Base):
    _object = CBMonthBegin

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = CBMonthBegin()
        offset2 = CBMonthBegin()
        offset2.normalize = True
        assert offset == offset2

    def test_repr(self):
        assert repr(self.offset) == '<CustomBusinessMonthBegin>'
        assert repr(self.offset2) == '<2 * CustomBusinessMonthBegins>'

    def testCall(self):
        assert self.offset2(self.d) == datetime(2008, 3, 3)

    def testRollback1(self):
        assert (CDay(10).rollback(datetime(2007, 12, 31)) ==
                datetime(2007, 12, 31))

    def testRollback2(self):
        assert CBMonthBegin(10).rollback(self.d) == datetime(2008, 1, 1)

    def testRollforward1(self):
        assert CBMonthBegin(10).rollforward(self.d) == datetime(2008, 1, 1)

    def test_roll_date_object(self):
        offset = CBMonthBegin()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 3)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 10, 1)

        offset = offsets.Day()
        result = offset.rollback(dt)
        assert result == datetime(2012, 9, 15)

        result = offset.rollforward(dt)
        assert result == datetime(2012, 9, 15)

    on_offset_cases = [(CBMonthBegin(), datetime(2008, 1, 1), True),
                       (CBMonthBegin(), datetime(2008, 1, 31), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, dt, expected = case
        assert_onOffset(offset, dt, expected)

    apply_cases = []
    apply_cases.append((CBMonthBegin(), {
        datetime(2008, 1, 1): datetime(2008, 2, 1),
        datetime(2008, 2, 7): datetime(2008, 3, 3)}))

    apply_cases.append((2 * CBMonthBegin(), {
        datetime(2008, 1, 1): datetime(2008, 3, 3),
        datetime(2008, 2, 7): datetime(2008, 4, 1)}))

    apply_cases.append((-CBMonthBegin(), {
        datetime(2008, 1, 1): datetime(2007, 12, 3),
        datetime(2008, 2, 8): datetime(2008, 2, 1)}))

    apply_cases.append((-2 * CBMonthBegin(), {
        datetime(2008, 1, 1): datetime(2007, 11, 1),
        datetime(2008, 2, 9): datetime(2008, 1, 1)}))

    apply_cases.append((CBMonthBegin(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2008, 1, 7): datetime(2008, 2, 1)}))

    @pytest.mark.parametrize('case', apply_cases)
    def test_apply(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + CBMonthBegin(10)
        assert result == datetime(2013, 8, 1)

        result = dt + CDay(100) - CDay(100)
        assert result == dt

        off = CBMonthBegin() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 7, 1)
        assert rs == xp

        st = datetime(2011, 12, 18)
        rs = st + off

        xp = datetime(2012, 6, 1)
        assert rs == xp

    def test_holidays(self):
        # Define a TradingDay offset
        holidays = ['2012-02-01', datetime(2012, 2, 2),
                    np.datetime64('2012-03-01')]
        bm_offset = CBMonthBegin(holidays=holidays)
        dt = datetime(2012, 1, 1)

        assert dt + bm_offset == datetime(2012, 1, 2)
        assert dt + 2 * bm_offset == datetime(2012, 2, 3)

    def test_datetimeindex(self):
        hcal = USFederalHolidayCalendar()
        cbmb = CBMonthBegin(calendar=hcal)
        assert (DatetimeIndex(start='20120101', end='20130101',
                              freq=cbmb).tolist()[0] == datetime(2012, 1, 3))


class TestWeek(Base):
    _offset = Week

    def test_repr(self):
        assert repr(Week(weekday=0)) == "<Week: weekday=0>"
        assert repr(Week(n=-1, weekday=0)) == "<-1 * Week: weekday=0>"
        assert repr(Week(n=-2, weekday=0)) == "<-2 * Weeks: weekday=0>"

    def test_corner(self):
        pytest.raises(ValueError, Week, weekday=7)
        tm.assert_raises_regex(
            ValueError, "Day must be", Week, weekday=-1)

    def test_isAnchored(self):
        assert Week(weekday=0).isAnchored()
        assert not Week().isAnchored()
        assert not Week(2, weekday=2).isAnchored()
        assert not Week(2).isAnchored()

    offset_cases = []
    # not business week
    offset_cases.append((Week(), {
        datetime(2008, 1, 1): datetime(2008, 1, 8),
        datetime(2008, 1, 4): datetime(2008, 1, 11),
        datetime(2008, 1, 5): datetime(2008, 1, 12),
        datetime(2008, 1, 6): datetime(2008, 1, 13),
        datetime(2008, 1, 7): datetime(2008, 1, 14)}))

    # Mon
    offset_cases.append((Week(weekday=0), {
        datetime(2007, 12, 31): datetime(2008, 1, 7),
        datetime(2008, 1, 4): datetime(2008, 1, 7),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 14)}))

    # n=0 -> roll forward. Mon
    offset_cases.append((Week(0, weekday=0), {
        datetime(2007, 12, 31): datetime(2007, 12, 31),
        datetime(2008, 1, 4): datetime(2008, 1, 7),
        datetime(2008, 1, 5): datetime(2008, 1, 7),
        datetime(2008, 1, 6): datetime(2008, 1, 7),
        datetime(2008, 1, 7): datetime(2008, 1, 7)}))

    # n=0 -> roll forward. Mon
    offset_cases.append((Week(-2, weekday=1), {
        datetime(2010, 4, 6): datetime(2010, 3, 23),
        datetime(2010, 4, 8): datetime(2010, 3, 30),
        datetime(2010, 4, 5): datetime(2010, 3, 23)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_onOffset(self):
        for weekday in range(7):
            offset = Week(weekday=weekday)

            for day in range(1, 8):
                date = datetime(2008, 1, day)

                if day % 7 == weekday:
                    expected = True
                else:
                    expected = False
            assert_onOffset(offset, date, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = Week()
        offset2 = Week()
        assert not offset1 != offset2


class TestWeekOfMonth(Base):
    _offset = WeekOfMonth

    def test_constructor(self):
        tm.assert_raises_regex(ValueError, "^Week", WeekOfMonth,
                               n=1, week=4, weekday=0)
        tm.assert_raises_regex(ValueError, "^Week", WeekOfMonth,
                               n=1, week=-1, weekday=0)
        tm.assert_raises_regex(ValueError, "^Day", WeekOfMonth,
                               n=1, week=0, weekday=-1)
        tm.assert_raises_regex(ValueError, "^Day", WeekOfMonth,
                               n=1, week=0, weekday=7)

    def test_repr(self):
        assert (repr(WeekOfMonth(weekday=1, week=2)) ==
                "<WeekOfMonth: week=2, weekday=1>")

    def test_offset(self):
        date1 = datetime(2011, 1, 4)  # 1st Tuesday of Month
        date2 = datetime(2011, 1, 11)  # 2nd Tuesday of Month
        date3 = datetime(2011, 1, 18)  # 3rd Tuesday of Month
        date4 = datetime(2011, 1, 25)  # 4th Tuesday of Month

        # see for loop for structure
        test_cases = [
            (-2, 2, 1, date1, datetime(2010, 11, 16)),
            (-2, 2, 1, date2, datetime(2010, 11, 16)),
            (-2, 2, 1, date3, datetime(2010, 11, 16)),
            (-2, 2, 1, date4, datetime(2010, 12, 21)),

            (-1, 2, 1, date1, datetime(2010, 12, 21)),
            (-1, 2, 1, date2, datetime(2010, 12, 21)),
            (-1, 2, 1, date3, datetime(2010, 12, 21)),
            (-1, 2, 1, date4, datetime(2011, 1, 18)),

            (0, 0, 1, date1, datetime(2011, 1, 4)),
            (0, 0, 1, date2, datetime(2011, 2, 1)),
            (0, 0, 1, date3, datetime(2011, 2, 1)),
            (0, 0, 1, date4, datetime(2011, 2, 1)),
            (0, 1, 1, date1, datetime(2011, 1, 11)),
            (0, 1, 1, date2, datetime(2011, 1, 11)),
            (0, 1, 1, date3, datetime(2011, 2, 8)),
            (0, 1, 1, date4, datetime(2011, 2, 8)),
            (0, 0, 1, date1, datetime(2011, 1, 4)),
            (0, 1, 1, date2, datetime(2011, 1, 11)),
            (0, 2, 1, date3, datetime(2011, 1, 18)),
            (0, 3, 1, date4, datetime(2011, 1, 25)),

            (1, 0, 0, date1, datetime(2011, 2, 7)),
            (1, 0, 0, date2, datetime(2011, 2, 7)),
            (1, 0, 0, date3, datetime(2011, 2, 7)),
            (1, 0, 0, date4, datetime(2011, 2, 7)),
            (1, 0, 1, date1, datetime(2011, 2, 1)),
            (1, 0, 1, date2, datetime(2011, 2, 1)),
            (1, 0, 1, date3, datetime(2011, 2, 1)),
            (1, 0, 1, date4, datetime(2011, 2, 1)),
            (1, 0, 2, date1, datetime(2011, 1, 5)),
            (1, 0, 2, date2, datetime(2011, 2, 2)),
            (1, 0, 2, date3, datetime(2011, 2, 2)),
            (1, 0, 2, date4, datetime(2011, 2, 2)),

            (1, 2, 1, date1, datetime(2011, 1, 18)),
            (1, 2, 1, date2, datetime(2011, 1, 18)),
            (1, 2, 1, date3, datetime(2011, 2, 15)),
            (1, 2, 1, date4, datetime(2011, 2, 15)),

            (2, 2, 1, date1, datetime(2011, 2, 15)),
            (2, 2, 1, date2, datetime(2011, 2, 15)),
            (2, 2, 1, date3, datetime(2011, 3, 15)),
            (2, 2, 1, date4, datetime(2011, 3, 15))]

        for n, week, weekday, dt, expected in test_cases:
            offset = WeekOfMonth(n, week=week, weekday=weekday)
            assert_offset_equal(offset, dt, expected)

        # try subtracting
        result = datetime(2011, 2, 1) - WeekOfMonth(week=1, weekday=2)
        assert result == datetime(2011, 1, 12)

        result = datetime(2011, 2, 3) - WeekOfMonth(week=0, weekday=2)
        assert result == datetime(2011, 2, 2)

    on_offset_cases = [(0, 0, datetime(2011, 2, 7), True),
                       (0, 0, datetime(2011, 2, 6), False),
                       (0, 0, datetime(2011, 2, 14), False),
                       (1, 0, datetime(2011, 2, 14), True),
                       (0, 1, datetime(2011, 2, 1), True),
                       (0, 1, datetime(2011, 2, 8), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        week, weekday, dt, expected = case
        offset = WeekOfMonth(week=week, weekday=weekday)
        assert offset.onOffset(dt) == expected


class TestLastWeekOfMonth(Base):
    _offset = LastWeekOfMonth

    def test_constructor(self):
        tm.assert_raises_regex(ValueError, "^N cannot be 0",
                               LastWeekOfMonth, n=0, weekday=1)

        tm.assert_raises_regex(ValueError, "^Day", LastWeekOfMonth, n=1,
                               weekday=-1)
        tm.assert_raises_regex(
            ValueError, "^Day", LastWeekOfMonth, n=1, weekday=7)

    def test_offset(self):
        # Saturday
        last_sat = datetime(2013, 8, 31)
        next_sat = datetime(2013, 9, 28)
        offset_sat = LastWeekOfMonth(n=1, weekday=5)

        one_day_before = (last_sat + timedelta(days=-1))
        assert one_day_before + offset_sat == last_sat

        one_day_after = (last_sat + timedelta(days=+1))
        assert one_day_after + offset_sat == next_sat

        # Test On that day
        assert last_sat + offset_sat == next_sat

        # Thursday

        offset_thur = LastWeekOfMonth(n=1, weekday=3)
        last_thurs = datetime(2013, 1, 31)
        next_thurs = datetime(2013, 2, 28)

        one_day_before = last_thurs + timedelta(days=-1)
        assert one_day_before + offset_thur == last_thurs

        one_day_after = last_thurs + timedelta(days=+1)
        assert one_day_after + offset_thur == next_thurs

        # Test on that day
        assert last_thurs + offset_thur == next_thurs

        three_before = last_thurs + timedelta(days=-3)
        assert three_before + offset_thur == last_thurs

        two_after = last_thurs + timedelta(days=+2)
        assert two_after + offset_thur == next_thurs

        offset_sunday = LastWeekOfMonth(n=1, weekday=WeekDay.SUN)
        assert datetime(2013, 7, 31) + offset_sunday == datetime(2013, 8, 25)

    on_offset_cases = [
        (WeekDay.SUN, datetime(2013, 1, 27), True),
        (WeekDay.SAT, datetime(2013, 3, 30), True),
        (WeekDay.MON, datetime(2013, 2, 18), False),  # Not the last Mon
        (WeekDay.SUN, datetime(2013, 2, 25), False),  # Not a SUN
        (WeekDay.MON, datetime(2013, 2, 25), True),
        (WeekDay.SAT, datetime(2013, 11, 30), True),

        (WeekDay.SAT, datetime(2006, 8, 26), True),
        (WeekDay.SAT, datetime(2007, 8, 25), True),
        (WeekDay.SAT, datetime(2008, 8, 30), True),
        (WeekDay.SAT, datetime(2009, 8, 29), True),
        (WeekDay.SAT, datetime(2010, 8, 28), True),
        (WeekDay.SAT, datetime(2011, 8, 27), True),
        (WeekDay.SAT, datetime(2019, 8, 31), True)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        weekday, dt, expected = case
        offset = LastWeekOfMonth(weekday=weekday)
        assert offset.onOffset(dt) == expected


class TestSemiMonthEnd(Base):
    _offset = SemiMonthEnd

    def test_offset_whole_year(self):
        dates = (datetime(2007, 12, 31),
                 datetime(2008, 1, 15),
                 datetime(2008, 1, 31),
                 datetime(2008, 2, 15),
                 datetime(2008, 2, 29),
                 datetime(2008, 3, 15),
                 datetime(2008, 3, 31),
                 datetime(2008, 4, 15),
                 datetime(2008, 4, 30),
                 datetime(2008, 5, 15),
                 datetime(2008, 5, 31),
                 datetime(2008, 6, 15),
                 datetime(2008, 6, 30),
                 datetime(2008, 7, 15),
                 datetime(2008, 7, 31),
                 datetime(2008, 8, 15),
                 datetime(2008, 8, 31),
                 datetime(2008, 9, 15),
                 datetime(2008, 9, 30),
                 datetime(2008, 10, 15),
                 datetime(2008, 10, 31),
                 datetime(2008, 11, 15),
                 datetime(2008, 11, 30),
                 datetime(2008, 12, 15),
                 datetime(2008, 12, 31))

        for base, exp_date in zip(dates[:-1], dates[1:]):
            assert_offset_equal(SemiMonthEnd(), base, exp_date)

        # ensure .apply_index works as expected
        s = DatetimeIndex(dates[:-1])
        result = SemiMonthEnd().apply_index(s)
        exp = DatetimeIndex(dates[1:])
        tm.assert_index_equal(result, exp)

        # ensure generating a range with DatetimeIndex gives same result
        result = DatetimeIndex(start=dates[0], end=dates[-1], freq='SM')
        exp = DatetimeIndex(dates)
        tm.assert_index_equal(result, exp)

    offset_cases = []
    offset_cases.append((SemiMonthEnd(), {
        datetime(2008, 1, 1): datetime(2008, 1, 15),
        datetime(2008, 1, 15): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 2, 15),
        datetime(2006, 12, 14): datetime(2006, 12, 15),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2007, 1, 15),
        datetime(2007, 1, 1): datetime(2007, 1, 15),
        datetime(2006, 12, 1): datetime(2006, 12, 15),
        datetime(2006, 12, 15): datetime(2006, 12, 31)}))

    offset_cases.append((SemiMonthEnd(day_of_month=20), {
        datetime(2008, 1, 1): datetime(2008, 1, 20),
        datetime(2008, 1, 15): datetime(2008, 1, 20),
        datetime(2008, 1, 21): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 2, 20),
        datetime(2006, 12, 14): datetime(2006, 12, 20),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2007, 1, 20),
        datetime(2007, 1, 1): datetime(2007, 1, 20),
        datetime(2006, 12, 1): datetime(2006, 12, 20),
        datetime(2006, 12, 15): datetime(2006, 12, 20)}))

    offset_cases.append((SemiMonthEnd(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 15),
        datetime(2008, 1, 16): datetime(2008, 1, 31),
        datetime(2008, 1, 15): datetime(2008, 1, 15),
        datetime(2008, 1, 31): datetime(2008, 1, 31),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2006, 12, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 15)}))

    offset_cases.append((SemiMonthEnd(0, day_of_month=16), {
        datetime(2008, 1, 1): datetime(2008, 1, 16),
        datetime(2008, 1, 16): datetime(2008, 1, 16),
        datetime(2008, 1, 15): datetime(2008, 1, 16),
        datetime(2008, 1, 31): datetime(2008, 1, 31),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2006, 12, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 16)}))

    offset_cases.append((SemiMonthEnd(2), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 2, 29),
        datetime(2006, 12, 29): datetime(2007, 1, 15),
        datetime(2006, 12, 31): datetime(2007, 1, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 31),
        datetime(2007, 1, 16): datetime(2007, 2, 15),
        datetime(2006, 11, 1): datetime(2006, 11, 30)}))

    offset_cases.append((SemiMonthEnd(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 31),
        datetime(2008, 6, 30): datetime(2008, 6, 15),
        datetime(2008, 12, 31): datetime(2008, 12, 15),
        datetime(2006, 12, 29): datetime(2006, 12, 15),
        datetime(2006, 12, 30): datetime(2006, 12, 15),
        datetime(2007, 1, 1): datetime(2006, 12, 31)}))

    offset_cases.append((SemiMonthEnd(-1, day_of_month=4), {
        datetime(2007, 1, 1): datetime(2006, 12, 31),
        datetime(2007, 1, 4): datetime(2006, 12, 31),
        datetime(2008, 6, 30): datetime(2008, 6, 4),
        datetime(2008, 12, 31): datetime(2008, 12, 4),
        datetime(2006, 12, 5): datetime(2006, 12, 4),
        datetime(2006, 12, 30): datetime(2006, 12, 4),
        datetime(2007, 1, 1): datetime(2006, 12, 31)}))

    offset_cases.append((SemiMonthEnd(-2), {
        datetime(2007, 1, 1): datetime(2006, 12, 15),
        datetime(2008, 6, 30): datetime(2008, 5, 31),
        datetime(2008, 3, 15): datetime(2008, 2, 15),
        datetime(2008, 12, 31): datetime(2008, 11, 30),
        datetime(2006, 12, 29): datetime(2006, 11, 30),
        datetime(2006, 12, 14): datetime(2006, 11, 15),
        datetime(2007, 1, 1): datetime(2006, 12, 15)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    @pytest.mark.parametrize('case', offset_cases)
    def test_apply_index(self, case):
        offset, cases = case
        s = DatetimeIndex(cases.keys())
        result = offset.apply_index(s)
        exp = DatetimeIndex(cases.values())
        tm.assert_index_equal(result, exp)

    on_offset_cases = [(datetime(2007, 12, 31), True),
                       (datetime(2007, 12, 15), True),
                       (datetime(2007, 12, 14), False),
                       (datetime(2007, 12, 1), False),
                       (datetime(2008, 2, 29), True)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        dt, expected = case
        assert_onOffset(SemiMonthEnd(), dt, expected)

    @pytest.mark.parametrize('klass,assert_func',
                             [(Series, tm.assert_series_equal),
                              (DatetimeIndex, tm.assert_index_equal)])
    def test_vectorized_offset_addition(self, klass, assert_func):
        s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                   Timestamp('2000-02-15', tz='US/Central')], name='a')

        result = s + SemiMonthEnd()
        result2 = SemiMonthEnd() + s
        exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                     Timestamp('2000-02-29', tz='US/Central')], name='a')
        assert_func(result, exp)
        assert_func(result2, exp)

        s = klass([Timestamp('2000-01-01 00:15:00', tz='US/Central'),
                   Timestamp('2000-02-01', tz='US/Central')], name='a')
        result = s + SemiMonthEnd()
        result2 = SemiMonthEnd() + s
        exp = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                     Timestamp('2000-02-15', tz='US/Central')], name='a')
        assert_func(result, exp)
        assert_func(result2, exp)


class TestSemiMonthBegin(Base):
    _offset = SemiMonthBegin

    def test_offset_whole_year(self):
        dates = (datetime(2007, 12, 15),
                 datetime(2008, 1, 1),
                 datetime(2008, 1, 15),
                 datetime(2008, 2, 1),
                 datetime(2008, 2, 15),
                 datetime(2008, 3, 1),
                 datetime(2008, 3, 15),
                 datetime(2008, 4, 1),
                 datetime(2008, 4, 15),
                 datetime(2008, 5, 1),
                 datetime(2008, 5, 15),
                 datetime(2008, 6, 1),
                 datetime(2008, 6, 15),
                 datetime(2008, 7, 1),
                 datetime(2008, 7, 15),
                 datetime(2008, 8, 1),
                 datetime(2008, 8, 15),
                 datetime(2008, 9, 1),
                 datetime(2008, 9, 15),
                 datetime(2008, 10, 1),
                 datetime(2008, 10, 15),
                 datetime(2008, 11, 1),
                 datetime(2008, 11, 15),
                 datetime(2008, 12, 1),
                 datetime(2008, 12, 15))

        for base, exp_date in zip(dates[:-1], dates[1:]):
            assert_offset_equal(SemiMonthBegin(), base, exp_date)

        # ensure .apply_index works as expected
        s = DatetimeIndex(dates[:-1])
        result = SemiMonthBegin().apply_index(s)
        exp = DatetimeIndex(dates[1:])
        tm.assert_index_equal(result, exp)

        # ensure generating a range with DatetimeIndex gives same result
        result = DatetimeIndex(start=dates[0], end=dates[-1], freq='SMS')
        exp = DatetimeIndex(dates)
        tm.assert_index_equal(result, exp)

    offset_cases = []
    offset_cases.append((SemiMonthBegin(), {
        datetime(2008, 1, 1): datetime(2008, 1, 15),
        datetime(2008, 1, 15): datetime(2008, 2, 1),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 14): datetime(2006, 12, 15),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2007, 1, 1): datetime(2007, 1, 15),
        datetime(2006, 12, 1): datetime(2006, 12, 15),
        datetime(2006, 12, 15): datetime(2007, 1, 1)}))

    offset_cases.append((SemiMonthBegin(day_of_month=20), {
        datetime(2008, 1, 1): datetime(2008, 1, 20),
        datetime(2008, 1, 15): datetime(2008, 1, 20),
        datetime(2008, 1, 21): datetime(2008, 2, 1),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 14): datetime(2006, 12, 20),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2007, 1, 1): datetime(2007, 1, 20),
        datetime(2006, 12, 1): datetime(2006, 12, 20),
        datetime(2006, 12, 15): datetime(2006, 12, 20)}))

    offset_cases.append((SemiMonthBegin(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2008, 1, 16): datetime(2008, 2, 1),
        datetime(2008, 1, 15): datetime(2008, 1, 15),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 2): datetime(2006, 12, 15),
        datetime(2007, 1, 1): datetime(2007, 1, 1)}))

    offset_cases.append((SemiMonthBegin(0, day_of_month=16), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2008, 1, 16): datetime(2008, 1, 16),
        datetime(2008, 1, 15): datetime(2008, 1, 16),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2007, 1, 5): datetime(2007, 1, 16),
        datetime(2007, 1, 1): datetime(2007, 1, 1)}))

    offset_cases.append((SemiMonthBegin(2), {
        datetime(2008, 1, 1): datetime(2008, 2, 1),
        datetime(2008, 1, 31): datetime(2008, 2, 15),
        datetime(2006, 12, 1): datetime(2007, 1, 1),
        datetime(2006, 12, 29): datetime(2007, 1, 15),
        datetime(2006, 12, 15): datetime(2007, 1, 15),
        datetime(2007, 1, 1): datetime(2007, 2, 1),
        datetime(2007, 1, 16): datetime(2007, 2, 15),
        datetime(2006, 11, 1): datetime(2006, 12, 1)}))

    offset_cases.append((SemiMonthBegin(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 15),
        datetime(2008, 6, 30): datetime(2008, 6, 15),
        datetime(2008, 6, 14): datetime(2008, 6, 1),
        datetime(2008, 12, 31): datetime(2008, 12, 15),
        datetime(2006, 12, 29): datetime(2006, 12, 15),
        datetime(2006, 12, 15): datetime(2006, 12, 1),
        datetime(2007, 1, 1): datetime(2006, 12, 15)}))

    offset_cases.append((SemiMonthBegin(-1, day_of_month=4), {
        datetime(2007, 1, 1): datetime(2006, 12, 4),
        datetime(2007, 1, 4): datetime(2007, 1, 1),
        datetime(2008, 6, 30): datetime(2008, 6, 4),
        datetime(2008, 12, 31): datetime(2008, 12, 4),
        datetime(2006, 12, 5): datetime(2006, 12, 4),
        datetime(2006, 12, 30): datetime(2006, 12, 4),
        datetime(2006, 12, 2): datetime(2006, 12, 1),
        datetime(2007, 1, 1): datetime(2006, 12, 4)}))

    offset_cases.append((SemiMonthBegin(-2), {
        datetime(2007, 1, 1): datetime(2006, 12, 1),
        datetime(2008, 6, 30): datetime(2008, 6, 1),
        datetime(2008, 6, 14): datetime(2008, 5, 15),
        datetime(2008, 12, 31): datetime(2008, 12, 1),
        datetime(2006, 12, 29): datetime(2006, 12, 1),
        datetime(2006, 12, 15): datetime(2006, 11, 15),
        datetime(2007, 1, 1): datetime(2006, 12, 1)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    @pytest.mark.parametrize('case', offset_cases)
    def test_apply_index(self, case):
        offset, cases = case
        s = DatetimeIndex(cases.keys())
        result = offset.apply_index(s)
        exp = DatetimeIndex(cases.values())
        tm.assert_index_equal(result, exp)

    on_offset_cases = [(datetime(2007, 12, 1), True),
                       (datetime(2007, 12, 15), True),
                       (datetime(2007, 12, 14), False),
                       (datetime(2007, 12, 31), False),
                       (datetime(2008, 2, 15), True)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        dt, expected = case
        assert_onOffset(SemiMonthBegin(), dt, expected)

    @pytest.mark.parametrize('klass,assert_func',
                             [(Series, tm.assert_series_equal),
                              (DatetimeIndex, tm.assert_index_equal)])
    def test_vectorized_offset_addition(self, klass, assert_func):
        s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                   Timestamp('2000-02-15', tz='US/Central')], name='a')
        result = s + SemiMonthBegin()
        result2 = SemiMonthBegin() + s
        exp = klass([Timestamp('2000-02-01 00:15:00', tz='US/Central'),
                     Timestamp('2000-03-01', tz='US/Central')], name='a')
        assert_func(result, exp)
        assert_func(result2, exp)

        s = klass([Timestamp('2000-01-01 00:15:00', tz='US/Central'),
                   Timestamp('2000-02-01', tz='US/Central')], name='a')
        result = s + SemiMonthBegin()
        result2 = SemiMonthBegin() + s
        exp = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                     Timestamp('2000-02-15', tz='US/Central')], name='a')
        assert_func(result, exp)
        assert_func(result2, exp)


def test_Easter():
    assert_offset_equal(Easter(), datetime(2010, 1, 1), datetime(2010, 4, 4))
    assert_offset_equal(Easter(), datetime(2010, 4, 5), datetime(2011, 4, 24))
    assert_offset_equal(Easter(2), datetime(2010, 1, 1), datetime(2011, 4, 24))

    assert_offset_equal(Easter(), datetime(2010, 4, 4), datetime(2011, 4, 24))
    assert_offset_equal(Easter(2), datetime(2010, 4, 4), datetime(2012, 4, 8))

    assert_offset_equal(-Easter(), datetime(2011, 1, 1), datetime(2010, 4, 4))
    assert_offset_equal(-Easter(), datetime(2010, 4, 5), datetime(2010, 4, 4))
    assert_offset_equal(-Easter(2),
                        datetime(2011, 1, 1),
                        datetime(2009, 4, 12))

    assert_offset_equal(-Easter(), datetime(2010, 4, 4), datetime(2009, 4, 12))
    assert_offset_equal(-Easter(2),
                        datetime(2010, 4, 4),
                        datetime(2008, 3, 23))


class TestOffsetNames(object):

    def test_get_offset_name(self):
        assert BDay().freqstr == 'B'
        assert BDay(2).freqstr == '2B'
        assert BMonthEnd().freqstr == 'BM'
        assert Week(weekday=0).freqstr == 'W-MON'
        assert Week(weekday=1).freqstr == 'W-TUE'
        assert Week(weekday=2).freqstr == 'W-WED'
        assert Week(weekday=3).freqstr == 'W-THU'
        assert Week(weekday=4).freqstr == 'W-FRI'

        assert LastWeekOfMonth(weekday=WeekDay.SUN).freqstr == "LWOM-SUN"


def test_get_offset():
    with tm.assert_raises_regex(ValueError, _INVALID_FREQ_ERROR):
        get_offset('gibberish')
    with tm.assert_raises_regex(ValueError, _INVALID_FREQ_ERROR):
        get_offset('QS-JAN-B')

    pairs = [
        ('B', BDay()), ('b', BDay()), ('bm', BMonthEnd()),
        ('Bm', BMonthEnd()), ('W-MON', Week(weekday=0)),
        ('W-TUE', Week(weekday=1)), ('W-WED', Week(weekday=2)),
        ('W-THU', Week(weekday=3)), ('W-FRI', Week(weekday=4))]

    for name, expected in pairs:
        offset = get_offset(name)
        assert offset == expected, ("Expected %r to yield %r (actual: %r)" %
                                    (name, expected, offset))


def test_get_offset_legacy():
    pairs = [('w@Sat', Week(weekday=5))]
    for name, expected in pairs:
        with tm.assert_raises_regex(ValueError, _INVALID_FREQ_ERROR):
            get_offset(name)


class TestOffsetAliases(object):

    def setup_method(self, method):
        _offset_map.clear()

    def test_alias_equality(self):
        for k, v in compat.iteritems(_offset_map):
            if v is None:
                continue
            assert k == v.copy()

    def test_rule_code(self):
        lst = ['M', 'MS', 'BM', 'BMS', 'D', 'B', 'H', 'T', 'S', 'L', 'U']
        for k in lst:
            assert k == get_offset(k).rule_code
            # should be cached - this is kind of an internals test...
            assert k in _offset_map
            assert k == (get_offset(k) * 3).rule_code

        suffix_lst = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
        base = 'W'
        for v in suffix_lst:
            alias = '-'.join([base, v])
            assert alias == get_offset(alias).rule_code
            assert alias == (get_offset(alias) * 5).rule_code

        suffix_lst = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG',
                      'SEP', 'OCT', 'NOV', 'DEC']
        base_lst = ['A', 'AS', 'BA', 'BAS', 'Q', 'QS', 'BQ', 'BQS']
        for base in base_lst:
            for v in suffix_lst:
                alias = '-'.join([base, v])
                assert alias == get_offset(alias).rule_code
                assert alias == (get_offset(alias) * 5).rule_code

        lst = ['M', 'D', 'B', 'H', 'T', 'S', 'L', 'U']
        for k in lst:
            code, stride = get_freq_code('3' + k)
            assert isinstance(code, int)
            assert stride == 3
            assert k == get_freq_str(code)


def test_dateoffset_misc():
    oset = offsets.DateOffset(months=2, days=4)
    # it works
    oset.freqstr

    assert (not offsets.DateOffset(months=2) == 2)


def test_freq_offsets():
    off = BDay(1, offset=timedelta(0, 1800))
    assert (off.freqstr == 'B+30Min')

    off = BDay(1, offset=timedelta(0, -1800))
    assert (off.freqstr == 'B-30Min')


def get_all_subclasses(cls):
    ret = set()
    this_subclasses = cls.__subclasses__()
    ret = ret | set(this_subclasses)
    for this_subclass in this_subclasses:
        ret | get_all_subclasses(this_subclass)
    return ret


class TestCaching(object):

    # as of GH 6479 (in 0.14.0), offset caching is turned off
    # as of v0.12.0 only BusinessMonth/Quarter were actually caching

    def setup_method(self, method):
        _daterange_cache.clear()
        _offset_map.clear()

    def run_X_index_creation(self, cls):
        inst1 = cls()
        if not inst1.isAnchored():
            assert not inst1._should_cache(), cls
            return

        assert inst1._should_cache(), cls

        DatetimeIndex(start=datetime(2013, 1, 31), end=datetime(2013, 3, 31),
                      freq=inst1, normalize=True)
        assert cls() in _daterange_cache, cls

    def test_should_cache_month_end(self):
        assert not MonthEnd()._should_cache()

    def test_should_cache_bmonth_end(self):
        assert not BusinessMonthEnd()._should_cache()

    def test_should_cache_week_month(self):
        assert not WeekOfMonth(weekday=1, week=2)._should_cache()

    def test_all_cacheableoffsets(self):
        for subclass in get_all_subclasses(CacheableOffset):
            if subclass.__name__[0] == "_" \
                    or subclass in TestCaching.no_simple_ctr:
                continue
            self.run_X_index_creation(subclass)

    def test_month_end_index_creation(self):
        DatetimeIndex(start=datetime(2013, 1, 31), end=datetime(2013, 3, 31),
                      freq=MonthEnd(), normalize=True)
        assert not MonthEnd() in _daterange_cache

    def test_bmonth_end_index_creation(self):
        DatetimeIndex(start=datetime(2013, 1, 31), end=datetime(2013, 3, 29),
                      freq=BusinessMonthEnd(), normalize=True)
        assert not BusinessMonthEnd() in _daterange_cache

    def test_week_of_month_index_creation(self):
        inst1 = WeekOfMonth(weekday=1, week=2)
        DatetimeIndex(start=datetime(2013, 1, 31), end=datetime(2013, 3, 29),
                      freq=inst1, normalize=True)
        inst2 = WeekOfMonth(weekday=1, week=2)
        assert inst2 not in _daterange_cache


class TestReprNames(object):

    def test_str_for_named_is_name(self):
        # look at all the amazing combinations!
        month_prefixes = ['A', 'AS', 'BA', 'BAS', 'Q', 'BQ', 'BQS', 'QS']
        names = [prefix + '-' + month
                 for prefix in month_prefixes
                 for month in ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
                               'AUG', 'SEP', 'OCT', 'NOV', 'DEC']]
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
        names += ['W-' + day for day in days]
        names += ['WOM-' + week + day
                  for week in ('1', '2', '3', '4') for day in days]
        _offset_map.clear()
        for name in names:
            offset = get_offset(name)
            assert offset.freqstr == name


def get_utc_offset_hours(ts):
    # take a Timestamp and compute total hours of utc offset
    o = ts.utcoffset()
    return (o.days * 24 * 3600 + o.seconds) / 3600.0


class TestDST(object):
    """
    test DateOffset additions over Daylight Savings Time
    """
    # one microsecond before the DST transition
    ts_pre_fallback = "2013-11-03 01:59:59.999999"
    ts_pre_springfwd = "2013-03-10 01:59:59.999999"

    # test both basic names and dateutil timezones
    timezone_utc_offsets = {
        'US/Eastern': dict(utc_offset_daylight=-4,
                           utc_offset_standard=-5, ),
        'dateutil/US/Pacific': dict(utc_offset_daylight=-7,
                                    utc_offset_standard=-8, )
    }
    valid_date_offsets_singular = [
        'weekday', 'day', 'hour', 'minute', 'second', 'microsecond'
    ]
    valid_date_offsets_plural = [
        'weeks', 'days',
        'hours', 'minutes', 'seconds',
        'milliseconds', 'microseconds'
    ]

    def _test_all_offsets(self, n, **kwds):
        valid_offsets = self.valid_date_offsets_plural if n > 1 \
            else self.valid_date_offsets_singular

        for name in valid_offsets:
            self._test_offset(offset_name=name, offset_n=n, **kwds)

    def _test_offset(self, offset_name, offset_n, tstart, expected_utc_offset):
        offset = DateOffset(**{offset_name: offset_n})

        t = tstart + offset
        if expected_utc_offset is not None:
            assert get_utc_offset_hours(t) == expected_utc_offset

        if offset_name == 'weeks':
            # dates should match
            assert t.date() == timedelta(days=7 * offset.kwds[
                'weeks']) + tstart.date()
            # expect the same day of week, hour of day, minute, second, ...
            assert (t.dayofweek == tstart.dayofweek and
                    t.hour == tstart.hour and
                    t.minute == tstart.minute and
                    t.second == tstart.second)
        elif offset_name == 'days':
            # dates should match
            assert timedelta(offset.kwds['days']) + tstart.date() == t.date()
            # expect the same hour of day, minute, second, ...
            assert (t.hour == tstart.hour and
                    t.minute == tstart.minute and
                    t.second == tstart.second)
        elif offset_name in self.valid_date_offsets_singular:
            # expect the singular offset value to match between tstart and t
            datepart_offset = getattr(t, offset_name
                                      if offset_name != 'weekday' else
                                      'dayofweek')
            assert datepart_offset == offset.kwds[offset_name]
        else:
            # the offset should be the same as if it was done in UTC
            assert (t == (tstart.tz_convert('UTC') + offset)
                    .tz_convert('US/Pacific'))

    def _make_timestamp(self, string, hrs_offset, tz):
        if hrs_offset >= 0:
            offset_string = '{hrs:02d}00'.format(hrs=hrs_offset)
        else:
            offset_string = '-{hrs:02d}00'.format(hrs=-1 * hrs_offset)
        return Timestamp(string + offset_string).tz_convert(tz)

    def test_fallback_plural(self):
        # test moving from daylight savings to standard time
        import dateutil
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets['utc_offset_daylight']
            hrs_post = utc_offsets['utc_offset_standard']

            if LooseVersion(dateutil.__version__) < LooseVersion('2.6.0'):
                # buggy ambiguous behavior in 2.6.0
                # GH 14621
                # https://github.com/dateutil/dateutil/issues/321
                self._test_all_offsets(
                    n=3, tstart=self._make_timestamp(self.ts_pre_fallback,
                                                     hrs_pre, tz),
                    expected_utc_offset=hrs_post)
            elif LooseVersion(dateutil.__version__) > LooseVersion('2.6.0'):
                # fixed, but skip the test
                continue

    def test_springforward_plural(self):
        # test moving from standard to daylight savings
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets['utc_offset_standard']
            hrs_post = utc_offsets['utc_offset_daylight']
            self._test_all_offsets(
                n=3, tstart=self._make_timestamp(self.ts_pre_springfwd,
                                                 hrs_pre, tz),
                expected_utc_offset=hrs_post)

    def test_fallback_singular(self):
        # in the case of singular offsets, we don't necessarily know which utc
        # offset the new Timestamp will wind up in (the tz for 1 month may be
        # different from 1 second) so we don't specify an expected_utc_offset
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets['utc_offset_standard']
            self._test_all_offsets(n=1, tstart=self._make_timestamp(
                self.ts_pre_fallback, hrs_pre, tz), expected_utc_offset=None)

    def test_springforward_singular(self):
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets['utc_offset_standard']
            self._test_all_offsets(n=1, tstart=self._make_timestamp(
                self.ts_pre_springfwd, hrs_pre, tz), expected_utc_offset=None)

    offset_classes = {MonthBegin: ['11/2/2012', '12/1/2012'],
                      MonthEnd: ['11/2/2012', '11/30/2012'],
                      BMonthBegin: ['11/2/2012', '12/3/2012'],
                      BMonthEnd: ['11/2/2012', '11/30/2012'],
                      CBMonthBegin: ['11/2/2012', '12/3/2012'],
                      CBMonthEnd: ['11/2/2012', '11/30/2012'],
                      SemiMonthBegin: ['11/2/2012', '11/15/2012'],
                      SemiMonthEnd: ['11/2/2012', '11/15/2012'],
                      Week: ['11/2/2012', '11/9/2012'],
                      YearBegin: ['11/2/2012', '1/1/2013'],
                      YearEnd: ['11/2/2012', '12/31/2012'],
                      BYearBegin: ['11/2/2012', '1/1/2013'],
                      BYearEnd: ['11/2/2012', '12/31/2012'],
                      QuarterBegin: ['11/2/2012', '12/1/2012'],
                      QuarterEnd: ['11/2/2012', '12/31/2012'],
                      BQuarterBegin: ['11/2/2012', '12/3/2012'],
                      BQuarterEnd: ['11/2/2012', '12/31/2012'],
                      Day: ['11/4/2012', '11/4/2012 23:00']}.items()

    @pytest.mark.parametrize('tup', offset_classes)
    def test_all_offset_classes(self, tup):
        offset, test_values = tup

        first = Timestamp(test_values[0], tz='US/Eastern') + offset()
        second = Timestamp(test_values[1], tz='US/Eastern')
        assert first == second


# ---------------------------------------------------------------------
def test_get_offset_day_error():
    # subclass of _BaseOffset must override _day_opt attribute, or we should
    # get a NotImplementedError

    with pytest.raises(NotImplementedError):
        DateOffset()._get_offset_day(datetime.now())


def test_valid_default_arguments(offset_types):
    # GH#19142 check that the calling the constructors without passing
    # any keyword arguments produce valid offsets
    cls = offset_types
    cls()


@pytest.mark.parametrize('kwd', sorted(list(liboffsets.relativedelta_kwds)))
def test_valid_month_attributes(kwd, month_classes):
    # GH#18226
    cls = month_classes
    # check that we cannot create e.g. MonthEnd(weeks=3)
    with pytest.raises(TypeError):
        cls(**{kwd: 3})


@pytest.mark.parametrize('kwd', sorted(list(liboffsets.relativedelta_kwds)))
def test_valid_tick_attributes(kwd, tick_classes):
    # GH#18226
    cls = tick_classes
    # check that we cannot create e.g. Hour(weeks=3)
    with pytest.raises(TypeError):
        cls(**{kwd: 3})


def test_validate_n_error():
    with pytest.raises(TypeError):
        DateOffset(n='Doh!')

    with pytest.raises(TypeError):
        MonthBegin(n=timedelta(1))

    with pytest.raises(TypeError):
        BDay(n=np.array([1, 2], dtype=np.int64))


def test_require_integers(offset_types):
    cls = offset_types
    with pytest.raises(ValueError):
        cls(n=1.5)


def test_weeks_onoffset():
    # GH#18510 Week with weekday = None, normalize = False should always
    # be onOffset
    offset = Week(n=2, weekday=None)
    ts = Timestamp('1862-01-13 09:03:34.873477378+0210', tz='Africa/Lusaka')
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow

    # negative n
    offset = Week(n=2, weekday=None)
    ts = Timestamp('1856-10-24 16:18:36.556360110-0717', tz='Pacific/Easter')
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow


def test_weekofmonth_onoffset():
    # GH#18864
    # Make sure that nanoseconds don't trip up onOffset (and with it apply)
    offset = WeekOfMonth(n=2, week=2, weekday=0)
    ts = Timestamp('1916-05-15 01:14:49.583410462+0422', tz='Asia/Qyzylorda')
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow

    # negative n
    offset = WeekOfMonth(n=-3, week=1, weekday=0)
    ts = Timestamp('1980-12-08 03:38:52.878321185+0500', tz='Asia/Oral')
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow


def test_last_week_of_month_on_offset():
    # GH#19036, GH#18977 _adjust_dst was incorrect for LastWeekOfMonth
    offset = LastWeekOfMonth(n=4, weekday=6)
    ts = Timestamp('1917-05-27 20:55:27.084284178+0200',
                   tz='Europe/Warsaw')
    slow = (ts + offset) - offset == ts
    fast = offset.onOffset(ts)
    assert fast == slow

    # negative n
    offset = LastWeekOfMonth(n=-4, weekday=5)
    ts = Timestamp('2005-08-27 05:01:42.799392561-0500',
                   tz='America/Rainy_River')
    slow = (ts + offset) - offset == ts
    fast = offset.onOffset(ts)
    assert fast == slow
