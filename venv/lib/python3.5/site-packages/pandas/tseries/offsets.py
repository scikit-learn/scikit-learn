# -*- coding: utf-8 -*-
from datetime import date, datetime, timedelta
import functools
import operator

from pandas.compat import range
from pandas import compat
import numpy as np

from pandas.core.dtypes.generic import ABCSeries, ABCDatetimeIndex, ABCPeriod
from pandas.core.tools.datetimes import to_datetime
import pandas.core.common as com

# import after tools, dateutil check
from dateutil.easter import easter
from pandas._libs import tslib, Timestamp, OutOfBoundsDatetime, Timedelta
from pandas.util._decorators import cache_readonly

from pandas._libs.tslibs import ccalendar, frequencies as libfrequencies
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds
import pandas._libs.tslibs.offsets as liboffsets
from pandas._libs.tslibs.offsets import (
    ApplyTypeError,
    as_datetime, _is_normalized,
    _get_calendar, _to_dt64,
    _determine_offset,
    apply_index_wraps,
    roll_yearday,
    shift_month,
    BaseOffset)


__all__ = ['Day', 'BusinessDay', 'BDay', 'CustomBusinessDay', 'CDay',
           'CBMonthEnd', 'CBMonthBegin',
           'MonthBegin', 'BMonthBegin', 'MonthEnd', 'BMonthEnd',
           'SemiMonthEnd', 'SemiMonthBegin',
           'BusinessHour', 'CustomBusinessHour',
           'YearBegin', 'BYearBegin', 'YearEnd', 'BYearEnd',
           'QuarterBegin', 'BQuarterBegin', 'QuarterEnd', 'BQuarterEnd',
           'LastWeekOfMonth', 'FY5253Quarter', 'FY5253',
           'Week', 'WeekOfMonth', 'Easter',
           'Hour', 'Minute', 'Second', 'Milli', 'Micro', 'Nano',
           'DateOffset']

# convert to/from datetime/timestamp to allow invalid Timestamp ranges to
# pass thru


def as_timestamp(obj):
    if isinstance(obj, Timestamp):
        return obj
    try:
        return Timestamp(obj)
    except (OutOfBoundsDatetime):
        pass
    return obj


def apply_wraps(func):
    @functools.wraps(func)
    def wrapper(self, other):
        if other is tslib.NaT:
            return tslib.NaT
        elif isinstance(other, (timedelta, Tick, DateOffset)):
            # timedelta path
            return func(self, other)
        elif isinstance(other, (np.datetime64, datetime, date)):
            other = as_timestamp(other)

        tz = getattr(other, 'tzinfo', None)
        nano = getattr(other, 'nanosecond', 0)

        try:
            if self._adjust_dst and isinstance(other, Timestamp):
                other = other.tz_localize(None)

            result = func(self, other)

            if self._adjust_dst:
                result = tslib._localize_pydatetime(result, tz)

            result = Timestamp(result)
            if self.normalize:
                result = result.normalize()

            # nanosecond may be deleted depending on offset process
            if not self.normalize and nano != 0:
                if not isinstance(self, Nano) and result.nanosecond != nano:
                    if result.tz is not None:
                        # convert to UTC
                        value = tslib.tz_convert_single(
                            result.value, 'UTC', result.tz)
                    else:
                        value = result.value
                    result = Timestamp(value + nano)

            if tz is not None and result.tzinfo is None:
                result = tslib._localize_pydatetime(result, tz)

        except OutOfBoundsDatetime:
            result = func(self, as_datetime(other))

            if self.normalize:
                # normalize_date returns normal datetime
                result = tslib.normalize_date(result)

            if tz is not None and result.tzinfo is None:
                result = tslib._localize_pydatetime(result, tz)

        return result
    return wrapper


def shift_day(other, days):
    """
    Increment the datetime `other` by the given number of days, retaining
    the time-portion of the datetime.  For tz-naive datetimes this is
    equivalent to adding a timedelta.  For tz-aware datetimes it is similar to
    dateutil's relativedelta.__add__, but handles pytz tzinfo objects.

    Parameters
    ----------
    other : datetime or Timestamp
    days : int

    Returns
    -------
    shifted: datetime or Timestamp
    """
    if other.tzinfo is None:
        return other + timedelta(days=days)

    tz = other.tzinfo
    naive = other.replace(tzinfo=None)
    shifted = naive + timedelta(days=days)
    return tslib._localize_pydatetime(shifted, tz)


# ---------------------------------------------------------------------
# DateOffset


class DateOffset(BaseOffset):
    """
    Standard kind of date increment used for a date range.

    Works exactly like relativedelta in terms of the keyword args you
    pass in, use of the keyword n is discouraged-- you would be better
    off specifying n in the keywords you use, but regardless it is
    there for you. n is needed for DateOffset subclasses.

    DateOffets work as follows.  Each offset specify a set of dates
    that conform to the DateOffset.  For example, Bday defines this
    set to be the set of dates that are weekdays (M-F).  To test if a
    date is in the set of a DateOffset dateOffset we can use the
    onOffset method: dateOffset.onOffset(date).

    If a date is not on a valid date, the rollback and rollforward
    methods can be used to roll the date to the nearest valid date
    before/after the date.

    DateOffsets can be created to move dates forward a given number of
    valid dates.  For example, Bday(2) can be added to a date to move
    it two business days forward.  If the date does not start on a
    valid date, first it is moved to a valid date.  Thus pseudo code
    is:

    def __add__(date):
      date = rollback(date) # does nothing if date is valid
      return date + <n number of periods>

    When a date offset is created for a negative number of periods,
    the date is first rolled forward.  The pseudo code is:

    def __add__(date):
      date = rollforward(date) # does nothing is date is valid
      return date + <n number of periods>

    Zero presents a problem.  Should it roll forward or back?  We
    arbitrarily have it rollforward:

    date + BDay(0) == BDay.rollforward(date)

    Since 0 is a bit weird, we suggest avoiding its use.
    """
    _use_relativedelta = False
    _adjust_dst = False
    _attributes = frozenset(['n', 'normalize'] +
                            list(liboffsets.relativedelta_kwds))

    # default for prior pickles
    normalize = False

    def __init__(self, n=1, normalize=False, **kwds):
        self.n = self._validate_n(n)
        self.normalize = normalize

        self._offset, self._use_relativedelta = _determine_offset(kwds)
        self.__dict__.update(kwds)

    @apply_wraps
    def apply(self, other):
        if self._use_relativedelta:
            other = as_datetime(other)

        if len(self.kwds) > 0:
            tzinfo = getattr(other, 'tzinfo', None)
            if tzinfo is not None and self._use_relativedelta:
                # perform calculation in UTC
                other = other.replace(tzinfo=None)

            if self.n > 0:
                for i in range(self.n):
                    other = other + self._offset
            else:
                for i in range(-self.n):
                    other = other - self._offset

            if tzinfo is not None and self._use_relativedelta:
                # bring tz back from UTC calculation
                other = tslib._localize_pydatetime(other, tzinfo)

            return as_timestamp(other)
        else:
            return other + timedelta(self.n)

    @apply_index_wraps
    def apply_index(self, i):
        """
        Vectorized apply of DateOffset to DatetimeIndex,
        raises NotImplentedError for offsets without a
        vectorized implementation

        Parameters
        ----------
        i : DatetimeIndex

        Returns
        -------
        y : DatetimeIndex
        """

        if type(self) is not DateOffset:
            raise NotImplementedError("DateOffset subclass {name} "
                                      "does not have a vectorized "
                                      "implementation".format(
                                          name=self.__class__.__name__))
        kwds = self.kwds
        relativedelta_fast = set(['years', 'months', 'weeks',
                                  'days', 'hours', 'minutes',
                                  'seconds', 'microseconds'])
        # relativedelta/_offset path only valid for base DateOffset
        if (self._use_relativedelta and
                set(kwds).issubset(relativedelta_fast)):

            months = ((kwds.get('years', 0) * 12 +
                       kwds.get('months', 0)) * self.n)
            if months:
                shifted = liboffsets.shift_months(i.asi8, months)
                i = i._shallow_copy(shifted)

            weeks = (kwds.get('weeks', 0)) * self.n
            if weeks:
                i = (i.to_period('W') + weeks).to_timestamp() + \
                    i.to_perioddelta('W')

            timedelta_kwds = {k: v for k, v in kwds.items()
                              if k in ['days', 'hours', 'minutes',
                                       'seconds', 'microseconds']}
            if timedelta_kwds:
                delta = Timedelta(**timedelta_kwds)
                i = i + (self.n * delta)
            return i
        elif not self._use_relativedelta and hasattr(self, '_offset'):
            # timedelta
            return i + (self._offset * self.n)
        else:
            # relativedelta with other keywords
            kwd = set(kwds) - relativedelta_fast
            raise NotImplementedError("DateOffset with relativedelta "
                                      "keyword(s) {kwd} not able to be "
                                      "applied vectorized".format(kwd=kwd))

    def isAnchored(self):
        # TODO: Does this make sense for the general case?  It would help
        # if there were a canonical docstring for what isAnchored means.
        return (self.n == 1)

    def _params(self):
        all_paras = self.__dict__.copy()
        if 'holidays' in all_paras and not all_paras['holidays']:
            all_paras.pop('holidays')
        exclude = ['kwds', 'name', 'normalize', 'calendar']
        attrs = [(k, v) for k, v in all_paras.items()
                 if (k not in exclude) and (k[0] != '_')]
        attrs = sorted(set(attrs))
        params = tuple([str(self.__class__)] + attrs)
        return params

    # TODO: Combine this with BusinessMixin version by defining a whitelisted
    # set of attributes on each object rather than the existing behavior of
    # iterating over internal ``__dict__``
    def _repr_attrs(self):
        exclude = set(['n', 'inc', 'normalize'])
        attrs = []
        for attr in sorted(self.__dict__):
            if attr.startswith('_') or attr == 'kwds':
                continue
            elif attr not in exclude:
                value = getattr(self, attr)
                attrs.append('{attr}={value}'.format(attr=attr, value=value))

        out = ''
        if attrs:
            out += ': ' + ', '.join(attrs)
        return out

    @property
    def name(self):
        return self.rule_code

    def __eq__(self, other):
        if other is None:
            return False

        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if not isinstance(other, DateOffset):
            return False

        return self._params() == other._params()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._params())

    def __add__(self, other):
        if isinstance(other, (ABCDatetimeIndex, ABCSeries)):
            return other + self
        elif isinstance(other, ABCPeriod):
            return other + self
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, datetime):
            raise TypeError('Cannot subtract datetime from offset.')
        elif type(other) == type(self):
            return self.__class__(self.n - other.n, normalize=self.normalize,
                                  **self.kwds)
        else:  # pragma: no cover
            return NotImplemented

    def rollback(self, dt):
        """Roll provided date backward to next offset only if not on offset"""
        dt = as_timestamp(dt)
        if not self.onOffset(dt):
            dt = dt - self.__class__(1, normalize=self.normalize, **self.kwds)
        return dt

    def rollforward(self, dt):
        """Roll provided date forward to next offset only if not on offset"""
        dt = as_timestamp(dt)
        if not self.onOffset(dt):
            dt = dt + self.__class__(1, normalize=self.normalize, **self.kwds)
        return dt

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        # XXX, see #1395
        if type(self) == DateOffset or isinstance(self, Tick):
            return True

        # Default (slow) method for determining if some date is a member of the
        # date range generated by this offset. Subclasses may have this
        # re-implemented in a nicer way.
        a = dt
        b = ((dt + self) - self)
        return a == b

    # way to get around weirdness with rule_code
    @property
    def _prefix(self):
        raise NotImplementedError('Prefix not defined')

    @property
    def rule_code(self):
        return self._prefix

    @property
    def freqstr(self):
        try:
            code = self.rule_code
        except NotImplementedError:
            return repr(self)

        if self.n != 1:
            fstr = '{n}{code}'.format(n=self.n, code=code)
        else:
            fstr = code

        try:
            if self._offset:
                fstr += self._offset_str()
        except AttributeError:
            # TODO: standardize `_offset` vs `offset` naming convention
            pass

        return fstr

    def _offset_str(self):
        return ''

    @property
    def nanos(self):
        raise ValueError("{name} is a non-fixed frequency".format(name=self))

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        if 'offset' in state:
            # Older (<0.22.0) versions have offset attribute instead of _offset
            if '_offset' in state:  # pragma: no cover
                raise AssertionError('Unexpected key `_offset`')
            state['_offset'] = state.pop('offset')
            state['kwds']['offset'] = state['_offset']

        if '_offset' in state and not isinstance(state['_offset'], timedelta):
            # relativedelta, we need to populate using its kwds
            offset = state['_offset']
            odict = offset.__dict__
            kwds = {key: odict[key] for key in odict if odict[key]}
            state.update(kwds)

        self.__dict__ = state
        if 'weekmask' in state and 'holidays' in state:
            calendar, holidays = _get_calendar(weekmask=self.weekmask,
                                               holidays=self.holidays,
                                               calendar=None)
            self.calendar = calendar
            self.holidays = holidays


class SingleConstructorOffset(DateOffset):
    @classmethod
    def _from_name(cls, suffix=None):
        # default _from_name calls cls with no args
        if suffix:
            raise ValueError("Bad freq suffix {suffix}".format(suffix=suffix))
        return cls()


class _CustomMixin(object):
    """
    Mixin for classes that define and validate calendar, holidays,
    and weekdays attributes
    """
    def __init__(self, weekmask, holidays, calendar):
        calendar, holidays = _get_calendar(weekmask=weekmask,
                                           holidays=holidays,
                                           calendar=calendar)
        # Custom offset instances are identified by the
        # following two attributes. See DateOffset._params()
        # holidays, weekmask

        self.weekmask = weekmask
        self.holidays = holidays
        self.calendar = calendar


class BusinessMixin(object):
    """ Mixin to business types to provide related functions """

    @property
    def offset(self):
        """Alias for self._offset"""
        # Alias for backward compat
        return self._offset

    def _repr_attrs(self):
        if self.offset:
            attrs = ['offset={offset!r}'.format(offset=self.offset)]
        else:
            attrs = None
        out = ''
        if attrs:
            out += ': ' + ', '.join(attrs)
        return out

    def __getstate__(self):
        """Return a pickleable state"""
        state = self.__dict__.copy()

        # we don't want to actually pickle the calendar object
        # as its a np.busyday; we recreate on deserilization
        if 'calendar' in state:
            del state['calendar']
        try:
            state['kwds'].pop('calendar')
        except KeyError:
            pass

        return state


class BusinessDay(BusinessMixin, SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n business days
    """
    _prefix = 'B'
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize', 'offset'])

    def __init__(self, n=1, normalize=False, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset

    def _offset_str(self):
        def get_str(td):
            off_str = ''
            if td.days > 0:
                off_str += str(td.days) + 'D'
            if td.seconds > 0:
                s = td.seconds
                hrs = int(s / 3600)
                if hrs != 0:
                    off_str += str(hrs) + 'H'
                    s -= hrs * 3600
                mts = int(s / 60)
                if mts != 0:
                    off_str += str(mts) + 'Min'
                    s -= mts * 60
                if s != 0:
                    off_str += str(s) + 's'
            if td.microseconds > 0:
                off_str += str(td.microseconds) + 'us'
            return off_str

        if isinstance(self.offset, timedelta):
            zero = timedelta(0, 0, 0)
            if self.offset >= zero:
                off_str = '+' + get_str(self.offset)
            else:
                off_str = '-' + get_str(-self.offset)
            return off_str
        else:
            return '+' + repr(self.offset)

    @apply_wraps
    def apply(self, other):
        if isinstance(other, datetime):
            n = self.n
            wday = other.weekday()

            # avoid slowness below by operating on weeks first
            weeks = n // 5
            if n <= 0 and wday > 4:
                # roll forward
                n += 1

            n -= 5 * weeks

            # n is always >= 0 at this point
            if n == 0 and wday > 4:
                # roll back
                days = 4 - wday
            elif wday > 4:
                # roll forward
                days = (7 - wday) + (n - 1)
            elif wday + n <= 4:
                # shift by n days without leaving the current week
                days = n
            else:
                # shift by n days plus 2 to get past the weekend
                days = n + 2

            result = other + timedelta(days=7 * weeks + days)
            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise ApplyTypeError('Only know how to combine business day with '
                                 'datetime or timedelta.')

    @apply_index_wraps
    def apply_index(self, i):
        time = i.to_perioddelta('D')
        # to_period rolls forward to next BDay; track and
        # reduce n where it does when rolling forward
        shifted = (i.to_perioddelta('B') - time).asi8 != 0
        if self.n > 0:
            roll = np.where(shifted, self.n - 1, self.n)
        else:
            roll = self.n

        return (i.to_period('B') + roll).to_timestamp() + time

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.weekday() < 5


class BusinessHourMixin(BusinessMixin):

    def __init__(self, start='09:00', end='17:00', offset=timedelta(0)):
        # must be validated here to equality check
        self.start = liboffsets._validate_business_time(start)
        self.end = liboffsets._validate_business_time(end)
        self._offset = offset

    @cache_readonly
    def next_bday(self):
        """used for moving to next businessday"""
        if self.n >= 0:
            nb_offset = 1
        else:
            nb_offset = -1
        if self._prefix.startswith('C'):
            # CustomBusinessHour
            return CustomBusinessDay(n=nb_offset,
                                     weekmask=self.weekmask,
                                     holidays=self.holidays,
                                     calendar=self.calendar)
        else:
            return BusinessDay(n=nb_offset)

    # TODO: Cache this once offsets are immutable
    def _get_daytime_flag(self):
        if self.start == self.end:
            raise ValueError('start and end must not be the same')
        elif self.start < self.end:
            return True
        else:
            return False

    def _next_opening_time(self, other):
        """
        If n is positive, return tomorrow's business day opening time.
        Otherwise yesterday's business day's opening time.

        Opening time always locates on BusinessDay.
        Otherwise, closing time may not if business hour extends over midnight.
        """
        if not self.next_bday.onOffset(other):
            other = other + self.next_bday
        else:
            if self.n >= 0 and self.start < other.time():
                other = other + self.next_bday
            elif self.n < 0 and other.time() < self.start:
                other = other + self.next_bday
        return datetime(other.year, other.month, other.day,
                        self.start.hour, self.start.minute)

    def _prev_opening_time(self, other):
        """
        If n is positive, return yesterday's business day opening time.
        Otherwise yesterday business day's opening time.
        """
        if not self.next_bday.onOffset(other):
            other = other - self.next_bday
        else:
            if self.n >= 0 and other.time() < self.start:
                other = other - self.next_bday
            elif self.n < 0 and other.time() > self.start:
                other = other - self.next_bday
        return datetime(other.year, other.month, other.day,
                        self.start.hour, self.start.minute)

    # TODO: cache this once offsets are immutable
    def _get_business_hours_by_sec(self):
        """
        Return business hours in a day by seconds.
        """
        if self._get_daytime_flag():
            # create dummy datetime to calculate businesshours in a day
            dtstart = datetime(2014, 4, 1, self.start.hour, self.start.minute)
            until = datetime(2014, 4, 1, self.end.hour, self.end.minute)
            return (until - dtstart).total_seconds()
        else:
            self.daytime = False
            dtstart = datetime(2014, 4, 1, self.start.hour, self.start.minute)
            until = datetime(2014, 4, 2, self.end.hour, self.end.minute)
            return (until - dtstart).total_seconds()

    @apply_wraps
    def rollback(self, dt):
        """Roll provided date backward to next offset only if not on offset"""
        if not self.onOffset(dt):
            businesshours = self._get_business_hours_by_sec()
            if self.n >= 0:
                dt = self._prev_opening_time(
                    dt) + timedelta(seconds=businesshours)
            else:
                dt = self._next_opening_time(
                    dt) + timedelta(seconds=businesshours)
        return dt

    @apply_wraps
    def rollforward(self, dt):
        """Roll provided date forward to next offset only if not on offset"""
        if not self.onOffset(dt):
            if self.n >= 0:
                return self._next_opening_time(dt)
            else:
                return self._prev_opening_time(dt)
        return dt

    @apply_wraps
    def apply(self, other):
        # calculate here because offset is not immutable
        daytime = self._get_daytime_flag()
        businesshours = self._get_business_hours_by_sec()
        bhdelta = timedelta(seconds=businesshours)

        if isinstance(other, datetime):
            # used for detecting edge condition
            nanosecond = getattr(other, 'nanosecond', 0)
            # reset timezone and nanosecond
            # other may be a Timestamp, thus not use replace
            other = datetime(other.year, other.month, other.day,
                             other.hour, other.minute,
                             other.second, other.microsecond)
            n = self.n
            if n >= 0:
                if (other.time() == self.end or
                        not self._onOffset(other, businesshours)):
                    other = self._next_opening_time(other)
            else:
                if other.time() == self.start:
                    # adjustment to move to previous business day
                    other = other - timedelta(seconds=1)
                if not self._onOffset(other, businesshours):
                    other = self._next_opening_time(other)
                    other = other + bhdelta

            bd, r = divmod(abs(n * 60), businesshours // 60)
            if n < 0:
                bd, r = -bd, -r

            if bd != 0:
                skip_bd = BusinessDay(n=bd)
                # midnight business hour may not on BusinessDay
                if not self.next_bday.onOffset(other):
                    remain = other - self._prev_opening_time(other)
                    other = self._next_opening_time(other + skip_bd) + remain
                else:
                    other = other + skip_bd

            hours, minutes = divmod(r, 60)
            result = other + timedelta(hours=hours, minutes=minutes)

            # because of previous adjustment, time will be larger than start
            if ((daytime and (result.time() < self.start or
                              self.end < result.time())) or
                    not daytime and (self.end < result.time() < self.start)):
                if n >= 0:
                    bday_edge = self._prev_opening_time(other)
                    bday_edge = bday_edge + bhdelta
                    # calculate remainder
                    bday_remain = result - bday_edge
                    result = self._next_opening_time(other)
                    result += bday_remain
                else:
                    bday_edge = self._next_opening_time(other)
                    bday_remain = result - bday_edge
                    result = self._next_opening_time(result) + bhdelta
                    result += bday_remain
            # edge handling
            if n >= 0:
                if result.time() == self.end:
                    result = self._next_opening_time(result)
            else:
                if result.time() == self.start and nanosecond == 0:
                    # adjustment to move to previous business day
                    result = self._next_opening_time(
                        result - timedelta(seconds=1)) + bhdelta

            return result
        else:
            # TODO: Figure out the end of this sente
            raise ApplyTypeError(
                'Only know how to combine business hour with ')

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False

        if dt.tzinfo is not None:
            dt = datetime(dt.year, dt.month, dt.day, dt.hour,
                          dt.minute, dt.second, dt.microsecond)
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        businesshours = self._get_business_hours_by_sec()
        return self._onOffset(dt, businesshours)

    def _onOffset(self, dt, businesshours):
        """
        Slight speedups using calculated values
        """
        # if self.normalize and not _is_normalized(dt):
        #     return False
        # Valid BH can be on the different BusinessDay during midnight
        # Distinguish by the time spent from previous opening time
        if self.n >= 0:
            op = self._prev_opening_time(dt)
        else:
            op = self._next_opening_time(dt)
        span = (dt - op).total_seconds()
        if span <= businesshours:
            return True
        else:
            return False

    def _repr_attrs(self):
        out = super(BusinessHourMixin, self)._repr_attrs()
        start = self.start.strftime('%H:%M')
        end = self.end.strftime('%H:%M')
        attrs = ['{prefix}={start}-{end}'.format(prefix=self._prefix,
                                                 start=start, end=end)]
        out += ': ' + ', '.join(attrs)
        return out


class BusinessHour(BusinessHourMixin, SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n business days

    .. versionadded:: 0.16.1

    """
    _prefix = 'BH'
    _anchor = 0
    _attributes = frozenset(['n', 'normalize', 'start', 'end', 'offset'])

    def __init__(self, n=1, normalize=False, start='09:00',
                 end='17:00', offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        super(BusinessHour, self).__init__(start=start, end=end, offset=offset)


class CustomBusinessDay(_CustomMixin, BusinessDay):
    """
    DateOffset subclass representing possibly n custom business days,
    excluding holidays

    Parameters
    ----------
    n : int, default 1
    offset : timedelta, default timedelta(0)
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    calendar : pd.HolidayCalendar or np.busdaycalendar
    """
    _cacheable = False
    _prefix = 'C'
    _attributes = frozenset(['n', 'normalize',
                             'weekmask', 'holidays', 'calendar', 'offset'])

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset

        _CustomMixin.__init__(self, weekmask, holidays, calendar)

    @apply_wraps
    def apply(self, other):
        if self.n <= 0:
            roll = 'forward'
        else:
            roll = 'backward'

        if isinstance(other, datetime):
            date_in = other
            np_dt = np.datetime64(date_in.date())

            np_incr_dt = np.busday_offset(np_dt, self.n, roll=roll,
                                          busdaycal=self.calendar)

            dt_date = np_incr_dt.astype(datetime)
            result = datetime.combine(dt_date, date_in.time())

            if self.offset:
                result = result + self.offset
            return result

        elif isinstance(other, (timedelta, Tick)):
            return BDay(self.n, offset=self.offset + other,
                        normalize=self.normalize)
        else:
            raise ApplyTypeError('Only know how to combine trading day with '
                                 'datetime, datetime64 or timedelta.')

    def apply_index(self, i):
        raise NotImplementedError

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        day64 = _to_dt64(dt, 'datetime64[D]')
        return np.is_busday(day64, busdaycal=self.calendar)


class CustomBusinessHour(_CustomMixin, BusinessHourMixin,
                         SingleConstructorOffset):
    """
    DateOffset subclass representing possibly n custom business days

    .. versionadded:: 0.18.1

    """
    _prefix = 'CBH'
    _anchor = 0
    _attributes = frozenset(['n', 'normalize',
                             'weekmask', 'holidays', 'calendar',
                             'start', 'end', 'offset'])

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None,
                 start='09:00', end='17:00', offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset

        _CustomMixin.__init__(self, weekmask, holidays, calendar)
        BusinessHourMixin.__init__(self, start=start, end=end, offset=offset)


# ---------------------------------------------------------------------
# Month-Based Offset Classes


class MonthOffset(SingleConstructorOffset):
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize'])

    def __init__(self, n=1, normalize=False):
        self.n = self._validate_n(n)
        self.normalize = normalize

    @property
    def name(self):
        if self.isAnchored:
            return self.rule_code
        else:
            month = ccalendar.MONTH_ALIASES[self.n]
            return "{code}-{month}".format(code=self.rule_code,
                                           month=month)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.day == self._get_offset_day(dt)

    @apply_wraps
    def apply(self, other):
        compare_day = self._get_offset_day(other)
        n = liboffsets.roll_convention(other.day, self.n, compare_day)
        return shift_month(other, n, self._day_opt)

    @apply_index_wraps
    def apply_index(self, i):
        shifted = liboffsets.shift_months(i.asi8, self.n, self._day_opt)
        return i._shallow_copy(shifted)


class MonthEnd(MonthOffset):
    """DateOffset of one month end"""
    _prefix = 'M'
    _day_opt = 'end'


class MonthBegin(MonthOffset):
    """DateOffset of one month at beginning"""
    _prefix = 'MS'
    _day_opt = 'start'


class BusinessMonthEnd(MonthOffset):
    """DateOffset increments between business EOM dates"""
    _prefix = 'BM'
    _day_opt = 'business_end'


class BusinessMonthBegin(MonthOffset):
    """DateOffset of one business month at beginning"""
    _prefix = 'BMS'
    _day_opt = 'business_start'


class _CustomBusinessMonth(_CustomMixin, BusinessMixin, MonthOffset):
    """
    DateOffset subclass representing one custom business month, incrementing
    between [BEGIN/END] of month dates

    Parameters
    ----------
    n : int, default 1
    offset : timedelta, default timedelta(0)
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    weekmask : str, Default 'Mon Tue Wed Thu Fri'
        weekmask of valid business days, passed to ``numpy.busdaycalendar``
    holidays : list
        list/array of dates to exclude from the set of valid business days,
        passed to ``numpy.busdaycalendar``
    calendar : pd.HolidayCalendar or np.busdaycalendar
    """
    _cacheable = False
    _attributes = frozenset(['n', 'normalize',
                             'weekmask', 'holidays', 'calendar', 'offset'])

    onOffset = DateOffset.onOffset        # override MonthOffset method
    apply_index = DateOffset.apply_index  # override MonthOffset method

    def __init__(self, n=1, normalize=False, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, calendar=None, offset=timedelta(0)):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self._offset = offset

        _CustomMixin.__init__(self, weekmask, holidays, calendar)

    @cache_readonly
    def cbday_roll(self):
        """Define default roll function to be called in apply method"""
        cbday = CustomBusinessDay(n=self.n, normalize=False, **self.kwds)

        if self._prefix.endswith('S'):
            # MonthBegin
            roll_func = cbday.rollforward
        else:
            # MonthEnd
            roll_func = cbday.rollback
        return roll_func

    @cache_readonly
    def m_offset(self):
        if self._prefix.endswith('S'):
            # MonthBegin
            moff = MonthBegin(n=1, normalize=False)
        else:
            # MonthEnd
            moff = MonthEnd(n=1, normalize=False)
        return moff

    @cache_readonly
    def month_roll(self):
        """Define default roll function to be called in apply method"""
        if self._prefix.endswith('S'):
            # MonthBegin
            roll_func = self.m_offset.rollback
        else:
            # MonthEnd
            roll_func = self.m_offset.rollforward
        return roll_func

    @apply_wraps
    def apply(self, other):
        # First move to month offset
        cur_month_offset_date = self.month_roll(other)

        # Find this custom month offset
        compare_date = self.cbday_roll(cur_month_offset_date)
        n = liboffsets.roll_convention(other.day, self.n, compare_date.day)

        new = cur_month_offset_date + n * self.m_offset
        result = self.cbday_roll(new)
        return result


class CustomBusinessMonthEnd(_CustomBusinessMonth):
    __doc__ = _CustomBusinessMonth.__doc__.replace('[BEGIN/END]', 'end')
    _prefix = 'CBM'


class CustomBusinessMonthBegin(_CustomBusinessMonth):
    __doc__ = _CustomBusinessMonth.__doc__.replace('[BEGIN/END]', 'beginning')
    _prefix = 'CBMS'


# ---------------------------------------------------------------------
# Semi-Month Based Offset Classes

class SemiMonthOffset(DateOffset):
    _adjust_dst = True
    _default_day_of_month = 15
    _min_day_of_month = 2
    _attributes = frozenset(['n', 'normalize', 'day_of_month'])

    def __init__(self, n=1, normalize=False, day_of_month=None):
        if day_of_month is None:
            self.day_of_month = self._default_day_of_month
        else:
            self.day_of_month = int(day_of_month)
        if not self._min_day_of_month <= self.day_of_month <= 27:
            msg = 'day_of_month must be {min}<=day_of_month<=27, got {day}'
            raise ValueError(msg.format(min=self._min_day_of_month,
                                        day=self.day_of_month))

        self.n = self._validate_n(n)
        self.normalize = normalize

    @classmethod
    def _from_name(cls, suffix=None):
        return cls(day_of_month=suffix)

    @property
    def rule_code(self):
        suffix = '-{day_of_month}'.format(day_of_month=self.day_of_month)
        return self._prefix + suffix

    @apply_wraps
    def apply(self, other):
        # shift `other` to self.day_of_month, incrementing `n` if necessary
        n = liboffsets.roll_convention(other.day, self.n, self.day_of_month)

        days_in_month = tslib.monthrange(other.year, other.month)[1]

        # For SemiMonthBegin on other.day == 1 and
        # SemiMonthEnd on other.day == days_in_month,
        # shifting `other` to `self.day_of_month` _always_ requires
        # incrementing/decrementing `n`, regardless of whether it is
        # initially positive.
        if type(self) is SemiMonthBegin and (self.n <= 0 and other.day == 1):
            n -= 1
        elif type(self) is SemiMonthEnd and (self.n > 0 and
                                             other.day == days_in_month):
            n += 1

        return self._apply(n, other)

    def _apply(self, n, other):
        """Handle specific apply logic for child classes"""
        raise com.AbstractMethodError(self)

    @apply_index_wraps
    def apply_index(self, i):
        # determine how many days away from the 1st of the month we are
        days_from_start = i.to_perioddelta('M').asi8
        delta = Timedelta(days=self.day_of_month - 1).value

        # get boolean array for each element before the day_of_month
        before_day_of_month = days_from_start < delta

        # get boolean array for each element after the day_of_month
        after_day_of_month = days_from_start > delta

        # determine the correct n for each date in i
        roll = self._get_roll(i, before_day_of_month, after_day_of_month)

        # isolate the time since it will be striped away one the next line
        time = i.to_perioddelta('D')

        # apply the correct number of months
        i = (i.to_period('M') + (roll // 2)).to_timestamp()

        # apply the correct day
        i = self._apply_index_days(i, roll)

        return i + time

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        """Return an array with the correct n for each date in i.

        The roll array is based on the fact that i gets rolled back to
        the first day of the month.
        """
        raise com.AbstractMethodError(self)

    def _apply_index_days(self, i, roll):
        """Apply the correct day for each date in i"""
        raise com.AbstractMethodError(self)


class SemiMonthEnd(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the last
    day of the month and day_of_month.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    n: int
    normalize : bool, default False
    day_of_month: int, {1, 3,...,27}, default 15
    """
    _prefix = 'SM'
    _min_day_of_month = 1

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        _, days_in_month = tslib.monthrange(dt.year, dt.month)
        return dt.day in (self.day_of_month, days_in_month)

    def _apply(self, n, other):
        months = n // 2
        day = 31 if n % 2 else self.day_of_month
        return shift_month(other, months, day)

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        n = self.n
        is_month_end = i.is_month_end
        if n > 0:
            roll_end = np.where(is_month_end, 1, 0)
            roll_before = np.where(before_day_of_month, n, n + 1)
            roll = roll_end + roll_before
        elif n == 0:
            roll_after = np.where(after_day_of_month, 2, 0)
            roll_before = np.where(~after_day_of_month, 1, 0)
            roll = roll_before + roll_after
        else:
            roll = np.where(after_day_of_month, n + 2, n + 1)
        return roll

    def _apply_index_days(self, i, roll):
        """Add days portion of offset to DatetimeIndex i

        Parameters
        ----------
        i : DatetimeIndex
        roll : ndarray[int64_t]

        Returns
        -------
        result : DatetimeIndex
        """
        nanos = (roll % 2) * Timedelta(days=self.day_of_month).value
        i += nanos.astype('timedelta64[ns]')
        return i + Timedelta(days=-1)


class SemiMonthBegin(SemiMonthOffset):
    """
    Two DateOffset's per month repeating on the first
    day of the month and day_of_month.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    n: int
    normalize : bool, default False
    day_of_month: int, {2, 3,...,27}, default 15
    """
    _prefix = 'SMS'

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.day in (1, self.day_of_month)

    def _apply(self, n, other):
        months = n // 2 + n % 2
        day = 1 if n % 2 else self.day_of_month
        return shift_month(other, months, day)

    def _get_roll(self, i, before_day_of_month, after_day_of_month):
        n = self.n
        is_month_start = i.is_month_start
        if n > 0:
            roll = np.where(before_day_of_month, n, n + 1)
        elif n == 0:
            roll_start = np.where(is_month_start, 0, 1)
            roll_after = np.where(after_day_of_month, 1, 0)
            roll = roll_start + roll_after
        else:
            roll_after = np.where(after_day_of_month, n + 2, n + 1)
            roll_start = np.where(is_month_start, -1, 0)
            roll = roll_after + roll_start
        return roll

    def _apply_index_days(self, i, roll):
        """Add days portion of offset to DatetimeIndex i

        Parameters
        ----------
        i : DatetimeIndex
        roll : ndarray[int64_t]

        Returns
        -------
        result : DatetimeIndex
        """
        nanos = (roll % 2) * Timedelta(days=self.day_of_month - 1).value
        return i + nanos.astype('timedelta64[ns]')


# ---------------------------------------------------------------------
# Week-Based Offset Classes

class Week(DateOffset):
    """
    Weekly offset

    Parameters
    ----------
    weekday : int, default None
        Always generate specific day of week. 0 for Monday
    """
    _adjust_dst = True
    _inc = timedelta(weeks=1)
    _prefix = 'W'
    _attributes = frozenset(['n', 'normalize', 'weekday'])

    def __init__(self, n=1, normalize=False, weekday=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday

        if self.weekday is not None:
            if self.weekday < 0 or self.weekday > 6:
                raise ValueError('Day must be 0<=day<=6, got {day}'
                                 .format(day=self.weekday))

    def isAnchored(self):
        return (self.n == 1 and self.weekday is not None)

    @apply_wraps
    def apply(self, other):
        if self.weekday is None:
            return other + self.n * self._inc

        k = self.n
        otherDay = other.weekday()
        if otherDay != self.weekday:
            other = other + timedelta((self.weekday - otherDay) % 7)
            if k > 0:
                k -= 1

        return other + timedelta(weeks=k)

    @apply_index_wraps
    def apply_index(self, i):
        if self.weekday is None:
            return ((i.to_period('W') + self.n).to_timestamp() +
                    i.to_perioddelta('W'))
        else:
            return self._end_apply_index(i)

    def _end_apply_index(self, dtindex):
        """Add self to the given DatetimeIndex, specialized for case where
        self.weekday is non-null.

        Parameters
        ----------
        dtindex : DatetimeIndex

        Returns
        -------
        result : DatetimeIndex
        """
        off = dtindex.to_perioddelta('D')

        base, mult = libfrequencies.get_freq_code(self.freqstr)
        base_period = dtindex.to_period(base)
        if self.n > 0:
            # when adding, dates on end roll to next
            normed = dtindex - off
            roll = np.where(base_period.to_timestamp(how='end') == normed,
                            self.n, self.n - 1)
        else:
            roll = self.n

        base = (base_period + roll).to_timestamp(how='end')
        return base + off

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        elif self.weekday is None:
            return True
        return dt.weekday() == self.weekday

    @property
    def rule_code(self):
        suffix = ''
        if self.weekday is not None:
            weekday = ccalendar.int_to_weekday[self.weekday]
            suffix = '-{weekday}'.format(weekday=weekday)
        return self._prefix + suffix

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            weekday = None
        else:
            weekday = ccalendar.weekday_to_int[suffix]
        return cls(weekday=weekday)


class _WeekOfMonthMixin(object):
    """Mixin for methods common to WeekOfMonth and LastWeekOfMonth"""
    @apply_wraps
    def apply(self, other):
        compare_day = self._get_offset_day(other)

        months = self.n
        if months > 0 and compare_day > other.day:
            months -= 1
        elif months <= 0 and compare_day < other.day:
            months += 1

        shifted = shift_month(other, months, 'start')
        to_day = self._get_offset_day(shifted)
        return shift_day(shifted, to_day - shifted.day)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.day == self._get_offset_day(dt)


class WeekOfMonth(_WeekOfMonthMixin, DateOffset):
    """
    Describes monthly dates like "the Tuesday of the 2nd week of each month"

    Parameters
    ----------
    n : int
    week : {0, 1, 2, 3, ...}, default 0
        0 is 1st week of month, 1 2nd week, etc.
    weekday : {0, 1, ..., 6}, default 0
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    """
    _prefix = 'WOM'
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize', 'week', 'weekday'])

    def __init__(self, n=1, normalize=False, week=0, weekday=0):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday
        self.week = week

        if self.weekday < 0 or self.weekday > 6:
            raise ValueError('Day must be 0<=day<=6, got {day}'
                             .format(day=self.weekday))
        if self.week < 0 or self.week > 3:
            raise ValueError('Week must be 0<=week<=3, got {week}'
                             .format(week=self.week))

    def _get_offset_day(self, other):
        """
        Find the day in the same month as other that has the same
        weekday as self.weekday and is the self.week'th such day in the month.

        Parameters
        ----------
        other: datetime

        Returns
        -------
        day: int
        """
        mstart = datetime(other.year, other.month, 1)
        wday = mstart.weekday()
        shift_days = (self.weekday - wday) % 7
        return 1 + shift_days + self.week * 7

    @property
    def rule_code(self):
        weekday = ccalendar.int_to_weekday.get(self.weekday, '')
        return '{prefix}-{week}{weekday}'.format(prefix=self._prefix,
                                                 week=self.week + 1,
                                                 weekday=weekday)

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError("Prefix {prefix!r} requires a suffix."
                             .format(prefix=cls._prefix))
        # TODO: handle n here...
        # only one digit weeks (1 --> week 0, 2 --> week 1, etc.)
        week = int(suffix[0]) - 1
        weekday = ccalendar.weekday_to_int[suffix[1:]]
        return cls(week=week, weekday=weekday)


class LastWeekOfMonth(_WeekOfMonthMixin, DateOffset):
    """
    Describes monthly dates in last week of month like "the last Tuesday of
    each month"

    Parameters
    ----------
    n : int, default 1
    weekday : {0, 1, ..., 6}, default 0
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays

    """
    _prefix = 'LWOM'
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize', 'weekday'])

    def __init__(self, n=1, normalize=False, weekday=0):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.weekday = weekday

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.weekday < 0 or self.weekday > 6:
            raise ValueError('Day must be 0<=day<=6, got {day}'
                             .format(day=self.weekday))

    def _get_offset_day(self, other):
        """
        Find the day in the same month as other that has the same
        weekday as self.weekday and is the last such day in the month.

        Parameters
        ----------
        other: datetime

        Returns
        -------
        day: int
        """
        dim = ccalendar.get_days_in_month(other.year, other.month)
        mend = datetime(other.year, other.month, dim)
        wday = mend.weekday()
        shift_days = (wday - self.weekday) % 7
        return dim - shift_days

    @property
    def rule_code(self):
        weekday = ccalendar.int_to_weekday.get(self.weekday, '')
        return '{prefix}-{weekday}'.format(prefix=self._prefix,
                                           weekday=weekday)

    @classmethod
    def _from_name(cls, suffix=None):
        if not suffix:
            raise ValueError("Prefix {prefix!r} requires a suffix."
                             .format(prefix=cls._prefix))
        # TODO: handle n here...
        weekday = ccalendar.weekday_to_int[suffix]
        return cls(weekday=weekday)

# ---------------------------------------------------------------------
# Quarter-Based Offset Classes


class QuarterOffset(DateOffset):
    """Quarter representation - doesn't call super"""
    _default_startingMonth = None
    _from_name_startingMonth = None
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize', 'startingMonth'])
    # TODO: Consider combining QuarterOffset and YearOffset __init__ at some
    #       point.  Also apply_index, onOffset, rule_code if
    #       startingMonth vs month attr names are resolved

    def __init__(self, n=1, normalize=False, startingMonth=None):
        self.n = self._validate_n(n)
        self.normalize = normalize
        if startingMonth is None:
            startingMonth = self._default_startingMonth
        self.startingMonth = startingMonth

    def isAnchored(self):
        return (self.n == 1 and self.startingMonth is not None)

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs['startingMonth'] = ccalendar.MONTH_TO_CAL_NUM[suffix]
        else:
            if cls._from_name_startingMonth is not None:
                kwargs['startingMonth'] = cls._from_name_startingMonth
        return cls(**kwargs)

    @property
    def rule_code(self):
        month = ccalendar.MONTH_ALIASES[self.startingMonth]
        return '{prefix}-{month}'.format(prefix=self._prefix, month=month)

    @apply_wraps
    def apply(self, other):
        # months_since: find the calendar quarter containing other.month,
        # e.g. if other.month == 8, the calendar quarter is [Jul, Aug, Sep].
        # Then find the month in that quarter containing an onOffset date for
        # self.  `months_since` is the number of months to shift other.month
        # to get to this on-offset month.
        months_since = other.month % 3 - self.startingMonth % 3
        qtrs = liboffsets.roll_qtrday(other, self.n, self.startingMonth,
                                      day_opt=self._day_opt, modby=3)
        months = qtrs * 3 - months_since
        return shift_month(other, months, self._day_opt)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        mod_month = (dt.month - self.startingMonth) % 3
        return mod_month == 0 and dt.day == self._get_offset_day(dt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = liboffsets.shift_quarters(dtindex.asi8, self.n,
                                            self.startingMonth, self._day_opt)
        return dtindex._shallow_copy(shifted)


class BQuarterEnd(QuarterOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...
    """
    _outputName = 'BusinessQuarterEnd'
    _default_startingMonth = 3
    _from_name_startingMonth = 12
    _prefix = 'BQ'
    _day_opt = 'business_end'


# TODO: This is basically the same as BQuarterEnd
class BQuarterBegin(QuarterOffset):
    _outputName = "BusinessQuarterBegin"
    # I suspect this is wrong for *all* of them.
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = 'BQS'
    _day_opt = 'business_start'


class QuarterEnd(QuarterOffset):
    """DateOffset increments between business Quarter dates
    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/31/2007, 6/30/2007, ...
    """
    _outputName = 'QuarterEnd'
    _default_startingMonth = 3
    _prefix = 'Q'
    _day_opt = 'end'


class QuarterBegin(QuarterOffset):
    _outputName = 'QuarterBegin'
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = 'QS'
    _day_opt = 'start'


# ---------------------------------------------------------------------
# Year-Based Offset Classes

class YearOffset(DateOffset):
    """DateOffset that just needs a month"""
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize', 'month'])

    def _get_offset_day(self, other):
        # override BaseOffset method to use self.month instead of other.month
        # TODO: there may be a more performant way to do this
        return liboffsets.get_day_of_month(other.replace(month=self.month),
                                           self._day_opt)

    @apply_wraps
    def apply(self, other):
        years = roll_yearday(other, self.n, self.month, self._day_opt)
        months = years * 12 + (self.month - other.month)
        return shift_month(other, months, self._day_opt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = liboffsets.shift_quarters(dtindex.asi8, self.n,
                                            self.month, self._day_opt,
                                            modby=12)
        return dtindex._shallow_copy(shifted)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return dt.month == self.month and dt.day == self._get_offset_day(dt)

    def __init__(self, n=1, normalize=False, month=None):
        self.n = self._validate_n(n)
        self.normalize = normalize

        month = month if month is not None else self._default_month
        self.month = month

        if self.month < 1 or self.month > 12:
            raise ValueError('Month must go from 1 to 12')

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs['month'] = ccalendar.MONTH_TO_CAL_NUM[suffix]
        return cls(**kwargs)

    @property
    def rule_code(self):
        month = ccalendar.MONTH_ALIASES[self.month]
        return '{prefix}-{month}'.format(prefix=self._prefix, month=month)


class BYearEnd(YearOffset):
    """DateOffset increments between business EOM dates"""
    _outputName = 'BusinessYearEnd'
    _default_month = 12
    _prefix = 'BA'
    _day_opt = 'business_end'


class BYearBegin(YearOffset):
    """DateOffset increments between business year begin dates"""
    _outputName = 'BusinessYearBegin'
    _default_month = 1
    _prefix = 'BAS'
    _day_opt = 'business_start'


class YearEnd(YearOffset):
    """DateOffset increments between calendar year ends"""
    _default_month = 12
    _prefix = 'A'
    _day_opt = 'end'


class YearBegin(YearOffset):
    """DateOffset increments between calendar year begin dates"""
    _default_month = 1
    _prefix = 'AS'
    _day_opt = 'start'


# ---------------------------------------------------------------------
# Special Offset Classes

class FY5253(DateOffset):
    """
    Describes 52-53 week fiscal year. This is also known as a 4-4-5 calendar.

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    http://en.wikipedia.org/wiki/4-4-5_calendar

    The year may either:
    - end on the last X day of the Y month.
    - end on the last X day closest to the last day of the Y month.

    X is a specific day of the week.
    Y is a certain month of the year

    Parameters
    ----------
    n : int
    weekday : {0, 1, ..., 6}
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    startingMonth : The month in which fiscal years end. {1, 2, ... 12}
    variation : str
        {"nearest", "last"} for "LastOfMonth" or "NearestEndMonth"
    """
    _prefix = 'RE'
    _adjust_dst = True
    _attributes = frozenset(['weekday', 'startingMonth', 'variation'])

    def __init__(self, n=1, normalize=False, weekday=0, startingMonth=1,
                 variation="nearest"):
        self.n = self._validate_n(n)
        self.normalize = normalize
        self.startingMonth = startingMonth
        self.weekday = weekday

        self.variation = variation

        if self.n == 0:
            raise ValueError('N cannot be 0')

        if self.variation not in ["nearest", "last"]:
            raise ValueError('{variation} is not a valid variation'
                             .format(variation=self.variation))

    def isAnchored(self):
        return (self.n == 1 and
                self.startingMonth is not None and
                self.weekday is not None)

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        dt = datetime(dt.year, dt.month, dt.day)
        year_end = self.get_year_end(dt)

        if self.variation == "nearest":
            # We have to check the year end of "this" cal year AND the previous
            return (year_end == dt or
                    self.get_year_end(shift_month(dt, -1, None)) == dt)
        else:
            return year_end == dt

    @apply_wraps
    def apply(self, other):
        norm = Timestamp(other).normalize()

        n = self.n
        prev_year = self.get_year_end(
            datetime(other.year - 1, self.startingMonth, 1))
        cur_year = self.get_year_end(
            datetime(other.year, self.startingMonth, 1))
        next_year = self.get_year_end(
            datetime(other.year + 1, self.startingMonth, 1))

        prev_year = tslib._localize_pydatetime(prev_year, other.tzinfo)
        cur_year = tslib._localize_pydatetime(cur_year, other.tzinfo)
        next_year = tslib._localize_pydatetime(next_year, other.tzinfo)

        # Note: next_year.year == other.year + 1, so we will always
        # have other < next_year
        if norm == prev_year:
            n -= 1
        elif norm == cur_year:
            pass
        elif n > 0:
            if norm < prev_year:
                n -= 2
            elif prev_year < norm < cur_year:
                n -= 1
            elif cur_year < norm < next_year:
                pass
        else:
            if cur_year < norm < next_year:
                n += 1
            elif prev_year < norm < cur_year:
                pass
            elif (norm.year == prev_year.year and norm < prev_year and
                  prev_year - norm <= timedelta(6)):
                # GH#14774, error when next_year.year == cur_year.year
                # e.g. prev_year == datetime(2004, 1, 3),
                # other == datetime(2004, 1, 1)
                n -= 1
            else:
                assert False

        shifted = datetime(other.year + n, self.startingMonth, 1)
        result = self.get_year_end(shifted)
        result = datetime(result.year, result.month, result.day,
                          other.hour, other.minute, other.second,
                          other.microsecond)
        return result

    def get_year_end(self, dt):
        assert dt.tzinfo is None

        dim = ccalendar.get_days_in_month(dt.year, self.startingMonth)
        target_date = datetime(dt.year, self.startingMonth, dim)
        wkday_diff = self.weekday - target_date.weekday()
        if wkday_diff == 0:
            # year_end is the same for "last" and "nearest" cases
            return target_date

        if self.variation == "last":
            days_forward = (wkday_diff % 7) - 7

            # days_forward is always negative, so we always end up
            # in the same year as dt
            return target_date + timedelta(days=days_forward)
        else:
            # variation == "nearest":
            days_forward = wkday_diff % 7
            if days_forward <= 3:
                # The upcoming self.weekday is closer than the previous one
                return target_date + timedelta(days_forward)
            else:
                # The previous self.weekday is closer than the upcoming one
                return target_date + timedelta(days_forward - 7)

    @property
    def rule_code(self):
        prefix = self._prefix
        suffix = self.get_rule_code_suffix()
        return "{prefix}-{suffix}".format(prefix=prefix, suffix=suffix)

    def _get_suffix_prefix(self):
        if self.variation == "nearest":
            return 'N'
        else:
            return 'L'

    def get_rule_code_suffix(self):
        prefix = self._get_suffix_prefix()
        month = ccalendar.MONTH_ALIASES[self.startingMonth]
        weekday = ccalendar.int_to_weekday[self.weekday]
        return '{prefix}-{month}-{weekday}'.format(prefix=prefix, month=month,
                                                   weekday=weekday)

    @classmethod
    def _parse_suffix(cls, varion_code, startingMonth_code, weekday_code):
        if varion_code == "N":
            variation = "nearest"
        elif varion_code == "L":
            variation = "last"
        else:
            raise ValueError("Unable to parse varion_code: "
                             "{code}".format(code=varion_code))

        startingMonth = ccalendar.MONTH_TO_CAL_NUM[startingMonth_code]
        weekday = ccalendar.weekday_to_int[weekday_code]

        return {"weekday": weekday,
                "startingMonth": startingMonth,
                "variation": variation}

    @classmethod
    def _from_name(cls, *args):
        return cls(**cls._parse_suffix(*args))


class FY5253Quarter(DateOffset):
    """
    DateOffset increments between business quarter dates
    for 52-53 week fiscal year (also known as a 4-4-5 calendar).

    It is used by companies that desire that their
    fiscal year always end on the same day of the week.

    It is a method of managing accounting periods.
    It is a common calendar structure for some industries,
    such as retail, manufacturing and parking industry.

    For more information see:
    http://en.wikipedia.org/wiki/4-4-5_calendar

    The year may either:
    - end on the last X day of the Y month.
    - end on the last X day closest to the last day of the Y month.

    X is a specific day of the week.
    Y is a certain month of the year

    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...

    Parameters
    ----------
    n : int
    weekday : {0, 1, ..., 6}
        0: Mondays
        1: Tuesdays
        2: Wednesdays
        3: Thursdays
        4: Fridays
        5: Saturdays
        6: Sundays
    startingMonth : The month in which fiscal years end. {1, 2, ... 12}
    qtr_with_extra_week : The quarter number that has the leap
        or 14 week when needed. {1,2,3,4}
    variation : str
        {"nearest", "last"} for "LastOfMonth" or "NearestEndMonth"
    """

    _prefix = 'REQ'
    _adjust_dst = True
    _attributes = frozenset(['weekday', 'startingMonth', 'qtr_with_extra_week',
                             'variation'])

    def __init__(self, n=1, normalize=False, weekday=0, startingMonth=1,
                 qtr_with_extra_week=1, variation="nearest"):
        self.n = self._validate_n(n)
        self.normalize = normalize

        self.weekday = weekday
        self.startingMonth = startingMonth
        self.qtr_with_extra_week = qtr_with_extra_week
        self.variation = variation

        if self.n == 0:
            raise ValueError('N cannot be 0')

    @cache_readonly
    def _offset(self):
        return FY5253(startingMonth=self.startingMonth,
                      weekday=self.weekday,
                      variation=self.variation)

    def isAnchored(self):
        return self.n == 1 and self._offset.isAnchored()

    def _rollback_to_year(self, other):
        """roll `other` back to the most recent date that was on a fiscal year
        end.  Return the date of that year-end, the number of full quarters
        elapsed between that year-end and other, and the remaining Timedelta
        since the most recent quarter-end.

        Parameters
        ----------
        other : datetime or Timestamp

        Returns
        -------
        tuple of
        prev_year_end : Timestamp giving most recent fiscal year end
        num_qtrs : int
        tdelta : Timedelta
        """
        num_qtrs = 0

        norm = Timestamp(other).tz_localize(None)
        start = self._offset.rollback(norm)
        # Note: start <= norm and self._offset.onOffset(start)

        if start < norm:
            # roll adjustment
            qtr_lens = self.get_weeks(norm)

            # check thet qtr_lens is consistent with self._offset addition
            end = shift_day(start, days=7 * sum(qtr_lens))
            assert self._offset.onOffset(end), (start, end, qtr_lens)

            tdelta = norm - start
            for qlen in qtr_lens:
                if qlen * 7 <= tdelta.days:
                    num_qtrs += 1
                    tdelta -= Timedelta(days=qlen * 7)
                else:
                    break
        else:
            tdelta = Timedelta(0)

        # Note: we always have tdelta.value >= 0
        return start, num_qtrs, tdelta

    @apply_wraps
    def apply(self, other):
        # Note: self.n == 0 is not allowed.
        n = self.n

        prev_year_end, num_qtrs, tdelta = self._rollback_to_year(other)
        res = prev_year_end
        n += num_qtrs
        if self.n <= 0 and tdelta.value > 0:
            n += 1

        # Possible speedup by handling years first.
        years = n // 4
        if years:
            res += self._offset * years
            n -= years * 4

        # Add an extra day to make *sure* we are getting the quarter lengths
        # for the upcoming year, not the previous year
        qtr_lens = self.get_weeks(res + Timedelta(days=1))

        # Note: we always have 0 <= n < 4
        weeks = sum(qtr_lens[:n])
        if weeks:
            res = shift_day(res, days=weeks * 7)

        return res

    def get_weeks(self, dt):
        ret = [13] * 4

        year_has_extra_week = self.year_has_extra_week(dt)

        if year_has_extra_week:
            ret[self.qtr_with_extra_week - 1] = 14

        return ret

    def year_has_extra_week(self, dt):
        # Avoid round-down errors --> normalize to get
        # e.g. '370D' instead of '360D23H'
        norm = Timestamp(dt).normalize().tz_localize(None)

        next_year_end = self._offset.rollforward(norm)
        prev_year_end = norm - self._offset
        weeks_in_year = (next_year_end - prev_year_end).days / 7
        assert weeks_in_year in [52, 53], weeks_in_year
        return weeks_in_year == 53

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        if self._offset.onOffset(dt):
            return True

        next_year_end = dt - self._offset

        qtr_lens = self.get_weeks(dt)

        current = next_year_end
        for qtr_len in qtr_lens:
            current = shift_day(current, days=qtr_len * 7)
            if dt == current:
                return True
        return False

    @property
    def rule_code(self):
        suffix = self._offset.get_rule_code_suffix()
        qtr = self.qtr_with_extra_week
        return "{prefix}-{suffix}-{qtr}".format(prefix=self._prefix,
                                                suffix=suffix, qtr=qtr)

    @classmethod
    def _from_name(cls, *args):
        return cls(**dict(FY5253._parse_suffix(*args[:-1]),
                          qtr_with_extra_week=int(args[-1])))


class Easter(DateOffset):
    """
    DateOffset for the Easter holiday using
    logic defined in dateutil.  Right now uses
    the revised method which is valid in years
    1583-4099.
    """
    _adjust_dst = True
    _attributes = frozenset(['n', 'normalize'])

    def __init__(self, n=1, normalize=False):
        self.n = self._validate_n(n)
        self.normalize = normalize

    @apply_wraps
    def apply(self, other):
        current_easter = easter(other.year)
        current_easter = datetime(current_easter.year,
                                  current_easter.month, current_easter.day)
        current_easter = tslib._localize_pydatetime(current_easter,
                                                    other.tzinfo)

        n = self.n
        if n >= 0 and other < current_easter:
            n -= 1
        elif n < 0 and other > current_easter:
            n += 1
        # TODO: Why does this handle the 0 case the opposite of others?

        # NOTE: easter returns a datetime.date so we have to convert to type of
        # other
        new = easter(other.year + n)
        new = datetime(new.year, new.month, new.day, other.hour,
                       other.minute, other.second, other.microsecond)
        return new

    def onOffset(self, dt):
        if self.normalize and not _is_normalized(dt):
            return False
        return date(dt.year, dt.month, dt.day) == easter(dt.year)

# ---------------------------------------------------------------------
# Ticks


def _tick_comp(op):
    def f(self, other):
        return op(self.delta, other.delta)

    return f


class Tick(SingleConstructorOffset):
    _inc = Timedelta(microseconds=1000)
    _prefix = 'undefined'
    _attributes = frozenset(['n', 'normalize'])

    def __init__(self, n=1, normalize=False):
        # TODO: do Tick classes with normalize=True make sense?
        self.n = self._validate_n(n)
        self.normalize = normalize

    __gt__ = _tick_comp(operator.gt)
    __ge__ = _tick_comp(operator.ge)
    __lt__ = _tick_comp(operator.lt)
    __le__ = _tick_comp(operator.le)
    __eq__ = _tick_comp(operator.eq)
    __ne__ = _tick_comp(operator.ne)

    def __add__(self, other):
        if isinstance(other, Tick):
            if type(self) == type(other):
                return type(self)(self.n + other.n)
            else:
                return _delta_to_tick(self.delta + other.delta)
        elif isinstance(other, ABCPeriod):
            return other + self
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented
        except OverflowError:
            raise OverflowError("the add operation between {self} and {other} "
                                "will overflow".format(self=self, other=other))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta == other.delta
        else:
            # TODO: Are there cases where this should raise TypeError?
            return False

    # This is identical to DateOffset.__hash__, but has to be redefined here
    # for Python 3, because we've redefined __eq__.
    def __hash__(self):
        return hash(self._params())

    def __ne__(self, other):
        if isinstance(other, compat.string_types):
            from pandas.tseries.frequencies import to_offset

            other = to_offset(other)

        if isinstance(other, Tick):
            return self.delta != other.delta
        else:
            # TODO: Are there cases where this should raise TypeError?
            return True

    @property
    def delta(self):
        return self.n * self._inc

    @property
    def nanos(self):
        return delta_to_nanoseconds(self.delta)

    # TODO: Should Tick have its own apply_index?
    def apply(self, other):
        # Timestamp can handle tz and nano sec, thus no need to use apply_wraps
        if isinstance(other, Timestamp):

            # GH 15126
            # in order to avoid a recursive
            # call of __add__ and __radd__ if there is
            # an exception, when we call using the + operator,
            # we directly call the known method
            result = other.__add__(self)
            if result == NotImplemented:
                raise OverflowError
            return result
        elif isinstance(other, (datetime, np.datetime64, date)):
            return as_timestamp(other) + self

        if isinstance(other, timedelta):
            return other + self.delta
        elif isinstance(other, type(self)):
            return type(self)(self.n + other.n)

        raise ApplyTypeError('Unhandled type: {type_str}'
                             .format(type_str=type(other).__name__))

    def isAnchored(self):
        return False


def _delta_to_tick(delta):
    if delta.microseconds == 0:
        if delta.seconds == 0:
            return Day(delta.days)
        else:
            seconds = delta.days * 86400 + delta.seconds
            if seconds % 3600 == 0:
                return Hour(seconds / 3600)
            elif seconds % 60 == 0:
                return Minute(seconds / 60)
            else:
                return Second(seconds)
    else:
        nanos = delta_to_nanoseconds(delta)
        if nanos % 1000000 == 0:
            return Milli(nanos // 1000000)
        elif nanos % 1000 == 0:
            return Micro(nanos // 1000)
        else:  # pragma: no cover
            return Nano(nanos)


class Day(Tick):
    _inc = Timedelta(days=1)
    _prefix = 'D'


class Hour(Tick):
    _inc = Timedelta(hours=1)
    _prefix = 'H'


class Minute(Tick):
    _inc = Timedelta(minutes=1)
    _prefix = 'T'


class Second(Tick):
    _inc = Timedelta(seconds=1)
    _prefix = 'S'


class Milli(Tick):
    _inc = Timedelta(milliseconds=1)
    _prefix = 'L'


class Micro(Tick):
    _inc = Timedelta(microseconds=1)
    _prefix = 'U'


class Nano(Tick):
    _inc = Timedelta(nanoseconds=1)
    _prefix = 'N'


BDay = BusinessDay
BMonthEnd = BusinessMonthEnd
BMonthBegin = BusinessMonthBegin
CBMonthEnd = CustomBusinessMonthEnd
CBMonthBegin = CustomBusinessMonthBegin
CDay = CustomBusinessDay

# ---------------------------------------------------------------------


def generate_range(start=None, end=None, periods=None,
                   offset=BDay(), time_rule=None):
    """
    Generates a sequence of dates corresponding to the specified time
    offset. Similar to dateutil.rrule except uses pandas DateOffset
    objects to represent time increments

    Parameters
    ----------
    start : datetime (default None)
    end : datetime (default None)
    periods : int, optional
    time_rule : (legacy) name of DateOffset object to be used, optional
        Corresponds with names expected by tseries.frequencies.get_offset

    Notes
    -----
    * This method is faster for generating weekdays than dateutil.rrule
    * At least two of (start, end, periods) must be specified.
    * If both start and end are specified, the returned dates will
    satisfy start <= date <= end.
    * If both time_rule and offset are specified, time_rule supersedes offset.

    Returns
    -------
    dates : generator object

    """
    if time_rule is not None:
        from pandas.tseries.frequencies import get_offset

        offset = get_offset(time_rule)

    start = to_datetime(start)
    end = to_datetime(end)

    if start and not offset.onOffset(start):
        start = offset.rollforward(start)

    elif end and not offset.onOffset(end):
        end = offset.rollback(end)

    if periods is None and end < start:
        end = None
        periods = 0

    if end is None:
        end = start + (periods - 1) * offset

    if start is None:
        start = end - (periods - 1) * offset

    cur = start
    if offset.n >= 0:
        while cur <= end:
            yield cur

            # faster than cur + offset
            next_date = offset.apply(cur)
            if next_date <= cur:
                raise ValueError('Offset {offset} did not increment date'
                                 .format(offset=offset))
            cur = next_date
    else:
        while cur >= end:
            yield cur

            # faster than cur + offset
            next_date = offset.apply(cur)
            if next_date >= cur:
                raise ValueError('Offset {offset} did not decrement date'
                                 .format(offset=offset))
            cur = next_date


prefix_mapping = dict((offset._prefix, offset) for offset in [
    YearBegin,                 # 'AS'
    YearEnd,                   # 'A'
    BYearBegin,                # 'BAS'
    BYearEnd,                  # 'BA'
    BusinessDay,               # 'B'
    BusinessMonthBegin,        # 'BMS'
    BusinessMonthEnd,          # 'BM'
    BQuarterEnd,               # 'BQ'
    BQuarterBegin,             # 'BQS'
    BusinessHour,              # 'BH'
    CustomBusinessDay,         # 'C'
    CustomBusinessMonthEnd,    # 'CBM'
    CustomBusinessMonthBegin,  # 'CBMS'
    CustomBusinessHour,        # 'CBH'
    MonthEnd,                  # 'M'
    MonthBegin,                # 'MS'
    Nano,                      # 'N'
    SemiMonthEnd,              # 'SM'
    SemiMonthBegin,            # 'SMS'
    Week,                      # 'W'
    Second,                    # 'S'
    Minute,                    # 'T'
    Micro,                     # 'U'
    QuarterEnd,                # 'Q'
    QuarterBegin,              # 'QS'
    Milli,                     # 'L'
    Hour,                      # 'H'
    Day,                       # 'D'
    WeekOfMonth,               # 'WOM'
    FY5253,
    FY5253Quarter,
])
