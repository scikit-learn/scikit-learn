"""
Matplotlib provides sophisticated date plotting capabilities, standing on the
shoulders of python :mod:`datetime`, the add-on modules :mod:`pytz` and
:mod:`dateutil`.


.. _date-format:

Matplotlib date format
----------------------
Matplotlib represents dates using floating point numbers specifying the number
of days since 0001-01-01 UTC, plus 1.  For example, 0001-01-01, 06:00 is 1.25,
not 0.25. Values < 1, i.e. dates before 0001-01-01 UTC are not supported.

There are a number of helper functions to convert between :mod:`datetime`
objects and Matplotlib dates:

.. currentmodule:: matplotlib.dates

.. autosummary::
   :nosignatures:

   date2num
   num2date
   num2timedelta
   epoch2num
   num2epoch
   mx2num
   drange

.. note::

   Like Python's datetime, mpl uses the Gregorian calendar for all
   conversions between dates and floating point numbers. This practice
   is not universal, and calendar differences can cause confusing
   differences between what Python and mpl give as the number of days
   since 0001-01-01 and what other software and databases yield.  For
   example, the US Naval Observatory uses a calendar that switches
   from Julian to Gregorian in October, 1582.  Hence, using their
   calculator, the number of days between 0001-01-01 and 2006-04-01 is
   732403, whereas using the Gregorian calendar via the datetime
   module we find::

     In [1]: date(2006, 4, 1).toordinal() - date(1, 1, 1).toordinal()
     Out[1]: 732401

All the Matplotlib date converters, tickers and formatters are timezone aware.
If no explicit timezone is provided, the rcParam ``timezone`` is assumend.  If
you want to use a custom time zone, pass a :class:`pytz.timezone` instance
with the tz keyword argument to :func:`num2date`, :func:`.plot_date`, and any
custom date tickers or locators you create.
See `pytz <http://pythonhosted.org/pytz/>`_ for information on :mod:`pytz` and
timezone handling.

A wide range of specific and general purpose date tick locators and
formatters are provided in this module.  See
:mod:`matplotlib.ticker` for general information on tick locators
and formatters.  These are described below.


The `dateutil module <https://dateutil.readthedocs.io/en/stable/>`_ provides
additional code to handle date ticking, making it easy to place ticks
on any kinds of dates.  See examples below.

Date tickers
------------

Most of the date tickers can locate single or multiple values.  For
example::

    # import constants for the days of the week
    from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU

    # tick on mondays every week
    loc = WeekdayLocator(byweekday=MO, tz=tz)

    # tick on mondays and saturdays
    loc = WeekdayLocator(byweekday=(MO, SA))

In addition, most of the constructors take an interval argument::

    # tick on mondays every second week
    loc = WeekdayLocator(byweekday=MO, interval=2)

The rrule locator allows completely general date ticking::

    # tick every 5th easter
    rule = rrulewrapper(YEARLY, byeaster=1, interval=5)
    loc = RRuleLocator(rule)

Here are all the date tickers:

    * :class:`MicrosecondLocator`: locate microseconds

    * :class:`SecondLocator`: locate seconds

    * :class:`MinuteLocator`: locate minutes

    * :class:`HourLocator`: locate hours

    * :class:`DayLocator`: locate specified days of the month

    * :class:`WeekdayLocator`: Locate days of the week, e.g., MO, TU

    * :class:`MonthLocator`: locate months, e.g., 7 for july

    * :class:`YearLocator`: locate years that are multiples of base

    * :class:`RRuleLocator`: locate using a
      :class:`matplotlib.dates.rrulewrapper`.  The
      :class:`rrulewrapper` is a simple wrapper around a
      :class:`dateutil.rrule` (`dateutil
      <https://dateutil.readthedocs.io/en/stable/>`_) which allow almost
      arbitrary date tick specifications.  See `rrule example
      <../gallery/ticks_and_spines/date_demo_rrule.html>`_.

    * :class:`AutoDateLocator`: On autoscale, this class picks the best
      :class:`DateLocator` (e.g., :class:`RRuleLocator`)
      to set the view limits and the tick
      locations.  If called with ``interval_multiples=True`` it will
      make ticks line up with sensible multiples of the tick intervals.  E.g.
      if the interval is 4 hours, it will pick hours 0, 4, 8, etc as ticks.
      This behaviour is not guaranteed by default.

Date formatters
---------------

Here all all the date formatters:

    * :class:`AutoDateFormatter`: attempts to figure out the best format
      to use.  This is most useful when used with the :class:`AutoDateLocator`.

    * :class:`DateFormatter`: use :func:`strftime` format strings

    * :class:`IndexDateFormatter`: date plots with implicit *x*
      indexing.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip
import re
import time
import math
import datetime
import functools

import warnings
import logging

from dateutil.rrule import (rrule, MO, TU, WE, TH, FR, SA, SU, YEARLY,
                            MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY,
                            SECONDLY)
from dateutil.relativedelta import relativedelta
import dateutil.parser
import logging
import numpy as np


import matplotlib
from matplotlib import rcParams
import matplotlib.units as units
import matplotlib.cbook as cbook
import matplotlib.ticker as ticker

_log = logging.getLogger(__name__)

__all__ = ('date2num', 'num2date', 'num2timedelta', 'drange', 'epoch2num',
           'num2epoch', 'mx2num', 'DateFormatter',
           'IndexDateFormatter', 'AutoDateFormatter', 'DateLocator',
           'RRuleLocator', 'AutoDateLocator', 'YearLocator',
           'MonthLocator', 'WeekdayLocator',
           'DayLocator', 'HourLocator', 'MinuteLocator',
           'SecondLocator', 'MicrosecondLocator',
           'rrule', 'MO', 'TU', 'WE', 'TH', 'FR', 'SA', 'SU',
           'YEARLY', 'MONTHLY', 'WEEKLY', 'DAILY',
           'HOURLY', 'MINUTELY', 'SECONDLY', 'MICROSECONDLY', 'relativedelta',
           'seconds', 'minutes', 'hours', 'weeks')


_log = logging.getLogger(__name__)


# Make a simple UTC instance so we don't always have to import
# pytz.  From the python datetime library docs:

class _UTC(datetime.tzinfo):
    """UTC"""

    def utcoffset(self, dt):
        return datetime.timedelta(0)

    def tzname(self, dt):
        return str("UTC")

    def dst(self, dt):
        return datetime.timedelta(0)


UTC = _UTC()


def _get_rc_timezone():
    """
    Retrieve the preferred timeszone from the rcParams dictionary.
    """
    s = matplotlib.rcParams['timezone']
    if s == 'UTC':
        return UTC
    import pytz
    return pytz.timezone(s)


"""
Time-related constants.
"""
EPOCH_OFFSET = float(datetime.datetime(1970, 1, 1).toordinal())
JULIAN_OFFSET = 1721424.5                         # Julian date at 0001-01-01
MICROSECONDLY = SECONDLY + 1
HOURS_PER_DAY = 24.
MIN_PER_HOUR = 60.
SEC_PER_MIN = 60.
MONTHS_PER_YEAR = 12.

DAYS_PER_WEEK = 7.
DAYS_PER_MONTH = 30.
DAYS_PER_YEAR = 365.0

MINUTES_PER_DAY = MIN_PER_HOUR * HOURS_PER_DAY

SEC_PER_HOUR = SEC_PER_MIN * MIN_PER_HOUR
SEC_PER_DAY = SEC_PER_HOUR * HOURS_PER_DAY
SEC_PER_WEEK = SEC_PER_DAY * DAYS_PER_WEEK

MUSECONDS_PER_DAY = 1e6 * SEC_PER_DAY

MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY = (
    MO, TU, WE, TH, FR, SA, SU)
WEEKDAYS = (MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY)


def _to_ordinalf(dt):
    """
    Convert :mod:`datetime` or :mod:`date` to the Gregorian date as UTC float
    days, preserving hours, minutes, seconds and microseconds.  Return value
    is a :func:`float`.
    """
    # Convert to UTC
    tzi = getattr(dt, 'tzinfo', None)
    if tzi is not None:
        dt = dt.astimezone(UTC)
        tzi = UTC

    base = float(dt.toordinal())

    # If it's sufficiently datetime-like, it will have a `date()` method
    cdate = getattr(dt, 'date', lambda: None)()
    if cdate is not None:
        # Get a datetime object at midnight UTC
        midnight_time = datetime.time(0, tzinfo=tzi)

        rdt = datetime.datetime.combine(cdate, midnight_time)

        # Append the seconds as a fraction of a day
        base += (dt - rdt).total_seconds() / SEC_PER_DAY

    return base


# a version of _to_ordinalf that can operate on numpy arrays
_to_ordinalf_np_vectorized = np.vectorize(_to_ordinalf)


def _dt64_to_ordinalf(d):
    """
    Convert `numpy.datetime64` or an ndarray of those types to Gregorian
    date as UTC float.  Roundoff is via float64 precision.  Practically:
    microseconds for dates between 290301 BC, 294241 AD, milliseconds for
    larger dates (see `numpy.datetime64`).  Nanoseconds aren't possible
    because we do times compared to ``0001-01-01T00:00:00`` (plus one day).
    """

    # the "extra" ensures that we at least allow the dynamic range out to
    # seconds.  That should get out to +/-2e11 years.
    extra = d - d.astype('datetime64[s]')
    extra = extra.astype('timedelta64[ns]')
    t0 = np.datetime64('0001-01-01T00:00:00').astype('datetime64[s]')
    dt = (d.astype('datetime64[s]') - t0).astype(np.float64)
    dt += extra.astype(np.float64) / 1.0e9
    dt = dt / SEC_PER_DAY + 1.0

    NaT_int = np.datetime64('NaT').astype(np.int64)
    d_int = d.astype(np.int64)
    try:
        dt[d_int == NaT_int] = np.nan
    except TypeError:
        if d_int == NaT_int:
            dt = np.nan
    return dt


def _from_ordinalf(x, tz=None):
    """
    Convert Gregorian float of the date, preserving hours, minutes,
    seconds and microseconds.  Return value is a `.datetime`.

    The input date *x* is a float in ordinal days at UTC, and the output will
    be the specified `.datetime` object corresponding to that time in
    timezone *tz*, or if *tz* is ``None``, in the timezone specified in
    :rc:`timezone`.
    """
    if tz is None:
        tz = _get_rc_timezone()

    ix, remainder = divmod(x, 1)
    ix = int(ix)
    if ix < 1:
        raise ValueError('Cannot convert {} to a date.  This often happens if '
                         'non-datetime values are passed to an axis that '
                         'expects datetime objects.'.format(ix))
    dt = datetime.datetime.fromordinal(ix).replace(tzinfo=UTC)

    # Since the input date `x` float is unable to preserve microsecond
    # precision of time representation in non-antique years, the
    # resulting datetime is rounded to the nearest multiple of
    # `musec_prec`. A value of 20 is appropriate for current dates.
    musec_prec = 20
    remainder_musec = int(round(remainder * MUSECONDS_PER_DAY / musec_prec)
                          * musec_prec)

    # For people trying to plot with full microsecond precision, enable
    # an early-year workaround
    if x < 30 * 365:
        remainder_musec = int(round(remainder * MUSECONDS_PER_DAY))

    # add hours, minutes, seconds, microseconds
    dt += datetime.timedelta(microseconds=remainder_musec)

    return dt.astimezone(tz)


# a version of _from_ordinalf that can operate on numpy arrays
_from_ordinalf_np_vectorized = np.vectorize(_from_ordinalf)


class strpdate2num(object):
    """
    Use this class to parse date strings to matplotlib datenums when
    you know the date format string of the date you are parsing.
    """
    def __init__(self, fmt):
        """ fmt: any valid strptime format is supported """
        self.fmt = fmt

    def __call__(self, s):
        """s : string to be converted
           return value: a date2num float
        """
        return date2num(datetime.datetime(*time.strptime(s, self.fmt)[:6]))


class bytespdate2num(strpdate2num):
    """
    Use this class to parse date strings to matplotlib datenums when
    you know the date format string of the date you are parsing.  See
    :file:`examples/misc/load_converter.py`.
    """
    def __init__(self, fmt, encoding='utf-8'):
        """
        Args:
            fmt: any valid strptime format is supported
            encoding: encoding to use on byte input (default: 'utf-8')
        """
        super(bytespdate2num, self).__init__(fmt)
        self.encoding = encoding

    def __call__(self, b):
        """
        Args:
            b: byte input to be converted
        Returns:
            A date2num float
        """
        s = b.decode(self.encoding)
        return super(bytespdate2num, self).__call__(s)


# a version of dateutil.parser.parse that can operate on nump0y arrays
_dateutil_parser_parse_np_vectorized = np.vectorize(dateutil.parser.parse)


def datestr2num(d, default=None):
    """
    Convert a date string to a datenum using
    :func:`dateutil.parser.parse`.

    Parameters
    ----------
    d : string or sequence of strings
        The dates to convert.

    default : datetime instance, optional
        The default date to use when fields are missing in *d*.
    """
    if isinstance(d, six.string_types):
        dt = dateutil.parser.parse(d, default=default)
        return date2num(dt)
    else:
        if default is not None:
            d = [dateutil.parser.parse(s, default=default) for s in d]
        d = np.asarray(d)
        if not d.size:
            return d
        return date2num(_dateutil_parser_parse_np_vectorized(d))


def date2num(d):
    """
    Convert datetime objects to Matplotlib dates.

    Parameters
    ----------
    d : `datetime.datetime` or `numpy.datetime64` or sequences of these

    Returns
    -------
    float or sequence of floats
        Number of days (fraction part represents hours, minutes, seconds, ms)
        since 0001-01-01 00:00:00 UTC, plus one.

    Notes
    -----
    The addition of one here is a historical artifact. Also, note that the
    Gregorian calendar is assumed; this is not universal practice.
    For details see the module docstring.
    """

    if hasattr(d, "values"):
        # this unpacks pandas series or dataframes...
        d = d.values

    if ((isinstance(d, np.ndarray) and np.issubdtype(d.dtype, np.datetime64))
            or isinstance(d, np.datetime64)):
        return _dt64_to_ordinalf(d)
    if not cbook.iterable(d):
        return _to_ordinalf(d)
    else:
        d = np.asarray(d)
        if not d.size:
            return d
        return _to_ordinalf_np_vectorized(d)


def julian2num(j):
    """
    Convert a Julian date (or sequence) to a Matplotlib date (or sequence).

    Parameters
    ----------
    j : float or sequence of floats
        Julian date(s)

    Returns
    -------
    float or sequence of floats
        Matplotlib date(s)
    """
    if cbook.iterable(j):
        j = np.asarray(j)
    return j - JULIAN_OFFSET


def num2julian(n):
    """
    Convert a Matplotlib date (or sequence) to a Julian date (or sequence).

    Parameters
    ----------
    n : float or sequence of floats
        Matplotlib date(s)

    Returns
    -------
    float or sequence of floats
        Julian date(s)
    """
    if cbook.iterable(n):
        n = np.asarray(n)
    return n + JULIAN_OFFSET


def num2date(x, tz=None):
    """
    Convert Matplotlib dates to `~datetime.datetime` objects.

    Parameters
    ----------
    x : float or sequence of floats
        Number of days (fraction part represents hours, minutes, seconds)
        since 0001-01-01 00:00:00 UTC, plus one.
    tz : string, optional
        Timezone of *x* (defaults to rcparams ``timezone``).

    Returns
    -------
    `~datetime.datetime` or sequence of `~datetime.datetime`
        Dates are returned in timezone *tz*.

        If *x* is a sequence, a sequence of :class:`datetime` objects will
        be returned.

    Notes
    -----
    The addition of one here is a historical artifact. Also, note that the
    Gregorian calendar is assumed; this is not universal practice.
    For details, see the module docstring.
    """
    if tz is None:
        tz = _get_rc_timezone()
    if not cbook.iterable(x):
        return _from_ordinalf(x, tz)
    else:
        x = np.asarray(x)
        if not x.size:
            return x
        return _from_ordinalf_np_vectorized(x, tz).tolist()


def _ordinalf_to_timedelta(x):
    return datetime.timedelta(days=x)


_ordinalf_to_timedelta_np_vectorized = np.vectorize(_ordinalf_to_timedelta)


def num2timedelta(x):
    """
    Convert number of days to a `~datetime.timedelta` object.

    If *x* is a sequence, a sequence of `~datetime.timedelta` objects will
    be returned.

    Parameters
    ----------
    x : float, sequence of floats
        Number of days. The fraction part represents hours, minutes, seconds.

    Returns
    -------
    `datetime.timedelta` or list[`datetime.timedelta`]

    """
    if not cbook.iterable(x):
        return _ordinalf_to_timedelta(x)
    else:
        x = np.asarray(x)
        if not x.size:
            return x
        return _ordinalf_to_timedelta_np_vectorized(x).tolist()


def drange(dstart, dend, delta):
    """
    Return a sequence of equally spaced Matplotlib dates.

    The dates start at *dstart* and reach up to, but not including *dend*.
    They are spaced by *delta*.

    Parameters
    ----------
    dstart, dend : `~datetime.datetime`
        The date limits.
    delta : `datetime.timedelta`
        Spacing of the dates.

    Returns
    -------
    drange : `numpy.array`
        A list floats representing Matplotlib dates.

    """
    f1 = date2num(dstart)
    f2 = date2num(dend)
    step = delta.total_seconds() / SEC_PER_DAY

    # calculate the difference between dend and dstart in times of delta
    num = int(np.ceil((f2 - f1) / step))

    # calculate end of the interval which will be generated
    dinterval_end = dstart + num * delta

    # ensure, that an half open interval will be generated [dstart, dend)
    if dinterval_end >= dend:
        # if the endpoint is greated than dend, just subtract one delta
        dinterval_end -= delta
        num -= 1

    f2 = date2num(dinterval_end)  # new float-endpoint
    return np.linspace(f1, f2, num + 1)

### date tickers and formatters ###


class DateFormatter(ticker.Formatter):
    """
    Tick location is seconds since the epoch.  Use a :func:`strftime`
    format string.

    Python only supports :mod:`datetime` :func:`strftime` formatting
    for years greater than 1900.  Thanks to Andrew Dalke, Dalke
    Scientific Software who contributed the :func:`strftime` code
    below to include dates earlier than this year.
    """

    illegal_s = re.compile(r"((^|[^%])(%%)*%s)")

    def __init__(self, fmt, tz=None):
        """
        *fmt* is a :func:`strftime` format string; *tz* is the
         :class:`tzinfo` instance.
        """
        if tz is None:
            tz = _get_rc_timezone()
        self.fmt = fmt
        self.tz = tz

    def __call__(self, x, pos=0):
        if x == 0:
            raise ValueError('DateFormatter found a value of x=0, which is '
                             'an illegal date.  This usually occurs because '
                             'you have not informed the axis that it is '
                             'plotting dates, e.g., with ax.xaxis_date()')
        dt = num2date(x, self.tz)
        return self.strftime(dt, self.fmt)

    def set_tzinfo(self, tz):
        self.tz = tz

    def _replace_common_substr(self, s1, s2, sub1, sub2, replacement):
        """Helper function for replacing substrings sub1 and sub2
        located at the same indexes in strings s1 and s2 respectively,
        with the string replacement.  It is expected that sub1 and sub2
        have the same length.  Returns the pair s1, s2 after the
        substitutions.
        """
        # Find common indexes of substrings sub1 in s1 and sub2 in s2
        # and make substitutions inplace. Because this is inplace,
        # it is okay if len(replacement) != len(sub1), len(sub2).
        i = 0
        while True:
            j = s1.find(sub1, i)
            if j == -1:
                break

            i = j + 1
            if s2[j:j + len(sub2)] != sub2:
                continue

            s1 = s1[:j] + replacement + s1[j + len(sub1):]
            s2 = s2[:j] + replacement + s2[j + len(sub2):]

        return s1, s2

    def strftime_pre_1900(self, dt, fmt=None):
        """Call time.strftime for years before 1900 by rolling
        forward a multiple of 28 years.

        *fmt* is a :func:`strftime` format string.

        Dalke: I hope I did this math right.  Every 28 years the
        calendar repeats, except through century leap years excepting
        the 400 year leap years.  But only if you're using the Gregorian
        calendar.
        """
        if fmt is None:
            fmt = self.fmt

        # Since python's time module's strftime implementation does not
        # support %f microsecond (but the datetime module does), use a
        # regular expression substitution to replace instances of %f.
        # Note that this can be useful since python's floating-point
        # precision representation for datetime causes precision to be
        # more accurate closer to year 0 (around the year 2000, precision
        # can be at 10s of microseconds).
        fmt = re.sub(r'((^|[^%])(%%)*)%f',
                     r'\g<1>{0:06d}'.format(dt.microsecond), fmt)

        year = dt.year
        # For every non-leap year century, advance by
        # 6 years to get into the 28-year repeat cycle
        delta = 2000 - year
        off = 6 * (delta // 100 + delta // 400)
        year = year + off

        # Move to between the years 1973 and 2000
        year1 = year + ((2000 - year) // 28) * 28
        year2 = year1 + 28
        timetuple = dt.timetuple()
        # Generate timestamp string for year and year+28
        s1 = time.strftime(fmt, (year1,) + timetuple[1:])
        s2 = time.strftime(fmt, (year2,) + timetuple[1:])

        # Replace instances of respective years (both 2-digit and 4-digit)
        # that are located at the same indexes of s1, s2 with dt's year.
        # Note that C++'s strftime implementation does not use padded
        # zeros or padded whitespace for %y or %Y for years before 100, but
        # uses padded zeros for %x. (For example, try the runnable examples
        # with .tm_year in the interval [-1900, -1800] on
        # http://en.cppreference.com/w/c/chrono/strftime.) For ease of
        # implementation, we always use padded zeros for %y, %Y, and %x.
        s1, s2 = self._replace_common_substr(s1, s2,
                                             "{0:04d}".format(year1),
                                             "{0:04d}".format(year2),
                                             "{0:04d}".format(dt.year))
        s1, s2 = self._replace_common_substr(s1, s2,
                                             "{0:02d}".format(year1 % 100),
                                             "{0:02d}".format(year2 % 100),
                                             "{0:02d}".format(dt.year % 100))
        return cbook.unicode_safe(s1)

    def strftime(self, dt, fmt=None):
        """
        Refer to documentation for :meth:`datetime.datetime.strftime`

        *fmt* is a :meth:`datetime.datetime.strftime` format string.

        Warning: For years before 1900, depending upon the current
        locale it is possible that the year displayed with %x might
        be incorrect. For years before 100, %y and %Y will yield
        zero-padded strings.
        """
        if fmt is None:
            fmt = self.fmt
        fmt = self.illegal_s.sub(r"\1", fmt)
        fmt = fmt.replace("%s", "s")
        if dt.year >= 1900:
            # Note: in python 3.3 this is okay for years >= 1000,
            # refer to http://bugs.python.org/issue1777412
            return cbook.unicode_safe(dt.strftime(fmt))

        return self.strftime_pre_1900(dt, fmt)


class IndexDateFormatter(ticker.Formatter):
    """
    Use with :class:`~matplotlib.ticker.IndexLocator` to cycle format
    strings by index.
    """
    def __init__(self, t, fmt, tz=None):
        """
        *t* is a sequence of dates (floating point days).  *fmt* is a
        :func:`strftime` format string.
        """
        if tz is None:
            tz = _get_rc_timezone()
        self.t = t
        self.fmt = fmt
        self.tz = tz

    def __call__(self, x, pos=0):
        'Return the label for time *x* at position *pos*'
        ind = int(np.round(x))
        if ind >= len(self.t) or ind <= 0:
            return ''

        dt = num2date(self.t[ind], self.tz)

        return cbook.unicode_safe(dt.strftime(self.fmt))


class AutoDateFormatter(ticker.Formatter):
    """
    This class attempts to figure out the best format to use.  This is
    most useful when used with the :class:`AutoDateLocator`.


    The AutoDateFormatter has a scale dictionary that maps the scale
    of the tick (the distance in days between one major tick) and a
    format string.  The default looks like this::

        self.scaled = {
            DAYS_PER_YEAR: rcParams['date.autoformat.year'],
            DAYS_PER_MONTH: rcParams['date.autoformat.month'],
            1.0: rcParams['date.autoformat.day'],
            1. / HOURS_PER_DAY: rcParams['date.autoformat.hour'],
            1. / (MINUTES_PER_DAY): rcParams['date.autoformat.minute'],
            1. / (SEC_PER_DAY): rcParams['date.autoformat.second'],
            1. / (MUSECONDS_PER_DAY): rcParams['date.autoformat.microsecond'],
            }


    The algorithm picks the key in the dictionary that is >= the
    current scale and uses that format string.  You can customize this
    dictionary by doing::


    >>> locator = AutoDateLocator()
    >>> formatter = AutoDateFormatter(locator)
    >>> formatter.scaled[1/(24.*60.)] = '%M:%S' # only show min and sec

    A custom :class:`~matplotlib.ticker.FuncFormatter` can also be used.
    The following example shows how to use a custom format function to strip
    trailing zeros from decimal seconds and adds the date to the first
    ticklabel::

        >>> def my_format_function(x, pos=None):
        ...     x = matplotlib.dates.num2date(x)
        ...     if pos == 0:
        ...         fmt = '%D %H:%M:%S.%f'
        ...     else:
        ...         fmt = '%H:%M:%S.%f'
        ...     label = x.strftime(fmt)
        ...     label = label.rstrip("0")
        ...     label = label.rstrip(".")
        ...     return label
        >>> from matplotlib.ticker import FuncFormatter
        >>> formatter.scaled[1/(24.*60.)] = FuncFormatter(my_format_function)
    """

    # This can be improved by providing some user-level direction on
    # how to choose the best format (precedence, etc...)

    # Perhaps a 'struct' that has a field for each time-type where a
    # zero would indicate "don't show" and a number would indicate
    # "show" with some sort of priority.  Same priorities could mean
    # show all with the same priority.

    # Or more simply, perhaps just a format string for each
    # possibility...

    def __init__(self, locator, tz=None, defaultfmt='%Y-%m-%d'):
        """
        Autoformat the date labels.  The default format is the one to use
        if none of the values in ``self.scaled`` are greater than the unit
        returned by ``locator._get_unit()``.
        """
        self._locator = locator
        self._tz = tz
        self.defaultfmt = defaultfmt
        self._formatter = DateFormatter(self.defaultfmt, tz)
        self.scaled = {DAYS_PER_YEAR: rcParams['date.autoformatter.year'],
                       DAYS_PER_MONTH: rcParams['date.autoformatter.month'],
                       1.0: rcParams['date.autoformatter.day'],
                       1. / HOURS_PER_DAY: rcParams['date.autoformatter.hour'],
                       1. / (MINUTES_PER_DAY):
                           rcParams['date.autoformatter.minute'],
                       1. / (SEC_PER_DAY):
                           rcParams['date.autoformatter.second'],
                       1. / (MUSECONDS_PER_DAY):
                           rcParams['date.autoformatter.microsecond']}

    def __call__(self, x, pos=None):
        locator_unit_scale = float(self._locator._get_unit())
        # Pick the first scale which is greater than the locator unit.
        fmt = next((fmt for scale, fmt in sorted(self.scaled.items())
                    if scale >= locator_unit_scale),
                   self.defaultfmt)

        if isinstance(fmt, six.string_types):
            self._formatter = DateFormatter(fmt, self._tz)
            result = self._formatter(x, pos)
        elif callable(fmt):
            result = fmt(x, pos)
        else:
            raise TypeError('Unexpected type passed to {0!r}.'.format(self))

        return result


class rrulewrapper(object):
    def __init__(self, freq, tzinfo=None, **kwargs):
        kwargs['freq'] = freq
        self._base_tzinfo = tzinfo

        self._update_rrule(**kwargs)

    def set(self, **kwargs):
        self._construct.update(kwargs)

        self._update_rrule(**self._construct)

    def _update_rrule(self, **kwargs):
        tzinfo = self._base_tzinfo

        # rrule does not play nicely with time zones - especially pytz time
        # zones, it's best to use naive zones and attach timezones once the
        # datetimes are returned
        if 'dtstart' in kwargs:
            dtstart = kwargs['dtstart']
            if dtstart.tzinfo is not None:
                if tzinfo is None:
                    tzinfo = dtstart.tzinfo
                else:
                    dtstart = dtstart.astimezone(tzinfo)

                kwargs['dtstart'] = dtstart.replace(tzinfo=None)

        if 'until' in kwargs:
            until = kwargs['until']
            if until.tzinfo is not None:
                if tzinfo is not None:
                    until = until.astimezone(tzinfo)
                else:
                    raise ValueError('until cannot be aware if dtstart '
                                     'is naive and tzinfo is None')

                kwargs['until'] = until.replace(tzinfo=None)

        self._construct = kwargs.copy()
        self._tzinfo = tzinfo
        self._rrule = rrule(**self._construct)

    def _attach_tzinfo(self, dt, tzinfo):
        # pytz zones are attached by "localizing" the datetime
        if hasattr(tzinfo, 'localize'):
            return tzinfo.localize(dt, is_dst=True)

        return dt.replace(tzinfo=tzinfo)

    def _aware_return_wrapper(self, f, returns_list=False):
        """Decorator function that allows rrule methods to handle tzinfo."""
        # This is only necessary if we're actually attaching a tzinfo
        if self._tzinfo is None:
            return f

        # All datetime arguments must be naive. If they are not naive, they are
        # converted to the _tzinfo zone before dropping the zone.
        def normalize_arg(arg):
            if isinstance(arg, datetime.datetime) and arg.tzinfo is not None:
                if arg.tzinfo is not self._tzinfo:
                    arg = arg.astimezone(self._tzinfo)

                return arg.replace(tzinfo=None)

            return arg

        def normalize_args(args, kwargs):
            args = tuple(normalize_arg(arg) for arg in args)
            kwargs = {kw: normalize_arg(arg) for kw, arg in kwargs.items()}

            return args, kwargs

        # There are two kinds of functions we care about - ones that return
        # dates and ones that return lists of dates.
        if not returns_list:
            def inner_func(*args, **kwargs):
                args, kwargs = normalize_args(args, kwargs)
                dt = f(*args, **kwargs)
                return self._attach_tzinfo(dt, self._tzinfo)
        else:
            def inner_func(*args, **kwargs):
                args, kwargs = normalize_args(args, kwargs)
                dts = f(*args, **kwargs)
                return [self._attach_tzinfo(dt, self._tzinfo) for dt in dts]

        return functools.wraps(f)(inner_func)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]

        f = getattr(self._rrule, name)

        if name in {'after', 'before'}:
            return self._aware_return_wrapper(f)
        elif name in {'xafter', 'xbefore', 'between'}:
            return self._aware_return_wrapper(f, returns_list=True)
        else:
            return f

    def __setstate__(self, state):
        self.__dict__.update(state)


class DateLocator(ticker.Locator):
    """
    Determines the tick locations when plotting dates.

    This class is subclassed by other Locators and
    is not meant to be used on its own.
    """
    hms0d = {'byhour': 0, 'byminute': 0, 'bysecond': 0}

    def __init__(self, tz=None):
        """
        *tz* is a :class:`tzinfo` instance.
        """
        if tz is None:
            tz = _get_rc_timezone()
        self.tz = tz

    def set_tzinfo(self, tz):
        """
        Set time zone info.
        """
        self.tz = tz

    def datalim_to_dt(self):
        """
        Convert axis data interval to datetime objects.
        """
        dmin, dmax = self.axis.get_data_interval()
        if dmin > dmax:
            dmin, dmax = dmax, dmin
        if dmin < 1:
            raise ValueError('datalim minimum {} is less than 1 and '
                             'is an invalid Matplotlib date value. This often '
                             'happens if you pass a non-datetime '
                             'value to an axis that has datetime units'
                             .format(dmin))
        return num2date(dmin, self.tz), num2date(dmax, self.tz)

    def viewlim_to_dt(self):
        """
        Converts the view interval to datetime objects.
        """
        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if vmin < 1:
            raise ValueError('view limit minimum {} is less than 1 and '
                             'is an invalid Matplotlib date value. This '
                             'often happens if you pass a non-datetime '
                             'value to an axis that has datetime units'
                             .format(vmin))
        return num2date(vmin, self.tz), num2date(vmax, self.tz)

    def _get_unit(self):
        """
        Return how many days a unit of the locator is; used for
        intelligent autoscaling.
        """
        return 1

    def _get_interval(self):
        """
        Return the number of units for each tick.
        """
        return 1

    def nonsingular(self, vmin, vmax):
        """
        Given the proposed upper and lower extent, adjust the range
        if it is too close to being singular (i.e. a range of ~0).

        """
        unit = self._get_unit()
        interval = self._get_interval()
        if abs(vmax - vmin) < 1e-6:
            vmin -= 2 * unit * interval
            vmax += 2 * unit * interval
        return vmin, vmax


class RRuleLocator(DateLocator):
    # use the dateutil rrule instance

    def __init__(self, o, tz=None):
        DateLocator.__init__(self, tz)
        self.rule = o

    def __call__(self):
        # if no data have been set, this will tank with a ValueError
        try:
            dmin, dmax = self.viewlim_to_dt()
        except ValueError:
            return []

        return self.tick_values(dmin, dmax)

    def tick_values(self, vmin, vmax):
        delta = relativedelta(vmax, vmin)

        # We need to cap at the endpoints of valid datetime
        try:
            start = vmin - delta
        except (ValueError, OverflowError):
            start = _from_ordinalf(1.0)

        try:
            stop = vmax + delta
        except (ValueError, OverflowError):
            # The magic number!
            stop = _from_ordinalf(3652059.9999999)

        self.rule.set(dtstart=start, until=stop)

        dates = self.rule.between(vmin, vmax, True)
        if len(dates) == 0:
            return date2num([vmin, vmax])
        return self.raise_if_exceeds(date2num(dates))

    def _get_unit(self):
        """
        Return how many days a unit of the locator is; used for
        intelligent autoscaling.
        """
        freq = self.rule._rrule._freq
        return self.get_unit_generic(freq)

    @staticmethod
    def get_unit_generic(freq):
        if freq == YEARLY:
            return DAYS_PER_YEAR
        elif freq == MONTHLY:
            return DAYS_PER_MONTH
        elif freq == WEEKLY:
            return DAYS_PER_WEEK
        elif freq == DAILY:
            return 1.0
        elif freq == HOURLY:
            return 1.0 / HOURS_PER_DAY
        elif freq == MINUTELY:
            return 1.0 / MINUTES_PER_DAY
        elif freq == SECONDLY:
            return 1.0 / SEC_PER_DAY
        else:
            # error
            return -1   # or should this just return '1'?

    def _get_interval(self):
        return self.rule._rrule._interval

    def autoscale(self):
        """
        Set the view limits to include the data range.
        """
        dmin, dmax = self.datalim_to_dt()
        delta = relativedelta(dmax, dmin)

        # We need to cap at the endpoints of valid datetime
        try:
            start = dmin - delta
        except ValueError:
            start = _from_ordinalf(1.0)

        try:
            stop = dmax + delta
        except ValueError:
            # The magic number!
            stop = _from_ordinalf(3652059.9999999)

        self.rule.set(dtstart=start, until=stop)
        dmin, dmax = self.datalim_to_dt()

        vmin = self.rule.before(dmin, True)
        if not vmin:
            vmin = dmin

        vmax = self.rule.after(dmax, True)
        if not vmax:
            vmax = dmax

        vmin = date2num(vmin)
        vmax = date2num(vmax)

        return self.nonsingular(vmin, vmax)


class AutoDateLocator(DateLocator):
    """
    On autoscale, this class picks the best
    :class:`DateLocator` to set the view limits and the tick
    locations.
    """
    def __init__(self, tz=None, minticks=5, maxticks=None,
                 interval_multiples=False):
        """
        *minticks* is the minimum number of ticks desired, which is used to
        select the type of ticking (yearly, monthly, etc.).

        *maxticks* is the maximum number of ticks desired, which controls
        any interval between ticks (ticking every other, every 3, etc.).
        For really fine-grained control, this can be a dictionary mapping
        individual rrule frequency constants (YEARLY, MONTHLY, etc.)
        to their own maximum number of ticks.  This can be used to keep
        the number of ticks appropriate to the format chosen in
        :class:`AutoDateFormatter`. Any frequency not specified in this
        dictionary is given a default value.

        *tz* is a :class:`tzinfo` instance.

        *interval_multiples* is a boolean that indicates whether ticks
        should be chosen to be multiple of the interval. This will lock
        ticks to 'nicer' locations. For example, this will force the
        ticks to be at hours 0,6,12,18 when hourly ticking is done at
        6 hour intervals.

        The AutoDateLocator has an interval dictionary that maps the
        frequency of the tick (a constant from dateutil.rrule) and a
        multiple allowed for that ticking.  The default looks like this::

          self.intervald = {
            YEARLY  : [1, 2, 4, 5, 10, 20, 40, 50, 100, 200, 400, 500,
                      1000, 2000, 4000, 5000, 10000],
            MONTHLY : [1, 2, 3, 4, 6],
            DAILY   : [1, 2, 3, 7, 14],
            HOURLY  : [1, 2, 3, 4, 6, 12],
            MINUTELY: [1, 5, 10, 15, 30],
            SECONDLY: [1, 5, 10, 15, 30],
            MICROSECONDLY: [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,
                           5000, 10000, 20000, 50000, 100000, 200000, 500000,
                           1000000],
            }

        The interval is used to specify multiples that are appropriate for
        the frequency of ticking. For instance, every 7 days is sensible
        for daily ticks, but for minutes/seconds, 15 or 30 make sense.
        You can customize this dictionary by doing::

          locator = AutoDateLocator()
          locator.intervald[HOURLY] = [3] # only show every 3 hours
        """
        DateLocator.__init__(self, tz)
        self._locator = YearLocator()
        self._freq = YEARLY
        self._freqs = [YEARLY, MONTHLY, DAILY, HOURLY, MINUTELY,
                       SECONDLY, MICROSECONDLY]
        self.minticks = minticks

        self.maxticks = {YEARLY: 11, MONTHLY: 12, DAILY: 11, HOURLY: 12,
                         MINUTELY: 11, SECONDLY: 11, MICROSECONDLY: 8}
        if maxticks is not None:
            try:
                self.maxticks.update(maxticks)
            except TypeError:
                # Assume we were given an integer. Use this as the maximum
                # number of ticks for every frequency and create a
                # dictionary for this
                self.maxticks = dict.fromkeys(self._freqs, maxticks)
        self.interval_multiples = interval_multiples
        self.intervald = {
            YEARLY:   [1, 2, 4, 5, 10, 20, 40, 50, 100, 200, 400, 500,
                       1000, 2000, 4000, 5000, 10000],
            MONTHLY:  [1, 2, 3, 4, 6],
            DAILY:    [1, 2, 3, 7, 14, 21],
            HOURLY:   [1, 2, 3, 4, 6, 12],
            MINUTELY: [1, 5, 10, 15, 30],
            SECONDLY: [1, 5, 10, 15, 30],
            MICROSECONDLY: [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,
                            5000, 10000, 20000, 50000, 100000, 200000, 500000,
                            1000000]}
        self._byranges = [None, range(1, 13), range(1, 32),
                          range(0, 24), range(0, 60), range(0, 60), None]

    def __call__(self):
        'Return the locations of the ticks'
        self.refresh()
        return self._locator()

    def tick_values(self, vmin, vmax):
        return self.get_locator(vmin, vmax).tick_values(vmin, vmax)

    def nonsingular(self, vmin, vmax):
        # whatever is thrown at us, we can scale the unit.
        # But default nonsingular date plots at an ~4 year period.
        if vmin == vmax:
            vmin = vmin - DAYS_PER_YEAR * 2
            vmax = vmax + DAYS_PER_YEAR * 2
        return vmin, vmax

    def set_axis(self, axis):
        DateLocator.set_axis(self, axis)
        self._locator.set_axis(axis)

    def refresh(self):
        'Refresh internal information based on current limits.'
        dmin, dmax = self.viewlim_to_dt()
        self._locator = self.get_locator(dmin, dmax)

    def _get_unit(self):
        if self._freq in [MICROSECONDLY]:
            return 1. / MUSECONDS_PER_DAY
        else:
            return RRuleLocator.get_unit_generic(self._freq)

    def autoscale(self):
        'Try to choose the view limits intelligently.'
        dmin, dmax = self.datalim_to_dt()
        self._locator = self.get_locator(dmin, dmax)
        return self._locator.autoscale()

    def get_locator(self, dmin, dmax):
        'Pick the best locator based on a distance.'
        delta = relativedelta(dmax, dmin)
        tdelta = dmax - dmin

        # take absolute difference
        if dmin > dmax:
            delta = -delta
            tdelta = -tdelta

        # The following uses a mix of calls to relativedelta and timedelta
        # methods because there is incomplete overlap in the functionality of
        # these similar functions, and it's best to avoid doing our own math
        # whenever possible.
        numYears = float(delta.years)
        numMonths = numYears * MONTHS_PER_YEAR + delta.months
        numDays = tdelta.days   # Avoids estimates of days/month, days/year
        numHours = numDays * HOURS_PER_DAY + delta.hours
        numMinutes = numHours * MIN_PER_HOUR + delta.minutes
        numSeconds = np.floor(tdelta.total_seconds())
        numMicroseconds = np.floor(tdelta.total_seconds() * 1e6)

        nums = [numYears, numMonths, numDays, numHours, numMinutes,
                numSeconds, numMicroseconds]

        use_rrule_locator = [True] * 6 + [False]

        # Default setting of bymonth, etc. to pass to rrule
        # [unused (for year), bymonth, bymonthday, byhour, byminute,
        #  bysecond, unused (for microseconds)]
        byranges = [None, 1, 1, 0, 0, 0, None]

        # Loop over all the frequencies and try to find one that gives at
        # least a minticks tick positions.  Once this is found, look for
        # an interval from an list specific to that frequency that gives no
        # more than maxticks tick positions. Also, set up some ranges
        # (bymonth, etc.) as appropriate to be passed to rrulewrapper.
        for i, (freq, num) in enumerate(zip(self._freqs, nums)):
            # If this particular frequency doesn't give enough ticks, continue
            if num < self.minticks:
                # Since we're not using this particular frequency, set
                # the corresponding by_ to None so the rrule can act as
                # appropriate
                byranges[i] = None
                continue

            # Find the first available interval that doesn't give too many
            # ticks
            for interval in self.intervald[freq]:
                if num <= interval * (self.maxticks[freq] - 1):
                    break
            else:
                # We went through the whole loop without breaking, default to
                # the last interval in the list and raise a warning
                warnings.warn('AutoDateLocator was unable to pick an '
                              'appropriate interval for this date range. '
                              'It may be necessary to add an interval value '
                              "to the AutoDateLocator's intervald dictionary."
                              ' Defaulting to {0}.'.format(interval))

            # Set some parameters as appropriate
            self._freq = freq

            if self._byranges[i] and self.interval_multiples:
                byranges[i] = self._byranges[i][::interval]
                interval = 1
            else:
                byranges[i] = self._byranges[i]

            break
        else:
            raise ValueError('No sensible date limit could be found in the '
                             'AutoDateLocator.')

        if (freq == YEARLY) and self.interval_multiples:
            locator = YearLocator(interval)
        elif use_rrule_locator[i]:
            _, bymonth, bymonthday, byhour, byminute, bysecond, _ = byranges
            rrule = rrulewrapper(self._freq, interval=interval,
                                 dtstart=dmin, until=dmax,
                                 bymonth=bymonth, bymonthday=bymonthday,
                                 byhour=byhour, byminute=byminute,
                                 bysecond=bysecond)

            locator = RRuleLocator(rrule, self.tz)
        else:
            locator = MicrosecondLocator(interval, tz=self.tz)
            if dmin.year > 20 and interval < 1000:
                _log.warn('Plotting microsecond time intervals is not'
                          ' well supported. Please see the'
                          ' MicrosecondLocator documentation'
                          ' for details.')

        locator.set_axis(self.axis)

        if self.axis is not None:
            locator.set_view_interval(*self.axis.get_view_interval())
            locator.set_data_interval(*self.axis.get_data_interval())
        return locator


class YearLocator(DateLocator):
    """
    Make ticks on a given day of each year that is a multiple of base.

    Examples::

      # Tick every year on Jan 1st
      locator = YearLocator()

      # Tick every 5 years on July 4th
      locator = YearLocator(5, month=7, day=4)
    """
    def __init__(self, base=1, month=1, day=1, tz=None):
        """
        Mark years that are multiple of base on a given month and day
        (default jan 1).
        """
        DateLocator.__init__(self, tz)
        self.base = ticker.Base(base)
        self.replaced = {'month':  month,
                         'day':    day,
                         'hour':   0,
                         'minute': 0,
                         'second': 0,
                         'tzinfo': tz
                         }

    def __call__(self):
        # if no data have been set, this will tank with a ValueError
        try:
            dmin, dmax = self.viewlim_to_dt()
        except ValueError:
            return []

        return self.tick_values(dmin, dmax)

    def tick_values(self, vmin, vmax):
        ymin = self.base.le(vmin.year)
        ymax = self.base.ge(vmax.year)

        ticks = [vmin.replace(year=ymin, **self.replaced)]
        while True:
            dt = ticks[-1]
            if dt.year >= ymax:
                return date2num(ticks)
            year = dt.year + self.base.get_base()
            ticks.append(dt.replace(year=year, **self.replaced))

    def autoscale(self):
        """
        Set the view limits to include the data range.
        """
        dmin, dmax = self.datalim_to_dt()

        ymin = self.base.le(dmin.year)
        ymax = self.base.ge(dmax.year)
        vmin = dmin.replace(year=ymin, **self.replaced)
        vmax = dmax.replace(year=ymax, **self.replaced)

        vmin = date2num(vmin)
        vmax = date2num(vmax)
        return self.nonsingular(vmin, vmax)


class MonthLocator(RRuleLocator):
    """
    Make ticks on occurrences of each month, e.g., 1, 3, 12.
    """
    def __init__(self, bymonth=None, bymonthday=1, interval=1, tz=None):
        """
        Mark every month in *bymonth*; *bymonth* can be an int or
        sequence.  Default is ``range(1,13)``, i.e. every month.

        *interval* is the interval between each iteration.  For
        example, if ``interval=2``, mark every second occurrence.
        """
        if bymonth is None:
            bymonth = range(1, 13)
        elif isinstance(bymonth, np.ndarray):
            # This fixes a bug in dateutil <= 2.3 which prevents the use of
            # numpy arrays in (among other things) the bymonthday, byweekday
            # and bymonth parameters.
            bymonth = [x.item() for x in bymonth.astype(int)]

        rule = rrulewrapper(MONTHLY, bymonth=bymonth, bymonthday=bymonthday,
                            interval=interval, **self.hms0d)
        RRuleLocator.__init__(self, rule, tz)


class WeekdayLocator(RRuleLocator):
    """
    Make ticks on occurrences of each weekday.
    """

    def __init__(self, byweekday=1, interval=1, tz=None):
        """
        Mark every weekday in *byweekday*; *byweekday* can be a number or
        sequence.

        Elements of *byweekday* must be one of MO, TU, WE, TH, FR, SA,
        SU, the constants from :mod:`dateutil.rrule`, which have been
        imported into the :mod:`matplotlib.dates` namespace.

        *interval* specifies the number of weeks to skip.  For example,
        ``interval=2`` plots every second week.
        """
        if isinstance(byweekday, np.ndarray):
            # This fixes a bug in dateutil <= 2.3 which prevents the use of
            # numpy arrays in (among other things) the bymonthday, byweekday
            # and bymonth parameters.
            [x.item() for x in byweekday.astype(int)]

        rule = rrulewrapper(DAILY, byweekday=byweekday,
                            interval=interval, **self.hms0d)
        RRuleLocator.__init__(self, rule, tz)


class DayLocator(RRuleLocator):
    """
    Make ticks on occurrences of each day of the month.  For example,
    1, 15, 30.
    """
    def __init__(self, bymonthday=None, interval=1, tz=None):
        """
        Mark every day in *bymonthday*; *bymonthday* can be an int or
        sequence.

        Default is to tick every day of the month: ``bymonthday=range(1,32)``
        """
        if not interval == int(interval) or interval < 1:
            raise ValueError("interval must be an integer greater than 0")
        if bymonthday is None:
            bymonthday = range(1, 32)
        elif isinstance(bymonthday, np.ndarray):
            # This fixes a bug in dateutil <= 2.3 which prevents the use of
            # numpy arrays in (among other things) the bymonthday, byweekday
            # and bymonth parameters.
            bymonthday = [x.item() for x in bymonthday.astype(int)]

        rule = rrulewrapper(DAILY, bymonthday=bymonthday,
                            interval=interval, **self.hms0d)
        RRuleLocator.__init__(self, rule, tz)


class HourLocator(RRuleLocator):
    """
    Make ticks on occurrences of each hour.
    """
    def __init__(self, byhour=None, interval=1, tz=None):
        """
        Mark every hour in *byhour*; *byhour* can be an int or sequence.
        Default is to tick every hour: ``byhour=range(24)``

        *interval* is the interval between each iteration.  For
        example, if ``interval=2``, mark every second occurrence.
        """
        if byhour is None:
            byhour = range(24)

        rule = rrulewrapper(HOURLY, byhour=byhour, interval=interval,
                            byminute=0, bysecond=0)
        RRuleLocator.__init__(self, rule, tz)


class MinuteLocator(RRuleLocator):
    """
    Make ticks on occurrences of each minute.
    """
    def __init__(self, byminute=None, interval=1, tz=None):
        """
        Mark every minute in *byminute*; *byminute* can be an int or
        sequence.  Default is to tick every minute: ``byminute=range(60)``

        *interval* is the interval between each iteration.  For
        example, if ``interval=2``, mark every second occurrence.
        """
        if byminute is None:
            byminute = range(60)

        rule = rrulewrapper(MINUTELY, byminute=byminute, interval=interval,
                            bysecond=0)
        RRuleLocator.__init__(self, rule, tz)


class SecondLocator(RRuleLocator):
    """
    Make ticks on occurrences of each second.
    """
    def __init__(self, bysecond=None, interval=1, tz=None):
        """
        Mark every second in *bysecond*; *bysecond* can be an int or
        sequence.  Default is to tick every second: ``bysecond = range(60)``

        *interval* is the interval between each iteration.  For
        example, if ``interval=2``, mark every second occurrence.

        """
        if bysecond is None:
            bysecond = range(60)

        rule = rrulewrapper(SECONDLY, bysecond=bysecond, interval=interval)
        RRuleLocator.__init__(self, rule, tz)


class MicrosecondLocator(DateLocator):
    """
    Make ticks on regular intervals of one or more microsecond(s).

    .. note::

        Due to the floating point representation of time in days since
        0001-01-01 UTC (plus 1), plotting data with microsecond time
        resolution does not work well with current dates.

        If you want microsecond resolution time plots, it is strongly
        recommended to use floating point seconds, not datetime-like
        time representation.

        If you really must use datetime.datetime() or similar and still
        need microsecond precision, your only chance is to use very
        early years; using year 0001 is recommended.

    """
    def __init__(self, interval=1, tz=None):
        """
        *interval* is the interval between each iteration.  For
        example, if ``interval=2``, mark every second microsecond.

        """
        self._interval = interval
        self._wrapped_locator = ticker.MultipleLocator(interval)
        self.tz = tz

    def set_axis(self, axis):
        self._wrapped_locator.set_axis(axis)
        return DateLocator.set_axis(self, axis)

    def set_view_interval(self, vmin, vmax):
        self._wrapped_locator.set_view_interval(vmin, vmax)
        return DateLocator.set_view_interval(self, vmin, vmax)

    def set_data_interval(self, vmin, vmax):
        self._wrapped_locator.set_data_interval(vmin, vmax)
        return DateLocator.set_data_interval(self, vmin, vmax)

    def __call__(self):
        # if no data have been set, this will tank with a ValueError
        try:
            dmin, dmax = self.viewlim_to_dt()
        except ValueError:
            return []

        return self.tick_values(dmin, dmax)

    def tick_values(self, vmin, vmax):
        nmin, nmax = date2num((vmin, vmax))
        nmin *= MUSECONDS_PER_DAY
        nmax *= MUSECONDS_PER_DAY
        ticks = self._wrapped_locator.tick_values(nmin, nmax)
        ticks = [tick / MUSECONDS_PER_DAY for tick in ticks]
        return ticks

    def _get_unit(self):
        """
        Return how many days a unit of the locator is; used for
        intelligent autoscaling.
        """
        return 1. / MUSECONDS_PER_DAY

    def _get_interval(self):
        """
        Return the number of units for each tick.
        """
        return self._interval


def _close_to_dt(d1, d2, epsilon=5):
    """
    Assert that datetimes *d1* and *d2* are within *epsilon* microseconds.
    """
    delta = d2 - d1
    mus = abs(delta.total_seconds() * 1e6)
    assert mus < epsilon


def _close_to_num(o1, o2, epsilon=5):
    """
    Assert that float ordinals *o1* and *o2* are within *epsilon*
    microseconds.
    """
    delta = abs((o2 - o1) * MUSECONDS_PER_DAY)
    assert delta < epsilon


def epoch2num(e):
    """
    Convert an epoch or sequence of epochs to the new date format,
    that is days since 0001.
    """
    return EPOCH_OFFSET + np.asarray(e) / SEC_PER_DAY


def num2epoch(d):
    """
    Convert days since 0001 to epoch.  *d* can be a number or sequence.
    """
    return (np.asarray(d) - EPOCH_OFFSET) * SEC_PER_DAY


def mx2num(mxdates):
    """
    Convert mx :class:`datetime` instance (or sequence of mx
    instances) to the new date format.
    """
    scalar = False
    if not cbook.iterable(mxdates):
        scalar = True
        mxdates = [mxdates]
    ret = epoch2num([m.ticks() for m in mxdates])
    if scalar:
        return ret[0]
    else:
        return ret


def date_ticker_factory(span, tz=None, numticks=5):
    """
    Create a date locator with *numticks* (approx) and a date formatter
    for *span* in days.  Return value is (locator, formatter).
    """

    if span == 0:
        span = 1 / HOURS_PER_DAY

    mins = span * MINUTES_PER_DAY
    hrs = span * HOURS_PER_DAY
    days = span
    wks = span / DAYS_PER_WEEK
    months = span / DAYS_PER_MONTH      # Approx
    years = span / DAYS_PER_YEAR        # Approx

    if years > numticks:
        locator = YearLocator(int(years / numticks), tz=tz)  # define
        fmt = '%Y'
    elif months > numticks:
        locator = MonthLocator(tz=tz)
        fmt = '%b %Y'
    elif wks > numticks:
        locator = WeekdayLocator(tz=tz)
        fmt = '%a, %b %d'
    elif days > numticks:
        locator = DayLocator(interval=int(math.ceil(days / numticks)), tz=tz)
        fmt = '%b %d'
    elif hrs > numticks:
        locator = HourLocator(interval=int(math.ceil(hrs / numticks)), tz=tz)
        fmt = '%H:%M\n%b %d'
    elif mins > numticks:
        locator = MinuteLocator(interval=int(math.ceil(mins / numticks)),
                                tz=tz)
        fmt = '%H:%M:%S'
    else:
        locator = MinuteLocator(tz=tz)
        fmt = '%H:%M:%S'

    formatter = DateFormatter(fmt, tz=tz)
    return locator, formatter


def seconds(s):
    """
    Return seconds as days.
    """
    return s / SEC_PER_DAY


def minutes(m):
    """
    Return minutes as days.
    """
    return m / MINUTES_PER_DAY


def hours(h):
    """
    Return hours as days.
    """
    return h / HOURS_PER_DAY


def weeks(w):
    """
    Return weeks as days.
    """
    return w * DAYS_PER_WEEK


class DateConverter(units.ConversionInterface):
    """
    Converter for datetime.date and datetime.datetime data,
    or for date/time data represented as it would be converted
    by :func:`date2num`.

    The 'unit' tag for such data is None or a tzinfo instance.
    """

    @staticmethod
    def axisinfo(unit, axis):
        """
        Return the :class:`~matplotlib.units.AxisInfo` for *unit*.

        *unit* is a tzinfo instance or None.
        The *axis* argument is required but not used.
        """
        tz = unit

        majloc = AutoDateLocator(tz=tz)
        majfmt = AutoDateFormatter(majloc, tz=tz)
        datemin = datetime.date(2000, 1, 1)
        datemax = datetime.date(2010, 1, 1)

        return units.AxisInfo(majloc=majloc, majfmt=majfmt, label='',
                              default_limits=(datemin, datemax))

    @staticmethod
    def convert(value, unit, axis):
        """
        If *value* is not already a number or sequence of numbers,
        convert it with :func:`date2num`.

        The *unit* and *axis* arguments are not used.
        """
        return date2num(value)

    @staticmethod
    def default_units(x, axis):
        """
        Return the tzinfo instance of *x* or of its first element, or None
        """
        if isinstance(x, np.ndarray):
            x = x.ravel()

        try:
            x = cbook.safe_first_element(x)
        except (TypeError, StopIteration):
            pass

        try:
            return x.tzinfo
        except AttributeError:
            pass
        return None


units.registry[np.datetime64] = DateConverter()
units.registry[datetime.date] = DateConverter()
units.registry[datetime.datetime] = DateConverter()
