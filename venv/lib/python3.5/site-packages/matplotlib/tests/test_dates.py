from __future__ import absolute_import, division, print_function

from six.moves import map


import datetime
import dateutil
import tempfile

import numpy as np
import pytest
import pytz

try:
    # mock in python 3.3+
    from unittest import mock
except ImportError:
    import mock

from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def test_date_numpyx():
    # test that numpy dates work properly...
    base = datetime.datetime(2017, 1, 1)
    time = [base + datetime.timedelta(days=x) for x in range(0, 3)]
    timenp = np.array(time, dtype='datetime64[ns]')
    data = np.array([0., 2., 1.])
    fig = plt.figure(figsize=(10, 2))
    ax = fig.add_subplot(1, 1, 1)
    h, = ax.plot(time, data)
    hnp, = ax.plot(timenp, data)
    assert np.array_equal(h.get_xdata(orig=False), hnp.get_xdata(orig=False))
    fig = plt.figure(figsize=(10, 2))
    ax = fig.add_subplot(1, 1, 1)
    h, = ax.plot(data, time)
    hnp, = ax.plot(data, timenp)
    assert np.array_equal(h.get_ydata(orig=False), hnp.get_ydata(orig=False))


@pytest.mark.parametrize('t0', [datetime.datetime(2017, 1, 1, 0, 1, 1),

                                [datetime.datetime(2017, 1, 1, 0, 1, 1),
                                 datetime.datetime(2017, 1, 1, 1, 1, 1)],

                                [[datetime.datetime(2017, 1, 1, 0, 1, 1),
                                  datetime.datetime(2017, 1, 1, 1, 1, 1)],
                                 [datetime.datetime(2017, 1, 1, 2, 1, 1),
                                  datetime.datetime(2017, 1, 1, 3, 1, 1)]]])
@pytest.mark.parametrize('dtype', ['datetime64[s]',
                                    'datetime64[us]',
                                    'datetime64[ms]',
                                    'datetime64[ns]'])
def test_date_date2num_numpy(t0, dtype):
    time = mdates.date2num(t0)
    tnp = np.array(t0, dtype=dtype)
    nptime = mdates.date2num(tnp)
    assert np.array_equal(time, nptime)


@pytest.mark.parametrize('dtype', ['datetime64[s]',
                                    'datetime64[us]',
                                    'datetime64[ms]',
                                    'datetime64[ns]'])
def test_date2num_NaT(dtype):
    t0 = datetime.datetime(2017, 1, 1, 0, 1, 1)
    tmpl = [mdates.date2num(t0), np.nan]
    tnp = np.array([t0, 'NaT'], dtype=dtype)
    nptime = mdates.date2num(tnp)
    np.testing.assert_array_equal(tmpl, nptime)


@pytest.mark.parametrize('units', ['s', 'ms', 'us', 'ns'])
def test_date2num_NaT_scalar(units):
    tmpl = mdates.date2num(np.datetime64('NaT', units))
    assert np.isnan(tmpl)


@image_comparison(baseline_images=['date_empty'], extensions=['png'])
def test_date_empty():
    # make sure mpl does the right thing when told to plot dates even
    # if no date data has been presented, cf
    # http://sourceforge.net/tracker/?func=detail&aid=2850075&group_id=80706&atid=560720
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.xaxis_date()


@image_comparison(baseline_images=['date_axhspan'], extensions=['png'])
def test_date_axhspan():
    # test ax hspan with date inputs
    t0 = datetime.datetime(2009, 1, 20)
    tf = datetime.datetime(2009, 1, 21)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axhspan(t0, tf, facecolor="blue", alpha=0.25)
    ax.set_ylim(t0 - datetime.timedelta(days=5),
                tf + datetime.timedelta(days=5))
    fig.subplots_adjust(left=0.25)


@image_comparison(baseline_images=['date_axvspan'], extensions=['png'])
def test_date_axvspan():
    # test ax hspan with date inputs
    t0 = datetime.datetime(2000, 1, 20)
    tf = datetime.datetime(2010, 1, 21)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axvspan(t0, tf, facecolor="blue", alpha=0.25)
    ax.set_xlim(t0 - datetime.timedelta(days=720),
                tf + datetime.timedelta(days=720))
    fig.autofmt_xdate()


@image_comparison(baseline_images=['date_axhline'],
                  extensions=['png'])
def test_date_axhline():
    # test ax hline with date inputs
    t0 = datetime.datetime(2009, 1, 20)
    tf = datetime.datetime(2009, 1, 31)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axhline(t0, color="blue", lw=3)
    ax.set_ylim(t0 - datetime.timedelta(days=5),
                tf + datetime.timedelta(days=5))
    fig.subplots_adjust(left=0.25)


@image_comparison(baseline_images=['date_axvline'],
                  extensions=['png'])
def test_date_axvline():
    # test ax hline with date inputs
    t0 = datetime.datetime(2000, 1, 20)
    tf = datetime.datetime(2000, 1, 21)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axvline(t0, color="red", lw=3)
    ax.set_xlim(t0 - datetime.timedelta(days=5),
                tf + datetime.timedelta(days=5))
    fig.autofmt_xdate()


def test_too_many_date_ticks():
    # Attempt to test SF 2715172, see
    # https://sourceforge.net/tracker/?func=detail&aid=2715172&group_id=80706&atid=560720
    # setting equal datetimes triggers and expander call in
    # transforms.nonsingular which results in too many ticks in the
    # DayLocator.  This should trigger a Locator.MAXTICKS RuntimeError
    t0 = datetime.datetime(2000, 1, 20)
    tf = datetime.datetime(2000, 1, 20)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    with pytest.warns(UserWarning) as rec:
        ax.set_xlim((t0, tf), auto=True)
        assert len(rec) == 1
        assert 'Attempting to set identical left==right' in str(rec[0].message)
    ax.plot([], [])
    ax.xaxis.set_major_locator(mdates.DayLocator())
    with pytest.raises(RuntimeError):
        fig.savefig('junk.png')


@image_comparison(baseline_images=['RRuleLocator_bounds'], extensions=['png'])
def test_RRuleLocator():
    import matplotlib.testing.jpl_units as units
    units.register()

    # This will cause the RRuleLocator to go out of bounds when it tries
    # to add padding to the limits, so we make sure it caps at the correct
    # boundary values.
    t0 = datetime.datetime(1000, 1, 1)
    tf = datetime.datetime(6000, 1, 1)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_autoscale_on(True)
    ax.plot([t0, tf], [0.0, 1.0], marker='o')

    rrule = mdates.rrulewrapper(dateutil.rrule.YEARLY, interval=500)
    locator = mdates.RRuleLocator(rrule)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))

    ax.autoscale_view()
    fig.autofmt_xdate()


def test_RRuleLocator_dayrange():
    loc = mdates.DayLocator()
    x1 = datetime.datetime(year=1, month=1, day=1, tzinfo=pytz.UTC)
    y1 = datetime.datetime(year=1, month=1, day=16, tzinfo=pytz.UTC)
    loc.tick_values(x1, y1)
    # On success, no overflow error shall be thrown


@image_comparison(baseline_images=['DateFormatter_fractionalSeconds'],
                  extensions=['png'])
def test_DateFormatter():
    import matplotlib.testing.jpl_units as units
    units.register()

    # Lets make sure that DateFormatter will allow us to have tick marks
    # at intervals of fractional seconds.

    t0 = datetime.datetime(2001, 1, 1, 0, 0, 0)
    tf = datetime.datetime(2001, 1, 1, 0, 0, 1)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_autoscale_on(True)
    ax.plot([t0, tf], [0.0, 1.0], marker='o')

    # rrule = mpldates.rrulewrapper( dateutil.rrule.YEARLY, interval=500 )
    # locator = mpldates.RRuleLocator( rrule )
    # ax.xaxis.set_major_locator( locator )
    # ax.xaxis.set_major_formatter( mpldates.AutoDateFormatter(locator) )

    ax.autoscale_view()
    fig.autofmt_xdate()


def test_date_formatter_strftime():
    """
    Tests that DateFormatter matches datetime.strftime,
    check microseconds for years before 1900 for bug #3179
    as well as a few related issues for years before 1900.
    """
    def test_strftime_fields(dt):
        """For datetime object dt, check DateFormatter fields"""
        # Note: the last couple of %%s are to check multiple %s are handled
        # properly; %% should get replaced by %.
        formatter = mdates.DateFormatter("%w %d %m %y %Y %H %I %M %S %%%f %%x")
        # Compute date fields without using datetime.strftime,
        # since datetime.strftime does not work before year 1900
        formatted_date_str = (
            "{weekday} {day:02d} {month:02d} {year:02d} {full_year:04d} "
            "{hour24:02d} {hour12:02d} {minute:02d} {second:02d} "
            "%{microsecond:06d} %x"
            .format(
                weekday=str((dt.weekday() + 1) % 7),
                day=dt.day,
                month=dt.month,
                year=dt.year % 100,
                full_year=dt.year,
                hour24=dt.hour,
                hour12=((dt.hour-1) % 12) + 1,
                minute=dt.minute,
                second=dt.second,
                microsecond=dt.microsecond))
        assert formatter.strftime(dt) == formatted_date_str

        try:
            # Test strftime("%x") with the current locale.
            import locale  # Might not exist on some platforms, such as Windows
            locale_formatter = mdates.DateFormatter("%x")
            locale_d_fmt = locale.nl_langinfo(locale.D_FMT)
            expanded_formatter = mdates.DateFormatter(locale_d_fmt)
            assert locale_formatter.strftime(dt) == \
                expanded_formatter.strftime(dt)
        except (ImportError, AttributeError):
            pass

    for year in range(1, 3000, 71):
        # Iterate through random set of years
        test_strftime_fields(datetime.datetime(year, 1, 1))
        test_strftime_fields(datetime.datetime(year, 2, 3, 4, 5, 6, 12345))


def test_date_formatter_callable():
    scale = -11
    locator = mock.Mock(_get_unit=mock.Mock(return_value=scale))
    callable_formatting_function = (lambda dates, _:
                                    [dt.strftime('%d-%m//%Y') for dt in dates])

    formatter = mdates.AutoDateFormatter(locator)
    formatter.scaled[-10] = callable_formatting_function
    assert formatter([datetime.datetime(2014, 12, 25)]) == ['25-12//2014']


def test_drange():
    """
    This test should check if drange works as expected, and if all the
    rounding errors are fixed
    """
    start = datetime.datetime(2011, 1, 1, tzinfo=mdates.UTC)
    end = datetime.datetime(2011, 1, 2, tzinfo=mdates.UTC)
    delta = datetime.timedelta(hours=1)
    # We expect 24 values in drange(start, end, delta), because drange returns
    # dates from an half open interval [start, end)
    assert len(mdates.drange(start, end, delta)) == 24

    # if end is a little bit later, we expect the range to contain one element
    # more
    end = end + datetime.timedelta(microseconds=1)
    assert len(mdates.drange(start, end, delta)) == 25

    # reset end
    end = datetime.datetime(2011, 1, 2, tzinfo=mdates.UTC)

    # and tst drange with "complicated" floats:
    # 4 hours = 1/6 day, this is an "dangerous" float
    delta = datetime.timedelta(hours=4)
    daterange = mdates.drange(start, end, delta)
    assert len(daterange) == 6
    assert mdates.num2date(daterange[-1]) == (end - delta)


def test_empty_date_with_year_formatter():
    # exposes sf bug 2861426:
    # https://sourceforge.net/tracker/?func=detail&aid=2861426&group_id=80706&atid=560720

    # update: I am no longer believe this is a bug, as I commented on
    # the tracker.  The question is now: what to do with this test

    import matplotlib.dates as dates

    fig = plt.figure()
    ax = fig.add_subplot(111)

    yearFmt = dates.DateFormatter('%Y')
    ax.xaxis.set_major_formatter(yearFmt)

    with tempfile.TemporaryFile() as fh:
        with pytest.raises(ValueError):
            fig.savefig(fh)


def test_auto_date_locator():
    def _create_auto_date_locator(date1, date2):
        locator = mdates.AutoDateLocator()
        locator.create_dummy_axis()
        locator.set_view_interval(mdates.date2num(date1),
                                  mdates.date2num(date2))
        return locator

    d1 = datetime.datetime(1990, 1, 1)
    results = ([datetime.timedelta(weeks=52 * 200),
                ['1990-01-01 00:00:00+00:00', '2010-01-01 00:00:00+00:00',
                 '2030-01-01 00:00:00+00:00', '2050-01-01 00:00:00+00:00',
                 '2070-01-01 00:00:00+00:00', '2090-01-01 00:00:00+00:00',
                 '2110-01-01 00:00:00+00:00', '2130-01-01 00:00:00+00:00',
                 '2150-01-01 00:00:00+00:00', '2170-01-01 00:00:00+00:00']
                ],
               [datetime.timedelta(weeks=52),
                ['1990-01-01 00:00:00+00:00', '1990-02-01 00:00:00+00:00',
                 '1990-03-01 00:00:00+00:00', '1990-04-01 00:00:00+00:00',
                 '1990-05-01 00:00:00+00:00', '1990-06-01 00:00:00+00:00',
                 '1990-07-01 00:00:00+00:00', '1990-08-01 00:00:00+00:00',
                 '1990-09-01 00:00:00+00:00', '1990-10-01 00:00:00+00:00',
                 '1990-11-01 00:00:00+00:00', '1990-12-01 00:00:00+00:00']
                ],
               [datetime.timedelta(days=141),
                ['1990-01-05 00:00:00+00:00', '1990-01-26 00:00:00+00:00',
                 '1990-02-16 00:00:00+00:00', '1990-03-09 00:00:00+00:00',
                 '1990-03-30 00:00:00+00:00', '1990-04-20 00:00:00+00:00',
                 '1990-05-11 00:00:00+00:00']
                ],
               [datetime.timedelta(days=40),
                ['1990-01-03 00:00:00+00:00', '1990-01-10 00:00:00+00:00',
                 '1990-01-17 00:00:00+00:00', '1990-01-24 00:00:00+00:00',
                 '1990-01-31 00:00:00+00:00', '1990-02-07 00:00:00+00:00']
                ],
               [datetime.timedelta(hours=40),
                ['1990-01-01 00:00:00+00:00', '1990-01-01 04:00:00+00:00',
                 '1990-01-01 08:00:00+00:00', '1990-01-01 12:00:00+00:00',
                 '1990-01-01 16:00:00+00:00', '1990-01-01 20:00:00+00:00',
                 '1990-01-02 00:00:00+00:00', '1990-01-02 04:00:00+00:00',
                 '1990-01-02 08:00:00+00:00', '1990-01-02 12:00:00+00:00',
                 '1990-01-02 16:00:00+00:00']
                ],
               [datetime.timedelta(minutes=20),
                ['1990-01-01 00:00:00+00:00', '1990-01-01 00:05:00+00:00',
                 '1990-01-01 00:10:00+00:00', '1990-01-01 00:15:00+00:00',
                 '1990-01-01 00:20:00+00:00']
                ],
               [datetime.timedelta(seconds=40),
                ['1990-01-01 00:00:00+00:00', '1990-01-01 00:00:05+00:00',
                 '1990-01-01 00:00:10+00:00', '1990-01-01 00:00:15+00:00',
                 '1990-01-01 00:00:20+00:00', '1990-01-01 00:00:25+00:00',
                 '1990-01-01 00:00:30+00:00', '1990-01-01 00:00:35+00:00',
                 '1990-01-01 00:00:40+00:00']
                ],
               [datetime.timedelta(microseconds=1500),
                ['1989-12-31 23:59:59.999500+00:00',
                 '1990-01-01 00:00:00+00:00',
                 '1990-01-01 00:00:00.000500+00:00',
                 '1990-01-01 00:00:00.001000+00:00',
                 '1990-01-01 00:00:00.001500+00:00']
                ],
               )

    for t_delta, expected in results:
        d2 = d1 + t_delta
        locator = _create_auto_date_locator(d1, d2)
        assert list(map(str, mdates.num2date(locator()))) == expected


def test_auto_date_locator_intmult():
    def _create_auto_date_locator(date1, date2):
        locator = mdates.AutoDateLocator(interval_multiples=True)
        locator.create_dummy_axis()
        locator.set_view_interval(mdates.date2num(date1),
                                  mdates.date2num(date2))
        return locator

    d1 = datetime.datetime(1997, 1, 1)
    results = ([datetime.timedelta(weeks=52 * 200),
                ['1980-01-01 00:00:00+00:00', '2000-01-01 00:00:00+00:00',
                 '2020-01-01 00:00:00+00:00', '2040-01-01 00:00:00+00:00',
                 '2060-01-01 00:00:00+00:00', '2080-01-01 00:00:00+00:00',
                 '2100-01-01 00:00:00+00:00', '2120-01-01 00:00:00+00:00',
                 '2140-01-01 00:00:00+00:00', '2160-01-01 00:00:00+00:00',
                 '2180-01-01 00:00:00+00:00', '2200-01-01 00:00:00+00:00']
                ],
               [datetime.timedelta(weeks=52),
                ['1997-01-01 00:00:00+00:00', '1997-02-01 00:00:00+00:00',
                 '1997-03-01 00:00:00+00:00', '1997-04-01 00:00:00+00:00',
                 '1997-05-01 00:00:00+00:00', '1997-06-01 00:00:00+00:00',
                 '1997-07-01 00:00:00+00:00', '1997-08-01 00:00:00+00:00',
                 '1997-09-01 00:00:00+00:00', '1997-10-01 00:00:00+00:00',
                 '1997-11-01 00:00:00+00:00', '1997-12-01 00:00:00+00:00']
                ],
               [datetime.timedelta(days=141),
                ['1997-01-01 00:00:00+00:00', '1997-01-22 00:00:00+00:00',
                 '1997-02-01 00:00:00+00:00', '1997-02-22 00:00:00+00:00',
                 '1997-03-01 00:00:00+00:00', '1997-03-22 00:00:00+00:00',
                 '1997-04-01 00:00:00+00:00', '1997-04-22 00:00:00+00:00',
                 '1997-05-01 00:00:00+00:00', '1997-05-22 00:00:00+00:00']
                ],
               [datetime.timedelta(days=40),
                ['1997-01-01 00:00:00+00:00', '1997-01-08 00:00:00+00:00',
                 '1997-01-15 00:00:00+00:00', '1997-01-22 00:00:00+00:00',
                 '1997-01-29 00:00:00+00:00', '1997-02-01 00:00:00+00:00',
                 '1997-02-08 00:00:00+00:00']
                ],
               [datetime.timedelta(hours=40),
                ['1997-01-01 00:00:00+00:00', '1997-01-01 04:00:00+00:00',
                 '1997-01-01 08:00:00+00:00', '1997-01-01 12:00:00+00:00',
                 '1997-01-01 16:00:00+00:00', '1997-01-01 20:00:00+00:00',
                 '1997-01-02 00:00:00+00:00', '1997-01-02 04:00:00+00:00',
                 '1997-01-02 08:00:00+00:00', '1997-01-02 12:00:00+00:00',
                 '1997-01-02 16:00:00+00:00']
                ],
               [datetime.timedelta(minutes=20),
                ['1997-01-01 00:00:00+00:00', '1997-01-01 00:05:00+00:00',
                 '1997-01-01 00:10:00+00:00', '1997-01-01 00:15:00+00:00',
                 '1997-01-01 00:20:00+00:00']
                ],
               [datetime.timedelta(seconds=40),
                ['1997-01-01 00:00:00+00:00', '1997-01-01 00:00:05+00:00',
                 '1997-01-01 00:00:10+00:00', '1997-01-01 00:00:15+00:00',
                 '1997-01-01 00:00:20+00:00', '1997-01-01 00:00:25+00:00',
                 '1997-01-01 00:00:30+00:00', '1997-01-01 00:00:35+00:00',
                 '1997-01-01 00:00:40+00:00']
                ],
               [datetime.timedelta(microseconds=1500),
                ['1996-12-31 23:59:59.999500+00:00',
                 '1997-01-01 00:00:00+00:00',
                 '1997-01-01 00:00:00.000500+00:00',
                 '1997-01-01 00:00:00.001000+00:00',
                 '1997-01-01 00:00:00.001500+00:00']
                ],
               )

    for t_delta, expected in results:
        d2 = d1 + t_delta
        locator = _create_auto_date_locator(d1, d2)
        assert list(map(str, mdates.num2date(locator()))) == expected


@image_comparison(baseline_images=['date_inverted_limit'],
                  extensions=['png'])
def test_date_inverted_limit():
    # test ax hline with date inputs
    t0 = datetime.datetime(2009, 1, 20)
    tf = datetime.datetime(2009, 1, 31)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axhline(t0, color="blue", lw=3)
    ax.set_ylim(t0 - datetime.timedelta(days=5),
                tf + datetime.timedelta(days=5))
    ax.invert_yaxis()
    fig.subplots_adjust(left=0.25)


def _test_date2num_dst(date_range, tz_convert):
    # Timezones
    BRUSSELS = pytz.timezone('Europe/Brussels')
    UTC = pytz.UTC

    # Create a list of timezone-aware datetime objects in UTC
    # Interval is 0b0.0000011 days, to prevent float rounding issues
    dtstart = datetime.datetime(2014, 3, 30, 0, 0, tzinfo=UTC)
    interval = datetime.timedelta(minutes=33, seconds=45)
    interval_days = 0.0234375   # 2025 / 86400 seconds
    N = 8

    dt_utc = date_range(start=dtstart, freq=interval, periods=N)
    dt_bxl = tz_convert(dt_utc, BRUSSELS)

    expected_ordinalf = [735322.0 + (i * interval_days) for i in range(N)]
    actual_ordinalf = list(mdates.date2num(dt_bxl))

    assert actual_ordinalf == expected_ordinalf


def test_date2num_dst():
    # Test for github issue #3896, but in date2num around DST transitions
    # with a timezone-aware pandas date_range object.

    class dt_tzaware(datetime.datetime):
        """
        This bug specifically occurs because of the normalization behavior of
        pandas Timestamp objects, so in order to replicate it, we need a
        datetime-like object that applies timezone normalization after
        subtraction.
        """
        def __sub__(self, other):
            r = super(dt_tzaware, self).__sub__(other)
            tzinfo = getattr(r, 'tzinfo', None)

            if tzinfo is not None:
                localizer = getattr(tzinfo, 'normalize', None)
                if localizer is not None:
                    r = tzinfo.normalize(r)

            if isinstance(r, datetime.datetime):
                r = self.mk_tzaware(r)

            return r

        def __add__(self, other):
            return self.mk_tzaware(super(dt_tzaware, self).__add__(other))

        def astimezone(self, tzinfo):
            dt = super(dt_tzaware, self).astimezone(tzinfo)
            return self.mk_tzaware(dt)

        @classmethod
        def mk_tzaware(cls, datetime_obj):
            kwargs = {}
            attrs = ('year',
                     'month',
                     'day',
                     'hour',
                     'minute',
                     'second',
                     'microsecond',
                     'tzinfo')

            for attr in attrs:
                val = getattr(datetime_obj, attr, None)
                if val is not None:
                    kwargs[attr] = val

            return cls(**kwargs)

    # Define a date_range function similar to pandas.date_range
    def date_range(start, freq, periods):
        dtstart = dt_tzaware.mk_tzaware(start)

        return [dtstart + (i * freq) for i in range(periods)]

    # Define a tz_convert function that converts a list to a new time zone.
    def tz_convert(dt_list, tzinfo):
        return [d.astimezone(tzinfo) for d in dt_list]

    _test_date2num_dst(date_range, tz_convert)


def test_date2num_dst_pandas(pd):
    # Test for github issue #3896, but in date2num around DST transitions
    # with a timezone-aware pandas date_range object.

    def tz_convert(*args):
        return pd.DatetimeIndex.tz_convert(*args).astype(object)

    _test_date2num_dst(pd.date_range, tz_convert)


@pytest.mark.parametrize("attach_tz, get_tz", [
    (lambda dt, zi: zi.localize(dt), lambda n: pytz.timezone(n)),
    (lambda dt, zi: dt.replace(tzinfo=zi), lambda n: dateutil.tz.gettz(n))])
def test_rrulewrapper(attach_tz, get_tz):
    SYD = get_tz('Australia/Sydney')

    dtstart = attach_tz(datetime.datetime(2017, 4, 1, 0), SYD)
    dtend = attach_tz(datetime.datetime(2017, 4, 4, 0), SYD)

    rule = mdates.rrulewrapper(freq=dateutil.rrule.DAILY, dtstart=dtstart)

    act = rule.between(dtstart, dtend)
    exp = [datetime.datetime(2017, 4, 1, 13, tzinfo=dateutil.tz.tzutc()),
           datetime.datetime(2017, 4, 2, 14, tzinfo=dateutil.tz.tzutc())]

    assert act == exp


def test_DayLocator():
    with pytest.raises(ValueError):
        mdates.DayLocator(interval=-1)
    with pytest.raises(ValueError):
        mdates.DayLocator(interval=-1.5)
    with pytest.raises(ValueError):
        mdates.DayLocator(interval=0)
    with pytest.raises(ValueError):
        mdates.DayLocator(interval=1.3)
    mdates.DayLocator(interval=1.0)


def test_tz_utc():
    dt = datetime.datetime(1970, 1, 1, tzinfo=mdates.UTC)
    dt.tzname()


@pytest.mark.parametrize("x, tdelta",
                         [(1, datetime.timedelta(days=1)),
                          ([1, 1.5], [datetime.timedelta(days=1),
                                      datetime.timedelta(days=1.5)])])
def test_num2timedelta(x, tdelta):
    dt = mdates.num2timedelta(x)
    assert dt == tdelta
