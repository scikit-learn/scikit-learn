# -*- coding: utf-8 -*-
from datetime import timedelta

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import Timedelta


def test_construction():
    expected = np.timedelta64(10, 'D').astype('m8[ns]').view('i8')
    assert Timedelta(10, unit='d').value == expected
    assert Timedelta(10.0, unit='d').value == expected
    assert Timedelta('10 days').value == expected
    assert Timedelta(days=10).value == expected
    assert Timedelta(days=10.0).value == expected

    expected += np.timedelta64(10, 's').astype('m8[ns]').view('i8')
    assert Timedelta('10 days 00:00:10').value == expected
    assert Timedelta(days=10, seconds=10).value == expected
    assert Timedelta(days=10, milliseconds=10 * 1000).value == expected
    assert Timedelta(days=10,
                     microseconds=10 * 1000 * 1000).value == expected

    # rounding cases
    assert Timedelta(82739999850000).value == 82739999850000
    assert ('0 days 22:58:59.999850' in str(Timedelta(82739999850000)))
    assert Timedelta(123072001000000).value == 123072001000000
    assert ('1 days 10:11:12.001' in str(Timedelta(123072001000000)))

    # string conversion with/without leading zero
    # GH#9570
    assert Timedelta('0:00:00') == timedelta(hours=0)
    assert Timedelta('00:00:00') == timedelta(hours=0)
    assert Timedelta('-1:00:00') == -timedelta(hours=1)
    assert Timedelta('-01:00:00') == -timedelta(hours=1)

    # more strings & abbrevs
    # GH#8190
    assert Timedelta('1 h') == timedelta(hours=1)
    assert Timedelta('1 hour') == timedelta(hours=1)
    assert Timedelta('1 hr') == timedelta(hours=1)
    assert Timedelta('1 hours') == timedelta(hours=1)
    assert Timedelta('-1 hours') == -timedelta(hours=1)
    assert Timedelta('1 m') == timedelta(minutes=1)
    assert Timedelta('1.5 m') == timedelta(seconds=90)
    assert Timedelta('1 minute') == timedelta(minutes=1)
    assert Timedelta('1 minutes') == timedelta(minutes=1)
    assert Timedelta('1 s') == timedelta(seconds=1)
    assert Timedelta('1 second') == timedelta(seconds=1)
    assert Timedelta('1 seconds') == timedelta(seconds=1)
    assert Timedelta('1 ms') == timedelta(milliseconds=1)
    assert Timedelta('1 milli') == timedelta(milliseconds=1)
    assert Timedelta('1 millisecond') == timedelta(milliseconds=1)
    assert Timedelta('1 us') == timedelta(microseconds=1)
    assert Timedelta('1 micros') == timedelta(microseconds=1)
    assert Timedelta('1 microsecond') == timedelta(microseconds=1)
    assert Timedelta('1.5 microsecond') == Timedelta('00:00:00.000001500')
    assert Timedelta('1 ns') == Timedelta('00:00:00.000000001')
    assert Timedelta('1 nano') == Timedelta('00:00:00.000000001')
    assert Timedelta('1 nanosecond') == Timedelta('00:00:00.000000001')

    # combos
    assert Timedelta('10 days 1 hour') == timedelta(days=10, hours=1)
    assert Timedelta('10 days 1 h') == timedelta(days=10, hours=1)
    assert Timedelta('10 days 1 h 1m 1s') == timedelta(
        days=10, hours=1, minutes=1, seconds=1)
    assert Timedelta('-10 days 1 h 1m 1s') == -timedelta(
        days=10, hours=1, minutes=1, seconds=1)
    assert Timedelta('-10 days 1 h 1m 1s') == -timedelta(
        days=10, hours=1, minutes=1, seconds=1)
    assert Timedelta('-10 days 1 h 1m 1s 3us') == -timedelta(
        days=10, hours=1, minutes=1, seconds=1, microseconds=3)
    assert Timedelta('-10 days 1 h 1.5m 1s 3us') == -timedelta(
        days=10, hours=1, minutes=1, seconds=31, microseconds=3)

    # Currently invalid as it has a - on the hh:mm:dd part
    # (only allowed on the days)
    with pytest.raises(ValueError):
        Timedelta('-10 days -1 h 1.5m 1s 3us')

    # only leading neg signs are allowed
    with pytest.raises(ValueError):
        Timedelta('10 days -1 h 1.5m 1s 3us')

    # no units specified
    with pytest.raises(ValueError):
        Timedelta('3.1415')

    # invalid construction
    tm.assert_raises_regex(ValueError, "cannot construct a Timedelta",
                           lambda: Timedelta())
    tm.assert_raises_regex(ValueError,
                           "unit abbreviation w/o a number",
                           lambda: Timedelta('foo'))
    tm.assert_raises_regex(ValueError,
                           "cannot construct a Timedelta from the "
                           "passed arguments, allowed keywords are ",
                           lambda: Timedelta(day=10))

    # floats
    expected = np.timedelta64(
        10, 's').astype('m8[ns]').view('i8') + np.timedelta64(
            500, 'ms').astype('m8[ns]').view('i8')
    assert Timedelta(10.5, unit='s').value == expected

    # offset
    assert pd.to_timedelta(pd.offsets.Hour(2)) == Timedelta(hours=2)
    assert Timedelta(pd.offsets.Hour(2)) == Timedelta(hours=2)
    assert Timedelta(pd.offsets.Second(2)) == Timedelta(seconds=2)

    # GH#11995: unicode
    expected = Timedelta('1H')
    result = pd.Timedelta(u'1H')
    assert result == expected
    assert (pd.to_timedelta(pd.offsets.Hour(2)) ==
            Timedelta(u'0 days, 02:00:00'))

    with pytest.raises(ValueError):
        Timedelta(u'foo bar')


@pytest.mark.parametrize('item', list({'days': 'D',
                                       'seconds': 's',
                                       'microseconds': 'us',
                                       'milliseconds': 'ms',
                                       'minutes': 'm',
                                       'hours': 'h',
                                       'weeks': 'W'}.items()))
@pytest.mark.parametrize('npdtype', [np.int64, np.int32, np.int16,
                                     np.float64, np.float32, np.float16])
def test_td_construction_with_np_dtypes(npdtype, item):
    # GH#8757: test construction with np dtypes
    pykwarg, npkwarg = item
    expected = np.timedelta64(1, npkwarg).astype('m8[ns]').view('i8')
    assert Timedelta(**{pykwarg: npdtype(1)}).value == expected


@pytest.mark.parametrize('val', [
    '1s', '-1s', '1us', '-1us', '1 day', '-1 day',
    '-23:59:59.999999', '-1 days +23:59:59.999999', '-1ns',
    '1ns', '-23:59:59.999999999'])
def test_td_from_repr_roundtrip(val):
    # round-trip both for string and value
    td = Timedelta(val)
    assert Timedelta(td.value) == td

    # str does not normally display nanos
    if not td.nanoseconds:
        assert Timedelta(str(td)) == td
    assert Timedelta(td._repr_base(format='all')) == td


def test_overflow_on_construction():
    # xref https://github.com/statsmodels/statsmodels/issues/3374
    value = pd.Timedelta('1day').value * 20169940
    with pytest.raises(OverflowError):
        pd.Timedelta(value)

    # xref GH#17637
    with pytest.raises(OverflowError):
        pd.Timedelta(7 * 19999, unit='D')

    with pytest.raises(OverflowError):
        pd.Timedelta(timedelta(days=13 * 19999))


@pytest.mark.parametrize('fmt,exp', [
    ('P6DT0H50M3.010010012S', Timedelta(days=6, minutes=50, seconds=3,
                                        milliseconds=10, microseconds=10,
                                        nanoseconds=12)),
    ('P-6DT0H50M3.010010012S', Timedelta(days=-6, minutes=50, seconds=3,
                                         milliseconds=10, microseconds=10,
                                         nanoseconds=12)),
    ('P4DT12H30M5S', Timedelta(days=4, hours=12, minutes=30, seconds=5)),
    ('P0DT0H0M0.000000123S', Timedelta(nanoseconds=123)),
    ('P0DT0H0M0.00001S', Timedelta(microseconds=10)),
    ('P0DT0H0M0.001S', Timedelta(milliseconds=1)),
    ('P0DT0H1M0S', Timedelta(minutes=1)),
    ('P1DT25H61M61S', Timedelta(days=1, hours=25, minutes=61, seconds=61))
])
def test_iso_constructor(fmt, exp):
    assert Timedelta(fmt) == exp


@pytest.mark.parametrize('fmt', [
    'PPPPPPPPPPPP', 'PDTHMS', 'P0DT999H999M999S',
    'P1DT0H0M0.0000000000000S', 'P1DT0H0M00000000000S',
    'P1DT0H0M0.S'])
def test_iso_constructor_raises(fmt):
    with tm.assert_raises_regex(ValueError, 'Invalid ISO 8601 Duration '
                                'format - {}'.format(fmt)):
        Timedelta(fmt)


@pytest.mark.parametrize('constructed_td, conversion', [
    (Timedelta(nanoseconds=100), '100ns'),
    (Timedelta(days=1, hours=1, minutes=1, weeks=1, seconds=1, milliseconds=1,
               microseconds=1, nanoseconds=1), 694861001001001),
    (Timedelta(microseconds=1) + Timedelta(nanoseconds=1), '1us1ns'),
    (Timedelta(microseconds=1) - Timedelta(nanoseconds=1), '999ns'),
    (Timedelta(microseconds=1) + 5 * Timedelta(nanoseconds=-2), '990ns')])
def test_td_constructor_on_nanoseconds(constructed_td, conversion):
    # GH#9273
    assert constructed_td == Timedelta(conversion)


def test_td_constructor_value_error():
    with pytest.raises(TypeError):
        Timedelta(nanoseconds='abc')
