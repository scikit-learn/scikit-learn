""" test the scalar Timedelta """
import pytest

import numpy as np
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm
from pandas.core.tools.timedeltas import _coerce_scalar_to_timedelta_type as ct
from pandas import (Timedelta, TimedeltaIndex, timedelta_range, Series,
                    to_timedelta, compat)
from pandas._libs.tslib import iNaT, NaT


class TestTimedeltaArithmetic(object):

    def test_arithmetic_overflow(self):
        with pytest.raises(OverflowError):
            pd.Timestamp('1700-01-01') + pd.Timedelta(13 * 19999, unit='D')

        with pytest.raises(OverflowError):
            pd.Timestamp('1700-01-01') + timedelta(days=13 * 19999)

    def test_array_timedelta_floordiv(self):
        # https://github.com/pandas-dev/pandas/issues/19761
        ints = pd.date_range('2012-10-08', periods=4, freq='D').view('i8')
        msg = r"Use 'array // timedelta.value'"
        with tm.assert_produces_warning(FutureWarning) as m:
            result = ints // pd.Timedelta(1, unit='s')

        assert msg in str(m[0].message)
        expected = np.array([1349654400, 1349740800, 1349827200, 1349913600],
                            dtype='i8')
        tm.assert_numpy_array_equal(result, expected)

    def test_ops_error_str(self):
        # GH 13624
        td = Timedelta('1 day')

        for left, right in [(td, 'a'), ('a', td)]:

            with pytest.raises(TypeError):
                left + right

            with pytest.raises(TypeError):
                left > right

            assert not left == right
            assert left != right

    def test_ops_notimplemented(self):
        class Other(object):
            pass

        other = Other()

        td = Timedelta('1 day')
        assert td.__add__(other) is NotImplemented
        assert td.__sub__(other) is NotImplemented
        assert td.__truediv__(other) is NotImplemented
        assert td.__mul__(other) is NotImplemented
        assert td.__floordiv__(other) is NotImplemented

    def test_unary_ops(self):
        td = Timedelta(10, unit='d')

        # __neg__, __pos__
        assert -td == Timedelta(-10, unit='d')
        assert -td == Timedelta('-10d')
        assert +td == Timedelta(10, unit='d')

        # __abs__, __abs__(__neg__)
        assert abs(td) == td
        assert abs(-td) == td
        assert abs(-td) == Timedelta('10d')


class TestTimedeltaComparison(object):
    def test_comparison_object_array(self):
        # analogous to GH#15183
        td = Timedelta('2 days')
        other = Timedelta('3 hours')

        arr = np.array([other, td], dtype=object)
        res = arr == td
        expected = np.array([False, True], dtype=bool)
        assert (res == expected).all()

        # 2D case
        arr = np.array([[other, td],
                        [td, other]],
                       dtype=object)
        res = arr != td
        expected = np.array([[True, False], [False, True]], dtype=bool)
        assert res.shape == expected.shape
        assert (res == expected).all()

    def test_compare_timedelta_ndarray(self):
        # GH11835
        periods = [Timedelta('0 days 01:00:00'), Timedelta('0 days 01:00:00')]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)


class TestTimedeltas(object):

    def test_total_seconds_scalar(self):
        # see gh-10939
        rng = Timedelta('1 days, 10:11:12.100123456')
        expt = 1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9
        tm.assert_almost_equal(rng.total_seconds(), expt)

        rng = Timedelta(np.nan)
        assert np.isnan(rng.total_seconds())

    def test_conversion(self):

        for td in [Timedelta(10, unit='d'),
                   Timedelta('1 days, 10:11:12.012345')]:
            pydt = td.to_pytimedelta()
            assert td == Timedelta(pydt)
            assert td == pydt
            assert (isinstance(pydt, timedelta) and not isinstance(
                pydt, Timedelta))

            assert td == np.timedelta64(td.value, 'ns')
            td64 = td.to_timedelta64()

            assert td64 == np.timedelta64(td.value, 'ns')
            assert td == td64

            assert isinstance(td64, np.timedelta64)

        # this is NOT equal and cannot be roundtriped (because of the nanos)
        td = Timedelta('1 days, 10:11:12.012345678')
        assert td != td.to_pytimedelta()

    def test_freq_conversion(self):

        # truediv
        td = Timedelta('1 days 2 hours 3 ns')
        result = td / np.timedelta64(1, 'D')
        assert result == td.value / float(86400 * 1e9)
        result = td / np.timedelta64(1, 's')
        assert result == td.value / float(1e9)
        result = td / np.timedelta64(1, 'ns')
        assert result == td.value

        # floordiv
        td = Timedelta('1 days 2 hours 3 ns')
        result = td // np.timedelta64(1, 'D')
        assert result == 1
        result = td // np.timedelta64(1, 's')
        assert result == 93600
        result = td // np.timedelta64(1, 'ns')
        assert result == td.value

    def test_fields(self):
        def check(value):
            # that we are int/long like
            assert isinstance(value, (int, compat.long))

        # compat to datetime.timedelta
        rng = to_timedelta('1 days, 10:11:12')
        assert rng.days == 1
        assert rng.seconds == 10 * 3600 + 11 * 60 + 12
        assert rng.microseconds == 0
        assert rng.nanoseconds == 0

        pytest.raises(AttributeError, lambda: rng.hours)
        pytest.raises(AttributeError, lambda: rng.minutes)
        pytest.raises(AttributeError, lambda: rng.milliseconds)

        # GH 10050
        check(rng.days)
        check(rng.seconds)
        check(rng.microseconds)
        check(rng.nanoseconds)

        td = Timedelta('-1 days, 10:11:12')
        assert abs(td) == Timedelta('13:48:48')
        assert str(td) == "-1 days +10:11:12"
        assert -td == Timedelta('0 days 13:48:48')
        assert -Timedelta('-1 days, 10:11:12').value == 49728000000000
        assert Timedelta('-1 days, 10:11:12').value == -49728000000000

        rng = to_timedelta('-1 days, 10:11:12.100123456')
        assert rng.days == -1
        assert rng.seconds == 10 * 3600 + 11 * 60 + 12
        assert rng.microseconds == 100 * 1000 + 123
        assert rng.nanoseconds == 456
        pytest.raises(AttributeError, lambda: rng.hours)
        pytest.raises(AttributeError, lambda: rng.minutes)
        pytest.raises(AttributeError, lambda: rng.milliseconds)

        # components
        tup = pd.to_timedelta(-1, 'us').components
        assert tup.days == -1
        assert tup.hours == 23
        assert tup.minutes == 59
        assert tup.seconds == 59
        assert tup.milliseconds == 999
        assert tup.microseconds == 999
        assert tup.nanoseconds == 0

        # GH 10050
        check(tup.days)
        check(tup.hours)
        check(tup.minutes)
        check(tup.seconds)
        check(tup.milliseconds)
        check(tup.microseconds)
        check(tup.nanoseconds)

        tup = Timedelta('-1 days 1 us').components
        assert tup.days == -2
        assert tup.hours == 23
        assert tup.minutes == 59
        assert tup.seconds == 59
        assert tup.milliseconds == 999
        assert tup.microseconds == 999
        assert tup.nanoseconds == 0

    def test_nat_converters(self):
        assert to_timedelta('nat', box=False).astype('int64') == iNaT
        assert to_timedelta('nan', box=False).astype('int64') == iNaT

        def testit(unit, transform):

            # array
            result = to_timedelta(np.arange(5), unit=unit)
            expected = TimedeltaIndex([np.timedelta64(i, transform(unit))
                                       for i in np.arange(5).tolist()])
            tm.assert_index_equal(result, expected)

            # scalar
            result = to_timedelta(2, unit=unit)
            expected = Timedelta(np.timedelta64(2, transform(unit)).astype(
                'timedelta64[ns]'))
            assert result == expected

        # validate all units
        # GH 6855
        for unit in ['Y', 'M', 'W', 'D', 'y', 'w', 'd']:
            testit(unit, lambda x: x.upper())
        for unit in ['days', 'day', 'Day', 'Days']:
            testit(unit, lambda x: 'D')
        for unit in ['h', 'm', 's', 'ms', 'us', 'ns', 'H', 'S', 'MS', 'US',
                     'NS']:
            testit(unit, lambda x: x.lower())

        # offsets

        # m
        testit('T', lambda x: 'm')

        # ms
        testit('L', lambda x: 'ms')

    def test_numeric_conversions(self):
        assert ct(0) == np.timedelta64(0, 'ns')
        assert ct(10) == np.timedelta64(10, 'ns')
        assert ct(10, unit='ns') == np.timedelta64(10, 'ns').astype('m8[ns]')

        assert ct(10, unit='us') == np.timedelta64(10, 'us').astype('m8[ns]')
        assert ct(10, unit='ms') == np.timedelta64(10, 'ms').astype('m8[ns]')
        assert ct(10, unit='s') == np.timedelta64(10, 's').astype('m8[ns]')
        assert ct(10, unit='d') == np.timedelta64(10, 'D').astype('m8[ns]')

    def test_timedelta_conversions(self):
        assert (ct(timedelta(seconds=1)) ==
                np.timedelta64(1, 's').astype('m8[ns]'))
        assert (ct(timedelta(microseconds=1)) ==
                np.timedelta64(1, 'us').astype('m8[ns]'))
        assert (ct(timedelta(days=1)) ==
                np.timedelta64(1, 'D').astype('m8[ns]'))

    def test_round(self):

        t1 = Timedelta('1 days 02:34:56.789123456')
        t2 = Timedelta('-1 days 02:34:56.789123456')

        for (freq, s1, s2) in [('N', t1, t2),
                               ('U', Timedelta('1 days 02:34:56.789123000'),
                                Timedelta('-1 days 02:34:56.789123000')),
                               ('L', Timedelta('1 days 02:34:56.789000000'),
                                Timedelta('-1 days 02:34:56.789000000')),
                               ('S', Timedelta('1 days 02:34:57'),
                                Timedelta('-1 days 02:34:57')),
                               ('2S', Timedelta('1 days 02:34:56'),
                                Timedelta('-1 days 02:34:56')),
                               ('5S', Timedelta('1 days 02:34:55'),
                                Timedelta('-1 days 02:34:55')),
                               ('T', Timedelta('1 days 02:35:00'),
                                Timedelta('-1 days 02:35:00')),
                               ('12T', Timedelta('1 days 02:36:00'),
                                Timedelta('-1 days 02:36:00')),
                               ('H', Timedelta('1 days 03:00:00'),
                                Timedelta('-1 days 03:00:00')),
                               ('d', Timedelta('1 days'),
                                Timedelta('-1 days'))]:
            r1 = t1.round(freq)
            assert r1 == s1
            r2 = t2.round(freq)
            assert r2 == s2

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            pytest.raises(ValueError, lambda: t1.round(freq))

        t1 = timedelta_range('1 days', periods=3, freq='1 min 2 s 3 us')
        t2 = -1 * t1
        t1a = timedelta_range('1 days', periods=3, freq='1 min 2 s')
        t1c = pd.TimedeltaIndex([1, 1, 1], unit='D')

        # note that negative times round DOWN! so don't give whole numbers
        for (freq, s1, s2) in [('N', t1, t2),
                               ('U', t1, t2),
                               ('L', t1a,
                                TimedeltaIndex(['-1 days +00:00:00',
                                                '-2 days +23:58:58',
                                                '-2 days +23:57:56'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('S', t1a,
                                TimedeltaIndex(['-1 days +00:00:00',
                                                '-2 days +23:58:58',
                                                '-2 days +23:57:56'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('12T', t1c,
                                TimedeltaIndex(['-1 days',
                                                '-1 days',
                                                '-1 days'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('H', t1c,
                                TimedeltaIndex(['-1 days',
                                                '-1 days',
                                                '-1 days'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('d', t1c,
                                pd.TimedeltaIndex([-1, -1, -1], unit='D')
                                )]:

            r1 = t1.round(freq)
            tm.assert_index_equal(r1, s1)
            r2 = t2.round(freq)
        tm.assert_index_equal(r2, s2)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            pytest.raises(ValueError, lambda: t1.round(freq))

    def test_contains(self):
        # Checking for any NaT-like objects
        # GH 13603
        td = to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)
        for v in [pd.NaT, None, float('nan'), np.nan]:
            assert not (v in td)

        td = to_timedelta([pd.NaT])
        for v in [pd.NaT, None, float('nan'), np.nan]:
            assert (v in td)

    def test_identity(self):

        td = Timedelta(10, unit='d')
        assert isinstance(td, Timedelta)
        assert isinstance(td, timedelta)

    def test_short_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        assert ct('10') == np.timedelta64(10, 'ns')
        assert ct('10ns') == np.timedelta64(10, 'ns')
        assert ct('100') == np.timedelta64(100, 'ns')
        assert ct('100ns') == np.timedelta64(100, 'ns')

        assert ct('1000') == np.timedelta64(1000, 'ns')
        assert ct('1000ns') == np.timedelta64(1000, 'ns')
        assert ct('1000NS') == np.timedelta64(1000, 'ns')

        assert ct('10us') == np.timedelta64(10000, 'ns')
        assert ct('100us') == np.timedelta64(100000, 'ns')
        assert ct('1000us') == np.timedelta64(1000000, 'ns')
        assert ct('1000Us') == np.timedelta64(1000000, 'ns')
        assert ct('1000uS') == np.timedelta64(1000000, 'ns')

        assert ct('1ms') == np.timedelta64(1000000, 'ns')
        assert ct('10ms') == np.timedelta64(10000000, 'ns')
        assert ct('100ms') == np.timedelta64(100000000, 'ns')
        assert ct('1000ms') == np.timedelta64(1000000000, 'ns')

        assert ct('-1s') == -np.timedelta64(1000000000, 'ns')
        assert ct('1s') == np.timedelta64(1000000000, 'ns')
        assert ct('10s') == np.timedelta64(10000000000, 'ns')
        assert ct('100s') == np.timedelta64(100000000000, 'ns')
        assert ct('1000s') == np.timedelta64(1000000000000, 'ns')

        assert ct('1d') == conv(np.timedelta64(1, 'D'))
        assert ct('-1d') == -conv(np.timedelta64(1, 'D'))
        assert ct('1D') == conv(np.timedelta64(1, 'D'))
        assert ct('10D') == conv(np.timedelta64(10, 'D'))
        assert ct('100D') == conv(np.timedelta64(100, 'D'))
        assert ct('1000D') == conv(np.timedelta64(1000, 'D'))
        assert ct('10000D') == conv(np.timedelta64(10000, 'D'))

        # space
        assert ct(' 10000D ') == conv(np.timedelta64(10000, 'D'))
        assert ct(' - 10000D ') == -conv(np.timedelta64(10000, 'D'))

        # invalid
        pytest.raises(ValueError, ct, '1foo')
        pytest.raises(ValueError, ct, 'foo')

    def test_full_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        d1 = np.timedelta64(1, 'D')

        assert ct('1days') == conv(d1)
        assert ct('1days,') == conv(d1)
        assert ct('- 1days,') == -conv(d1)

        assert ct('00:00:01') == conv(np.timedelta64(1, 's'))
        assert ct('06:00:01') == conv(np.timedelta64(6 * 3600 + 1, 's'))
        assert ct('06:00:01.0') == conv(np.timedelta64(6 * 3600 + 1, 's'))
        assert ct('06:00:01.01') == conv(np.timedelta64(
            1000 * (6 * 3600 + 1) + 10, 'ms'))

        assert (ct('- 1days, 00:00:01') ==
                conv(-d1 + np.timedelta64(1, 's')))
        assert (ct('1days, 06:00:01') ==
                conv(d1 + np.timedelta64(6 * 3600 + 1, 's')))
        assert (ct('1days, 06:00:01.01') ==
                conv(d1 + np.timedelta64(1000 * (6 * 3600 + 1) + 10, 'ms')))

        # invalid
        pytest.raises(ValueError, ct, '- 1days, 00')

    def test_overflow(self):
        # GH 9442
        s = Series(pd.date_range('20130101', periods=100000, freq='H'))
        s[0] += pd.Timedelta('1s 1ms')

        # mean
        result = (s - s.min()).mean()
        expected = pd.Timedelta((pd.DatetimeIndex((s - s.min())).asi8 / len(s)
                                 ).sum())

        # the computation is converted to float so
        # might be some loss of precision
        assert np.allclose(result.value / 1000, expected.value / 1000)

        # sum
        pytest.raises(ValueError, lambda: (s - s.min()).sum())
        s1 = s[0:10000]
        pytest.raises(ValueError, lambda: (s1 - s1.min()).sum())
        s2 = s[0:1000]
        result = (s2 - s2.min()).sum()

    def test_pickle(self):

        v = Timedelta('1 days 10:11:12.0123456')
        v_p = tm.round_trip_pickle(v)
        assert v == v_p

    def test_timedelta_hash_equality(self):
        # GH 11129
        v = Timedelta(1, 'D')
        td = timedelta(days=1)
        assert hash(v) == hash(td)

        d = {td: 2}
        assert d[v] == 2

        tds = timedelta_range('1 second', periods=20)
        assert all(hash(td) == hash(td.to_pytimedelta()) for td in tds)

        # python timedeltas drop ns resolution
        ns_td = Timedelta(1, 'ns')
        assert hash(ns_td) != hash(ns_td.to_pytimedelta())

    def test_implementation_limits(self):
        min_td = Timedelta(Timedelta.min)
        max_td = Timedelta(Timedelta.max)

        # GH 12727
        # timedelta limits correspond to int64 boundaries
        assert min_td.value == np.iinfo(np.int64).min + 1
        assert max_td.value == np.iinfo(np.int64).max

        # Beyond lower limit, a NAT before the Overflow
        assert (min_td - Timedelta(1, 'ns')) is NaT

        with pytest.raises(OverflowError):
            min_td - Timedelta(2, 'ns')

        with pytest.raises(OverflowError):
            max_td + Timedelta(1, 'ns')

        # Same tests using the internal nanosecond values
        td = Timedelta(min_td.value - 1, 'ns')
        assert td is NaT

        with pytest.raises(OverflowError):
            Timedelta(min_td.value - 2, 'ns')

        with pytest.raises(OverflowError):
            Timedelta(max_td.value + 1, 'ns')

    def test_total_seconds_precision(self):
        # GH 19458
        assert Timedelta('30S').total_seconds() == 30.0
        assert Timedelta('0').total_seconds() == 0.0
        assert Timedelta('-2S').total_seconds() == -2.0
        assert Timedelta('5.324S').total_seconds() == 5.324
        assert (Timedelta('30S').total_seconds() - 30.0) < 1e-20
        assert (30.0 - Timedelta('30S').total_seconds()) < 1e-20

    def test_timedelta_arithmetic(self):
        data = pd.Series(['nat', '32 days'], dtype='timedelta64[ns]')
        deltas = [timedelta(days=1), Timedelta(1, unit='D')]
        for delta in deltas:
            result_method = data.add(delta)
            result_operator = data + delta
            expected = pd.Series(['nat', '33 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

            result_method = data.sub(delta)
            result_operator = data - delta
            expected = pd.Series(['nat', '31 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)
            # GH 9396
            result_method = data.div(delta)
            result_operator = data / delta
            expected = pd.Series([np.nan, 32.], dtype='float64')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

    def test_apply_to_timedelta(self):
        timedelta_NaT = pd.to_timedelta('NaT')

        list_of_valid_strings = ['00:00:01', '00:00:02']
        a = pd.to_timedelta(list_of_valid_strings)
        b = Series(list_of_valid_strings).apply(pd.to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

        list_of_strings = ['00:00:01', np.nan, pd.NaT, timedelta_NaT]

        # TODO: unused?
        a = pd.to_timedelta(list_of_strings)  # noqa
        b = Series(list_of_strings).apply(pd.to_timedelta)  # noqa
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

    def test_components(self):
        rng = timedelta_range('1 days, 10:11:12', periods=2, freq='s')
        rng.components

        # with nat
        s = Series(rng)
        s[1] = np.nan

        result = s.dt.components
        assert not result.iloc[0].isna().all()
        assert result.iloc[1].isna().all()
