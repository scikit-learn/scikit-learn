from __future__ import division

import numpy as np
from pandas import Interval, Timestamp, Timedelta
import pandas.core.common as com

import pytest
import pandas.util.testing as tm


@pytest.fixture
def interval():
    return Interval(0, 1)


class TestInterval(object):

    def test_properties(self, interval):
        assert interval.closed == 'right'
        assert interval.left == 0
        assert interval.right == 1
        assert interval.mid == 0.5

    def test_repr(self, interval):
        assert repr(interval) == "Interval(0, 1, closed='right')"
        assert str(interval) == "(0, 1]"

        interval_left = Interval(0, 1, closed='left')
        assert repr(interval_left) == "Interval(0, 1, closed='left')"
        assert str(interval_left) == "[0, 1)"

    def test_contains(self, interval):
        assert 0.5 in interval
        assert 1 in interval
        assert 0 not in interval

        msg = "__contains__ not defined for two intervals"
        with tm.assert_raises_regex(TypeError, msg):
            interval in interval

        interval_both = Interval(0, 1, closed='both')
        assert 0 in interval_both
        assert 1 in interval_both

        interval_neither = Interval(0, 1, closed='neither')
        assert 0 not in interval_neither
        assert 0.5 in interval_neither
        assert 1 not in interval_neither

    def test_equal(self):
        assert Interval(0, 1) == Interval(0, 1, closed='right')
        assert Interval(0, 1) != Interval(0, 1, closed='left')
        assert Interval(0, 1) != 0

    def test_comparison(self):
        with tm.assert_raises_regex(TypeError, 'unorderable types'):
            Interval(0, 1) < 2

        assert Interval(0, 1) < Interval(1, 2)
        assert Interval(0, 1) < Interval(0, 2)
        assert Interval(0, 1) < Interval(0.5, 1.5)
        assert Interval(0, 1) <= Interval(0, 1)
        assert Interval(0, 1) > Interval(-1, 2)
        assert Interval(0, 1) >= Interval(0, 1)

    def test_hash(self, interval):
        # should not raise
        hash(interval)

    @pytest.mark.parametrize('left, right, expected', [
        (0, 5, 5),
        (-2, 5.5, 7.5),
        (10, 10, 0),
        (10, np.inf, np.inf),
        (-np.inf, -5, np.inf),
        (-np.inf, np.inf, np.inf),
        (Timedelta('0 days'), Timedelta('5 days'), Timedelta('5 days')),
        (Timedelta('10 days'), Timedelta('10 days'), Timedelta('0 days')),
        (Timedelta('1H10M'), Timedelta('5H5M'), Timedelta('3H55M')),
        (Timedelta('5S'), Timedelta('1H'), Timedelta('59M55S'))])
    def test_length(self, left, right, expected):
        # GH 18789
        iv = Interval(left, right)
        result = iv.length
        assert result == expected

    @pytest.mark.parametrize('left, right, expected', [
        ('2017-01-01', '2017-01-06', '5 days'),
        ('2017-01-01', '2017-01-01 12:00:00', '12 hours'),
        ('2017-01-01 12:00', '2017-01-01 12:00:00', '0 days'),
        ('2017-01-01 12:01', '2017-01-05 17:31:00', '4 days 5 hours 30 min')])
    @pytest.mark.parametrize('tz', (None, 'UTC', 'CET', 'US/Eastern'))
    def test_length_timestamp(self, tz, left, right, expected):
        # GH 18789
        iv = Interval(Timestamp(left, tz=tz), Timestamp(right, tz=tz))
        result = iv.length
        expected = Timedelta(expected)
        assert result == expected

    @pytest.mark.parametrize('left, right', [
        ('a', 'z'),
        (('a', 'b'), ('c', 'd')),
        (list('AB'), list('ab')),
        (Interval(0, 1), Interval(1, 2))])
    def test_length_errors(self, left, right):
        # GH 18789
        iv = Interval(left, right)
        msg = 'cannot compute length between .* and .*'
        with tm.assert_raises_regex(TypeError, msg):
            iv.length

    def test_math_add(self, interval):
        expected = Interval(1, 2)
        actual = interval + 1
        assert expected == actual

        expected = Interval(1, 2)
        actual = 1 + interval
        assert expected == actual

        actual = interval
        actual += 1
        assert expected == actual

        msg = r"unsupported operand type\(s\) for \+"
        with tm.assert_raises_regex(TypeError, msg):
            interval + Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval + 'foo'

    def test_math_sub(self, interval):
        expected = Interval(-1, 0)
        actual = interval - 1
        assert expected == actual

        actual = interval
        actual -= 1
        assert expected == actual

        msg = r"unsupported operand type\(s\) for -"
        with tm.assert_raises_regex(TypeError, msg):
            interval - Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval - 'foo'

    def test_math_mult(self, interval):
        expected = Interval(0, 2)
        actual = interval * 2
        assert expected == actual

        expected = Interval(0, 2)
        actual = 2 * interval
        assert expected == actual

        actual = interval
        actual *= 2
        assert expected == actual

        msg = r"unsupported operand type\(s\) for \*"
        with tm.assert_raises_regex(TypeError, msg):
            interval * Interval(1, 2)

        msg = r"can\'t multiply sequence by non-int"
        with tm.assert_raises_regex(TypeError, msg):
            interval * 'foo'

    def test_math_div(self, interval):
        expected = Interval(0, 0.5)
        actual = interval / 2.0
        assert expected == actual

        actual = interval
        actual /= 2.0
        assert expected == actual

        msg = r"unsupported operand type\(s\) for /"
        with tm.assert_raises_regex(TypeError, msg):
            interval / Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval / 'foo'

    def test_constructor_errors(self):
        msg = "invalid option for 'closed': foo"
        with tm.assert_raises_regex(ValueError, msg):
            Interval(0, 1, closed='foo')

        msg = 'left side of interval must be <= right side'
        with tm.assert_raises_regex(ValueError, msg):
            Interval(1, 0)

    @pytest.mark.parametrize('tz_left, tz_right', [
        (None, 'UTC'), ('UTC', None), ('UTC', 'US/Eastern')])
    def test_constructor_errors_tz(self, tz_left, tz_right):
        # GH 18538
        left = Timestamp('2017-01-01', tz=tz_left)
        right = Timestamp('2017-01-02', tz=tz_right)
        error = TypeError if com._any_none(tz_left, tz_right) else ValueError
        with pytest.raises(error):
            Interval(left, right)
