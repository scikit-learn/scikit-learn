# -*- coding: utf-8 -*-
"""
Tests for scalar Timedelta arithmetic ops
"""
from datetime import datetime, timedelta
import operator

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.core import ops
from pandas import Timedelta, Timestamp, NaT


class TestTimedeltaAdditionSubtraction(object):
    """
    Tests for Timedelta methods:

        __add__, __radd__,
        __sub__, __rsub__
    """
    @pytest.mark.parametrize('ten_seconds', [
        Timedelta(10, unit='s'),
        timedelta(seconds=10),
        np.timedelta64(10, 's'),
        np.timedelta64(10000000000, 'ns'),
        pd.offsets.Second(10)])
    def test_td_add_sub_ten_seconds(self, ten_seconds):
        # GH#6808
        base = Timestamp('20130101 09:01:12.123456')
        expected_add = Timestamp('20130101 09:01:22.123456')
        expected_sub = Timestamp('20130101 09:01:02.123456')

        result = base + ten_seconds
        assert result == expected_add

        result = base - ten_seconds
        assert result == expected_sub

    @pytest.mark.parametrize('one_day_ten_secs', [
        Timedelta('1 day, 00:00:10'),
        Timedelta('1 days, 00:00:10'),
        timedelta(days=1, seconds=10),
        np.timedelta64(1, 'D') + np.timedelta64(10, 's'),
        pd.offsets.Day() + pd.offsets.Second(10)])
    def test_td_add_sub_one_day_ten_seconds(self, one_day_ten_secs):
        # GH#6808
        base = Timestamp('20130102 09:01:12.123456')
        expected_add = Timestamp('20130103 09:01:22.123456')
        expected_sub = Timestamp('20130101 09:01:02.123456')

        result = base + one_day_ten_secs
        assert result == expected_add

        result = base - one_day_ten_secs
        assert result == expected_sub

    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_td_add_datetimelike_scalar(self, op):
        # GH#19738
        td = Timedelta(10, unit='d')

        result = op(td, datetime(2016, 1, 1))
        if op is operator.add:
            # datetime + Timedelta does _not_ call Timedelta.__radd__,
            # so we get a datetime back instead of a Timestamp
            assert isinstance(result, Timestamp)
        assert result == Timestamp(2016, 1, 11)

        result = op(td, Timestamp('2018-01-12 18:09'))
        assert isinstance(result, Timestamp)
        assert result == Timestamp('2018-01-22 18:09')

        result = op(td, np.datetime64('2018-01-12'))
        assert isinstance(result, Timestamp)
        assert result == Timestamp('2018-01-22')

        result = op(td, NaT)
        assert result is NaT

        with pytest.raises(TypeError):
            op(td, 2)
        with pytest.raises(TypeError):
            op(td, 2.0)

    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_td_add_td(self, op):
        td = Timedelta(10, unit='d')

        result = op(td, Timedelta(days=10))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=20)

    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_td_add_pytimedelta(self, op):
        td = Timedelta(10, unit='d')
        result = op(td, timedelta(days=9))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=19)

    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_td_add_timedelta64(self, op):
        td = Timedelta(10, unit='d')
        result = op(td, np.timedelta64(-4, 'D'))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=6)

    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_td_add_offset(self, op):
        td = Timedelta(10, unit='d')

        result = op(td, pd.offsets.Hour(6))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=10, hours=6)

    def test_td_sub_td(self):
        td = Timedelta(10, unit='d')
        expected = Timedelta(0, unit='ns')
        result = td - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_pytimedelta(self):
        td = Timedelta(10, unit='d')
        expected = Timedelta(0, unit='ns')
        result = td - td.to_pytimedelta()
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_timedelta64(self):
        td = Timedelta(10, unit='d')
        expected = Timedelta(0, unit='ns')
        result = td - td.to_timedelta64()
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_nat(self):
        td = Timedelta(10, unit='d')
        result = td - NaT
        assert result is NaT

    def test_td_sub_td64_nat(self):
        td = Timedelta(10, unit='d')
        result = td - np.timedelta64('NaT')
        assert result is NaT

    def test_td_sub_offset(self):
        td = Timedelta(10, unit='d')
        result = td - pd.offsets.Hour(1)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(239, unit='h')

    def test_td_sub_numeric_raises(self):
        td = td = Timedelta(10, unit='d')
        with pytest.raises(TypeError):
            td - 2
        with pytest.raises(TypeError):
            td - 2.0

    def test_td_rsub_pytimedelta(self):
        td = Timedelta(10, unit='d')
        expected = Timedelta(0, unit='ns')

        result = td.to_pytimedelta() - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_rsub_timedelta64(self):
        td = Timedelta(10, unit='d')
        expected = Timedelta(0, unit='ns')

        result = td.to_timedelta64() - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_rsub_nat(self):
        td = Timedelta(10, unit='d')
        result = NaT - td
        assert result is NaT

        result = np.datetime64('NaT') - td
        assert result is NaT

    def test_td_rsub_td64_nat(self):
        td = Timedelta(10, unit='d')
        result = np.timedelta64('NaT') - td
        assert result is NaT

    def test_td_rsub_offset(self):
        result = pd.offsets.Hour(1) - Timedelta(10, unit='d')
        assert isinstance(result, Timedelta)
        assert result == Timedelta(-239, unit='h')

    def test_td_rsub_numeric_raises(self):
        td = td = Timedelta(10, unit='d')
        with pytest.raises(TypeError):
            2 - td
        with pytest.raises(TypeError):
            2.0 - td


class TestTimedeltaMultiplicationDivision(object):
    """
    Tests for Timedelta methods:

        __mul__, __rmul__,
        __div__, __rdiv__,
        __truediv__, __rtruediv__,
        __floordiv__, __rfloordiv__,
        __mod__, __rmod__,
        __divmod__, __rdivmod__
    """

    # ---------------------------------------------------------------
    # Timedelta.__mul__, __rmul__

    @pytest.mark.parametrize('td_nat', [pd.NaT,
                                        np.timedelta64('NaT', 'ns'),
                                        np.timedelta64('NaT')])
    @pytest.mark.parametrize('op', [operator.mul, ops.rmul])
    def test_td_mul_nat(self, op, td_nat):
        # GH#19819
        td = Timedelta(10, unit='d')
        with pytest.raises(TypeError):
            op(td, td_nat)

    @pytest.mark.parametrize('op', [operator.mul, ops.rmul])
    def test_td_mul_scalar(self, op):
        # GH#19738
        td = Timedelta(minutes=3)

        result = op(td, 2)
        assert result == Timedelta(minutes=6)

        result = op(td, 1.5)
        assert result == Timedelta(minutes=4, seconds=30)

        assert op(td, np.nan) is NaT

        assert op(-1, td).value == -1 * td.value
        assert op(-1.0, td).value == -1.0 * td.value

        with pytest.raises(TypeError):
            # timedelta * datetime is gibberish
            op(td, Timestamp(2016, 1, 2))

        with pytest.raises(TypeError):
            # invalid multiply with another timedelta
            op(td, td)

    # ---------------------------------------------------------------
    # Timedelta.__div__, __truediv__

    def test_td_div_timedeltalike_scalar(self):
        # GH#19738
        td = Timedelta(10, unit='d')

        result = td / pd.offsets.Hour(1)
        assert result == 240

        assert td / td == 1
        assert td / np.timedelta64(60, 'h') == 4

        assert np.isnan(td / NaT)

    def test_td_div_numeric_scalar(self):
        # GH#19738
        td = Timedelta(10, unit='d')

        result = td / 2
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=5)

        result = td / 5.0
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=2)

    # ---------------------------------------------------------------
    # Timedelta.__rdiv__

    def test_td_rdiv_timedeltalike_scalar(self):
        # GH#19738
        td = Timedelta(10, unit='d')
        result = pd.offsets.Hour(1) / td
        assert result == 1 / 240.0

        assert np.timedelta64(60, 'h') / td == 0.25

    # ---------------------------------------------------------------
    # Timedelta.__floordiv__

    def test_td_floordiv_timedeltalike_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        scalar = Timedelta(hours=3, minutes=3)

        assert td // scalar == 1
        assert -td // scalar.to_pytimedelta() == -2
        assert (2 * td) // scalar.to_timedelta64() == 2

    def test_td_floordiv_null_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        assert td // np.nan is NaT
        assert np.isnan(td // NaT)
        assert np.isnan(td // np.timedelta64('NaT'))

    def test_td_floordiv_offsets(self):
        # GH#19738
        td = Timedelta(hours=3, minutes=4)
        assert td // pd.offsets.Hour(1) == 3
        assert td // pd.offsets.Minute(2) == 92

    def test_td_floordiv_invalid_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        with pytest.raises(TypeError):
            td // np.datetime64('2016-01-01', dtype='datetime64[us]')

    def test_td_floordiv_numeric_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        expected = Timedelta(hours=1, minutes=32)
        assert td // 2 == expected
        assert td // 2.0 == expected
        assert td // np.float64(2.0) == expected
        assert td // np.int32(2.0) == expected
        assert td // np.uint8(2.0) == expected

    def test_td_floordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        scalar = Timedelta(hours=3, minutes=3)

        # Array-like others
        assert td // np.array(scalar.to_timedelta64()) == 1

        res = (3 * td) // np.array([scalar.to_timedelta64()])
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        res = (10 * td) // np.array([scalar.to_timedelta64(),
                                     np.timedelta64('NaT')])
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_td_floordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        ser = pd.Series([1], dtype=np.int64)
        res = td // ser
        assert res.dtype.kind == 'm'

    # ---------------------------------------------------------------
    # Timedelta.__rfloordiv__

    def test_td_rfloordiv_timedeltalike_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        scalar = Timedelta(hours=3, minutes=4)

        # scalar others
        # x // Timedelta is defined only for timedelta-like x. int-like,
        # float-like, and date-like, in particular, should all either
        # a) raise TypeError directly or
        # b) return NotImplemented, following which the reversed
        #    operation will raise TypeError.
        assert td.__rfloordiv__(scalar) == 1
        assert (-td).__rfloordiv__(scalar.to_pytimedelta()) == -2
        assert (2 * td).__rfloordiv__(scalar.to_timedelta64()) == 0

    def test_td_rfloordiv_null_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        assert np.isnan(td.__rfloordiv__(NaT))
        assert np.isnan(td.__rfloordiv__(np.timedelta64('NaT')))

    def test_td_rfloordiv_offsets(self):
        # GH#19738
        assert pd.offsets.Hour(1) // Timedelta(minutes=25) == 2

    def test_td_rfloordiv_invalid_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        dt64 = np.datetime64('2016-01-01', dtype='datetime64[us]')
        with pytest.raises(TypeError):
            td.__rfloordiv__(dt64)

    def test_td_rfloordiv_numeric_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        assert td.__rfloordiv__(np.nan) is NotImplemented
        assert td.__rfloordiv__(3.5) is NotImplemented
        assert td.__rfloordiv__(2) is NotImplemented

        with pytest.raises(TypeError):
            td.__rfloordiv__(np.float64(2.0))
        with pytest.raises(TypeError):
            td.__rfloordiv__(np.uint8(9))
        with tm.assert_produces_warning(FutureWarning):
            # GH-19761: Change to TypeError.
            td.__rfloordiv__(np.int32(2.0))

    def test_td_rfloordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        scalar = Timedelta(hours=3, minutes=4)

        # Array-like others
        assert td.__rfloordiv__(np.array(scalar.to_timedelta64())) == 1

        res = td.__rfloordiv__(np.array([(3 * scalar).to_timedelta64()]))
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        arr = np.array([(10 * scalar).to_timedelta64(),
                        np.timedelta64('NaT')])
        res = td.__rfloordiv__(arr)
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_td_rfloordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        ser = pd.Series([1], dtype=np.int64)
        res = td.__rfloordiv__(ser)
        assert res is NotImplemented
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # TODO: GH-19761. Change to TypeError.
            ser // td

    def test_mod_timedeltalike(self):
        # GH#19365
        td = Timedelta(hours=37)

        # Timedelta-like others
        result = td % Timedelta(hours=6)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=1)

        result = td % timedelta(minutes=60)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(0)

        result = td % NaT
        assert result is NaT

    def test_mod_timedelta64_nat(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % np.timedelta64('NaT', 'ns')
        assert result is NaT

    def test_mod_timedelta64(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % np.timedelta64(2, 'h')
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=1)

    def test_mod_offset(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % pd.offsets.Hour(5)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=2)

    # ----------------------------------------------------------------
    # Timedelta.__mod__, __rmod__

    def test_mod_numeric(self):
        # GH#19365
        td = Timedelta(hours=37)

        # Numeric Others
        result = td % 2
        assert isinstance(result, Timedelta)
        assert result == Timedelta(0)

        result = td % 1e12
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

        result = td % int(1e12)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

    def test_mod_invalid(self):
        # GH#19365
        td = Timedelta(hours=37)

        with pytest.raises(TypeError):
            td % pd.Timestamp('2018-01-22')

        with pytest.raises(TypeError):
            td % []

    def test_rmod_pytimedelta(self):
        # GH#19365
        td = Timedelta(minutes=3)

        result = timedelta(minutes=4) % td
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=1)

    def test_rmod_timedelta64(self):
        # GH#19365
        td = Timedelta(minutes=3)
        result = np.timedelta64(5, 'm') % td
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=2)

    def test_rmod_invalid(self):
        # GH#19365
        td = Timedelta(minutes=3)

        with pytest.raises(TypeError):
            pd.Timestamp('2018-01-22') % td

        with pytest.raises(TypeError):
            15 % td

        with pytest.raises(TypeError):
            16.0 % td

        with pytest.raises(TypeError):
            np.array([22, 24]) % td

    # ----------------------------------------------------------------
    # Timedelta.__divmod__, __rdivmod__

    def test_divmod_numeric(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, 53 * 3600 * 1e9)
        assert result[0] == Timedelta(1, unit='ns')
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=1)

        assert result
        result = divmod(td, np.nan)
        assert result[0] is pd.NaT
        assert result[1] is pd.NaT

    def test_divmod(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, timedelta(days=1))
        assert result[0] == 2
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=6)

        result = divmod(td, 54)
        assert result[0] == Timedelta(hours=1)
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(0)

        result = divmod(td, pd.NaT)
        assert np.isnan(result[0])
        assert result[1] is pd.NaT

    def test_divmod_offset(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, pd.offsets.Hour(-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_divmod_invalid(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        with pytest.raises(TypeError):
            divmod(td, pd.Timestamp('2018-01-22'))

    def test_rdivmod_pytimedelta(self):
        # GH#19365
        result = divmod(timedelta(days=2, hours=6), Timedelta(days=1))
        assert result[0] == 2
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=6)

    def test_rdivmod_offset(self):
        result = divmod(pd.offsets.Hour(54), Timedelta(hours=-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_rdivmod_invalid(self):
        # GH#19365
        td = Timedelta(minutes=3)

        with pytest.raises(TypeError):
            divmod(pd.Timestamp('2018-01-22'), td)

        with pytest.raises(TypeError):
            divmod(15, td)

        with pytest.raises(TypeError):
            divmod(16.0, td)

        with pytest.raises(TypeError):
            divmod(np.array([22, 24]), td)
