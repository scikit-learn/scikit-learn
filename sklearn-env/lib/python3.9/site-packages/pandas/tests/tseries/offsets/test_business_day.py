"""
Tests for offsets.BDay
"""
from datetime import (
    date,
    datetime,
    timedelta,
)

import numpy as np
import pytest

from pandas._libs.tslibs.offsets import (
    ApplyTypeError,
    BDay,
    BMonthEnd,
)

from pandas import (
    DatetimeIndex,
    Timedelta,
    _testing as tm,
)
from pandas.tests.tseries.offsets.common import (
    Base,
    assert_is_on_offset,
    assert_offset_equal,
)
from pandas.tests.tseries.offsets.test_offsets import _ApplyCases

from pandas.tseries import offsets as offsets


class TestBusinessDay(Base):
    _offset = BDay

    def setup_method(self, method):
        self.d = datetime(2008, 1, 1)
        self.nd = np.datetime64("2008-01-01 00:00:00")

        self.offset = self._offset()
        self.offset1 = self.offset
        self.offset2 = self._offset(2)

    def test_different_normalize_equals(self):
        # GH#21404 changed __eq__ to return False when `normalize` does not match
        offset = self._offset()
        offset2 = self._offset(normalize=True)
        assert offset != offset2

    def test_repr(self):
        assert repr(self.offset) == "<BusinessDay>"
        assert repr(self.offset2) == "<2 * BusinessDays>"

        expected = "<BusinessDay: offset=datetime.timedelta(days=1)>"
        assert repr(self.offset + timedelta(1)) == expected

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    @pytest.mark.parametrize(
        "td",
        [
            Timedelta(hours=2),
            Timedelta(hours=2).to_pytimedelta(),
            Timedelta(hours=2).to_timedelta64(),
        ],
        ids=lambda x: type(x),
    )
    def test_with_offset_index(self, td):

        dti = DatetimeIndex([self.d])
        expected = DatetimeIndex([datetime(2008, 1, 2, 2)])

        result = dti + (td + self.offset)
        tm.assert_index_equal(result, expected)

        result = dti + (self.offset + td)
        tm.assert_index_equal(result, expected)

    def test_eq(self):
        assert self.offset2 == self.offset2

    def test_mul(self):
        pass

    def test_hash(self):
        assert hash(self.offset2) == hash(self.offset2)

    def test_call(self):
        with tm.assert_produces_warning(FutureWarning):
            # GH#34171 DateOffset.__call__ is deprecated
            assert self.offset2(self.d) == datetime(2008, 1, 3)
            assert self.offset2(self.nd) == datetime(2008, 1, 3)

    def testRollback1(self):
        assert self._offset(10).rollback(self.d) == self.d

    def testRollback2(self):
        assert self._offset(10).rollback(datetime(2008, 1, 5)) == datetime(2008, 1, 4)

    def testRollforward1(self):
        assert self._offset(10).rollforward(self.d) == self.d

    def testRollforward2(self):
        assert self._offset(10).rollforward(datetime(2008, 1, 5)) == datetime(
            2008, 1, 7
        )

    def test_roll_date_object(self):
        offset = self._offset()

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

    def test_is_on_offset(self):
        tests = [
            (self._offset(), datetime(2008, 1, 1), True),
            (self._offset(), datetime(2008, 1, 5), False),
        ]

        for offset, d, expected in tests:
            assert_is_on_offset(offset, d, expected)

    apply_cases: _ApplyCases = [
        (
            1,
            {
                datetime(2008, 1, 1): datetime(2008, 1, 2),
                datetime(2008, 1, 4): datetime(2008, 1, 7),
                datetime(2008, 1, 5): datetime(2008, 1, 7),
                datetime(2008, 1, 6): datetime(2008, 1, 7),
                datetime(2008, 1, 7): datetime(2008, 1, 8),
            },
        ),
        (
            2,
            {
                datetime(2008, 1, 1): datetime(2008, 1, 3),
                datetime(2008, 1, 4): datetime(2008, 1, 8),
                datetime(2008, 1, 5): datetime(2008, 1, 8),
                datetime(2008, 1, 6): datetime(2008, 1, 8),
                datetime(2008, 1, 7): datetime(2008, 1, 9),
            },
        ),
        (
            -1,
            {
                datetime(2008, 1, 1): datetime(2007, 12, 31),
                datetime(2008, 1, 4): datetime(2008, 1, 3),
                datetime(2008, 1, 5): datetime(2008, 1, 4),
                datetime(2008, 1, 6): datetime(2008, 1, 4),
                datetime(2008, 1, 7): datetime(2008, 1, 4),
                datetime(2008, 1, 8): datetime(2008, 1, 7),
            },
        ),
        (
            -2,
            {
                datetime(2008, 1, 1): datetime(2007, 12, 28),
                datetime(2008, 1, 4): datetime(2008, 1, 2),
                datetime(2008, 1, 5): datetime(2008, 1, 3),
                datetime(2008, 1, 6): datetime(2008, 1, 3),
                datetime(2008, 1, 7): datetime(2008, 1, 3),
                datetime(2008, 1, 8): datetime(2008, 1, 4),
                datetime(2008, 1, 9): datetime(2008, 1, 7),
            },
        ),
        (
            0,
            {
                datetime(2008, 1, 1): datetime(2008, 1, 1),
                datetime(2008, 1, 4): datetime(2008, 1, 4),
                datetime(2008, 1, 5): datetime(2008, 1, 7),
                datetime(2008, 1, 6): datetime(2008, 1, 7),
                datetime(2008, 1, 7): datetime(2008, 1, 7),
            },
        ),
    ]

    @pytest.mark.parametrize("case", apply_cases)
    def test_apply(self, case):
        n, cases = case
        offset = self._offset(n)
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + self._offset(10)
        assert result == datetime(2012, 11, 6)

        result = dt + self._offset(100) - self._offset(100)
        assert result == dt

        off = self._offset() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 12, 23)
        assert rs == xp

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2011, 12, 26)
        assert rs == xp

        off = self._offset() * 10
        rs = datetime(2014, 1, 5) + off  # see #5890
        xp = datetime(2014, 1, 17)
        assert rs == xp

    def test_apply_corner(self):
        if self._offset is BDay:
            msg = "Only know how to combine business day with datetime or timedelta"
        else:
            msg = (
                "Only know how to combine trading day "
                "with datetime, datetime64 or timedelta"
            )
        with pytest.raises(ApplyTypeError, match=msg):
            self._offset()._apply(BMonthEnd())
