"""
Tests for the following offsets:
- YearBegin
- YearEnd
"""
from datetime import datetime

import pytest

from pandas.tests.tseries.offsets.common import (
    Base,
    assert_is_on_offset,
    assert_offset_equal,
)

from pandas.tseries.offsets import (
    YearBegin,
    YearEnd,
)


class TestYearBegin(Base):
    _offset = YearBegin

    def test_misspecified(self):
        with pytest.raises(ValueError, match="Month must go from 1 to 12"):
            YearBegin(month=13)

    offset_cases = []
    offset_cases.append(
        (
            YearBegin(),
            {
                datetime(2008, 1, 1): datetime(2009, 1, 1),
                datetime(2008, 6, 30): datetime(2009, 1, 1),
                datetime(2008, 12, 31): datetime(2009, 1, 1),
                datetime(2005, 12, 30): datetime(2006, 1, 1),
                datetime(2005, 12, 31): datetime(2006, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(0),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 1),
                datetime(2008, 6, 30): datetime(2009, 1, 1),
                datetime(2008, 12, 31): datetime(2009, 1, 1),
                datetime(2005, 12, 30): datetime(2006, 1, 1),
                datetime(2005, 12, 31): datetime(2006, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(3),
            {
                datetime(2008, 1, 1): datetime(2011, 1, 1),
                datetime(2008, 6, 30): datetime(2011, 1, 1),
                datetime(2008, 12, 31): datetime(2011, 1, 1),
                datetime(2005, 12, 30): datetime(2008, 1, 1),
                datetime(2005, 12, 31): datetime(2008, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(-1),
            {
                datetime(2007, 1, 1): datetime(2006, 1, 1),
                datetime(2007, 1, 15): datetime(2007, 1, 1),
                datetime(2008, 6, 30): datetime(2008, 1, 1),
                datetime(2008, 12, 31): datetime(2008, 1, 1),
                datetime(2006, 12, 29): datetime(2006, 1, 1),
                datetime(2006, 12, 30): datetime(2006, 1, 1),
                datetime(2007, 1, 1): datetime(2006, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(-2),
            {
                datetime(2007, 1, 1): datetime(2005, 1, 1),
                datetime(2008, 6, 30): datetime(2007, 1, 1),
                datetime(2008, 12, 31): datetime(2007, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(month=4),
            {
                datetime(2007, 4, 1): datetime(2008, 4, 1),
                datetime(2007, 4, 15): datetime(2008, 4, 1),
                datetime(2007, 3, 1): datetime(2007, 4, 1),
                datetime(2007, 12, 15): datetime(2008, 4, 1),
                datetime(2012, 1, 31): datetime(2012, 4, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(0, month=4),
            {
                datetime(2007, 4, 1): datetime(2007, 4, 1),
                datetime(2007, 3, 1): datetime(2007, 4, 1),
                datetime(2007, 12, 15): datetime(2008, 4, 1),
                datetime(2012, 1, 31): datetime(2012, 4, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(4, month=4),
            {
                datetime(2007, 4, 1): datetime(2011, 4, 1),
                datetime(2007, 4, 15): datetime(2011, 4, 1),
                datetime(2007, 3, 1): datetime(2010, 4, 1),
                datetime(2007, 12, 15): datetime(2011, 4, 1),
                datetime(2012, 1, 31): datetime(2015, 4, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(-1, month=4),
            {
                datetime(2007, 4, 1): datetime(2006, 4, 1),
                datetime(2007, 3, 1): datetime(2006, 4, 1),
                datetime(2007, 12, 15): datetime(2007, 4, 1),
                datetime(2012, 1, 31): datetime(2011, 4, 1),
            },
        )
    )

    offset_cases.append(
        (
            YearBegin(-3, month=4),
            {
                datetime(2007, 4, 1): datetime(2004, 4, 1),
                datetime(2007, 3, 1): datetime(2004, 4, 1),
                datetime(2007, 12, 15): datetime(2005, 4, 1),
                datetime(2012, 1, 31): datetime(2009, 4, 1),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [
        (YearBegin(), datetime(2007, 1, 3), False),
        (YearBegin(), datetime(2008, 1, 1), True),
        (YearBegin(), datetime(2006, 12, 31), False),
        (YearBegin(), datetime(2006, 1, 2), False),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)


class TestYearEnd(Base):
    _offset = YearEnd

    def test_misspecified(self):
        with pytest.raises(ValueError, match="Month must go from 1 to 12"):
            YearEnd(month=13)

    offset_cases = []
    offset_cases.append(
        (
            YearEnd(),
            {
                datetime(2008, 1, 1): datetime(2008, 12, 31),
                datetime(2008, 6, 30): datetime(2008, 12, 31),
                datetime(2008, 12, 31): datetime(2009, 12, 31),
                datetime(2005, 12, 30): datetime(2005, 12, 31),
                datetime(2005, 12, 31): datetime(2006, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(0),
            {
                datetime(2008, 1, 1): datetime(2008, 12, 31),
                datetime(2008, 6, 30): datetime(2008, 12, 31),
                datetime(2008, 12, 31): datetime(2008, 12, 31),
                datetime(2005, 12, 30): datetime(2005, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(-1),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 31),
                datetime(2008, 6, 30): datetime(2007, 12, 31),
                datetime(2008, 12, 31): datetime(2007, 12, 31),
                datetime(2006, 12, 29): datetime(2005, 12, 31),
                datetime(2006, 12, 30): datetime(2005, 12, 31),
                datetime(2007, 1, 1): datetime(2006, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(-2),
            {
                datetime(2007, 1, 1): datetime(2005, 12, 31),
                datetime(2008, 6, 30): datetime(2006, 12, 31),
                datetime(2008, 12, 31): datetime(2006, 12, 31),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [
        (YearEnd(), datetime(2007, 12, 31), True),
        (YearEnd(), datetime(2008, 1, 1), False),
        (YearEnd(), datetime(2006, 12, 31), True),
        (YearEnd(), datetime(2006, 12, 29), False),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)


class TestYearEndDiffMonth(Base):
    offset_cases = []
    offset_cases.append(
        (
            YearEnd(month=3),
            {
                datetime(2008, 1, 1): datetime(2008, 3, 31),
                datetime(2008, 2, 15): datetime(2008, 3, 31),
                datetime(2008, 3, 31): datetime(2009, 3, 31),
                datetime(2008, 3, 30): datetime(2008, 3, 31),
                datetime(2005, 3, 31): datetime(2006, 3, 31),
                datetime(2006, 7, 30): datetime(2007, 3, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(0, month=3),
            {
                datetime(2008, 1, 1): datetime(2008, 3, 31),
                datetime(2008, 2, 28): datetime(2008, 3, 31),
                datetime(2008, 3, 31): datetime(2008, 3, 31),
                datetime(2005, 3, 30): datetime(2005, 3, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(-1, month=3),
            {
                datetime(2007, 1, 1): datetime(2006, 3, 31),
                datetime(2008, 2, 28): datetime(2007, 3, 31),
                datetime(2008, 3, 31): datetime(2007, 3, 31),
                datetime(2006, 3, 29): datetime(2005, 3, 31),
                datetime(2006, 3, 30): datetime(2005, 3, 31),
                datetime(2007, 3, 1): datetime(2006, 3, 31),
            },
        )
    )

    offset_cases.append(
        (
            YearEnd(-2, month=3),
            {
                datetime(2007, 1, 1): datetime(2005, 3, 31),
                datetime(2008, 6, 30): datetime(2007, 3, 31),
                datetime(2008, 3, 31): datetime(2006, 3, 31),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [
        (YearEnd(month=3), datetime(2007, 3, 31), True),
        (YearEnd(month=3), datetime(2008, 1, 1), False),
        (YearEnd(month=3), datetime(2006, 3, 31), True),
        (YearEnd(month=3), datetime(2006, 3, 29), False),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)
