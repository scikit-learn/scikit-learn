# -*- coding: utf-8 -*-
"""
Assertion helpers for offsets tests
"""


def assert_offset_equal(offset, base, expected):
    actual = offset + base
    actual_swapped = base + offset
    actual_apply = offset.apply(base)
    try:
        assert actual == expected
        assert actual_swapped == expected
        assert actual_apply == expected
    except AssertionError:
        raise AssertionError("\nExpected: %s\nActual: %s\nFor Offset: %s)"
                             "\nAt Date: %s" %
                             (expected, actual, offset, base))


def assert_onOffset(offset, date, expected):
    actual = offset.onOffset(date)
    assert actual == expected, ("\nExpected: %s\nActual: %s\nFor Offset: %s)"
                                "\nAt Date: %s" %
                                (expected, actual, offset, date))
