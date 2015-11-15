# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

from sklearn.utils import testing as testing
from sklearn import basic_checks as bc

# Test private methods


def test__get_wrapper():
    """
    Tests if the function `__get_wrapper` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `__get_wrapper` failed the test.
    """
    wrap_fct = bc.__get_wrapper
    wrap_raises = wrap_fct(
        launch_exception=True, exception_type=AssertionError)
    wrap_no_raises = wrap_fct(launch_exception=False)
    testing.assert_raise_message(
        AssertionError, "Error", wrap_raises, False, "Error")
    testing.assert_true(wrap_raises(True))
    testing.assert_true(wrap_no_raises(True))
    testing.assert_false(wrap_no_raises(False))

# Tests public methods


def test_check_all_in():
    """
    Tests if the function `check_all_in` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_all_in` failed the test.
    """
    testing.assert_true(bc.check_all_in(42, 42, both=True))
    testing.assert_true(bc.check_all_in("test", ("test",), both=True))
    testing.assert_true(
        bc.check_all_in((42, "test"), ("test", 42,), both=True))
    testing.assert_false(bc.check_all_in(42, (42, 41), both=True))
    testing.assert_false(bc.check_all_in(42, (41, 42), both=True))
    testing.assert_false(bc.check_all_in((41, 42), 42, both=True))
    testing.assert_false(bc.check_all_in((42, 41), 42, both=True))
    testing.assert_true(
        bc.check_all_in(42, 42, both=False,
                        first_in_second=True, second_in_first=True))
    testing.assert_true(bc.check_all_in(
        "test", ("test",), both=False, first_in_second=True,
        second_in_first=True))
    testing.assert_true(bc.check_all_in(
        (42, "test"), ("test", 42,), both=False,
        first_in_second=True, second_in_first=True))
    testing.assert_false(bc.check_all_in(
        42, (42, 41), both=False, first_in_second=True, second_in_first=True))
    testing.assert_false(bc.check_all_in(
        42, (41, 42), both=False, first_in_second=True, second_in_first=True))
    testing.assert_false(bc.check_all_in(
        (41, 42), 42, both=False, first_in_second=True, second_in_first=True))
    testing.assert_false(bc.check_all_in(
        (42, 41), 42, both=False, first_in_second=True, second_in_first=True))

    testing.assert_true(bc.check_all_in(42, 42, first_in_second=True))
    testing.assert_true(
        bc.check_all_in("test", ("test",), first_in_second=True))
    testing.assert_true(
        bc.check_all_in((42, "test"), ("test", 42,), first_in_second=True))
    testing.assert_true(bc.check_all_in(42, (42, 41), first_in_second=True))
    testing.assert_true(bc.check_all_in(42, (41, 42), first_in_second=True))
    testing.assert_false(bc.check_all_in((41, 42), 42, first_in_second=True))
    testing.assert_false(bc.check_all_in((42, 41), 42, first_in_second=True))

    testing.assert_true(bc.check_all_in(42, 42, second_in_first=True))
    testing.assert_true(
        bc.check_all_in("test", ("test",), second_in_first=True))
    testing.assert_true(
        bc.check_all_in((42, "test"), ("test", 42,), second_in_first=True))
    testing.assert_false(bc.check_all_in(42, (42, 41), second_in_first=True))
    testing.assert_false(bc.check_all_in(42, (41, 42), second_in_first=True))
    testing.assert_true(bc.check_all_in((41, 42), 42, second_in_first=True))
    testing.assert_true(bc.check_all_in((42, 41), 42, second_in_first=True))

    testing.assert_true(bc.check_all_in(
        42, 42, both=False, first_in_second=False, second_in_first=False))
    testing.assert_true(bc.check_all_in([], 42, first_in_second=True))
    testing.assert_false(bc.check_all_in([], 42, second_in_first=True))


def test_check_is_binary_01():
    """
    Tests if the function `check_is_binary_01` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_binary_01` failed the test.
    """
    testing.assert_true(bc.check_is_binary_01(
        (0,), launch_exception=False))
    testing.assert_true(bc.check_is_binary_01(
        [1], launch_exception=False))
    testing.assert_true(bc.check_is_binary_01(
        [[0, 0], [0, 0]], launch_exception=False))
    testing.assert_true(bc.check_is_binary_01(
        [[1, 0], [0, 1]], launch_exception=False))
    testing.assert_false(bc.check_is_binary_01(
        0, launch_exception=False))
    testing.assert_false(bc.check_is_binary_01(
        2, launch_exception=False))
    testing.assert_false(bc.check_is_binary_01(
        [0, 2], launch_exception=False))
    testing.assert_false(bc.check_is_binary_01(
        [[[0, 2], [2, 0]], [0, 0]], launch_exception=False))
    testing.assert_false(bc.check_is_binary_01(
        [[[0, 2], [1, 0]], [0, 0]], launch_exception=False))
    testing.assert_raise_message(
        bc.CheckFalseError, "Error:",
        bc.check_is_binary_01, None, launch_exception=True)


def test_check_is_in():
    """
    Tests if the function `check_is_in` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_in` failed the test.
    """
    testing.assert_true(bc.check_is_in(None, None))
    testing.assert_true(bc.check_is_in(42, 42))
    testing.assert_true(bc.check_is_in(42, [42]))
    testing.assert_true(bc.check_is_in([42, 30], [[42, 30], [54, 65]]))
    testing.assert_true(bc.check_is_in(42, [5, "rnd", 42]))
    testing.assert_false(bc.check_is_in(42, []))
    testing.assert_false(bc.check_is_in(42, "0024156482"))
    testing.assert_true(bc.check_is_in(42, "42", stringed=True))
    testing.assert_true(bc.check_is_in(42, "0256420468", stringed=True))
    testing.assert_false(bc.check_is_in(None, [42, "rnd", "", False]))
    testing.assert_false(bc.check_is_in([], 42))
    testing.assert_false(bc.check_is_in([], [42]))
    testing.assert_raise_message(
        bc.CheckFalseError, "Error:",
        bc.check_is_in, None, [], launch_exception=True)


def test_check_is_in_list():
    """
    Tests if the function `check_is_in_list` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_in_list` failed the test.
    """
    testing.assert_false(bc.check_is_in_list(42, None, False))
    testing.assert_raise_message(
        bc.CheckFalseError, "Error:", bc.check_is_in_list, 42, (41, 43), True)
    testing.assert_true(bc.check_is_in_list(42, (42,), False))


def test_check_is_in_range():
    """
    Tests if the function `check_is_in_range` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_in_range` failed the test.
    """
    testing.assert_false(bc.check_is_in_range(
        42, lower_bound=42, low_exclusive=True, launch_exception=False))
    testing.assert_true(bc.check_is_in_range(
        42, lower_bound=42, low_exclusive=False, launch_exception=False))
    testing.assert_false(bc.check_is_in_range(
        42, higher_bound=42, high_exclusive=True, launch_exception=False))
    testing.assert_true(bc.check_is_in_range(
        42, higher_bound=42, high_exclusive=False, launch_exception=False))
    testing.assert_raise_message(
        bc.CheckFalseError, "Error 42",
        bc.check_is_in_range, elem=42, lower_bound=43,
        higher_bound=41, msg="Error 42", launch_exception=True)


def test_check_is_iterable():
    """
    Tests if the function `check_is_iterable` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_iterable` failed the test.
    """
    testing.assert_true(bc.check_is_iterable([42]))
    testing.assert_false(bc.check_is_iterable(42))
    testing.assert_true(bc.check_is_iterable(iter(range(10))))
    testing.assert_true(bc.check_is_iterable(["42"]))
    testing.assert_false(bc.check_is_iterable("42", string_allowed=False))
    testing.assert_true(bc.check_is_iterable("42", string_allowed=True))
    testing.assert_raise_message(
        bc.CheckFalseError, "Error",
        bc.check_is_iterable, "test",
        string_allowed=False, launch_exception=True)
    itl = [10, 20, 30]
    it = iter(itl)
    testing.assert_true(bc.check_is_iterable(it))
    testing.assert_array_equal(list(it), itl)


def test_check_is_list_or_tuple():
    """
    Tests if the function `check_is_list_or_tuple` behaves correctly.

    Raises
    ------
    AssertionError
        If the method `check_is_list_or_tuple` failed the test.
    """
    testing.assert_false(bc.check_is_list_or_tuple(None, False, None))
    testing.assert_raise_message(
        bc.CheckFalseError, "not a list", bc.check_is_list_or_tuple,
        42, True, None)
    testing.assert_false(
        bc.check_is_list_or_tuple((42,), False, (lambda x: (x != 1))))
    testing.assert_raise_message(
        bc.CheckFalseError, "incorrect size", bc.check_is_list_or_tuple,
        (42,), True, (lambda x: (x != 1)))
    testing.assert_true(bc.check_is_list_or_tuple(list((42,)), False, None))
    testing.assert_true(bc.check_is_list_or_tuple((42,), False, None))
    testing.assert_true(
        bc.check_is_list_or_tuple(list((42,)), False, (lambda x: (x == 1))))
