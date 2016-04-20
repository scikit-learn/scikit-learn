import warnings
import unittest
import sys

from nose.tools import assert_raises

from sklearn.utils.testing import (
    _assert_less,
    _assert_greater,
    assert_less_equal,
    assert_greater_equal,
    assert_warns,
    assert_no_warnings,
    assert_equal,
    set_random_state,
    assert_raise_message,
    ignore_warnings)

from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

try:
    from nose.tools import assert_less

    def test_assert_less():
        # Check that the nose implementation of assert_less gives the
        # same thing as the scikit's
        assert_less(0, 1)
        _assert_less(0, 1)
        assert_raises(AssertionError, assert_less, 1, 0)
        assert_raises(AssertionError, _assert_less, 1, 0)

except ImportError:
    pass

try:
    from nose.tools import assert_greater

    def test_assert_greater():
        # Check that the nose implementation of assert_less gives the
        # same thing as the scikit's
        assert_greater(1, 0)
        _assert_greater(1, 0)
        assert_raises(AssertionError, assert_greater, 0, 1)
        assert_raises(AssertionError, _assert_greater, 0, 1)

except ImportError:
    pass


def test_assert_less_equal():
    assert_less_equal(0, 1)
    assert_less_equal(1, 1)
    assert_raises(AssertionError, assert_less_equal, 1, 0)


def test_assert_greater_equal():
    assert_greater_equal(1, 0)
    assert_greater_equal(1, 1)
    assert_raises(AssertionError, assert_greater_equal, 0, 1)


def test_set_random_state():
    lda = LinearDiscriminantAnalysis()
    tree = DecisionTreeClassifier()
    # Linear Discriminant Analysis doesn't have random state: smoke test
    set_random_state(lda, 3)
    set_random_state(tree, 3)
    assert_equal(tree.random_state, 3)


def test_assert_raise_message():
    def _raise_ValueError(message):
        raise ValueError(message)

    def _no_raise():
        pass

    assert_raise_message(ValueError, "test",
                         _raise_ValueError, "test")

    assert_raises(AssertionError,
                  assert_raise_message, ValueError, "something else",
                  _raise_ValueError, "test")

    assert_raises(ValueError,
                  assert_raise_message, TypeError, "something else",
                  _raise_ValueError, "test")

    assert_raises(AssertionError,
                  assert_raise_message, ValueError, "test",
                  _no_raise)

    # multiple exceptions in a tuple
    assert_raises(AssertionError,
                  assert_raise_message, (ValueError, AttributeError),
                  "test", _no_raise)


def test_ignore_warning():
    # This check that ignore_warning decorateur and context manager are working
    # as expected
    def _warning_function():
        warnings.warn("deprecation warning", DeprecationWarning)

    def _multiple_warning_function():
        warnings.warn("deprecation warning", DeprecationWarning)
        warnings.warn("deprecation warning")

    # Check the function directly
    assert_no_warnings(ignore_warnings(_warning_function))
    assert_no_warnings(ignore_warnings(_warning_function,
                                       category=DeprecationWarning))
    assert_warns(DeprecationWarning, ignore_warnings(_warning_function,
                                                     category=UserWarning))
    assert_warns(UserWarning,
                 ignore_warnings(_multiple_warning_function,
                                 category=DeprecationWarning))
    assert_warns(DeprecationWarning,
                 ignore_warnings(_multiple_warning_function,
                                 category=UserWarning))
    assert_no_warnings(ignore_warnings(_warning_function,
                                       category=(DeprecationWarning,
                                                 UserWarning)))

    # Check the decorator
    @ignore_warnings
    def decorator_no_warning():
        _warning_function()
        _multiple_warning_function()

    @ignore_warnings(category=(DeprecationWarning, UserWarning))
    def decorator_no_warning_multiple():
        _multiple_warning_function()

    @ignore_warnings(category=DeprecationWarning)
    def decorator_no_deprecation_warning():
        _warning_function()

    @ignore_warnings(category=UserWarning)
    def decorator_no_user_warning():
        _warning_function()

    @ignore_warnings(category=DeprecationWarning)
    def decorator_no_deprecation_multiple_warning():
        _multiple_warning_function()

    @ignore_warnings(category=UserWarning)
    def decorator_no_user_multiple_warning():
        _multiple_warning_function()

    assert_no_warnings(decorator_no_warning)
    assert_no_warnings(decorator_no_warning_multiple)
    assert_no_warnings(decorator_no_deprecation_warning)
    assert_warns(DeprecationWarning, decorator_no_user_warning)
    assert_warns(UserWarning, decorator_no_deprecation_multiple_warning)
    assert_warns(DeprecationWarning, decorator_no_user_multiple_warning)

    # Check the context manager
    def context_manager_no_warning():
        with ignore_warnings():
            _warning_function()

    def context_manager_no_warning_multiple():
        with ignore_warnings(category=(DeprecationWarning, UserWarning)):
            _multiple_warning_function()

    def context_manager_no_deprecation_warning():
        with ignore_warnings(category=DeprecationWarning):
            _warning_function()

    def context_manager_no_user_warning():
        with ignore_warnings(category=UserWarning):
            _warning_function()

    def context_manager_no_deprecation_multiple_warning():
        with ignore_warnings(category=DeprecationWarning):
            _multiple_warning_function()

    def context_manager_no_user_multiple_warning():
        with ignore_warnings(category=UserWarning):
            _multiple_warning_function()

    assert_no_warnings(context_manager_no_warning)
    assert_no_warnings(context_manager_no_warning_multiple)
    assert_no_warnings(context_manager_no_deprecation_warning)
    assert_warns(DeprecationWarning, context_manager_no_user_warning)
    assert_warns(UserWarning, context_manager_no_deprecation_multiple_warning)
    assert_warns(DeprecationWarning, context_manager_no_user_multiple_warning)


# This class is inspired from numpy 1.7 with an alteration to check
# the reset warning filters after calls to assert_warns.
# This assert_warns behavior is specific to scikit-learn because
#`clean_warning_registry()` is called internally by assert_warns
# and clears all previous filters.
class TestWarns(unittest.TestCase):
    def test_warn(self):
        def f():
            warnings.warn("yo")
            return 3

        # Test that assert_warns is not impacted by externally set
        # filters and is reset internally.
        # This is because `clean_warning_registry()` is called internally by
        # assert_warns and clears all previous filters.
        warnings.simplefilter("ignore", UserWarning)
        assert_equal(assert_warns(UserWarning, f), 3)

        # Test that the warning registry is empty after assert_warns
        assert_equal(sys.modules['warnings'].filters, [])

        assert_raises(AssertionError, assert_no_warnings, f)
        assert_equal(assert_no_warnings(lambda x: x, 1), 1)

    def test_warn_wrong_warning(self):
        def f():
            warnings.warn("yo", DeprecationWarning)

        failed = False
        filters = sys.modules['warnings'].filters[:]
        try:
            try:
                # Should raise an AssertionError
                assert_warns(UserWarning, f)
                failed = True
            except AssertionError:
                pass
        finally:
            sys.modules['warnings'].filters = filters

        if failed:
            raise AssertionError("wrong warning caught by assert_warn")
