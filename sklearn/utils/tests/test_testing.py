import warnings
import unittest
import sys
import numpy as np
from scipy import sparse

from sklearn.utils.deprecation import deprecated
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.testing import (
    assert_true,
    assert_raises,
    assert_less,
    assert_greater,
    assert_less_equal,
    assert_greater_equal,
    assert_warns,
    assert_no_warnings,
    assert_equal,
    set_random_state,
    assert_raise_message,
    ignore_warnings,
    check_docstring_parameters,
    assert_allclose_dense_sparse)

from sklearn.utils.testing import SkipTest
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def test_assert_less():
    assert_less(0, 1)
    assert_raises(AssertionError, assert_less, 1, 0)


def test_assert_greater():
    assert_greater(1, 0)
    assert_raises(AssertionError, assert_greater, 0, 1)


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


def test_assert_allclose_dense_sparse():
    x = np.arange(9).reshape(3, 3)
    msg = "Not equal to tolerance "
    y = sparse.csc_matrix(x)
    for X in [x, y]:
        # basic compare
        assert_raise_message(AssertionError, msg, assert_allclose_dense_sparse,
                             X, X * 2)
        assert_allclose_dense_sparse(X, X)

    assert_raise_message(ValueError, "Can only compare two sparse",
                         assert_allclose_dense_sparse, x, y)

    A = sparse.diags(np.ones(5), offsets=0).tocsr()
    B = sparse.csr_matrix(np.ones((1, 5)))

    assert_raise_message(AssertionError, "Arrays are not equal",
                         assert_allclose_dense_sparse, B, A)


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
# `clean_warning_registry()` is called internally by assert_warns
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


# Tests for docstrings:

def f_ok(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_bad_sections(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Results
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_bad_order(b, a):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_missing(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_check_param_definition(a, b, c, d):
    """Function f

    Parameters
    ----------
    a: int
        Parameter a
    b:
        Parameter b
    c :
        Parameter c
    d:int
        Parameter d
    """
    return a + b + c + d


class Klass(object):
    def f_missing(self, X, y):
        pass

    def f_bad_sections(self, X, y):
        """Function f

        Parameter
        ----------
        a : int
            Parameter a
        b : float
            Parameter b

        Results
        -------
        c : list
            Parameter c
        """
        pass


class MockEst(object):
    def __init__(self):
        """MockEstimator"""
    def fit(self, X, y):
        return X

    def predict(self, X):
        return X

    def predict_proba(self, X):
        return X

    def score(self, X):
        return 1.


class MockMetaEstimator(object):
    def __init__(self, delegate):
        """MetaEstimator to check if doctest on delegated methods work.

        Parameters
        ---------
        delegate : estimator
            Delegated estimator.
        """
        self.delegate = delegate

    @if_delegate_has_method(delegate=('delegate'))
    def predict(self, X):
        """This is available only if delegate has predict.

        Parameters
        ----------
        y : ndarray
            Parameter y
        """
        return self.delegate.predict(X)

    @deprecated("Testing a deprecated delegated method")
    @if_delegate_has_method(delegate=('delegate'))
    def score(self, X):
        """This is available only if delegate has score.

        Parameters
        ---------
        y : ndarray
            Parameter y
        """

    @if_delegate_has_method(delegate=('delegate'))
    def predict_proba(self, X):
        """This is available only if delegate has predict_proba.

        Parameters
        ---------
        X : ndarray
            Parameter X
        """
        return X

    @deprecated('Testing deprecated function with incorrect params')
    @if_delegate_has_method(delegate=('delegate'))
    def predict_log_proba(self, X):
        """This is available only if delegate has predict_proba.

        Parameters
        ---------
        y : ndarray
            Parameter X
        """
        return X

    @deprecated('Testing deprecated function with wrong params')
    @if_delegate_has_method(delegate=('delegate'))
    def fit(self, X, y):
        """Incorrect docstring but should not be tested"""


def test_check_docstring_parameters():
    try:
        import numpydoc  # noqa
        assert sys.version_info >= (3, 5)
    except (ImportError, AssertionError):
        raise SkipTest(
            "numpydoc is required to test the docstrings")

    incorrect = check_docstring_parameters(f_ok)
    assert_equal(incorrect, [])
    incorrect = check_docstring_parameters(f_ok, ignore=['b'])
    assert_equal(incorrect, [])
    incorrect = check_docstring_parameters(f_missing, ignore=['b'])
    assert_equal(incorrect, [])
    assert_raise_message(RuntimeError, 'Unknown section Results',
                         check_docstring_parameters, f_bad_sections)
    assert_raise_message(RuntimeError, 'Unknown section Parameter',
                         check_docstring_parameters, Klass.f_bad_sections)

    messages = ["a != b", "arg mismatch: ['b']", "arg mismatch: ['X', 'y']",
                "predict y != X",
                "predict_proba arg mismatch: ['X']",
                "predict_log_proba arg mismatch: ['X']",
                "score arg mismatch: ['X']",
                ".fit arg mismatch: ['X', 'y']"]

    mock_meta = MockMetaEstimator(delegate=MockEst())

    for mess, f in zip(messages,
                       [f_bad_order, f_missing, Klass.f_missing,
                        mock_meta.predict, mock_meta.predict_proba,
                        mock_meta.predict_log_proba,
                        mock_meta.score, mock_meta.fit]):
        incorrect = check_docstring_parameters(f)
        assert_true(len(incorrect) >= 1)
        assert_true(mess in incorrect[0],
                    '"%s" not in "%s"' % (mess, incorrect[0]))

    incorrect = check_docstring_parameters(f_check_param_definition)
    assert_equal(
        incorrect,
        ['sklearn.utils.tests.test_testing.f_check_param_definition There was '
         'no space between the param name and colon ("a: int")',
         'sklearn.utils.tests.test_testing.f_check_param_definition There was '
         'no space between the param name and colon ("b:")',
         'sklearn.utils.tests.test_testing.f_check_param_definition Incorrect '
         'type definition for param: "c " (type definition was "")',
         'sklearn.utils.tests.test_testing.f_check_param_definition There was '
         'no space between the param name and colon ("d:int")'])
