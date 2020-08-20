import warnings
import unittest
import sys
import os
import atexit

import numpy as np

from scipy import sparse

import pytest

from sklearn.utils.deprecation import deprecated
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils._testing import (
    assert_raises,
    assert_warns,
    assert_no_warnings,
    set_random_state,
    assert_raise_message,
    ignore_warnings,
    check_docstring_parameters,
    assert_allclose_dense_sparse,
    assert_raises_regex,
    TempMemmap,
    create_memmap_backed_data,
    _delete_folder,
    _convert_container)

from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def test_set_random_state():
    lda = LinearDiscriminantAnalysis()
    tree = DecisionTreeClassifier()
    # Linear Discriminant Analysis doesn't have random state: smoke test
    set_random_state(lda, 3)
    set_random_state(tree, 3)
    assert tree.random_state == 3


def test_assert_allclose_dense_sparse():
    x = np.arange(9).reshape(3, 3)
    msg = "Not equal to tolerance "
    y = sparse.csc_matrix(x)
    for X in [x, y]:
        # basic compare
        with pytest.raises(AssertionError, match=msg):
            assert_allclose_dense_sparse(X, X*2)
        assert_allclose_dense_sparse(X, X)

    with pytest.raises(ValueError, match="Can only compare two sparse"):
        assert_allclose_dense_sparse(x, y)

    A = sparse.diags(np.ones(5), offsets=0).tocsr()
    B = sparse.csr_matrix(np.ones((1, 5)))
    with pytest.raises(AssertionError, match="Arrays are not equal"):
        assert_allclose_dense_sparse(B, A)


def test_assert_raises_msg():
    with assert_raises_regex(AssertionError, 'Hello world'):
        with assert_raises(ValueError, msg='Hello world'):
            pass


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
                                 category=FutureWarning))
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

    # Check that passing warning class as first positional argument
    warning_class = UserWarning
    match = "'obj' should be a callable.+you should use 'category=UserWarning'"

    with pytest.raises(ValueError, match=match):
        silence_warnings_func = ignore_warnings(warning_class)(
            _warning_function)
        silence_warnings_func()

    with pytest.raises(ValueError, match=match):
        @ignore_warnings(warning_class)
        def test():
            pass


class TestWarns(unittest.TestCase):
    def test_warn(self):
        def f():
            warnings.warn("yo")
            return 3

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            filters_orig = warnings.filters[:]
            assert assert_warns(UserWarning, f) == 3
            # test that assert_warns doesn't have side effects on warnings
            # filters
            assert warnings.filters == filters_orig
        with pytest.raises(AssertionError):
            assert_no_warnings(f)
        assert assert_no_warnings(lambda x: x, 1) == 1

    def test_warn_wrong_warning(self):
        def f():
            warnings.warn("yo", FutureWarning)

        failed = False
        filters = sys.modules['warnings'].filters[:]
        try:
            try:
                # Should raise an AssertionError

                # assert_warns has a special handling of "FutureWarning" that
                # pytest.warns does not have
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


def f_too_many_param_docstring(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : int
        Parameter b
    c : int
        Parameter c

    Returns
    -------
    d : list
        Parameter c
    """
    d = a + b
    return d


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


def f_check_param_definition(a, b, c, d, e):
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
    e
        No typespec is allowed without colon
    """
    return a + b + c + d


class Klass:
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


class MockEst:
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


class MockMetaEstimator:
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

    @if_delegate_has_method(delegate=('delegate'))
    @deprecated("Testing a deprecated delegated method")
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

    @deprecated('Testing deprecated function with wrong params')
    def fit(self, X, y):
        """Incorrect docstring but should not be tested"""


def test_check_docstring_parameters():
    pytest.importorskip('numpydoc',
                        reason="numpydoc is required to test the docstrings")

    incorrect = check_docstring_parameters(f_ok)
    assert incorrect == []
    incorrect = check_docstring_parameters(f_ok, ignore=['b'])
    assert incorrect == []
    incorrect = check_docstring_parameters(f_missing, ignore=['b'])
    assert incorrect == []
    with pytest.raises(RuntimeError, match="Unknown section Results"):
        check_docstring_parameters(f_bad_sections)
    with pytest.raises(RuntimeError, match="Unknown section Parameter"):
        check_docstring_parameters(Klass.f_bad_sections)

    incorrect = check_docstring_parameters(f_check_param_definition)
    assert (
        incorrect == [
            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('a: int')",

            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('b:')",

            "sklearn.utils.tests.test_testing.f_check_param_definition "
            "Parameter 'c :' has an empty type spec. Remove the colon",

            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('d:int')",
        ])

    messages = [
            ["In function: sklearn.utils.tests.test_testing.f_bad_order",
             "There's a parameter name mismatch in function docstring w.r.t."
             " function signature, at index 0 diff: 'b' != 'a'",
             "Full diff:",
             "- ['b', 'a']",
             "+ ['a', 'b']"],

            ["In function: " +
                "sklearn.utils.tests.test_testing.f_too_many_param_docstring",
             "Parameters in function docstring have more items w.r.t. function"
             " signature, first extra item: c",
             "Full diff:",
             "- ['a', 'b']",
             "+ ['a', 'b', 'c']",
             "?          +++++"],

            ["In function: sklearn.utils.tests.test_testing.f_missing",
             "Parameters in function docstring have less items w.r.t. function"
             " signature, first missing item: b",
             "Full diff:",
             "- ['a', 'b']",
             "+ ['a']"],

            ["In function: sklearn.utils.tests.test_testing.Klass.f_missing",
             "Parameters in function docstring have less items w.r.t. function"
             " signature, first missing item: X",
             "Full diff:",
             "- ['X', 'y']",
             "+ []"],

            ["In function: " +
             "sklearn.utils.tests.test_testing.MockMetaEstimator.predict",
             "There's a parameter name mismatch in function docstring w.r.t."
             " function signature, at index 0 diff: 'X' != 'y'",
             "Full diff:",
             "- ['X']",
             "?   ^",
             "+ ['y']",
             "?   ^"],

            ["In function: " +
             "sklearn.utils.tests.test_testing.MockMetaEstimator."
             + "predict_proba",
             "Parameters in function docstring have less items w.r.t. function"
             " signature, first missing item: X",
             "Full diff:",
             "- ['X']",
             "+ []"],

            ["In function: " +
                "sklearn.utils.tests.test_testing.MockMetaEstimator.score",
             "Parameters in function docstring have less items w.r.t. function"
             " signature, first missing item: X",
             "Full diff:",
             "- ['X']",
             "+ []"],

            ["In function: " +
                "sklearn.utils.tests.test_testing.MockMetaEstimator.fit",
             "Parameters in function docstring have less items w.r.t. function"
             " signature, first missing item: X",
             "Full diff:",
             "- ['X', 'y']",
             "+ []"],

            ]

    mock_meta = MockMetaEstimator(delegate=MockEst())

    for msg, f in zip(messages,
                      [f_bad_order,
                       f_too_many_param_docstring,
                       f_missing,
                       Klass.f_missing,
                       mock_meta.predict,
                       mock_meta.predict_proba,
                       mock_meta.score,
                       mock_meta.fit]):
        incorrect = check_docstring_parameters(f)
        assert msg == incorrect, ('\n"%s"\n not in \n"%s"' % (msg, incorrect))


class RegistrationCounter:
    def __init__(self):
        self.nb_calls = 0

    def __call__(self, to_register_func):
        self.nb_calls += 1
        assert to_register_func.func is _delete_folder


def check_memmap(input_array, mmap_data, mmap_mode='r'):
    assert isinstance(mmap_data, np.memmap)
    writeable = mmap_mode != 'r'
    assert mmap_data.flags.writeable is writeable
    np.testing.assert_array_equal(input_array, mmap_data)


def test_tempmemmap(monkeypatch):
    registration_counter = RegistrationCounter()
    monkeypatch.setattr(atexit, 'register', registration_counter)

    input_array = np.ones(3)
    with TempMemmap(input_array) as data:
        check_memmap(input_array, data)
        temp_folder = os.path.dirname(data.filename)
    if os.name != 'nt':
        assert not os.path.exists(temp_folder)
    assert registration_counter.nb_calls == 1

    mmap_mode = 'r+'
    with TempMemmap(input_array, mmap_mode=mmap_mode) as data:
        check_memmap(input_array, data, mmap_mode=mmap_mode)
        temp_folder = os.path.dirname(data.filename)
    if os.name != 'nt':
        assert not os.path.exists(temp_folder)
    assert registration_counter.nb_calls == 2


def test_create_memmap_backed_data(monkeypatch):
    registration_counter = RegistrationCounter()
    monkeypatch.setattr(atexit, 'register', registration_counter)

    input_array = np.ones(3)
    data = create_memmap_backed_data(input_array)
    check_memmap(input_array, data)
    assert registration_counter.nb_calls == 1

    data, folder = create_memmap_backed_data(input_array,
                                             return_folder=True)
    check_memmap(input_array, data)
    assert folder == os.path.dirname(data.filename)
    assert registration_counter.nb_calls == 2

    mmap_mode = 'r+'
    data = create_memmap_backed_data(input_array, mmap_mode=mmap_mode)
    check_memmap(input_array, data, mmap_mode)
    assert registration_counter.nb_calls == 3

    input_list = [input_array, input_array + 1, input_array + 2]
    mmap_data_list = create_memmap_backed_data(input_list)
    for input_array, data in zip(input_list, mmap_data_list):
        check_memmap(input_array, data)
    assert registration_counter.nb_calls == 4


@pytest.mark.parametrize(
    "constructor_name, container_type",
    [('list', list),
     ('tuple', tuple),
     ('array', np.ndarray),
     ('sparse', sparse.csr_matrix),
     ('dataframe', pytest.importorskip('pandas').DataFrame),
     ('series', pytest.importorskip('pandas').Series),
     ('index', pytest.importorskip('pandas').Index),
     ('slice', slice)]
)
def test_convert_container(constructor_name, container_type):
    container = [0, 1]
    assert isinstance(_convert_container(container, constructor_name),
                      container_type)
