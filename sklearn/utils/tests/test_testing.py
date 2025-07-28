import atexit
import os
import re
import warnings

import numpy as np
import pytest
from scipy import sparse

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.tree import DecisionTreeClassifier
from sklearn.utils._testing import (
    TempMemmap,
    _convert_container,
    _delete_folder,
    _get_warnings_filters_info_list,
    assert_allclose,
    assert_allclose_dense_sparse,
    assert_docstring_consistency,
    assert_run_python_script_without_output,
    check_docstring_parameters,
    create_memmap_backed_data,
    ignore_warnings,
    raises,
    set_random_state,
    skip_if_no_numpydoc,
    turn_warnings_into_errors,
)
from sklearn.utils.deprecation import deprecated
from sklearn.utils.fixes import (
    _IS_WASM,
    CSC_CONTAINERS,
    CSR_CONTAINERS,
)
from sklearn.utils.metaestimators import available_if


def test_set_random_state():
    lda = LinearDiscriminantAnalysis()
    tree = DecisionTreeClassifier()
    # Linear Discriminant Analysis doesn't have random state: smoke test
    set_random_state(lda, 3)
    set_random_state(tree, 3)
    assert tree.random_state == 3


@pytest.mark.parametrize("csr_container", CSC_CONTAINERS)
def test_assert_allclose_dense_sparse(csr_container):
    x = np.arange(9).reshape(3, 3)
    msg = "Not equal to tolerance "
    y = csr_container(x)
    for X in [x, y]:
        # basic compare
        with pytest.raises(AssertionError, match=msg):
            assert_allclose_dense_sparse(X, X * 2)
        assert_allclose_dense_sparse(X, X)

    with pytest.raises(ValueError, match="Can only compare two sparse"):
        assert_allclose_dense_sparse(x, y)

    A = sparse.diags(np.ones(5), offsets=0).tocsr()
    B = csr_container(np.ones((1, 5)))
    with pytest.raises(AssertionError, match="Arrays are not equal"):
        assert_allclose_dense_sparse(B, A)


def test_ignore_warning():
    # This check that ignore_warning decorator and context manager are working
    # as expected
    def _warning_function():
        warnings.warn("deprecation warning", DeprecationWarning)

    def _multiple_warning_function():
        warnings.warn("deprecation warning", DeprecationWarning)
        warnings.warn("deprecation warning")

    # Check the function directly
    with warnings.catch_warnings():
        warnings.simplefilter("error")

        ignore_warnings(_warning_function)
        ignore_warnings(_warning_function, category=DeprecationWarning)

    with pytest.warns(DeprecationWarning):
        ignore_warnings(_warning_function, category=UserWarning)()

    with pytest.warns() as record:
        ignore_warnings(_multiple_warning_function, category=FutureWarning)()
    assert len(record) == 2
    assert isinstance(record[0].message, DeprecationWarning)
    assert isinstance(record[1].message, UserWarning)

    with pytest.warns() as record:
        ignore_warnings(_multiple_warning_function, category=UserWarning)()
    assert len(record) == 1
    assert isinstance(record[0].message, DeprecationWarning)

    with warnings.catch_warnings():
        warnings.simplefilter("error")

        ignore_warnings(_warning_function, category=(DeprecationWarning, UserWarning))

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

    with warnings.catch_warnings():
        warnings.simplefilter("error")

        decorator_no_warning()
        decorator_no_warning_multiple()
        decorator_no_deprecation_warning()

    with pytest.warns(DeprecationWarning):
        decorator_no_user_warning()
    with pytest.warns(UserWarning):
        decorator_no_deprecation_multiple_warning()
    with pytest.warns(DeprecationWarning):
        decorator_no_user_multiple_warning()

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

    with warnings.catch_warnings():
        warnings.simplefilter("error")

        context_manager_no_warning()
        context_manager_no_warning_multiple()
        context_manager_no_deprecation_warning()

    with pytest.warns(DeprecationWarning):
        context_manager_no_user_warning()
    with pytest.warns(UserWarning):
        context_manager_no_deprecation_multiple_warning()
    with pytest.warns(DeprecationWarning):
        context_manager_no_user_multiple_warning()

    # Check that passing warning class as first positional argument
    warning_class = UserWarning
    match = "'obj' should be a callable.+you should use 'category=UserWarning'"

    with pytest.raises(ValueError, match=match):
        silence_warnings_func = ignore_warnings(warning_class)(_warning_function)
        silence_warnings_func()

    with pytest.raises(ValueError, match=match):

        @ignore_warnings(warning_class)
        def test():
            pass


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
        This is parsed correctly in numpydoc 1.2
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
        ---------
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
        return 1.0


class MockMetaEstimator:
    def __init__(self, delegate):
        """MetaEstimator to check if doctest on delegated methods work.

        Parameters
        ---------
        delegate : estimator
            Delegated estimator.
        """
        self.delegate = delegate

    @available_if(lambda self: hasattr(self.delegate, "predict"))
    def predict(self, X):
        """This is available only if delegate has predict.

        Parameters
        ----------
        y : ndarray
            Parameter y
        """
        return self.delegate.predict(X)

    @available_if(lambda self: hasattr(self.delegate, "score"))
    @deprecated("Testing a deprecated delegated method")
    def score(self, X):
        """This is available only if delegate has score.

        Parameters
        ---------
        y : ndarray
            Parameter y
        """

    @available_if(lambda self: hasattr(self.delegate, "predict_proba"))
    def predict_proba(self, X):
        """This is available only if delegate has predict_proba.

        Parameters
        ---------
        X : ndarray
            Parameter X
        """
        return X

    @deprecated("Testing deprecated function with wrong params")
    def fit(self, X, y):
        """Incorrect docstring but should not be tested"""


@skip_if_no_numpydoc
def test_check_docstring_parameters():
    incorrect = check_docstring_parameters(f_ok)
    assert incorrect == []
    incorrect = check_docstring_parameters(f_ok, ignore=["b"])
    assert incorrect == []
    incorrect = check_docstring_parameters(f_missing, ignore=["b"])
    assert incorrect == []
    with pytest.raises(RuntimeError, match="Unknown section Results"):
        check_docstring_parameters(f_bad_sections)
    with pytest.raises(RuntimeError, match="Unknown section Parameter"):
        check_docstring_parameters(Klass.f_bad_sections)

    incorrect = check_docstring_parameters(f_check_param_definition)
    mock_meta = MockMetaEstimator(delegate=MockEst())
    mock_meta_name = mock_meta.__class__.__name__
    assert incorrect == [
        (
            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('a: int')"
        ),
        (
            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('b:')"
        ),
        (
            "sklearn.utils.tests.test_testing.f_check_param_definition There "
            "was no space between the param name and colon ('d:int')"
        ),
    ]

    messages = [
        [
            "In function: sklearn.utils.tests.test_testing.f_bad_order",
            (
                "There's a parameter name mismatch in function docstring w.r.t."
                " function signature, at index 0 diff: 'b' != 'a'"
            ),
            "Full diff:",
            "- ['b', 'a']",
            "+ ['a', 'b']",
        ],
        [
            "In function: sklearn.utils.tests.test_testing.f_too_many_param_docstring",
            (
                "Parameters in function docstring have more items w.r.t. function"
                " signature, first extra item: c"
            ),
            "Full diff:",
            "- ['a', 'b']",
            "+ ['a', 'b', 'c']",
            "?          +++++",
        ],
        [
            "In function: sklearn.utils.tests.test_testing.f_missing",
            (
                "Parameters in function docstring have less items w.r.t. function"
                " signature, first missing item: b"
            ),
            "Full diff:",
            "- ['a', 'b']",
            "+ ['a']",
        ],
        [
            "In function: sklearn.utils.tests.test_testing.Klass.f_missing",
            (
                "Parameters in function docstring have less items w.r.t. function"
                " signature, first missing item: X"
            ),
            "Full diff:",
            "- ['X', 'y']",
            "+ []",
        ],
        [
            f"In function: sklearn.utils.tests.test_testing.{mock_meta_name}.predict",
            (
                "There's a parameter name mismatch in function docstring w.r.t."
                " function signature, at index 0 diff: 'X' != 'y'"
            ),
            "Full diff:",
            "- ['X']",
            "?   ^",
            "+ ['y']",
            "?   ^",
        ],
        [
            "In function: "
            f"sklearn.utils.tests.test_testing.{mock_meta_name}."
            "predict_proba",
            "potentially wrong underline length... ",
            "Parameters ",
            "--------- in ",
        ],
        [
            f"In function: sklearn.utils.tests.test_testing.{mock_meta_name}.score",
            "potentially wrong underline length... ",
            "Parameters ",
            "--------- in ",
        ],
        [
            f"In function: sklearn.utils.tests.test_testing.{mock_meta_name}.fit",
            (
                "Parameters in function docstring have less items w.r.t. function"
                " signature, first missing item: X"
            ),
            "Full diff:",
            "- ['X', 'y']",
            "+ []",
        ],
    ]

    for msg, f in zip(
        messages,
        [
            f_bad_order,
            f_too_many_param_docstring,
            f_missing,
            Klass.f_missing,
            mock_meta.predict,
            mock_meta.predict_proba,
            mock_meta.score,
            mock_meta.fit,
        ],
    ):
        incorrect = check_docstring_parameters(f)
        assert msg == incorrect, '\n"%s"\n not in \n"%s"' % (msg, incorrect)


def f_one(a, b):  # pragma: no cover
    """Function one.

    Parameters
    ----------
    a : int,   float
        Parameter a.
        Second    line.

    b : str
        Parameter b.

    Returns
    -------
    c : int
       Returning

    d : int
       Returning
    """
    pass


def f_two(a, b):  # pragma: no cover
    """Function two.

    Parameters
    ----------
    a :   int, float
        Parameter a.
          Second line.

    b : str
        Parameter bb.

    e : int
        Extra parameter.

    Returns
    -------
    c : int
       Returning

    d : int
       Returning
    """
    pass


def f_three(a, b):  # pragma: no cover
    """Function two.

    Parameters
    ----------
    a :   int, float
        Parameter a.

    b : str
        Parameter B!

    e :
        Extra parameter.

    Returns
    -------
    c : int
       Returning.

    d : int
       Returning
    """
    pass


@skip_if_no_numpydoc
def test_assert_docstring_consistency_object_type():
    """Check error raised when `objects` incorrect type."""
    with pytest.raises(TypeError, match="All 'objects' must be one of"):
        assert_docstring_consistency(["string", f_one])


@skip_if_no_numpydoc
@pytest.mark.parametrize(
    "objects, kwargs, error",
    [
        (
            [f_one, f_two],
            {"include_params": ["a"], "exclude_params": ["b"]},
            "The 'exclude_params' argument",
        ),
        (
            [f_one, f_two],
            {"include_returns": False, "exclude_returns": ["c"]},
            "The 'exclude_returns' argument",
        ),
    ],
)
def test_assert_docstring_consistency_arg_checks(objects, kwargs, error):
    """Check `assert_docstring_consistency` argument checking correct."""
    with pytest.raises(TypeError, match=error):
        assert_docstring_consistency(objects, **kwargs)


@skip_if_no_numpydoc
@pytest.mark.parametrize(
    "objects, kwargs, error, warn",
    [
        pytest.param(
            [f_one, f_two], {"include_params": ["a"]}, "", "", id="whitespace"
        ),
        pytest.param([f_one, f_two], {"include_returns": True}, "", "", id="incl_all"),
        pytest.param(
            [f_one, f_two, f_three],
            {"include_params": ["a"]},
            (
                r"The description of Parameter 'a' is inconsistent between "
                r"\['f_one',\n'f_two'\]"
            ),
            "",
            id="2-1 group",
        ),
        pytest.param(
            [f_one, f_two, f_three],
            {"include_params": ["b"]},
            (
                r"The description of Parameter 'b' is inconsistent between "
                r"\['f_one'\] and\n\['f_two'\] and"
            ),
            "",
            id="1-1-1 group",
        ),
        pytest.param(
            [f_two, f_three],
            {"include_params": ["e"]},
            (
                r"The type specification of Parameter 'e' is inconsistent between\n"
                r"\['f_two'\] and"
            ),
            "",
            id="empty type",
        ),
        pytest.param(
            [f_one, f_two],
            {"include_params": True, "exclude_params": ["b"]},
            "",
            r"Checking was skipped for Parameters: \['e'\]",
            id="skip warn",
        ),
    ],
)
def test_assert_docstring_consistency(objects, kwargs, error, warn):
    """Check `assert_docstring_consistency` gives correct results."""
    if error:
        with pytest.raises(AssertionError, match=error):
            assert_docstring_consistency(objects, **kwargs)
    elif warn:
        with pytest.warns(UserWarning, match=warn):
            assert_docstring_consistency(objects, **kwargs)
    else:
        assert_docstring_consistency(objects, **kwargs)


def f_four(labels):  # pragma: no cover
    """Function four.

    Parameters
    ----------

    labels : array-like, default=None
        The set of labels to include when `average != 'binary'`, and their
        order if `average is None`. Labels present in the data can be excluded.
    """
    pass


def f_five(labels):  # pragma: no cover
    """Function five.

    Parameters
    ----------

    labels : array-like, default=None
        The set of labels to include when `average != 'binary'`, and their
        order if `average is None`. This is an extra line. Labels present in the
        data can be excluded.
    """
    pass


def f_six(labels):  # pragma: no cover
    """Function six.

    Parameters
    ----------

    labels : array-like, default=None
        The group of labels to add when `average != 'binary'`, and the
        order if `average is None`. Labels present on them datas can be excluded.
    """
    pass


@skip_if_no_numpydoc
def test_assert_docstring_consistency_error_msg():
    """Check `assert_docstring_consistency` difference message."""
    msg = r"""The description of Parameter 'labels' is inconsistent between
\['f_four'\] and \['f_five'\] and \['f_six'\]:

\*\*\* \['f_four'\]
--- \['f_five'\]
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*\*\* 10,25 \*\*\*\*

--- 10,30 ----

  'binary'`, and their order if `average is None`.
\+ This is an extra line.
  Labels present in the data can be excluded.

\*\*\* \['f_four'\]
--- \['f_six'\]
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*\*\* 1,25 \*\*\*\*

  The
! set
  of labels to
! include
  when `average != 'binary'`, and
! their
  order if `average is None`. Labels present
! in the data
  can be excluded.
--- 1,25 ----

  The
! group
  of labels to
! add
  when `average != 'binary'`, and
! the
  order if `average is None`. Labels present
! on them datas
  can be excluded."""

    with pytest.raises(AssertionError, match=msg):
        assert_docstring_consistency([f_four, f_five, f_six], include_params=True)


@skip_if_no_numpydoc
def test_assert_docstring_consistency_descr_regex_pattern():
    """Check `assert_docstring_consistency` `descr_regex_pattern` works."""
    # Check regex that matches full parameter descriptions
    regex_full = (
        r"The (set|group) "  # match 'set' or 'group'
        r"of labels to (include|add) "  # match 'include' or 'add'
        r"when `average \!\= 'binary'`, and (their|the) "  #  match 'their' or 'the'
        r"order if `average is None`\."
        r"[\s\w]*\.* "  # optionally match additional sentence
        r"Labels present (on|in) "  # match 'on' or 'in'
        r"(them|the) "  # match 'them' or 'the'
        r"datas? can be excluded\."  # match 'data' or 'datas'
    )

    assert_docstring_consistency(
        [f_four, f_five, f_six],
        include_params=True,
        descr_regex_pattern=" ".join(regex_full.split()),
    )
    # Check we can just match a few alternate words
    regex_words = r"(labels|average|binary)"  # match any of these 3 words
    assert_docstring_consistency(
        [f_four, f_five, f_six],
        include_params=True,
        descr_regex_pattern=" ".join(regex_words.split()),
    )
    # Check error raised when regex doesn't match
    regex_error = r"The set of labels to include when.+"
    msg = r"The description of Parameter 'labels' in \['f_six'\] does not match"
    with pytest.raises(AssertionError, match=msg):
        assert_docstring_consistency(
            [f_four, f_five, f_six],
            include_params=True,
            descr_regex_pattern=" ".join(regex_error.split()),
        )


class RegistrationCounter:
    def __init__(self):
        self.nb_calls = 0

    def __call__(self, to_register_func):
        self.nb_calls += 1
        assert to_register_func.func is _delete_folder


def check_memmap(input_array, mmap_data, mmap_mode="r"):
    assert isinstance(mmap_data, np.memmap)
    writeable = mmap_mode != "r"
    assert mmap_data.flags.writeable is writeable
    np.testing.assert_array_equal(input_array, mmap_data)


def test_tempmemmap(monkeypatch):
    registration_counter = RegistrationCounter()
    monkeypatch.setattr(atexit, "register", registration_counter)

    input_array = np.ones(3)
    with TempMemmap(input_array) as data:
        check_memmap(input_array, data)
        temp_folder = os.path.dirname(data.filename)
    if os.name != "nt":
        assert not os.path.exists(temp_folder)
    assert registration_counter.nb_calls == 1

    mmap_mode = "r+"
    with TempMemmap(input_array, mmap_mode=mmap_mode) as data:
        check_memmap(input_array, data, mmap_mode=mmap_mode)
        temp_folder = os.path.dirname(data.filename)
    if os.name != "nt":
        assert not os.path.exists(temp_folder)
    assert registration_counter.nb_calls == 2


def test_create_memmap_backed_data(monkeypatch):
    registration_counter = RegistrationCounter()
    monkeypatch.setattr(atexit, "register", registration_counter)

    input_array = np.ones(3)
    data = create_memmap_backed_data(input_array)
    check_memmap(input_array, data)
    assert registration_counter.nb_calls == 1

    data, folder = create_memmap_backed_data(input_array, return_folder=True)
    check_memmap(input_array, data)
    assert folder == os.path.dirname(data.filename)
    assert registration_counter.nb_calls == 2

    mmap_mode = "r+"
    data = create_memmap_backed_data(input_array, mmap_mode=mmap_mode)
    check_memmap(input_array, data, mmap_mode)
    assert registration_counter.nb_calls == 3

    input_list = [input_array, input_array + 1, input_array + 2]
    mmap_data_list = create_memmap_backed_data(input_list)
    for input_array, data in zip(input_list, mmap_data_list):
        check_memmap(input_array, data)
    assert registration_counter.nb_calls == 4

    output_data, other = create_memmap_backed_data([input_array, "not-an-array"])
    check_memmap(input_array, output_data)
    assert other == "not-an-array"


@pytest.mark.parametrize(
    "constructor_kwargs, expected_output_type",
    [
        ({"constructor_type": "list"}, list),
        ({"constructor_type": "tuple"}, tuple),
        ({"constructor_type": "array"}, np.ndarray),
        # Use `zip` to keep only the available sparse containers depending of the
        # installed scipy version
        *zip(
            (
                {
                    "constructor_type": "array",
                    "sparse_container": sparse_container,
                    "sparse_format": "csr",
                }
                for sparse_container in ("matrix", "array")
            ),
            CSR_CONTAINERS,
        ),
        *zip(
            (
                {
                    "constructor_type": "array",
                    "sparse_container": sparse_container,
                    "sparse_format": "csc",
                }
                for sparse_container in ("matrix", "array")
            ),
            CSC_CONTAINERS,
        ),
        # Use a lambda function to delay the import of optional dependencies to within
        # the function to only skip a specific test instead of the whole file
        (
            {"constructor_type": "dataframe", "constructor_lib": "pandas"},
            lambda: pytest.importorskip("pandas").DataFrame,
        ),
        (
            {"constructor_type": "dataframe", "constructor_lib": "polars"},
            lambda: pytest.importorskip("polars").DataFrame,
        ),
        (
            {"constructor_type": "dataframe", "constructor_lib": "pyarrow"},
            lambda: pytest.importorskip("pyarrow").Table,
        ),
        (
            {"constructor_type": "series", "constructor_lib": "pandas"},
            lambda: pytest.importorskip("pandas").Series,
        ),
        (
            {"constructor_type": "series", "constructor_lib": "polars"},
            lambda: pytest.importorskip("polars").Series,
        ),
        (
            {"constructor_type": "series", "constructor_lib": "pyarrow"},
            lambda: pytest.importorskip("pyarrow").Array,
        ),
        (
            {"constructor_type": "index", "constructor_lib": "pandas"},
            lambda: pytest.importorskip("pandas").Index,
        ),
    ],
)
@pytest.mark.parametrize(
    "input_type, input_shape",
    [
        (list, (-1,)),
        (list, (2, -1)),
        (np.array, (-1,)),
        (np.array, (2, -1)),
        # sparse containers are always 2-dimensional
        *[(csr_container, (2, -1)) for csr_container in CSR_CONTAINERS],
    ],
)
@pytest.mark.parametrize(
    "dtype, superdtype, polars_dtype, pyarrow_dtype",
    [
        (
            np.int32,
            np.integer,
            lambda: pytest.importorskip("polars").Int32,
            lambda: pytest.importorskip("pyarrow").int32(),
        ),
        (
            np.int64,
            np.integer,
            lambda: pytest.importorskip("polars").Int64,
            lambda: pytest.importorskip("pyarrow").int64(),
        ),
        (
            np.float32,
            np.floating,
            lambda: pytest.importorskip("polars").Float32,
            lambda: pytest.importorskip("pyarrow").float32(),
        ),
        (
            np.float64,
            np.floating,
            lambda: pytest.importorskip("polars").Float64,
            lambda: pytest.importorskip("pyarrow").float64(),
        ),
    ],
)
def test_convert_container(
    constructor_kwargs,
    expected_output_type,
    input_type,
    input_shape,
    dtype,
    superdtype,
    polars_dtype,
    pyarrow_dtype,
):
    """
    Check that we convert the container to the right type of array with the
    right data type and shape.
    """
    constructor_type = constructor_kwargs["constructor_type"]
    if constructor_type in ("dataframe", "series", "index"):
        # Run the lambda function to import the library or skip
        expected_output_type = expected_output_type()

    container = np.arange(12).reshape(input_shape)

    target_dtype = dtype
    constructor_lib = constructor_kwargs.get("constructor_lib")
    if constructor_lib == "polars":
        target_dtype = polars_dtype()
    elif constructor_lib == "pyarrow":
        target_dtype = pyarrow_dtype()

    # Certain constructor types require specific-dimensional input containers
    # None stands for no constraint
    desired_dimension = None
    if constructor_type in ("series", "index"):
        desired_dimension = 1
    elif constructor_type == "dataframe":
        desired_dimension = 2

    if desired_dimension is None or desired_dimension == len(input_shape):
        container_converted = _convert_container(
            input_type(container.tolist()),
            dtype=target_dtype,
            **constructor_kwargs,
        )
    else:  # dimension mismatch
        msg = (
            f"Only {desired_dimension}D containers can be converted to "
            f"{constructor_type}; got shape {container.shape} instead"
        )
        with pytest.raises(ValueError, match=re.escape(msg)):
            _convert_container(
                input_type(container.tolist()),
                dtype=target_dtype,
                **constructor_kwargs,
            )
        return

    # Check output data type and shape
    assert isinstance(container_converted, expected_output_type)
    if constructor_type in ("list", "tuple"):
        # list and tuple do not have a shape attribute
        if len(input_shape) == 1:
            assert len(container_converted) == container.size
        else:
            assert len(container_converted) == container.shape[0]
            assert len(container_converted[0]) == container.shape[1]
    elif constructor_type == "series" and constructor_lib == "pyarrow":
        # pyarrow.Array does not have a shape attribute
        assert len(container_converted) == container.size
    elif (
        constructor_type == "dataframe"
        or constructor_kwargs.get("sparse_container") is not None
    ) and len(input_shape) == 1:
        # dataframe and sparse container will convert to 2-dimensional
        assert container_converted.shape == (1, container.size)
    else:
        assert container_converted.shape == container.shape

    if constructor_type in ("list", "tuple", "index"):
        # list and tuple will use Python class dtype: int, float
        # pandas index will always use high precision: np.int64 and np.float64
        first_element = (
            container_converted[0]
            if len(input_shape) == 1
            else container_converted[0][0]
        )
        assert np.issubdtype(type(first_element), superdtype)
    else:
        # Determine the dtype of the converted container
        if constructor_type == "dataframe" and constructor_lib == "pyarrow":
            converted_dtype = container_converted.schema.types[0]
        elif constructor_type == "series" and constructor_lib == "pyarrow":
            converted_dtype = container_converted.type
        elif hasattr(container_converted, "dtype"):
            converted_dtype = container_converted.dtype
        elif hasattr(container_converted, "dtypes"):
            converted_dtype = container_converted.dtypes[0]
        else:
            assert False, (
                f"{type(container_converted).__name__} has no dtype"
            )  # pragma: no cover
        assert converted_dtype == target_dtype, (
            f"Expected {target_dtype}, got {converted_dtype}"
        )


def test_convert_container_slice():
    """Check that converting to slice works and raises correctly."""
    container = [0, 1]
    result = _convert_container(container, "slice")
    assert result == slice(*container, None)

    container = [0, 1, 2]
    msg = (
        "Only 1D containers with exactly 2 elements can be converted to slice; got "
        "shape (3,) instead"
    )
    with pytest.raises(ValueError, match=re.escape(msg)):
        _convert_container(container, "slice")


def test_convert_container_invalid_sparse_options():
    """Check that invalid sparse options raise the correct error."""
    container = np.arange(12).reshape(3, -1)

    msg = "Invalid sparse_format='dok'; expected one of ('csr', 'csc')"
    with pytest.raises(ValueError, match=re.escape(msg)):
        _convert_container(
            container,
            constructor_type="array",
            sparse_container="array",
            sparse_format="dok",
        )

    msg = "Invalid sparse_container='series'; expected one of ('matrix', 'array')"
    with pytest.raises(ValueError, match=re.escape(msg)):
        _convert_container(
            container, constructor_type="array", sparse_container="series"
        )


@pytest.mark.parametrize("constructor_type", ["dataframe", "series", "index"])
def test_convert_container_invalid_library(constructor_type):
    """Check that incompatible library raises the correct error."""
    if constructor_type == "dataframe":
        container = [[0, 1], [2, 3]]
    else:
        container = [0, 1]

    msg = f"constructor_lib='scipy' is incompatible with {constructor_type}"
    with pytest.raises(ValueError, match=re.escape(msg)):
        _convert_container(
            container, constructor_type=constructor_type, constructor_lib="scipy"
        )


@pytest.mark.parametrize("constructor_lib", ["pandas", "polars", "pyarrow"])
def test_convert_container_categories(constructor_lib):
    lib = pytest.importorskip(constructor_lib)
    df = _convert_container(
        [["x"]],
        constructor_type="dataframe",
        constructor_lib=constructor_lib,
        column_names=["A"],
        categorical_feature_names=["A"],
    )

    if constructor_lib == "pandas":
        assert isinstance(df, lib.DataFrame)
        assert df.dtypes.iloc[0] == "category"
    elif constructor_lib == "polars":
        assert isinstance(df, lib.DataFrame)
        assert df.schema["A"] == lib.Categorical()
    elif constructor_lib == "pyarrow":
        assert isinstance(df, lib.Table)
        assert type(df.schema[0].type) is lib.DictionaryType


@pytest.mark.parametrize("sparse_container", ["matrix", "array"])
@pytest.mark.parametrize("sparse_format", ["csc", "csr"])
@pytest.mark.parametrize("orig_csr_container", CSR_CONTAINERS)
def test_convert_container_sparse_to_sparse(
    sparse_container, sparse_format, orig_csr_container
):
    """
    Non-regression test to check that we can still convert a sparse container
    from a given format to another format.
    """
    X_sparse = orig_csr_container(np.random.randn(10, 10))
    X_converted = _convert_container(
        X_sparse,
        constructor_type="array",
        sparse_container=sparse_container,
        sparse_format=sparse_format,
    )
    assert X_converted.format == sparse_format


@pytest.mark.parametrize("constructor_lib_from", ["pandas", "polars", "pyarrow"])
@pytest.mark.parametrize("constructor_lib_to", ["pandas", "polars", "pyarrow"])
@pytest.mark.parametrize("column_names", [None, ["C", "D"]])
@pytest.mark.parametrize(
    "dtype, polars_dtype, pyarrow_dtype",
    [
        (
            np.int32,
            lambda: pytest.importorskip("polars").Int32,
            lambda: pytest.importorskip("pyarrow").int32(),
        ),
        (
            np.float64,
            lambda: pytest.importorskip("polars").Float64,
            lambda: pytest.importorskip("pyarrow").float64(),
        ),
    ],
)
def test_convert_container_df_to_df(
    constructor_lib_from,
    constructor_lib_to,
    column_names,
    dtype,
    polars_dtype,
    pyarrow_dtype,
):
    """Check that we can convert a DataFrame to another DataFrame."""
    pytest.importorskip(constructor_lib_from)
    lib_to = pytest.importorskip(constructor_lib_to)

    column_names_from = ["A", "B"]
    df_from = _convert_container(
        [[1, 2], [3, 4]],
        "dataframe",
        constructor_lib=constructor_lib_from,
        column_names=column_names_from,
    )

    target_dtype = dtype
    if constructor_lib_to == "polars":
        target_dtype = polars_dtype()
    elif constructor_lib_to == "pyarrow":
        target_dtype = pyarrow_dtype()

    df_to = _convert_container(
        df_from,
        "dataframe",
        column_names=column_names,
        constructor_lib=constructor_lib_to,
        dtype=target_dtype,
    )
    target_column_names = column_names_from if column_names is None else column_names
    if constructor_lib_to == "pandas":
        assert isinstance(df_to, lib_to.DataFrame)
        assert df_to.dtypes.iloc[0] == target_dtype
        assert df_to.columns.to_list() == target_column_names
    elif constructor_lib_to == "polars":
        assert isinstance(df_to, lib_to.DataFrame)
        assert df_to.dtypes[0] == target_dtype
        assert df_to.columns == target_column_names
    elif constructor_lib_to == "pyarrow":
        assert isinstance(df_to, lib_to.Table)
        assert df_to.schema[0].type == target_dtype
        assert df_to.column_names == target_column_names


@pytest.mark.parametrize("constructor_lib_from", ["pandas", "polars"])
@pytest.mark.parametrize("constructor_lib_to", ["pandas", "polars"])
@pytest.mark.parametrize(
    "dtype, polars_dtype",
    [
        (np.int32, lambda: pytest.importorskip("polars").Int32),
        (np.float64, lambda: pytest.importorskip("polars").Float64),
    ],
)
def test_convert_container_series_to_series(
    constructor_lib_from, constructor_lib_to, dtype, polars_dtype
):
    """Check that we can convert a Series to another Series."""
    lib_to = pytest.importorskip(constructor_lib_to)

    ser_from = _convert_container(
        [1, 2, 3, 4],
        "series",
        constructor_lib=constructor_lib_from,
    )

    target_dtype = dtype
    if constructor_lib_to == "polars":
        target_dtype = polars_dtype()

    ser_to = _convert_container(
        ser_from,
        "series",
        constructor_lib=constructor_lib_to,
        dtype=target_dtype,
    )
    assert isinstance(ser_to, lib_to.Series)
    assert ser_to.dtype == target_dtype


def test_raises():
    # Tests for the raises context manager

    # Proper type, no match
    with raises(TypeError):
        raise TypeError()

    # Proper type, proper match
    with raises(TypeError, match="how are you") as cm:
        raise TypeError("hello how are you")
    assert cm.raised_and_matched

    # Proper type, proper match with multiple patterns
    with raises(TypeError, match=["not this one", "how are you"]) as cm:
        raise TypeError("hello how are you")
    assert cm.raised_and_matched

    # bad type, no match
    with pytest.raises(ValueError, match="this will be raised"):
        with raises(TypeError) as cm:
            raise ValueError("this will be raised")
    assert not cm.raised_and_matched

    # Bad type, no match, with a err_msg
    with pytest.raises(AssertionError, match="the failure message"):
        with raises(TypeError, err_msg="the failure message") as cm:
            raise ValueError()
    assert not cm.raised_and_matched

    # bad type, with match (is ignored anyway)
    with pytest.raises(ValueError, match="this will be raised"):
        with raises(TypeError, match="this is ignored") as cm:
            raise ValueError("this will be raised")
    assert not cm.raised_and_matched

    # proper type but bad match
    with pytest.raises(
        AssertionError, match="should contain one of the following patterns"
    ):
        with raises(TypeError, match="hello") as cm:
            raise TypeError("Bad message")
    assert not cm.raised_and_matched

    # proper type but bad match, with err_msg
    with pytest.raises(AssertionError, match="the failure message"):
        with raises(TypeError, match="hello", err_msg="the failure message") as cm:
            raise TypeError("Bad message")
    assert not cm.raised_and_matched

    # no raise with default may_pass=False
    with pytest.raises(AssertionError, match="Did not raise"):
        with raises(TypeError) as cm:
            pass
    assert not cm.raised_and_matched

    # no raise with may_pass=True
    with raises(TypeError, match="hello", may_pass=True) as cm:
        pass  # still OK
    assert not cm.raised_and_matched

    # Multiple exception types:
    with raises((TypeError, ValueError)):
        raise TypeError()
    with raises((TypeError, ValueError)):
        raise ValueError()
    with pytest.raises(AssertionError):
        with raises((TypeError, ValueError)):
            pass


def test_float32_aware_assert_allclose():
    # The relative tolerance for float32 inputs is 1e-4
    assert_allclose(np.array([1.0 + 2e-5], dtype=np.float32), 1.0)
    with pytest.raises(AssertionError):
        assert_allclose(np.array([1.0 + 2e-4], dtype=np.float32), 1.0)

    # The relative tolerance for other inputs is left to 1e-7 as in
    # the original numpy version.
    assert_allclose(np.array([1.0 + 2e-8], dtype=np.float64), 1.0)
    with pytest.raises(AssertionError):
        assert_allclose(np.array([1.0 + 2e-7], dtype=np.float64), 1.0)

    # atol is left to 0.0 by default, even for float32
    with pytest.raises(AssertionError):
        assert_allclose(np.array([1e-5], dtype=np.float32), 0.0)
    assert_allclose(np.array([1e-5], dtype=np.float32), 0.0, atol=2e-5)


@pytest.mark.xfail(_IS_WASM, reason="cannot start subprocess")
def test_assert_run_python_script_without_output():
    code = "x = 1"
    assert_run_python_script_without_output(code)

    code = "print('something to stdout')"
    with pytest.raises(AssertionError, match="Expected no output"):
        assert_run_python_script_without_output(code)

    code = "print('something to stdout')"
    with pytest.raises(
        AssertionError,
        match="output was not supposed to match.+got.+something to stdout",
    ):
        assert_run_python_script_without_output(code, pattern="to.+stdout")

    code = "\n".join(["import sys", "print('something to stderr', file=sys.stderr)"])
    with pytest.raises(
        AssertionError,
        match="output was not supposed to match.+got.+something to stderr",
    ):
        assert_run_python_script_without_output(code, pattern="to.+stderr")


def check_warnings_as_errors(warning_info, warnings_as_errors):
    if warning_info.action == "error" and warnings_as_errors:
        with pytest.raises(warning_info.category, match=warning_info.message):
            warnings.warn(
                message=warning_info.message,
                category=warning_info.category,
            )
    if warning_info.action == "ignore":
        with warnings.catch_warnings(record=True) as record:
            message = warning_info.message
            # Special treatment when regex is used
            if "Pyarrow" in message:
                message = "\nPyarrow will become a required dependency"

            warnings.warn(
                message=message,
                category=warning_info.category,
            )
            assert len(record) == 0 if warnings_as_errors else 1
            if record:
                assert str(record[0].message) == message
                assert record[0].category == warning_info.category


@pytest.mark.parametrize("warning_info", _get_warnings_filters_info_list())
def test_sklearn_warnings_as_errors(warning_info):
    warnings_as_errors = os.environ.get("SKLEARN_WARNINGS_AS_ERRORS", "0") != "0"
    check_warnings_as_errors(warning_info, warnings_as_errors=warnings_as_errors)


@pytest.mark.parametrize("warning_info", _get_warnings_filters_info_list())
def test_turn_warnings_into_errors(warning_info):
    with warnings.catch_warnings():
        turn_warnings_into_errors()
        check_warnings_as_errors(warning_info, warnings_as_errors=True)
