from __future__ import annotations

import contextlib
from typing import Any

from polars._utils.deprecation import deprecate_renamed_parameter
from polars.series import Series
from polars.testing.asserts.utils import raise_assertion_error

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import assert_series_equal_py


def _assert_correct_input_type(left: Any, right: Any) -> bool:
    __tracebackhide__ = True

    if not (isinstance(left, Series) and isinstance(right, Series)):
        raise_assertion_error(
            "inputs",
            "unexpected input types",
            type(left).__name__,
            type(right).__name__,
        )
    return True


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
@deprecate_renamed_parameter("rtol", "rel_tol", version="1.32.3")
@deprecate_renamed_parameter("atol", "abs_tol", version="1.32.3")
def assert_series_equal(
    left: Series,
    right: Series,
    *,
    check_dtypes: bool = True,
    check_names: bool = True,
    check_order: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-5,
    abs_tol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right Series are equal.

    Raises a detailed `AssertionError` if the Series differ.
    This function is intended for use in unit tests.

    .. versionchanged:: 0.20.31
        The `check_dtype` parameter was renamed `check_dtypes`.

    .. versionchanged:: 1.32.3
        The `rtol` and `atol` parameters were renamed to `rel_tol` and `abs_tol`,
        respectively.

    Parameters
    ----------
    left
        The first Series to compare.
    right
        The second Series to compare.
    check_dtypes
        Requires data types to match.
    check_names
        Requires names to match.
    check_order
        Requires elements to appear in the same order.
    check_exact
        Requires float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rel_tol` and `abs_tol`).
        Only affects columns with a Float data type.
    rel_tol
        Relative tolerance for inexact checking, given as a fraction of the values in
        `right`.
    abs_tol
        Absolute tolerance for inexact checking.
    categorical_as_str
        Cast categorical columns to string before comparing. Enabling this helps
        compare columns that do not share the same string cache.

    See Also
    --------
    assert_frame_equal
    assert_series_not_equal

    Notes
    -----
    When using pytest, it may be worthwhile to shorten Python traceback printing
    by passing `--tb=short`. The default mode tends to be unhelpfully verbose.
    More information in the
    `pytest docs <https://docs.pytest.org/en/latest/how-to/output.html#modifying-python-traceback-printing>`_.

    Examples
    --------
    >>> from polars.testing import assert_series_equal
    >>> s1 = pl.Series([1, 2, 3])
    >>> s2 = pl.Series([1, 5, 3])
    >>> assert_series_equal(s1, s2)
    Traceback (most recent call last):
    ...
    AssertionError: Series are different (exact value mismatch)
    [left]: shape: (3,)
    Series: '' [i64]
    [
        1
        2
        3
    ]
    [right]: shape: (3,)
    Series: '' [i64]
    [
        1
        5
        3
    ]
    """
    __tracebackhide__ = True

    _assert_correct_input_type(left, right)

    assert_series_equal_py(
        left._s,
        right._s,
        check_dtypes=check_dtypes,
        check_names=check_names,
        check_order=check_order,
        check_exact=check_exact,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        categorical_as_str=categorical_as_str,
    )


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
@deprecate_renamed_parameter("rtol", "rel_tol", version="1.32.3")
@deprecate_renamed_parameter("atol", "abs_tol", version="1.32.3")
def assert_series_not_equal(
    left: Series,
    right: Series,
    *,
    check_dtypes: bool = True,
    check_names: bool = True,
    check_order: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-5,
    abs_tol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right Series are **not** equal.

    This function is intended for use in unit tests.

    .. versionchanged:: 0.20.31
        The `check_dtype` parameter was renamed `check_dtypes`.

    .. versionchanged:: 1.32.3
        The `rtol` and `atol` parameters were renamed to `rel_tol` and `abs_tol`,
        respectively.

    Parameters
    ----------
    left
        The first Series to compare.
    right
        The second Series to compare.
    check_dtypes
        Requires data types to match.
    check_names
        Requires names to match.
    check_order
        Requires elements to appear in the same order.
    check_exact
        Requires float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rel_tol` and `abs_tol`).
        Only affects columns with a Float data type.
    rel_tol
        Relative tolerance for inexact checking, given as a fraction of the values in
        `right`.
    abs_tol
        Absolute tolerance for inexact checking.
    categorical_as_str
        Cast categorical columns to string before comparing. Enabling this helps
        compare columns that do not share the same string cache.

    See Also
    --------
    assert_series_equal
    assert_frame_not_equal

    Examples
    --------
    >>> from polars.testing import assert_series_not_equal
    >>> s1 = pl.Series([1, 2, 3])
    >>> s2 = pl.Series([1, 2, 3])
    >>> assert_series_not_equal(s1, s2)
    Traceback (most recent call last):
    ...
    AssertionError: Series are equal (but are expected not to be)
    """
    __tracebackhide__ = True

    _assert_correct_input_type(left, right)
    try:
        assert_series_equal(
            left=left,
            right=right,
            check_dtypes=check_dtypes,
            check_names=check_names,
            check_order=check_order,
            check_exact=check_exact,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            categorical_as_str=categorical_as_str,
        )
    except AssertionError:
        return
    else:
        msg = "Series are equal (but are expected not to be)"
        raise AssertionError(msg)
