from __future__ import annotations

import contextlib
from typing import cast

from polars._utils.deprecation import deprecate_renamed_parameter
from polars.dataframe import DataFrame
from polars.lazyframe import LazyFrame
from polars.testing.asserts.utils import raise_assertion_error

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import assert_dataframe_equal_py


def _assert_correct_input_type(
    left: DataFrame | LazyFrame, right: DataFrame | LazyFrame
) -> bool:
    __tracebackhide__ = True

    if isinstance(left, DataFrame) and isinstance(right, DataFrame):
        return False
    elif isinstance(left, LazyFrame) and isinstance(right, LazyFrame):
        return True
    else:
        raise_assertion_error(
            "inputs",
            "unexpected input types",
            type(left).__name__,
            type(right).__name__,
        )


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
@deprecate_renamed_parameter("rtol", "rel_tol", version="1.32.3")
@deprecate_renamed_parameter("atol", "abs_tol", version="1.32.3")
def assert_frame_equal(
    left: DataFrame | LazyFrame,
    right: DataFrame | LazyFrame,
    *,
    check_row_order: bool = True,
    check_column_order: bool = True,
    check_dtypes: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-5,
    abs_tol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right frame are equal.

    Raises a detailed `AssertionError` if the frames differ.
    This function is intended for use in unit tests.

    .. versionchanged:: 0.20.31
        The `check_dtype` parameter was renamed `check_dtypes`.

    .. versionchanged:: 1.32.3
        The `rtol` and `atol` parameters were renamed to `rel_tol` and `abs_tol`,
        respectively.

    Parameters
    ----------
    left
        The first DataFrame or LazyFrame to compare.
    right
        The second DataFrame or LazyFrame to compare.
    check_row_order
        Requires row order to match.
    check_column_order
        Requires column order to match.
    check_dtypes
        Requires data types to match.
    check_exact
        Requires float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rel_tol` and `abs_tol`).
        Only affects columns with a Float data type.
    rel_tol
        Relative tolerance for inexact checking. Fraction of values in `right`.
    abs_tol
        Absolute tolerance for inexact checking.
    categorical_as_str
        Cast categorical columns to string before comparing. Enabling this helps
        compare columns that do not share the same string cache.

    See Also
    --------
    assert_series_equal
    assert_frame_not_equal

    Notes
    -----
    When using pytest, it may be worthwhile to shorten Python traceback printing
    by passing `--tb=short`. The default mode tends to be unhelpfully verbose.
    More information in the
    `pytest docs <https://docs.pytest.org/en/latest/how-to/output.html#modifying-python-traceback-printing>`_.

    Examples
    --------
    >>> from polars.testing import assert_frame_equal
    >>> df1 = pl.DataFrame({"a": [1, 2, 3]})
    >>> df2 = pl.DataFrame({"a": [1, 5, 3]})
    >>> assert_frame_equal(df1, df2)
    Traceback (most recent call last):
    ...
    AssertionError: DataFrames are different (value mismatch for column "a")
    [left]: shape: (3,)
    Series: 'a' [i64]
    [
        1
        2
        3
    ]
    [right]: shape: (3,)
    Series: 'a' [i64]
    [
        1
        5
        3
    ]
    """
    __tracebackhide__ = True

    lazy = _assert_correct_input_type(left, right)

    # Rust back-end function expects DataFrames so LazyFrames must be collected
    if lazy:
        left, right = left.collect(), right.collect()  # type: ignore[union-attr]

    # Tell type checker these are now DataFrames to prevent type errors
    left, right = cast("DataFrame", left), cast("DataFrame", right)

    assert_dataframe_equal_py(
        left._df,
        right._df,
        check_row_order=check_row_order,
        check_column_order=check_column_order,
        check_dtypes=check_dtypes,
        check_exact=check_exact,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        categorical_as_str=categorical_as_str,
    )


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
@deprecate_renamed_parameter("rtol", "rel_tol", version="1.32.3")
@deprecate_renamed_parameter("atol", "abs_tol", version="1.32.3")
def assert_frame_not_equal(
    left: DataFrame | LazyFrame,
    right: DataFrame | LazyFrame,
    *,
    check_row_order: bool = True,
    check_column_order: bool = True,
    check_dtypes: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-5,
    abs_tol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right frame are **not** equal.

    This function is intended for use in unit tests.

    .. versionchanged:: 0.20.31
        The `check_dtype` parameter was renamed `check_dtypes`.

    .. versionchanged:: 1.32.3
        The `rtol` and `atol` parameters were renamed to `rel_tol` and `abs_tol`,
        respectively.

    Parameters
    ----------
    left
        The first DataFrame or LazyFrame to compare.
    right
        The second DataFrame or LazyFrame to compare.
    check_row_order
        Requires row order to match.
    check_column_order
        Requires column order to match.
    check_dtypes
        Requires data types to match.
    check_exact
        Requires float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rel_tol` and `abs_tol`).
        Only affects columns with a Float data type.
    rel_tol
        Relative tolerance for inexact checking. Fraction of values in `right`.
    abs_tol
        Absolute tolerance for inexact checking.
    categorical_as_str
        Cast categorical columns to string before comparing. Enabling this helps
        compare columns that do not share the same string cache.

    See Also
    --------
    assert_frame_equal
    assert_series_not_equal

    Examples
    --------
    >>> from polars.testing import assert_frame_not_equal
    >>> df1 = pl.DataFrame({"a": [1, 2, 3]})
    >>> df2 = pl.DataFrame({"a": [1, 2, 3]})
    >>> assert_frame_not_equal(df1, df2)
    Traceback (most recent call last):
    ...
    AssertionError: DataFrames are equal (but are expected not to be)
    """
    __tracebackhide__ = True

    lazy = _assert_correct_input_type(left, right)
    try:
        assert_frame_equal(
            left=left,
            right=right,
            check_column_order=check_column_order,
            check_row_order=check_row_order,
            check_dtypes=check_dtypes,
            check_exact=check_exact,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            categorical_as_str=categorical_as_str,
        )
    except AssertionError:
        return
    else:
        objects = "LazyFrames" if lazy else "DataFrames"
        msg = f"{objects} are equal (but are expected not to be)"
        raise AssertionError(msg)
