from __future__ import annotations

from typing import cast

from polars._utils.deprecation import deprecate_renamed_parameter
from polars.dataframe import DataFrame
from polars.exceptions import InvalidOperationError
from polars.lazyframe import LazyFrame
from polars.testing.asserts.series import _assert_series_values_equal
from polars.testing.asserts.utils import raise_assertion_error


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
def assert_frame_equal(
    left: DataFrame | LazyFrame,
    right: DataFrame | LazyFrame,
    *,
    check_row_order: bool = True,
    check_column_order: bool = True,
    check_dtypes: bool = True,
    check_exact: bool = False,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right frame are equal.

    Raises a detailed `AssertionError` if the frames differ.
    This function is intended for use in unit tests.

    Parameters
    ----------
    left
        The first DataFrame or LazyFrame to compare.
    right
        The second DataFrame or LazyFrame to compare.
    check_row_order
        Require row order to match.
    check_column_order
        Require column order to match.
    check_dtypes
        Require data types to match.
    check_exact
        Require float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rtol` and `atol`).
        Only affects columns with a Float data type.
    rtol
        Relative tolerance for inexact checking. Fraction of values in `right`.
    atol
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
    AssertionError: DataFrames are different (value mismatch for column 'a')
    [left]:  [1, 2, 3]
    [right]: [1, 5, 3]
    """
    __tracebackhide__ = True

    lazy = _assert_correct_input_type(left, right)
    objects = "LazyFrames" if lazy else "DataFrames"

    _assert_frame_schema_equal(
        left,
        right,
        check_column_order=check_column_order,
        check_dtypes=check_dtypes,
        objects=objects,
    )

    if lazy:
        left, right = left.collect(), right.collect()  # type: ignore[union-attr]
    left, right = cast(DataFrame, left), cast(DataFrame, right)

    if left.height != right.height:
        raise_assertion_error(
            objects, "number of rows does not match", left.height, right.height
        )

    if not check_row_order:
        left, right = _sort_dataframes(left, right)

    for c in left.columns:
        s_left, s_right = left.get_column(c), right.get_column(c)
        try:
            _assert_series_values_equal(
                s_left,
                s_right,
                check_order=True,
                check_exact=check_exact,
                rtol=rtol,
                atol=atol,
                categorical_as_str=categorical_as_str,
            )
        except AssertionError as exc:
            raise_assertion_error(
                objects,
                f"value mismatch for column {c!r}",
                s_left.to_list(),
                s_right.to_list(),
                cause=exc,
            )


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


def _assert_frame_schema_equal(
    left: DataFrame | LazyFrame,
    right: DataFrame | LazyFrame,
    *,
    check_dtypes: bool,
    check_column_order: bool,
    objects: str,
) -> None:
    __tracebackhide__ = True

    left_schema, right_schema = left.collect_schema(), right.collect_schema()

    # Fast path for equal frames
    if left_schema == right_schema:
        return

    # Special error message for when column names do not match
    if left_schema.keys() != right_schema.keys():
        if left_not_right := [c for c in left_schema if c not in right_schema]:
            msg = f"columns {left_not_right!r} in left {objects[:-1]}, but not in right"
            raise AssertionError(msg)
        else:
            right_not_left = [c for c in right_schema if c not in left_schema]
            msg = f"columns {right_not_left!r} in right {objects[:-1]}, but not in left"
            raise AssertionError(msg)

    if check_column_order:
        left_columns, right_columns = list(left_schema), list(right_schema)
        if left_columns != right_columns:
            detail = "columns are not in the same order"
            raise_assertion_error(objects, detail, left_columns, right_columns)

    if check_dtypes:
        left_schema_dict, right_schema_dict = dict(left_schema), dict(right_schema)
        if check_column_order or left_schema_dict != right_schema_dict:
            print(left_schema_dict, right_schema_dict)
            detail = "dtypes do not match"
            raise_assertion_error(objects, detail, left_schema_dict, right_schema_dict)


def _sort_dataframes(left: DataFrame, right: DataFrame) -> tuple[DataFrame, DataFrame]:
    by = left.columns
    try:
        left = left.sort(by)
        right = right.sort(by)
    except InvalidOperationError as exc:
        msg = "cannot set `check_row_order=False` on frame with unsortable columns"
        raise TypeError(msg) from exc
    return left, right


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
def assert_frame_not_equal(
    left: DataFrame | LazyFrame,
    right: DataFrame | LazyFrame,
    *,
    check_row_order: bool = True,
    check_column_order: bool = True,
    check_dtypes: bool = True,
    check_exact: bool = False,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right frame are **not** equal.

    This function is intended for use in unit tests.

    Parameters
    ----------
    left
        The first DataFrame or LazyFrame to compare.
    right
        The second DataFrame or LazyFrame to compare.
    check_row_order
        Require row order to match.
    check_column_order
        Require column order to match.
    check_dtypes
        Require data types to match.
    check_exact
        Require float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rtol` and `atol`).
        Only affects columns with a Float data type.
    rtol
        Relative tolerance for inexact checking. Fraction of values in `right`.
    atol
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
            rtol=rtol,
            atol=atol,
            categorical_as_str=categorical_as_str,
        )
    except AssertionError:
        return
    else:
        objects = "LazyFrames" if lazy else "DataFrames"
        msg = f"{objects} are equal (but are expected not to be)"
        raise AssertionError(msg)
