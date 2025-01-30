from __future__ import annotations

from typing import TYPE_CHECKING, Any

from polars._utils.deprecation import deprecate_renamed_parameter
from polars.datatypes import (
    Array,
    Categorical,
    List,
    String,
    Struct,
    unpack_dtypes,
)
from polars.datatypes.group import FLOAT_DTYPES
from polars.exceptions import ComputeError, InvalidOperationError, ShapeError
from polars.series import Series
from polars.testing.asserts.utils import raise_assertion_error

if TYPE_CHECKING:
    from polars import DataType


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
def assert_series_equal(
    left: Series,
    right: Series,
    *,
    check_dtypes: bool = True,
    check_names: bool = True,
    check_order: bool = True,
    check_exact: bool = False,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right Series are equal.

    Raises a detailed `AssertionError` if the Series differ.
    This function is intended for use in unit tests.

    Parameters
    ----------
    left
        The first Series to compare.
    right
        The second Series to compare.
    check_dtypes
        Require data types to match.
    check_names
        Require names to match.
    check_order
        Require elements to appear in the same order.
    check_exact
        Require float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rtol` and `atol`).
        Only affects columns with a Float data type.
    rtol
        Relative tolerance for inexact checking, given as a fraction of the values in
        `right`.
    atol
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
    [left]:  [1, 2, 3]
    [right]: [1, 5, 3]
    """
    __tracebackhide__ = True

    _assert_correct_input_type(left, right)

    if left.len() != right.len():
        raise_assertion_error("Series", "length mismatch", left.len(), right.len())

    if check_names and left.name != right.name:
        raise_assertion_error("Series", "name mismatch", left.name, right.name)

    if check_dtypes and left.dtype != right.dtype:
        raise_assertion_error("Series", "dtype mismatch", left.dtype, right.dtype)

    _assert_series_values_equal(
        left,
        right,
        check_order=check_order,
        check_exact=check_exact,
        rtol=rtol,
        atol=atol,
        categorical_as_str=categorical_as_str,
    )


def _assert_series_values_equal(
    left: Series,
    right: Series,
    *,
    check_order: bool,
    check_exact: bool,
    rtol: float,
    atol: float,
    categorical_as_str: bool,
) -> None:
    """Assert that the values in both Series are equal."""
    __tracebackhide__ = True

    # Handle categoricals
    if categorical_as_str:
        left = _categorical_series_to_string(left)
        right = _categorical_series_to_string(right)

    if not check_order:
        left, right = _sort_series(left, right)

    # Determine unequal elements
    try:
        unequal = left.ne_missing(right)
    except ComputeError as exc:
        raise_assertion_error(
            "Series",
            "incompatible data types",
            left=left.dtype,
            right=right.dtype,
            cause=exc,
        )
    except ShapeError as exc:
        raise_assertion_error(
            "Series",
            "incompatible lengths",
            left=left,
            right=right,
            cause=exc,
        )

    # Check nested dtypes in separate function
    if _comparing_nested_floats(left.dtype, right.dtype):
        try:
            _assert_series_nested_values_equal(
                left=left.filter(unequal),
                right=right.filter(unequal),
                check_exact=check_exact,
                rtol=rtol,
                atol=atol,
                categorical_as_str=categorical_as_str,
            )
        except AssertionError as exc:
            raise_assertion_error(
                "Series",
                "nested value mismatch",
                left=left.to_list(),
                right=right.to_list(),
                cause=exc,
            )
        else:  # All nested values match
            return

    # If no differences found during exact checking, we're done
    if not unequal.any():
        return

    # Only do inexact checking for float types
    if check_exact or not left.dtype.is_float() or not right.dtype.is_float():
        raise_assertion_error(
            "Series", "exact value mismatch", left=left.to_list(), right=right.to_list()
        )

    _assert_series_null_values_match(left, right)
    _assert_series_nan_values_match(left, right)
    _assert_series_values_within_tolerance(
        left,
        right,
        unequal,
        rtol=rtol,
        atol=atol,
    )


def _sort_series(left: Series, right: Series) -> tuple[Series, Series]:
    try:
        left = left.sort()
        right = right.sort()
    except InvalidOperationError as exc:
        msg = "cannot set `check_order=False` on Series with unsortable data type"
        raise TypeError(msg) from exc
    return left, right


def _assert_series_nested_values_equal(
    left: Series,
    right: Series,
    *,
    check_exact: bool,
    rtol: float,
    atol: float,
    categorical_as_str: bool,
) -> None:
    __tracebackhide__ = True

    # compare nested lists element-wise
    if _comparing_lists(left.dtype, right.dtype):
        for s1, s2 in zip(left, right):
            if s1 is None or s2 is None:
                raise_assertion_error("Series", "nested value mismatch", s1, s2)

            _assert_series_values_equal(
                s1,
                s2,
                check_order=True,
                check_exact=check_exact,
                rtol=rtol,
                atol=atol,
                categorical_as_str=categorical_as_str,
            )

    # unnest structs as series and compare
    else:
        ls, rs = left.struct.unnest(), right.struct.unnest()
        for s1, s2 in zip(ls, rs):
            _assert_series_values_equal(
                s1,
                s2,
                check_order=True,
                check_exact=check_exact,
                rtol=rtol,
                atol=atol,
                categorical_as_str=categorical_as_str,
            )


def _assert_series_null_values_match(left: Series, right: Series) -> None:
    __tracebackhide__ = True
    null_value_mismatch = left.is_null() != right.is_null()
    if null_value_mismatch.any():
        raise_assertion_error(
            "Series", "null value mismatch", left.to_list(), right.to_list()
        )


def _assert_series_nan_values_match(left: Series, right: Series) -> None:
    __tracebackhide__ = True
    if not _comparing_floats(left.dtype, right.dtype):
        return
    nan_value_mismatch = left.is_nan() != right.is_nan()
    if nan_value_mismatch.any():
        raise_assertion_error(
            "Series",
            "nan value mismatch",
            left.to_list(),
            right.to_list(),
        )


def _comparing_floats(left: DataType, right: DataType) -> bool:
    return left.is_float() and right.is_float()


def _comparing_lists(left: DataType, right: DataType) -> bool:
    return left in (List, Array) and right in (List, Array)


def _comparing_structs(left: DataType, right: DataType) -> bool:
    return left == Struct and right == Struct


def _comparing_nested_floats(left: DataType, right: DataType) -> bool:
    if not (_comparing_lists(left, right) or _comparing_structs(left, right)):
        return False

    return bool(FLOAT_DTYPES & unpack_dtypes(left)) and bool(
        FLOAT_DTYPES & unpack_dtypes(right)
    )


def _assert_series_values_within_tolerance(
    left: Series,
    right: Series,
    unequal: Series,
    *,
    rtol: float,
    atol: float,
) -> None:
    __tracebackhide__ = True

    left_unequal, right_unequal = left.filter(unequal), right.filter(unequal)

    difference = (left_unequal - right_unequal).abs()
    tolerance = atol + rtol * right_unequal.abs()
    exceeds_tolerance = difference > tolerance

    if exceeds_tolerance.any():
        raise_assertion_error(
            "Series",
            "value mismatch",
            left.to_list(),
            right.to_list(),
        )


def _categorical_series_to_string(s: Series) -> Series:
    """Cast a (possibly nested) Categorical Series to a String Series."""
    dtype = s.dtype
    noncat_dtype = _categorical_dtype_to_string_dtype(dtype)
    if dtype != noncat_dtype:
        s = s.cast(noncat_dtype)
    return s


def _categorical_dtype_to_string_dtype(dtype: DataType) -> DataType:
    """Change a (possibly nested) Categorical data type to a String data type."""
    if isinstance(dtype, Categorical):
        return String()
    elif isinstance(dtype, List):
        inner_cast = _categorical_dtype_to_string_dtype(dtype.inner)  # type: ignore[arg-type]
        return List(inner_cast)
    elif isinstance(dtype, Array):
        inner_cast = _categorical_dtype_to_string_dtype(dtype.inner)  # type: ignore[arg-type]
        return Array(inner_cast, dtype.size)
    elif isinstance(dtype, Struct):
        fields = {
            f.name: _categorical_dtype_to_string_dtype(f.dtype)  # type: ignore[arg-type]
            for f in dtype.fields
        }
        return Struct(fields)
    else:
        return dtype


@deprecate_renamed_parameter("check_dtype", "check_dtypes", version="0.20.31")
def assert_series_not_equal(
    left: Series,
    right: Series,
    *,
    check_dtypes: bool = True,
    check_names: bool = True,
    check_order: bool = True,
    check_exact: bool = False,
    rtol: float = 1e-5,
    atol: float = 1e-8,
    categorical_as_str: bool = False,
) -> None:
    """
    Assert that the left and right Series are **not** equal.

    This function is intended for use in unit tests.

    Parameters
    ----------
    left
        The first Series to compare.
    right
        The second Series to compare.
    check_dtypes
        Require data types to match.
    check_names
        Require names to match.
    check_order
        Require elements to appear in the same order.
    check_exact
        Require float values to match exactly. If set to `False`, values are considered
        equal when within tolerance of each other (see `rtol` and `atol`).
        Only affects columns with a Float data type.
    rtol
        Relative tolerance for inexact checking, given as a fraction of the values in
        `right`.
    atol
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
            rtol=rtol,
            atol=atol,
            categorical_as_str=categorical_as_str,
        )
    except AssertionError:
        return
    else:
        msg = "Series are equal (but are expected not to be)"
        raise AssertionError(msg)
