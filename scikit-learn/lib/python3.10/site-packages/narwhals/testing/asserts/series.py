from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable

from narwhals._utils import qualified_type_name, zip_strict
from narwhals.dependencies import is_narwhals_series
from narwhals.dtypes import Array, Boolean, Categorical, List, String, Struct
from narwhals.functions import new_series
from narwhals.testing.asserts.utils import raise_series_assertion_error

if TYPE_CHECKING:
    from typing_extensions import TypeAlias

    from narwhals.series import Series
    from narwhals.typing import IntoSeriesT, SeriesT

    CheckFn: TypeAlias = Callable[[Series[Any], Series[Any]], None]


def assert_series_equal(
    left: Series[IntoSeriesT],
    right: Series[IntoSeriesT],
    *,
    check_dtypes: bool = True,
    check_names: bool = True,
    check_order: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-05,
    abs_tol: float = 1e-08,
    categorical_as_str: bool = False,
) -> None:
    """Assert that the left and right Series are equal.

    Raises a detailed `AssertionError` if the Series differ.
    This function is intended for use in unit tests.

    Arguments:
        left: The first Series to compare.
        right: The second Series to compare.
        check_dtypes: Requires data types to match.
        check_names: Requires names to match.
        check_order: Requires elements to appear in the same order.
        check_exact: Requires float values to match exactly. If set to `False`, values are
            considered equal when within tolerance of each other (see `rel_tol` and
            `abs_tol`). Only affects columns with a Float data type.
        rel_tol: Relative tolerance for inexact checking, given as a fraction of the
            values in `right`.
        abs_tol: Absolute tolerance for inexact checking.
        categorical_as_str: Cast categorical columns to string before comparing.
            Enabling this helps compare columns that do not share the same string cache.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>> from narwhals.testing import assert_series_equal
        >>> s1 = nw.from_native(pd.Series([1, 2, 3]), series_only=True)
        >>> s2 = nw.from_native(pd.Series([1, 5, 3]), series_only=True)
        >>> assert_series_equal(s1, s2)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        AssertionError: Series are different (exact value mismatch)
        [left]:
        ┌───────────────┐
        |Narwhals Series|
        |---------------|
        | 0    1        |
        | 1    2        |
        | 2    3        |
        | dtype: int64  |
        └───────────────┘
        [right]:
        ┌───────────────┐
        |Narwhals Series|
        |---------------|
        | 0    1        |
        | 1    5        |
        | 2    3        |
        | dtype: int64  |
        └───────────────┘
    """
    __tracebackhide__ = True

    if any(not is_narwhals_series(obj) for obj in (left, right)):
        msg = (
            "Expected `narwhals.Series` instance, found:\n"
            f"[left]: {qualified_type_name(type(left))}\n"
            f"[right]: {qualified_type_name(type(right))}\n\n"
            "Hint: Use `nw.from_native(obj, series_only=True) to convert each native "
            "object into a `narwhals.Series` first."
        )
        raise TypeError(msg)

    _check_metadata(left, right, check_dtypes=check_dtypes, check_names=check_names)

    if not check_order:
        if left.dtype.is_nested():
            msg = "`check_order=False` is not supported (yet) with nested data type."
            raise NotImplementedError(msg)
        left, right = left.sort(), right.sort()

    left_vals, right_vals = _check_null_values(left, right)

    if check_exact or not left.dtype.is_float():
        _check_exact_values(
            left_vals,
            right_vals,
            check_dtypes=check_dtypes,
            check_exact=check_exact,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            categorical_as_str=categorical_as_str,
        )
    else:
        _check_approximate_values(left_vals, right_vals, rel_tol=rel_tol, abs_tol=abs_tol)


def _check_metadata(
    left: SeriesT, right: SeriesT, *, check_dtypes: bool, check_names: bool
) -> None:
    """Check metadata information: implementation, length, dtype, and names."""
    left_impl, right_impl = left.implementation, right.implementation
    if left_impl != right_impl:
        raise_series_assertion_error("implementation mismatch", left_impl, right_impl)

    left_len, right_len = len(left), len(right)
    if left_len != right_len:
        raise_series_assertion_error("length mismatch", left_len, right_len)

    left_dtype, right_dtype = left.dtype, right.dtype
    if check_dtypes and left_dtype != right_dtype:
        raise_series_assertion_error("dtype mismatch", left_dtype, right_dtype)

    left_name, right_name = left.name, right.name
    if check_names and left_name != right_name:
        raise_series_assertion_error("name mismatch", left_name, right_name)


def _check_null_values(left: SeriesT, right: SeriesT) -> tuple[SeriesT, SeriesT]:
    """Check null value consistency and return non-null values."""
    left_null_count, right_null_count = left.null_count(), right.null_count()
    left_null_mask, right_null_mask = left.is_null(), right.is_null()

    if left_null_count != right_null_count or (left_null_mask != right_null_mask).any():
        raise_series_assertion_error(
            "null value mismatch", left_null_count, right_null_count
        )

    return left.filter(~left_null_mask), right.filter(~right_null_mask)


def _check_exact_values(
    left: SeriesT,
    right: SeriesT,
    *,
    check_dtypes: bool,
    check_exact: bool,
    rel_tol: float,
    abs_tol: float,
    categorical_as_str: bool,
) -> None:
    """Check exact value equality for various data types."""
    left_impl = left.implementation
    left_dtype, right_dtype = left.dtype, right.dtype

    is_not_equal_mask: Series[Any]
    if left_dtype.is_numeric():
        # For _all_ numeric dtypes, we can use `is_close` with 0-tolerances to handle
        # inf and nan values out of the box.
        is_not_equal_mask = ~left.is_close(right, rel_tol=0, abs_tol=0, nans_equal=True)
    elif (
        isinstance(left_dtype, (Array, List)) and isinstance(right_dtype, (Array, List))
    ) and left_dtype == right_dtype:
        check_fn = partial(
            assert_series_equal,
            check_dtypes=check_dtypes,
            check_names=False,
            check_order=True,
            check_exact=check_exact,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            categorical_as_str=categorical_as_str,
        )
        _check_list_like(left, right, left_dtype, right_dtype, check_fn=check_fn)
        # If `_check_list_like` didn't raise, then every nested element is equal
        is_not_equal_mask = new_series("", [False], dtype=Boolean(), backend=left_impl)
    elif isinstance(left_dtype, Struct) and isinstance(right_dtype, Struct):
        check_fn = partial(
            assert_series_equal,
            check_dtypes=True,
            check_names=True,
            check_order=True,
            check_exact=check_exact,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            categorical_as_str=categorical_as_str,
        )
        _check_struct(left, right, left_dtype, right_dtype, check_fn=check_fn)
        # If `_check_struct` didn't raise, then every nested element is equal
        is_not_equal_mask = new_series("", [False], dtype=Boolean(), backend=left_impl)
    elif isinstance(left_dtype, Categorical) and isinstance(right_dtype, Categorical):
        # If `_check_categorical` didn't raise, then the categories sources/encodings are
        # the same, and we can use equality
        _not_equal = _check_categorical(
            left, right, categorical_as_str=categorical_as_str
        )
        is_not_equal_mask = new_series(
            "", [_not_equal], dtype=Boolean(), backend=left_impl
        )
    else:
        is_not_equal_mask = left != right

    if is_not_equal_mask.any():
        raise_series_assertion_error("exact value mismatch", left, right)


def _check_approximate_values(
    left: SeriesT, right: SeriesT, *, rel_tol: float, abs_tol: float
) -> None:
    """Check approximate value equality with tolerance."""
    is_not_close_mask = ~left.is_close(
        right, rel_tol=rel_tol, abs_tol=abs_tol, nans_equal=True
    )

    if is_not_close_mask.any():
        raise_series_assertion_error(
            "values not within tolerance",
            left.filter(is_not_close_mask),
            right.filter(is_not_close_mask),
        )


def _check_list_like(
    left_vals: SeriesT,
    right_vals: SeriesT,
    left_dtype: List | Array,
    right_dtype: List | Array,
    check_fn: CheckFn,
) -> None:
    # Check row by row after transforming each array/list into a new series.
    # Notice that order within the array/list must be the same, regardless of
    # `check_order` value at the top level.
    impl = left_vals.implementation
    try:
        for left_val, right_val in zip_strict(left_vals, right_vals):
            check_fn(
                new_series("", values=left_val, dtype=left_dtype.inner, backend=impl),
                new_series("", values=right_val, dtype=right_dtype.inner, backend=impl),
            )
    except AssertionError:
        raise_series_assertion_error("nested value mismatch", left_vals, right_vals)


def _check_struct(
    left_vals: SeriesT,
    right_vals: SeriesT,
    left_dtype: Struct,
    right_dtype: Struct,
    check_fn: CheckFn,
) -> None:
    # Check field by field as a separate column.
    # Notice that for struct's polars raises if:
    #   * field names are different but values are equal
    #   * dtype differs, regardless of `check_dtypes=False`
    #   * order applies only at top level
    try:
        for left_field, right_field in zip_strict(left_dtype.fields, right_dtype.fields):
            check_fn(
                left_vals.struct.field(left_field.name),
                right_vals.struct.field(right_field.name),
            )
    except AssertionError:
        raise_series_assertion_error("exact value mismatch", left_vals, right_vals)


def _check_categorical(
    left_vals: SeriesT, right_vals: SeriesT, *, categorical_as_str: bool
) -> bool:
    """Try to compare if any element of categorical series' differ.

    Inability to compare means that the encoding is different, and an exception is raised.
    """
    if categorical_as_str:
        left_vals, right_vals = left_vals.cast(String()), right_vals.cast(String())

    try:
        return (left_vals != right_vals).any()
    except Exception as exc:
        msg = "Cannot compare categoricals coming from different sources."
        # TODO(FBruzzesi): Improve error message?
        raise AssertionError(msg) from exc
