from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, NoReturn, overload

import polars._reexport as pl
import polars.functions as F
from polars._dependencies import _check_for_numpy
from polars._dependencies import numpy as np
from polars._utils.constants import U32_MAX
from polars._utils.slice import PolarsSlice
from polars._utils.various import qualified_type_name, range_to_slice
from polars.datatypes.classes import (
    Boolean,
    Int8,
    Int16,
    Int32,
    Int64,
    String,
    UInt32,
    UInt64,
)
from polars.meta.index_type import get_index_type

if TYPE_CHECKING:
    from collections.abc import Iterable

    from polars import DataFrame, Series
    from polars._typing import (
        MultiColSelector,
        MultiIndexSelector,
        SingleColSelector,
        SingleIndexSelector,
    )

__all__ = [
    "get_df_item_by_key",
    "get_series_item_by_key",
]


@overload
def get_series_item_by_key(s: Series, key: SingleIndexSelector) -> Any: ...


@overload
def get_series_item_by_key(s: Series, key: MultiIndexSelector) -> Series: ...


def get_series_item_by_key(
    s: Series, key: SingleIndexSelector | MultiIndexSelector
) -> Any | Series:
    """Select one or more elements from the Series."""
    if isinstance(key, int):
        return s._s.get_index_signed(key)

    elif isinstance(key, slice):
        return _select_elements_by_slice(s, key)

    elif isinstance(key, range):
        key = range_to_slice(key)
        return _select_elements_by_slice(s, key)

    elif isinstance(key, Sequence):
        if not key:
            return s.clear()

        first = key[0]
        if isinstance(first, bool):
            _raise_on_boolean_mask()

        try:
            indices = pl.Series("", key, dtype=Int64)
        except TypeError:
            msg = f"cannot select elements using Sequence with elements of type {qualified_type_name(first)!r}"
            raise TypeError(msg) from None

        indices = _convert_series_to_indices(indices, s.len())
        return _select_elements_by_index(s, indices)

    elif isinstance(key, pl.Series):
        indices = _convert_series_to_indices(key, s.len())
        return _select_elements_by_index(s, indices)

    elif _check_for_numpy(key) and isinstance(key, np.ndarray):
        indices = _convert_np_ndarray_to_indices(key, s.len())
        return _select_elements_by_index(s, indices)

    msg = f"cannot select elements using key of type {qualified_type_name(key)!r}: {key!r}"
    raise TypeError(msg)


def _select_elements_by_slice(s: Series, key: slice) -> Series:
    return PolarsSlice(s).apply(key)  # type: ignore[return-value]


def _select_elements_by_index(s: Series, key: Series) -> Series:
    return s._from_pyseries(s._s.gather_with_series(key._s))


# `str` overlaps with `Sequence[str]`
# We can ignore this but we must keep this overload ordering
@overload
def get_df_item_by_key(
    df: DataFrame, key: tuple[SingleIndexSelector, SingleColSelector]
) -> Any: ...


@overload
def get_df_item_by_key(  # type: ignore[overload-overlap]
    df: DataFrame, key: str | tuple[MultiIndexSelector, SingleColSelector]
) -> Series: ...


@overload
def get_df_item_by_key(
    df: DataFrame,
    key: (
        SingleIndexSelector
        | MultiIndexSelector
        | MultiColSelector
        | tuple[SingleIndexSelector, MultiColSelector]
        | tuple[MultiIndexSelector, MultiColSelector]
    ),
) -> DataFrame: ...


def get_df_item_by_key(
    df: DataFrame,
    key: (
        SingleIndexSelector
        | SingleColSelector
        | MultiColSelector
        | MultiIndexSelector
        | tuple[SingleIndexSelector, SingleColSelector]
        | tuple[SingleIndexSelector, MultiColSelector]
        | tuple[MultiIndexSelector, SingleColSelector]
        | tuple[MultiIndexSelector, MultiColSelector]
    ),
) -> DataFrame | Series | Any:
    """Get part of the DataFrame as a new DataFrame, Series, or scalar."""
    # Two inputs, e.g. df[1, 2:5]
    if isinstance(key, tuple) and len(key) == 2:
        row_key, col_key = key

        # Support df[True, False] and df["a", "b"] as these are not ambiguous
        if isinstance(row_key, (bool, str)):
            return _select_columns(df, key)  # type: ignore[arg-type]

        selection = _select_columns(df, col_key)

        if selection.is_empty():
            return selection
        elif isinstance(selection, pl.Series):
            return get_series_item_by_key(selection, row_key)
        else:
            return _select_rows(selection, row_key)

    # Single string input, e.g. df["a"]
    if isinstance(key, str):
        # This case is required because empty strings are otherwise treated
        # as an empty Sequence in `_select_rows`
        return df.get_column(key)

    # Single input - df[1] - or multiple inputs - df["a", "b", "c"]
    try:
        return _select_rows(df, key)  # type: ignore[arg-type]
    except TypeError:
        return _select_columns(df, key)


# `str` overlaps with `Sequence[str]`
# We can ignore this but we must keep this overload ordering
@overload
def _select_columns(df: DataFrame, key: SingleColSelector) -> Series: ...  # type: ignore[overload-overlap]


@overload
def _select_columns(df: DataFrame, key: MultiColSelector) -> DataFrame: ...


def _select_columns(
    df: DataFrame, key: SingleColSelector | MultiColSelector
) -> DataFrame | Series:
    """Select one or more columns from the DataFrame."""
    if isinstance(key, int):
        return df.to_series(key)

    elif isinstance(key, str):
        return df.get_column(key)

    elif isinstance(key, slice):
        start, stop, step = key.start, key.stop, key.step
        # Fast path for common case: df[x, :]
        if start is None and stop is None and step is None:
            return df
        if isinstance(start, str):
            start = df.get_column_index(start)
        if isinstance(stop, str):
            stop = df.get_column_index(stop) + 1
        int_slice = slice(start, stop, step)
        rng = range(df.width)[int_slice]
        return _select_columns_by_index(df, rng)

    elif isinstance(key, range):
        return _select_columns_by_index(df, key)

    elif isinstance(key, Sequence):
        if not key:
            return df.__class__()
        first = key[0]
        if isinstance(first, bool):
            return _select_columns_by_mask(df, key)  # type: ignore[arg-type]
        elif isinstance(first, int):
            return _select_columns_by_index(df, key)  # type: ignore[arg-type]
        elif isinstance(first, str):
            return _select_columns_by_name(df, key)  # type: ignore[arg-type]
        else:
            msg = f"cannot select columns using Sequence with elements of type {qualified_type_name(first)!r}"
            raise TypeError(msg)

    elif isinstance(key, pl.Series):
        if key.is_empty():
            return df.__class__()
        dtype = key.dtype
        if dtype == String:
            return _select_columns_by_name(df, key)
        elif dtype.is_integer():
            return _select_columns_by_index(df, key)
        elif dtype == Boolean:
            return _select_columns_by_mask(df, key)
        else:
            msg = f"cannot select columns using Series of type {dtype}"
            raise TypeError(msg)

    elif _check_for_numpy(key) and isinstance(key, np.ndarray):
        if key.ndim == 0:
            key = np.atleast_1d(key)
        elif key.ndim != 1:
            msg = "multi-dimensional NumPy arrays not supported as index"
            raise TypeError(msg)

        if len(key) == 0:
            return df.__class__()

        dtype_kind = key.dtype.kind
        if dtype_kind in ("i", "u"):
            return _select_columns_by_index(df, key)
        elif dtype_kind == "b":
            return _select_columns_by_mask(df, key)
        elif isinstance(key[0], str):
            return _select_columns_by_name(df, key)
        else:
            msg = f"cannot select columns using NumPy array of type {key.dtype}"
            raise TypeError(msg)

    msg = (
        f"cannot select columns using key of type {qualified_type_name(key)!r}: {key!r}"
    )
    raise TypeError(msg)


def _select_columns_by_index(df: DataFrame, key: Iterable[int]) -> DataFrame:
    series = [df.to_series(i) for i in key]
    return df.__class__(series)


def _select_columns_by_name(df: DataFrame, key: Iterable[str]) -> DataFrame:
    return df._from_pydf(df._df.select(list(key)))


def _select_columns_by_mask(
    df: DataFrame, key: Sequence[bool] | Series | np.ndarray[Any, Any]
) -> DataFrame:
    if len(key) != df.width:
        msg = f"expected {df.width} values when selecting columns by boolean mask, got {len(key)}"
        raise ValueError(msg)

    indices = (i for i, val in enumerate(key) if val)
    return _select_columns_by_index(df, indices)


@overload
def _select_rows(df: DataFrame, key: SingleIndexSelector) -> Series: ...


@overload
def _select_rows(df: DataFrame, key: MultiIndexSelector) -> DataFrame: ...


def _select_rows(
    df: DataFrame, key: SingleIndexSelector | MultiIndexSelector
) -> DataFrame | Series:
    """Select one or more rows from the DataFrame."""
    if isinstance(key, int):
        num_rows = df.height
        if (key >= num_rows) or (key < -num_rows):
            msg = f"index {key} is out of bounds for DataFrame of height {num_rows}"
            raise IndexError(msg)
        return df.slice(key, 1)

    if isinstance(key, slice):
        return _select_rows_by_slice(df, key)

    elif isinstance(key, range):
        key = range_to_slice(key)
        return _select_rows_by_slice(df, key)

    elif isinstance(key, Sequence):
        if not key:
            return df.clear()
        if isinstance(key[0], bool):
            _raise_on_boolean_mask()
        s = pl.Series("", key, dtype=Int64)
        indices = _convert_series_to_indices(s, df.height)
        return _select_rows_by_index(df, indices)

    elif isinstance(key, pl.Series):
        indices = _convert_series_to_indices(key, df.height)
        return _select_rows_by_index(df, indices)

    elif _check_for_numpy(key) and isinstance(key, np.ndarray):
        indices = _convert_np_ndarray_to_indices(key, df.height)
        return _select_rows_by_index(df, indices)

    else:
        msg = f"cannot select rows using key of type {qualified_type_name(key)!r}: {key!r}"
        raise TypeError(msg)


def _select_rows_by_slice(df: DataFrame, key: slice) -> DataFrame:
    return PolarsSlice(df).apply(key)  # type: ignore[return-value]


def _select_rows_by_index(df: DataFrame, key: Series) -> DataFrame:
    return df._from_pydf(df._df.gather_with_series(key._s))


# UTILS


def _convert_series_to_indices(s: Series, size: int) -> Series:
    """Convert a Series to indices, taking into account negative values."""
    # Unsigned or signed Series (ordered from fastest to slowest).
    #   - pl.UInt32 (polars) or pl.UInt64 (polars_u64_idx) Series indexes.
    #   - Other unsigned Series indexes are converted to pl.UInt32 (polars)
    #     or pl.UInt64 (polars_u64_idx).
    #   - Signed Series indexes are converted pl.UInt32 (polars) or
    #     pl.UInt64 (polars_u64_idx) after negative indexes are converted
    #     to absolute indexes.

    # pl.UInt32 (polars) or pl.UInt64 (polars_u64_idx).
    idx_type = get_index_type()

    if s.dtype == idx_type:
        return s

    if not s.dtype.is_integer():
        if s.dtype == Boolean:
            _raise_on_boolean_mask()
        else:
            msg = f"cannot treat Series of type {s.dtype} as indices"
            raise TypeError(msg)

    if s.len() == 0:
        return pl.Series(s.name, [], dtype=idx_type)

    if idx_type == UInt32:
        if s.dtype in {Int64, UInt64} and s.max() >= U32_MAX:  # type: ignore[operator]
            msg = "index positions should be smaller than 2^32"
            raise ValueError(msg)
        if s.dtype == Int64 and s.min() < -U32_MAX:  # type: ignore[operator]
            msg = "index positions should be greater than or equal to -2^32"
            raise ValueError(msg)

    if s.dtype.is_signed_integer():
        if s.min() < 0:  # type: ignore[operator]
            if idx_type == UInt32:
                idxs = s.cast(Int32) if s.dtype in {Int8, Int16} else s
            else:
                idxs = s.cast(Int64) if s.dtype in {Int8, Int16, Int32} else s

            # Update negative indexes to absolute indexes.
            return (
                idxs.to_frame()
                .select(
                    F.when(F.col(idxs.name) < 0)
                    .then(size + F.col(idxs.name))
                    .otherwise(F.col(idxs.name))
                    .cast(idx_type)
                )
                .to_series(0)
            )

    return s.cast(idx_type)


def _convert_np_ndarray_to_indices(arr: np.ndarray[Any, Any], size: int) -> Series:
    """Convert a NumPy ndarray to indices, taking into account negative values."""
    # Unsigned or signed Numpy array (ordered from fastest to slowest).
    #   - np.uint32 (polars) or np.uint64 (polars_u64_idx) numpy array
    #     indexes.
    #   - Other unsigned numpy array indexes are converted to pl.UInt32
    #     (polars) or pl.UInt64 (polars_u64_idx).
    #   - Signed numpy array indexes are converted pl.UInt32 (polars) or
    #     pl.UInt64 (polars_u64_idx) after negative indexes are converted
    #     to absolute indexes.
    if arr.ndim == 0:
        arr = np.atleast_1d(arr)
    if arr.ndim != 1:
        msg = "only 1D NumPy arrays can be treated as indices"
        raise TypeError(msg)

    idx_type = get_index_type()

    if len(arr) == 0:
        return pl.Series("", [], dtype=idx_type)

    # Numpy array with signed or unsigned integers.
    if arr.dtype.kind not in ("i", "u"):
        if arr.dtype.kind == "b":
            _raise_on_boolean_mask()
        else:
            msg = f"cannot treat NumPy array of type {arr.dtype} as indices"
            raise TypeError(msg)

    if idx_type == UInt32:
        if arr.dtype in {np.int64, np.uint64} and arr.max() >= U32_MAX:
            msg = "index positions should be smaller than 2^32"
            raise ValueError(msg)
        if arr.dtype == np.int64 and arr.min() < -U32_MAX:
            msg = "index positions should be greater than or equal to -2^32"
            raise ValueError(msg)

    if arr.dtype.kind == "i" and arr.min() < 0:
        if idx_type == UInt32:
            if arr.dtype in (np.int8, np.int16):
                arr = arr.astype(np.int32)
        else:
            if arr.dtype in (np.int8, np.int16, np.int32):
                arr = arr.astype(np.int64)

        # Update negative indexes to absolute indexes.
        arr = np.where(arr < 0, size + arr, arr)

    # numpy conversion is much faster
    arr = arr.astype(np.uint32) if idx_type == UInt32 else arr.astype(np.uint64)

    return pl.Series("", arr, dtype=idx_type)


def _raise_on_boolean_mask() -> NoReturn:
    msg = (
        "selecting rows by passing a boolean mask to `__getitem__` is not supported"
        "\n\nHint: Use the `filter` method instead."
    )
    raise TypeError(msg)
