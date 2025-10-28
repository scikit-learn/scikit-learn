from __future__ import annotations

import functools
from decimal import Decimal as PyDecimal
from typing import TYPE_CHECKING, Any, Callable

from polars import datatypes as dt
from polars._dependencies import numpy as np

# Module not available when building docs
try:
    from polars._plr import PySeries

    _DOCUMENTING = False
except ImportError:
    _DOCUMENTING = True

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars._typing import PolarsDataType

if not _DOCUMENTING:
    _POLARS_TYPE_TO_CONSTRUCTOR: dict[
        PolarsDataType, Callable[[str, Sequence[Any], bool], PySeries]
    ] = {
        dt.Float32: PySeries.new_opt_f32,
        dt.Float64: PySeries.new_opt_f64,
        dt.Int8: PySeries.new_opt_i8,
        dt.Int16: PySeries.new_opt_i16,
        dt.Int32: PySeries.new_opt_i32,
        dt.Int64: PySeries.new_opt_i64,
        dt.Int128: PySeries.new_opt_i128,
        dt.UInt8: PySeries.new_opt_u8,
        dt.UInt16: PySeries.new_opt_u16,
        dt.UInt32: PySeries.new_opt_u32,
        dt.UInt64: PySeries.new_opt_u64,
        dt.UInt128: PySeries.new_opt_u128,
        dt.Decimal: PySeries.new_decimal,
        dt.Date: PySeries.new_opt_i32,
        dt.Datetime: PySeries.new_opt_i64,
        dt.Duration: PySeries.new_opt_i64,
        dt.Time: PySeries.new_opt_i64,
        dt.Boolean: PySeries.new_opt_bool,
        dt.String: PySeries.new_str,
        dt.Object: PySeries.new_object,
        dt.Categorical: PySeries.new_str,
        dt.Enum: PySeries.new_str,
        dt.Binary: PySeries.new_binary,
        dt.Null: PySeries.new_null,
    }


def polars_type_to_constructor(
    dtype: PolarsDataType,
) -> Callable[[str, Sequence[Any], bool], PySeries]:
    """Get the right PySeries constructor for the given Polars dtype."""
    # Special case for Array as it needs to pass the dtype argument on construction
    if isinstance(dtype, dt.Array):
        return functools.partial(PySeries.new_array, dtype=dtype)

    try:
        base_type = dtype.base_type()
        return _POLARS_TYPE_TO_CONSTRUCTOR[base_type]
    except KeyError:  # pragma: no cover
        msg = f"cannot construct PySeries for type {dtype!r}"
        raise ValueError(msg) from None


_NUMPY_TYPE_TO_CONSTRUCTOR = None


def _set_numpy_to_constructor() -> None:
    global _NUMPY_TYPE_TO_CONSTRUCTOR
    _NUMPY_TYPE_TO_CONSTRUCTOR = {
        np.float32: PySeries.new_f32,
        np.float64: PySeries.new_f64,
        np.int8: PySeries.new_i8,
        np.int16: PySeries.new_i16,
        np.int32: PySeries.new_i32,
        np.int64: PySeries.new_i64,
        np.uint8: PySeries.new_u8,
        np.uint16: PySeries.new_u16,
        np.uint32: PySeries.new_u32,
        np.uint64: PySeries.new_u64,
        np.str_: PySeries.new_str,
        np.bytes_: PySeries.new_binary,
        np.bool_: PySeries.new_bool,
        np.datetime64: PySeries.new_i64,
        np.timedelta64: PySeries.new_i64,
    }


@functools.lru_cache(maxsize=32)
def _normalise_numpy_dtype(dtype: Any) -> tuple[Any, Any]:
    normalised_dtype = (
        np.dtype(dtype.base.name) if dtype.kind in ("i", "u", "f") else dtype
    ).type
    cast_as: Any = None
    if normalised_dtype == np.float16:
        normalised_dtype = cast_as = np.float32
    elif normalised_dtype in (np.datetime64, np.timedelta64):
        time_unit = np.datetime_data(dtype)[0]
        if time_unit in dt.DTYPE_TEMPORAL_UNITS or (
            time_unit == "D" and normalised_dtype == np.datetime64
        ):
            cast_as = np.int64
        else:
            msg = (
                "incorrect NumPy datetime resolution"
                "\n\n'D' (datetime only), 'ms', 'us', and 'ns' resolutions are supported when converting from numpy.{datetime64,timedelta64}."
                " Please cast to the closest supported unit before converting."
            )
            raise ValueError(msg)
    return normalised_dtype, cast_as


def numpy_values_and_dtype(
    values: np.ndarray[Any, Any],
) -> tuple[np.ndarray[Any, Any], type]:
    """Return numpy values and their associated dtype, adjusting if required."""
    # Create new dtype object from dtype base name so architecture specific
    # dtypes (np.longlong np.ulonglong np.intc np.uintc np.longdouble, ...)
    # get converted to their normalized dtype (np.int*, np.uint*, np.float*).
    dtype, cast_as = _normalise_numpy_dtype(values.dtype)
    if cast_as:
        values = values.astype(cast_as)
    return values, dtype


def numpy_type_to_constructor(
    values: np.ndarray[Any, Any], dtype: type[np.dtype[Any]]
) -> Callable[..., PySeries]:
    """Get the right PySeries constructor for the given Polars dtype."""
    if _NUMPY_TYPE_TO_CONSTRUCTOR is None:
        _set_numpy_to_constructor()
    try:
        return _NUMPY_TYPE_TO_CONSTRUCTOR[dtype]  # type:ignore[index]
    except KeyError:
        if len(values) > 0:
            first_non_nan = next(
                (v for v in values if isinstance(v, np.ndarray) or v == v), None
            )
            if isinstance(first_non_nan, str):
                return PySeries.new_str
            if isinstance(first_non_nan, bytes):
                return PySeries.new_binary
        return PySeries.new_object
    except NameError:  # pragma: no cover
        msg = f"'numpy' is required to convert numpy dtype {dtype!r}"
        raise ModuleNotFoundError(msg) from None


if not _DOCUMENTING:
    _PY_TYPE_TO_CONSTRUCTOR = {
        float: PySeries.new_opt_f64,
        bool: PySeries.new_opt_bool,
        int: PySeries.new_opt_i64,
        str: PySeries.new_str,
        bytes: PySeries.new_binary,
        PyDecimal: PySeries.new_decimal,
    }


def py_type_to_constructor(py_type: type[Any]) -> Callable[..., PySeries]:
    """Get the right PySeries constructor for the given Python dtype."""
    py_type = (
        next((tp for tp in _PY_TYPE_TO_CONSTRUCTOR if issubclass(py_type, tp)), py_type)
        if py_type not in _PY_TYPE_TO_CONSTRUCTOR
        else py_type
    )
    return _PY_TYPE_TO_CONSTRUCTOR.get(py_type, PySeries.new_object)
