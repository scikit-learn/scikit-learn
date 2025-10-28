from __future__ import annotations

import contextlib
import functools
import re
import sys
from collections.abc import Collection
from datetime import date, datetime, time, timedelta
from decimal import Decimal as PyDecimal
from typing import TYPE_CHECKING, Any, Optional, Union

from polars._dependencies import numpy as np
from polars._dependencies import pyarrow as pa
from polars.datatypes.classes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    DataType,
    DataTypeClass,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Field,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Null,
    Object,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Unknown,
)

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import dtype_str_repr as _dtype_str_repr


OptionType = type(Optional[type])
if sys.version_info >= (3, 10):
    from types import NoneType, UnionType
else:
    # infer equivalent class
    NoneType = type(None)
    UnionType = type(Union[int, float])

if TYPE_CHECKING:
    from polars._typing import PolarsDataType, PythonDataType, TimeUnit

    if sys.version_info >= (3, 10):
        from typing import TypeGuard
    else:
        from typing_extensions import TypeGuard


def is_polars_dtype(
    dtype: Any,
    *,
    include_unknown: bool = False,
    require_instantiated: bool = False,
) -> TypeGuard[PolarsDataType]:
    """Indicate whether the given input is a Polars dtype, or dtype specialization."""
    check_classes = DataType if require_instantiated else (DataType, DataTypeClass)
    is_dtype = isinstance(dtype, check_classes)

    if not include_unknown:
        return is_dtype and dtype != Unknown
    else:
        return is_dtype


def unpack_dtypes(
    *dtypes: PolarsDataType | None,
    include_compound: bool = False,
) -> set[PolarsDataType]:
    """
    Return a set of unique dtypes found in one or more (potentially compound) dtypes.

    Parameters
    ----------
    *dtypes
        One or more Polars dtypes.
    include_compound
        * if True, any parent/compound dtypes (List, Struct) are included in the result.
        * if False, only the child/scalar dtypes are returned from these types.

    Examples
    --------
    >>> from polars.datatypes import unpack_dtypes
    >>> list_dtype = [pl.List(pl.Float64)]
    >>> struct_dtype = pl.Struct(
    ...     [
    ...         pl.Field("a", pl.Int64),
    ...         pl.Field("b", pl.String),
    ...         pl.Field("c", pl.List(pl.Float64)),
    ...     ]
    ... )
    >>> unpack_dtypes([struct_dtype, list_dtype])  # doctest: +IGNORE_RESULT
    {Float64, Int64, String}
    >>> unpack_dtypes(
    ...     [struct_dtype, list_dtype], include_compound=True
    ... )  # doctest: +IGNORE_RESULT
    {Float64, Int64, String, List(Float64), Struct([Field('a', Int64), Field('b', String), Field('c', List(Float64))])}
    """  # noqa: W505
    if not dtypes:
        return set()
    elif len(dtypes) == 1 and isinstance(dtypes[0], Collection):
        dtypes = dtypes[0]

    unpacked: set[PolarsDataType] = set()
    for tp in dtypes:
        if isinstance(tp, (List, Array)):
            if include_compound:
                unpacked.add(tp)
            unpacked.update(unpack_dtypes(tp.inner, include_compound=include_compound))
        elif isinstance(tp, Struct):
            if include_compound:
                unpacked.add(tp)
            unpacked.update(unpack_dtypes(tp.fields, include_compound=include_compound))  # type: ignore[arg-type]
        elif isinstance(tp, Field):
            unpacked.update(unpack_dtypes(tp.dtype, include_compound=include_compound))
        elif tp is not None and is_polars_dtype(tp):
            unpacked.add(tp)
    return unpacked


class _DataTypeMappings:
    @property
    @functools.lru_cache  # noqa: B019
    def DTYPE_TO_FFINAME(self) -> dict[PolarsDataType, str]:
        return {
            Binary: "binary",
            Boolean: "bool",
            Categorical: "categorical",
            Date: "date",
            Datetime: "datetime",
            Decimal: "decimal",
            Duration: "duration",
            Float32: "f32",
            Float64: "f64",
            Int8: "i8",
            Int16: "i16",
            Int32: "i32",
            Int64: "i64",
            Int128: "i128",
            List: "list",
            Object: "object",
            String: "str",
            Struct: "struct",
            Time: "time",
            UInt8: "u8",
            UInt16: "u16",
            UInt32: "u32",
            UInt64: "u64",
            UInt128: "u128",
        }

    @property
    @functools.lru_cache  # noqa: B019
    def DTYPE_TO_PY_TYPE(self) -> dict[PolarsDataType, PythonDataType]:
        return {
            Array: list,
            Binary: bytes,
            Boolean: bool,
            Date: date,
            Datetime: datetime,
            Decimal: PyDecimal,
            Duration: timedelta,
            Float32: float,
            Float64: float,
            Int8: int,
            Int16: int,
            Int32: int,
            Int64: int,
            Int128: int,
            List: list,
            Null: None.__class__,
            Object: object,
            String: str,
            Struct: dict,
            Time: time,
            UInt8: int,
            UInt16: int,
            UInt32: int,
            UInt64: int,
            UInt128: int,
            # the below mappings are appropriate as we restrict cat/enum to strings
            Enum: str,
            Categorical: str,
        }

    @property
    @functools.lru_cache  # noqa: B019
    def NUMPY_KIND_AND_ITEMSIZE_TO_DTYPE(self) -> dict[tuple[str, int], PolarsDataType]:
        return {
            # (np.dtype().kind, np.dtype().itemsize)
            ("M", 8): Datetime,
            ("b", 1): Boolean,
            ("f", 2): Float32,
            ("f", 4): Float32,
            ("f", 8): Float64,
            ("i", 1): Int8,
            ("i", 2): Int16,
            ("i", 4): Int32,
            ("i", 8): Int64,
            ("m", 8): Duration,
            ("u", 1): UInt8,
            ("u", 2): UInt16,
            ("u", 4): UInt32,
            ("u", 8): UInt64,
        }

    @property
    @functools.lru_cache  # noqa: B019
    def PY_TYPE_TO_ARROW_TYPE(self) -> dict[PythonDataType, pa.lib.DataType]:
        return {
            bool: pa.bool_(),
            date: pa.date32(),
            datetime: pa.timestamp("us"),
            float: pa.float64(),
            int: pa.int64(),
            str: pa.large_utf8(),
            time: pa.time64("us"),
            timedelta: pa.duration("us"),
            None.__class__: pa.null(),
        }

    @property
    @functools.lru_cache  # noqa: B019
    def REPR_TO_DTYPE(self) -> dict[str, PolarsDataType]:
        def _dtype_str_repr_safe(o: Any) -> PolarsDataType | None:
            try:
                return _dtype_str_repr(o.base_type()).split("[")[0]  # type: ignore[return-value]
            except TypeError:
                return None

        return {
            _dtype_str_repr_safe(obj): obj  # type: ignore[misc]
            for obj in globals().values()
            if is_polars_dtype(obj) and _dtype_str_repr_safe(obj) is not None
        }


# Initialize once (poor man's singleton :)
DataTypeMappings = _DataTypeMappings()


def dtype_to_ffiname(dtype: PolarsDataType) -> str:
    """Return FFI function name associated with the given Polars dtype."""
    try:
        dtype = dtype.base_type()
        return DataTypeMappings.DTYPE_TO_FFINAME[dtype]
    except KeyError:  # pragma: no cover
        msg = f"conversion of polars data type {dtype!r} to FFI not implemented"
        raise NotImplementedError(msg) from None


def dtype_to_py_type(dtype: PolarsDataType) -> PythonDataType:
    """Convert a Polars dtype to a Python dtype."""
    try:
        dtype = dtype.base_type()
        return DataTypeMappings.DTYPE_TO_PY_TYPE[dtype]
    except KeyError:  # pragma: no cover
        msg = f"conversion of polars data type {dtype!r} to Python type not implemented"
        raise NotImplementedError(msg) from None


def py_type_to_arrow_type(dtype: PythonDataType) -> pa.lib.DataType:
    """Convert a Python dtype to an Arrow dtype."""
    try:
        return DataTypeMappings.PY_TYPE_TO_ARROW_TYPE[dtype]
    except KeyError:  # pragma: no cover
        msg = f"cannot parse Python data type {dtype!r} into Arrow data type"
        raise ValueError(msg) from None


def dtype_short_repr_to_dtype(dtype_string: str | None) -> PolarsDataType | None:
    """Map a PolarsDataType short repr (eg: 'i64', 'list[str]') back into a dtype."""
    if dtype_string is None:
        return None

    m = re.match(r"^(\w+)(?:\[(.+)\])?$", dtype_string)
    if m is None:
        return None

    dtype_base, subtype = m.groups()
    dtype = DataTypeMappings.REPR_TO_DTYPE.get(dtype_base)
    if dtype and subtype:
        # TODO: further-improve handling for nested types (such as List,Struct)
        try:
            if dtype == Decimal:
                subtype = (None, int(subtype))
            else:
                subtype = (
                    s.strip("'\" ") for s in subtype.replace("μs", "us").split(",")
                )
            return dtype(*subtype)  # type: ignore[operator]
        except ValueError:
            pass
    return dtype


def supported_numpy_char_code(dtype_char: str) -> bool:
    """Check if the input can be mapped to a Polars dtype."""
    dtype = np.dtype(dtype_char)
    return (
        dtype.kind,
        dtype.itemsize,
    ) in DataTypeMappings.NUMPY_KIND_AND_ITEMSIZE_TO_DTYPE


def numpy_char_code_to_dtype(dtype_char: str) -> PolarsDataType:
    """Convert a numpy character dtype to a Polars dtype."""
    dtype = np.dtype(dtype_char)
    if dtype.kind == "U":
        return String
    elif dtype.kind == "S":
        return Binary
    try:
        return DataTypeMappings.NUMPY_KIND_AND_ITEMSIZE_TO_DTYPE[
            dtype.kind, dtype.itemsize
        ]
    except KeyError:  # pragma: no cover
        msg = f"cannot parse numpy data type {dtype!r} into Polars data type"
        raise ValueError(msg) from None


def maybe_cast(el: Any, dtype: PolarsDataType) -> Any:
    """Try casting a value to a value that is valid for the given Polars dtype."""
    # cast el if it doesn't match
    from polars._utils.convert import (
        datetime_to_int,
        timedelta_to_int,
    )

    time_unit: TimeUnit
    if isinstance(el, datetime):
        time_unit = getattr(dtype, "time_unit", "us")
        return datetime_to_int(el, time_unit)
    elif isinstance(el, timedelta):
        time_unit = getattr(dtype, "time_unit", "us")
        return timedelta_to_int(el, time_unit)

    py_type = dtype_to_py_type(dtype)
    if not isinstance(el, py_type):
        try:
            el = py_type(el)  # type: ignore[call-arg]
        except Exception:
            from polars._utils.various import qualified_type_name

            msg = f"cannot convert Python type {qualified_type_name(el)!r} to {dtype!r}"
            raise TypeError(msg) from None
    return el
