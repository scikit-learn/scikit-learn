from __future__ import annotations

import enum
import functools
import re
import sys
from datetime import date, datetime, time, timedelta
from decimal import Decimal as PyDecimal
from inspect import isclass
from typing import TYPE_CHECKING, Any, ForwardRef, NoReturn, Union, get_args

import polars._reexport as pl
from polars.datatypes.classes import (
    Binary,
    Boolean,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Float64,
    Int64,
    List,
    Null,
    Object,
    String,
    Time,
    Unknown,
)
from polars.datatypes.convert import is_polars_dtype

if TYPE_CHECKING:
    from polars._typing import PolarsDataType, PythonDataType, SchemaDict


UnionTypeOld = type(Union[int, str])
if sys.version_info >= (3, 10):
    from types import NoneType, UnionType
else:  # pragma: no cover
    # Define equivalent for older Python versions
    NoneType = type(None)
    UnionType = UnionTypeOld


def parse_into_datatype_expr(input: Any) -> pl.DataTypeExpr:
    """Parse an input into a DataTypeExpr."""
    if isinstance(input, pl.DataTypeExpr):
        return input
    else:
        return parse_into_dtype(input).to_dtype_expr()


def parse_into_dtype(input: Any) -> PolarsDataType:
    """
    Parse an input into a Polars data type.

    Raises
    ------
    TypeError
        If the input cannot be parsed into a Polars data type.
    """
    if is_polars_dtype(input):
        return input
    elif isinstance(input, ForwardRef):
        return _parse_forward_ref_into_dtype(input)
    elif isinstance(input, (UnionType, UnionTypeOld)):
        return _parse_union_type_into_dtype(input)
    else:
        return parse_py_type_into_dtype(input)


def try_parse_into_dtype(input: Any) -> PolarsDataType | None:
    """Try parsing an input into a Polars data type, returning None on failure."""
    try:
        return parse_into_dtype(input)
    except TypeError:
        return None


@functools.lru_cache(16)
def parse_py_type_into_dtype(input: PythonDataType | type[object]) -> PolarsDataType:
    """Convert Python data type to Polars data type."""
    if input is int:
        return Int64()
    elif input is float:
        return Float64()
    elif input is str:
        return String()
    elif input is bool:
        return Boolean()
    elif isinstance(input, type) and issubclass(input, datetime):  # type: ignore[redundant-expr]
        return Datetime("us")
    elif isinstance(input, type) and issubclass(input, date):  # type: ignore[redundant-expr]
        return Date()
    elif input is timedelta:
        return Duration
    elif input is time:
        return Time()
    elif input is PyDecimal:
        return Decimal
    elif input is bytes:
        return Binary()
    elif input is object:
        return Object()
    elif input is NoneType:
        return Null()
    elif input is list or input is tuple:
        return List
    elif isclass(input) and issubclass(input, enum.Enum):
        return Enum(input)
    # this is required as pass through. Don't remove
    elif input == Unknown:
        return Unknown
    elif hasattr(input, "__origin__") and hasattr(input, "__args__"):
        return _parse_generic_into_dtype(input)
    else:
        _raise_on_invalid_dtype(input)


def _parse_generic_into_dtype(input: Any) -> PolarsDataType:
    """Parse a generic type (from typing annotation) into a Polars data type."""
    base_type = input.__origin__
    if base_type not in (tuple, list):
        _raise_on_invalid_dtype(input)

    inner_types = input.__args__
    inner_type = inner_types[0]
    if len(inner_types) > 1:
        all_equal = all(t in (inner_type, ...) for t in inner_types)
        if not all_equal:
            _raise_on_invalid_dtype(input)

    inner_type = inner_types[0]
    inner_dtype = parse_py_type_into_dtype(inner_type)
    return List(inner_dtype)


PY_TYPE_STR_TO_DTYPE: SchemaDict = {
    "Decimal": Decimal,
    "NoneType": Null(),
    "bool": Boolean(),
    "bytes": Binary(),
    "date": Date(),
    "datetime": Datetime("us"),
    "float": Float64(),
    "int": Int64(),
    "list": List,
    "object": Object(),
    "str": String(),
    "time": Time(),
    "timedelta": Duration,
    "tuple": List,
}


def _parse_forward_ref_into_dtype(input: ForwardRef) -> PolarsDataType:
    """Parse a ForwardRef into a Polars data type."""
    annotation = input.__forward_arg__

    # Strip "optional" designation - Polars data types are always nullable
    formatted = re.sub(r"(^None \|)|(\| None$)", "", annotation).strip()

    try:
        return PY_TYPE_STR_TO_DTYPE[formatted]
    except KeyError:
        _raise_on_invalid_dtype(input)


def _parse_union_type_into_dtype(input: Any) -> PolarsDataType:
    """
    Parse a union of types into a Polars data type.

    Unions of multiple non-null types (e.g. `int | float`) are not supported.

    Parameters
    ----------
    input
        A union type, e.g. `str | None` (new syntax) or `Union[str, None]` (old syntax).
    """
    # Strip "optional" designation - Polars data types are always nullable
    inner_types = [tp for tp in get_args(input) if tp is not NoneType]

    if len(inner_types) != 1:
        _raise_on_invalid_dtype(input)

    input = inner_types[0]
    return parse_into_dtype(input)


def _raise_on_invalid_dtype(input: Any) -> NoReturn:
    """Raise an informative error if the input could not be parsed."""
    input_type = input if type(input) is type else f"of type {type(input).__name__!r}"
    input_detail = "" if type(input) is type else f" (given: {input!r})"
    msg = f"cannot parse input {input_type} into Polars data type{input_detail}"
    raise TypeError(msg) from None
