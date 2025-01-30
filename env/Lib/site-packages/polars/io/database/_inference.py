from __future__ import annotations

import functools
import re
from contextlib import suppress
from inspect import isclass
from typing import TYPE_CHECKING, Any

from polars.datatypes import (
    Binary,
    Boolean,
    Date,
    Datetime,
    Decimal,
    Duration,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    List,
    Null,
    String,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
)
from polars.datatypes._parse import parse_py_type_into_dtype
from polars.datatypes.group import (
    INTEGER_DTYPES,
    UNSIGNED_INTEGER_DTYPES,
)

if TYPE_CHECKING:
    from polars._typing import PolarsDataType


def _infer_dtype_from_database_typename(
    value: str,
    *,
    raise_unmatched: bool = True,
) -> PolarsDataType | None:
    """
    Attempt to infer Polars dtype from database cursor `type_code` string value.

    Examples
    --------
    >>> _infer_dtype_from_database_typename("INT2")
    Int16
    >>> _infer_dtype_from_database_typename("NVARCHAR")
    String
    >>> _infer_dtype_from_database_typename("NUMERIC(10,2)")
    Decimal(precision=10, scale=2)
    >>> _infer_dtype_from_database_typename("TIMESTAMP WITHOUT TZ")
    Datetime(time_unit='us', time_zone=None)
    """
    dtype: PolarsDataType | None = None

    # normalise string name/case (eg: 'IntegerType' -> 'INTEGER')
    original_value = value
    value = value.upper().replace("TYPE", "")

    # extract optional type modifier (eg: 'VARCHAR(64)' -> '64')
    if re.search(r"\([\w,: ]+\)$", value):
        modifier = value[value.find("(") + 1 : -1]
        value = value.split("(")[0]
    elif (
        not value.startswith(("<", ">")) and re.search(r"\[[\w,\]\[: ]+]$", value)
    ) or value.endswith(("[S]", "[MS]", "[US]", "[NS]")):
        modifier = value[value.find("[") + 1 : -1]
        value = value.split("[")[0]
    else:
        modifier = ""

    # array dtypes
    array_aliases = ("ARRAY", "LIST", "[]")
    if value.endswith(array_aliases) or value.startswith(array_aliases):
        for a in array_aliases:
            value = value.replace(a, "", 1) if value else ""

        nested: PolarsDataType | None = None
        if not value and modifier:
            nested = _infer_dtype_from_database_typename(
                value=modifier,
                raise_unmatched=False,
            )
        else:
            if inner_value := _infer_dtype_from_database_typename(
                value[1:-1]
                if (value[0], value[-1]) == ("<", ">")
                else re.sub(r"\W", "", re.sub(r"\WOF\W", "", value)),
                raise_unmatched=False,
            ):
                nested = inner_value
            elif modifier:
                nested = _infer_dtype_from_database_typename(
                    value=modifier,
                    raise_unmatched=False,
                )
        if nested:
            dtype = List(nested)

    # float dtypes
    elif value.startswith("FLOAT") or ("DOUBLE" in value) or (value == "REAL"):
        dtype = (
            Float32
            if value == "FLOAT4"
            or (value.endswith(("16", "32")) or (modifier in ("16", "32")))
            else Float64
        )

    # integer dtypes
    elif ("INTERVAL" not in value) and (
        value.startswith(("INT", "UINT", "UNSIGNED"))
        or value.endswith(("INT", "SERIAL"))
        or ("INTEGER" in value)
        or value == "ROWID"
    ):
        sz: Any
        if "LARGE" in value or value.startswith("BIG") or value == "INT8":
            sz = 64
        elif "MEDIUM" in value or value in ("INT4", "SERIAL"):
            sz = 32
        elif "SMALL" in value or value == "INT2":
            sz = 16
        elif "TINY" in value:
            sz = 8
        else:
            sz = None

        sz = modifier if (not sz and modifier) else sz
        if not isinstance(sz, int):
            sz = int(sz) if isinstance(sz, str) and sz.isdigit() else None
        if (
            ("U" in value and "MEDIUM" not in value)
            or ("UNSIGNED" in value)
            or value == "ROWID"
        ):
            dtype = _integer_dtype_from_nbits(sz, unsigned=True, default=UInt64)
        else:
            dtype = _integer_dtype_from_nbits(sz, unsigned=False, default=Int64)

    # number types (note: 'number' alone is not that helpful and requires refinement)
    elif "NUMBER" in value and "CARDINAL" in value:
        dtype = UInt64

    # decimal dtypes
    elif (is_dec := ("DECIMAL" in value)) or ("NUMERIC" in value):
        if "," in modifier:
            prec, scale = modifier.split(",")
            dtype = Decimal(int(prec), int(scale))
        else:
            dtype = Decimal if is_dec else Float64

    # string dtypes
    elif (
        any(tp in value for tp in ("VARCHAR", "STRING", "TEXT", "UNICODE"))
        or value.startswith(("STR", "CHAR", "BPCHAR", "NCHAR", "UTF"))
        or value.endswith(("_UTF8", "_UTF16", "_UTF32"))
    ):
        dtype = String

    # binary dtypes
    elif value in ("BYTEA", "BYTES", "BLOB", "CLOB", "BINARY"):
        dtype = Binary

    # boolean dtypes
    elif value.startswith("BOOL"):
        dtype = Boolean

    # null dtype; odd, but valid
    elif value == "NULL":
        dtype = Null

    # temporal dtypes
    elif value.startswith(("DATETIME", "TIMESTAMP")) and not (value.endswith("[D]")):
        if any((tz in value.replace(" ", "")) for tz in ("TZ", "TIMEZONE")):
            if "WITHOUT" not in value:
                return None  # there's a timezone, but we don't know what it is
        unit = _timeunit_from_precision(modifier) if modifier else "us"
        dtype = Datetime(time_unit=(unit or "us"))  # type: ignore[arg-type]
    else:
        value = re.sub(r"\d", "", value)
        if value in ("INTERVAL", "TIMEDELTA", "DURATION"):
            dtype = Duration
        elif value == "DATE":
            dtype = Date
        elif value == "TIME":
            dtype = Time

    if not dtype and raise_unmatched:
        msg = f"cannot infer dtype from {original_value!r} string value"
        raise ValueError(msg)

    return dtype


def _infer_dtype_from_cursor_description(
    cursor: Any,
    description: tuple[Any, ...],
) -> PolarsDataType | None:
    """Attempt to infer Polars dtype from database cursor description `type_code`."""
    type_code, _disp_size, internal_size, precision, scale, *_ = description
    dtype: PolarsDataType | None = None

    if isclass(type_code):
        # python types, eg: int, float, str, etc
        with suppress(TypeError):
            dtype = parse_py_type_into_dtype(type_code)  # type: ignore[arg-type]

    elif isinstance(type_code, str):
        # database/sql type names, eg: "VARCHAR", "NUMERIC", "BLOB", etc
        dtype = _infer_dtype_from_database_typename(
            value=type_code,
            raise_unmatched=False,
        )

    # check additional cursor attrs to refine dtype specification
    if dtype is not None:
        if dtype == Float64 and internal_size == 4:
            dtype = Float32

        elif dtype in INTEGER_DTYPES and internal_size in (2, 4, 8):
            bits = internal_size * 8
            dtype = _integer_dtype_from_nbits(
                bits,
                unsigned=(dtype in UNSIGNED_INTEGER_DTYPES),
                default=dtype,
            )
        elif (
            dtype == Decimal
            and isinstance(precision, int)
            and isinstance(scale, int)
            and precision <= 38
            and scale <= 38
        ):
            dtype = Decimal(precision, scale)

    return dtype


@functools.lru_cache(8)
def _integer_dtype_from_nbits(
    bits: int,
    *,
    unsigned: bool,
    default: PolarsDataType | None = None,
) -> PolarsDataType | None:
    """
    Return matching Polars integer dtype from num bits and signed/unsigned flag.

    Examples
    --------
    >>> _integer_dtype_from_nbits(8, unsigned=False)
    Int8
    >>> _integer_dtype_from_nbits(32, unsigned=True)
    UInt32
    """
    dtype = {
        (8, False): Int8,
        (8, True): UInt8,
        (16, False): Int16,
        (16, True): UInt16,
        (32, False): Int32,
        (32, True): UInt32,
        (64, False): Int64,
        (64, True): UInt64,
    }.get((bits, unsigned), None)

    if dtype is None and default is not None:
        return default
    return dtype


def _timeunit_from_precision(precision: int | str | None) -> str | None:
    """
    Return `time_unit` from integer precision value.

    Examples
    --------
    >>> _timeunit_from_precision(3)
    'ms'
    >>> _timeunit_from_precision(5)
    'us'
    >>> _timeunit_from_precision(7)
    'ns'
    """
    from math import ceil

    if not precision:
        return None
    elif isinstance(precision, str):
        if precision.isdigit():
            precision = int(precision)
        elif (precision := precision.lower()) in ("s", "ms", "us", "ns"):
            return "ms" if precision == "s" else precision
    try:
        n = min(max(3, int(ceil(precision / 3)) * 3), 9)  # type: ignore[operator]
        return {3: "ms", 6: "us", 9: "ns"}.get(n)
    except TypeError:
        return None
