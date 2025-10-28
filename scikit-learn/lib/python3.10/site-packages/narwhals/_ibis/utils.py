from __future__ import annotations

from functools import lru_cache, partial
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import ibis
import ibis.expr.datatypes as ibis_dtypes

from narwhals._utils import Version, isinstance_or_issubclass

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from datetime import timedelta

    import ibis.expr.types as ir
    from ibis.common.temporal import TimestampUnit
    from ibis.expr.datatypes import DataType as IbisDataType
    from typing_extensions import TypeAlias, TypeIs

    from narwhals._duration import IntervalUnit
    from narwhals._ibis.dataframe import IbisLazyFrame
    from narwhals._ibis.expr import IbisExpr
    from narwhals.dtypes import DType
    from narwhals.typing import IntoDType, PythonLiteral

IntoColumn: TypeAlias = "str | ir.Value | ir.Column"
SortFn: TypeAlias = "Callable[[IntoColumn], ir.Column]"
Incomplete: TypeAlias = Any
"""Marker for upstream issues."""


@overload
def lit(value: bool, dtype: None = ...) -> ir.BooleanScalar: ...  # noqa: FBT001
@overload
def lit(value: int, dtype: None = ...) -> ir.IntegerScalar: ...
@overload
def lit(value: float, dtype: None = ...) -> ir.FloatingScalar: ...
@overload
def lit(value: str, dtype: None = ...) -> ir.StringScalar: ...
@overload
def lit(value: PythonLiteral | ir.Value, dtype: None = ...) -> ir.Scalar: ...
@overload
def lit(value: Any, dtype: Any) -> Incomplete: ...
def lit(value: Any, dtype: Any | None = None) -> Incomplete:
    """Alias for `ibis.literal`."""
    literal: Incomplete = ibis.literal
    return literal(value, dtype)


asc_nulls_first = cast("SortFn", partial(ibis.asc, nulls_first=True))
asc_nulls_last = cast("SortFn", partial(ibis.asc, nulls_first=False))
desc_nulls_first = cast("SortFn", partial(ibis.desc, nulls_first=True))
desc_nulls_last = cast("SortFn", partial(ibis.desc, nulls_first=False))


BucketUnit: TypeAlias = Literal[
    "years",
    "quarters",
    "months",
    "days",
    "hours",
    "minutes",
    "seconds",
    "milliseconds",
    "microseconds",
    "nanoseconds",
]
TruncateUnit: TypeAlias = Literal[
    "Y", "Q", "M", "W", "D", "h", "m", "s", "ms", "us", "ns"
]

UNITS_DICT_BUCKET: Mapping[IntervalUnit, BucketUnit] = {
    "y": "years",
    "q": "quarters",
    "mo": "months",
    "d": "days",
    "h": "hours",
    "m": "minutes",
    "s": "seconds",
    "ms": "milliseconds",
    "us": "microseconds",
    "ns": "nanoseconds",
}

UNITS_DICT_TRUNCATE: Mapping[IntervalUnit, TruncateUnit] = {
    "y": "Y",
    "q": "Q",
    "mo": "M",
    "d": "D",
    "h": "h",
    "m": "m",
    "s": "s",
    "ms": "ms",
    "us": "us",
    "ns": "ns",
}

FUNCTION_REMAPPING = {
    "starts_with": "startswith",
    "ends_with": "endswith",
    "regexp_matches": "re_search",
    "str_split": "split",
    "dayofyear": "day_of_year",
    "to_date": "date",
}


def evaluate_exprs(df: IbisLazyFrame, /, *exprs: IbisExpr) -> list[tuple[str, ir.Value]]:
    native_results: list[tuple[str, ir.Value]] = []
    for expr in exprs:
        native_series_list = expr(df)
        output_names = expr._evaluate_output_names(df)
        if expr._alias_output_names is not None:
            output_names = expr._alias_output_names(output_names)
        if len(output_names) != len(native_series_list):  # pragma: no cover
            msg = f"Internal error: got output names {output_names}, but only got {len(native_series_list)} results"
            raise AssertionError(msg)
        native_results.extend(zip(output_names, native_series_list))
    return native_results


@lru_cache(maxsize=16)
def native_to_narwhals_dtype(ibis_dtype: IbisDataType, version: Version) -> DType:  # noqa: C901, PLR0912
    dtypes = version.dtypes
    if ibis_dtype.is_int64():
        return dtypes.Int64()
    if ibis_dtype.is_int32():
        return dtypes.Int32()
    if ibis_dtype.is_int16():
        return dtypes.Int16()
    if ibis_dtype.is_int8():
        return dtypes.Int8()
    if ibis_dtype.is_uint64():
        return dtypes.UInt64()
    if ibis_dtype.is_uint32():
        return dtypes.UInt32()
    if ibis_dtype.is_uint16():
        return dtypes.UInt16()
    if ibis_dtype.is_uint8():
        return dtypes.UInt8()
    if ibis_dtype.is_boolean():
        return dtypes.Boolean()
    if ibis_dtype.is_float64():
        return dtypes.Float64()
    if ibis_dtype.is_float32():
        return dtypes.Float32()
    if ibis_dtype.is_string():
        return dtypes.String()
    if ibis_dtype.is_date():
        return dtypes.Date()
    if is_timestamp(ibis_dtype):
        _unit = cast("TimestampUnit", ibis_dtype.unit)
        return dtypes.Datetime(time_unit=_unit.value, time_zone=ibis_dtype.timezone)
    if is_interval(ibis_dtype):
        _time_unit = ibis_dtype.unit.value
        if _time_unit not in {"ns", "us", "ms", "s"}:  # pragma: no cover
            msg = f"Unsupported interval unit: {_time_unit}"
            raise NotImplementedError(msg)
        return dtypes.Duration(_time_unit)
    if is_array(ibis_dtype):
        if ibis_dtype.length:
            return dtypes.Array(
                native_to_narwhals_dtype(ibis_dtype.value_type, version),
                ibis_dtype.length,
            )
        return dtypes.List(native_to_narwhals_dtype(ibis_dtype.value_type, version))
    if is_struct(ibis_dtype):
        return dtypes.Struct(
            [
                dtypes.Field(name, native_to_narwhals_dtype(dtype, version))
                for name, dtype in ibis_dtype.items()
            ]
        )
    if ibis_dtype.is_decimal():  # pragma: no cover
        return dtypes.Decimal()
    if ibis_dtype.is_time():
        return dtypes.Time()
    if ibis_dtype.is_binary():
        return dtypes.Binary()
    return dtypes.Unknown()  # pragma: no cover


def is_timestamp(obj: IbisDataType) -> TypeIs[ibis_dtypes.Timestamp]:
    return obj.is_timestamp()


def is_interval(obj: IbisDataType) -> TypeIs[ibis_dtypes.Interval]:
    return obj.is_interval()


def is_array(obj: IbisDataType) -> TypeIs[ibis_dtypes.Array[Any]]:
    return obj.is_array()


def is_struct(obj: IbisDataType) -> TypeIs[ibis_dtypes.Struct]:
    return obj.is_struct()


def is_floating(obj: IbisDataType) -> TypeIs[ibis_dtypes.Floating]:
    return obj.is_floating()


dtypes = Version.MAIN.dtypes
NW_TO_IBIS_DTYPES: Mapping[type[DType], IbisDataType] = {
    dtypes.Float64: ibis_dtypes.Float64(),
    dtypes.Float32: ibis_dtypes.Float32(),
    dtypes.Binary: ibis_dtypes.Binary(),
    dtypes.String: ibis_dtypes.String(),
    dtypes.Boolean: ibis_dtypes.Boolean(),
    dtypes.Date: ibis_dtypes.Date(),
    dtypes.Time: ibis_dtypes.Time(),
    dtypes.Int8: ibis_dtypes.Int8(),
    dtypes.Int16: ibis_dtypes.Int16(),
    dtypes.Int32: ibis_dtypes.Int32(),
    dtypes.Int64: ibis_dtypes.Int64(),
    dtypes.UInt8: ibis_dtypes.UInt8(),
    dtypes.UInt16: ibis_dtypes.UInt16(),
    dtypes.UInt32: ibis_dtypes.UInt32(),
    dtypes.UInt64: ibis_dtypes.UInt64(),
    dtypes.Decimal: ibis_dtypes.Decimal(),
}
# Enum support: https://github.com/ibis-project/ibis/issues/10991
UNSUPPORTED_DTYPES = (dtypes.Int128, dtypes.UInt128, dtypes.Categorical, dtypes.Enum)


def narwhals_to_native_dtype(dtype: IntoDType, version: Version) -> IbisDataType:
    dtypes = version.dtypes
    base_type = dtype.base_type()
    if ibis_type := NW_TO_IBIS_DTYPES.get(base_type):
        return ibis_type
    if isinstance_or_issubclass(dtype, dtypes.Datetime):
        return ibis_dtypes.Timestamp.from_unit(dtype.time_unit, timezone=dtype.time_zone)
    if isinstance_or_issubclass(dtype, dtypes.Duration):
        return ibis_dtypes.Interval(unit=dtype.time_unit)  # pyright: ignore[reportArgumentType]
    if isinstance_or_issubclass(dtype, dtypes.List):
        inner = narwhals_to_native_dtype(dtype.inner, version)
        return ibis_dtypes.Array(value_type=inner)
    if isinstance_or_issubclass(dtype, dtypes.Struct):
        fields = [
            (field.name, narwhals_to_native_dtype(field.dtype, version))
            for field in dtype.fields
        ]
        return ibis_dtypes.Struct.from_tuples(fields)
    if isinstance_or_issubclass(dtype, dtypes.Array):
        inner = narwhals_to_native_dtype(dtype.inner, version)
        return ibis_dtypes.Array(value_type=inner, length=dtype.size)
    if issubclass(base_type, UNSUPPORTED_DTYPES):
        msg = f"Converting to {base_type.__name__} dtype is not supported for Ibis."
        raise NotImplementedError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)


def timedelta_to_ibis_interval(td: timedelta) -> ibis.expr.types.temporal.IntervalScalar:
    return ibis.interval(days=td.days, seconds=td.seconds, microseconds=td.microseconds)


def function(name: str, *args: ir.Value | PythonLiteral) -> ir.Value:
    # Workaround SQL vs Ibis differences.
    if name == "row_number":
        return ibis.row_number() + lit(1)
    if name == "least":
        return ibis.least(*args)
    if name == "greatest":
        return ibis.greatest(*args)
    expr = args[0]
    if name == "var_pop":
        return cast("ir.NumericColumn", expr).var(how="pop")
    if name == "var_samp":
        return cast("ir.NumericColumn", expr).var(how="sample")
    if name == "stddev_pop":
        return cast("ir.NumericColumn", expr).std(how="pop")
    if name == "stddev_samp":
        return cast("ir.NumericColumn", expr).std(how="sample")
    if name == "substr":
        # Ibis is 0-indexed here, SQL is 1-indexed
        return cast("ir.StringColumn", expr).substr(args[1] - 1, *args[2:])  # type: ignore[operator]  # pyright: ignore[reportArgumentType]
    return getattr(expr, FUNCTION_REMAPPING.get(name, name))(*args[1:])
