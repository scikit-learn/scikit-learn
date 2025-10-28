from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING, Any, cast

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._compliant import EagerSeriesNamespace
from narwhals._utils import Implementation, Version, isinstance_or_issubclass

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping

    from typing_extensions import TypeAlias, TypeIs

    from narwhals._arrow.series import ArrowSeries
    from narwhals._arrow.typing import (
        ArrayAny,
        ArrayOrScalar,
        ArrayOrScalarT1,
        ArrayOrScalarT2,
        ChunkedArrayAny,
        Incomplete,
        NativeIntervalUnit,
        PromoteOptions,
        ScalarAny,
    )
    from narwhals._duration import IntervalUnit
    from narwhals.dtypes import DType
    from narwhals.typing import IntoDType, PythonLiteral

    # NOTE: stubs don't allow for `ChunkedArray[StructArray]`
    # Intended to represent the `.chunks` property storing `list[pa.StructArray]`
    ChunkedArrayStructArray: TypeAlias = ChunkedArrayAny

    def is_timestamp(t: Any) -> TypeIs[pa.TimestampType[Any, Any]]: ...
    def is_duration(t: Any) -> TypeIs[pa.DurationType[Any]]: ...
    def is_list(t: Any) -> TypeIs[pa.ListType[Any]]: ...
    def is_large_list(t: Any) -> TypeIs[pa.LargeListType[Any]]: ...
    def is_fixed_size_list(t: Any) -> TypeIs[pa.FixedSizeListType[Any, Any]]: ...
    def is_dictionary(t: Any) -> TypeIs[pa.DictionaryType[Any, Any, Any]]: ...
    def extract_regex(
        strings: ChunkedArrayAny,
        /,
        pattern: str,
        *,
        options: Any = None,
        memory_pool: Any = None,
    ) -> ChunkedArrayStructArray: ...
else:
    from pyarrow.compute import extract_regex
    from pyarrow.types import (
        is_dictionary,  # noqa: F401
        is_duration,
        is_fixed_size_list,
        is_large_list,
        is_list,
        is_timestamp,
    )

BACKEND_VERSION = Implementation.PYARROW._backend_version()
"""Static backend version for `pyarrow`."""

UNITS_DICT: Mapping[IntervalUnit, NativeIntervalUnit] = {
    "y": "year",
    "q": "quarter",
    "mo": "month",
    "d": "day",
    "h": "hour",
    "m": "minute",
    "s": "second",
    "ms": "millisecond",
    "us": "microsecond",
    "ns": "nanosecond",
}

lit = pa.scalar
"""Alias for `pyarrow.scalar`."""


def extract_py_scalar(value: Any, /) -> Any:
    from narwhals._arrow.series import maybe_extract_py_scalar

    return maybe_extract_py_scalar(value, return_py_scalar=True)


def is_array_or_scalar(obj: Any) -> TypeIs[ArrayOrScalar]:
    """Return True for any base `pyarrow` container."""
    return isinstance(obj, (pa.ChunkedArray, pa.Array, pa.Scalar))


def chunked_array(
    arr: ArrayOrScalar | list[Iterable[Any]], dtype: pa.DataType | None = None, /
) -> ChunkedArrayAny:
    if isinstance(arr, pa.ChunkedArray):
        return arr
    if isinstance(arr, list):
        return pa.chunked_array(arr, dtype)
    return pa.chunked_array([arr], dtype)


def nulls_like(n: int, series: ArrowSeries) -> ArrayAny:
    """Create a strongly-typed Array instance with all elements null.

    Uses the type of `series`, without upseting `mypy`.
    """
    return pa.nulls(n, series.native.type)


def repeat(
    value: PythonLiteral | ScalarAny, n: int, /, dtype: pa.DataType | None = None
) -> ArrayAny:
    """Create an Array instance whose slots are the given scalar.

    *Optionally*, casting to `dtype` **before** repeating `n` times.
    """
    lit_: Incomplete = lit
    return pa.repeat(lit_(value, type=dtype), n)


def zeros(n: int, /) -> pa.Int64Array:
    return pa.repeat(0, n)


def native_to_narwhals_dtype(dtype: pa.DataType, version: Version) -> DType:
    if isinstance(dtype, pa.ExtensionType):
        return version.dtypes.Unknown()
    return native_non_extension_to_narwhals_dtype(dtype, version)


@lru_cache(maxsize=16)
def native_non_extension_to_narwhals_dtype(dtype: pa.DataType, version: Version) -> DType:  # noqa: C901, PLR0912
    dtypes = version.dtypes
    if pa.types.is_int64(dtype):
        return dtypes.Int64()
    if pa.types.is_int32(dtype):
        return dtypes.Int32()
    if pa.types.is_int16(dtype):
        return dtypes.Int16()
    if pa.types.is_int8(dtype):
        return dtypes.Int8()
    if pa.types.is_uint64(dtype):
        return dtypes.UInt64()
    if pa.types.is_uint32(dtype):
        return dtypes.UInt32()
    if pa.types.is_uint16(dtype):
        return dtypes.UInt16()
    if pa.types.is_uint8(dtype):
        return dtypes.UInt8()
    if pa.types.is_boolean(dtype):
        return dtypes.Boolean()
    if pa.types.is_float64(dtype):
        return dtypes.Float64()
    if pa.types.is_float32(dtype):
        return dtypes.Float32()
    # bug in coverage? it shows `31->exit` (where `31` is currently the line number of
    # the next line), even though both when the if condition is true and false are covered
    if (  # pragma: no cover
        pa.types.is_string(dtype)
        or pa.types.is_large_string(dtype)
        or getattr(pa.types, "is_string_view", lambda _: False)(dtype)
    ):
        return dtypes.String()
    if pa.types.is_date32(dtype):
        return dtypes.Date()
    if is_timestamp(dtype):
        return dtypes.Datetime(time_unit=dtype.unit, time_zone=dtype.tz)
    if is_duration(dtype):
        return dtypes.Duration(time_unit=dtype.unit)
    if pa.types.is_dictionary(dtype):
        return dtypes.Categorical()
    if pa.types.is_struct(dtype):
        return dtypes.Struct(
            [
                dtypes.Field(
                    dtype.field(i).name,
                    native_to_narwhals_dtype(dtype.field(i).type, version),
                )
                for i in range(dtype.num_fields)
            ]
        )
    if is_list(dtype) or is_large_list(dtype):
        return dtypes.List(native_to_narwhals_dtype(dtype.value_type, version))
    if is_fixed_size_list(dtype):
        return dtypes.Array(
            native_to_narwhals_dtype(dtype.value_type, version), dtype.list_size
        )
    if pa.types.is_decimal(dtype):
        return dtypes.Decimal()
    if pa.types.is_time32(dtype) or pa.types.is_time64(dtype):
        return dtypes.Time()
    if pa.types.is_binary(dtype):
        return dtypes.Binary()
    return dtypes.Unknown()  # pragma: no cover


dtypes = Version.MAIN.dtypes
NW_TO_PA_DTYPES: Mapping[type[DType], pa.DataType] = {
    dtypes.Float64: pa.float64(),
    dtypes.Float32: pa.float32(),
    dtypes.Binary: pa.binary(),
    dtypes.String: pa.string(),
    dtypes.Boolean: pa.bool_(),
    dtypes.Categorical: pa.dictionary(pa.uint32(), pa.string()),
    dtypes.Date: pa.date32(),
    dtypes.Time: pa.time64("ns"),
    dtypes.Int8: pa.int8(),
    dtypes.Int16: pa.int16(),
    dtypes.Int32: pa.int32(),
    dtypes.Int64: pa.int64(),
    dtypes.UInt8: pa.uint8(),
    dtypes.UInt16: pa.uint16(),
    dtypes.UInt32: pa.uint32(),
    dtypes.UInt64: pa.uint64(),
}
UNSUPPORTED_DTYPES = (dtypes.Decimal, dtypes.Object)


def narwhals_to_native_dtype(dtype: IntoDType, version: Version) -> pa.DataType:
    dtypes = version.dtypes
    base_type = dtype.base_type()
    if pa_type := NW_TO_PA_DTYPES.get(base_type):
        return pa_type
    if isinstance_or_issubclass(dtype, dtypes.Datetime):
        unit = dtype.time_unit
        return pa.timestamp(unit, tz) if (tz := dtype.time_zone) else pa.timestamp(unit)
    if isinstance_or_issubclass(dtype, dtypes.Duration):
        return pa.duration(dtype.time_unit)
    if isinstance_or_issubclass(dtype, dtypes.List):
        return pa.list_(value_type=narwhals_to_native_dtype(dtype.inner, version=version))
    if isinstance_or_issubclass(dtype, dtypes.Struct):
        return pa.struct(
            [
                (field.name, narwhals_to_native_dtype(field.dtype, version=version))
                for field in dtype.fields
            ]
        )
    if isinstance_or_issubclass(dtype, dtypes.Array):  # pragma: no cover
        inner = narwhals_to_native_dtype(dtype.inner, version=version)
        list_size = dtype.size
        return pa.list_(inner, list_size=list_size)
    if issubclass(base_type, UNSUPPORTED_DTYPES):
        msg = f"Converting to {base_type.__name__} dtype is not supported for PyArrow."
        raise NotImplementedError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)


def extract_native(
    lhs: ArrowSeries, rhs: ArrowSeries | PythonLiteral | ScalarAny
) -> tuple[ChunkedArrayAny | ScalarAny, ChunkedArrayAny | ScalarAny]:
    """Extract native objects in binary  operation.

    If the comparison isn't supported, return `NotImplemented` so that the
    "right-hand-side" operation (e.g. `__radd__`) can be tried.

    If one of the two sides has a `_broadcast` flag, then extract the scalar
    underneath it so that PyArrow can do its own broadcasting.
    """
    from narwhals._arrow.series import ArrowSeries

    if rhs is None:  # pragma: no cover
        return lhs.native, lit(None, type=lhs._type)

    if isinstance(rhs, ArrowSeries):
        if lhs._broadcast and not rhs._broadcast:
            return lhs.native[0], rhs.native
        if rhs._broadcast:
            return lhs.native, rhs.native[0]
        return lhs.native, rhs.native

    if isinstance(rhs, list):
        msg = "Expected Series or scalar, got list."
        raise TypeError(msg)

    return lhs.native, rhs if isinstance(rhs, pa.Scalar) else lit(rhs)


def floordiv_compat(left: ArrayOrScalar, right: ArrayOrScalar, /) -> Any:
    # The following lines are adapted from pandas' pyarrow implementation.
    # Ref: https://github.com/pandas-dev/pandas/blob/262fcfbffcee5c3116e86a951d8b693f90411e68/pandas/core/arrays/arrow/array.py#L124-L154
    # We modify it to add `safe_mask` so as to align with Polars' behaviour when dividing by zero.
    safe_mask = pc.not_equal(right, lit(0, type=right.type))
    safe_right = pc.if_else(safe_mask, right, lit(1, type=right.type))  # Dummy value
    non_safe_result = lit(None, type=left.type)

    if pa.types.is_integer(left.type) and pa.types.is_integer(right.type):
        divided = pc.if_else(
            safe_mask, pc.divide_checked(left, safe_right), non_safe_result
        )
        # TODO @dangotbanned: Use a `TypeVar` in guards
        # Narrowing to a `Union` isn't interacting well with the rest of the stubs
        # https://github.com/zen-xu/pyarrow-stubs/pull/215
        if pa.types.is_signed_integer(divided.type):
            div_type = cast("pa._lib.Int64Type", divided.type)
            has_remainder = pc.not_equal(pc.multiply(divided, right), left)
            has_one_negative_operand = pc.less(
                pc.bit_wise_xor(left, right), lit(0, div_type)
            )
            result = pc.if_else(
                pc.and_(has_remainder, has_one_negative_operand),
                pc.subtract(divided, lit(1, div_type)),
                divided,
            )
        else:
            result = divided  # pragma: no cover
        result = result.cast(left.type)
    else:
        divided = pc.if_else(safe_mask, pc.divide(left, safe_right), non_safe_result)
        result = pc.floor(divided)
    return result


def cast_for_truediv(
    arrow_array: ArrayOrScalarT1, pa_object: ArrayOrScalarT2
) -> tuple[ArrayOrScalarT1, ArrayOrScalarT2]:
    # Lifted from:
    # https://github.com/pandas-dev/pandas/blob/262fcfbffcee5c3116e86a951d8b693f90411e68/pandas/core/arrays/arrow/array.py#L108-L122
    # Ensure int / int -> float mirroring Python/Numpy behavior
    # as pc.divide_checked(int, int) -> int
    if pa.types.is_integer(arrow_array.type) and pa.types.is_integer(pa_object.type):
        # GH: 56645.  # noqa: ERA001
        # https://github.com/apache/arrow/issues/35563
        return arrow_array.cast(pa.float64(), safe=False), pa_object.cast(
            pa.float64(), safe=False
        )

    return arrow_array, pa_object


# Regex for date, time, separator and timezone components
DATE_RE = r"(?P<date>\d{1,4}[-/.]\d{1,2}[-/.]\d{1,4}|\d{8})"
SEP_RE = r"(?P<sep>\s|T)"
TIME_RE = r"(?P<time>\d{2}:\d{2}(?::\d{2})?|\d{6}?)"  # \s*(?P<period>[AP]M)?)?
HMS_RE = r"^(?P<hms>\d{2}:\d{2}:\d{2})$"
HM_RE = r"^(?P<hm>\d{2}:\d{2})$"
HMS_RE_NO_SEP = r"^(?P<hms_no_sep>\d{6})$"
TZ_RE = r"(?P<tz>Z|[+-]\d{2}:?\d{2})"  # Matches 'Z', '+02:00', '+0200', '+02', etc.
FULL_RE = rf"{DATE_RE}{SEP_RE}?{TIME_RE}?{TZ_RE}?$"

# Separate regexes for different date formats
YMD_RE = r"^(?P<year>(?:[12][0-9])?[0-9]{2})(?P<sep1>[-/.])(?P<month>0[1-9]|1[0-2])(?P<sep2>[-/.])(?P<day>0[1-9]|[12][0-9]|3[01])$"
DMY_RE = r"^(?P<day>0[1-9]|[12][0-9]|3[01])(?P<sep1>[-/.])(?P<month>0[1-9]|1[0-2])(?P<sep2>[-/.])(?P<year>(?:[12][0-9])?[0-9]{2})$"
MDY_RE = r"^(?P<month>0[1-9]|1[0-2])(?P<sep1>[-/.])(?P<day>0[1-9]|[12][0-9]|3[01])(?P<sep2>[-/.])(?P<year>(?:[12][0-9])?[0-9]{2})$"
YMD_RE_NO_SEP = r"^(?P<year>(?:[12][0-9])?[0-9]{2})(?P<month>0[1-9]|1[0-2])(?P<day>0[1-9]|[12][0-9]|3[01])$"

DATE_FORMATS = (
    (YMD_RE_NO_SEP, "%Y%m%d"),
    (YMD_RE, "%Y-%m-%d"),
    (DMY_RE, "%d-%m-%Y"),
    (MDY_RE, "%m-%d-%Y"),
)
TIME_FORMATS = ((HMS_RE, "%H:%M:%S"), (HM_RE, "%H:%M"), (HMS_RE_NO_SEP, "%H%M%S"))


def _extract_regex_concat_arrays(
    strings: ChunkedArrayAny,
    /,
    pattern: str,
    *,
    options: Any = None,
    memory_pool: Any = None,
) -> pa.StructArray:
    r = pa.concat_arrays(
        extract_regex(strings, pattern, options=options, memory_pool=memory_pool).chunks
    )
    return cast("pa.StructArray", r)


def parse_datetime_format(arr: ChunkedArrayAny) -> str:
    """Try to infer datetime format from StringArray."""
    matches = _extract_regex_concat_arrays(arr.drop_null().slice(0, 10), pattern=FULL_RE)
    if not pc.all(matches.is_valid()).as_py():
        msg = (
            "Unable to infer datetime format, provided format is not supported. "
            "Please report a bug to https://github.com/narwhals-dev/narwhals/issues"
        )
        raise NotImplementedError(msg)

    separators = matches.field("sep")
    tz = matches.field("tz")

    # separators and time zones must be unique
    if pc.count(pc.unique(separators)).as_py() > 1:
        msg = "Found multiple separator values while inferring datetime format."
        raise ValueError(msg)

    if pc.count(pc.unique(tz)).as_py() > 1:
        msg = "Found multiple timezone values while inferring datetime format."
        raise ValueError(msg)

    date_value = _parse_date_format(cast("pc.StringArray", matches.field("date")))
    time_value = _parse_time_format(cast("pc.StringArray", matches.field("time")))

    sep_value = separators[0].as_py()
    tz_value = "%z" if tz[0].as_py() else ""

    return f"{date_value}{sep_value}{time_value}{tz_value}"


def _parse_date_format(arr: pc.StringArray) -> str:
    for date_rgx, date_fmt in DATE_FORMATS:
        matches = pc.extract_regex(arr, pattern=date_rgx)
        if date_fmt == "%Y%m%d" and pc.all(matches.is_valid()).as_py():
            return date_fmt
        if (
            pc.all(matches.is_valid()).as_py()
            and pc.count(pc.unique(sep1 := matches.field("sep1"))).as_py() == 1
            and pc.count(pc.unique(sep2 := matches.field("sep2"))).as_py() == 1
            and (date_sep_value := sep1[0].as_py()) == sep2[0].as_py()
        ):
            return date_fmt.replace("-", date_sep_value)

    msg = (
        "Unable to infer datetime format. "
        "Please report a bug to https://github.com/narwhals-dev/narwhals/issues"
    )
    raise ValueError(msg)


def _parse_time_format(arr: pc.StringArray) -> str:
    for time_rgx, time_fmt in TIME_FORMATS:
        matches = pc.extract_regex(arr, pattern=time_rgx)
        if pc.all(matches.is_valid()).as_py():
            return time_fmt
    return ""


def pad_series(
    series: ArrowSeries, *, window_size: int, center: bool
) -> tuple[ArrowSeries, int]:
    """Pad series with None values on the left and/or right side, depending on the specified parameters.

    Arguments:
        series: The input ArrowSeries to be padded.
        window_size: The desired size of the window.
        center: Specifies whether to center the padding or not.

    Returns:
        A tuple containing the padded ArrowSeries and the offset value.
    """
    if not center:
        return series, 0
    offset_left = window_size // 2
    # subtract one if window_size is even
    offset_right = offset_left - (window_size % 2 == 0)
    chunks = series.native.chunks
    arrays = nulls_like(offset_left, series), *chunks, nulls_like(offset_right, series)
    return series._with_native(pa.concat_arrays(arrays)), offset_left + offset_right


def cast_to_comparable_string_types(
    *chunked_arrays: ChunkedArrayAny | ScalarAny, separator: str
) -> tuple[Iterator[ChunkedArrayAny | ScalarAny], ScalarAny]:
    # Ensure `chunked_arrays` are either all `string` or all `large_string`.
    dtype = (
        pa.string()  # (PyArrow default)
        if not any(pa.types.is_large_string(ca.type) for ca in chunked_arrays)
        else pa.large_string()
    )
    return (ca.cast(dtype) for ca in chunked_arrays), lit(separator, dtype)


if BACKEND_VERSION >= (14,):
    # https://arrow.apache.org/docs/14.0/python/generated/pyarrow.concat_tables.html
    _PROMOTE: Mapping[PromoteOptions, Mapping[str, Any]] = {
        "default": {"promote_options": "default"},
        "permissive": {"promote_options": "permissive"},
        "none": {"promote_options": "none"},
    }
else:  # pragma: no cover
    # https://arrow.apache.org/docs/13.0/python/generated/pyarrow.concat_tables.html
    _PROMOTE = {
        "default": {"promote": True},
        "permissive": {"promote": True},
        "none": {"promote": False},
    }


def concat_tables(
    tables: Iterable[pa.Table], promote_options: PromoteOptions = "none"
) -> pa.Table:
    return pa.concat_tables(tables, **_PROMOTE[promote_options])


class ArrowSeriesNamespace(EagerSeriesNamespace["ArrowSeries", "ChunkedArrayAny"]): ...
