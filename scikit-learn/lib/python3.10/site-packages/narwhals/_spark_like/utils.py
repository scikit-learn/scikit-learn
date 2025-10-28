from __future__ import annotations

import operator
from collections.abc import Callable
from functools import lru_cache
from importlib import import_module
from operator import attrgetter
from types import ModuleType
from typing import TYPE_CHECKING, Any, overload

from narwhals._utils import Implementation, Version, isinstance_or_issubclass
from narwhals.exceptions import ColumnNotFoundError, UnsupportedDTypeError

if TYPE_CHECKING:
    from collections.abc import Mapping

    import sqlframe.base.types as sqlframe_types
    from sqlframe.base.column import Column
    from sqlframe.base.session import _BaseSession as Session
    from typing_extensions import TypeAlias

    from narwhals._compliant.typing import CompliantLazyFrameAny
    from narwhals._spark_like.dataframe import SparkLikeLazyFrame
    from narwhals._spark_like.expr import SparkLikeExpr
    from narwhals.dtypes import DType
    from narwhals.typing import IntoDType

    _NativeDType: TypeAlias = sqlframe_types.DataType
    SparkSession = Session[Any, Any, Any, Any, Any, Any, Any]

UNITS_DICT = {
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

# see https://spark.apache.org/docs/latest/sql-ref-datetime-pattern.html
# and https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior
DATETIME_PATTERNS_MAPPING = {
    "%Y": "yyyy",  # Year with century (4 digits)
    "%y": "yy",  # Year without century (2 digits)
    "%m": "MM",  # Month (01-12)
    "%d": "dd",  # Day of the month (01-31)
    "%H": "HH",  # Hour (24-hour clock) (00-23)
    "%I": "hh",  # Hour (12-hour clock) (01-12)
    "%M": "mm",  # Minute (00-59)
    "%S": "ss",  # Second (00-59)
    "%f": "S",  # Microseconds -> Milliseconds
    "%p": "a",  # AM/PM
    "%a": "E",  # Abbreviated weekday name
    "%A": "E",  # Full weekday name
    "%j": "D",  # Day of the year
    "%z": "Z",  # Timezone offset
    "%s": "X",  # Unix timestamp
}


# NOTE: don't lru_cache this as `ModuleType` isn't hashable
def native_to_narwhals_dtype(  # noqa: C901, PLR0912
    dtype: _NativeDType, version: Version, spark_types: ModuleType, session: SparkSession
) -> DType:
    dtypes = version.dtypes
    if TYPE_CHECKING:
        native = sqlframe_types
    else:
        native = spark_types

    if isinstance(dtype, native.DoubleType):
        return dtypes.Float64()
    if isinstance(dtype, native.FloatType):
        return dtypes.Float32()
    if isinstance(dtype, native.LongType):
        return dtypes.Int64()
    if isinstance(dtype, native.IntegerType):
        return dtypes.Int32()
    if isinstance(dtype, native.ShortType):
        return dtypes.Int16()
    if isinstance(dtype, native.ByteType):
        return dtypes.Int8()
    if isinstance(dtype, (native.StringType, native.VarcharType, native.CharType)):
        return dtypes.String()
    if isinstance(dtype, native.BooleanType):
        return dtypes.Boolean()
    if isinstance(dtype, native.DateType):
        return dtypes.Date()
    if isinstance(dtype, native.TimestampNTZType):
        # TODO(marco): cover this
        return dtypes.Datetime()  # pragma: no cover
    if isinstance(dtype, native.TimestampType):
        return dtypes.Datetime(time_zone=fetch_session_time_zone(session))
    if isinstance(dtype, native.DecimalType):
        # TODO(marco): cover this
        return dtypes.Decimal()  # pragma: no cover
    if isinstance(dtype, native.ArrayType):
        return dtypes.List(
            inner=native_to_narwhals_dtype(
                dtype.elementType, version, spark_types, session
            )
        )
    if isinstance(dtype, native.StructType):
        return dtypes.Struct(
            fields=[
                dtypes.Field(
                    name=field.name,
                    dtype=native_to_narwhals_dtype(
                        field.dataType, version, spark_types, session
                    ),
                )
                for field in dtype
            ]
        )
    if isinstance(dtype, native.BinaryType):
        return dtypes.Binary()
    return dtypes.Unknown()  # pragma: no cover


@lru_cache(maxsize=4)
def fetch_session_time_zone(session: SparkSession) -> str:
    # Timezone can't be changed in PySpark session, so this can be cached.
    try:
        return session.conf.get("spark.sql.session.timeZone")  # type: ignore[attr-defined]
    except Exception:  # noqa: BLE001
        # https://github.com/eakmanrq/sqlframe/issues/406
        return "<unknown>"


IntoSparkDType: TypeAlias = Callable[[ModuleType], Callable[[], "_NativeDType"]]
dtypes = Version.MAIN.dtypes
NW_TO_SPARK_DTYPES: Mapping[type[DType], IntoSparkDType] = {
    dtypes.Float64: attrgetter("DoubleType"),
    dtypes.Float32: attrgetter("FloatType"),
    dtypes.Binary: attrgetter("BinaryType"),
    dtypes.String: attrgetter("StringType"),
    dtypes.Boolean: attrgetter("BooleanType"),
    dtypes.Date: attrgetter("DateType"),
    dtypes.Int8: attrgetter("ByteType"),
    dtypes.Int16: attrgetter("ShortType"),
    dtypes.Int32: attrgetter("IntegerType"),
    dtypes.Int64: attrgetter("LongType"),
}
UNSUPPORTED_DTYPES = (
    dtypes.UInt64,
    dtypes.UInt32,
    dtypes.UInt16,
    dtypes.UInt8,
    dtypes.Enum,
    dtypes.Categorical,
    dtypes.Time,
)


def narwhals_to_native_dtype(
    dtype: IntoDType, version: Version, spark_types: ModuleType, session: SparkSession
) -> _NativeDType:
    dtypes = version.dtypes
    if TYPE_CHECKING:
        native = sqlframe_types
    else:
        native = spark_types
    base_type = dtype.base_type()
    if into_spark_type := NW_TO_SPARK_DTYPES.get(base_type):
        return into_spark_type(native)()
    if isinstance_or_issubclass(dtype, dtypes.Datetime):
        if (tu := dtype.time_unit) != "us":  # pragma: no cover
            msg = f"Only microsecond precision is supported for PySpark, got: {tu}."
            raise ValueError(msg)
        dt_time_zone = dtype.time_zone
        if dt_time_zone is None:
            return native.TimestampNTZType()
        if dt_time_zone != (tz := fetch_session_time_zone(session)):  # pragma: no cover
            msg = f"Only {tz} time zone is supported, as that's the connection time zone, got: {dt_time_zone}"
            raise ValueError(msg)
        # TODO(unassigned): cover once https://github.com/narwhals-dev/narwhals/issues/2742 addressed
        return native.TimestampType()  # pragma: no cover
    if isinstance_or_issubclass(dtype, (dtypes.List, dtypes.Array)):
        return native.ArrayType(
            elementType=narwhals_to_native_dtype(dtype.inner, version, native, session)
        )
    if isinstance_or_issubclass(dtype, dtypes.Struct):  # pragma: no cover
        return native.StructType(
            fields=[
                native.StructField(
                    name=field.name,
                    dataType=narwhals_to_native_dtype(
                        field.dtype, version, native, session
                    ),
                )
                for field in dtype.fields
            ]
        )
    if issubclass(base_type, UNSUPPORTED_DTYPES):  # pragma: no cover
        msg = f"Converting to {base_type.__name__} dtype is not supported for Spark-Like backend."
        raise UnsupportedDTypeError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)


def evaluate_exprs(
    df: SparkLikeLazyFrame, /, *exprs: SparkLikeExpr
) -> list[tuple[str, Column]]:
    native_results: list[tuple[str, Column]] = []

    for expr in exprs:
        native_series_list = expr._call(df)
        output_names = expr._evaluate_output_names(df)
        if expr._alias_output_names is not None:
            output_names = expr._alias_output_names(output_names)
        if len(output_names) != len(native_series_list):  # pragma: no cover
            msg = f"Internal error: got output names {output_names}, but only got {len(native_series_list)} results"
            raise AssertionError(msg)
        native_results.extend(zip(output_names, native_series_list))

    return native_results


def import_functions(implementation: Implementation, /) -> ModuleType:
    if implementation is Implementation.PYSPARK:
        from pyspark.sql import functions

        return functions
    if implementation is Implementation.PYSPARK_CONNECT:
        from pyspark.sql.connect import functions

        return functions
    from sqlframe.base.session import _BaseSession

    return import_module(f"sqlframe.{_BaseSession().execution_dialect_name}.functions")


def import_native_dtypes(implementation: Implementation, /) -> ModuleType:
    if implementation is Implementation.PYSPARK:
        from pyspark.sql import types

        return types
    if implementation is Implementation.PYSPARK_CONNECT:
        from pyspark.sql.connect import types

        return types
    from sqlframe.base.session import _BaseSession

    return import_module(f"sqlframe.{_BaseSession().execution_dialect_name}.types")


def import_window(implementation: Implementation, /) -> type[Any]:
    if implementation is Implementation.PYSPARK:
        from pyspark.sql import Window

        return Window

    if implementation is Implementation.PYSPARK_CONNECT:
        from pyspark.sql.connect.window import Window

        return Window
    from sqlframe.base.session import _BaseSession

    return import_module(
        f"sqlframe.{_BaseSession().execution_dialect_name}.window"
    ).Window


@overload
def strptime_to_pyspark_format(format: None) -> None: ...


@overload
def strptime_to_pyspark_format(format: str) -> str: ...


def strptime_to_pyspark_format(format: str | None) -> str | None:
    """Converts a Python strptime datetime format string to a PySpark datetime format string."""
    if format is None:  # pragma: no cover
        return None

    # Replace Python format specifiers with PySpark specifiers
    pyspark_format = format
    for py_format, spark_format in DATETIME_PATTERNS_MAPPING.items():
        pyspark_format = pyspark_format.replace(py_format, spark_format)
    return pyspark_format.replace("T", " ")


def true_divide(F: ModuleType, left: Column, right: Column) -> Column:
    # PySpark before 3.5 doesn't have `try_divide`, SQLFrame doesn't have it.
    divide = getattr(F, "try_divide", operator.truediv)
    zero = F.lit(0)

    safe_case = divide(left, right)  # Dividend is not zero
    unsafe_case = (  # Dividend is zero, the result depends on the numerator
        F.when(left == zero, F.lit(float("nan")))
        .when(left > zero, F.lit(float("inf")))
        .when(left < zero, F.lit(float("-inf")))
    )
    return F.when(right != zero, safe_case).when(right == zero, unsafe_case)


def catch_pyspark_sql_exception(
    exception: Exception, frame: CompliantLazyFrameAny, /
) -> ColumnNotFoundError | Exception:  # pragma: no cover
    from pyspark.errors import AnalysisException

    if isinstance(exception, AnalysisException) and str(exception).startswith(
        "[UNRESOLVED_COLUMN.WITH_SUGGESTION]"
    ):
        return ColumnNotFoundError.from_available_column_names(
            available_columns=frame.columns
        )
    # Just return exception as-is.
    return exception


def catch_pyspark_connect_exception(
    exception: Exception, /
) -> ColumnNotFoundError | Exception:  # pragma: no cover
    from pyspark.errors.exceptions.connect import AnalysisException

    if isinstance(exception, AnalysisException) and str(exception).startswith(
        "[UNRESOLVED_COLUMN.WITH_SUGGESTION]"
    ):
        return ColumnNotFoundError(str(exception))
    # Just return exception as-is.
    return exception
