from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import duckdb
from duckdb import Expression

try:
    import duckdb.sqltypes as duckdb_dtypes
except ModuleNotFoundError:
    # DuckDB pre 1.3
    import duckdb.typing as duckdb_dtypes

from narwhals._utils import Version, extend_bool, isinstance_or_issubclass, zip_strict
from narwhals.exceptions import ColumnNotFoundError

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from duckdb import DuckDBPyRelation

    from narwhals._compliant.typing import CompliantLazyFrameAny
    from narwhals._duckdb.dataframe import DuckDBLazyFrame
    from narwhals._duckdb.expr import DuckDBExpr
    from narwhals.dtypes import DType
    from narwhals.typing import IntoDType, TimeUnit


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
DESCENDING_TO_ORDER = {True: "desc", False: "asc"}
NULLS_LAST_TO_NULLS_POS = {True: "nulls last", False: "nulls first"}

col = duckdb.ColumnExpression
"""Alias for `duckdb.ColumnExpression`."""

lit = duckdb.ConstantExpression
"""Alias for `duckdb.ConstantExpression`."""

when = duckdb.CaseExpression
"""Alias for `duckdb.CaseExpression`."""

F = duckdb.FunctionExpression
"""Alias for `duckdb.FunctionExpression`."""


def lambda_expr(
    params: str | Expression | tuple[Expression, ...], expr: Expression, /
) -> Expression:
    """Wraps [`duckdb.LambdaExpression`].

    [`duckdb.LambdaExpression`]: https://duckdb.org/docs/stable/sql/functions/lambda
    """
    try:
        from duckdb import LambdaExpression
    except ModuleNotFoundError as exc:
        msg = f"DuckDB>=1.2.0 is required for this operation. Found: DuckDB {duckdb.__version__}"
        raise NotImplementedError(msg) from exc
    args = (params,) if isinstance(params, Expression) else params
    return LambdaExpression(args, expr)


def concat_str(*exprs: Expression, separator: str = "") -> Expression:
    """Concatenate many strings, NULL inputs are skipped.

    Wraps [concat] and [concat_ws] `FunctionExpression`(s).

    Arguments:
        exprs: Native columns.
        separator: String that will be used to separate the values of each column.


    [concat]: https://duckdb.org/docs/stable/sql/functions/char.html#concatstring-
    [concat_ws]: https://duckdb.org/docs/stable/sql/functions/char.html#concat_wsseparator-string-
    """
    return F("concat_ws", lit(separator), *exprs) if separator else F("concat", *exprs)


def evaluate_exprs(
    df: DuckDBLazyFrame, /, *exprs: DuckDBExpr
) -> list[tuple[str, Expression]]:
    native_results: list[tuple[str, Expression]] = []
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


class DeferredTimeZone:
    """Object which gets passed between `native_to_narwhals_dtype` calls.

    DuckDB stores the time zone in the connection, rather than in the dtypes, so
    this ensures that when calculating the schema of a dataframe with multiple
    timezone-aware columns, that the connection's time zone is only fetched once.

    Note: we cannot make the time zone a cached `DuckDBLazyFrame` property because
    the time zone can be modified after `DuckDBLazyFrame` creation:

    ```python
    df = nw.from_native(rel)
    print(df.collect_schema())
    rel.query("set timezone = 'Asia/Kolkata'")
    print(df.collect_schema())  # should change to reflect new time zone
    ```
    """

    _cached_time_zone: str | None = None

    def __init__(self, rel: DuckDBPyRelation) -> None:
        self._rel = rel

    @property
    def time_zone(self) -> str:
        """Fetch relation time zone (if it wasn't calculated already)."""
        if self._cached_time_zone is None:
            self._cached_time_zone = fetch_rel_time_zone(self._rel)
        return self._cached_time_zone


def native_to_narwhals_dtype(
    duckdb_dtype: duckdb_dtypes.DuckDBPyType,
    version: Version,
    deferred_time_zone: DeferredTimeZone,
) -> DType:
    duckdb_dtype_id = duckdb_dtype.id
    dtypes = version.dtypes

    # Handle nested data types first
    if duckdb_dtype_id == "list":
        return dtypes.List(
            native_to_narwhals_dtype(duckdb_dtype.child, version, deferred_time_zone)
        )

    if duckdb_dtype_id == "struct":
        children = duckdb_dtype.children
        return dtypes.Struct(
            [
                dtypes.Field(
                    name=child[0],
                    dtype=native_to_narwhals_dtype(child[1], version, deferred_time_zone),
                )
                for child in children
            ]
        )

    if duckdb_dtype_id == "array":
        child, size = duckdb_dtype.children
        shape: list[int] = [size[1]]

        while child[1].id == "array":
            child, size = child[1].children
            shape.insert(0, size[1])

        inner = native_to_narwhals_dtype(child[1], version, deferred_time_zone)
        return dtypes.Array(inner=inner, shape=tuple(shape))

    if duckdb_dtype_id == "enum":
        if version is Version.V1:
            return dtypes.Enum()  # type: ignore[call-arg]
        categories = duckdb_dtype.children[0][1]
        return dtypes.Enum(categories=categories)

    if duckdb_dtype_id == "timestamp with time zone":
        return dtypes.Datetime(time_zone=deferred_time_zone.time_zone)

    return _non_nested_native_to_narwhals_dtype(duckdb_dtype_id, version)


def fetch_rel_time_zone(rel: duckdb.DuckDBPyRelation) -> str:
    result = rel.query(
        "duckdb_settings()", "select value from duckdb_settings() where name = 'TimeZone'"
    ).fetchone()
    assert result is not None  # noqa: S101
    return result[0]  # type: ignore[no-any-return]


@lru_cache(maxsize=16)
def _non_nested_native_to_narwhals_dtype(duckdb_dtype_id: str, version: Version) -> DType:
    dtypes = version.dtypes
    return {
        "hugeint": dtypes.Int128(),
        "bigint": dtypes.Int64(),
        "integer": dtypes.Int32(),
        "smallint": dtypes.Int16(),
        "tinyint": dtypes.Int8(),
        "uhugeint": dtypes.UInt128(),
        "ubigint": dtypes.UInt64(),
        "uinteger": dtypes.UInt32(),
        "usmallint": dtypes.UInt16(),
        "utinyint": dtypes.UInt8(),
        "double": dtypes.Float64(),
        "float": dtypes.Float32(),
        "varchar": dtypes.String(),
        "date": dtypes.Date(),
        "timestamp_s": dtypes.Datetime("s"),
        "timestamp_ms": dtypes.Datetime("ms"),
        "timestamp": dtypes.Datetime(),
        "timestamp_ns": dtypes.Datetime("ns"),
        "boolean": dtypes.Boolean(),
        "interval": dtypes.Duration(),
        "decimal": dtypes.Decimal(),
        "time": dtypes.Time(),
        "blob": dtypes.Binary(),
    }.get(duckdb_dtype_id, dtypes.Unknown())


dtypes = Version.MAIN.dtypes
NW_TO_DUCKDB_DTYPES: Mapping[type[DType], duckdb_dtypes.DuckDBPyType] = {
    dtypes.Float64: duckdb_dtypes.DOUBLE,
    dtypes.Float32: duckdb_dtypes.FLOAT,
    dtypes.Binary: duckdb_dtypes.BLOB,
    dtypes.String: duckdb_dtypes.VARCHAR,
    dtypes.Boolean: duckdb_dtypes.BOOLEAN,
    dtypes.Date: duckdb_dtypes.DATE,
    dtypes.Time: duckdb_dtypes.TIME,
    dtypes.Int8: duckdb_dtypes.TINYINT,
    dtypes.Int16: duckdb_dtypes.SMALLINT,
    dtypes.Int32: duckdb_dtypes.INTEGER,
    dtypes.Int64: duckdb_dtypes.BIGINT,
    dtypes.Int128: duckdb_dtypes.HUGEINT,
    dtypes.UInt8: duckdb_dtypes.UTINYINT,
    dtypes.UInt16: duckdb_dtypes.USMALLINT,
    dtypes.UInt32: duckdb_dtypes.UINTEGER,
    dtypes.UInt64: duckdb_dtypes.UBIGINT,
    dtypes.UInt128: duckdb_dtypes.UHUGEINT,
}
TIME_UNIT_TO_TIMESTAMP: Mapping[TimeUnit, duckdb_dtypes.DuckDBPyType] = {
    "s": duckdb_dtypes.TIMESTAMP_S,
    "ms": duckdb_dtypes.TIMESTAMP_MS,
    "us": duckdb_dtypes.TIMESTAMP,
    "ns": duckdb_dtypes.TIMESTAMP_NS,
}
UNSUPPORTED_DTYPES = (dtypes.Decimal, dtypes.Categorical)


def narwhals_to_native_dtype(  # noqa: PLR0912, C901
    dtype: IntoDType, version: Version, deferred_time_zone: DeferredTimeZone
) -> duckdb_dtypes.DuckDBPyType:
    dtypes = version.dtypes
    base_type = dtype.base_type()
    if duckdb_type := NW_TO_DUCKDB_DTYPES.get(base_type):
        return duckdb_type
    if isinstance_or_issubclass(dtype, dtypes.Enum):
        if version is Version.V1:
            msg = "Converting to Enum is not supported in narwhals.stable.v1"
            raise NotImplementedError(msg)
        if isinstance(dtype, dtypes.Enum):
            return duckdb_dtypes.DuckDBPyType(f"ENUM{dtype.categories!r}")
        msg = "Can not cast / initialize Enum without categories present"
        raise ValueError(msg)
    if isinstance_or_issubclass(dtype, dtypes.Datetime):
        tu = dtype.time_unit
        tz = dtype.time_zone
        if not tz:
            return TIME_UNIT_TO_TIMESTAMP[tu]
        if tu != "us":
            msg = f"Only microsecond precision is supported for timezone-aware `Datetime` in DuckDB, got {tu} precision"
            raise ValueError(msg)
        if tz != (rel_tz := deferred_time_zone.time_zone):  # pragma: no cover
            msg = f"Only the connection time zone {rel_tz} is supported, got: {tz}."
            raise ValueError(msg)
        # TODO(unassigned): cover once https://github.com/narwhals-dev/narwhals/issues/2742 addressed
        return duckdb_dtypes.TIMESTAMP_TZ  # pragma: no cover
    if isinstance_or_issubclass(dtype, dtypes.Duration):
        if (tu := dtype.time_unit) != "us":  # pragma: no cover
            msg = f"Only microsecond-precision Duration is supported, got {tu} precision"
        return duckdb_dtypes.INTERVAL
    if isinstance_or_issubclass(dtype, dtypes.List):
        inner = narwhals_to_native_dtype(dtype.inner, version, deferred_time_zone)
        return duckdb.list_type(inner)
    if isinstance_or_issubclass(dtype, dtypes.Struct):
        fields = {
            field.name: narwhals_to_native_dtype(field.dtype, version, deferred_time_zone)
            for field in dtype.fields
        }
        return duckdb.struct_type(fields)
    if isinstance(dtype, dtypes.Array):
        nw_inner: IntoDType = dtype
        while isinstance(nw_inner, dtypes.Array):
            nw_inner = nw_inner.inner
        duckdb_inner = narwhals_to_native_dtype(nw_inner, version, deferred_time_zone)
        duckdb_shape_fmt = "".join(f"[{item}]" for item in dtype.shape)
        return duckdb_dtypes.DuckDBPyType(f"{duckdb_inner}{duckdb_shape_fmt}")
    if issubclass(base_type, UNSUPPORTED_DTYPES):
        msg = f"Converting to {base_type.__name__} dtype is not supported for DuckDB."
        raise NotImplementedError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)


def parse_into_expression(into_expression: str | Expression) -> Expression:
    return col(into_expression) if isinstance(into_expression, str) else into_expression


def generate_partition_by_sql(*partition_by: str | Expression) -> str:
    if not partition_by:
        return ""
    by_sql = ", ".join([f"{parse_into_expression(x)}" for x in partition_by])
    return f"partition by {by_sql}"


def join_column_names(*names: str) -> str:
    return ", ".join(str(col(name)) for name in names)


def generate_order_by_sql(
    *order_by: str | Expression, descending: Sequence[bool], nulls_last: Sequence[bool]
) -> str:
    if not order_by:
        return ""
    by_sql = ",".join(
        f"{parse_into_expression(x)} {DESCENDING_TO_ORDER[_descending]} {NULLS_LAST_TO_NULLS_POS[_nulls_last]}"
        for x, _descending, _nulls_last in zip_strict(order_by, descending, nulls_last)
    )
    return f"order by {by_sql}"


def window_expression(
    expr: Expression,
    partition_by: Sequence[str | Expression] = (),
    order_by: Sequence[str | Expression] = (),
    rows_start: int | None = None,
    rows_end: int | None = None,
    *,
    descending: Sequence[bool] | None = None,
    nulls_last: Sequence[bool] | None = None,
    ignore_nulls: bool = False,
) -> Expression:
    # TODO(unassigned): Replace with `duckdb.WindowExpression` when they release it.
    # https://github.com/duckdb/duckdb/discussions/14725#discussioncomment-11200348
    pb = generate_partition_by_sql(*partition_by)
    flags = extend_bool(False, len(order_by))
    descending = descending or flags
    nulls_last = nulls_last or flags
    ob = generate_order_by_sql(*order_by, descending=descending, nulls_last=nulls_last)

    if rows_start is not None and rows_end is not None:
        rows = f"rows between {-rows_start} preceding and {rows_end} following"
    elif rows_end is not None:
        rows = f"rows between unbounded preceding and {rows_end} following"
    elif rows_start is not None:
        rows = f"rows between {-rows_start} preceding and unbounded following"
    else:
        rows = ""

    func = f"{str(expr).removesuffix(')')} ignore nulls)" if ignore_nulls else str(expr)
    return sql_expression(f"{func} over ({pb} {ob} {rows})")


def catch_duckdb_exception(
    exception: Exception, frame: CompliantLazyFrameAny, /
) -> ColumnNotFoundError | Exception:
    if isinstance(exception, duckdb.BinderException) and any(
        msg in str(exception)
        for msg in (
            "not found in FROM clause",
            "this column cannot be referenced before it is defined",
        )
    ):
        return ColumnNotFoundError.from_available_column_names(
            available_columns=frame.columns
        )
    # Just return exception as-is.
    return exception


def function(name: str, *args: Expression) -> Expression:
    if name == "isnull":
        return args[0].isnull()
    if name == "count_distinct":
        return sql_expression(f"count(distinct {args[0]})")
    return F(name, *args)


def sql_expression(expr: str) -> Expression:
    try:
        from duckdb import SQLExpression
    except ImportError as exc:  # pragma: no cover
        msg = f"DuckDB>=1.3.0 is required for this operation. Found: DuckDB {duckdb.__version__}"
        raise NotImplementedError(msg) from exc
    return SQLExpression(expr)


__all__ = [
    "UNITS_DICT",
    "DeferredTimeZone",
    "F",
    "catch_duckdb_exception",
    "col",
    "concat_str",
    "duckdb_dtypes",
    "evaluate_exprs",
    "fetch_rel_time_zone",
    "function",
    "generate_order_by_sql",
    "generate_partition_by_sql",
    "join_column_names",
    "lambda_expr",
    "lit",
    "narwhals_to_native_dtype",
    "native_to_narwhals_dtype",
    "parse_into_expression",
    "sql_expression",
    "when",
    "window_expression",
]
