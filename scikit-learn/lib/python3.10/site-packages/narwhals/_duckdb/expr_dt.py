from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._constants import (
    MS_PER_MINUTE,
    MS_PER_SECOND,
    NS_PER_SECOND,
    SECONDS_PER_MINUTE,
    US_PER_MINUTE,
    US_PER_SECOND,
)
from narwhals._duckdb.utils import UNITS_DICT, F, fetch_rel_time_zone, lit
from narwhals._duration import Interval
from narwhals._sql.expr_dt import SQLExprDateTimeNamesSpace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    from collections.abc import Sequence

    from duckdb import Expression

    from narwhals._duckdb.dataframe import DuckDBLazyFrame
    from narwhals._duckdb.expr import DuckDBExpr


class DuckDBExprDateTimeNamespace(SQLExprDateTimeNamesSpace["DuckDBExpr"]):
    def millisecond(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("millisecond", expr) - F("second", expr) * lit(MS_PER_SECOND)
        )

    def microsecond(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("microsecond", expr) - F("second", expr) * lit(US_PER_SECOND)
        )

    def nanosecond(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("nanosecond", expr) - F("second", expr) * lit(NS_PER_SECOND)
        )

    def to_string(self, format: str) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("strftime", expr, lit(format))
        )

    def weekday(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("isodow", expr))

    def date(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: expr.cast("date"))

    def total_minutes(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("datepart", lit("minute"), expr)
        )

    def total_seconds(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: lit(SECONDS_PER_MINUTE) * F("datepart", lit("minute"), expr)
            + F("datepart", lit("second"), expr)
        )

    def total_milliseconds(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: lit(MS_PER_MINUTE) * F("datepart", lit("minute"), expr)
            + F("datepart", lit("millisecond"), expr)
        )

    def total_microseconds(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: lit(US_PER_MINUTE) * F("datepart", lit("minute"), expr)
            + F("datepart", lit("microsecond"), expr)
        )

    def truncate(self, every: str) -> DuckDBExpr:
        interval = Interval.parse(every)
        multiple, unit = interval.multiple, interval.unit
        if multiple != 1:
            # https://github.com/duckdb/duckdb/issues/17554
            msg = f"Only multiple 1 is currently supported for DuckDB.\nGot {multiple!s}."
            raise ValueError(msg)
        if unit == "ns":
            msg = "Truncating to nanoseconds is not yet supported for DuckDB."
            raise NotImplementedError(msg)
        format = lit(UNITS_DICT[unit])

        def _truncate(expr: Expression) -> Expression:
            return F("date_trunc", format, expr)

        return self.compliant._with_elementwise(_truncate)

    def offset_by(self, by: str) -> DuckDBExpr:
        interval = Interval.parse_no_constraints(by)
        format = lit(f"{interval.multiple!s} {UNITS_DICT[interval.unit]}")

        def _offset_by(expr: Expression) -> Expression:
            return F("date_add", format, expr)

        return self.compliant._with_callable(_offset_by)

    def _no_op_time_zone(self, time_zone: str) -> DuckDBExpr:
        def func(df: DuckDBLazyFrame) -> Sequence[Expression]:
            native_series_list = self.compliant(df)
            conn_time_zone = fetch_rel_time_zone(df.native)
            if conn_time_zone != time_zone:
                msg = (
                    "DuckDB stores the time zone in the connection, rather than in the "
                    f"data type, so changing the timezone to anything other than {conn_time_zone} "
                    " (the current connection time zone) is not supported."
                )
                raise NotImplementedError(msg)
            return native_series_list

        return self.compliant.__class__(
            func,
            evaluate_output_names=self.compliant._evaluate_output_names,
            alias_output_names=self.compliant._alias_output_names,
            version=self.compliant._version,
        )

    def convert_time_zone(self, time_zone: str) -> DuckDBExpr:
        return self._no_op_time_zone(time_zone)

    def replace_time_zone(self, time_zone: str | None) -> DuckDBExpr:
        if time_zone is None:
            return self.compliant._with_elementwise(lambda expr: expr.cast("timestamp"))
        return self._no_op_time_zone(time_zone)

    total_nanoseconds = not_implemented()
    timestamp = not_implemented()
