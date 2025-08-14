from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._constants import US_PER_SECOND
from narwhals._duration import Interval
from narwhals._spark_like.utils import (
    UNITS_DICT,
    fetch_session_time_zone,
    strptime_to_pyspark_format,
)
from narwhals._sql.expr_dt import SQLExprDateTimeNamesSpace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    from collections.abc import Sequence

    from sqlframe.base.column import Column

    from narwhals._spark_like.dataframe import SparkLikeLazyFrame
    from narwhals._spark_like.expr import SparkLikeExpr


class SparkLikeExprDateTimeNamespace(SQLExprDateTimeNamesSpace["SparkLikeExpr"]):
    def _weekday(self, expr: Column) -> Column:
        # PySpark's dayofweek returns 1-7 for Sunday-Saturday
        return (self.compliant._F.dayofweek(expr) + 6) % 7

    def to_string(self, format: str) -> SparkLikeExpr:
        F = self.compliant._F

        def _to_string(expr: Column) -> Column:
            # Handle special formats
            if format == "%G-W%V":
                return self._format_iso_week(expr)
            if format == "%G-W%V-%u":
                return self._format_iso_week_with_day(expr)

            format_, suffix = self._format_microseconds(expr, format)

            # Convert Python format to PySpark format
            pyspark_fmt = strptime_to_pyspark_format(format_)

            result = F.date_format(expr, pyspark_fmt)
            if "T" in format_:
                # `strptime_to_pyspark_format` replaces "T" with " " since pyspark
                # does not support the literal "T" in `date_format`.
                # If no other spaces are in the given format, then we can revert this
                # operation, otherwise we raise an exception.
                if " " not in format_:
                    result = F.replace(result, F.lit(" "), F.lit("T"))
                else:  # pragma: no cover
                    msg = (
                        "`dt.to_string` with a format that contains both spaces and "
                        " the literal 'T' is not supported for spark-like backends."
                    )
                    raise NotImplementedError(msg)

            return F.concat(result, *suffix)

        return self.compliant._with_elementwise(_to_string)

    def millisecond(self) -> SparkLikeExpr:
        def _millisecond(expr: Column) -> Column:
            return self.compliant._F.floor(
                (self.compliant._F.unix_micros(expr) % US_PER_SECOND) / 1000
            )

        return self.compliant._with_elementwise(_millisecond)

    def microsecond(self) -> SparkLikeExpr:
        def _microsecond(expr: Column) -> Column:
            return self.compliant._F.unix_micros(expr) % US_PER_SECOND

        return self.compliant._with_elementwise(_microsecond)

    def nanosecond(self) -> SparkLikeExpr:
        def _nanosecond(expr: Column) -> Column:
            return (self.compliant._F.unix_micros(expr) % US_PER_SECOND) * 1000

        return self.compliant._with_elementwise(_nanosecond)

    def weekday(self) -> SparkLikeExpr:
        return self.compliant._with_elementwise(self._weekday)

    def truncate(self, every: str) -> SparkLikeExpr:
        interval = Interval.parse(every)
        multiple, unit = interval.multiple, interval.unit
        if multiple != 1:
            msg = f"Only multiple 1 is currently supported for Spark-like.\nGot {multiple!s}."
            raise ValueError(msg)
        if unit == "ns":
            msg = "Truncating to nanoseconds is not yet supported for Spark-like."
            raise NotImplementedError(msg)
        format = UNITS_DICT[unit]

        def _truncate(expr: Column) -> Column:
            return self.compliant._F.date_trunc(format, expr)

        return self.compliant._with_elementwise(_truncate)

    def offset_by(self, by: str) -> SparkLikeExpr:
        interval = Interval.parse_no_constraints(by)
        multiple, unit = interval.multiple, interval.unit
        if unit == "ns":  # pragma: no cover
            msg = "Offsetting by nanoseconds is not yet supported for Spark-like."
            raise NotImplementedError(msg)

        F = self.compliant._F

        def _offset_by(expr: Column) -> Column:
            # https://github.com/eakmanrq/sqlframe/issues/441
            return F.timestamp_add(  # pyright: ignore[reportAttributeAccessIssue]
                UNITS_DICT[unit], F.lit(multiple), expr
            )

        return self.compliant._with_callable(_offset_by)

    def _no_op_time_zone(self, time_zone: str) -> SparkLikeExpr:  # pragma: no cover
        def func(df: SparkLikeLazyFrame) -> Sequence[Column]:
            native_series_list = self.compliant(df)
            conn_time_zone = fetch_session_time_zone(df.native.sparkSession)
            if conn_time_zone != time_zone:
                msg = (
                    "PySpark stores the time zone in the session, rather than in the "
                    f"data type, so changing the timezone to anything other than {conn_time_zone} "
                    " (the current session time zone) is not supported."
                )
                raise NotImplementedError(msg)
            return native_series_list

        return self.compliant.__class__(
            func,
            evaluate_output_names=self.compliant._evaluate_output_names,
            alias_output_names=self.compliant._alias_output_names,
            version=self.compliant._version,
            implementation=self.compliant._implementation,
        )

    def convert_time_zone(self, time_zone: str) -> SparkLikeExpr:  # pragma: no cover
        return self._no_op_time_zone(time_zone)

    def replace_time_zone(
        self, time_zone: str | None
    ) -> SparkLikeExpr:  # pragma: no cover
        if time_zone is None:
            return self.compliant._with_elementwise(
                lambda expr: expr.cast("timestamp_ntz")
            )
        return self._no_op_time_zone(time_zone)

    def _format_iso_week_with_day(self, expr: Column) -> Column:
        """Format datetime as ISO week string with day."""
        F = self.compliant._F

        year = F.date_format(expr, "yyyy")
        week = F.lpad(F.weekofyear(expr).cast("string"), 2, "0")
        day = self._weekday(expr)
        return F.concat(year, F.lit("-W"), week, F.lit("-"), day.cast("string"))

    def _format_iso_week(self, expr: Column) -> Column:
        """Format datetime as ISO week string."""
        F = self.compliant._F

        year = F.date_format(expr, "yyyy")
        week = F.lpad(F.weekofyear(expr).cast("string"), 2, "0")
        return F.concat(year, F.lit("-W"), week)

    def _format_microseconds(
        self, expr: Column, format: str
    ) -> tuple[str, tuple[Column, ...]]:
        """Format microseconds if present in format, else it's a no-op."""
        F = self.compliant._F

        suffix: tuple[Column, ...]
        if format.endswith((".%f", "%.f")):
            import re

            micros = F.unix_micros(expr) % US_PER_SECOND
            micros_str = F.lpad(micros.cast("string"), 6, "0")
            suffix = (F.lit("."), micros_str)
            format_ = re.sub(r"(.%|%.)f$", "", format)
            return format_, suffix

        return format, ()

    timestamp = not_implemented()
    total_seconds = not_implemented()
    total_minutes = not_implemented()
    total_milliseconds = not_implemented()
    total_microseconds = not_implemented()
    total_nanoseconds = not_implemented()
