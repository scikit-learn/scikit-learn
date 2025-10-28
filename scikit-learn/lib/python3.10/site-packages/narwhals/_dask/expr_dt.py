from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import DateTimeNamespace
from narwhals._constants import MS_PER_SECOND, NS_PER_SECOND, US_PER_SECOND
from narwhals._duration import Interval
from narwhals._pandas_like.utils import (
    ALIAS_DICT,
    calculate_timestamp_date,
    calculate_timestamp_datetime,
    native_to_narwhals_dtype,
)
from narwhals._utils import Implementation

if TYPE_CHECKING:
    import dask.dataframe.dask_expr as dx

    from narwhals._dask.expr import DaskExpr
    from narwhals.typing import TimeUnit


class DaskExprDateTimeNamespace(
    LazyExprNamespace["DaskExpr"], DateTimeNamespace["DaskExpr"]
):
    def date(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.date, "date")

    def year(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.year, "year")

    def month(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.month, "month")

    def day(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.day, "day")

    def hour(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.hour, "hour")

    def minute(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.minute, "minute")

    def second(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.dt.second, "second")

    def millisecond(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.microsecond // 1000, "millisecond"
        )

    def microsecond(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.microsecond, "microsecond"
        )

    def nanosecond(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.microsecond * 1000 + expr.dt.nanosecond, "nanosecond"
        )

    def ordinal_day(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.dayofyear, "ordinal_day"
        )

    def weekday(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.weekday + 1,  # Dask is 0-6
            "weekday",
        )

    def to_string(self, format: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, format: expr.dt.strftime(format.replace("%.f", ".%f")),
            "strftime",
            format=format,
        )

    def replace_time_zone(self, time_zone: str | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, time_zone: expr.dt.tz_localize(None).dt.tz_localize(time_zone)
            if time_zone is not None
            else expr.dt.tz_localize(None),
            "tz_localize",
            time_zone=time_zone,
        )

    def convert_time_zone(self, time_zone: str) -> DaskExpr:
        def func(s: dx.Series, time_zone: str) -> dx.Series:
            dtype = native_to_narwhals_dtype(
                s.dtype, self.compliant._version, Implementation.DASK
            )
            if dtype.time_zone is None:  # type: ignore[attr-defined]
                return s.dt.tz_localize("UTC").dt.tz_convert(time_zone)  # pyright: ignore[reportAttributeAccessIssue]
            return s.dt.tz_convert(time_zone)  # pyright: ignore[reportAttributeAccessIssue]

        return self.compliant._with_callable(func, "tz_convert", time_zone=time_zone)

    # ignoring coverage due to https://github.com/narwhals-dev/narwhals/issues/2808.
    def timestamp(self, time_unit: TimeUnit) -> DaskExpr:  # pragma: no cover
        def func(s: dx.Series, time_unit: TimeUnit) -> dx.Series:
            dtype = native_to_narwhals_dtype(
                s.dtype, self.compliant._version, Implementation.DASK
            )
            is_pyarrow_dtype = "pyarrow" in str(dtype)
            mask_na = s.isna()
            dtypes = self.compliant._version.dtypes
            if dtype == dtypes.Date:
                # Date is only supported in pandas dtypes if pyarrow-backed
                s_cast = s.astype("Int32[pyarrow]")
                result = calculate_timestamp_date(s_cast, time_unit)
            elif isinstance(dtype, dtypes.Datetime):
                original_time_unit = dtype.time_unit
                s_cast = (
                    s.astype("Int64[pyarrow]") if is_pyarrow_dtype else s.astype("int64")
                )
                result = calculate_timestamp_datetime(
                    s_cast, original_time_unit, time_unit
                )
            else:
                msg = "Input should be either of Date or Datetime type"
                raise TypeError(msg)
            return result.where(~mask_na)  # pyright: ignore[reportReturnType]

        return self.compliant._with_callable(func, "datetime", time_unit=time_unit)

    def total_minutes(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.total_seconds() // 60, "total_minutes"
        )

    def total_seconds(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.total_seconds() // 1, "total_seconds"
        )

    def total_milliseconds(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.total_seconds() * MS_PER_SECOND // 1,
            "total_milliseconds",
        )

    def total_microseconds(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.total_seconds() * US_PER_SECOND // 1,
            "total_microseconds",
        )

    def total_nanoseconds(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.dt.total_seconds() * NS_PER_SECOND // 1, "total_nanoseconds"
        )

    def truncate(self, every: str) -> DaskExpr:
        interval = Interval.parse(every)
        unit = interval.unit
        if unit in {"mo", "q", "y"}:
            msg = f"Truncating to {unit} is not yet supported for dask."
            raise NotImplementedError(msg)
        freq = f"{interval.multiple}{ALIAS_DICT.get(unit, unit)}"
        return self.compliant._with_callable(lambda expr: expr.dt.floor(freq), "truncate")

    def offset_by(self, by: str) -> DaskExpr:
        def func(s: dx.Series, by: str) -> dx.Series:
            interval = Interval.parse_no_constraints(by)
            unit = interval.unit
            if unit in {"y", "q", "mo", "d", "ns"}:
                msg = f"Offsetting by {unit} is not yet supported for dask."
                raise NotImplementedError(msg)
            offset = interval.to_timedelta()
            return s.add(offset)

        return self.compliant._with_callable(func, "offset_by", by=by)
