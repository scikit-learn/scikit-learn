from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable

from narwhals._duration import Interval
from narwhals._ibis.utils import (
    UNITS_DICT_BUCKET,
    UNITS_DICT_TRUNCATE,
    timedelta_to_ibis_interval,
)
from narwhals._sql.expr_dt import SQLExprDateTimeNamesSpace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    import ibis.expr.types as ir

    from narwhals._ibis.expr import IbisExpr
    from narwhals._ibis.utils import BucketUnit, TruncateUnit


class IbisExprDateTimeNamespace(SQLExprDateTimeNamesSpace["IbisExpr"]):
    def millisecond(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.millisecond())

    def microsecond(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.microsecond())

    def to_string(self, format: str) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.strftime(format))

    def weekday(self) -> IbisExpr:
        # Ibis uses 0-6 for Monday-Sunday. Add 1 to match polars.
        return self.compliant._with_callable(lambda expr: expr.day_of_week.index() + 1)

    def _bucket(self, kwds: dict[BucketUnit, Any], /) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.TimestampValue) -> ir.TimestampValue:
            return expr.bucket(**kwds)

        return fn

    def _truncate(self, unit: TruncateUnit, /) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.TimestampValue) -> ir.TimestampValue:
            return expr.truncate(unit)

        return fn

    def truncate(self, every: str) -> IbisExpr:
        interval = Interval.parse(every)
        multiple, unit = interval.multiple, interval.unit
        if unit == "q":
            multiple, unit = 3 * multiple, "mo"
        if multiple != 1:
            if self.compliant._backend_version < (7, 1):  # pragma: no cover
                msg = "Truncating datetimes with multiples of the unit is only supported in Ibis >= 7.1."
                raise NotImplementedError(msg)
            fn = self._bucket({UNITS_DICT_BUCKET[unit]: multiple})
        else:
            fn = self._truncate(UNITS_DICT_TRUNCATE[unit])
        return self.compliant._with_callable(fn)

    def offset_by(self, by: str) -> IbisExpr:
        interval = Interval.parse_no_constraints(by)
        unit = interval.unit
        if unit in {"y", "q", "mo", "d", "ns"}:
            msg = f"Offsetting by {unit} is not yet supported for ibis."
            raise NotImplementedError(msg)
        offset = timedelta_to_ibis_interval(interval.to_timedelta())
        return self.compliant._with_callable(lambda expr: expr.add(offset))

    def replace_time_zone(self, time_zone: str | None) -> IbisExpr:
        if time_zone is None:
            return self.compliant._with_callable(lambda expr: expr.cast("timestamp"))
        msg = "`replace_time_zone` with non-null `time_zone` not yet implemented for Ibis"  # pragma: no cover
        raise NotImplementedError(msg)

    nanosecond = not_implemented()
    total_minutes = not_implemented()
    total_seconds = not_implemented()
    total_milliseconds = not_implemented()
    total_microseconds = not_implemented()
    total_nanoseconds = not_implemented()
    convert_time_zone = not_implemented()
    timestamp = not_implemented()
