from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, ClassVar, cast

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import UNITS_DICT, ArrowSeriesNamespace, floordiv_compat, lit
from narwhals._compliant.any_namespace import DateTimeNamespace
from narwhals._constants import (
    MS_PER_MINUTE,
    MS_PER_SECOND,
    NS_PER_MICROSECOND,
    NS_PER_MILLISECOND,
    NS_PER_MINUTE,
    NS_PER_SECOND,
    SECONDS_PER_DAY,
    SECONDS_PER_MINUTE,
    US_PER_MINUTE,
    US_PER_SECOND,
)
from narwhals._duration import Interval

if TYPE_CHECKING:
    from collections.abc import Mapping

    from typing_extensions import TypeAlias

    from narwhals._arrow.series import ArrowSeries
    from narwhals._arrow.typing import ChunkedArrayAny, ScalarAny
    from narwhals.dtypes import Datetime
    from narwhals.typing import TimeUnit

    UnitCurrent: TypeAlias = TimeUnit
    UnitTarget: TypeAlias = TimeUnit
    BinOpBroadcast: TypeAlias = Callable[[ChunkedArrayAny, ScalarAny], ChunkedArrayAny]
    IntoRhs: TypeAlias = int


class ArrowSeriesDateTimeNamespace(
    ArrowSeriesNamespace, DateTimeNamespace["ArrowSeries"]
):
    _TIMESTAMP_DATE_FACTOR: ClassVar[Mapping[TimeUnit, int]] = {
        "ns": NS_PER_SECOND,
        "us": US_PER_SECOND,
        "ms": MS_PER_SECOND,
        "s": 1,
    }
    _TIMESTAMP_DATETIME_OP_FACTOR: ClassVar[
        Mapping[tuple[UnitCurrent, UnitTarget], tuple[BinOpBroadcast, IntoRhs]]
    ] = {
        ("ns", "us"): (floordiv_compat, 1_000),
        ("ns", "ms"): (floordiv_compat, 1_000_000),
        ("us", "ns"): (pc.multiply, NS_PER_MICROSECOND),
        ("us", "ms"): (floordiv_compat, 1_000),
        ("ms", "ns"): (pc.multiply, NS_PER_MILLISECOND),
        ("ms", "us"): (pc.multiply, 1_000),
        ("s", "ns"): (pc.multiply, NS_PER_SECOND),
        ("s", "us"): (pc.multiply, US_PER_SECOND),
        ("s", "ms"): (pc.multiply, MS_PER_SECOND),
    }

    @property
    def unit(self) -> TimeUnit:  # NOTE: Unsafe (native).
        return cast("pa.TimestampType[TimeUnit, Any]", self.native.type).unit

    @property
    def time_zone(self) -> str | None:  # NOTE: Unsafe (narwhals).
        return cast("Datetime", self.compliant.dtype).time_zone

    def to_string(self, format: str) -> ArrowSeries:
        # PyArrow differs from other libraries in that %S also prints out
        # the fractional part of the second...:'(
        # https://arrow.apache.org/docs/python/generated/pyarrow.compute.strftime.html
        format = format.replace("%S.%f", "%S").replace("%S%.f", "%S")
        return self.with_native(pc.strftime(self.native, format))

    def replace_time_zone(self, time_zone: str | None) -> ArrowSeries:
        if time_zone is not None:
            result = pc.assume_timezone(pc.local_timestamp(self.native), time_zone)
        else:
            result = pc.local_timestamp(self.native)
        return self.with_native(result)

    def convert_time_zone(self, time_zone: str) -> ArrowSeries:
        ser = self.replace_time_zone("UTC") if self.time_zone is None else self.compliant
        return self.with_native(ser.native.cast(pa.timestamp(self.unit, time_zone)))

    def timestamp(self, time_unit: TimeUnit) -> ArrowSeries:
        ser = self.compliant
        dtypes = ser._version.dtypes
        if isinstance(ser.dtype, dtypes.Datetime):
            current = ser.dtype.time_unit
            s_cast = self.native.cast(pa.int64())
            if current == time_unit:
                result = s_cast
            elif item := self._TIMESTAMP_DATETIME_OP_FACTOR.get((current, time_unit)):
                fn, factor = item
                result = fn(s_cast, lit(factor))
            else:  # pragma: no cover
                msg = f"unexpected time unit {current}, please report an issue at https://github.com/narwhals-dev/narwhals"
                raise AssertionError(msg)
            return self.with_native(result)
        if isinstance(ser.dtype, dtypes.Date):
            time_s = pc.multiply(self.native.cast(pa.int32()), lit(SECONDS_PER_DAY))
            factor = self._TIMESTAMP_DATE_FACTOR[time_unit]
            return self.with_native(pc.multiply(time_s, lit(factor)))
        msg = "Input should be either of Date or Datetime type"
        raise TypeError(msg)

    def date(self) -> ArrowSeries:
        return self.with_native(self.native.cast(pa.date32()))

    def year(self) -> ArrowSeries:
        return self.with_native(pc.year(self.native))

    def month(self) -> ArrowSeries:
        return self.with_native(pc.month(self.native))

    def day(self) -> ArrowSeries:
        return self.with_native(pc.day(self.native))

    def hour(self) -> ArrowSeries:
        return self.with_native(pc.hour(self.native))

    def minute(self) -> ArrowSeries:
        return self.with_native(pc.minute(self.native))

    def second(self) -> ArrowSeries:
        return self.with_native(pc.second(self.native))

    def millisecond(self) -> ArrowSeries:
        return self.with_native(pc.millisecond(self.native))

    def microsecond(self) -> ArrowSeries:
        arr = self.native
        result = pc.add(pc.multiply(pc.millisecond(arr), lit(1000)), pc.microsecond(arr))
        return self.with_native(result)

    def nanosecond(self) -> ArrowSeries:
        result = pc.add(
            pc.multiply(self.microsecond().native, lit(1000)), pc.nanosecond(self.native)
        )
        return self.with_native(result)

    def ordinal_day(self) -> ArrowSeries:
        return self.with_native(pc.day_of_year(self.native))

    def weekday(self) -> ArrowSeries:
        return self.with_native(pc.day_of_week(self.native, count_from_zero=False))

    def total_minutes(self) -> ArrowSeries:
        unit_to_minutes_factor = {
            "s": SECONDS_PER_MINUTE,
            "ms": MS_PER_MINUTE,
            "us": US_PER_MINUTE,
            "ns": NS_PER_MINUTE,
        }
        factor = lit(unit_to_minutes_factor[self.unit], type=pa.int64())
        return self.with_native(pc.divide(self.native, factor).cast(pa.int64()))

    def total_seconds(self) -> ArrowSeries:
        unit_to_seconds_factor = {
            "s": 1,
            "ms": MS_PER_SECOND,
            "us": US_PER_SECOND,
            "ns": NS_PER_SECOND,
        }
        factor = lit(unit_to_seconds_factor[self.unit], type=pa.int64())
        return self.with_native(pc.divide(self.native, factor).cast(pa.int64()))

    def total_milliseconds(self) -> ArrowSeries:
        unit_to_milli_factor = {
            "s": 1e3,  # seconds
            "ms": 1,  # milli
            "us": 1e3,  # micro
            "ns": 1e6,  # nano
        }
        factor = lit(unit_to_milli_factor[self.unit], type=pa.int64())
        if self.unit == "s":
            return self.with_native(pc.multiply(self.native, factor).cast(pa.int64()))
        return self.with_native(pc.divide(self.native, factor).cast(pa.int64()))

    def total_microseconds(self) -> ArrowSeries:
        unit_to_micro_factor = {
            "s": 1e6,  # seconds
            "ms": 1e3,  # milli
            "us": 1,  # micro
            "ns": 1e3,  # nano
        }
        factor = lit(unit_to_micro_factor[self.unit], type=pa.int64())
        if self.unit in {"s", "ms"}:
            return self.with_native(pc.multiply(self.native, factor).cast(pa.int64()))
        return self.with_native(pc.divide(self.native, factor).cast(pa.int64()))

    def total_nanoseconds(self) -> ArrowSeries:
        unit_to_nano_factor = {
            "s": NS_PER_SECOND,
            "ms": NS_PER_MILLISECOND,
            "us": NS_PER_MICROSECOND,
            "ns": 1,
        }
        factor = lit(unit_to_nano_factor[self.unit], type=pa.int64())
        return self.with_native(pc.multiply(self.native, factor).cast(pa.int64()))

    def truncate(self, every: str) -> ArrowSeries:
        interval = Interval.parse(every)
        return self.with_native(
            pc.floor_temporal(self.native, interval.multiple, UNITS_DICT[interval.unit])
        )

    def offset_by(self, by: str) -> ArrowSeries:
        interval = Interval.parse_no_constraints(by)
        native = self.native
        if interval.unit in {"y", "q", "mo"}:
            msg = f"Offsetting by {interval.unit} is not yet supported for pyarrow."
            raise NotImplementedError(msg)
        dtype = self.compliant.dtype
        datetime_dtype = self.version.dtypes.Datetime
        if interval.unit == "d" and isinstance(dtype, datetime_dtype) and dtype.time_zone:
            offset: pa.DurationScalar[Any] = lit(interval.to_timedelta())
            native_naive = pc.local_timestamp(native)
            result = pc.assume_timezone(pc.add(native_naive, offset), dtype.time_zone)
            return self.with_native(result)
        if interval.unit == "ns":  # pragma: no cover
            offset = lit(interval.multiple, pa.duration("ns"))  # type: ignore[assignment]
        else:
            offset = lit(interval.to_timedelta())
        return self.with_native(pc.add(native, offset))
