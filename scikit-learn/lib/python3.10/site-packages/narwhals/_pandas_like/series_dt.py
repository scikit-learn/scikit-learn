from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals._compliant.any_namespace import DateTimeNamespace
from narwhals._constants import (
    EPOCH_YEAR,
    MS_PER_SECOND,
    NS_PER_SECOND,
    SECONDS_PER_DAY,
    US_PER_SECOND,
)
from narwhals._duration import Interval
from narwhals._pandas_like.utils import (
    ALIAS_DICT,
    UNITS_DICT,
    PandasLikeSeriesNamespace,
    calculate_timestamp_date,
    calculate_timestamp_datetime,
    get_dtype_backend,
    int_dtype_mapper,
    is_dtype_pyarrow,
)

if TYPE_CHECKING:
    from datetime import timedelta

    import pandas as pd

    from narwhals._pandas_like.series import PandasLikeSeries
    from narwhals.typing import TimeUnit


class PandasLikeSeriesDateTimeNamespace(
    PandasLikeSeriesNamespace, DateTimeNamespace["PandasLikeSeries"]
):
    def date(self) -> PandasLikeSeries:
        result = self.with_native(self.native.dt.date)
        if str(result.dtype).lower() == "object":
            msg = (
                "Accessing `date` on the default pandas backend "
                "will return a Series of type `object`."
                "\nThis differs from polars API and will prevent `.dt` chaining. "
                "Please switch to the `pyarrow` backend:"
                '\ndf.convert_dtypes(dtype_backend="pyarrow")'
            )
            raise NotImplementedError(msg)
        return result

    def year(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.year)

    def month(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.month)

    def day(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.day)

    def hour(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.hour)

    def minute(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.minute)

    def second(self) -> PandasLikeSeries:
        return self.with_native(self.native.dt.second)

    def millisecond(self) -> PandasLikeSeries:
        return self.microsecond() // 1000

    def microsecond(self) -> PandasLikeSeries:
        if self.backend_version < (3, 0, 0) and self._is_pyarrow():
            # crazy workaround for https://github.com/pandas-dev/pandas/issues/59154
            import pyarrow.compute as pc  # ignore-banned-import()

            from narwhals._arrow.utils import lit

            arr_ns = self.native.array
            arr = arr_ns.__arrow_array__()
            result_arr = pc.add(
                pc.multiply(pc.millisecond(arr), lit(1_000)), pc.microsecond(arr)
            )
            result = type(self.native)(type(arr_ns)(result_arr), name=self.native.name)
            return self.with_native(result)

        return self.with_native(self.native.dt.microsecond)

    def nanosecond(self) -> PandasLikeSeries:
        return self.microsecond() * 1_000 + self.native.dt.nanosecond

    def ordinal_day(self) -> PandasLikeSeries:
        year_start = self.native.dt.year
        result = (
            self.native.to_numpy().astype("datetime64[D]")
            - (year_start.to_numpy() - EPOCH_YEAR).astype("datetime64[Y]")
        ).astype("int32") + 1
        dtype = "Int64[pyarrow]" if self._is_pyarrow() else "int32"
        return self.with_native(
            type(self.native)(result, dtype=dtype, name=year_start.name)
        )

    def weekday(self) -> PandasLikeSeries:
        # Pandas is 0-6 while Polars is 1-7
        return self.with_native(self.native.dt.weekday) + 1

    def _is_pyarrow(self) -> bool:
        return is_dtype_pyarrow(self.native.dtype)

    def _get_total_seconds(self) -> Any:
        if hasattr(self.native.dt, "total_seconds"):
            return self.native.dt.total_seconds()
        return (  # pragma: no cover
            self.native.dt.days * SECONDS_PER_DAY
            + self.native.dt.seconds
            + (self.native.dt.microseconds / US_PER_SECOND)
            + (self.native.dt.nanoseconds / NS_PER_SECOND)
        )

    def total_minutes(self) -> PandasLikeSeries:
        s = self._get_total_seconds()
        # this calculates the sign of each series element
        s_sign = 2 * (s > 0).astype(int_dtype_mapper(s.dtype)) - 1
        s_abs = s.abs() // 60
        if ~s.isna().any():
            s_abs = s_abs.astype(int_dtype_mapper(s.dtype))
        return self.with_native(s_abs * s_sign)

    def total_seconds(self) -> PandasLikeSeries:
        s = self._get_total_seconds()
        # this calculates the sign of each series element
        s_sign = 2 * (s > 0).astype(int_dtype_mapper(s.dtype)) - 1
        s_abs = s.abs() // 1
        if ~s.isna().any():
            s_abs = s_abs.astype(int_dtype_mapper(s.dtype))
        return self.with_native(s_abs * s_sign)

    def total_milliseconds(self) -> PandasLikeSeries:
        s = self._get_total_seconds() * MS_PER_SECOND
        # this calculates the sign of each series element
        s_sign = 2 * (s > 0).astype(int_dtype_mapper(s.dtype)) - 1
        s_abs = s.abs() // 1
        if ~s.isna().any():
            s_abs = s_abs.astype(int_dtype_mapper(s.dtype))
        return self.with_native(s_abs * s_sign)

    def total_microseconds(self) -> PandasLikeSeries:
        s = self._get_total_seconds() * US_PER_SECOND
        # this calculates the sign of each series element
        s_sign = 2 * (s > 0).astype(int_dtype_mapper(s.dtype)) - 1
        s_abs = s.abs() // 1
        if ~s.isna().any():
            s_abs = s_abs.astype(int_dtype_mapper(s.dtype))
        return self.with_native(s_abs * s_sign)

    def total_nanoseconds(self) -> PandasLikeSeries:
        s = self._get_total_seconds() * NS_PER_SECOND
        # this calculates the sign of each series element
        s_sign = 2 * (s > 0).astype(int_dtype_mapper(s.dtype)) - 1
        s_abs = s.abs() // 1
        if ~s.isna().any():
            s_abs = s_abs.astype(int_dtype_mapper(s.dtype))
        return self.with_native(s_abs * s_sign)

    def to_string(self, format: str) -> PandasLikeSeries:
        # Polars' parser treats `'%.f'` as pandas does `'.%f'`
        # PyArrow interprets `'%S'` as "seconds, plus fractional seconds"
        # and doesn't support `%f`
        if not self._is_pyarrow():
            format = format.replace("%S%.f", "%S.%f")
        else:
            format = format.replace("%S.%f", "%S").replace("%S%.f", "%S")
        return self.with_native(self.native.dt.strftime(format))

    def replace_time_zone(self, time_zone: str | None) -> PandasLikeSeries:
        de_zone = self.native.dt.tz_localize(None)
        result = de_zone.dt.tz_localize(time_zone) if time_zone is not None else de_zone
        return self.with_native(result)

    def convert_time_zone(self, time_zone: str) -> PandasLikeSeries:
        if self.compliant.dtype.time_zone is None:  # type: ignore[attr-defined]
            result = self.native.dt.tz_localize("UTC").dt.tz_convert(time_zone)
        else:
            result = self.native.dt.tz_convert(time_zone)
        return self.with_native(result)

    def timestamp(self, time_unit: TimeUnit) -> PandasLikeSeries:
        s = self.native
        dtype = self.compliant.dtype
        mask_na = s.isna()
        dtypes = self.version.dtypes
        if dtype == dtypes.Date:
            # Date is only supported in pandas dtypes if pyarrow-backed
            s_cast = s.astype("Int32[pyarrow]")
            result = calculate_timestamp_date(s_cast, time_unit)
        elif isinstance(dtype, dtypes.Datetime):
            fn = (
                s.view
                if (self.implementation.is_pandas() and self.backend_version < (2,))
                else s.astype
            )
            s_cast = fn("Int64[pyarrow]") if self._is_pyarrow() else fn("int64")
            result = calculate_timestamp_datetime(s_cast, dtype.time_unit, time_unit)
        else:
            msg = "Input should be either of Date or Datetime type"
            raise TypeError(msg)
        result[mask_na] = None
        return self.with_native(result)

    def truncate(self, every: str) -> PandasLikeSeries:
        interval = Interval.parse(every)
        multiple, unit = interval.multiple, interval.unit
        native = self.native
        if self.implementation.is_cudf():
            if multiple != 1:
                msg = f"Only multiple `1` is supported for cuDF, got: {multiple}."
                raise NotImplementedError(msg)
            return self.with_native(self.native.dt.floor(ALIAS_DICT.get(unit, unit)))
        dtype_backend = get_dtype_backend(native.dtype, self.compliant._implementation)
        if unit in {"mo", "q", "y"}:
            if self.implementation.is_cudf():
                msg = f"Truncating to {unit} is not supported yet for cuDF."
                raise NotImplementedError(msg)
            if dtype_backend == "pyarrow":
                import pyarrow.compute as pc  # ignore-banned-import

                ca = native.array._pa_array
                result_arr = pc.floor_temporal(ca, multiple, UNITS_DICT[unit])
            else:
                if unit == "q":
                    multiple *= 3
                    np_unit = "M"
                elif unit == "mo":
                    np_unit = "M"
                else:
                    np_unit = "Y"
                arr = native.values  # noqa: PD011
                arr_dtype = arr.dtype
                result_arr = arr.astype(f"datetime64[{multiple}{np_unit}]").astype(
                    arr_dtype
                )
            result_native = type(native)(
                result_arr, dtype=native.dtype, index=native.index, name=native.name
            )
            return self.with_native(result_native)
        return self.with_native(
            self.native.dt.floor(f"{multiple}{ALIAS_DICT.get(unit, unit)}")
        )

    def offset_by(self, by: str) -> PandasLikeSeries:
        native = self.native
        pdx = self.compliant.__native_namespace__()
        if self._is_pyarrow():
            import pyarrow as pa  # ignore-banned-import

            compliant = self.compliant
            ca = pa.chunked_array([compliant.to_arrow()])  # type: ignore[arg-type]
            result = (
                compliant._version.namespace.from_backend("pyarrow")
                .compliant.from_native(ca)
                .dt.offset_by(by)
                .native
            )
            result_pd = native.__class__(
                result, dtype=native.dtype, index=native.index, name=native.name
            )
        else:
            interval = Interval.parse_no_constraints(by)
            multiple, unit = interval.multiple, interval.unit
            if unit == "q":
                multiple *= 3
                unit = "mo"
            offset: pd.DateOffset | timedelta
            if unit == "y":
                offset = pdx.DateOffset(years=multiple)
            elif unit == "mo":
                offset = pdx.DateOffset(months=multiple)
            elif unit == "ns":
                offset = pdx.Timedelta(multiple, unit=UNITS_DICT[unit])
            else:
                offset = interval.to_timedelta()
            dtype = self.compliant.dtype
            datetime_dtype = self.version.dtypes.Datetime
            if unit == "d" and isinstance(dtype, datetime_dtype) and dtype.time_zone:
                native_without_timezone = native.dt.tz_localize(None)
                result_pd = native_without_timezone + offset
                result_pd = result_pd.dt.tz_localize(dtype.time_zone)
            else:
                result_pd = native + offset

        return self.with_native(result_pd)
