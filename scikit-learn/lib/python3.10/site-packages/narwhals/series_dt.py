from __future__ import annotations

from typing import TYPE_CHECKING, Generic

from narwhals.typing import SeriesT

if TYPE_CHECKING:
    from narwhals.typing import TimeUnit


class SeriesDateTimeNamespace(Generic[SeriesT]):
    def __init__(self, series: SeriesT) -> None:
        self._narwhals_series = series

    def date(self) -> SeriesT:
        """Get the date in a datetime series.

        Raises:
            NotImplementedError: If pandas default backend is being used.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [datetime(2012, 1, 7, 10, 20), datetime(2023, 3, 10, 11, 32)]
            ... ).convert_dtypes(dtype_backend="pyarrow")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.date().to_native()
            0    2012-01-07
            1    2023-03-10
            dtype: date32[day][pyarrow]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.date()
        )

    def year(self) -> SeriesT:
        """Get the year in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([datetime(2012, 1, 7), datetime(2023, 3, 10)])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.year().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i32]
            [
                    2012
                    2023
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.year()
        )

    def month(self) -> SeriesT:
        """Gets the month in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([datetime(2012, 1, 7), datetime(2023, 3, 10)])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.month().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i8]
            [
                    1
                    3
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.month()
        )

    def day(self) -> SeriesT:
        """Extracts the day in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array(
            ...     [[datetime(2022, 1, 1), datetime(2022, 1, 5)]]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.day().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                5
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.day()
        )

    def hour(self) -> SeriesT:
        """Extracts the hour in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array(
            ...     [[datetime(2022, 1, 1, 5, 3), datetime(2022, 1, 5, 9, 12)]]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.hour().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                5,
                9
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.hour()
        )

    def minute(self) -> SeriesT:
        """Extracts the minute in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [datetime(2022, 1, 1, 5, 3), datetime(2022, 1, 5, 9, 12)]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.minute().to_native()
            0     3
            1    12
            dtype: int32
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.minute()
        )

    def second(self) -> SeriesT:
        """Extracts the seconds in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [datetime(2022, 1, 1, 5, 3, 10), datetime(2022, 1, 5, 9, 12, 4)]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.second().to_native()
            0    10
            1     4
            dtype: int32
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.second()
        )

    def millisecond(self) -> SeriesT:
        """Extracts the milliseconds in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [
            ...         datetime(2022, 1, 1, 5, 3, 7, 400000),
            ...         datetime(2022, 1, 1, 5, 3, 7, 0),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.millisecond().alias("datetime").to_native()
            0    400
            1      0
            Name: datetime, dtype: int32
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.millisecond()
        )

    def microsecond(self) -> SeriesT:
        """Extracts the microseconds in a datetime series.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [
            ...         datetime(2022, 1, 1, 5, 3, 7, 400000),
            ...         datetime(2022, 1, 1, 5, 3, 7, 0),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.microsecond().alias("datetime").to_native()
            0    400000
            1         0
            Name: datetime, dtype: int32
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.microsecond()
        )

    def nanosecond(self) -> SeriesT:
        """Extract the nanoseconds in a date series.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [
            ...         datetime(2022, 1, 1, 5, 3, 7, 400000),
            ...         datetime(2022, 1, 1, 5, 3, 7, 0),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.nanosecond().alias("datetime").to_native()
            0    400000000
            1            0
            Name: datetime, dtype: int32
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.nanosecond()
        )

    def ordinal_day(self) -> SeriesT:
        """Get ordinal day.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array(
            ...     [[datetime(2020, 1, 1), datetime(2020, 8, 3)]]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.ordinal_day().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                1,
                216
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.ordinal_day()
        )

    def weekday(self) -> SeriesT:
        """Extract the week day in a datetime series.

        Note that Monday = 1 and Sunday = 7.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array(
            ...     [[datetime(2020, 1, 1), datetime(2020, 8, 3)]]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.weekday().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                3,
                1
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.weekday()
        )

    def total_minutes(self) -> SeriesT:
        """Get total minutes.

        Notes:
            The function outputs the total minutes in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [timedelta(minutes=10), timedelta(minutes=20, seconds=40)]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.total_minutes().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
                    10
                    20
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.total_minutes()
        )

    def total_seconds(self) -> SeriesT:
        """Get total seconds.

        Notes:
            The function outputs the total seconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [timedelta(minutes=10), timedelta(minutes=20, seconds=40)]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.total_seconds().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
                    600
                    1240
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.total_seconds()
        )

    def total_milliseconds(self) -> SeriesT:
        """Get total milliseconds.

        Notes:
            The function outputs the total milliseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [
            ...         timedelta(milliseconds=10),
            ...         timedelta(milliseconds=20, microseconds=40),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.total_milliseconds().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
                    10
                    20
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.total_milliseconds()
        )

    def total_microseconds(self) -> SeriesT:
        """Get total microseconds.

        Notes:
            The function outputs the total microseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [
            ...         timedelta(microseconds=10),
            ...         timedelta(milliseconds=1, microseconds=200),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.total_microseconds().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
                    10
                    1200
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.total_microseconds()
        )

    def total_nanoseconds(self) -> SeriesT:
        """Get total nanoseconds.

        Notes:
            The function outputs the total nanoseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` in this case.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     ["2024-01-01 00:00:00.000000001", "2024-01-01 00:00:00.000000002"]
            ... ).str.to_datetime(time_unit="ns")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.diff().dt.total_nanoseconds().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [i64]
            [
                    null
                    1
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.total_nanoseconds()
        )

    def to_string(self, format: str) -> SeriesT:
        """Convert a Date/Time/Datetime series into a String series with the given format.

        Arguments:
            format: Format string for converting the datetime to string.

        Notes:
            Unfortunately, different libraries interpret format directives a bit
            differently.

            - Chrono, the library used by Polars, uses `"%.f"` for fractional seconds,
              whereas pandas and Python stdlib use `".%f"`.
            - PyArrow interprets `"%S"` as "seconds, including fractional seconds"
              whereas most other tools interpret it as "just seconds, as 2 digits".
            ---
            Therefore, we make the following adjustments.

            - for pandas-like libraries, we replace `"%S.%f"` with `"%S%.f"`.
            - for PyArrow, we replace `"%S.%f"` with `"%S"`.
            ---
            Workarounds like these don't make us happy, and we try to avoid them as
            much as possible, but here we feel like it's the best compromise.

            If you just want to format a date/datetime Series as a local datetime
            string, and have it work as consistently as possible across libraries,
            we suggest using:

            - `"%Y-%m-%dT%H:%M:%S%.f"` for datetimes
            - `"%Y-%m-%d"` for dates
            ---
            Though note that, even then, different tools may return a different number
            of trailing zeros. Nonetheless, this is probably consistent enough for
            most applications.

            If you have an application where this is not enough, please open an issue
            and let us know.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([datetime(2020, 3, 1), datetime(2020, 4, 1)])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.to_string("%Y/%m/%d")
            ┌───────────────┐
            |Narwhals Series|
            |---------------|
            |0    2020/03/01|
            |1    2020/04/01|
            |dtype: object  |
            └───────────────┘
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.to_string(format)
        )

    def replace_time_zone(self, time_zone: str | None) -> SeriesT:
        """Replace time zone.

        Arguments:
            time_zone: Target time zone.

        Examples:
            >>> from datetime import datetime, timezone
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [
            ...         datetime(2024, 1, 1, tzinfo=timezone.utc),
            ...         datetime(2024, 1, 2, tzinfo=timezone.utc),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.replace_time_zone(
            ...     "Asia/Kathmandu"
            ... ).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [datetime[μs, Asia/Kathmandu]]
            [
                    2024-01-01 00:00:00 +0545
                    2024-01-02 00:00:00 +0545
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.replace_time_zone(time_zone)
        )

    def convert_time_zone(self, time_zone: str) -> SeriesT:
        """Convert time zone.

        If converting from a time-zone-naive column, then conversion happens
        as if converting from UTC.

        Arguments:
            time_zone: Target time zone.

        Examples:
            >>> from datetime import datetime, timezone
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [
            ...         datetime(2024, 1, 1, tzinfo=timezone.utc),
            ...         datetime(2024, 1, 2, tzinfo=timezone.utc),
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.convert_time_zone("Asia/Kathmandu").to_native()
            0   2024-01-01 05:45:00+05:45
            1   2024-01-02 05:45:00+05:45
            dtype: datetime64[ns, Asia/Kathmandu]
        """
        if time_zone is None:
            msg = "Target `time_zone` cannot be `None` in `convert_time_zone`. Please use `replace_time_zone(None)` if you want to remove the time zone."
            raise TypeError(msg)
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.convert_time_zone(time_zone)
        )

    def timestamp(self, time_unit: TimeUnit) -> SeriesT:
        """Return a timestamp in the given time unit.

        Arguments:
            time_unit: One of
                - 'ns': nanosecond.
                - 'us': microsecond.
                - 'ms': millisecond.

        Examples:
            >>> from datetime import date
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(
            ...     [date(2001, 1, 1), None, date(2001, 1, 3)], dtype="datetime64[ns]"
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.timestamp("ms").to_native()
            0    9.783072e+11
            1             NaN
            2    9.784800e+11
            dtype: float64
        """
        if time_unit not in {"ns", "us", "ms"}:
            msg = (
                "invalid `time_unit`"
                f"\n\nExpected one of {{'ns', 'us', 'ms'}}, got {time_unit!r}."
            )
            raise ValueError(msg)
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.timestamp(time_unit)
        )

    def truncate(self, every: str) -> SeriesT:
        """Divide the date/datetime range into buckets.

        Arguments:
            every: Length of bucket. Must be of form `<multiple><unit>`,
                where `multiple` is a positive integer and `unit` is one of

                - 'ns': nanosecond.
                - 'us': microsecond.
                - 'ms': millisecond.
                - 's': second.
                - 'm': minute.
                - 'h': hour.
                - 'd': day.
                - 'mo': month.
                - 'q': quarter.
                - 'y': year.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([datetime(2021, 3, 1, 12, 34)])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.truncate("1h").to_native()
            0   2021-03-01 12:00:00
            dtype: datetime64[ns]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.truncate(every)
        )

    def offset_by(self, by: str) -> SeriesT:
        """Offset this date by a relative time offset.

        Arguments:
            by: The offset. Must be of form `<multiple><unit>`,
                where `multiple` is a positive integer and `unit` is one of

                - 'ns': nanosecond.
                - 'us': microsecond.
                - 'ms': millisecond.
                - 's': second.
                - 'm': minute.
                - 'h': hour.
                - 'd': day.
                - 'mo': month.
                - 'q': quarter.
                - 'y': year.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series([datetime(2021, 3, 1, 12, 34)])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.dt.offset_by("1h").to_native()
            0   2021-03-01 13:34:00
            dtype: datetime64[ns]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.dt.offset_by(by)
        )
