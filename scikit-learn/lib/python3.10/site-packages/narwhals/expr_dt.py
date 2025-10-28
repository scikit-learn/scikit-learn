from __future__ import annotations

from typing import TYPE_CHECKING, Generic, TypeVar

if TYPE_CHECKING:
    from narwhals.expr import Expr
    from narwhals.typing import TimeUnit

ExprT = TypeVar("ExprT", bound="Expr")


class ExprDateTimeNamespace(Generic[ExprT]):
    def __init__(self, expr: ExprT) -> None:
        self._expr = expr

    def date(self) -> ExprT:
        """Extract the date from underlying DateTime representation.

        Raises:
            NotImplementedError: If pandas default backend is being used.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"a": [datetime(2012, 1, 7, 10), datetime(2027, 12, 13)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").dt.date()).to_native()
            shape: (2, 1)
            ┌────────────┐
            │ a          │
            │ ---        │
            │ date       │
            ╞════════════╡
            │ 2012-01-07 │
            │ 2027-12-13 │
            └────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.date()
        )

    def year(self) -> ExprT:
        """Extract year from underlying DateTime representation.

        Returns the year number in the calendar date.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [datetime(1978, 6, 1), datetime(2065, 1, 1)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.year().alias("year"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |           a  year|
            |0 1978-06-01  1978|
            |1 2065-01-01  2065|
            └──────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.year()
        )

    def month(self) -> ExprT:
        """Extract month from underlying DateTime representation.

        Returns the month number starting from 1. The return value ranges from 1 to 12.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [datetime(1978, 6, 1), datetime(2065, 1, 1)]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.month().alias("month")).to_native()
            pyarrow.Table
            a: timestamp[us]
            month: int64
            ----
            a: [[1978-06-01 00:00:00.000000,2065-01-01 00:00:00.000000]]
            month: [[6,1]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.month()
        )

    def day(self) -> ExprT:
        """Extract day from underlying DateTime representation.

        Returns the day of month starting from 1. The return value ranges from 1 to 31. (The last day of month differs by months.)

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [datetime(1978, 6, 1), datetime(2065, 1, 1)]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.day().alias("day")).to_native()
            pyarrow.Table
            a: timestamp[us]
            day: int64
            ----
            a: [[1978-06-01 00:00:00.000000,2065-01-01 00:00:00.000000]]
            day: [[1,1]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.day()
        )

    def hour(self) -> ExprT:
        """Extract hour from underlying DateTime representation.

        Returns the hour number from 0 to 23.

        Examples:
            >>> from datetime import datetime
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"a": [datetime(1978, 1, 1, 1), datetime(2065, 1, 1, 10)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.hour().alias("hour"))
            ┌──────────────────────────────┐
            |      Narwhals DataFrame      |
            |------------------------------|
            |shape: (2, 2)                 |
            |┌─────────────────────┬──────┐|
            |│ a                   ┆ hour │|
            |│ ---                 ┆ ---  │|
            |│ datetime[μs]        ┆ i8   │|
            |╞═════════════════════╪══════╡|
            |│ 1978-01-01 01:00:00 ┆ 1    │|
            |│ 2065-01-01 10:00:00 ┆ 10   │|
            |└─────────────────────┴──────┘|
            └──────────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.hour()
        )

    def minute(self) -> ExprT:
        """Extract minutes from underlying DateTime representation.

        Returns the minute number from 0 to 59.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [datetime(1978, 1, 1, 1, 1), datetime(2065, 1, 1, 10, 20)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.minute().alias("minute")).to_native()
                                a  minute
            0 1978-01-01 01:01:00       1
            1 2065-01-01 10:20:00      20
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.minute()
        )

    def second(self) -> ExprT:
        """Extract seconds from underlying DateTime representation.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table(
            ...     {
            ...         "a": [
            ...             datetime(1978, 1, 1, 1, 1, 1),
            ...             datetime(2065, 1, 1, 10, 20, 30),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a").dt.second().alias("second")).to_native()
            pyarrow.Table
            a: timestamp[us]
            second: int64
            ----
            a: [[1978-01-01 01:01:01.000000,2065-01-01 10:20:30.000000]]
            second: [[1,30]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.second()
        )

    def millisecond(self) -> ExprT:
        """Extract milliseconds from underlying DateTime representation.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table(
            ...     {
            ...         "a": [
            ...             datetime(1978, 1, 1, 1, 1, 1, 0),
            ...             datetime(2065, 1, 1, 10, 20, 30, 67000),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").dt.millisecond().alias("millisecond")
            ... ).to_native()
            pyarrow.Table
            a: timestamp[us]
            millisecond: int64
            ----
            a: [[1978-01-01 01:01:01.000000,2065-01-01 10:20:30.067000]]
            millisecond: [[0,67]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.millisecond()
        )

    def microsecond(self) -> ExprT:
        """Extract microseconds from underlying DateTime representation.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table(
            ...     {
            ...         "a": [
            ...             datetime(1978, 1, 1, 1, 1, 1, 0),
            ...             datetime(2065, 1, 1, 10, 20, 30, 67000),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").dt.microsecond().alias("microsecond")
            ... ).to_native()
            pyarrow.Table
            a: timestamp[us]
            microsecond: int64
            ----
            a: [[1978-01-01 01:01:01.000000,2065-01-01 10:20:30.067000]]
            microsecond: [[0,67000]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.microsecond()
        )

    def nanosecond(self) -> ExprT:
        """Extract Nanoseconds from underlying DateTime representation.

        Examples:
            >>> from datetime import datetime
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table(
            ...     {
            ...         "a": [
            ...             datetime(1978, 1, 1, 1, 1, 1, 0),
            ...             datetime(2065, 1, 1, 10, 20, 30, 67000),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").dt.nanosecond().alias("nanosecond")
            ... ).to_native()
            pyarrow.Table
            a: timestamp[us]
            nanosecond: int64
            ----
            a: [[1978-01-01 01:01:01.000000,2065-01-01 10:20:30.067000]]
            nanosecond: [[0,67000000]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.nanosecond()
        )

    def ordinal_day(self) -> ExprT:
        """Get ordinal day.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [datetime(2020, 1, 1), datetime(2020, 8, 3)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_ordinal_day=nw.col("a").dt.ordinal_day())
            ┌───────────────────────────┐
            |    Narwhals DataFrame     |
            |---------------------------|
            |           a  a_ordinal_day|
            |0 2020-01-01              1|
            |1 2020-08-03            216|
            └───────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.ordinal_day()
        )

    def weekday(self) -> ExprT:
        """Extract the week day from the underlying Date representation.

        Note that Monday = 1 and Sunday = 7.

        Examples:
            >>> from datetime import datetime
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [datetime(2020, 1, 1), datetime(2020, 8, 3)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_week_day=nw.col("a").dt.weekday())
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |           a  a_week_day|
            |0 2020-01-01           3|
            |1 2020-08-03           1|
            └────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.weekday()
        )

    def total_minutes(self) -> ExprT:
        """Get total minutes.

        Notes:
            The function outputs the total minutes in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` and `cast` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"a": [timedelta(minutes=10), timedelta(minutes=20, seconds=40)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_total_minutes=nw.col("a").dt.total_minutes()
            ... ).to_native()
            shape: (2, 2)
            ┌──────────────┬─────────────────┐
            │ a            ┆ a_total_minutes │
            │ ---          ┆ ---             │
            │ duration[μs] ┆ i64             │
            ╞══════════════╪═════════════════╡
            │ 10m          ┆ 10              │
            │ 20m 40s      ┆ 20              │
            └──────────────┴─────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.total_minutes()
        )

    def total_seconds(self) -> ExprT:
        """Get total seconds.

        Notes:
            The function outputs the total seconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` and `cast` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"a": [timedelta(seconds=10), timedelta(seconds=20, milliseconds=40)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_total_seconds=nw.col("a").dt.total_seconds()
            ... ).to_native()
            shape: (2, 2)
            ┌──────────────┬─────────────────┐
            │ a            ┆ a_total_seconds │
            │ ---          ┆ ---             │
            │ duration[μs] ┆ i64             │
            ╞══════════════╪═════════════════╡
            │ 10s          ┆ 10              │
            │ 20s 40ms     ┆ 20              │
            └──────────────┴─────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.total_seconds()
        )

    def total_milliseconds(self) -> ExprT:
        """Get total milliseconds.

        Notes:
            The function outputs the total milliseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` and `cast` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {
            ...         "a": [
            ...             timedelta(milliseconds=10),
            ...             timedelta(milliseconds=20, microseconds=40),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_total_milliseconds=nw.col("a").dt.total_milliseconds()
            ... ).to_native()
            shape: (2, 2)
            ┌──────────────┬──────────────────────┐
            │ a            ┆ a_total_milliseconds │
            │ ---          ┆ ---                  │
            │ duration[μs] ┆ i64                  │
            ╞══════════════╪══════════════════════╡
            │ 10ms         ┆ 10                   │
            │ 20040µs      ┆ 20                   │
            └──────────────┴──────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.total_milliseconds()
        )

    def total_microseconds(self) -> ExprT:
        """Get total microseconds.

        Notes:
            The function outputs the total microseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` and `cast` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table(
            ...     {
            ...         "a": [
            ...             timedelta(microseconds=10),
            ...             timedelta(milliseconds=1, microseconds=200),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_total_microseconds=nw.col("a").dt.total_microseconds()
            ... ).to_native()
            pyarrow.Table
            a: duration[us]
            a_total_microseconds: int64
            ----
            a: [[10,1200]]
            a_total_microseconds: [[10,1200]]
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.total_microseconds()
        )

    def total_nanoseconds(self) -> ExprT:
        """Get total nanoseconds.

        Notes:
            The function outputs the total nanoseconds in the int dtype by default,
            however, pandas may change the dtype to float when there are missing values,
            consider using `fill_null()` and `cast` in this case.

        Examples:
            >>> from datetime import timedelta
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {
            ...         "a": pd.to_datetime(
            ...             [
            ...                 "2024-01-01 00:00:00.000000001",
            ...                 "2024-01-01 00:00:00.000000002",
            ...             ]
            ...         )
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_diff_total_nanoseconds=nw.col("a").diff().dt.total_nanoseconds()
            ... ).to_native()
                                          a  a_diff_total_nanoseconds
            0 2024-01-01 00:00:00.000000001                       NaN
            1 2024-01-01 00:00:00.000000002                       1.0
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.total_nanoseconds()
        )

    def to_string(self, format: str) -> ExprT:
        """Convert a Date/Time/Datetime column into a String column with the given format.

        Arguments:
            format: Format to format temporal column with.

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
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"a": [datetime(2020, 3, 1), datetime(2020, 5, 1)]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").dt.to_string("%Y/%m/%d %H:%M:%S"))
            ┌───────────────────────┐
            |  Narwhals DataFrame   |
            |-----------------------|
            |shape: (2, 1)          |
            |┌─────────────────────┐|
            |│ a                   │|
            |│ ---                 │|
            |│ str                 │|
            |╞═════════════════════╡|
            |│ 2020/03/01 00:00:00 │|
            |│ 2020/05/01 00:00:00 │|
            |└─────────────────────┘|
            └───────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.to_string(format)
        )

    def replace_time_zone(self, time_zone: str | None) -> ExprT:
        """Replace time zone.

        Arguments:
            time_zone: Target time zone.

        Examples:
            >>> from datetime import datetime, timezone
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {
            ...         "a": [
            ...             datetime(2024, 1, 1, tzinfo=timezone.utc),
            ...             datetime(2024, 1, 2, tzinfo=timezone.utc),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").dt.replace_time_zone("Asia/Kathmandu")).to_native()
                                      a
            0 2024-01-01 00:00:00+05:45
            1 2024-01-02 00:00:00+05:45
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.replace_time_zone(time_zone)
        )

    def convert_time_zone(self, time_zone: str) -> ExprT:
        """Convert to a new time zone.

        If converting from a time-zone-naive column, then conversion happens
        as if converting from UTC.

        Arguments:
            time_zone: Target time zone.

        Examples:
            >>> from datetime import datetime, timezone
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {
            ...         "a": [
            ...             datetime(2024, 1, 1, tzinfo=timezone.utc),
            ...             datetime(2024, 1, 2, tzinfo=timezone.utc),
            ...         ]
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").dt.convert_time_zone("Asia/Kathmandu")).to_native()
                                      a
            0 2024-01-01 05:45:00+05:45
            1 2024-01-02 05:45:00+05:45
        """
        if time_zone is None:
            msg = "Target `time_zone` cannot be `None` in `convert_time_zone`. Please use `replace_time_zone(None)` if you want to remove the time zone."
            raise TypeError(msg)
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.convert_time_zone(time_zone)
        )

    def timestamp(self, time_unit: TimeUnit = "us") -> ExprT:
        """Return a timestamp in the given time unit.

        Arguments:
            time_unit: One of
                - 'ns': nanosecond.
                - 'us': microsecond.
                - 'ms': millisecond.

        Examples:
            >>> from datetime import date
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"date": [date(2001, 1, 1), None]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("date").dt.timestamp("ms").alias("timestamp_ms"))
            ┌─────────────────────────────┐
            |     Narwhals DataFrame      |
            |-----------------------------|
            |shape: (2, 2)                |
            |┌────────────┬──────────────┐|
            |│ date       ┆ timestamp_ms │|
            |│ ---        ┆ ---          │|
            |│ date       ┆ i64          │|
            |╞════════════╪══════════════╡|
            |│ 2001-01-01 ┆ 978307200000 │|
            |│ null       ┆ null         │|
            |└────────────┴──────────────┘|
            └─────────────────────────────┘
        """
        if time_unit not in {"ns", "us", "ms"}:
            msg = (
                "invalid `time_unit`"
                f"\n\nExpected one of {{'ns', 'us', 'ms'}}, got {time_unit!r}."
            )
            raise ValueError(msg)
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.timestamp(time_unit)
        )

    def truncate(self, every: str) -> ExprT:
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
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"datetime": [datetime(2021, 3, 1, 12, 34)]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("datetime").dt.truncate("1h").alias("datetime_trunc")
            ... )
            ┌─────────────────────────────────────────────┐
            |             Narwhals DataFrame              |
            |---------------------------------------------|
            |shape: (1, 2)                                |
            |┌─────────────────────┬─────────────────────┐|
            |│ datetime            ┆ datetime_trunc      │|
            |│ ---                 ┆ ---                 │|
            |│ datetime[μs]        ┆ datetime[μs]        │|
            |╞═════════════════════╪═════════════════════╡|
            |│ 2021-03-01 12:34:00 ┆ 2021-03-01 12:00:00 │|
            |└─────────────────────┴─────────────────────┘|
            └─────────────────────────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.truncate(every)
        )

    def offset_by(self, by: str) -> ExprT:
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
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"datetime": [datetime(2021, 3, 1, 12, 34)]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("datetime").dt.offset_by("1h").alias("datetime_offset_by_1h")
            ... )
            ┌───────────────────────────────────────────────┐
            |              Narwhals DataFrame               |
            |-----------------------------------------------|
            |shape: (1, 2)                                  |
            |┌─────────────────────┬───────────────────────┐|
            |│ datetime            ┆ datetime_offset_by_1h │|
            |│ ---                 ┆ ---                   │|
            |│ datetime[μs]        ┆ datetime[μs]          │|
            |╞═════════════════════╪═══════════════════════╡|
            |│ 2021-03-01 12:34:00 ┆ 2021-03-01 13:34:00   │|
            |└─────────────────────┴───────────────────────┘|
            └───────────────────────────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).dt.offset_by(by)
        )
