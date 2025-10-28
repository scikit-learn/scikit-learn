from __future__ import annotations

import datetime as dt
from typing import TYPE_CHECKING

import polars._reexport as pl
from polars import functions as F
from polars._utils.convert import parse_as_duration_string
from polars._utils.deprecation import deprecate_nonkeyword_arguments, deprecated
from polars._utils.parse import parse_into_expression, parse_into_list_of_expressions
from polars._utils.unstable import unstable
from polars._utils.various import qualified_type_name
from polars._utils.wrap import wrap_expr
from polars.datatypes import DTYPE_TEMPORAL_UNITS, Date, Int32, Int64

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterable

    from polars import Expr
    from polars._typing import (
        Ambiguous,
        EpochTimeUnit,
        IntoExpr,
        IntoExprColumn,
        NonExistent,
        Roll,
        TimeUnit,
    )

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


class ExprDateTimeNameSpace:
    """Namespace for datetime related expressions."""

    _accessor = "dt"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    @unstable()
    @deprecate_nonkeyword_arguments(allowed_args=["self", "n"], version="1.12.0")
    def add_business_days(
        self,
        n: int | IntoExpr,
        week_mask: Iterable[bool] = (True, True, True, True, True, False, False),
        holidays: Iterable[dt.date] = (),
        roll: Roll = "raise",
    ) -> Expr:
        """
        Offset by `n` business days.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.12.0
            Parameters after `n` should now be passed as keyword arguments.

        Parameters
        ----------
        n
            Number of business days to offset by. Can be a single number of an
            expression.
        week_mask
            Which days of the week to count. The default is Monday to Friday.
            If you wanted to count only Monday to Thursday, you would pass
            `(True, True, True, True, False, False, False)`.
        holidays
            Holidays to exclude from the count. The Python package
            `python-holidays <https://github.com/vacanza/python-holidays>`_
            may come in handy here. You can install it with ``pip install holidays``,
            and then, to get all Dutch holidays for years 2020-2024:

            .. code-block:: python

                import holidays

                my_holidays = holidays.country_holidays("NL", years=range(2020, 2025))

            and pass `holidays=my_holidays` when you call `add_business_days`.
        roll
            What to do when the start date lands on a non-business day. Options are:

            - `'raise'`: raise an error
            - `'forward'`: move to the next business day
            - `'backward'`: move to the previous business day

        Returns
        -------
        Expr
            Data type is preserved.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame({"start": [date(2020, 1, 1), date(2020, 1, 2)]})
        >>> df.with_columns(result=pl.col("start").dt.add_business_days(5))
        shape: (2, 2)
        ┌────────────┬────────────┐
        │ start      ┆ result     │
        │ ---        ┆ ---        │
        │ date       ┆ date       │
        ╞════════════╪════════════╡
        │ 2020-01-01 ┆ 2020-01-08 │
        │ 2020-01-02 ┆ 2020-01-09 │
        └────────────┴────────────┘

        You can pass a custom weekend - for example, if you only take Sunday off:

        >>> week_mask = (True, True, True, True, True, True, False)
        >>> df.with_columns(
        ...     result=pl.col("start").dt.add_business_days(5, week_mask=week_mask)
        ... )
        shape: (2, 2)
        ┌────────────┬────────────┐
        │ start      ┆ result     │
        │ ---        ┆ ---        │
        │ date       ┆ date       │
        ╞════════════╪════════════╡
        │ 2020-01-01 ┆ 2020-01-07 │
        │ 2020-01-02 ┆ 2020-01-08 │
        └────────────┴────────────┘

        You can also pass a list of holidays:

        >>> from datetime import date
        >>> holidays = [date(2020, 1, 3), date(2020, 1, 6)]
        >>> df.with_columns(
        ...     result=pl.col("start").dt.add_business_days(5, holidays=holidays)
        ... )
        shape: (2, 2)
        ┌────────────┬────────────┐
        │ start      ┆ result     │
        │ ---        ┆ ---        │
        │ date       ┆ date       │
        ╞════════════╪════════════╡
        │ 2020-01-01 ┆ 2020-01-10 │
        │ 2020-01-02 ┆ 2020-01-13 │
        └────────────┴────────────┘

        Roll all dates forwards to the next business day:

        >>> df = pl.DataFrame({"start": [date(2020, 1, 5), date(2020, 1, 6)]})
        >>> df.with_columns(
        ...     rolled_forwards=pl.col("start").dt.add_business_days(0, roll="forward")
        ... )
        shape: (2, 2)
        ┌────────────┬─────────────────┐
        │ start      ┆ rolled_forwards │
        │ ---        ┆ ---             │
        │ date       ┆ date            │
        ╞════════════╪═════════════════╡
        │ 2020-01-05 ┆ 2020-01-06      │
        │ 2020-01-06 ┆ 2020-01-06      │
        └────────────┴─────────────────┘
        """
        n_pyexpr = parse_into_expression(n)
        unix_epoch = dt.date(1970, 1, 1)
        return wrap_expr(
            self._pyexpr.dt_add_business_days(
                n_pyexpr,
                list(week_mask),
                [(holiday - unix_epoch).days for holiday in holidays],
                roll,
            )
        )

    def truncate(self, every: str | dt.timedelta | Expr) -> Expr:
        """
        Divide the date/datetime range into buckets.

        Each date/datetime is mapped to the start of its bucket using the corresponding
        local datetime. Note that:

        - Weekly buckets start on Monday.
        - All other buckets start on the Unix epoch (1970-01-01).
        - Ambiguous results are localised using the DST offset of the original
          timestamp - for example, truncating `'2022-11-06 01:30:00 CST'` by
          `'1h'` results in `'2022-11-06 01:00:00 CST'`, whereas truncating
          `'2022-11-06 01:30:00 CDT'` by `'1h'` results in
          `'2022-11-06 01:00:00 CDT'`.

        Parameters
        ----------
        every
            The size of each bucket.

        Notes
        -----
        The `every` argument is created with
        the following string language:

        - 1ns   (1 nanosecond)
        - 1us   (1 microsecond)
        - 1ms   (1 millisecond)
        - 1s    (1 second)
        - 1m    (1 minute)
        - 1h    (1 hour)
        - 1d    (1 calendar day)
        - 1w    (1 calendar week)
        - 1mo   (1 calendar month)
        - 1q    (1 calendar quarter)
        - 1y    (1 calendar year)

        By "calendar day", we mean the corresponding time on the next day (which may
        not be 24 hours, due to daylight savings). Similarly for "calendar week",
        "calendar month", "calendar quarter", and "calendar year".

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime`.

        Examples
        --------
        >>> from datetime import timedelta, datetime
        >>> df = (
        ...     pl.datetime_range(
        ...         datetime(2001, 1, 1),
        ...         datetime(2001, 1, 2),
        ...         timedelta(minutes=225),
        ...         eager=True,
        ...     )
        ...     .alias("datetime")
        ...     .to_frame()
        ... )
        >>> df
        shape: (7, 1)
        ┌─────────────────────┐
        │ datetime            │
        │ ---                 │
        │ datetime[μs]        │
        ╞═════════════════════╡
        │ 2001-01-01 00:00:00 │
        │ 2001-01-01 03:45:00 │
        │ 2001-01-01 07:30:00 │
        │ 2001-01-01 11:15:00 │
        │ 2001-01-01 15:00:00 │
        │ 2001-01-01 18:45:00 │
        │ 2001-01-01 22:30:00 │
        └─────────────────────┘
        >>> df.select(pl.col("datetime").dt.truncate("1h"))
        shape: (7, 1)
        ┌─────────────────────┐
        │ datetime            │
        │ ---                 │
        │ datetime[μs]        │
        ╞═════════════════════╡
        │ 2001-01-01 00:00:00 │
        │ 2001-01-01 03:00:00 │
        │ 2001-01-01 07:00:00 │
        │ 2001-01-01 11:00:00 │
        │ 2001-01-01 15:00:00 │
        │ 2001-01-01 18:00:00 │
        │ 2001-01-01 22:00:00 │
        └─────────────────────┘
        >>> truncate_str = df.select(pl.col("datetime").dt.truncate("1h"))
        >>> truncate_td = df.select(pl.col("datetime").dt.truncate(timedelta(hours=1)))
        >>> truncate_str.equals(truncate_td)
        True

        >>> df = (
        ...     pl.datetime_range(
        ...         datetime(2001, 1, 1), datetime(2001, 1, 1, 1), "10m", eager=True
        ...     )
        ...     .alias("datetime")
        ...     .to_frame()
        ... )
        >>> df.select(
        ...     "datetime", pl.col("datetime").dt.truncate("30m").alias("truncate")
        ... )
        shape: (7, 2)
        ┌─────────────────────┬─────────────────────┐
        │ datetime            ┆ truncate            │
        │ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ datetime[μs]        │
        ╞═════════════════════╪═════════════════════╡
        │ 2001-01-01 00:00:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 00:10:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 00:20:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 00:30:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 00:40:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 00:50:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 01:00:00 ┆ 2001-01-01 01:00:00 │
        └─────────────────────┴─────────────────────┘
        """
        if isinstance(every, dt.timedelta):
            every = parse_as_duration_string(every)
        every_pyexpr = parse_into_expression(every, str_as_lit=True)
        return wrap_expr(self._pyexpr.dt_truncate(every_pyexpr))

    def round(self, every: str | dt.timedelta | IntoExprColumn) -> Expr:
        """
        Divide the date/datetime range into buckets.

        - Each date/datetime in the first half of the interval
          is mapped to the start of its bucket.
        - Each date/datetime in the second half of the interval
          is mapped to the end of its bucket.
        - Half-way points are mapped to the start of their bucket.

        Ambiguous results are localised using the DST offset of the original timestamp -
        for example, rounding `'2022-11-06 01:20:00 CST'` by `'1h'` results in
        `'2022-11-06 01:00:00 CST'`, whereas rounding `'2022-11-06 01:20:00 CDT'` by
        `'1h'` results in `'2022-11-06 01:00:00 CDT'`.

        Parameters
        ----------
        every
            Every interval start and period length

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        The `every` argument is created with
        the following small string formatting language:

        - 1ns   (1 nanosecond)
        - 1us   (1 microsecond)
        - 1ms   (1 millisecond)
        - 1s    (1 second)
        - 1m    (1 minute)
        - 1h    (1 hour)
        - 1d    (1 calendar day)
        - 1w    (1 calendar week)
        - 1mo   (1 calendar month)
        - 1q    (1 calendar quarter)
        - 1y    (1 calendar year)

        By "calendar day", we mean the corresponding time on the next day (which may
        not be 24 hours, due to daylight savings). Similarly for "calendar week",
        "calendar month", "calendar quarter", and "calendar year".

        Examples
        --------
        >>> from datetime import timedelta, datetime
        >>> df = (
        ...     pl.datetime_range(
        ...         datetime(2001, 1, 1),
        ...         datetime(2001, 1, 2),
        ...         timedelta(minutes=225),
        ...         eager=True,
        ...     )
        ...     .alias("datetime")
        ...     .to_frame()
        ... )
        >>> df.with_columns(pl.col("datetime").dt.round("1h").alias("round"))
        shape: (7, 2)
        ┌─────────────────────┬─────────────────────┐
        │ datetime            ┆ round               │
        │ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ datetime[μs]        │
        ╞═════════════════════╪═════════════════════╡
        │ 2001-01-01 00:00:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 03:45:00 ┆ 2001-01-01 04:00:00 │
        │ 2001-01-01 07:30:00 ┆ 2001-01-01 08:00:00 │
        │ 2001-01-01 11:15:00 ┆ 2001-01-01 11:00:00 │
        │ 2001-01-01 15:00:00 ┆ 2001-01-01 15:00:00 │
        │ 2001-01-01 18:45:00 ┆ 2001-01-01 19:00:00 │
        │ 2001-01-01 22:30:00 ┆ 2001-01-01 23:00:00 │
        └─────────────────────┴─────────────────────┘

        >>> df = (
        ...     pl.datetime_range(
        ...         datetime(2001, 1, 1), datetime(2001, 1, 1, 1), "10m", eager=True
        ...     )
        ...     .alias("datetime")
        ...     .to_frame()
        ... )
        >>> df.with_columns(pl.col("datetime").dt.round("30m").alias("round"))
        shape: (7, 2)
        ┌─────────────────────┬─────────────────────┐
        │ datetime            ┆ round               │
        │ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ datetime[μs]        │
        ╞═════════════════════╪═════════════════════╡
        │ 2001-01-01 00:00:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 00:10:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-01 00:20:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 00:30:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 00:40:00 ┆ 2001-01-01 00:30:00 │
        │ 2001-01-01 00:50:00 ┆ 2001-01-01 01:00:00 │
        │ 2001-01-01 01:00:00 ┆ 2001-01-01 01:00:00 │
        └─────────────────────┴─────────────────────┘
        """
        if isinstance(every, dt.timedelta):
            every = parse_as_duration_string(every)
        every_pyexpr = parse_into_expression(every, str_as_lit=True)
        return wrap_expr(self._pyexpr.dt_round(every_pyexpr))

    def replace(
        self,
        *,
        year: int | IntoExpr | None = None,
        month: int | IntoExpr | None = None,
        day: int | IntoExpr | None = None,
        hour: int | IntoExpr | None = None,
        minute: int | IntoExpr | None = None,
        second: int | IntoExpr | None = None,
        microsecond: int | IntoExpr | None = None,
        ambiguous: Ambiguous | Expr = "raise",
    ) -> Expr:
        """
        Replace time unit.

        Parameters
        ----------
        year
            Column or literal.
        month
            Column or literal, ranging from 1-12.
        day
            Column or literal, ranging from 1-31.
        hour
            Column or literal, ranging from 0-23.
        minute
            Column or literal, ranging from 0-59.
        second
            Column or literal, ranging from 0-59.
        microsecond
            Column or literal, ranging from 0-999999.
        ambiguous
            Determine how to deal with ambiguous datetimes:

            - `'raise'` (default): raise
            - `'earliest'`: use the earliest datetime
            - `'latest'`: use the latest datetime
            - `'null'`: set to null

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime` with the
            specified time units replaced.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": [date(2024, 4, 1), date(2025, 3, 16)],
        ...         "new_day": [10, 15],
        ...     }
        ... )
        >>> df.with_columns(pl.col("date").dt.replace(day="new_day").alias("replaced"))
        shape: (2, 3)
        ┌────────────┬─────────┬────────────┐
        │ date       ┆ new_day ┆ replaced   │
        │ ---        ┆ ---     ┆ ---        │
        │ date       ┆ i64     ┆ date       │
        ╞════════════╪═════════╪════════════╡
        │ 2024-04-01 ┆ 10      ┆ 2024-04-10 │
        │ 2025-03-16 ┆ 15      ┆ 2025-03-15 │
        └────────────┴─────────┴────────────┘
        >>> df.with_columns(pl.col("date").dt.replace(year=1800).alias("replaced"))
        shape: (2, 3)
        ┌────────────┬─────────┬────────────┐
        │ date       ┆ new_day ┆ replaced   │
        │ ---        ┆ ---     ┆ ---        │
        │ date       ┆ i64     ┆ date       │
        ╞════════════╪═════════╪════════════╡
        │ 2024-04-01 ┆ 10      ┆ 1800-04-01 │
        │ 2025-03-16 ┆ 15      ┆ 1800-03-16 │
        └────────────┴─────────┴────────────┘
        """
        (
            day_pyexpr,
            month_pyexpr,
            year_pyexpr,
            hour_pyexpr,
            minute_pyexpr,
            second_pyexpr,
            microsecond_pyexpr,
        ) = parse_into_list_of_expressions(
            day, month, year, hour, minute, second, microsecond
        )
        ambiguous_expr = parse_into_expression(ambiguous, str_as_lit=True)
        return wrap_expr(
            self._pyexpr.dt_replace(
                year_pyexpr,
                month_pyexpr,
                day_pyexpr,
                hour_pyexpr,
                minute_pyexpr,
                second_pyexpr,
                microsecond_pyexpr,
                ambiguous_expr,
            )
        )

    def combine(self, time: dt.time | Expr, time_unit: TimeUnit = "us") -> Expr:
        """
        Create a naive Datetime from an existing Date/Datetime expression and a Time.

        If the underlying expression is a Datetime then its time component is replaced,
        and if it is a Date then a new Datetime is created by combining the two values.

        Parameters
        ----------
        time
            A python time literal or polars expression/column that resolves to a time.
        time_unit : {'ns', 'us', 'ms'}
            Unit of time.

        Examples
        --------
        >>> from datetime import datetime, date, time
        >>> df = pl.DataFrame(
        ...     {
        ...         "dtm": [
        ...             datetime(2022, 12, 31, 10, 30, 45),
        ...             datetime(2023, 7, 5, 23, 59, 59),
        ...         ],
        ...         "dt": [date(2022, 10, 10), date(2022, 7, 5)],
        ...         "tm": [time(1, 2, 3, 456000), time(7, 8, 9, 101000)],
        ...     }
        ... )
        >>> df
        shape: (2, 3)
        ┌─────────────────────┬────────────┬──────────────┐
        │ dtm                 ┆ dt         ┆ tm           │
        │ ---                 ┆ ---        ┆ ---          │
        │ datetime[μs]        ┆ date       ┆ time         │
        ╞═════════════════════╪════════════╪══════════════╡
        │ 2022-12-31 10:30:45 ┆ 2022-10-10 ┆ 01:02:03.456 │
        │ 2023-07-05 23:59:59 ┆ 2022-07-05 ┆ 07:08:09.101 │
        └─────────────────────┴────────────┴──────────────┘
        >>> df.select(
        ...     [
        ...         pl.col("dtm").dt.combine(pl.col("tm")).alias("d1"),
        ...         pl.col("dt").dt.combine(pl.col("tm")).alias("d2"),
        ...         pl.col("dt").dt.combine(time(4, 5, 6)).alias("d3"),
        ...     ]
        ... )
        shape: (2, 3)
        ┌─────────────────────────┬─────────────────────────┬─────────────────────┐
        │ d1                      ┆ d2                      ┆ d3                  │
        │ ---                     ┆ ---                     ┆ ---                 │
        │ datetime[μs]            ┆ datetime[μs]            ┆ datetime[μs]        │
        ╞═════════════════════════╪═════════════════════════╪═════════════════════╡
        │ 2022-12-31 01:02:03.456 ┆ 2022-10-10 01:02:03.456 ┆ 2022-10-10 04:05:06 │
        │ 2023-07-05 07:08:09.101 ┆ 2022-07-05 07:08:09.101 ┆ 2022-07-05 04:05:06 │
        └─────────────────────────┴─────────────────────────┴─────────────────────┘
        """
        if not isinstance(time, (dt.time, pl.Expr)):
            msg = f"expected 'time' to be a Python time or Polars expression, found {qualified_type_name(time)!r}"
            raise TypeError(msg)
        time_pyexpr = parse_into_expression(time)
        return wrap_expr(self._pyexpr.dt_combine(time_pyexpr, time_unit))

    def to_string(self, format: str | None = None) -> Expr:
        """
        Convert a Date/Time/Datetime column into a String column with the given format.

        .. versionchanged:: 1.15.0
            Added support for the use of "iso:strict" as a format string.
        .. versionchanged:: 1.14.0
            Added support for the `Duration` dtype, and use of "iso" as a format string.

        Parameters
        ----------
        format
            * Format to use, refer to the `chrono strftime documentation
              <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
              for specification. Example: `"%y-%m-%d"`.

            * If no format is provided, the appropriate ISO format for the underlying
              data type is used. This can be made explicit by passing `"iso"` or
              `"iso:strict"` as the format string (see notes below for details).

        Notes
        -----
        * Similar to `cast(pl.String)`, but this method allows you to customize
          the formatting of the resulting string; if no format is provided, the
          appropriate ISO format for the underlying data type is used.

        * Datetime dtype expressions distinguish between "iso" and "iso:strict"
          format strings. The difference is in the inclusion of a "T" separator
          between the date and time components ("iso" results in ISO compliant
          date and time components, separated with a space; "iso:strict" returns
          the same components separated with a "T"). All other temporal types
          return the same value for both format strings.

        * Duration dtype expressions cannot be formatted with `strftime`. Instead,
          only "iso" and "polars" are supported as format strings. The "iso" format
          string results in ISO8601 duration string output, and "polars" results
          in the same form seen in the frame `repr`.

        Examples
        --------
        >>> from datetime import datetime, date, timedelta, time
        >>> df = pl.DataFrame(
        ...     {
        ...         "dt": [
        ...             date(1999, 3, 1),
        ...             date(2020, 5, 3),
        ...             date(2077, 7, 5),
        ...         ],
        ...         "dtm": [
        ...             datetime(1980, 8, 10, 0, 10, 20),
        ...             datetime(2010, 10, 20, 8, 25, 35),
        ...             datetime(2040, 12, 30, 16, 40, 50),
        ...         ],
        ...         "tm": [
        ...             time(1, 2, 3, 456789),
        ...             time(23, 59, 9, 101),
        ...             time(0, 0, 0, 100),
        ...         ],
        ...         "td": [
        ...             timedelta(days=-1, seconds=-42),
        ...             timedelta(days=14, hours=-10, microseconds=100),
        ...             timedelta(seconds=0),
        ...         ],
        ...     }
        ... )

        Default format for temporal dtypes is ISO8601:

        >>> import polars.selectors as cs
        >>> df.select(cs.temporal().dt.to_string().name.prefix("s_"))
        shape: (3, 4)
        ┌────────────┬────────────────────────────┬─────────────────┬─────────────────┐
        │ s_dt       ┆ s_dtm                      ┆ s_tm            ┆ s_td            │
        │ ---        ┆ ---                        ┆ ---             ┆ ---             │
        │ str        ┆ str                        ┆ str             ┆ str             │
        ╞════════════╪════════════════════════════╪═════════════════╪═════════════════╡
        │ 1999-03-01 ┆ 1980-08-10 00:10:20.000000 ┆ 01:02:03.456789 ┆ -P1DT42S        │
        │ 2020-05-03 ┆ 2010-10-20 08:25:35.000000 ┆ 23:59:09.000101 ┆ P13DT14H0.0001S │
        │ 2077-07-05 ┆ 2040-12-30 16:40:50.000000 ┆ 00:00:00.000100 ┆ PT0S            │
        └────────────┴────────────────────────────┴─────────────────┴─────────────────┘

        For `Datetime` specifically you can choose between "iso" (where the date and
        time components are ISO, separated by a space) and "iso:strict" (where these
        components are separated by a "T"):

        >>> df.select(
        ...     pl.col("dtm").dt.to_string("iso").alias("dtm_iso"),
        ...     pl.col("dtm").dt.to_string("iso:strict").alias("dtm_iso_strict"),
        ... )
        shape: (3, 2)
        ┌────────────────────────────┬────────────────────────────┐
        │ dtm_iso                    ┆ dtm_iso_strict             │
        │ ---                        ┆ ---                        │
        │ str                        ┆ str                        │
        ╞════════════════════════════╪════════════════════════════╡
        │ 1980-08-10 00:10:20.000000 ┆ 1980-08-10T00:10:20.000000 │
        │ 2010-10-20 08:25:35.000000 ┆ 2010-10-20T08:25:35.000000 │
        │ 2040-12-30 16:40:50.000000 ┆ 2040-12-30T16:40:50.000000 │
        └────────────────────────────┴────────────────────────────┘

        All temporal types (aside from `Duration`) support strftime formatting:

        >>> df.select(
        ...     pl.col("dtm"),
        ...     s_dtm=pl.col("dtm").dt.to_string("%Y/%m/%d (%H.%M.%S)"),
        ... )
        shape: (3, 2)
        ┌─────────────────────┬───────────────────────┐
        │ dtm                 ┆ s_dtm                 │
        │ ---                 ┆ ---                   │
        │ datetime[μs]        ┆ str                   │
        ╞═════════════════════╪═══════════════════════╡
        │ 1980-08-10 00:10:20 ┆ 1980/08/10 (00.10.20) │
        │ 2010-10-20 08:25:35 ┆ 2010/10/20 (08.25.35) │
        │ 2040-12-30 16:40:50 ┆ 2040/12/30 (16.40.50) │
        └─────────────────────┴───────────────────────┘

        The Polars Duration string format (as seen in the frame repr) is also available:

        >>> df.select(
        ...     pl.col("td"),
        ...     s_td=pl.col("td").dt.to_string("polars"),
        ... )
        shape: (3, 2)
        ┌───────────────┬───────────────┐
        │ td            ┆ s_td          │
        │ ---           ┆ ---           │
        │ duration[μs]  ┆ str           │
        ╞═══════════════╪═══════════════╡
        │ -1d -42s      ┆ -1d -42s      │
        │ 13d 14h 100µs ┆ 13d 14h 100µs │
        │ 0µs           ┆ 0µs           │
        └───────────────┴───────────────┘

        If you're interested in extracting the day or month names, you can use
        the `'%A'` and `'%B'` strftime specifiers:

        >>> df.select(
        ...     pl.col("dt"),
        ...     day_name=pl.col("dtm").dt.to_string("%A"),
        ...     month_name=pl.col("dtm").dt.to_string("%B"),
        ... )
        shape: (3, 3)
        ┌────────────┬───────────┬────────────┐
        │ dt         ┆ day_name  ┆ month_name │
        │ ---        ┆ ---       ┆ ---        │
        │ date       ┆ str       ┆ str        │
        ╞════════════╪═══════════╪════════════╡
        │ 1999-03-01 ┆ Sunday    ┆ August     │
        │ 2020-05-03 ┆ Wednesday ┆ October    │
        │ 2077-07-05 ┆ Sunday    ┆ December   │
        └────────────┴───────────┴────────────┘
        """
        if format is None:
            format = "iso"
        return wrap_expr(self._pyexpr.dt_to_string(format))

    def strftime(self, format: str) -> Expr:
        """
        Convert a Date/Time/Datetime column into a String column with the given format.

        Similar to `cast(pl.String)`, but this method allows you to customize the
        formatting of the resulting string.

        Alias for :func:`to_string`.

        Parameters
        ----------
        format
            Format to use, refer to the `chrono strftime documentation
            <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            for specification. Example: `"%y-%m-%d"`.

        See Also
        --------
        to_string : The identical expression for which `strftime` is an alias.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(2020, 3, 1),
        ...             datetime(2020, 4, 1),
        ...             datetime(2020, 5, 1),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime")
        ...     .dt.strftime("%Y/%m/%d %H:%M:%S")
        ...     .alias("datetime_string")
        ... )
        shape: (3, 2)
        ┌─────────────────────┬─────────────────────┐
        │ datetime            ┆ datetime_string     │
        │ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ str                 │
        ╞═════════════════════╪═════════════════════╡
        │ 2020-03-01 00:00:00 ┆ 2020/03/01 00:00:00 │
        │ 2020-04-01 00:00:00 ┆ 2020/04/01 00:00:00 │
        │ 2020-05-01 00:00:00 ┆ 2020/05/01 00:00:00 │
        └─────────────────────┴─────────────────────┘

        If you're interested in the day name / month name, you can use
        `'%A'` / `'%B'`:

        >>> df.with_columns(
        ...     day_name=pl.col("datetime").dt.strftime("%A"),
        ...     month_name=pl.col("datetime").dt.strftime("%B"),
        ... )
        shape: (3, 3)
        ┌─────────────────────┬───────────┬────────────┐
        │ datetime            ┆ day_name  ┆ month_name │
        │ ---                 ┆ ---       ┆ ---        │
        │ datetime[μs]        ┆ str       ┆ str        │
        ╞═════════════════════╪═══════════╪════════════╡
        │ 2020-03-01 00:00:00 ┆ Sunday    ┆ March      │
        │ 2020-04-01 00:00:00 ┆ Wednesday ┆ April      │
        │ 2020-05-01 00:00:00 ┆ Friday    ┆ May        │
        └─────────────────────┴───────────┴────────────┘
        """
        return self.to_string(format)

    def millennium(self) -> Expr:
        """
        Extract the millennium from underlying representation.

        Applies to Date and Datetime columns.

        Returns the millennium number in the calendar date.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": [
        ...             date(999, 12, 31),
        ...             date(1897, 5, 7),
        ...             date(2000, 1, 1),
        ...             date(2001, 7, 5),
        ...             date(3002, 10, 20),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(mlnm=pl.col("date").dt.millennium())
        shape: (5, 2)
        ┌────────────┬──────┐
        │ date       ┆ mlnm │
        │ ---        ┆ ---  │
        │ date       ┆ i32  │
        ╞════════════╪══════╡
        │ 0999-12-31 ┆ 1    │
        │ 1897-05-07 ┆ 2    │
        │ 2000-01-01 ┆ 2    │
        │ 2001-07-05 ┆ 3    │
        │ 3002-10-20 ┆ 4    │
        └────────────┴──────┘
        """
        return wrap_expr(self._pyexpr.dt_millennium())

    def century(self) -> Expr:
        """
        Extract the century from underlying representation.

        Applies to Date and Datetime columns.

        Returns the century number in the calendar date.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": [
        ...             date(999, 12, 31),
        ...             date(1897, 5, 7),
        ...             date(2000, 1, 1),
        ...             date(2001, 7, 5),
        ...             date(3002, 10, 20),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(cent=pl.col("date").dt.century())
        shape: (5, 2)
        ┌────────────┬──────┐
        │ date       ┆ cent │
        │ ---        ┆ ---  │
        │ date       ┆ i32  │
        ╞════════════╪══════╡
        │ 0999-12-31 ┆ 10   │
        │ 1897-05-07 ┆ 19   │
        │ 2000-01-01 ┆ 20   │
        │ 2001-07-05 ┆ 21   │
        │ 3002-10-20 ┆ 31   │
        └────────────┴──────┘
        """
        return wrap_expr(self._pyexpr.dt_century())

    def year(self) -> Expr:
        """
        Extract year from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the year number in the calendar date.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(1977, 1, 1), date(1978, 1, 1), date(1979, 1, 1)]}
        ... )
        >>> df.with_columns(
        ...     calendar_year=pl.col("date").dt.year(),
        ...     iso_year=pl.col("date").dt.iso_year(),
        ... )
        shape: (3, 3)
        ┌────────────┬───────────────┬──────────┐
        │ date       ┆ calendar_year ┆ iso_year │
        │ ---        ┆ ---           ┆ ---      │
        │ date       ┆ i32           ┆ i32      │
        ╞════════════╪═══════════════╪══════════╡
        │ 1977-01-01 ┆ 1977          ┆ 1976     │
        │ 1978-01-01 ┆ 1978          ┆ 1977     │
        │ 1979-01-01 ┆ 1979          ┆ 1979     │
        └────────────┴───────────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.dt_year())

    @unstable()
    def is_business_day(
        self,
        *,
        week_mask: Iterable[bool] = (True, True, True, True, True, False, False),
        holidays: Iterable[dt.date] = (),
    ) -> Expr:
        """
        Determine whether each day lands on a business day.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        week_mask
            Which days of the week to count. The default is Monday to Friday.
            If you wanted to count only Monday to Thursday, you would pass
            `(True, True, True, True, False, False, False)`.
        holidays
            Holidays to exclude from the count. The Python package
            `python-holidays <https://github.com/vacanza/python-holidays>`_
            may come in handy here. You can install it with ``pip install holidays``,
            and then, to get all Dutch holidays for years 2020-2024:

            .. code-block:: python

                import holidays

                my_holidays = holidays.country_holidays("NL", years=range(2020, 2025))

            and pass `holidays=my_holidays` when you call `is_business_day`.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame({"start": [date(2020, 1, 3), date(2020, 1, 5)]})
        >>> df.with_columns(is_business_day=pl.col("start").dt.is_business_day())
        shape: (2, 2)
        ┌────────────┬─────────────────┐
        │ start      ┆ is_business_day │
        │ ---        ┆ ---             │
        │ date       ┆ bool            │
        ╞════════════╪═════════════════╡
        │ 2020-01-03 ┆ true            │
        │ 2020-01-05 ┆ false           │
        └────────────┴─────────────────┘

        You can pass a custom weekend - for example, if you only take Sunday off:

        >>> week_mask = (True, True, True, True, True, True, False)
        >>> df.with_columns(
        ...     is_business_day=pl.col("start").dt.is_business_day(week_mask=week_mask)
        ... )
        shape: (2, 2)
        ┌────────────┬─────────────────┐
        │ start      ┆ is_business_day │
        │ ---        ┆ ---             │
        │ date       ┆ bool            │
        ╞════════════╪═════════════════╡
        │ 2020-01-03 ┆ true            │
        │ 2020-01-05 ┆ false           │
        └────────────┴─────────────────┘

        You can also pass a list of holidays:

        >>> from datetime import date
        >>> holidays = [date(2020, 1, 3), date(2020, 1, 6)]
        >>> df.with_columns(
        ...     is_business_day=pl.col("start").dt.is_business_day(holidays=holidays)
        ... )
        shape: (2, 2)
        ┌────────────┬─────────────────┐
        │ start      ┆ is_business_day │
        │ ---        ┆ ---             │
        │ date       ┆ bool            │
        ╞════════════╪═════════════════╡
        │ 2020-01-03 ┆ false           │
        │ 2020-01-05 ┆ false           │
        └────────────┴─────────────────┘
        """
        unix_epoch = dt.date(1970, 1, 1)
        return wrap_expr(
            self._pyexpr.dt_is_business_day(
                list(week_mask),
                [(holiday - unix_epoch).days for holiday in holidays],
            )
        )

    def is_leap_year(self) -> Expr:
        """
        Determine whether the year of the underlying date is a leap year.

        Applies to Date and Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(2000, 1, 1), date(2001, 1, 1), date(2002, 1, 1)]}
        ... )
        >>> df.with_columns(
        ...     leap_year=pl.col("date").dt.is_leap_year(),
        ... )
        shape: (3, 2)
        ┌────────────┬───────────┐
        │ date       ┆ leap_year │
        │ ---        ┆ ---       │
        │ date       ┆ bool      │
        ╞════════════╪═══════════╡
        │ 2000-01-01 ┆ true      │
        │ 2001-01-01 ┆ false     │
        │ 2002-01-01 ┆ false     │
        └────────────┴───────────┘
        """
        return wrap_expr(self._pyexpr.dt_is_leap_year())

    def iso_year(self) -> Expr:
        """
        Extract ISO year from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the year number in the ISO standard.
        This may not correspond with the calendar year.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(1977, 1, 1), date(1978, 1, 1), date(1979, 1, 1)]}
        ... )
        >>> df.select(
        ...     "date",
        ...     pl.col("date").dt.year().alias("calendar_year"),
        ...     pl.col("date").dt.iso_year().alias("iso_year"),
        ... )
        shape: (3, 3)
        ┌────────────┬───────────────┬──────────┐
        │ date       ┆ calendar_year ┆ iso_year │
        │ ---        ┆ ---           ┆ ---      │
        │ date       ┆ i32           ┆ i32      │
        ╞════════════╪═══════════════╪══════════╡
        │ 1977-01-01 ┆ 1977          ┆ 1976     │
        │ 1978-01-01 ┆ 1978          ┆ 1977     │
        │ 1979-01-01 ┆ 1979          ┆ 1979     │
        └────────────┴───────────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.dt_iso_year())

    def quarter(self) -> Expr:
        """
        Extract quarter from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the quarter ranging from 1 to 4.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(2001, 1, 1), date(2001, 6, 30), date(2001, 12, 27)]}
        ... )
        >>> df.with_columns(pl.col("date").dt.quarter().alias("quarter"))
        shape: (3, 2)
        ┌────────────┬─────────┐
        │ date       ┆ quarter │
        │ ---        ┆ ---     │
        │ date       ┆ i8      │
        ╞════════════╪═════════╡
        │ 2001-01-01 ┆ 1       │
        │ 2001-06-30 ┆ 2       │
        │ 2001-12-27 ┆ 4       │
        └────────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.dt_quarter())

    def month(self) -> Expr:
        """
        Extract month from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the month number starting from 1.
        The return value ranges from 1 to 12.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(2001, 1, 1), date(2001, 6, 30), date(2001, 12, 27)]}
        ... )
        >>> df.with_columns(pl.col("date").dt.month().alias("month"))
        shape: (3, 2)
        ┌────────────┬───────┐
        │ date       ┆ month │
        │ ---        ┆ ---   │
        │ date       ┆ i8    │
        ╞════════════╪═══════╡
        │ 2001-01-01 ┆ 1     │
        │ 2001-06-30 ┆ 6     │
        │ 2001-12-27 ┆ 12    │
        └────────────┴───────┘
        """
        return wrap_expr(self._pyexpr.dt_month())

    def days_in_month(self) -> Expr:
        """
        Extract the number of days in the month from the underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the number of days in the month.
        The return value ranges from 28 to 31.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        See Also
        --------
        month
        is_leap_year

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(2001, 1, 1), date(2001, 2, 1), date(2000, 2, 1)]}
        ... )
        >>> df.with_columns(pl.col("date").dt.days_in_month().alias("days_in_month"))
        shape: (3, 2)
        ┌────────────┬───────────────┐
        │ date       ┆ days_in_month │
        │ ---        ┆ ---           │
        │ date       ┆ i8            │
        ╞════════════╪═══════════════╡
        │ 2001-01-01 ┆ 31            │
        │ 2001-02-01 ┆ 28            │
        │ 2000-02-01 ┆ 29            │
        └────────────┴───────────────┘
        """
        return wrap_expr(self._pyexpr.dt_days_in_month())

    def week(self) -> Expr:
        """
        Extract the week from the underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the ISO week number starting from 1.
        The return value ranges from 1 to 53. (The last week of year differs by years.)

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {"date": [date(2001, 1, 1), date(2001, 6, 30), date(2001, 12, 27)]}
        ... )
        >>> df.with_columns(pl.col("date").dt.week().alias("week"))
        shape: (3, 2)
        ┌────────────┬──────┐
        │ date       ┆ week │
        │ ---        ┆ ---  │
        │ date       ┆ i8   │
        ╞════════════╪══════╡
        │ 2001-01-01 ┆ 1    │
        │ 2001-06-30 ┆ 26   │
        │ 2001-12-27 ┆ 52   │
        └────────────┴──────┘
        """
        return wrap_expr(self._pyexpr.dt_week())

    def weekday(self) -> Expr:
        """
        Extract the week day from the underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the ISO weekday number where monday = 1 and sunday = 7

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        See Also
        --------
        day
        ordinal_day

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.date_range(
        ...             date(2001, 12, 22), date(2001, 12, 25), eager=True
        ...         )
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("date").dt.weekday().alias("weekday"),
        ...     pl.col("date").dt.day().alias("day_of_month"),
        ...     pl.col("date").dt.ordinal_day().alias("day_of_year"),
        ... )
        shape: (4, 4)
        ┌────────────┬─────────┬──────────────┬─────────────┐
        │ date       ┆ weekday ┆ day_of_month ┆ day_of_year │
        │ ---        ┆ ---     ┆ ---          ┆ ---         │
        │ date       ┆ i8      ┆ i8           ┆ i16         │
        ╞════════════╪═════════╪══════════════╪═════════════╡
        │ 2001-12-22 ┆ 6       ┆ 22           ┆ 356         │
        │ 2001-12-23 ┆ 7       ┆ 23           ┆ 357         │
        │ 2001-12-24 ┆ 1       ┆ 24           ┆ 358         │
        │ 2001-12-25 ┆ 2       ┆ 25           ┆ 359         │
        └────────────┴─────────┴──────────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_weekday())

    def day(self) -> Expr:
        """
        Extract day from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the day of month starting from 1.
        The return value ranges from 1 to 31. (The last day of month differs by months.)

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        See Also
        --------
        weekday
        ordinal_day

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.date_range(
        ...             date(2001, 12, 22), date(2001, 12, 25), eager=True
        ...         )
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("date").dt.weekday().alias("weekday"),
        ...     pl.col("date").dt.day().alias("day_of_month"),
        ...     pl.col("date").dt.ordinal_day().alias("day_of_year"),
        ... )
        shape: (4, 4)
        ┌────────────┬─────────┬──────────────┬─────────────┐
        │ date       ┆ weekday ┆ day_of_month ┆ day_of_year │
        │ ---        ┆ ---     ┆ ---          ┆ ---         │
        │ date       ┆ i8      ┆ i8           ┆ i16         │
        ╞════════════╪═════════╪══════════════╪═════════════╡
        │ 2001-12-22 ┆ 6       ┆ 22           ┆ 356         │
        │ 2001-12-23 ┆ 7       ┆ 23           ┆ 357         │
        │ 2001-12-24 ┆ 1       ┆ 24           ┆ 358         │
        │ 2001-12-25 ┆ 2       ┆ 25           ┆ 359         │
        └────────────┴─────────┴──────────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_day())

    def ordinal_day(self) -> Expr:
        """
        Extract ordinal day from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the day of year starting from 1.
        The return value ranges from 1 to 366. (The last day of year differs by years.)

        Returns
        -------
        Expr
            Expression of data type :class:`Int16`.

        See Also
        --------
        weekday
        day

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.date_range(
        ...             date(2001, 12, 22), date(2001, 12, 25), eager=True
        ...         )
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("date").dt.weekday().alias("weekday"),
        ...     pl.col("date").dt.day().alias("day_of_month"),
        ...     pl.col("date").dt.ordinal_day().alias("day_of_year"),
        ... )
        shape: (4, 4)
        ┌────────────┬─────────┬──────────────┬─────────────┐
        │ date       ┆ weekday ┆ day_of_month ┆ day_of_year │
        │ ---        ┆ ---     ┆ ---          ┆ ---         │
        │ date       ┆ i8      ┆ i8           ┆ i16         │
        ╞════════════╪═════════╪══════════════╪═════════════╡
        │ 2001-12-22 ┆ 6       ┆ 22           ┆ 356         │
        │ 2001-12-23 ┆ 7       ┆ 23           ┆ 357         │
        │ 2001-12-24 ┆ 1       ┆ 24           ┆ 358         │
        │ 2001-12-25 ┆ 2       ┆ 25           ┆ 359         │
        └────────────┴─────────┴──────────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_ordinal_day())

    def time(self) -> Expr:
        """
        Extract time.

        Applies to Datetime columns only; fails on Date.

        Returns
        -------
        Expr
            Expression of data type :class:`Time`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(pl.col("datetime").dt.time().alias("time"))
        shape: (3, 2)
        ┌─────────────────────────┬──────────────┐
        │ datetime                ┆ time         │
        │ ---                     ┆ ---          │
        │ datetime[μs]            ┆ time         │
        ╞═════════════════════════╪══════════════╡
        │ 1978-01-01 01:01:01     ┆ 01:01:01     │
        │ 2024-10-13 05:30:14.500 ┆ 05:30:14.500 │
        │ 2065-01-01 10:20:30.060 ┆ 10:20:30.060 │
        └─────────────────────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.dt_time())

    def date(self) -> Expr:
        """
        Extract date from date(time).

        Applies to Date and Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Date`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(pl.col("datetime").dt.date().alias("date"))
        shape: (3, 2)
        ┌─────────────────────────┬────────────┐
        │ datetime                ┆ date       │
        │ ---                     ┆ ---        │
        │ datetime[μs]            ┆ date       │
        ╞═════════════════════════╪════════════╡
        │ 1978-01-01 01:01:01     ┆ 1978-01-01 │
        │ 2024-10-13 05:30:14.500 ┆ 2024-10-13 │
        │ 2065-01-01 10:20:30.060 ┆ 2065-01-01 │
        └─────────────────────────┴────────────┘
        """
        return wrap_expr(self._pyexpr.dt_date())

    @deprecated(
        "`dt.datetime` is deprecated; use `dt.replace_time_zone(None)` instead."
    )
    def datetime(self) -> Expr:
        """
        Return datetime.

        .. deprecated:: 0.20.4
            Use the `dt.replace_time_zone(None)` method instead.

        Applies to Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Datetime`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime UTC": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     },
        ...     schema={"datetime UTC": pl.Datetime(time_zone="UTC")},
        ... )
        >>> df.with_columns(  # doctest: +SKIP
        ...     pl.col("datetime UTC").dt.datetime().alias("datetime (no timezone)"),
        ... )
        shape: (3, 2)
        ┌─────────────────────────────┬─────────────────────────┐
        │ datetime UTC                ┆ datetime (no timezone)  │
        │ ---                         ┆ ---                     │
        │ datetime[μs, UTC]           ┆ datetime[μs]            │
        ╞═════════════════════════════╪═════════════════════════╡
        │ 1978-01-01 01:01:01 UTC     ┆ 1978-01-01 01:01:01     │
        │ 2024-10-13 05:30:14.500 UTC ┆ 2024-10-13 05:30:14.500 │
        │ 2065-01-01 10:20:30.060 UTC ┆ 2065-01-01 10:20:30.060 │
        └─────────────────────────────┴─────────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_datetime())

    def hour(self) -> Expr:
        """
        Extract hour from underlying DateTime representation.

        Applies to Datetime columns.

        Returns the hour number from 0 to 23.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ...     pl.col("datetime").dt.millisecond().alias("millisecond"),
        ... )
        shape: (3, 5)
        ┌─────────────────────────┬──────┬────────┬────────┬─────────────┐
        │ datetime                ┆ hour ┆ minute ┆ second ┆ millisecond │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    ┆ ---         │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     ┆ i32         │
        ╞═════════════════════════╪══════╪════════╪════════╪═════════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      ┆ 0           │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     ┆ 500         │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     ┆ 60          │
        └─────────────────────────┴──────┴────────┴────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_hour())

    def minute(self) -> Expr:
        """
        Extract minutes from underlying DateTime representation.

        Applies to Datetime columns.

        Returns the minute number from 0 to 59.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ...     pl.col("datetime").dt.millisecond().alias("millisecond"),
        ... )
        shape: (3, 5)
        ┌─────────────────────────┬──────┬────────┬────────┬─────────────┐
        │ datetime                ┆ hour ┆ minute ┆ second ┆ millisecond │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    ┆ ---         │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     ┆ i32         │
        ╞═════════════════════════╪══════╪════════╪════════╪═════════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      ┆ 0           │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     ┆ 500         │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     ┆ 60          │
        └─────────────────────────┴──────┴────────┴────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_minute())

    def second(self, *, fractional: bool = False) -> Expr:
        """
        Extract seconds from underlying DateTime representation.

        Applies to Datetime columns.

        Returns the integer second number from 0 to 59, or a floating
        point number from 0 < 60 if `fractional=True` that includes
        any milli/micro/nanosecond component.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the second.

        Returns
        -------
        Expr
            Expression of data type :class:`Int8` or :class:`Float64`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ... )
        shape: (3, 4)
        ┌─────────────────────────┬──────┬────────┬────────┐
        │ datetime                ┆ hour ┆ minute ┆ second │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     │
        ╞═════════════════════════╪══════╪════════╪════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     │
        └─────────────────────────┴──────┴────────┴────────┘
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second(fractional=True).alias("second"),
        ... )
        shape: (3, 4)
        ┌─────────────────────────┬──────┬────────┬────────┐
        │ datetime                ┆ hour ┆ minute ┆ second │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ f64    │
        ╞═════════════════════════╪══════╪════════╪════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1.0    │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14.5   │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30.06  │
        └─────────────────────────┴──────┴────────┴────────┘
        """
        sec = wrap_expr(self._pyexpr.dt_second())
        return (
            sec + (wrap_expr(self._pyexpr.dt_nanosecond()) / F.lit(1_000_000_000.0))
            if fractional
            else sec
        )

    def millisecond(self) -> Expr:
        """
        Extract milliseconds from underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ...     pl.col("datetime").dt.millisecond().alias("millisecond"),
        ... )
        shape: (3, 5)
        ┌─────────────────────────┬──────┬────────┬────────┬─────────────┐
        │ datetime                ┆ hour ┆ minute ┆ second ┆ millisecond │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    ┆ ---         │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     ┆ i32         │
        ╞═════════════════════════╪══════╪════════╪════════╪═════════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      ┆ 0           │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     ┆ 500         │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     ┆ 60          │
        └─────────────────────────┴──────┴────────┴────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_millisecond())

    def microsecond(self) -> Expr:
        """
        Extract microseconds from underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ...     pl.col("datetime").dt.microsecond().alias("microsecond"),
        ... )
        shape: (3, 5)
        ┌─────────────────────────┬──────┬────────┬────────┬─────────────┐
        │ datetime                ┆ hour ┆ minute ┆ second ┆ microsecond │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    ┆ ---         │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     ┆ i32         │
        ╞═════════════════════════╪══════╪════════╪════════╪═════════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      ┆ 0           │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     ┆ 500000      │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     ┆ 60000       │
        └─────────────────────────┴──────┴────────┴────────┴─────────────┘
        """
        return wrap_expr(self._pyexpr.dt_microsecond())

    def nanosecond(self) -> Expr:
        """
        Extract nanoseconds from underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Expr
            Expression of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "datetime": [
        ...             datetime(1978, 1, 1, 1, 1, 1, 0),
        ...             datetime(2024, 10, 13, 5, 30, 14, 500_000),
        ...             datetime(2065, 1, 1, 10, 20, 30, 60_000),
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("datetime").dt.hour().alias("hour"),
        ...     pl.col("datetime").dt.minute().alias("minute"),
        ...     pl.col("datetime").dt.second().alias("second"),
        ...     pl.col("datetime").dt.nanosecond().alias("nanosecond"),
        ... )
        shape: (3, 5)
        ┌─────────────────────────┬──────┬────────┬────────┬────────────┐
        │ datetime                ┆ hour ┆ minute ┆ second ┆ nanosecond │
        │ ---                     ┆ ---  ┆ ---    ┆ ---    ┆ ---        │
        │ datetime[μs]            ┆ i8   ┆ i8     ┆ i8     ┆ i32        │
        ╞═════════════════════════╪══════╪════════╪════════╪════════════╡
        │ 1978-01-01 01:01:01     ┆ 1    ┆ 1      ┆ 1      ┆ 0          │
        │ 2024-10-13 05:30:14.500 ┆ 5    ┆ 30     ┆ 14     ┆ 500000000  │
        │ 2065-01-01 10:20:30.060 ┆ 10   ┆ 20     ┆ 30     ┆ 60000000   │
        └─────────────────────────┴──────┴────────┴────────┴────────────┘
        """
        return wrap_expr(self._pyexpr.dt_nanosecond())

    def epoch(self, time_unit: EpochTimeUnit = "us") -> Expr:
        """
        Get the time passed since the Unix EPOCH in the give time unit.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms', 's', 'd'}
            Time unit.

        Examples
        --------
        >>> from datetime import date
        >>> df = (
        ...     pl.date_range(date(2001, 1, 1), date(2001, 1, 3), eager=True)
        ...     .alias("date")
        ...     .to_frame()
        ... )
        >>> df.with_columns(
        ...     pl.col("date").dt.epoch().alias("epoch_ns"),
        ...     pl.col("date").dt.epoch(time_unit="s").alias("epoch_s"),
        ... )
        shape: (3, 3)
        ┌────────────┬─────────────────┬───────────┐
        │ date       ┆ epoch_ns        ┆ epoch_s   │
        │ ---        ┆ ---             ┆ ---       │
        │ date       ┆ i64             ┆ i64       │
        ╞════════════╪═════════════════╪═══════════╡
        │ 2001-01-01 ┆ 978307200000000 ┆ 978307200 │
        │ 2001-01-02 ┆ 978393600000000 ┆ 978393600 │
        │ 2001-01-03 ┆ 978480000000000 ┆ 978480000 │
        └────────────┴─────────────────┴───────────┘
        """
        if time_unit in DTYPE_TEMPORAL_UNITS:
            return self.timestamp(time_unit)  # type: ignore[arg-type]
        elif time_unit == "s":
            return self.timestamp("ms") // F.lit(1000, Int64)
        elif time_unit == "d":
            return wrap_expr(self._pyexpr).cast(Date).cast(Int32)
        else:
            msg = f"`time_unit` must be one of {{'ns', 'us', 'ms', 's', 'd'}}, got {time_unit!r}"
            raise ValueError(msg)

    def timestamp(self, time_unit: TimeUnit = "us") -> Expr:
        """
        Return a timestamp in the given time unit.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms'}
            Time unit.

        Examples
        --------
        >>> from datetime import date
        >>> df = (
        ...     pl.date_range(date(2001, 1, 1), date(2001, 1, 3), eager=True)
        ...     .alias("date")
        ...     .to_frame()
        ... )
        >>> df.with_columns(
        ...     pl.col("date").dt.timestamp().alias("timestamp_us"),
        ...     pl.col("date").dt.timestamp("ms").alias("timestamp_ms"),
        ... )
        shape: (3, 3)
        ┌────────────┬─────────────────┬──────────────┐
        │ date       ┆ timestamp_us    ┆ timestamp_ms │
        │ ---        ┆ ---             ┆ ---          │
        │ date       ┆ i64             ┆ i64          │
        ╞════════════╪═════════════════╪══════════════╡
        │ 2001-01-01 ┆ 978307200000000 ┆ 978307200000 │
        │ 2001-01-02 ┆ 978393600000000 ┆ 978393600000 │
        │ 2001-01-03 ┆ 978480000000000 ┆ 978480000000 │
        └────────────┴─────────────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.dt_timestamp(time_unit))

    @deprecated(
        "`dt.with_time_unit` is deprecated; instead, first cast "
        "to `Int64` and then cast to the desired data type."
    )
    def with_time_unit(self, time_unit: TimeUnit) -> Expr:
        """
        Set time unit of an expression of dtype Datetime or Duration.

        .. deprecated:: 0.20.5
            First cast to `Int64` and then cast to the desired data type.

        This does not modify underlying data, and should be used to fix an incorrect
        time unit.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms'}
            Unit of time for the `Datetime` or `Duration` expression.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2001, 1, 1),
        ...             datetime(2001, 1, 3),
        ...             "1d",
        ...             time_unit="ns",
        ...             eager=True,
        ...         )
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("date"),
        ...     pl.col("date").dt.with_time_unit("us").alias("time_unit_us"),
        ... )  # doctest: +SKIP
        shape: (3, 2)
        ┌─────────────────────┬───────────────────────┐
        │ date                ┆ time_unit_us          │
        │ ---                 ┆ ---                   │
        │ datetime[ns]        ┆ datetime[μs]          │
        ╞═════════════════════╪═══════════════════════╡
        │ 2001-01-01 00:00:00 ┆ +32971-04-28 00:00:00 │
        │ 2001-01-02 00:00:00 ┆ +32974-01-22 00:00:00 │
        │ 2001-01-03 00:00:00 ┆ +32976-10-18 00:00:00 │
        └─────────────────────┴───────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_with_time_unit(time_unit))

    def cast_time_unit(self, time_unit: TimeUnit) -> Expr:
        """
        Cast the underlying data to another time unit. This may lose precision.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms'}
            Time unit for the `Datetime` expression.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2001, 1, 1), datetime(2001, 1, 3), "1d", eager=True
        ...         )
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("date"),
        ...         pl.col("date").dt.cast_time_unit("ms").alias("time_unit_ms"),
        ...         pl.col("date").dt.cast_time_unit("ns").alias("time_unit_ns"),
        ...     ]
        ... )
        shape: (3, 3)
        ┌─────────────────────┬─────────────────────┬─────────────────────┐
        │ date                ┆ time_unit_ms        ┆ time_unit_ns        │
        │ ---                 ┆ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ datetime[ms]        ┆ datetime[ns]        │
        ╞═════════════════════╪═════════════════════╪═════════════════════╡
        │ 2001-01-01 00:00:00 ┆ 2001-01-01 00:00:00 ┆ 2001-01-01 00:00:00 │
        │ 2001-01-02 00:00:00 ┆ 2001-01-02 00:00:00 ┆ 2001-01-02 00:00:00 │
        │ 2001-01-03 00:00:00 ┆ 2001-01-03 00:00:00 ┆ 2001-01-03 00:00:00 │
        └─────────────────────┴─────────────────────┴─────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_cast_time_unit(time_unit))

    def convert_time_zone(self, time_zone: str) -> Expr:
        """
        Convert to given time zone for an expression of type Datetime.

        Parameters
        ----------
        time_zone
            Time zone for the `Datetime` expression.

        Notes
        -----
        If converting from a time-zone-naive datetime, then conversion will happen
        as if converting from UTC, regardless of your system's time zone.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 3, 1),
        ...             datetime(2020, 5, 1),
        ...             "1mo",
        ...             time_zone="UTC",
        ...             eager=True,
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("date"),
        ...         pl.col("date")
        ...         .dt.convert_time_zone(time_zone="Europe/London")
        ...         .alias("London"),
        ...     ]
        ... )
        shape: (3, 2)
        ┌─────────────────────────┬─────────────────────────────┐
        │ date                    ┆ London                      │
        │ ---                     ┆ ---                         │
        │ datetime[μs, UTC]       ┆ datetime[μs, Europe/London] │
        ╞═════════════════════════╪═════════════════════════════╡
        │ 2020-03-01 00:00:00 UTC ┆ 2020-03-01 00:00:00 GMT     │
        │ 2020-04-01 00:00:00 UTC ┆ 2020-04-01 01:00:00 BST     │
        │ 2020-05-01 00:00:00 UTC ┆ 2020-05-01 01:00:00 BST     │
        └─────────────────────────┴─────────────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_convert_time_zone(time_zone))

    def replace_time_zone(
        self,
        time_zone: str | None,
        *,
        ambiguous: Ambiguous | Expr = "raise",
        non_existent: NonExistent = "raise",
    ) -> Expr:
        """
        Replace time zone for an expression of type Datetime.

        Different from `convert_time_zone`, this will also modify
        the underlying timestamp and will ignore the original time zone.

        Parameters
        ----------
        time_zone
            Time zone for the `Datetime` expression. Pass `None` to unset time zone.
        ambiguous
            Determine how to deal with ambiguous datetimes:

            - `'raise'` (default): raise
            - `'earliest'`: use the earliest datetime
            - `'latest'`: use the latest datetime
            - `'null'`: set to null
        non_existent
            Determine how to deal with non-existent datetimes:

            - `'raise'` (default): raise
            - `'null'`: set to null

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "london_timezone": pl.datetime_range(
        ...             datetime(2020, 3, 1),
        ...             datetime(2020, 7, 1),
        ...             "1mo",
        ...             time_zone="UTC",
        ...             eager=True,
        ...         ).dt.convert_time_zone(time_zone="Europe/London"),
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("london_timezone"),
        ...         pl.col("london_timezone")
        ...         .dt.replace_time_zone(time_zone="Europe/Amsterdam")
        ...         .alias("London_to_Amsterdam"),
        ...     ]
        ... )
        shape: (5, 2)
        ┌─────────────────────────────┬────────────────────────────────┐
        │ london_timezone             ┆ London_to_Amsterdam            │
        │ ---                         ┆ ---                            │
        │ datetime[μs, Europe/London] ┆ datetime[μs, Europe/Amsterdam] │
        ╞═════════════════════════════╪════════════════════════════════╡
        │ 2020-03-01 00:00:00 GMT     ┆ 2020-03-01 00:00:00 CET        │
        │ 2020-04-01 01:00:00 BST     ┆ 2020-04-01 01:00:00 CEST       │
        │ 2020-05-01 01:00:00 BST     ┆ 2020-05-01 01:00:00 CEST       │
        │ 2020-06-01 01:00:00 BST     ┆ 2020-06-01 01:00:00 CEST       │
        │ 2020-07-01 01:00:00 BST     ┆ 2020-07-01 01:00:00 CEST       │
        └─────────────────────────────┴────────────────────────────────┘

        You can use `ambiguous` to deal with ambiguous datetimes:

        >>> dates = [
        ...     "2018-10-28 01:30",
        ...     "2018-10-28 02:00",
        ...     "2018-10-28 02:30",
        ...     "2018-10-28 02:00",
        ... ]
        >>> df = pl.DataFrame(
        ...     {
        ...         "ts": pl.Series(dates).str.strptime(pl.Datetime),
        ...         "ambiguous": ["earliest", "earliest", "latest", "latest"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     ts_localized=pl.col("ts").dt.replace_time_zone(
        ...         "Europe/Brussels", ambiguous=pl.col("ambiguous")
        ...     )
        ... )
        shape: (4, 3)
        ┌─────────────────────┬───────────┬───────────────────────────────┐
        │ ts                  ┆ ambiguous ┆ ts_localized                  │
        │ ---                 ┆ ---       ┆ ---                           │
        │ datetime[μs]        ┆ str       ┆ datetime[μs, Europe/Brussels] │
        ╞═════════════════════╪═══════════╪═══════════════════════════════╡
        │ 2018-10-28 01:30:00 ┆ earliest  ┆ 2018-10-28 01:30:00 CEST      │
        │ 2018-10-28 02:00:00 ┆ earliest  ┆ 2018-10-28 02:00:00 CEST      │
        │ 2018-10-28 02:30:00 ┆ latest    ┆ 2018-10-28 02:30:00 CET       │
        │ 2018-10-28 02:00:00 ┆ latest    ┆ 2018-10-28 02:00:00 CET       │
        └─────────────────────┴───────────┴───────────────────────────────┘
        """
        if not isinstance(ambiguous, pl.Expr):
            ambiguous = F.lit(ambiguous)
        return wrap_expr(
            self._pyexpr.dt_replace_time_zone(
                time_zone, ambiguous._pyexpr, non_existent
            )
        )

    def total_days(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total days from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the day.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 3, 1), datetime(2020, 5, 1), "1mo", eager=True
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("date"),
        ...         pl.col("date").diff().dt.total_days().alias("days_diff"),
        ...     ]
        ... )
        shape: (3, 2)
        ┌─────────────────────┬───────────┐
        │ date                ┆ days_diff │
        │ ---                 ┆ ---       │
        │ datetime[μs]        ┆ i64       │
        ╞═════════════════════╪═══════════╡
        │ 2020-03-01 00:00:00 ┆ null      │
        │ 2020-04-01 00:00:00 ┆ 31        │
        │ 2020-05-01 00:00:00 ┆ 30        │
        └─────────────────────┴───────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_days(fractional))

    def total_hours(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total hours from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the hour.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1), datetime(2020, 1, 4), "1d", eager=True
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("date"),
        ...         pl.col("date").diff().dt.total_hours().alias("hours_diff"),
        ...     ]
        ... )
        shape: (4, 2)
        ┌─────────────────────┬────────────┐
        │ date                ┆ hours_diff │
        │ ---                 ┆ ---        │
        │ datetime[μs]        ┆ i64        │
        ╞═════════════════════╪════════════╡
        │ 2020-01-01 00:00:00 ┆ null       │
        │ 2020-01-02 00:00:00 ┆ 24         │
        │ 2020-01-03 00:00:00 ┆ 24         │
        │ 2020-01-04 00:00:00 ┆ 24         │
        └─────────────────────┴────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_hours(fractional))

    def total_minutes(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total minutes from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the minute.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1), datetime(2020, 1, 4), "1d", eager=True
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("date"),
        ...         pl.col("date").diff().dt.total_minutes().alias("minutes_diff"),
        ...     ]
        ... )
        shape: (4, 2)
        ┌─────────────────────┬──────────────┐
        │ date                ┆ minutes_diff │
        │ ---                 ┆ ---          │
        │ datetime[μs]        ┆ i64          │
        ╞═════════════════════╪══════════════╡
        │ 2020-01-01 00:00:00 ┆ null         │
        │ 2020-01-02 00:00:00 ┆ 1440         │
        │ 2020-01-03 00:00:00 ┆ 1440         │
        │ 2020-01-04 00:00:00 ┆ 1440         │
        └─────────────────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_minutes(fractional))

    def total_seconds(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total seconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the second.

        Returns
        -------
        Expr
            Expression of data type :py:class:`.Int64` or :py:class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1),
        ...             datetime(2020, 1, 1, 0, 4, 0),
        ...             "1m",
        ...             eager=True,
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("date"),
        ...     pl.col("date").diff().dt.total_seconds().alias("seconds_diff"),
        ... )
        shape: (5, 2)
        ┌─────────────────────┬──────────────┐
        │ date                ┆ seconds_diff │
        │ ---                 ┆ ---          │
        │ datetime[μs]        ┆ i64          │
        ╞═════════════════════╪══════════════╡
        │ 2020-01-01 00:00:00 ┆ null         │
        │ 2020-01-01 00:01:00 ┆ 60           │
        │ 2020-01-01 00:02:00 ┆ 60           │
        │ 2020-01-01 00:03:00 ┆ 60           │
        │ 2020-01-01 00:04:00 ┆ 60           │
        └─────────────────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_seconds(fractional))

    def total_milliseconds(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total milliseconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the millisecond.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1),
        ...             datetime(2020, 1, 1, 0, 0, 1, 0),
        ...             "200ms",
        ...             eager=True,
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("date"),
        ...     milliseconds_diff=pl.col("date").diff().dt.total_milliseconds(),
        ... )
        shape: (6, 2)
        ┌─────────────────────────┬───────────────────┐
        │ date                    ┆ milliseconds_diff │
        │ ---                     ┆ ---               │
        │ datetime[μs]            ┆ i64               │
        ╞═════════════════════════╪═══════════════════╡
        │ 2020-01-01 00:00:00     ┆ null              │
        │ 2020-01-01 00:00:00.200 ┆ 200               │
        │ 2020-01-01 00:00:00.400 ┆ 200               │
        │ 2020-01-01 00:00:00.600 ┆ 200               │
        │ 2020-01-01 00:00:00.800 ┆ 200               │
        │ 2020-01-01 00:00:01     ┆ 200               │
        └─────────────────────────┴───────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_milliseconds(fractional))

    def total_microseconds(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total microseconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the microsecond.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1),
        ...             datetime(2020, 1, 1, 0, 0, 1, 0),
        ...             "200ms",
        ...             eager=True,
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("date"),
        ...     milliseconds_diff=pl.col("date").diff().dt.total_microseconds(),
        ... )
        shape: (6, 2)
        ┌─────────────────────────┬───────────────────┐
        │ date                    ┆ milliseconds_diff │
        │ ---                     ┆ ---               │
        │ datetime[μs]            ┆ i64               │
        ╞═════════════════════════╪═══════════════════╡
        │ 2020-01-01 00:00:00     ┆ null              │
        │ 2020-01-01 00:00:00.200 ┆ 200000            │
        │ 2020-01-01 00:00:00.400 ┆ 200000            │
        │ 2020-01-01 00:00:00.600 ┆ 200000            │
        │ 2020-01-01 00:00:00.800 ┆ 200000            │
        │ 2020-01-01 00:00:01     ┆ 200000            │
        └─────────────────────────┴───────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_microseconds(fractional))

    def total_nanoseconds(self, *, fractional: bool = False) -> Expr:
        """
        Extract the total nanoseconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include return the result as a :class:`.Float64`.
            Because the smallest :type:`.TimeUnit` is `'ns'`, the
            fractional component will always be zero.

        Returns
        -------
        Expr
            Expression of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": pl.datetime_range(
        ...             datetime(2020, 1, 1),
        ...             datetime(2020, 1, 1, 0, 0, 1, 0),
        ...             "200ms",
        ...             eager=True,
        ...         ),
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("date"),
        ...     milliseconds_diff=pl.col("date").diff().dt.total_nanoseconds(),
        ... )
        shape: (6, 2)
        ┌─────────────────────────┬───────────────────┐
        │ date                    ┆ milliseconds_diff │
        │ ---                     ┆ ---               │
        │ datetime[μs]            ┆ i64               │
        ╞═════════════════════════╪═══════════════════╡
        │ 2020-01-01 00:00:00     ┆ null              │
        │ 2020-01-01 00:00:00.200 ┆ 200000000         │
        │ 2020-01-01 00:00:00.400 ┆ 200000000         │
        │ 2020-01-01 00:00:00.600 ┆ 200000000         │
        │ 2020-01-01 00:00:00.800 ┆ 200000000         │
        │ 2020-01-01 00:00:01     ┆ 200000000         │
        └─────────────────────────┴───────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_total_nanoseconds(fractional))

    def offset_by(self, by: str | Expr) -> Expr:
        """
        Offset this date by a relative time offset.

        This differs from `pl.col("foo") + timedelta` in that it can
        take months and leap years into account. Note that only a single minus
        sign is allowed in the `by` string, as the first character.

        Parameters
        ----------
        by
            The offset is dictated by the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)

            By "calendar day", we mean the corresponding time on the next day (which may
            not be 24 hours, due to daylight savings). Similarly for "calendar week",
            "calendar month", "calendar quarter", and "calendar year".

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "dates": pl.datetime_range(
        ...             datetime(2000, 1, 1), datetime(2005, 1, 1), "1y", eager=True
        ...         ),
        ...         "offset": ["1d", "2d", "-1d", "1mo", None, "1y"],
        ...     }
        ... )
        >>> df.select(
        ...     [
        ...         pl.col("dates").dt.offset_by("1y").alias("date_plus_1y"),
        ...         pl.col("dates").dt.offset_by("-1y2mo").alias("date_min"),
        ...     ]
        ... )
        shape: (6, 2)
        ┌─────────────────────┬─────────────────────┐
        │ date_plus_1y        ┆ date_min            │
        │ ---                 ┆ ---                 │
        │ datetime[μs]        ┆ datetime[μs]        │
        ╞═════════════════════╪═════════════════════╡
        │ 2001-01-01 00:00:00 ┆ 1998-11-01 00:00:00 │
        │ 2002-01-01 00:00:00 ┆ 1999-11-01 00:00:00 │
        │ 2003-01-01 00:00:00 ┆ 2000-11-01 00:00:00 │
        │ 2004-01-01 00:00:00 ┆ 2001-11-01 00:00:00 │
        │ 2005-01-01 00:00:00 ┆ 2002-11-01 00:00:00 │
        │ 2006-01-01 00:00:00 ┆ 2003-11-01 00:00:00 │
        └─────────────────────┴─────────────────────┘

        You can also pass the relative offset as an expression:

        >>> df.with_columns(new_dates=pl.col("dates").dt.offset_by(pl.col("offset")))
        shape: (6, 3)
        ┌─────────────────────┬────────┬─────────────────────┐
        │ dates               ┆ offset ┆ new_dates           │
        │ ---                 ┆ ---    ┆ ---                 │
        │ datetime[μs]        ┆ str    ┆ datetime[μs]        │
        ╞═════════════════════╪════════╪═════════════════════╡
        │ 2000-01-01 00:00:00 ┆ 1d     ┆ 2000-01-02 00:00:00 │
        │ 2001-01-01 00:00:00 ┆ 2d     ┆ 2001-01-03 00:00:00 │
        │ 2002-01-01 00:00:00 ┆ -1d    ┆ 2001-12-31 00:00:00 │
        │ 2003-01-01 00:00:00 ┆ 1mo    ┆ 2003-02-01 00:00:00 │
        │ 2004-01-01 00:00:00 ┆ null   ┆ null                │
        │ 2005-01-01 00:00:00 ┆ 1y     ┆ 2006-01-01 00:00:00 │
        └─────────────────────┴────────┴─────────────────────┘
        """
        by_pyexpr = parse_into_expression(by, str_as_lit=True)
        return wrap_expr(self._pyexpr.dt_offset_by(by_pyexpr))

    def month_start(self) -> Expr:
        """
        Roll backward to the first day of the month.

        For datetimes, the time-of-day is preserved.

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        If you're coming from pandas, you can think of this as a vectorised version
        of `pandas.tseries.offsets.MonthBegin().rollback(datetime)`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "dates": pl.datetime_range(
        ...             datetime(2000, 1, 15, 2),
        ...             datetime(2000, 12, 15, 2),
        ...             "1mo",
        ...             eager=True,
        ...         )
        ...     }
        ... )
        >>> df.select(pl.col("dates").dt.month_start())
        shape: (12, 1)
        ┌─────────────────────┐
        │ dates               │
        │ ---                 │
        │ datetime[μs]        │
        ╞═════════════════════╡
        │ 2000-01-01 02:00:00 │
        │ 2000-02-01 02:00:00 │
        │ 2000-03-01 02:00:00 │
        │ 2000-04-01 02:00:00 │
        │ 2000-05-01 02:00:00 │
        │ …                   │
        │ 2000-08-01 02:00:00 │
        │ 2000-09-01 02:00:00 │
        │ 2000-10-01 02:00:00 │
        │ 2000-11-01 02:00:00 │
        │ 2000-12-01 02:00:00 │
        └─────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_month_start())

    def month_end(self) -> Expr:
        """
        Roll forward to the last day of the month.

        For datetimes, the time-of-day is preserved.

        Returns
        -------
        Expr
            Expression of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        If you're coming from pandas, you can think of this as a vectorised version
        of `pandas.tseries.offsets.MonthEnd().rollforward(datetime)`.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "dates": pl.datetime_range(
        ...             datetime(2000, 1, 1, 2),
        ...             datetime(2000, 12, 1, 2),
        ...             "1mo",
        ...             eager=True,
        ...         )
        ...     }
        ... )
        >>> df.select(pl.col("dates").dt.month_end())
        shape: (12, 1)
        ┌─────────────────────┐
        │ dates               │
        │ ---                 │
        │ datetime[μs]        │
        ╞═════════════════════╡
        │ 2000-01-31 02:00:00 │
        │ 2000-02-29 02:00:00 │
        │ 2000-03-31 02:00:00 │
        │ 2000-04-30 02:00:00 │
        │ 2000-05-31 02:00:00 │
        │ …                   │
        │ 2000-08-31 02:00:00 │
        │ 2000-09-30 02:00:00 │
        │ 2000-10-31 02:00:00 │
        │ 2000-11-30 02:00:00 │
        │ 2000-12-31 02:00:00 │
        └─────────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_month_end())

    def base_utc_offset(self) -> Expr:
        """
        Base offset from UTC.

        This is usually constant for all datetimes in a given time zone, but
        may vary in the rare case that a country switches time zone, like
        Samoa (Apia) did at the end of 2011.

        Returns
        -------
        Expr
            Expression of data type :class:`Duration`.

        See Also
        --------
        Expr.dt.dst_offset : Daylight savings offset from UTC.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "ts": [datetime(2011, 12, 29), datetime(2012, 1, 1)],
        ...     }
        ... )
        >>> df = df.with_columns(pl.col("ts").dt.replace_time_zone("Pacific/Apia"))
        >>> df.with_columns(pl.col("ts").dt.base_utc_offset().alias("base_utc_offset"))
        shape: (2, 2)
        ┌────────────────────────────┬─────────────────┐
        │ ts                         ┆ base_utc_offset │
        │ ---                        ┆ ---             │
        │ datetime[μs, Pacific/Apia] ┆ duration[ms]    │
        ╞════════════════════════════╪═════════════════╡
        │ 2011-12-29 00:00:00 -10    ┆ -11h            │
        │ 2012-01-01 00:00:00 +14    ┆ 13h             │
        └────────────────────────────┴─────────────────┘
        """
        return wrap_expr(self._pyexpr.dt_base_utc_offset())

    def dst_offset(self) -> Expr:
        """
        Additional offset currently in effect (typically due to daylight saving time).

        Returns
        -------
        Expr
            Expression of data type :class:`Duration`.

        See Also
        --------
        Expr.dt.base_utc_offset : Base offset from UTC.

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "ts": [datetime(2020, 10, 25), datetime(2020, 10, 26)],
        ...     }
        ... )
        >>> df = df.with_columns(pl.col("ts").dt.replace_time_zone("Europe/London"))
        >>> df.with_columns(pl.col("ts").dt.dst_offset().alias("dst_offset"))
        shape: (2, 2)
        ┌─────────────────────────────┬──────────────┐
        │ ts                          ┆ dst_offset   │
        │ ---                         ┆ ---          │
        │ datetime[μs, Europe/London] ┆ duration[ms] │
        ╞═════════════════════════════╪══════════════╡
        │ 2020-10-25 00:00:00 BST     ┆ 1h           │
        │ 2020-10-26 00:00:00 GMT     ┆ 0ms          │
        └─────────────────────────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.dt_dst_offset())
