from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, overload

from polars import functions as F
from polars._utils.parse import (
    parse_into_expression,
    parse_into_list_of_expressions,
)
from polars._utils.unstable import issue_unstable_warning
from polars._utils.wrap import wrap_expr
from polars.datatypes import Date, Struct, Time

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr


if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from polars import Expr, Series
    from polars._typing import Ambiguous, IntoExpr, SchemaDict, TimeUnit


def datetime_(
    year: int | IntoExpr,
    month: int | IntoExpr,
    day: int | IntoExpr,
    hour: int | IntoExpr | None = None,
    minute: int | IntoExpr | None = None,
    second: int | IntoExpr | None = None,
    microsecond: int | IntoExpr | None = None,
    *,
    time_unit: TimeUnit = "us",
    time_zone: str | None = None,
    ambiguous: Ambiguous | Expr = "raise",
) -> Expr:
    """
    Create a Polars literal expression of type Datetime.

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
    time_unit : {'us', 'ms', 'ns'}
        Time unit of the resulting expression.
    time_zone
        Time zone of the resulting expression.
    ambiguous
        Determine how to deal with ambiguous datetimes:

        - `'raise'` (default): raise
        - `'earliest'`: use the earliest datetime
        - `'latest'`: use the latest datetime
        - `'null'`: set to null

    Returns
    -------
    Expr
        Expression of data type :class:`Datetime`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "month": [1, 2, 3],
    ...         "day": [4, 5, 6],
    ...         "hour": [12, 13, 14],
    ...         "minute": [15, 30, 45],
    ...     }
    ... )
    >>> df.with_columns(
    ...     pl.datetime(
    ...         2024,
    ...         pl.col("month"),
    ...         pl.col("day"),
    ...         pl.col("hour"),
    ...         pl.col("minute"),
    ...         time_zone="Australia/Sydney",
    ...     )
    ... )
    shape: (3, 5)
    ┌───────┬─────┬──────┬────────┬────────────────────────────────┐
    │ month ┆ day ┆ hour ┆ minute ┆ datetime                       │
    │ ---   ┆ --- ┆ ---  ┆ ---    ┆ ---                            │
    │ i64   ┆ i64 ┆ i64  ┆ i64    ┆ datetime[μs, Australia/Sydney] │
    ╞═══════╪═════╪══════╪════════╪════════════════════════════════╡
    │ 1     ┆ 4   ┆ 12   ┆ 15     ┆ 2024-01-04 12:15:00 AEDT       │
    │ 2     ┆ 5   ┆ 13   ┆ 30     ┆ 2024-02-05 13:30:00 AEDT       │
    │ 3     ┆ 6   ┆ 14   ┆ 45     ┆ 2024-03-06 14:45:00 AEDT       │
    └───────┴─────┴──────┴────────┴────────────────────────────────┘

    We can also use `pl.datetime` for filtering:

    >>> from datetime import datetime
    >>> df = pl.DataFrame(
    ...     {
    ...         "start": [
    ...             datetime(2024, 1, 1, 0, 0, 0),
    ...             datetime(2024, 1, 1, 0, 0, 0),
    ...             datetime(2024, 1, 1, 0, 0, 0),
    ...         ],
    ...         "end": [
    ...             datetime(2024, 5, 1, 20, 15, 10),
    ...             datetime(2024, 7, 1, 21, 25, 20),
    ...             datetime(2024, 9, 1, 22, 35, 30),
    ...         ],
    ...     }
    ... )
    >>> df.filter(pl.col("end") > pl.datetime(2024, 6, 1))
        shape: (2, 2)
    ┌─────────────────────┬─────────────────────┐
    │ start               ┆ end                 │
    │ ---                 ┆ ---                 │
    │ datetime[μs]        ┆ datetime[μs]        │
    ╞═════════════════════╪═════════════════════╡
    │ 2024-01-01 00:00:00 ┆ 2024-07-01 21:25:20 │
    │ 2024-01-01 00:00:00 ┆ 2024-09-01 22:35:30 │
    └─────────────────────┴─────────────────────┘
    """
    ambiguous_expr = parse_into_expression(ambiguous, str_as_lit=True)
    year_expr = parse_into_expression(year)
    month_expr = parse_into_expression(month)
    day_expr = parse_into_expression(day)

    hour_expr = parse_into_expression(hour) if hour is not None else None
    minute_expr = parse_into_expression(minute) if minute is not None else None
    second_expr = parse_into_expression(second) if second is not None else None
    microsecond_expr = (
        parse_into_expression(microsecond) if microsecond is not None else None
    )

    return wrap_expr(
        plr.datetime(
            year_expr,
            month_expr,
            day_expr,
            hour_expr,
            minute_expr,
            second_expr,
            microsecond_expr,
            time_unit,
            time_zone,
            ambiguous_expr,
        )
    )


def date_(
    year: Expr | str | int,
    month: Expr | str | int,
    day: Expr | str | int,
) -> Expr:
    """
    Create a Polars literal expression of type Date.

    Parameters
    ----------
    year
        column or literal.
    month
        column or literal, ranging from 1-12.
    day
        column or literal, ranging from 1-31.

    Returns
    -------
    Expr
        Expression of data type :class:`Date`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "month": [1, 2, 3],
    ...         "day": [4, 5, 6],
    ...     }
    ... )
    >>> df.with_columns(pl.date(2024, pl.col("month"), pl.col("day")))
    shape: (3, 3)
    ┌───────┬─────┬────────────┐
    │ month ┆ day ┆ date       │
    │ ---   ┆ --- ┆ ---        │
    │ i64   ┆ i64 ┆ date       │
    ╞═══════╪═════╪════════════╡
    │ 1     ┆ 4   ┆ 2024-01-04 │
    │ 2     ┆ 5   ┆ 2024-02-05 │
    │ 3     ┆ 6   ┆ 2024-03-06 │
    └───────┴─────┴────────────┘

    We can also use `pl.date` for filtering:

    >>> from datetime import date
    >>> df = pl.DataFrame(
    ...     {
    ...         "start": [date(2024, 1, 1), date(2024, 1, 1), date(2024, 1, 1)],
    ...         "end": [date(2024, 5, 1), date(2024, 7, 1), date(2024, 9, 1)],
    ...     }
    ... )
    >>> df.filter(pl.col("end") > pl.date(2024, 6, 1))
    shape: (2, 2)
    ┌────────────┬────────────┐
    │ start      ┆ end        │
    │ ---        ┆ ---        │
    │ date       ┆ date       │
    ╞════════════╪════════════╡
    │ 2024-01-01 ┆ 2024-07-01 │
    │ 2024-01-01 ┆ 2024-09-01 │
    └────────────┴────────────┘
    """
    return datetime_(year, month, day).cast(Date).alias("date")


def time_(
    hour: Expr | str | int | None = None,
    minute: Expr | str | int | None = None,
    second: Expr | str | int | None = None,
    microsecond: Expr | str | int | None = None,
) -> Expr:
    """
    Create a Polars literal expression of type Time.

    Parameters
    ----------
    hour
        column or literal, ranging from 0-23.
    minute
        column or literal, ranging from 0-59.
    second
        column or literal, ranging from 0-59.
    microsecond
        column or literal, ranging from 0-999999.

    Returns
    -------
    Expr
        Expression of data type :class:`Date`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "hour": [12, 13, 14],
    ...         "minute": [15, 30, 45],
    ...     }
    ... )

    >>> df.with_columns(pl.time(pl.col("hour"), pl.col("minute")))
    shape: (3, 3)
    ┌──────┬────────┬──────────┐
    │ hour ┆ minute ┆ time     │
    │ ---  ┆ ---    ┆ ---      │
    │ i64  ┆ i64    ┆ time     │
    ╞══════╪════════╪══════════╡
    │ 12   ┆ 15     ┆ 12:15:00 │
    │ 13   ┆ 30     ┆ 13:30:00 │
    │ 14   ┆ 45     ┆ 14:45:00 │
    └──────┴────────┴──────────┘
    """
    epoch_start = (1970, 1, 1)
    return (
        datetime_(*epoch_start, hour, minute, second, microsecond)
        .cast(Time)
        .alias("time")
    )


def duration(
    *,
    weeks: Expr | str | int | float | None = None,
    days: Expr | str | int | float | None = None,
    hours: Expr | str | int | float | None = None,
    minutes: Expr | str | int | float | None = None,
    seconds: Expr | str | int | float | None = None,
    milliseconds: Expr | str | int | float | None = None,
    microseconds: Expr | str | int | float | None = None,
    nanoseconds: Expr | str | int | float | None = None,
    time_unit: TimeUnit | None = None,
) -> Expr:
    """
    Create polars `Duration` from distinct time components.

    Parameters
    ----------
    weeks
        Number of weeks.
    days
        Number of days.
    hours
        Number of hours.
    minutes
        Number of minutes.
    seconds
        Number of seconds.
    milliseconds
        Number of milliseconds.
    microseconds
        Number of microseconds.
    nanoseconds
        Number of nanoseconds.
    time_unit : {None, 'us', 'ms', 'ns'}
        Time unit of the resulting expression. If set to `None` (default), the time
        unit will be inferred from the other inputs: `'ns'` if `nanoseconds` was
        specified, `'us'` otherwise.

    Returns
    -------
    Expr
        Expression of data type :class:`Duration`.

    Notes
    -----
    A `duration` represents a fixed amount of time. For example,
    `pl.duration(days=1)` means "exactly 24 hours". By contrast,
    `Expr.dt.offset_by('1d')` means "1 calendar day", which could sometimes be
    23 hours or 25 hours depending on Daylight Savings Time.
    For non-fixed durations such as "calendar month" or "calendar day",
    please use :meth:`polars.Expr.dt.offset_by` instead.

    Examples
    --------
    >>> from datetime import datetime
    >>> df = pl.DataFrame(
    ...     {
    ...         "dt": [datetime(2022, 1, 1), datetime(2022, 1, 2)],
    ...         "add": [1, 2],
    ...     }
    ... )
    >>> df
    shape: (2, 2)
    ┌─────────────────────┬─────┐
    │ dt                  ┆ add │
    │ ---                 ┆ --- │
    │ datetime[μs]        ┆ i64 │
    ╞═════════════════════╪═════╡
    │ 2022-01-01 00:00:00 ┆ 1   │
    │ 2022-01-02 00:00:00 ┆ 2   │
    └─────────────────────┴─────┘
    >>> with pl.Config(tbl_width_chars=120):
    ...     df.select(
    ...         (pl.col("dt") + pl.duration(weeks="add")).alias("add_weeks"),
    ...         (pl.col("dt") + pl.duration(days="add")).alias("add_days"),
    ...         (pl.col("dt") + pl.duration(seconds="add")).alias("add_seconds"),
    ...         (pl.col("dt") + pl.duration(milliseconds="add")).alias("add_millis"),
    ...         (pl.col("dt") + pl.duration(hours="add")).alias("add_hours"),
    ...     )
    shape: (2, 5)
    ┌─────────────────────┬─────────────────────┬─────────────────────┬─────────────────────────┬─────────────────────┐
    │ add_weeks           ┆ add_days            ┆ add_seconds         ┆ add_millis              ┆ add_hours           │
    │ ---                 ┆ ---                 ┆ ---                 ┆ ---                     ┆ ---                 │
    │ datetime[μs]        ┆ datetime[μs]        ┆ datetime[μs]        ┆ datetime[μs]            ┆ datetime[μs]        │
    ╞═════════════════════╪═════════════════════╪═════════════════════╪═════════════════════════╪═════════════════════╡
    │ 2022-01-08 00:00:00 ┆ 2022-01-02 00:00:00 ┆ 2022-01-01 00:00:01 ┆ 2022-01-01 00:00:00.001 ┆ 2022-01-01 01:00:00 │
    │ 2022-01-16 00:00:00 ┆ 2022-01-04 00:00:00 ┆ 2022-01-02 00:00:02 ┆ 2022-01-02 00:00:00.002 ┆ 2022-01-02 02:00:00 │
    └─────────────────────┴─────────────────────┴─────────────────────┴─────────────────────────┴─────────────────────┘

    If you need to add non-fixed durations, you should use :meth:`polars.Expr.dt.offset_by` instead:

    >>> with pl.Config(tbl_width_chars=120):
    ...     df.select(
    ...         add_calendar_days=pl.col("dt").dt.offset_by(
    ...             pl.format("{}d", pl.col("add"))
    ...         ),
    ...         add_calendar_months=pl.col("dt").dt.offset_by(
    ...             pl.format("{}mo", pl.col("add"))
    ...         ),
    ...         add_calendar_years=pl.col("dt").dt.offset_by(
    ...             pl.format("{}y", pl.col("add"))
    ...         ),
    ...     )
    shape: (2, 3)
    ┌─────────────────────┬─────────────────────┬─────────────────────┐
    │ add_calendar_days   ┆ add_calendar_months ┆ add_calendar_years  │
    │ ---                 ┆ ---                 ┆ ---                 │
    │ datetime[μs]        ┆ datetime[μs]        ┆ datetime[μs]        │
    ╞═════════════════════╪═════════════════════╪═════════════════════╡
    │ 2022-01-02 00:00:00 ┆ 2022-02-01 00:00:00 ┆ 2023-01-01 00:00:00 │
    │ 2022-01-04 00:00:00 ┆ 2022-03-02 00:00:00 ┆ 2024-01-02 00:00:00 │
    └─────────────────────┴─────────────────────┴─────────────────────┘
    """  # noqa: W505
    if nanoseconds is not None and time_unit is None:
        time_unit = "ns"

    weeks_expr = parse_into_expression(weeks) if weeks is not None else None
    days_expr = parse_into_expression(days) if days is not None else None
    hours_expr = parse_into_expression(hours) if hours is not None else None
    minutes_expr = parse_into_expression(minutes) if minutes is not None else None
    seconds_expr = parse_into_expression(seconds) if seconds is not None else None
    milliseconds_expr = (
        parse_into_expression(milliseconds) if milliseconds is not None else None
    )
    microseconds_expr = (
        parse_into_expression(microseconds) if microseconds is not None else None
    )
    nanoseconds_expr = (
        parse_into_expression(nanoseconds) if nanoseconds is not None else None
    )

    if time_unit is None:
        time_unit = "us"

    return wrap_expr(
        plr.duration(
            weeks_expr,
            days_expr,
            hours_expr,
            minutes_expr,
            seconds_expr,
            milliseconds_expr,
            microseconds_expr,
            nanoseconds_expr,
            time_unit,
        )
    )


def concat_list(exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr) -> Expr:
    """
    Horizontally concatenate columns into a single list column.

    Operates in linear time.

    Parameters
    ----------
    exprs
        Columns to concatenate into a single list column. Accepts expression input.
        Strings are parsed as column names, other non-expression inputs are parsed as
        literals.
    *more_exprs
        Additional columns to concatenate into a single list column, specified as
        positional arguments.

    Examples
    --------
    Concatenate two existing list columns. Null values are propagated.

    >>> df = pl.DataFrame({"a": [[1, 2], [3], [4, 5]], "b": [[4], [], None]})
    >>> df.with_columns(concat_list=pl.concat_list("a", "b"))
    shape: (3, 3)
    ┌───────────┬───────────┬─────────────┐
    │ a         ┆ b         ┆ concat_list │
    │ ---       ┆ ---       ┆ ---         │
    │ list[i64] ┆ list[i64] ┆ list[i64]   │
    ╞═══════════╪═══════════╪═════════════╡
    │ [1, 2]    ┆ [4]       ┆ [1, 2, 4]   │
    │ [3]       ┆ []        ┆ [3]         │
    │ [4, 5]    ┆ null      ┆ null        │
    └───────────┴───────────┴─────────────┘

    Non-list columns are cast to a list before concatenation. The output data type
    is the supertype of the concatenated columns.

    >>> df.select("a", concat_list=pl.concat_list("a", pl.lit("x")))
    shape: (3, 2)
    ┌───────────┬─────────────────┐
    │ a         ┆ concat_list     │
    │ ---       ┆ ---             │
    │ list[i64] ┆ list[str]       │
    ╞═══════════╪═════════════════╡
    │ [1, 2]    ┆ ["1", "2", "x"] │
    │ [3]       ┆ ["3", "x"]      │
    │ [4, 5]    ┆ ["4", "5", "x"] │
    └───────────┴─────────────────┘

    Create lagged columns and collect them into a list. This mimics a rolling window.

    >>> df = pl.DataFrame({"A": [1.0, 2.0, 9.0, 2.0, 13.0]})
    >>> df = df.select([pl.col("A").shift(i).alias(f"A_lag_{i}") for i in range(3)])
    >>> df.select(
    ...     pl.concat_list([f"A_lag_{i}" for i in range(3)][::-1]).alias("A_rolling")
    ... )
    shape: (5, 1)
    ┌───────────────────┐
    │ A_rolling         │
    │ ---               │
    │ list[f64]         │
    ╞═══════════════════╡
    │ [null, null, 1.0] │
    │ [null, 1.0, 2.0]  │
    │ [1.0, 2.0, 9.0]   │
    │ [2.0, 9.0, 2.0]   │
    │ [9.0, 2.0, 13.0]  │
    └───────────────────┘
    """
    exprs = parse_into_list_of_expressions(exprs, *more_exprs)
    return wrap_expr(plr.concat_list(exprs))


def concat_arr(exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr) -> Expr:
    """
    Horizontally concatenate columns into a single array column.

    Non-array columns are reshaped to a unit-width array. All columns must have
    a dtype of either `pl.Array(<DataType>, width)` or `pl.<DataType>`.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

    Parameters
    ----------
    exprs
        Columns to concatenate into a single array column. Accepts expression input.
        Strings are parsed as column names, other non-expression inputs are parsed as
        literals.
    *more_exprs
        Additional columns to concatenate into a single array column, specified as
        positional arguments.

    Examples
    --------
    Concatenate 2 array columns:

    >>> (
    ...     pl.select(
    ...         a=pl.Series([[1], [3], None], dtype=pl.Array(pl.Int64, 1)),
    ...         b=pl.Series([[3], [None], [5]], dtype=pl.Array(pl.Int64, 1)),
    ...     ).with_columns(
    ...         pl.concat_arr("a", "b").alias("concat_arr(a, b)"),
    ...         pl.concat_arr("a", pl.first("b")).alias("concat_arr(a, first(b))"),
    ...     )
    ... )
    shape: (3, 4)
    ┌───────────────┬───────────────┬──────────────────┬─────────────────────────┐
    │ a             ┆ b             ┆ concat_arr(a, b) ┆ concat_arr(a, first(b)) │
    │ ---           ┆ ---           ┆ ---              ┆ ---                     │
    │ array[i64, 1] ┆ array[i64, 1] ┆ array[i64, 2]    ┆ array[i64, 2]           │
    ╞═══════════════╪═══════════════╪══════════════════╪═════════════════════════╡
    │ [1]           ┆ [3]           ┆ [1, 3]           ┆ [1, 3]                  │
    │ [3]           ┆ [null]        ┆ [3, null]        ┆ [3, 3]                  │
    │ null          ┆ [5]           ┆ null             ┆ null                    │
    └───────────────┴───────────────┴──────────────────┴─────────────────────────┘

    Concatenate non-array columns:

    >>> (
    ...     pl.select(
    ...         c=pl.Series([None, 5, 6], dtype=pl.Int64),
    ...     )
    ...     .with_columns(d=pl.col("c").reverse())
    ...     .with_columns(
    ...         pl.concat_arr("c", "d").alias("concat_arr(c, d)"),
    ...     )
    ... )
    shape: (3, 3)
    ┌──────┬──────┬──────────────────┐
    │ c    ┆ d    ┆ concat_arr(c, d) │
    │ ---  ┆ ---  ┆ ---              │
    │ i64  ┆ i64  ┆ array[i64, 2]    │
    ╞══════╪══════╪══════════════════╡
    │ null ┆ 6    ┆ [null, 6]        │
    │ 5    ┆ 5    ┆ [5, 5]           │
    │ 6    ┆ null ┆ [6, null]        │
    └──────┴──────┴──────────────────┘

    Concatenate mixed array and non-array columns:

    >>> (
    ...     pl.select(
    ...         a=pl.Series([[1], [3], None], dtype=pl.Array(pl.Int64, 1)),
    ...         b=pl.Series([[3], [None], [5]], dtype=pl.Array(pl.Int64, 1)),
    ...         c=pl.Series([None, 5, 6], dtype=pl.Int64),
    ...     ).with_columns(
    ...         pl.concat_arr("a", "b", "c").alias("concat_arr(a, b, c)"),
    ...     )
    ... )
    shape: (3, 4)
    ┌───────────────┬───────────────┬──────┬─────────────────────┐
    │ a             ┆ b             ┆ c    ┆ concat_arr(a, b, c) │
    │ ---           ┆ ---           ┆ ---  ┆ ---                 │
    │ array[i64, 1] ┆ array[i64, 1] ┆ i64  ┆ array[i64, 3]       │
    ╞═══════════════╪═══════════════╪══════╪═════════════════════╡
    │ [1]           ┆ [3]           ┆ null ┆ [1, 3, null]        │
    │ [3]           ┆ [null]        ┆ 5    ┆ [3, null, 5]        │
    │ null          ┆ [5]           ┆ 6    ┆ null                │
    └───────────────┴───────────────┴──────┴─────────────────────┘

    Unit-length columns are broadcasted:

    >>> (
    ...     pl.select(
    ...         a=pl.Series([1, 3, None]),
    ...     ).with_columns(
    ...         pl.concat_arr("a", pl.lit(0, dtype=pl.Int64)).alias("concat_arr(a, 0)"),
    ...         pl.concat_arr("a", pl.sum("a")).alias("concat_arr(a, sum(a))"),
    ...         pl.concat_arr("a", pl.max("a")).alias("concat_arr(a, max(a))"),
    ...     )
    ... )
    shape: (3, 4)
    ┌──────┬──────────────────┬───────────────────────┬───────────────────────┐
    │ a    ┆ concat_arr(a, 0) ┆ concat_arr(a, sum(a)) ┆ concat_arr(a, max(a)) │
    │ ---  ┆ ---              ┆ ---                   ┆ ---                   │
    │ i64  ┆ array[i64, 2]    ┆ array[i64, 2]         ┆ array[i64, 2]         │
    ╞══════╪══════════════════╪═══════════════════════╪═══════════════════════╡
    │ 1    ┆ [1, 0]           ┆ [1, 4]                ┆ [1, 3]                │
    │ 3    ┆ [3, 0]           ┆ [3, 4]                ┆ [3, 3]                │
    │ null ┆ [null, 0]        ┆ [null, 4]             ┆ [null, 3]             │
    └──────┴──────────────────┴───────────────────────┴───────────────────────┘
    """
    msg = "`concat_arr` functionality is considered unstable"
    issue_unstable_warning(msg)

    exprs = parse_into_list_of_expressions(exprs, *more_exprs)
    return wrap_expr(plr.concat_arr(exprs))


@overload
def struct(
    *exprs: IntoExpr | Iterable[IntoExpr],
    schema: SchemaDict | None = ...,
    eager: Literal[False] = ...,
    **named_exprs: IntoExpr,
) -> Expr: ...


@overload
def struct(
    *exprs: IntoExpr | Iterable[IntoExpr],
    schema: SchemaDict | None = ...,
    eager: Literal[True],
    **named_exprs: IntoExpr,
) -> Series: ...


@overload
def struct(
    *exprs: IntoExpr | Iterable[IntoExpr],
    schema: SchemaDict | None = ...,
    eager: bool,
    **named_exprs: IntoExpr,
) -> Expr | Series: ...


def struct(
    *exprs: IntoExpr | Iterable[IntoExpr],
    schema: SchemaDict | None = None,
    eager: bool = False,
    **named_exprs: IntoExpr,
) -> Expr | Series:
    """
    Collect columns into a struct column.

    Parameters
    ----------
    *exprs
        Column(s) to collect into a struct column, specified as positional arguments.
        Accepts expression input. Strings are parsed as column names,
        other non-expression inputs are parsed as literals.
    schema
        Optional schema that explicitly defines the struct field dtypes. If no columns
        or expressions are provided, schema keys are used to define columns.
    eager
        Evaluate immediately and return a `Series`. If set to `False` (default),
        return an expression instead.
    **named_exprs
        Additional columns to collect into the struct column, specified as keyword
        arguments. The columns will be renamed to the keyword used.

    Examples
    --------
    Collect all columns of a dataframe into a struct by passing `pl.all()`.

    >>> df = pl.DataFrame(
    ...     {
    ...         "int": [1, 2],
    ...         "str": ["a", "b"],
    ...         "bool": [True, None],
    ...         "list": [[1, 2], [3]],
    ...     }
    ... )
    >>> df.select(pl.struct(pl.all()).alias("my_struct"))
    shape: (2, 1)
    ┌─────────────────────┐
    │ my_struct           │
    │ ---                 │
    │ struct[4]           │
    ╞═════════════════════╡
    │ {1,"a",true,[1, 2]} │
    │ {2,"b",null,[3]}    │
    └─────────────────────┘

    Collect selected columns into a struct by either passing a list of columns, or by
    specifying each column as a positional argument.

    >>> df.select(pl.struct("int", False).alias("my_struct"))
    shape: (2, 1)
    ┌───────────┐
    │ my_struct │
    │ ---       │
    │ struct[2] │
    ╞═══════════╡
    │ {1,false} │
    │ {2,false} │
    └───────────┘

    Use keyword arguments to easily name each struct field.

    >>> df.select(pl.struct(p="int", q="bool").alias("my_struct")).schema
    Schema({'my_struct': Struct({'p': Int64, 'q': Boolean})})
    """
    pyexprs = parse_into_list_of_expressions(*exprs, **named_exprs)

    if schema:
        if not exprs and not named_exprs:
            # no columns or expressions provided; create one from schema keys
            expr = wrap_expr(
                plr.as_struct(parse_into_list_of_expressions(list(schema.keys())))
            )
        else:
            expr = wrap_expr(plr.as_struct(pyexprs))
        expr = expr.cast(Struct(schema), strict=False)
    else:
        expr = wrap_expr(plr.as_struct(pyexprs))

    if eager:
        return F.select(expr).to_series()
    else:
        return expr


def concat_str(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    separator: str = "",
    ignore_nulls: bool = False,
) -> Expr:
    """
    Horizontally concatenate columns into a single string column.

    Operates in linear time.

    Parameters
    ----------
    exprs
        Columns to concatenate into a single string column. Accepts expression input.
        Strings are parsed as column names, other non-expression inputs are parsed as
        literals. Non-`String` columns are cast to `String`.
    *more_exprs
        Additional columns to concatenate into a single string column, specified as
        positional arguments.
    separator
        String that will be used to separate the values of each column.
    ignore_nulls
        Ignore null values (default is ``False``).

        If set to ``False``, null values will be propagated.
        if the row contains any null values, the output is null.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": ["dogs", "cats", None],
    ...         "c": ["play", "swim", "walk"],
    ...     }
    ... )
    >>> df.with_columns(
    ...     pl.concat_str(
    ...         [
    ...             pl.col("a") * 2,
    ...             pl.col("b"),
    ...             pl.col("c"),
    ...         ],
    ...         separator=" ",
    ...     ).alias("full_sentence"),
    ... )
    shape: (3, 4)
    ┌─────┬──────┬──────┬───────────────┐
    │ a   ┆ b    ┆ c    ┆ full_sentence │
    │ --- ┆ ---  ┆ ---  ┆ ---           │
    │ i64 ┆ str  ┆ str  ┆ str           │
    ╞═════╪══════╪══════╪═══════════════╡
    │ 1   ┆ dogs ┆ play ┆ 2 dogs play   │
    │ 2   ┆ cats ┆ swim ┆ 4 cats swim   │
    │ 3   ┆ null ┆ walk ┆ null          │
    └─────┴──────┴──────┴───────────────┘
    """
    exprs = parse_into_list_of_expressions(exprs, *more_exprs)
    return wrap_expr(plr.concat_str(exprs, separator, ignore_nulls))


def format(f_string: str, *args: Expr | str) -> Expr:
    """
    Format expressions as a string.

    Parameters
    ----------
    f_string
        A string that with placeholders.
        For example: "hello_{}" or "{}_world
    args
        Expression(s) that fill the placeholders

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": ["a", "b", "c"],
    ...         "b": [1, 2, 3],
    ...     }
    ... )
    >>> df.select(
    ...     [
    ...         pl.format("foo_{}_bar_{}", pl.col("a"), "b").alias("fmt"),
    ...     ]
    ... )
    shape: (3, 1)
    ┌─────────────┐
    │ fmt         │
    │ ---         │
    │ str         │
    ╞═════════════╡
    │ foo_a_bar_1 │
    │ foo_b_bar_2 │
    │ foo_c_bar_3 │
    └─────────────┘
    """
    exprs = [parse_into_expression(arg) for arg in args]
    return wrap_expr(plr.PyExpr.str_format(f_string, exprs))
