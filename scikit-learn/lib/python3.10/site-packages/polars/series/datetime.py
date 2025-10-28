from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.deprecation import deprecate_nonkeyword_arguments, deprecated
from polars._utils.unstable import unstable
from polars._utils.wrap import wrap_s
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    import datetime as dt
    import sys
    from collections.abc import Iterable

    from polars import Series
    from polars._plr import PySeries
    from polars._typing import (
        Ambiguous,
        EpochTimeUnit,
        IntoExpr,
        IntoExprColumn,
        NonExistent,
        Roll,
        TemporalLiteral,
        TimeUnit,
    )

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


@expr_dispatch
class DateTimeNameSpace:
    """Series.dt namespace."""

    _accessor = "dt"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def __getitem__(self, item: int) -> dt.date | dt.datetime | dt.timedelta:
        s = wrap_s(self._s)
        return s[item]

    @unstable()
    @deprecate_nonkeyword_arguments(allowed_args=["self", "n"], version="1.27.0")
    def add_business_days(
        self,
        n: int | IntoExpr,
        week_mask: Iterable[bool] = (True, True, True, True, True, False, False),
        holidays: Iterable[dt.date] = (),
        roll: Roll = "raise",
    ) -> Series:
        """
        Offset by `n` business days.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.27.0
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
        Series
            Data type is preserved.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series("start", [date(2020, 1, 1), date(2020, 1, 2)])
        >>> s.dt.add_business_days(5)
        shape: (2,)
        Series: 'start' [date]
        [
                2020-01-08
                2020-01-09
        ]

        You can pass a custom weekend - for example, if you only take Sunday off:

        >>> week_mask = (True, True, True, True, True, True, False)
        >>> s.dt.add_business_days(5, week_mask=week_mask)
        shape: (2,)
        Series: 'start' [date]
        [
                2020-01-07
                2020-01-08
        ]

        You can also pass a list of holidays:

        >>> from datetime import date
        >>> holidays = [date(2020, 1, 3), date(2020, 1, 6)]
        >>> s.dt.add_business_days(5, holidays=holidays)
        shape: (2,)
        Series: 'start' [date]
        [
                2020-01-10
                2020-01-13
        ]

        Roll all dates forwards to the next business day:

        >>> s = pl.Series("start", [date(2020, 1, 5), date(2020, 1, 6)])
        >>> s.dt.add_business_days(0, roll="forward")
        shape: (2,)
        Series: 'start' [date]
        [
                2020-01-06
                2020-01-06
        ]
        """

    def min(self) -> dt.date | dt.datetime | dt.timedelta | None:
        """
        Return minimum as Python datetime.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series([date(2001, 1, 1), date(2001, 1, 2), date(2001, 1, 3)])
        >>> s.dt.min()
        datetime.date(2001, 1, 1)
        """
        return wrap_s(self._s).min()  # type: ignore[return-value]

    def max(self) -> dt.date | dt.datetime | dt.timedelta | None:
        """
        Return maximum as Python datetime.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series([date(2001, 1, 1), date(2001, 1, 2), date(2001, 1, 3)])
        >>> s.dt.max()
        datetime.date(2001, 1, 3)
        """
        return wrap_s(self._s).max()  # type: ignore[return-value]

    @deprecated("`Series.dt.median` is deprecated; use `Series.median` instead.")
    def median(self) -> TemporalLiteral | None:
        """
        Return median as python DateTime.

        .. deprecated:: 1.0.0
            Use the `Series.median` method instead.

        Examples
        --------
        >>> from datetime import date, datetime
        >>> s = pl.Series([date(2001, 1, 1), date(2001, 1, 2)])
        >>> s.dt.median()  # doctest: +SKIP
        datetime.datetime(2001, 1, 1, 12, 0)
        >>> date = pl.datetime_range(
        ...     datetime(2001, 1, 1), datetime(2001, 1, 3), "1d", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        >>> date.dt.median()  # doctest: +SKIP
        datetime.datetime(2001, 1, 2, 0, 0)
        """
        return self._s.median()

    @deprecated("`Series.dt.mean` is deprecated; use `Series.mean` instead.")
    def mean(self) -> TemporalLiteral | None:
        """
        Return mean as python DateTime.

        .. deprecated:: 1.0.0
            Use the `Series.mean` method instead.

        Examples
        --------
        >>> from datetime import date, datetime
        >>> s = pl.Series([date(2001, 1, 1), date(2001, 1, 2)])
        >>> s.dt.mean()  # doctest: +SKIP
        datetime.datetime(2001, 1, 1, 12, 0)
        >>> s = pl.Series(
        ...     [datetime(2001, 1, 1), datetime(2001, 1, 2), datetime(2001, 1, 3)]
        ... )
        >>> s.dt.mean()  # doctest: +SKIP
        datetime.datetime(2001, 1, 2, 0, 0)
        """
        return self._s.mean()

    def to_string(self, format: str | None = None) -> Series:
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
        >>> from datetime import datetime
        >>> s = pl.Series(
        ...     "dtm",
        ...     [
        ...         datetime(1999, 12, 31, 6, 12, 30, 800),
        ...         datetime(2020, 7, 5, 10, 20, 45, 12345),
        ...         datetime(2077, 10, 20, 18, 25, 10, 999999),
        ...     ],
        ... )

        Default for temporal dtypes (if not specifying a format string) is ISO8601:

        >>> s.dt.to_string()  # or s.dt.to_string("iso")
        shape: (3,)
        Series: 'dtm' [str]
        [
            "1999-12-31 06:12:30.000800"
            "2020-07-05 10:20:45.012345"
            "2077-10-20 18:25:10.999999"
        ]

        For `Datetime` specifically you can choose between "iso" (where the date and
        time components are ISO, separated by a space) and "iso:strict" (where these
        components are separated by a "T"):

        >>> s.dt.to_string("iso:strict")
        shape: (3,)
        Series: 'dtm' [str]
        [
            "1999-12-31T06:12:30.000800"
            "2020-07-05T10:20:45.012345"
            "2077-10-20T18:25:10.999999"
        ]

        The output can be customized by using a strftime-compatible format string:

        >>> s.dt.to_string("%d/%m/%y")
        shape: (3,)
        Series: 'dtm' [str]
        [
            "31/12/99"
            "05/07/20"
            "20/10/77"
        ]

        If you're interested in using day or month names, you can use
        the `'%A'` and/or `'%B'` format strings:

        >>> s.dt.to_string("%A")
        shape: (3,)
        Series: 'dtm' [str]
        [
            "Friday"
            "Sunday"
            "Wednesday"
        ]

        >>> s.dt.to_string("%B")
        shape: (3,)
        Series: 'dtm' [str]
        [
            "December"
            "July"
            "October"
        ]
        """

    def strftime(self, format: str) -> Series:
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
        to_string : The identical Series method for which `strftime` is an alias.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.Series(
        ...     "datetime",
        ...     [datetime(2020, 3, 1), datetime(2020, 4, 1), datetime(2020, 5, 1)],
        ... )
        >>> s.dt.strftime("%Y/%m/%d")
        shape: (3,)
        Series: 'datetime' [str]
        [
            "2020/03/01"
            "2020/04/01"
            "2020/05/01"
        ]

        If you're interested in the day name / month name, you can use
        `'%A'` / `'%B'`:

        >>> s.dt.strftime("%A")
        shape: (3,)
        Series: 'datetime' [str]
        [
                "Sunday"
                "Wednesday"
                "Friday"
        ]

        >>> s.dt.strftime("%B")
        shape: (3,)
        Series: 'datetime' [str]
        [
                "March"
                "April"
                "May"
        ]
        """
        return self.to_string(format)

    def millennium(self) -> Series:
        """
        Extract the millennium from underlying representation.

        Applies to Date and Datetime columns.

        Returns the millennium number in the calendar date.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series(
        ...     "dt",
        ...     [
        ...         date(999, 12, 31),
        ...         date(1897, 5, 7),
        ...         date(2000, 1, 1),
        ...         date(2001, 7, 5),
        ...         date(3002, 10, 20),
        ...     ],
        ... )
        >>> s.dt.millennium()
        shape: (5,)
        Series: 'dt' [i32]
        [
            1
            2
            2
            3
            4
        ]
        """

    def century(self) -> Series:
        """
        Extract the century from underlying representation.

        Applies to Date and Datetime columns.

        Returns the century number in the calendar date.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series(
        ...     "dt",
        ...     [
        ...         date(999, 12, 31),
        ...         date(1897, 5, 7),
        ...         date(2000, 1, 1),
        ...         date(2001, 7, 5),
        ...         date(3002, 10, 20),
        ...     ],
        ... )
        >>> s.dt.century()
        shape: (5,)
        Series: 'dt' [i32]
        [
            10
            19
            20
            21
            31
        ]
        """

    def year(self) -> Series:
        """
        Extract the year from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the year number in the calendar date.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series("date", [date(2001, 1, 1), date(2002, 1, 1)])
        >>> s.dt.year()
        shape: (2,)
        Series: 'date' [i32]
        [
                2001
                2002
        ]
        """

    @unstable()
    def is_business_day(
        self,
        *,
        week_mask: Iterable[bool] = (True, True, True, True, True, False, False),
        holidays: Iterable[dt.date] = (),
    ) -> Series:
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
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series([date(2020, 1, 3), date(2020, 1, 5)])
        >>> s.dt.is_business_day()
        shape: (2,)
        Series: '' [bool]
        [
            true
            false
        ]

        You can pass a custom weekend - for example, if you only take Sunday off:

        >>> week_mask = (True, True, True, True, True, True, False)
        >>> s.dt.is_business_day(week_mask=week_mask)
        shape: (2,)
        Series: '' [bool]
        [
            true
            false
        ]

        You can also pass a list of holidays:

        >>> from datetime import date
        >>> holidays = [date(2020, 1, 3), date(2020, 1, 6)]
        >>> s.dt.is_business_day(holidays=holidays)
        shape: (2,)
        Series: '' [bool]
        [
            false
            false
        ]
        """

    def is_leap_year(self) -> Series:
        """
        Determine whether the year of the underlying date representation is a leap year.

        Applies to Date and Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series(
        ...     "date", [date(2000, 1, 1), date(2001, 1, 1), date(2002, 1, 1)]
        ... )
        >>> s.dt.is_leap_year()
        shape: (3,)
        Series: 'date' [bool]
        [
                true
                false
                false
        ]
        """

    def iso_year(self) -> Series:
        """
        Extract ISO year from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the year number according to the ISO standard.
        This may not correspond with the calendar year.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> dt = datetime(2022, 1, 1, 7, 8, 40)
        >>> pl.Series([dt]).dt.iso_year()
        shape: (1,)
        Series: '' [i32]
        [
                2021
        ]
        """

    def quarter(self) -> Series:
        """
        Extract quarter from underlying Date representation.

        Applies to Date and Datetime columns.

        Returns the quarter ranging from 1 to 4.

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> date = pl.date_range(
        ...     date(2001, 1, 1), date(2001, 4, 1), interval="1mo", eager=True
        ... ).alias("date")
        >>> date.dt.quarter()
        shape: (4,)
        Series: 'date' [i8]
        [
                1
                1
                1
                2
        ]
        """

    def month(self) -> Series:
        """
        Extract the month from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the month number starting from 1.
        The return value ranges from 1 to 12.

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> date = pl.date_range(
        ...     date(2001, 1, 1), date(2001, 4, 1), interval="1mo", eager=True
        ... ).alias("date")
        >>> date.dt.month()
        shape: (4,)
        Series: 'date' [i8]
        [
                1
                2
                3
                4
        ]
        """

    def days_in_month(self) -> Series:
        """
        Extract the number of days in the month from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the number of days in the month.
        The return value ranges from 28 to 31.

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        See Also
        --------
        month
        is_leap_year

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series(
        ...     "date", [date(2001, 1, 1), date(2001, 2, 1), date(2000, 2, 1)]
        ... )
        >>> s.dt.days_in_month()
        shape: (3,)
        Series: 'date' [i8]
        [
                31
                28
                29
        ]
        """

    def week(self) -> Series:
        """
        Extract the week from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the ISO week number starting from 1.
        The return value ranges from 1 to 53. (The last week of year differs by years.)

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> date = pl.date_range(
        ...     date(2001, 1, 1), date(2001, 4, 1), interval="1mo", eager=True
        ... ).alias("date")
        >>> date.dt.week()
        shape: (4,)
        Series: 'date' [i8]
        [
                1
                5
                9
                13
        ]
        """

    def weekday(self) -> Series:
        """
        Extract the week day from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the ISO weekday number where monday = 1 and sunday = 7

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.date_range(date(2001, 1, 1), date(2001, 1, 7), eager=True).alias(
        ...     "date"
        ... )
        >>> s.dt.weekday()
        shape: (7,)
        Series: 'date' [i8]
        [
                1
                2
                3
                4
                5
                6
                7
        ]
        """

    def day(self) -> Series:
        """
        Extract the day from the underlying date representation.

        Applies to Date and Datetime columns.

        Returns the day of month starting from 1.
        The return value ranges from 1 to 31. (The last day of month differs by months.)

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.date_range(
        ...     date(2001, 1, 1), date(2001, 1, 9), interval="2d", eager=True
        ... ).alias("date")
        >>> s.dt.day()
        shape: (5,)
        Series: 'date' [i8]
        [
                1
                3
                5
                7
                9
        ]
        """

    def ordinal_day(self) -> Series:
        """
        Extract ordinal day from underlying date representation.

        Applies to Date and Datetime columns.

        Returns the day of year starting from 1.
        The return value ranges from 1 to 366. (The last day of year differs by years.)

        Returns
        -------
        Series
            Series of data type :class:`Int16`.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.date_range(
        ...     date(2001, 1, 1), date(2001, 3, 1), interval="1mo", eager=True
        ... ).alias("date")
        >>> s.dt.ordinal_day()
        shape: (3,)
        Series: 'date' [i16]
        [
                1
                32
                60
        ]
        """

    def time(self) -> Series:
        """
        Extract (local) time.

        Applies to Date/Datetime/Time columns.

        Returns
        -------
        Series
            Series of data type :class:`Time`.

        Examples
        --------
        >>> from datetime import datetime
        >>> ser = pl.Series([datetime(2021, 1, 2, 5)]).dt.replace_time_zone(
        ...     "Asia/Kathmandu"
        ... )
        >>> ser
        shape: (1,)
        Series: '' [datetime[μs, Asia/Kathmandu]]
        [
                2021-01-02 05:00:00 +0545
        ]
        >>> ser.dt.time()
        shape: (1,)
        Series: '' [time]
        [
                05:00:00
        ]
        """

    def date(self) -> Series:
        """
        Extract (local) date.

        Applies to Date/Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Date`.

        Examples
        --------
        >>> from datetime import datetime
        >>> ser = pl.Series([datetime(2021, 1, 2, 5)]).dt.replace_time_zone(
        ...     "Asia/Kathmandu"
        ... )
        >>> ser
        shape: (1,)
        Series: '' [datetime[μs, Asia/Kathmandu]]
        [
                2021-01-02 05:00:00 +0545
        ]
        >>> ser.dt.date()
        shape: (1,)
        Series: '' [date]
        [
                2021-01-02
        ]
        """

    @deprecated(
        "`Series.dt.datetime` is deprecated; "
        "use `Series.dt.replace_time_zone(None)` instead."
    )
    def datetime(self) -> Series:
        """
        Extract (local) datetime.

        .. deprecated:: 0.20.4
            Use `dt.replace_time_zone(None)` instead.

        Applies to Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Datetime`.

        Examples
        --------
        >>> from datetime import datetime
        >>> ser = pl.Series([datetime(2021, 1, 2, 5)]).dt.replace_time_zone(
        ...     "Asia/Kathmandu"
        ... )
        >>> ser
        shape: (1,)
        Series: '' [datetime[μs, Asia/Kathmandu]]
        [
                2021-01-02 05:00:00 +0545
        ]
        >>> ser.dt.datetime()  # doctest: +SKIP
        shape: (1,)
        Series: '' [datetime[μs]]
        [
                2021-01-02 05:00:00
        ]
        """

    def hour(self) -> Series:
        """
        Extract the hour from the underlying DateTime representation.

        Applies to Datetime columns.

        Returns the hour number from 0 to 23.

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 3)
        >>> date = pl.datetime_range(start, stop, interval="1h", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (4,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 01:00:00
                2001-01-01 02:00:00
                2001-01-01 03:00:00
        ]
        >>> date.dt.hour()
        shape: (4,)
        Series: 'datetime' [i8]
        [
                0
                1
                2
                3
        ]
        """

    def minute(self) -> Series:
        """
        Extract the minutes from the underlying DateTime representation.

        Applies to Datetime columns.

        Returns the minute number from 0 to 59.

        Returns
        -------
        Series
            Series of data type :class:`Int8`.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 0, 4, 0)
        >>> date = pl.datetime_range(start, stop, interval="2m", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:02:00
                2001-01-01 00:04:00
        ]
        >>> date.dt.minute()
        shape: (3,)
        Series: 'datetime' [i8]
        [
                0
                2
                4
        ]
        """

    def second(self, *, fractional: bool = False) -> Series:
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
        Series
            Series of data type :class:`Int8` or :class:`Float64`.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.Series(
        ...     "datetime",
        ...     [
        ...         datetime(2000, 1, 1, 0, 0, 0, 456789),
        ...         datetime(2000, 1, 1, 0, 0, 3, 111110),
        ...         datetime(2000, 1, 1, 0, 0, 5, 765431),
        ...     ],
        ... )
        >>> s.dt.second()
        shape: (3,)
        Series: 'datetime' [i8]
        [
                0
                3
                5
        ]
        >>> s.dt.second(fractional=True)
        shape: (3,)
        Series: 'datetime' [f64]
        [
                0.456789
                3.11111
                5.765431
        ]
        """

    def millisecond(self) -> Series:
        """
        Extract the milliseconds from the underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 0, 0, 4)
        >>> s = pl.datetime_range(start, stop, interval="500ms", eager=True).alias(
        ...     "datetime"
        ... )
        >>> s.dt.millisecond()
        shape: (9,)
        Series: 'datetime' [i32]
        [
                0
                500
                0
                500
                0
                500
                0
                500
                0
        ]
        """

    def microsecond(self) -> Series:
        """
        Extract the microseconds from the underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 0, 0, 4)
        >>> date = pl.datetime_range(start, stop, interval="500ms", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:00:00.500
                2001-01-01 00:00:01
                2001-01-01 00:00:01.500
                2001-01-01 00:00:02
                2001-01-01 00:00:02.500
                2001-01-01 00:00:03
                2001-01-01 00:00:03.500
                2001-01-01 00:00:04
        ]
        >>> date.dt.microsecond()
        shape: (9,)
        Series: 'datetime' [i32]
        [
                0
                500000
                0
                500000
                0
                500000
                0
                500000
                0
        ]
        """

    def nanosecond(self) -> Series:
        """
        Extract the nanoseconds from the underlying DateTime representation.

        Applies to Datetime columns.

        Returns
        -------
        Series
            Series of data type :class:`Int32`.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 0, 0, 4)
        >>> date = pl.datetime_range(start, stop, interval="500ms", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:00:00.500
                2001-01-01 00:00:01
                2001-01-01 00:00:01.500
                2001-01-01 00:00:02
                2001-01-01 00:00:02.500
                2001-01-01 00:00:03
                2001-01-01 00:00:03.500
                2001-01-01 00:00:04
        ]
        >>> date.dt.nanosecond()
        shape: (9,)
        Series: 'datetime' [i32]
        [
                0
                500000000
                0
                500000000
                0
                500000000
                0
                500000000
                0
        ]
        """

    def timestamp(self, time_unit: TimeUnit = "us") -> Series:
        """
        Return a timestamp in the given time unit.

        Parameters
        ----------
        time_unit : {'us', 'ns', 'ms'}
            Time unit.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 3)
        >>> date = pl.datetime_range(start, stop, interval="1d", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        >>> date.dt.timestamp().alias("timestamp_us")
        shape: (3,)
        Series: 'timestamp_us' [i64]
        [
                978307200000000
                978393600000000
                978480000000000
        ]
        >>> date.dt.timestamp("ns").alias("timestamp_ns")
        shape: (3,)
        Series: 'timestamp_ns' [i64]
        [
                978307200000000000
                978393600000000000
                978480000000000000
        ]
        """

    def epoch(self, time_unit: EpochTimeUnit = "us") -> Series:
        """
        Get the time passed since the Unix EPOCH in the give time unit.

        Parameters
        ----------
        time_unit : {'us', 'ns', 'ms', 's', 'd'}
            Unit of time.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 3)
        >>> date = pl.datetime_range(start, stop, interval="1d", eager=True).alias(
        ...     "datetime"
        ... )
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        >>> date.dt.epoch().alias("epoch_ns")
        shape: (3,)
        Series: 'epoch_ns' [i64]
        [
                978307200000000
                978393600000000
                978480000000000
        ]
        >>> date.dt.epoch(time_unit="s").alias("epoch_s")
        shape: (3,)
        Series: 'epoch_s' [i64]
        [
                978307200
                978393600
                978480000
        ]
        """

    def with_time_unit(self, time_unit: TimeUnit) -> Series:
        """
        Set time unit a Series of dtype Datetime or Duration.

        .. deprecated:: 0.20.5
            First cast to `Int64` and then cast to the desired data type.

        This does not modify underlying data, and should be used to fix an incorrect
        time unit.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms'}
            Unit of time for the `Datetime` or `Duration` Series.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.Series(
        ...     "datetime",
        ...     [datetime(2001, 1, 1), datetime(2001, 1, 2), datetime(2001, 1, 3)],
        ...     dtype=pl.Datetime(time_unit="ns"),
        ... )
        >>> s.dt.with_time_unit("us")  # doctest: +SKIP
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                +32971-04-28 00:00:00
                +32974-01-22 00:00:00
                +32976-10-18 00:00:00
        ]
        """

    def cast_time_unit(self, time_unit: TimeUnit) -> Series:
        """
        Cast the underlying data to another time unit. This may lose precision.

        Parameters
        ----------
        time_unit : {'ns', 'us', 'ms'}
            Unit of time for the `Datetime` Series.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 3)
        >>> date = pl.datetime_range(start, stop, "1d", eager=True).alias("datetime")
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        >>> date.dt.cast_time_unit("ms").alias("time_unit_ms")
        shape: (3,)
        Series: 'time_unit_ms' [datetime[ms]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        >>> date.dt.cast_time_unit("ns").alias("time_unit_ns")
        shape: (3,)
        Series: 'time_unit_ns' [datetime[ns]]
        [
                2001-01-01 00:00:00
                2001-01-02 00:00:00
                2001-01-03 00:00:00
        ]
        """

    def convert_time_zone(self, time_zone: str) -> Series:
        """
        Convert to given time zone for a Series of type Datetime.

        Parameters
        ----------
        time_zone
            Time zone for the `Datetime` Series.

        Notes
        -----
        If converting from a time-zone-naive datetime, then conversion will happen
        as if converting from UTC, regardless of your system's time zone.

        Examples
        --------
        >>> from datetime import datetime
        >>> start = datetime(2020, 3, 1)
        >>> stop = datetime(2020, 5, 1)
        >>> date = pl.datetime_range(
        ...     start, stop, "1mo", time_zone="UTC", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs, UTC]]
        [
                2020-03-01 00:00:00 UTC
                2020-04-01 00:00:00 UTC
                2020-05-01 00:00:00 UTC
        ]
        >>> date = date.dt.convert_time_zone("Europe/London").alias("London")
        >>> date
        shape: (3,)
        Series: 'London' [datetime[μs, Europe/London]]
        [
            2020-03-01 00:00:00 GMT
            2020-04-01 01:00:00 BST
            2020-05-01 01:00:00 BST
        ]
        """

    def replace_time_zone(
        self,
        time_zone: str | None,
        *,
        ambiguous: Ambiguous | Series = "raise",
        non_existent: NonExistent = "raise",
    ) -> Series:
        """
        Replace time zone for a Series of type Datetime.

        Different from `convert_time_zone`, this will also modify
        the underlying timestamp and will ignore the original time zone.

        Parameters
        ----------
        time_zone
            Time zone for the `Datetime` Series. Pass `None` to unset time zone.
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
        ...         )
        ...         .alias("datetime")
        ...         .dt.convert_time_zone(time_zone="Europe/London"),
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
        ...         "ambiguous": ["earliest", "earliest", "earliest", "latest"],
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
        │ 2018-10-28 02:30:00 ┆ earliest  ┆ 2018-10-28 02:30:00 CEST      │
        │ 2018-10-28 02:00:00 ┆ latest    ┆ 2018-10-28 02:00:00 CET       │
        └─────────────────────┴───────────┴───────────────────────────────┘
        """

    def total_days(self, *, fractional: bool = False) -> Series:
        """
        Extract the total days from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the day.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 3, 1), datetime(2020, 5, 1), "1mo", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-03-01 00:00:00
                2020-04-01 00:00:00
                2020-05-01 00:00:00
        ]
        >>> date.diff().dt.total_days()
        shape: (3,)
        Series: 'datetime' [i64]
        [
                null
                31
                30
        ]
        """

    def total_hours(self, *, fractional: bool = False) -> Series:
        """
        Extract the total hours from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the hour.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1), datetime(2020, 1, 4), "1d", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (4,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-02 00:00:00
                2020-01-03 00:00:00
                2020-01-04 00:00:00
        ]
        >>> date.diff().dt.total_hours()
        shape: (4,)
        Series: 'datetime' [i64]
        [
                null
                24
                24
                24
        ]
        """

    def total_minutes(self, *, fractional: bool = False) -> Series:
        """
        Extract the total minutes from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the minute.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1), datetime(2020, 1, 4), "1d", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (4,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-02 00:00:00
                2020-01-03 00:00:00
                2020-01-04 00:00:00
        ]
        >>> date.diff().dt.total_minutes()
        shape: (4,)
        Series: 'datetime' [i64]
        [
                null
                1440
                1440
                1440
        ]
        """

    def total_seconds(self, *, fractional: bool = False) -> Series:
        """
        Extract the total seconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the second.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1), datetime(2020, 1, 1, 0, 4, 0), "1m", eager=True
        ... ).alias("datetime")
        >>> date
        shape: (5,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-01 00:01:00
                2020-01-01 00:02:00
                2020-01-01 00:03:00
                2020-01-01 00:04:00
        ]
        >>> date.diff().dt.total_seconds()
        shape: (5,)
        Series: 'datetime' [i64]
        [
                null
                60
                60
                60
                60
        ]
        """

    def total_milliseconds(self, *, fractional: bool = False) -> Series:
        """
        Extract the total milliseconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the millisecond.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1),
        ...     datetime(2020, 1, 1, 0, 0, 1, 0),
        ...     "1ms",
        ...     eager=True,
        ... ).alias("datetime")[:3]
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-01 00:00:00.001
                2020-01-01 00:00:00.002
        ]
        >>> date.diff().dt.total_milliseconds()
        shape: (3,)
        Series: 'datetime' [i64]
        [
                null
                1
                1
        ]
        """

    def total_microseconds(self, *, fractional: bool = False) -> Series:
        """
        Extract the total microseconds from a Duration type.

        Parameters
        ----------
        fractional
            Whether to include the fractional component of the microsecond.

        Returns
        -------
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1),
        ...     datetime(2020, 1, 1, 0, 0, 1, 0),
        ...     "1ms",
        ...     eager=True,
        ... ).alias("datetime")[:3]
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-01 00:00:00.001
                2020-01-01 00:00:00.002
        ]
        >>> date.diff().dt.total_microseconds()
        shape: (3,)
        Series: 'datetime' [i64]
        [
                null
                1000
                1000
        ]
        """

    def total_nanoseconds(self, *, fractional: bool = False) -> Series:
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
        Series
            Series of data type :class:`.Int64` or :class:`.Float64` if
            `fractional` is set.

        Examples
        --------
        >>> from datetime import datetime
        >>> date = pl.datetime_range(
        ...     datetime(2020, 1, 1),
        ...     datetime(2020, 1, 1, 0, 0, 1, 0),
        ...     "1ms",
        ...     eager=True,
        ... ).alias("datetime")[:3]
        >>> date
        shape: (3,)
        Series: 'datetime' [datetime[μs]]
        [
                2020-01-01 00:00:00
                2020-01-01 00:00:00.001
                2020-01-01 00:00:00.002
        ]
        >>> date.diff().dt.total_nanoseconds()
        shape: (3,)
        Series: 'datetime' [i64]
        [
                null
                1000000
                1000000
        ]
        """

    def offset_by(self, by: str | IntoExprColumn) -> Series:
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

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".

        Returns
        -------
        Series
            Series of data type :class:`Date` or :class:`Datetime`.

        Examples
        --------
        >>> from datetime import datetime
        >>> dates = pl.datetime_range(
        ...     datetime(2000, 1, 1), datetime(2005, 1, 1), "1y", eager=True
        ... ).alias("datetime")
        >>> dates
        shape: (6,)
        Series: 'datetime' [datetime[μs]]
        [
                2000-01-01 00:00:00
                2001-01-01 00:00:00
                2002-01-01 00:00:00
                2003-01-01 00:00:00
                2004-01-01 00:00:00
                2005-01-01 00:00:00
        ]
        >>> dates.dt.offset_by("1y").alias("date_plus_1y")
        shape: (6,)
        Series: 'date_plus_1y' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2002-01-01 00:00:00
                2003-01-01 00:00:00
                2004-01-01 00:00:00
                2005-01-01 00:00:00
                2006-01-01 00:00:00
        ]
        >>> dates.dt.offset_by("-1y2mo").alias("date_minus_1y_2mon")
        shape: (6,)
        Series: 'date_minus_1y_2mon' [datetime[μs]]
        [
                1998-11-01 00:00:00
                1999-11-01 00:00:00
                2000-11-01 00:00:00
                2001-11-01 00:00:00
                2002-11-01 00:00:00
                2003-11-01 00:00:00
        ]
        """

    def truncate(self, every: str | dt.timedelta | IntoExprColumn) -> Series:
        """
        Divide the date/ datetime range into buckets.

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
        The `every` argument is created with the
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
        Series
            Series of data type :class:`Date` or :class:`Datetime`.

        Examples
        --------
        >>> from datetime import timedelta, datetime
        >>> s = pl.datetime_range(
        ...     datetime(2001, 1, 1),
        ...     datetime(2001, 1, 2),
        ...     timedelta(minutes=165),
        ...     eager=True,
        ... ).alias("datetime")
        >>> s
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 02:45:00
            2001-01-01 05:30:00
            2001-01-01 08:15:00
            2001-01-01 11:00:00
            2001-01-01 13:45:00
            2001-01-01 16:30:00
            2001-01-01 19:15:00
            2001-01-01 22:00:00
        ]
        >>> s.dt.truncate("1h")
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 02:00:00
            2001-01-01 05:00:00
            2001-01-01 08:00:00
            2001-01-01 11:00:00
            2001-01-01 13:00:00
            2001-01-01 16:00:00
            2001-01-01 19:00:00
            2001-01-01 22:00:00
        ]

        >>> s = pl.datetime_range(
        ...     datetime(2001, 1, 1), datetime(2001, 1, 1, 1), "10m", eager=True
        ... ).alias("datetime")
        >>> s
        shape: (7,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:10:00
                2001-01-01 00:20:00
                2001-01-01 00:30:00
                2001-01-01 00:40:00
                2001-01-01 00:50:00
                2001-01-01 01:00:00
        ]
        >>> s.dt.truncate("30m")
        shape: (7,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:00:00
                2001-01-01 00:00:00
                2001-01-01 00:30:00
                2001-01-01 00:30:00
                2001-01-01 00:30:00
                2001-01-01 01:00:00
        ]
        """

    def round(self, every: str | dt.timedelta | IntoExprColumn) -> Series:
        """
        Divide the date/ datetime range into buckets.

        - Each date/datetime in the first half of the interval
          is mapped to the start of its bucket.
        - Each date/datetime in the second half of the interval
          is mapped to the end of its bucket.
        - Half-way points are mapped to the start of their bucket.

        Ambiguous results are localized using the DST offset of the original timestamp -
        for example, rounding `'2022-11-06 01:20:00 CST'` by `'1h'` results in
        `'2022-11-06 01:00:00 CST'`, whereas rounding `'2022-11-06 01:20:00 CDT'` by
        `'1h'` results in `'2022-11-06 01:00:00 CDT'`.

        Parameters
        ----------
        every
            Every interval start and period length

        Returns
        -------
        Series
            Series of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        The `every` argument is created with the
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

        Examples
        --------
        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.datetime_range(
        ...     start, stop, timedelta(minutes=165), eager=True
        ... ).alias("datetime")
        >>> s
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 02:45:00
            2001-01-01 05:30:00
            2001-01-01 08:15:00
            2001-01-01 11:00:00
            2001-01-01 13:45:00
            2001-01-01 16:30:00
            2001-01-01 19:15:00
            2001-01-01 22:00:00
        ]
        >>> s.dt.round("1h")
        shape: (9,)
        Series: 'datetime' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 03:00:00
            2001-01-01 06:00:00
            2001-01-01 08:00:00
            2001-01-01 11:00:00
            2001-01-01 14:00:00
            2001-01-01 17:00:00
            2001-01-01 19:00:00
            2001-01-01 22:00:00
        ]
        >>> round_str = s.dt.round("1h")
        >>> round_td = s.dt.round(timedelta(hours=1))
        >>> round_str.equals(round_td)
        True

        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 1, 1)
        >>> s = pl.datetime_range(start, stop, "10m", eager=True).alias("datetime")
        >>> s.dt.round("30m")
        shape: (7,)
        Series: 'datetime' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 00:00:00
                2001-01-01 00:30:00
                2001-01-01 00:30:00
                2001-01-01 00:30:00
                2001-01-01 01:00:00
                2001-01-01 01:00:00
        ]
        """

    def combine(self, time: dt.time | Series, time_unit: TimeUnit = "us") -> Series:
        """
        Create a naive Datetime from an existing Date/Datetime expression and a Time.

        If the underlying expression is a Datetime then its time component is replaced,
        and if it is a Date then a new Datetime is created by combining the two values.

        Parameters
        ----------
        time
            A python time literal or Series of the same length as this Series.
        time_unit : {'ns', 'us', 'ms'}
            Unit of time.

        Examples
        --------
        >>> from datetime import datetime, time
        >>> s = pl.Series(
        ...     "dtm",
        ...     [datetime(2022, 12, 31, 10, 30, 45), datetime(2023, 7, 5, 23, 59, 59)],
        ... )
        >>> s.dt.combine(time(1, 2, 3, 456000))
        shape: (2,)
        Series: 'dtm' [datetime[μs]]
        [
            2022-12-31 01:02:03.456
            2023-07-05 01:02:03.456
        ]
        """

    def month_start(self) -> Series:
        """
        Roll backward to the first day of the month.

        Returns
        -------
        Series
            Series of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        If you're coming from pandas, you can think of this as a vectorised version
        of `pandas.tseries.offsets.MonthBegin().rollback(datetime)`.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.datetime_range(
        ...     datetime(2000, 1, 2, 2), datetime(2000, 4, 2, 2), "1mo", eager=True
        ... ).alias("datetime")
        >>> s.dt.month_start()
        shape: (4,)
        Series: 'datetime' [datetime[μs]]
        [
                2000-01-01 02:00:00
                2000-02-01 02:00:00
                2000-03-01 02:00:00
                2000-04-01 02:00:00
        ]
        """

    def month_end(self) -> Series:
        """
        Roll forward to the last day of the month.

        Returns
        -------
        Series
            Series of data type :class:`Date` or :class:`Datetime`.

        Notes
        -----
        If you're coming from pandas, you can think of this as a vectorised version
        of `pandas.tseries.offsets.MonthEnd().rollforward(datetime)`.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.datetime_range(
        ...     datetime(2000, 1, 2, 2), datetime(2000, 4, 2, 2), "1mo", eager=True
        ... ).alias("datetime")
        >>> s.dt.month_end()
        shape: (4,)
        Series: 'datetime' [datetime[μs]]
        [
                2000-01-31 02:00:00
                2000-02-29 02:00:00
                2000-03-31 02:00:00
                2000-04-30 02:00:00
        ]
        """

    def base_utc_offset(self) -> Series:
        """
        Base offset from UTC.

        This is usually constant for all datetimes in a given time zone, but
        may vary in the rare case that a country switches time zone, like
        Samoa (Apia) did at the end of 2011.

        Returns
        -------
        Series
            Series of data type :class:`Duration`.

        See Also
        --------
        Series.dt.dst_offset : Additional offset currently in effect.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.datetime_range(
        ...     datetime(2011, 12, 29),
        ...     datetime(2012, 1, 1),
        ...     "2d",
        ...     time_zone="Pacific/Apia",
        ...     eager=True,
        ... ).alias("datetime")
        >>> s
        shape: (2,)
        Series: 'datetime' [datetime[μs, Pacific/Apia]]
        [
                2011-12-29 00:00:00 -10
                2011-12-31 00:00:00 +14
        ]
        >>> s.dt.base_utc_offset()
        shape: (2,)
        Series: 'datetime' [duration[ms]]
        [
                -11h
                13h
        ]
        """

    def dst_offset(self) -> Series:
        """
        Additional offset currently in effect (typically due to daylight saving time).

        Returns
        -------
        Series
            Series of data type :class:`Duration`.

        See Also
        --------
        Series.dt.base_utc_offset : Base offset from UTC.

        Examples
        --------
        >>> from datetime import datetime
        >>> s = pl.datetime_range(
        ...     datetime(2020, 10, 25),
        ...     datetime(2020, 10, 26),
        ...     time_zone="Europe/London",
        ...     eager=True,
        ... ).alias("datetime")
        >>> s
        shape: (2,)
        Series: 'datetime' [datetime[μs, Europe/London]]
        [
                2020-10-25 00:00:00 BST
                2020-10-26 00:00:00 GMT
        ]
        >>> s.dt.dst_offset()
        shape: (2,)
        Series: 'datetime' [duration[ms]]
        [
                1h
                0ms
        ]
        """

    def replace(
        self,
        *,
        year: int | Series | None = None,
        month: int | Series | None = None,
        day: int | Series | None = None,
        hour: int | Series | None = None,
        minute: int | Series | None = None,
        second: int | Series | None = None,
        microsecond: int | Series | None = None,
        ambiguous: Ambiguous | Series = "raise",
    ) -> Series:
        """
        Replace time unit.

        Parameters
        ----------
        year
            Literal or Series.
        month
            Literal or Series, ranging from 1-12.
        day
            Literal or Series, ranging from 1-31.
        hour
            Literal or Series, ranging from 0-23.
        minute
            Literal or Series, ranging from 0-59.
        second
            Literal or Series, ranging from 0-59.
        microsecond
            Literal or Series, ranging from 0-999999.
        ambiguous
            Determine how to deal with ambiguous datetimes:

            - `'raise'` (default): raise
            - `'earliest'`: use the earliest datetime
            - `'latest'`: use the latest datetime
            - `'null'`: set to null

        Returns
        -------
        Series
            Series of data type :class:`Date` or :class:`Datetime` with the specified
            time units replaced.

        Examples
        --------
        >>> from datetime import date
        >>> s = pl.Series("date", [date(2013, 1, 1), date(2024, 1, 2)])
        >>> s.dt.replace(year=1800)
        shape: (2,)
        Series: 'date' [date]
        [
                1800-01-01
                1800-01-02
        ]
        """
