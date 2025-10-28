from __future__ import annotations

import warnings
from collections.abc import Mapping
from typing import TYPE_CHECKING

import polars._reexport as pl
from polars import functions as F
from polars._utils.deprecation import deprecate_nonkeyword_arguments, deprecated
from polars._utils.parse import parse_into_expression
from polars._utils.unstable import unstable
from polars._utils.various import (
    find_stacklevel,
    issue_warning,
    no_default,
    qualified_type_name,
)
from polars._utils.wrap import wrap_expr
from polars.datatypes import Date, Datetime, Int64, Time, parse_into_datatype_expr
from polars.exceptions import ChronoFormatWarning

if TYPE_CHECKING:
    import sys

    from polars import Expr
    from polars._typing import (
        Ambiguous,
        IntoExpr,
        IntoExprColumn,
        PolarsDataType,
        PolarsIntegerType,
        PolarsTemporalType,
        TimeUnit,
        TransferEncoding,
        UnicodeForm,
    )
    from polars._utils.various import NoDefault

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


class ExprStringNameSpace:
    """Namespace for string related expressions."""

    _accessor = "str"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def to_date(
        self,
        format: str | None = None,
        *,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
    ) -> Expr:
        """
        Convert a String column into a Date column.

        Parameters
        ----------
        format
            Format to use for conversion. Refer to the `chrono crate documentation
            <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            for the full specification. Example: `"%Y-%m-%d"`.
            If set to None (default), the format is inferred from the data.
        strict
            Raise an error if any conversion fails.
        exact
            Require an exact format match. If False, allow the format to match anywhere
            in the target string.

            .. note::
                Using `exact=False` introduces a performance penalty - cleaning your
                data beforehand will almost certainly be more performant.
        cache
            Use a cache of unique, converted dates to apply the conversion.

        Examples
        --------
        >>> s = pl.Series(["2020/01/01", "2020/02/01", "2020/03/01"])
        >>> s.str.to_date()
        shape: (3,)
        Series: '' [date]
        [
                2020-01-01
                2020-02-01
                2020-03-01
        ]
        """
        _validate_format_argument(format)
        return wrap_expr(self._pyexpr.str_to_date(format, strict, exact, cache))

    def to_datetime(
        self,
        format: str | None = None,
        *,
        time_unit: TimeUnit | None = None,
        time_zone: str | None = None,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
        ambiguous: Ambiguous | Expr = "raise",
    ) -> Expr:
        """
        Convert a String column into a Datetime column.

        Parameters
        ----------
        format
            Format to use for conversion. Refer to the `chrono crate documentation
            <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            for the full specification. Example: `"%Y-%m-%d %H:%M:%S"`.
            If set to None (default), the format is inferred from the data.
        time_unit : {None, 'us', 'ns', 'ms'}
            Unit of time for the resulting Datetime column. If set to None (default),
            the time unit is inferred from the format string if given, eg:
            `"%F %T%.3f"` => `Datetime("ms")`. If no fractional second component is
            found, the default is `"us"`.
        time_zone
            Time zone for the resulting Datetime column. Rules are:

            - If inputs are tz-naive and `time_zone` is None, the result time zone is
              `None`.
            - If inputs are offset-aware and `time_zone` is None, inputs are converted
              to `'UTC'` and the result time zone is `'UTC'`.
            - If inputs are offset-aware and `time_zone` is given, inputs are converted
              to `time_zone` and the result time zone is `time_zone`.
            - If inputs are tz-naive and `time_zone` is given, input time zones are
              replaced with (not converted to!) `time_zone`, and the result time zone
              is `time_zone`.
        strict
            Raise an error if any conversion fails.
        exact
            Require an exact format match. If False, allow the format to match anywhere
            in the target string.

            .. note::
                Using `exact=False` introduces a performance penalty - cleaning your
                data beforehand will almost certainly be more performant.
        cache
            Use a cache of unique, converted datetimes to apply the conversion.
        ambiguous
            Determine how to deal with ambiguous datetimes:

            - `'raise'` (default): raise
            - `'earliest'`: use the earliest datetime
            - `'latest'`: use the latest datetime
            - `'null'`: set to null

        Examples
        --------
        >>> s = pl.Series(["2020-01-01 01:00Z", "2020-01-01 02:00Z"])
        >>> s.str.to_datetime("%Y-%m-%d %H:%M%#z")
        shape: (2,)
        Series: '' [datetime[μs, UTC]]
        [
                2020-01-01 01:00:00 UTC
                2020-01-01 02:00:00 UTC
        ]
        """
        _validate_format_argument(format)
        if not isinstance(ambiguous, pl.Expr):
            ambiguous = F.lit(ambiguous)
        return wrap_expr(
            self._pyexpr.str_to_datetime(
                format,
                time_unit,
                time_zone,
                strict,
                exact,
                cache,
                ambiguous._pyexpr,
            )
        )

    def to_time(
        self,
        format: str | None = None,
        *,
        strict: bool = True,
        cache: bool = True,
    ) -> Expr:
        """
        Convert a String column into a Time column.

        Parameters
        ----------
        format
            Format to use for conversion. Refer to the `chrono crate documentation
            <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            for the full specification. Example: `"%H:%M:%S"`.
            If set to None (default), the format is inferred from the data.
        strict
            Raise an error if any conversion fails.
        cache
            Use a cache of unique, converted times to apply the conversion.

        Examples
        --------
        >>> s = pl.Series(["01:00", "02:00", "03:00"])
        >>> s.str.to_time("%H:%M")
        shape: (3,)
        Series: '' [time]
        [
                01:00:00
                02:00:00
                03:00:00
        ]
        """
        _validate_format_argument(format)
        return wrap_expr(self._pyexpr.str_to_time(format, strict, cache))

    def strptime(
        self,
        dtype: PolarsTemporalType,
        format: str | None = None,
        *,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
        ambiguous: Ambiguous | Expr = "raise",
    ) -> Expr:
        """
        Convert a String column into a Date/Datetime/Time column.

        Parameters
        ----------
        dtype
            The data type to convert into. Can be either Date, Datetime, or Time.
        format
            Format to use for conversion. Refer to the `chrono crate documentation
            <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            for the full specification. Example: `"%Y-%m-%d %H:%M:%S"`.
            If set to None (default), the format is inferred from the data.
        strict
            Raise an error if any conversion fails.
        exact
            Require an exact format match. If False, allow the format to match anywhere
            in the target string. Conversion to the Time type is always exact.

            .. note::
                Using `exact=False` introduces a performance penalty - cleaning your
                data beforehand will almost certainly be more performant.
        cache
            Use a cache of unique, converted dates to apply the datetime conversion.
        ambiguous
            Determine how to deal with ambiguous datetimes:

            - `'raise'` (default): raise
            - `'earliest'`: use the earliest datetime
            - `'latest'`: use the latest datetime
            - `'null'`: set to null

        Notes
        -----
        When converting to a Datetime type, the time unit is inferred from the format
        string if given, eg: `"%F %T%.3f"` => `Datetime("ms")`. If no fractional
        second component is found, the default is `"us"`.

        Examples
        --------
        Dealing with a consistent format:

        >>> s = pl.Series(["2020-01-01 01:00Z", "2020-01-01 02:00Z"])
        >>> s.str.strptime(pl.Datetime, "%Y-%m-%d %H:%M%#z")
        shape: (2,)
        Series: '' [datetime[μs, UTC]]
        [
                2020-01-01 01:00:00 UTC
                2020-01-01 02:00:00 UTC
        ]

        Dealing with different formats.

        >>> s = pl.Series(
        ...     "date",
        ...     [
        ...         "2021-04-22",
        ...         "2022-01-04 00:00:00",
        ...         "01/31/22",
        ...         "Sun Jul  8 00:34:60 2001",
        ...     ],
        ... )
        >>> s.to_frame().select(
        ...     pl.coalesce(
        ...         pl.col("date").str.strptime(pl.Date, "%F", strict=False),
        ...         pl.col("date").str.strptime(pl.Date, "%F %T", strict=False),
        ...         pl.col("date").str.strptime(pl.Date, "%D", strict=False),
        ...         pl.col("date").str.strptime(pl.Date, "%c", strict=False),
        ...     )
        ... ).to_series()
        shape: (4,)
        Series: 'date' [date]
        [
                2021-04-22
                2022-01-04
                2022-01-31
                2001-07-08
        ]
        """
        if dtype == Date:
            return self.to_date(format, strict=strict, exact=exact, cache=cache)
        elif dtype == Datetime:
            time_unit = getattr(dtype, "time_unit", None)
            time_zone = getattr(dtype, "time_zone", None)
            return self.to_datetime(
                format,
                time_unit=time_unit,
                time_zone=time_zone,
                strict=strict,
                exact=exact,
                cache=cache,
                ambiguous=ambiguous,
            )
        elif dtype == Time:
            return self.to_time(format, strict=strict, cache=cache)
        else:
            msg = "`dtype` must be of type {Date, Datetime, Time}"
            raise ValueError(msg)

    @deprecate_nonkeyword_arguments(allowed_args=["self"], version="1.20.0")
    @unstable()
    def to_decimal(self, *, scale: int) -> Expr:
        """
        Convert a String column into a Decimal column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.20.0
            Parameter `inference_length` should now be passed as a keyword argument.

        .. versionchanged:: 1.33.0
            Parameter `inference_length` was removed and `scale` was made non-optional.

        Parameters
        ----------
        scale
            Number of digits after the comma to use for the decimals.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "numbers": [
        ...             "40.12",
        ...             "3420.13",
        ...             "120134.19",
        ...             "3212.98",
        ...             "12.90",
        ...             "143.09",
        ...             "143.9",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(numbers_decimal=pl.col("numbers").str.to_decimal(scale=2))
        shape: (7, 2)
        ┌───────────┬─────────────────┐
        │ numbers   ┆ numbers_decimal │
        │ ---       ┆ ---             │
        │ str       ┆ decimal[38,2]   │
        ╞═══════════╪═════════════════╡
        │ 40.12     ┆ 40.12           │
        │ 3420.13   ┆ 3420.13         │
        │ 120134.19 ┆ 120134.19       │
        │ 3212.98   ┆ 3212.98         │
        │ 12.90     ┆ 12.90           │
        │ 143.09    ┆ 143.09          │
        │ 143.9     ┆ 143.90          │
        └───────────┴─────────────────┘
        """
        return wrap_expr(self._pyexpr.str_to_decimal(scale=scale))

    def len_bytes(self) -> Expr:
        """
        Return the length of each string as the number of bytes.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`.

        See Also
        --------
        len_chars

        Notes
        -----
        When working with non-ASCII text, the length in bytes is not the same as the
        length in characters. You may want to use :func:`len_chars` instead.
        Note that :func:`len_bytes` is much more performant (_O(1)_) than
        :func:`len_chars` (_O(n)_).

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["Café", "345", "東京", None]})
        >>> df.with_columns(
        ...     pl.col("a").str.len_bytes().alias("n_bytes"),
        ...     pl.col("a").str.len_chars().alias("n_chars"),
        ... )
        shape: (4, 3)
        ┌──────┬─────────┬─────────┐
        │ a    ┆ n_bytes ┆ n_chars │
        │ ---  ┆ ---     ┆ ---     │
        │ str  ┆ u32     ┆ u32     │
        ╞══════╪═════════╪═════════╡
        │ Café ┆ 5       ┆ 4       │
        │ 345  ┆ 3       ┆ 3       │
        │ 東京 ┆ 6       ┆ 2       │
        │ null ┆ null    ┆ null    │
        └──────┴─────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.str_len_bytes())

    def len_chars(self) -> Expr:
        """
        Return the length of each string as the number of characters.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`.

        See Also
        --------
        len_bytes

        Notes
        -----
        When working with ASCII text, use :func:`len_bytes` instead to achieve
        equivalent output with much better performance:
        :func:`len_bytes` runs in _O(1)_, while :func:`len_chars` runs in (_O(n)_).

        A character is defined as a `Unicode scalar value`_. A single character is
        represented by a single byte when working with ASCII text, and a maximum of
        4 bytes otherwise.

        .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["Café", "345", "東京", None]})
        >>> df.with_columns(
        ...     pl.col("a").str.len_chars().alias("n_chars"),
        ...     pl.col("a").str.len_bytes().alias("n_bytes"),
        ... )
        shape: (4, 3)
        ┌──────┬─────────┬─────────┐
        │ a    ┆ n_chars ┆ n_bytes │
        │ ---  ┆ ---     ┆ ---     │
        │ str  ┆ u32     ┆ u32     │
        ╞══════╪═════════╪═════════╡
        │ Café ┆ 4       ┆ 5       │
        │ 345  ┆ 3       ┆ 3       │
        │ 東京 ┆ 2       ┆ 6       │
        │ null ┆ null    ┆ null    │
        └──────┴─────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.str_len_chars())

    def to_uppercase(self) -> Expr:
        """
        Modify strings to their uppercase equivalent.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": ["cat", "dog"]})
        >>> df.with_columns(foo_upper=pl.col("foo").str.to_uppercase())
        shape: (2, 2)
        ┌─────┬───────────┐
        │ foo ┆ foo_upper │
        │ --- ┆ ---       │
        │ str ┆ str       │
        ╞═════╪═══════════╡
        │ cat ┆ CAT       │
        │ dog ┆ DOG       │
        └─────┴───────────┘
        """
        return wrap_expr(self._pyexpr.str_to_uppercase())

    def to_lowercase(self) -> Expr:
        """
        Modify strings to their lowercase equivalent.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": ["CAT", "DOG"]})
        >>> df.with_columns(foo_lower=pl.col("foo").str.to_lowercase())
        shape: (2, 2)
        ┌─────┬───────────┐
        │ foo ┆ foo_lower │
        │ --- ┆ ---       │
        │ str ┆ str       │
        ╞═════╪═══════════╡
        │ CAT ┆ cat       │
        │ DOG ┆ dog       │
        └─────┴───────────┘
        """
        return wrap_expr(self._pyexpr.str_to_lowercase())

    def to_titlecase(self) -> Expr:
        """
        Modify strings to their titlecase equivalent.

        Notes
        -----
        This is a form of case transform where the first letter of each word is
        capitalized, with the rest of the word in lowercase. Non-alphanumeric
        characters define the word boundaries.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "quotes": [
        ...             "'e.t. phone home'",
        ...             "you talkin' to me?",
        ...             "to infinity,and BEYOND!",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     quotes_title=pl.col("quotes").str.to_titlecase(),
        ... )
        shape: (3, 2)
        ┌─────────────────────────┬─────────────────────────┐
        │ quotes                  ┆ quotes_title            │
        │ ---                     ┆ ---                     │
        │ str                     ┆ str                     │
        ╞═════════════════════════╪═════════════════════════╡
        │ 'e.t. phone home'       ┆ 'E.T. Phone Home'       │
        │ you talkin' to me?      ┆ You Talkin' To Me?      │
        │ to infinity,and BEYOND! ┆ To Infinity,And Beyond! │
        └─────────────────────────┴─────────────────────────┘
        """
        return wrap_expr(self._pyexpr.str_to_titlecase())

    def strip_chars(self, characters: IntoExpr = None) -> Expr:
        r"""
        Remove leading and trailing characters.

        Parameters
        ----------
        characters
            The set of characters to be removed. All combinations of this set of
            characters will be stripped from the start and end of the string. If set to
            None (default), all leading and trailing whitespace is removed instead.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [" hello", "\nworld"]})
        >>> df
        shape: (2, 1)
        ┌────────┐
        │ foo    │
        │ ---    │
        │ str    │
        ╞════════╡
        │  hello │
        │        │
        │ world  │
        └────────┘

        >>> df.with_columns(foo_stripped=pl.col("foo").str.strip_chars())
        shape: (2, 2)
        ┌────────┬──────────────┐
        │ foo    ┆ foo_stripped │
        │ ---    ┆ ---          │
        │ str    ┆ str          │
        ╞════════╪══════════════╡
        │  hello ┆ hello        │
        │        ┆ world        │
        │ world  ┆              │
        └────────┴──────────────┘

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so, unless that whitespace is
        also included in the string.

        >>> df.with_columns(foo_stripped=pl.col("foo").str.strip_chars("ow\n"))
        shape: (2, 2)
        ┌────────┬──────────────┐
        │ foo    ┆ foo_stripped │
        │ ---    ┆ ---          │
        │ str    ┆ str          │
        ╞════════╪══════════════╡
        │  hello ┆  hell        │
        │        ┆ rld          │
        │ world  ┆              │
        └────────┴──────────────┘
        """
        characters_pyexpr = parse_into_expression(characters, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_strip_chars(characters_pyexpr))

    def strip_chars_start(self, characters: IntoExpr = None) -> Expr:
        r"""
        Remove leading characters.

        .. note::
            This method strips any characters present in `characters` from the
            start of the input, no matter their order. To strip a prefix (i.e.
            a "word" of characters in a certain order), use
            :func:`strip_prefix` instead.

        Parameters
        ----------
        characters
            The set of characters to be removed. All combinations of this set of
            characters will be stripped from the start of the string. If set to None
            (default), all leading whitespace is removed instead.

        See Also
        --------
        strip_prefix
        strip_chars_end

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [" hello ", "\tworld"]})
        >>> df.with_columns(foo_strip_start=pl.col("foo").str.strip_chars_start())
        shape: (2, 2)
        ┌─────────┬─────────────────┐
        │ foo     ┆ foo_strip_start │
        │ ---     ┆ ---             │
        │ str     ┆ str             │
        ╞═════════╪═════════════════╡
        │  hello  ┆ hello           │
        │   world   ┆ world           │
        └─────────┴─────────────────┘

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so.

        >>> df.with_columns(
        ...     foo_strip_start=pl.col("foo").str.strip_chars_start("wod\t"),
        ... )
        shape: (2, 2)
        ┌─────────┬─────────────────┐
        │ foo     ┆ foo_strip_start │
        │ ---     ┆ ---             │
        │ str     ┆ str             │
        ╞═════════╪═════════════════╡
        │  hello  ┆  hello          │
        │   world   ┆ rld             │
        └─────────┴─────────────────┘

        The order of the provided characters does not matter, they behave like a set.

        >>> pl.DataFrame({"foo": ["aabcdef"]}).with_columns(
        ...     foo_strip_start=pl.col("foo").str.strip_chars_start("cba")
        ... )
        shape: (1, 2)
        ┌─────────┬─────────────────┐
        │ foo     ┆ foo_strip_start │
        │ ---     ┆ ---             │
        │ str     ┆ str             │
        ╞═════════╪═════════════════╡
        │ aabcdef ┆ def             │
        └─────────┴─────────────────┘
        """
        characters_pyexpr = parse_into_expression(characters, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_strip_chars_start(characters_pyexpr))

    def strip_chars_end(self, characters: IntoExpr = None) -> Expr:
        r"""
        Remove trailing characters.

        .. note::
            This method strips any characters present in `characters` from the
            end of the input, no matter their order. To strip a suffix (i.e.
            a "word" of characters in a certain order), use
            :func:`strip_suffix` instead.

        Parameters
        ----------
        characters
            The set of characters to be removed. All combinations of this set of
            characters will be stripped from the end of the string. If set to None
            (default), all trailing whitespace is removed instead.

        See Also
        --------
        strip_suffix
        strip_chars_start

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [" hello", "world\n"]})
        >>> df
        shape: (2, 1)
        ┌────────┐
        │ foo    │
        │ ---    │
        │ str    │
        ╞════════╡
        │  hello │
        │ world  │
        │        │
        └────────┘
        >>> df.with_columns(foo_strip_end=pl.col("foo").str.strip_chars_end())
        shape: (2, 2)
        ┌────────┬───────────────┐
        │ foo    ┆ foo_strip_end │
        │ ---    ┆ ---           │
        │ str    ┆ str           │
        ╞════════╪═══════════════╡
        │  hello ┆  hello        │
        │ world  ┆ world         │
        │        ┆               │
        └────────┴───────────────┘

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so, unless that whitespace is
        also included in the string.

        >>> df.with_columns(foo_strip_end=pl.col("foo").str.strip_chars_end("oldw "))
        shape: (2, 2)
        ┌────────┬───────────────┐
        │ foo    ┆ foo_strip_end │
        │ ---    ┆ ---           │
        │ str    ┆ str           │
        ╞════════╪═══════════════╡
        │  hello ┆  he           │
        │ world  ┆ world         │
        │        ┆               │
        └────────┴───────────────┘

        The order of the provided characters does not matter, they behave like a set.

        >>> pl.DataFrame({"foo": ["abcdeff"]}).with_columns(
        ...     foo_strip_end=pl.col("foo").str.strip_chars_end("fed")
        ... )
        shape: (1, 2)
        ┌─────────┬───────────────┐
        │ foo     ┆ foo_strip_end │
        │ ---     ┆ ---           │
        │ str     ┆ str           │
        ╞═════════╪═══════════════╡
        │ abcdeff ┆ abc           │
        └─────────┴───────────────┘
        """
        characters_pyexpr = parse_into_expression(characters, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_strip_chars_end(characters_pyexpr))

    def strip_prefix(self, prefix: IntoExpr) -> Expr:
        """
        Remove prefix.

        The prefix will be removed from the string exactly once, if found.

        .. note::
            This method strips the exact character sequence provided in
            `prefix` from the start of the input. To strip a set of characters
            in any order, use :func:`strip_chars_start` instead.

        Parameters
        ----------
        prefix
            The prefix to be removed.

        See Also
        --------
        strip_chars_start
        strip_suffix

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["foobar", "foofoobar", "foo", "bar"]})
        >>> df.with_columns(pl.col("a").str.strip_prefix("foo").alias("stripped"))
        shape: (4, 2)
        ┌───────────┬──────────┐
        │ a         ┆ stripped │
        │ ---       ┆ ---      │
        │ str       ┆ str      │
        ╞═══════════╪══════════╡
        │ foobar    ┆ bar      │
        │ foofoobar ┆ foobar   │
        │ foo       ┆          │
        │ bar       ┆ bar      │
        └───────────┴──────────┘
        """
        prefix_pyexpr = parse_into_expression(prefix, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_strip_prefix(prefix_pyexpr))

    def strip_suffix(self, suffix: IntoExpr) -> Expr:
        """
        Remove suffix.

        The suffix will be removed from the string exactly once, if found.

        .. note::
            This method strips the exact character sequence provided in
            `suffix` from the end of the input. To strip a set of characters
            in any order, use :func:`strip_chars_end` instead.

        Parameters
        ----------
        suffix
            The suffix to be removed.

        See Also
        --------
        strip_chars_end
        strip_prefix

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["foobar", "foobarbar", "foo", "bar"]})
        >>> df.with_columns(pl.col("a").str.strip_suffix("bar").alias("stripped"))
        shape: (4, 2)
        ┌───────────┬──────────┐
        │ a         ┆ stripped │
        │ ---       ┆ ---      │
        │ str       ┆ str      │
        ╞═══════════╪══════════╡
        │ foobar    ┆ foo      │
        │ foobarbar ┆ foobar   │
        │ foo       ┆ foo      │
        │ bar       ┆          │
        └───────────┴──────────┘
        """
        suffix_pyexpr = parse_into_expression(suffix, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_strip_suffix(suffix_pyexpr))

    def pad_start(self, length: int | IntoExprColumn, fill_char: str = " ") -> Expr:
        """
        Pad the start of the string until it reaches the given length.

        Parameters
        ----------
        length
            Pad the string until it reaches this length. Strings with length equal to or
            greater than this value are returned as-is.
        fill_char
            The character to pad the string with.

        See Also
        --------
        pad_end
        zfill

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["cow", "monkey", "hippopotamus", None]})
        >>> df.with_columns(padded=pl.col("a").str.pad_start(8, "*"))
        shape: (4, 2)
        ┌──────────────┬──────────────┐
        │ a            ┆ padded       │
        │ ---          ┆ ---          │
        │ str          ┆ str          │
        ╞══════════════╪══════════════╡
        │ cow          ┆ *****cow     │
        │ monkey       ┆ **monkey     │
        │ hippopotamus ┆ hippopotamus │
        │ null         ┆ null         │
        └──────────────┴──────────────┘
        """
        length_pyexpr = parse_into_expression(length)
        if not isinstance(fill_char, str):
            msg = f'"pad_start" expects a `str`, given a {qualified_type_name(fill_char)!r}'
            raise TypeError(msg)
        return wrap_expr(self._pyexpr.str_pad_start(length_pyexpr, fill_char))

    def pad_end(self, length: int | IntoExprColumn, fill_char: str = " ") -> Expr:
        """
        Pad the end of the string until it reaches the given length.

        Parameters
        ----------
        length
            Pad the string until it reaches this length. Strings with length equal to or
            greater than this value are returned as-is. Can be int or expression.
        fill_char
            The character to pad the string with.

        See Also
        --------
        pad_start

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["cow", "monkey", "hippopotamus", None]})
        >>> df.with_columns(padded=pl.col("a").str.pad_end(8, "*"))
        shape: (4, 2)
        ┌──────────────┬──────────────┐
        │ a            ┆ padded       │
        │ ---          ┆ ---          │
        │ str          ┆ str          │
        ╞══════════════╪══════════════╡
        │ cow          ┆ cow*****     │
        │ monkey       ┆ monkey**     │
        │ hippopotamus ┆ hippopotamus │
        │ null         ┆ null         │
        └──────────────┴──────────────┘
        """
        length_pyexpr = parse_into_expression(length)
        if not isinstance(fill_char, str):
            msg = (
                f'"pad_end" expects a `str`, given a {qualified_type_name(fill_char)!r}'
            )
            raise TypeError(msg)
        return wrap_expr(self._pyexpr.str_pad_end(length_pyexpr, fill_char))

    def zfill(self, length: int | IntoExprColumn) -> Expr:
        """
        Pad the start of the string with zeros until it reaches the given length.

        A sign prefix (`-`) is handled by inserting the padding after the sign
        character rather than before.

        Parameters
        ----------
        length
            Pad the string until it reaches this length. Strings with length equal to
            or greater than this value are returned as-is.

        See Also
        --------
        pad_start

        Notes
        -----
        This method is intended for padding numeric strings. If your data contains
        non-ASCII characters, use :func:`pad_start` instead.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [-1, 123, 999999, None]})
        >>> df.with_columns(zfill=pl.col("a").cast(pl.String).str.zfill(4))
        shape: (4, 2)
        ┌────────┬────────┐
        │ a      ┆ zfill  │
        │ ---    ┆ ---    │
        │ i64    ┆ str    │
        ╞════════╪════════╡
        │ -1     ┆ -001   │
        │ 123    ┆ 0123   │
        │ 999999 ┆ 999999 │
        │ null   ┆ null   │
        └────────┴────────┘
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [-1, 123, 999999, None],
        ...         "length": [8, 4, 1, 2],
        ...     }
        ... )
        >>> df.with_columns(zfill=pl.col("a").cast(pl.String).str.zfill("length"))
        shape: (4, 3)
        ┌────────┬────────┬──────────┐
        │ a      ┆ length ┆ zfill    │
        │ ---    ┆ ---    ┆ ---      │
        │ i64    ┆ i64    ┆ str      │
        ╞════════╪════════╪══════════╡
        │ -1     ┆ 8      ┆ -0000001 │
        │ 123    ┆ 4      ┆ 0123     │
        │ 999999 ┆ 1      ┆ 999999   │
        │ null   ┆ 2      ┆ null     │
        └────────┴────────┴──────────┘
        """
        length_pyexpr = parse_into_expression(length)
        return wrap_expr(self._pyexpr.str_zfill(length_pyexpr))

    def contains(
        self, pattern: str | Expr, *, literal: bool = False, strict: bool = True
    ) -> Expr:
        """
        Check if the string contains a substring that matches a pattern.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        literal
            Treat `pattern` as a literal string, not as a regular expression.
        strict
            Raise an error if the underlying pattern is not a valid regex,
            otherwise mask out with a null value.

        Notes
        -----
        To modify regular expression behaviour (such as case-sensitivity) with
        flags, use the inline `(?iLmsuxU)` syntax. For example:

        >>> pl.DataFrame({"s": ["AAA", "aAa", "aaa"]}).with_columns(
        ...     default_match=pl.col("s").str.contains("AA"),
        ...     insensitive_match=pl.col("s").str.contains("(?i)AA"),
        ... )
        shape: (3, 3)
        ┌─────┬───────────────┬───────────────────┐
        │ s   ┆ default_match ┆ insensitive_match │
        │ --- ┆ ---           ┆ ---               │
        │ str ┆ bool          ┆ bool              │
        ╞═════╪═══════════════╪═══════════════════╡
        │ AAA ┆ true          ┆ true              │
        │ aAa ┆ false         ┆ true              │
        │ aaa ┆ false         ┆ true              │
        └─────┴───────────────┴───────────────────┘

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        See Also
        --------
        starts_with : Check if string values start with a substring.
        ends_with : Check if string values end with a substring.
        find: Return the index of the first substring matching a pattern.

        Examples
        --------
        >>> df = pl.DataFrame({"txt": ["Crab", "cat and dog", "rab$bit", None]})
        >>> df.select(
        ...     pl.col("txt"),
        ...     pl.col("txt").str.contains("cat|bit").alias("regex"),
        ...     pl.col("txt").str.contains("rab$", literal=True).alias("literal"),
        ... )
        shape: (4, 3)
        ┌─────────────┬───────┬─────────┐
        │ txt         ┆ regex ┆ literal │
        │ ---         ┆ ---   ┆ ---     │
        │ str         ┆ bool  ┆ bool    │
        ╞═════════════╪═══════╪═════════╡
        │ Crab        ┆ false ┆ false   │
        │ cat and dog ┆ true  ┆ false   │
        │ rab$bit     ┆ true  ┆ true    │
        │ null        ┆ null  ┆ null    │
        └─────────────┴───────┴─────────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_contains(pattern_pyexpr, literal, strict))

    def find(
        self, pattern: str | Expr, *, literal: bool = False, strict: bool = True
    ) -> Expr:
        """
        Return the bytes offset of the first substring matching a pattern.

        If the pattern is not found, returns None.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        literal
            Treat `pattern` as a literal string, not as a regular expression.
        strict
            Raise an error if the underlying pattern is not a valid regex,
            otherwise mask out with a null value.

        Notes
        -----
        To modify regular expression behaviour (such as case-sensitivity) with
        flags, use the inline `(?iLmsuxU)` syntax. For example:

        >>> pl.DataFrame({"s": ["AAA", "aAa", "aaa"]}).with_columns(
        ...     default_match=pl.col("s").str.find("Aa"),
        ...     insensitive_match=pl.col("s").str.find("(?i)Aa"),
        ... )
        shape: (3, 3)
        ┌─────┬───────────────┬───────────────────┐
        │ s   ┆ default_match ┆ insensitive_match │
        │ --- ┆ ---           ┆ ---               │
        │ str ┆ u32           ┆ u32               │
        ╞═════╪═══════════════╪═══════════════════╡
        │ AAA ┆ null          ┆ 0                 │
        │ aAa ┆ 1             ┆ 0                 │
        │ aaa ┆ null          ┆ 0                 │
        └─────┴───────────────┴───────────────────┘

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        See Also
        --------
        contains : Check if the string contains a substring that matches a pattern.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "txt": ["Crab", "Lobster", None, "Crustacean"],
        ...         "pat": ["a[bc]", "b.t", "[aeiuo]", "(?i)A[BC]"],
        ...     }
        ... )

        Find the index of the first substring matching a regex or literal pattern:

        >>> df.select(
        ...     pl.col("txt"),
        ...     pl.col("txt").str.find("a|e").alias("a|e (regex)"),
        ...     pl.col("txt").str.find("e", literal=True).alias("e (lit)"),
        ... )
        shape: (4, 3)
        ┌────────────┬─────────────┬─────────┐
        │ txt        ┆ a|e (regex) ┆ e (lit) │
        │ ---        ┆ ---         ┆ ---     │
        │ str        ┆ u32         ┆ u32     │
        ╞════════════╪═════════════╪═════════╡
        │ Crab       ┆ 2           ┆ null    │
        │ Lobster    ┆ 5           ┆ 5       │
        │ null       ┆ null        ┆ null    │
        │ Crustacean ┆ 5           ┆ 7       │
        └────────────┴─────────────┴─────────┘

        Match against a pattern found in another column or (expression):

        >>> df.with_columns(pl.col("txt").str.find(pl.col("pat")).alias("find_pat"))
        shape: (4, 3)
        ┌────────────┬───────────┬──────────┐
        │ txt        ┆ pat       ┆ find_pat │
        │ ---        ┆ ---       ┆ ---      │
        │ str        ┆ str       ┆ u32      │
        ╞════════════╪═══════════╪══════════╡
        │ Crab       ┆ a[bc]     ┆ 2        │
        │ Lobster    ┆ b.t       ┆ 2        │
        │ null       ┆ [aeiuo]   ┆ null     │
        │ Crustacean ┆ (?i)A[BC] ┆ 5        │
        └────────────┴───────────┴──────────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_find(pattern_pyexpr, literal, strict))

    def ends_with(self, suffix: str | Expr) -> Expr:
        """
        Check if string values end with a substring.

        Parameters
        ----------
        suffix
            Suffix substring.

        See Also
        --------
        contains : Check if the string contains a substring that matches a pattern.
        starts_with : Check if string values start with a substring.

        Examples
        --------
        >>> df = pl.DataFrame({"fruits": ["apple", "mango", None]})
        >>> df.with_columns(
        ...     pl.col("fruits").str.ends_with("go").alias("has_suffix"),
        ... )
        shape: (3, 2)
        ┌────────┬────────────┐
        │ fruits ┆ has_suffix │
        │ ---    ┆ ---        │
        │ str    ┆ bool       │
        ╞════════╪════════════╡
        │ apple  ┆ false      │
        │ mango  ┆ true       │
        │ null   ┆ null       │
        └────────┴────────────┘

        >>> df = pl.DataFrame(
        ...     {"fruits": ["apple", "mango", "banana"], "suffix": ["le", "go", "nu"]}
        ... )
        >>> df.with_columns(
        ...     pl.col("fruits").str.ends_with(pl.col("suffix")).alias("has_suffix"),
        ... )
        shape: (3, 3)
        ┌────────┬────────┬────────────┐
        │ fruits ┆ suffix ┆ has_suffix │
        │ ---    ┆ ---    ┆ ---        │
        │ str    ┆ str    ┆ bool       │
        ╞════════╪════════╪════════════╡
        │ apple  ┆ le     ┆ true       │
        │ mango  ┆ go     ┆ true       │
        │ banana ┆ nu     ┆ false      │
        └────────┴────────┴────────────┘

        Using `ends_with` as a filter condition:

        >>> df.filter(pl.col("fruits").str.ends_with("go"))
        shape: (1, 2)
        ┌────────┬────────┐
        │ fruits ┆ suffix │
        │ ---    ┆ ---    │
        │ str    ┆ str    │
        ╞════════╪════════╡
        │ mango  ┆ go     │
        └────────┴────────┘
        """
        suffix_pyexpr = parse_into_expression(suffix, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_ends_with(suffix_pyexpr))

    def starts_with(self, prefix: str | Expr) -> Expr:
        """
        Check if string values start with a substring.

        Parameters
        ----------
        prefix
            Prefix substring.

        See Also
        --------
        contains : Check if the string contains a substring that matches a pattern.
        ends_with : Check if string values end with a substring.

        Examples
        --------
        >>> df = pl.DataFrame({"fruits": ["apple", "mango", None]})
        >>> df.with_columns(
        ...     pl.col("fruits").str.starts_with("app").alias("has_prefix"),
        ... )
        shape: (3, 2)
        ┌────────┬────────────┐
        │ fruits ┆ has_prefix │
        │ ---    ┆ ---        │
        │ str    ┆ bool       │
        ╞════════╪════════════╡
        │ apple  ┆ true       │
        │ mango  ┆ false      │
        │ null   ┆ null       │
        └────────┴────────────┘

        >>> df = pl.DataFrame(
        ...     {"fruits": ["apple", "mango", "banana"], "prefix": ["app", "na", "ba"]}
        ... )
        >>> df.with_columns(
        ...     pl.col("fruits").str.starts_with(pl.col("prefix")).alias("has_prefix"),
        ... )
        shape: (3, 3)
        ┌────────┬────────┬────────────┐
        │ fruits ┆ prefix ┆ has_prefix │
        │ ---    ┆ ---    ┆ ---        │
        │ str    ┆ str    ┆ bool       │
        ╞════════╪════════╪════════════╡
        │ apple  ┆ app    ┆ true       │
        │ mango  ┆ na     ┆ false      │
        │ banana ┆ ba     ┆ true       │
        └────────┴────────┴────────────┘

        Using `starts_with` as a filter condition:

        >>> df.filter(pl.col("fruits").str.starts_with("app"))
        shape: (1, 2)
        ┌────────┬────────┐
        │ fruits ┆ prefix │
        │ ---    ┆ ---    │
        │ str    ┆ str    │
        ╞════════╪════════╡
        │ apple  ┆ app    │
        └────────┴────────┘
        """
        prefix_pyexpr = parse_into_expression(prefix, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_starts_with(prefix_pyexpr))

    def json_decode(
        self,
        dtype: PolarsDataType | pl.DataTypeExpr,
        *,
        infer_schema_length: int | None = None,
    ) -> Expr:
        """
        Parse string values as JSON.

        Throws an error if invalid JSON strings are encountered.

        Parameters
        ----------
        dtype
            The dtype to cast the extracted value to.
        infer_schema_length
            Deprecated and ignored.

            .. versionchanged:: 1.33.0
                Deprecate `infer_schema_length` and make `dtype` non-optional to
                ensure that the planner can determine the output datatype.

        See Also
        --------
        json_path_match : Extract the first match from a JSON string using the provided
            JSONPath.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"json": ['{"a":1, "b": true}', None, '{"a":2, "b": false}']}
        ... )
        >>> dtype = pl.Struct([pl.Field("a", pl.Int64), pl.Field("b", pl.Boolean)])
        >>> df.with_columns(decoded=pl.col("json").str.json_decode(dtype))
        shape: (3, 2)
        ┌─────────────────────┬───────────┐
        │ json                ┆ decoded   │
        │ ---                 ┆ ---       │
        │ str                 ┆ struct[2] │
        ╞═════════════════════╪═══════════╡
        │ {"a":1, "b": true}  ┆ {1,true}  │
        │ null                ┆ null      │
        │ {"a":2, "b": false} ┆ {2,false} │
        └─────────────────────┴───────────┘
        """
        if dtype is None:
            msg = "`Expr.str.json_decode` needs an explicitly given `dtype` otherwise Polars is not able to determine the output type. If you want to eagerly infer datatype you can use `Series.str.json_decode`."
            raise TypeError(msg)

        if infer_schema_length is not None:
            issue_warning(
                "`Expr.str.json_decode` with `infer_schema_length` is deprecated and has no effect on execution.",
                DeprecationWarning,
            )

        dtype_expr = parse_into_datatype_expr(dtype)._pydatatype_expr
        return wrap_expr(self._pyexpr.str_json_decode(dtype_expr))

    def json_path_match(self, json_path: IntoExprColumn) -> Expr:
        """
        Extract the first match from a JSON string using the provided JSONPath.

        Throws errors if invalid JSON strings are encountered. All return values
        are cast to :class:`String`, regardless of the original value.

        Documentation on the JSONPath standard can be found
        `here <https://goessner.net/articles/JsonPath/>`_.

        Parameters
        ----------
        json_path
            A valid JSONPath query string.

        Returns
        -------
        Expr
            Expression of data type :class:`String`. Contains null values if original
            value is null or the json_path returns nothing.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"json_val": ['{"a":"1"}', None, '{"a":2}', '{"a":2.1}', '{"a":true}']}
        ... )
        >>> df.with_columns(matched=pl.col("json_val").str.json_path_match("$.a"))
        shape: (5, 2)
        ┌────────────┬─────────┐
        │ json_val   ┆ matched │
        │ ---        ┆ ---     │
        │ str        ┆ str     │
        ╞════════════╪═════════╡
        │ {"a":"1"}  ┆ 1       │
        │ null       ┆ null    │
        │ {"a":2}    ┆ 2       │
        │ {"a":2.1}  ┆ 2.1     │
        │ {"a":true} ┆ true    │
        └────────────┴─────────┘
        """
        json_path_pyexpr = parse_into_expression(json_path, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_json_path_match(json_path_pyexpr))

    def decode(self, encoding: TransferEncoding, *, strict: bool = True) -> Expr:
        r"""
        Decode values using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.
        strict
            Raise an error if the underlying value cannot be decoded,
            otherwise mask out with a null value.

        Returns
        -------
        Expr
            Expression of data type :class:`Binary`.

        Examples
        --------
        >>> df = pl.DataFrame({"color": ["000000", "ffff00", "0000ff"]})
        >>> df.with_columns(pl.col("color").str.decode("hex").alias("decoded"))
        shape: (3, 2)
        ┌────────┬─────────────────┐
        │ color  ┆ decoded         │
        │ ---    ┆ ---             │
        │ str    ┆ binary          │
        ╞════════╪═════════════════╡
        │ 000000 ┆ b"\x00\x00\x00" │
        │ ffff00 ┆ b"\xff\xff\x00" │
        │ 0000ff ┆ b"\x00\x00\xff" │
        └────────┴─────────────────┘
        """
        if encoding == "hex":
            return wrap_expr(self._pyexpr.str_hex_decode(strict))
        elif encoding == "base64":
            return wrap_expr(self._pyexpr.str_base64_decode(strict))
        else:
            msg = f"`encoding` must be one of {{'hex', 'base64'}}, got {encoding!r}"
            raise ValueError(msg)

    def encode(self, encoding: TransferEncoding) -> Expr:
        """
        Encode values using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"strings": ["foo", "bar", None]})
        >>> df.with_columns(strings_hex=pl.col("strings").str.encode("hex"))
        shape: (3, 2)
        ┌─────────┬─────────────┐
        │ strings ┆ strings_hex │
        │ ---     ┆ ---         │
        │ str     ┆ str         │
        ╞═════════╪═════════════╡
        │ foo     ┆ 666f6f      │
        │ bar     ┆ 626172      │
        │ null    ┆ null        │
        └─────────┴─────────────┘
        """
        if encoding == "hex":
            return wrap_expr(self._pyexpr.str_hex_encode())
        elif encoding == "base64":
            return wrap_expr(self._pyexpr.str_base64_encode())
        else:
            msg = f"`encoding` must be one of {{'hex', 'base64'}}, got {encoding!r}"
            raise ValueError(msg)

    def extract(self, pattern: IntoExprColumn, group_index: int = 1) -> Expr:
        r"""
        Extract the target capture group from provided patterns.

        Parameters
        ----------
        pattern
            A valid regular expression pattern containing at least one capture group,
            compatible with the `regex crate <https://docs.rs/regex/latest/regex/>`_.
        group_index
            Index of the targeted capture group.
            Group 0 means the whole pattern, the first group begins at index 1.
            Defaults to the first capture group.

        Notes
        -----
        To modify regular expression behaviour (such as multi-line matching)
        with flags, use the inline `(?iLmsuxU)` syntax. For example:

        >>> df = pl.DataFrame(
        ...     data={
        ...         "lines": [
        ...             "I Like\nThose\nOdds",
        ...             "This is\nThe Way",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("lines").str.extract(r"(?m)^(T\w+)", 1).alias("matches"),
        ... )
        shape: (2, 2)
        ┌─────────┬─────────┐
        │ lines   ┆ matches │
        │ ---     ┆ ---     │
        │ str     ┆ str     │
        ╞═════════╪═════════╡
        │ I Like  ┆ Those   │
        │ Those   ┆         │
        │ Odds    ┆         │
        │ This is ┆ This    │
        │ The Way ┆         │
        └─────────┴─────────┘

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        Returns
        -------
        Expr
            Expression of data type :class:`String`. Contains null values if original
            value is null or the regex captures nothing.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "url": [
        ...             "http://vote.com/ballon_dor?error=404&ref=unknown",
        ...             "http://vote.com/ballon_dor?ref=polars&candidate=messi",
        ...             "http://vote.com/ballon_dor?candidate=ronaldo&ref=polars",
        ...         ]
        ...     }
        ... )
        >>> df.select(
        ...     pl.col("url").str.extract(r"candidate=(\w+)", 1).alias("candidate"),
        ...     pl.col("url").str.extract(r"ref=(\w+)", 1).alias("referer"),
        ...     pl.col("url").str.extract(r"error=(\w+)", 1).alias("error"),
        ... )
        shape: (3, 3)
        ┌───────────┬─────────┬───────┐
        │ candidate ┆ referer ┆ error │
        │ ---       ┆ ---     ┆ ---   │
        │ str       ┆ str     ┆ str   │
        ╞═══════════╪═════════╪═══════╡
        │ null      ┆ unknown ┆ 404   │
        │ messi     ┆ polars  ┆ null  │
        │ ronaldo   ┆ polars  ┆ null  │
        └───────────┴─────────┴───────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_extract(pattern_pyexpr, group_index))

    def extract_all(self, pattern: str | Expr) -> Expr:
        r'''
        Extract all matches for the given regex pattern.

        Extract each successive non-overlapping regex match in an individual string
        as a list. If the haystack string is `null`, `null` is returned.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.

        Notes
        -----
        To modify regular expression behaviour (such as "verbose" mode and/or
        case-sensitive matching) with flags, use the inline `(?iLmsuxU)` syntax.
        For example:

        >>> df = pl.DataFrame(
        ...     data={
        ...         "email": [
        ...             "real.email@spam.com",
        ...             "some_account@somewhere.net",
        ...             "abc.def.ghi.jkl@uvw.xyz.co.uk",
        ...         ]
        ...     }
        ... )
        >>> # extract name/domain parts from the addresses, using verbose regex
        >>> df.with_columns(
        ...     pl.col("email")
        ...     .str.extract_all(
        ...         r"""(?xi)   # activate 'verbose' and 'case-insensitive' flags
        ...         [           # (start character group)
        ...           A-Z       # letters
        ...           0-9       # digits
        ...           ._%+\-    # special chars
        ...         ]           # (end character group)
        ...         +           # 'one or more' quantifier
        ...         """
        ...     )
        ...     .list.to_struct(fields=["name", "domain"])
        ...     .alias("email_parts")
        ... ).unnest("email_parts")
        shape: (3, 3)
        ┌───────────────────────────────┬─────────────────┬───────────────┐
        │ email                         ┆ name            ┆ domain        │
        │ ---                           ┆ ---             ┆ ---           │
        │ str                           ┆ str             ┆ str           │
        ╞═══════════════════════════════╪═════════════════╪═══════════════╡
        │ real.email@spam.com           ┆ real.email      ┆ spam.com      │
        │ some_account@somewhere.net    ┆ some_account    ┆ somewhere.net │
        │ abc.def.ghi.jkl@uvw.xyz.co.uk ┆ abc.def.ghi.jkl ┆ uvw.xyz.co.uk │
        └───────────────────────────────┴─────────────────┴───────────────┘

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        Returns
        -------
        Expr
            Expression of data type `List(String)`.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": ["123 bla 45 asd", "xyz 678 910t", "bar", None]})
        >>> df.select(
        ...     pl.col("foo").str.extract_all(r"\d+").alias("extracted_nrs"),
        ... )
        shape: (4, 1)
        ┌────────────────┐
        │ extracted_nrs  │
        │ ---            │
        │ list[str]      │
        ╞════════════════╡
        │ ["123", "45"]  │
        │ ["678", "910"] │
        │ []             │
        │ null           │
        └────────────────┘

        '''
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_extract_all(pattern_pyexpr))

    def extract_groups(self, pattern: str) -> Expr:
        r"""
        Extract all capture groups for the given regex pattern.

        Parameters
        ----------
        pattern
            A valid regular expression pattern containing at least one capture group,
            compatible with the `regex crate <https://docs.rs/regex/latest/regex/>`_.

        Notes
        -----
        All group names are **strings**.

        If your pattern contains unnamed groups, their numerical position is converted
        to a string.

        For example, here we access groups 2 and 3 via the names `"2"` and `"3"`::

            >>> df = pl.DataFrame({"col": ["foo bar baz"]})
            >>> (
            ...     df.with_columns(
            ...         pl.col("col").str.extract_groups(r"(\S+) (\S+) (.+)")
            ...     ).select(pl.col("col").struct["2"], pl.col("col").struct["3"])
            ... )
            shape: (1, 2)
            ┌─────┬─────┐
            │ 2   ┆ 3   │
            │ --- ┆ --- │
            │ str ┆ str │
            ╞═════╪═════╡
            │ bar ┆ baz │
            └─────┴─────┘

        Returns
        -------
        Expr
            Expression of data type :class:`Struct` with fields of data type
            :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "url": [
        ...             "http://vote.com/ballon_dor?candidate=messi&ref=python",
        ...             "http://vote.com/ballon_dor?candidate=weghorst&ref=polars",
        ...             "http://vote.com/ballon_dor?error=404&ref=rust",
        ...         ]
        ...     }
        ... )
        >>> pattern = r"candidate=(?<candidate>\w+)&ref=(?<ref>\w+)"
        >>> df.select(captures=pl.col("url").str.extract_groups(pattern)).unnest(
        ...     "captures"
        ... )
        shape: (3, 2)
        ┌───────────┬────────┐
        │ candidate ┆ ref    │
        │ ---       ┆ ---    │
        │ str       ┆ str    │
        ╞═══════════╪════════╡
        │ messi     ┆ python │
        │ weghorst  ┆ polars │
        │ null      ┆ null   │
        └───────────┴────────┘

        Unnamed groups have their numerical position converted to a string:

        >>> pattern = r"candidate=(\w+)&ref=(\w+)"
        >>> (
        ...     df.with_columns(
        ...         captures=pl.col("url").str.extract_groups(pattern)
        ...     ).with_columns(name=pl.col("captures").struct["1"].str.to_uppercase())
        ... )
        shape: (3, 3)
        ┌─────────────────────────────────┬───────────────────────┬──────────┐
        │ url                             ┆ captures              ┆ name     │
        │ ---                             ┆ ---                   ┆ ---      │
        │ str                             ┆ struct[2]             ┆ str      │
        ╞═════════════════════════════════╪═══════════════════════╪══════════╡
        │ http://vote.com/ballon_dor?can… ┆ {"messi","python"}    ┆ MESSI    │
        │ http://vote.com/ballon_dor?can… ┆ {"weghorst","polars"} ┆ WEGHORST │
        │ http://vote.com/ballon_dor?err… ┆ {null,null}           ┆ null     │
        └─────────────────────────────────┴───────────────────────┴──────────┘
        """
        if not isinstance(pattern, str):
            msg = f'"extract_groups" expects a `str`, given a {qualified_type_name(pattern)!r}'
            raise TypeError(msg)
        return wrap_expr(self._pyexpr.str_extract_groups(pattern))

    def count_matches(self, pattern: str | Expr, *, literal: bool = False) -> Expr:
        r"""
        Count all successive non-overlapping regex matches.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        literal
            Treat `pattern` as a literal string, not as a regular expression.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32`. Returns null if the
            original value is null.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": ["123 bla 45 asd", "xyz 678 910t", "bar", None]})
        >>> df.with_columns(
        ...     pl.col("foo").str.count_matches(r"\d").alias("count_digits"),
        ... )
        shape: (4, 2)
        ┌────────────────┬──────────────┐
        │ foo            ┆ count_digits │
        │ ---            ┆ ---          │
        │ str            ┆ u32          │
        ╞════════════════╪══════════════╡
        │ 123 bla 45 asd ┆ 5            │
        │ xyz 678 910t   ┆ 6            │
        │ bar            ┆ 0            │
        │ null           ┆ null         │
        └────────────────┴──────────────┘

        >>> df = pl.DataFrame({"bar": ["12 dbc 3xy", "cat\\w", "1zy3\\d\\d", None]})
        >>> df.with_columns(
        ...     pl.col("bar")
        ...     .str.count_matches(r"\d", literal=True)
        ...     .alias("count_digits"),
        ... )
        shape: (4, 2)
        ┌────────────┬──────────────┐
        │ bar        ┆ count_digits │
        │ ---        ┆ ---          │
        │ str        ┆ u32          │
        ╞════════════╪══════════════╡
        │ 12 dbc 3xy ┆ 0            │
        │ cat\w      ┆ 0            │
        │ 1zy3\d\d   ┆ 2            │
        │ null       ┆ null         │
        └────────────┴──────────────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_count_matches(pattern_pyexpr, literal))

    def split(self, by: IntoExpr, *, inclusive: bool = False) -> Expr:
        """
        Split the string by a substring.

        Parameters
        ----------
        by
            Substring to split by.
        inclusive
            If True, include the split character/string in the results.

        Examples
        --------
        >>> df = pl.DataFrame({"s": ["foo bar", "foo_bar", "foo_bar_baz"]})
        >>> df.with_columns(
        ...     pl.col("s").str.split(by="_").alias("split"),
        ...     pl.col("s").str.split(by="_", inclusive=True).alias("split_inclusive"),
        ... )
        shape: (3, 3)
        ┌─────────────┬───────────────────────┬─────────────────────────┐
        │ s           ┆ split                 ┆ split_inclusive         │
        │ ---         ┆ ---                   ┆ ---                     │
        │ str         ┆ list[str]             ┆ list[str]               │
        ╞═════════════╪═══════════════════════╪═════════════════════════╡
        │ foo bar     ┆ ["foo bar"]           ┆ ["foo bar"]             │
        │ foo_bar     ┆ ["foo", "bar"]        ┆ ["foo_", "bar"]         │
        │ foo_bar_baz ┆ ["foo", "bar", "baz"] ┆ ["foo_", "bar_", "baz"] │
        └─────────────┴───────────────────────┴─────────────────────────┘

        >>> df = pl.DataFrame(
        ...     {"s": ["foo^bar", "foo_bar", "foo*bar*baz"], "by": ["_", "_", "*"]}
        ... )
        >>> df.with_columns(
        ...     pl.col("s").str.split(by=pl.col("by")).alias("split"),
        ...     pl.col("s")
        ...     .str.split(by=pl.col("by"), inclusive=True)
        ...     .alias("split_inclusive"),
        ... )
        shape: (3, 4)
        ┌─────────────┬─────┬───────────────────────┬─────────────────────────┐
        │ s           ┆ by  ┆ split                 ┆ split_inclusive         │
        │ ---         ┆ --- ┆ ---                   ┆ ---                     │
        │ str         ┆ str ┆ list[str]             ┆ list[str]               │
        ╞═════════════╪═════╪═══════════════════════╪═════════════════════════╡
        │ foo^bar     ┆ _   ┆ ["foo^bar"]           ┆ ["foo^bar"]             │
        │ foo_bar     ┆ _   ┆ ["foo", "bar"]        ┆ ["foo_", "bar"]         │
        │ foo*bar*baz ┆ *   ┆ ["foo", "bar", "baz"] ┆ ["foo*", "bar*", "baz"] │
        └─────────────┴─────┴───────────────────────┴─────────────────────────┘

        Returns
        -------
        Expr
            Expression of data type :class:`String`.
        """
        by_pyexpr = parse_into_expression(by, str_as_lit=True)
        if inclusive:
            return wrap_expr(self._pyexpr.str_split_inclusive(by_pyexpr))
        return wrap_expr(self._pyexpr.str_split(by_pyexpr))

    def split_exact(self, by: IntoExpr, n: int, *, inclusive: bool = False) -> Expr:
        """
        Split the string by a substring using `n` splits.

        Results in a struct of `n+1` fields.

        If it cannot make `n` splits, the remaining field elements will be null.

        Parameters
        ----------
        by
            Substring to split by.
        n
            Number of splits to make.
        inclusive
            If True, include the split character/string in the results.

        Returns
        -------
        Expr
            Expression of data type :class:`Struct` with fields of data type
            :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"x": ["a_1", None, "c", "d_4"]})
        >>> df.with_columns(
        ...     extracted=pl.col("x").str.split_exact("_", 1).alias("fields"),
        ... )
        shape: (4, 2)
        ┌──────┬─────────────┐
        │ x    ┆ extracted   │
        │ ---  ┆ ---         │
        │ str  ┆ struct[2]   │
        ╞══════╪═════════════╡
        │ a_1  ┆ {"a","1"}   │
        │ null ┆ {null,null} │
        │ c    ┆ {"c",null}  │
        │ d_4  ┆ {"d","4"}   │
        └──────┴─────────────┘


        Split string values in column x in exactly 2 parts and assign
        each part to a new column.

        >>> df.with_columns(
        ...     pl.col("x")
        ...     .str.split_exact("_", 1)
        ...     .struct.rename_fields(["first_part", "second_part"])
        ...     .alias("fields")
        ... ).unnest("fields")
        shape: (4, 3)
        ┌──────┬────────────┬─────────────┐
        │ x    ┆ first_part ┆ second_part │
        │ ---  ┆ ---        ┆ ---         │
        │ str  ┆ str        ┆ str         │
        ╞══════╪════════════╪═════════════╡
        │ a_1  ┆ a          ┆ 1           │
        │ null ┆ null       ┆ null        │
        │ c    ┆ c          ┆ null        │
        │ d_4  ┆ d          ┆ 4           │
        └──────┴────────────┴─────────────┘
        """
        by_pyexpr = parse_into_expression(by, str_as_lit=True)
        if inclusive:
            return wrap_expr(self._pyexpr.str_split_exact_inclusive(by_pyexpr, n))
        return wrap_expr(self._pyexpr.str_split_exact(by_pyexpr, n))

    def splitn(self, by: IntoExpr, n: int) -> Expr:
        """
        Split the string by a substring, restricted to returning at most `n` items.

        If the number of possible splits is less than `n-1`, the remaining field
        elements will be null. If the number of possible splits is `n-1` or greater,
        the last (nth) substring will contain the remainder of the string.

        Parameters
        ----------
        by
            Substring to split by.
        n
            Max number of items to return.

        Returns
        -------
        Expr
            Expression of data type :class:`Struct` with fields of data type
            :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"s": ["foo bar", None, "foo-bar", "foo bar baz"]})
        >>> df.with_columns(pl.col("s").str.splitn(" ", 2).alias("fields"))
        shape: (4, 2)
        ┌─────────────┬───────────────────┐
        │ s           ┆ fields            │
        │ ---         ┆ ---               │
        │ str         ┆ struct[2]         │
        ╞═════════════╪═══════════════════╡
        │ foo bar     ┆ {"foo","bar"}     │
        │ null        ┆ {null,null}       │
        │ foo-bar     ┆ {"foo-bar",null}  │
        │ foo bar baz ┆ {"foo","bar baz"} │
        └─────────────┴───────────────────┘

        Split string values in column s in exactly 2 parts and assign
        each part to a new column.

        >>> df.with_columns(
        ...     pl.col("s")
        ...     .str.splitn(" ", 2)
        ...     .struct.rename_fields(["first_part", "second_part"])
        ...     .alias("fields")
        ... ).unnest("fields")
        shape: (4, 3)
        ┌─────────────┬────────────┬─────────────┐
        │ s           ┆ first_part ┆ second_part │
        │ ---         ┆ ---        ┆ ---         │
        │ str         ┆ str        ┆ str         │
        ╞═════════════╪════════════╪═════════════╡
        │ foo bar     ┆ foo        ┆ bar         │
        │ null        ┆ null       ┆ null        │
        │ foo-bar     ┆ foo-bar    ┆ null        │
        │ foo bar baz ┆ foo        ┆ bar baz     │
        └─────────────┴────────────┴─────────────┘
        """
        by_pyexpr = parse_into_expression(by, str_as_lit=True)
        return wrap_expr(self._pyexpr.str_splitn(by_pyexpr, n))

    def replace(
        self,
        pattern: str | Expr,
        value: str | Expr,
        *,
        literal: bool = False,
        n: int = 1,
    ) -> Expr:
        r"""
        Replace first matching regex/literal substring with a new string value.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        value
            String that will replace the matched substring.
        literal
            Treat `pattern` as a literal string, not a regex.
        n
            Number of matches to replace.

        See Also
        --------
        replace_all

        Notes
        -----
        * To modify regular expression behaviour (such as case-sensitivity) with flags,
          use the inline `(?iLmsuxU)` syntax. See the regex crate's section on
          `grouping and flags <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_
          for additional information about the use of inline expression modifiers.

        * The dollar sign (`$`) is a special character related to capture groups; if you
          want to replace some target pattern with characters that include a literal `$`
          you should escape it by doubling it up as `$$`, or set `literal=True` if you
          do not need a full regular expression pattern match. Otherwise, you will be
          referencing a (potentially non-existent) capture group.

          In the example below we need to double up `$` (to represent a literal dollar
          sign, and then refer to the capture group using `$n` or `${n}`, hence the
          three consecutive `$` characters in the replacement value:

          .. code-block:: python

              >>> df = pl.DataFrame({"cost": ["#12.34", "#56.78"]})
              >>> df.with_columns(
              ...     cost_usd=pl.col("cost").str.replace(r"#(\d+)", "$$${1}")
              ... )
              shape: (2, 2)
              ┌────────┬──────────┐
              │ cost   ┆ cost_usd │
              │ ---    ┆ ---      │
              │ str    ┆ str      │
              ╞════════╪══════════╡
              │ #12.34 ┆ $12.34   │
              │ #56.78 ┆ $56.78   │
              └────────┴──────────┘

        Examples
        --------
        >>> df = pl.DataFrame({"id": [1, 2], "text": ["123abc", "abc456"]})
        >>> df.with_columns(pl.col("text").str.replace(r"abc\b", "ABC"))
        shape: (2, 2)
        ┌─────┬────────┐
        │ id  ┆ text   │
        │ --- ┆ ---    │
        │ i64 ┆ str    │
        ╞═════╪════════╡
        │ 1   ┆ 123ABC │
        │ 2   ┆ abc456 │
        └─────┴────────┘

        Capture groups are supported. Use `$1` or `${1}` in the `value` string to refer
        to the first capture group in the `pattern`, `$2` or `${2}` to refer to the
        second capture group, and so on. You can also use *named* capture groups.

        >>> df = pl.DataFrame({"word": ["hat", "hut"]})
        >>> df.with_columns(
        ...     positional=pl.col.word.str.replace("h(.)t", "b${1}d"),
        ...     named=pl.col.word.str.replace("h(?<vowel>.)t", "b${vowel}d"),
        ... )
        shape: (2, 3)
        ┌──────┬────────────┬───────┐
        │ word ┆ positional ┆ named │
        │ ---  ┆ ---        ┆ ---   │
        │ str  ┆ str        ┆ str   │
        ╞══════╪════════════╪═══════╡
        │ hat  ┆ bad        ┆ bad   │
        │ hut  ┆ bud        ┆ bud   │
        └──────┴────────────┴───────┘

        Apply case-insensitive string replacement using the `(?i)` flag.

        >>> df = pl.DataFrame(
        ...     {
        ...         "city": "Philadelphia",
        ...         "season": ["Spring", "Summer", "Autumn", "Winter"],
        ...         "weather": ["Rainy", "Sunny", "Cloudy", "Snowy"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("weather").str.replace(r"(?i)foggy|rainy|cloudy|snowy", "Sunny")
        ... )
        shape: (4, 3)
        ┌──────────────┬────────┬─────────┐
        │ city         ┆ season ┆ weather │
        │ ---          ┆ ---    ┆ ---     │
        │ str          ┆ str    ┆ str     │
        ╞══════════════╪════════╪═════════╡
        │ Philadelphia ┆ Spring ┆ Sunny   │
        │ Philadelphia ┆ Summer ┆ Sunny   │
        │ Philadelphia ┆ Autumn ┆ Sunny   │
        │ Philadelphia ┆ Winter ┆ Sunny   │
        └──────────────┴────────┴─────────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        value_pyexpr = parse_into_expression(value, str_as_lit=True)
        return wrap_expr(
            self._pyexpr.str_replace_n(pattern_pyexpr, value_pyexpr, literal, n)
        )

    def replace_all(
        self, pattern: str | Expr, value: str | Expr, *, literal: bool = False
    ) -> Expr:
        r"""
        Replace all matching regex/literal substrings with a new string value.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        value
            String that will replace the matched substring.
        literal
            Treat `pattern` as a literal string, not a regex.

        See Also
        --------
        replace

        Notes
        -----
        * To modify regular expression behaviour (such as case-sensitivity) with flags,
          use the inline `(?iLmsuxU)` syntax. See the regex crate's section on
          `grouping and flags <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_
          for additional information about the use of inline expression modifiers.

        * The dollar sign (`$`) is a special character related to capture groups; if you
          want to replace some target pattern with characters that include a literal `$`
          you should escape it by doubling it up as `$$`, or set `literal=True` if you
          do not need a full regular expression pattern match. Otherwise, you will be
          referencing a (potentially non-existent) capture group.

          In the example below we need to double up `$` to represent a literal dollar
          sign, otherwise we are referring to a capture group (which may or may not
          exist):

          .. code-block:: python

              >>> df = pl.DataFrame({"text": ["ab12cd34ef", "gh45ij67kl"]})
              >>> df.with_columns(
              ...     # the replacement pattern refers back to the capture group
              ...     text1=pl.col("text").str.replace_all(r"(?<N>\d{2,})", "$N$"),
              ...     # doubling-up the `$` results in it appearing as a literal value
              ...     text2=pl.col("text").str.replace_all(r"(?<N>\d{2,})", "$$N$$"),
              ... )
              shape: (2, 3)
              ┌────────────┬──────────────┬──────────────┐
              │ text       ┆ text1        ┆ text2        │
              │ ---        ┆ ---          ┆ ---          │
              │ str        ┆ str          ┆ str          │
              ╞════════════╪══════════════╪══════════════╡
              │ ab12cd34ef ┆ ab12$cd34$ef ┆ ab$N$cd$N$ef │
              │ gh45ij67kl ┆ gh45$ij67$kl ┆ gh$N$ij$N$kl │
              └────────────┴──────────────┴──────────────┘

        Examples
        --------
        >>> df = pl.DataFrame({"id": [1, 2], "text": ["abcabc", "123a123"]})
        >>> df.with_columns(pl.col("text").str.replace_all("a", "-"))
        shape: (2, 2)
        ┌─────┬─────────┐
        │ id  ┆ text    │
        │ --- ┆ ---     │
        │ i64 ┆ str     │
        ╞═════╪═════════╡
        │ 1   ┆ -bc-bc  │
        │ 2   ┆ 123-123 │
        └─────┴─────────┘

        Capture groups are supported. Use `$1` or `${1}` in the `value` string to refer
        to the first capture group in the `pattern`, `$2` or `${2}` to refer to the
        second capture group, and so on. You can also use *named* capture groups.

        >>> df = pl.DataFrame({"word": ["hat", "hut"]})
        >>> df.with_columns(
        ...     positional=pl.col.word.str.replace_all("h(.)t", "b${1}d"),
        ...     named=pl.col.word.str.replace_all("h(?<vowel>.)t", "b${vowel}d"),
        ... )
        shape: (2, 3)
        ┌──────┬────────────┬───────┐
        │ word ┆ positional ┆ named │
        │ ---  ┆ ---        ┆ ---   │
        │ str  ┆ str        ┆ str   │
        ╞══════╪════════════╪═══════╡
        │ hat  ┆ bad        ┆ bad   │
        │ hut  ┆ bud        ┆ bud   │
        └──────┴────────────┴───────┘

        Apply case-insensitive string replacement using the `(?i)` flag.

        >>> df = pl.DataFrame(
        ...     {
        ...         "city": "Philadelphia",
        ...         "season": ["Spring", "Summer", "Autumn", "Winter"],
        ...         "weather": ["Rainy", "Sunny", "Cloudy", "Snowy"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     # apply case-insensitive string replacement
        ...     pl.col("weather").str.replace_all(
        ...         r"(?i)foggy|rainy|cloudy|snowy", "Sunny"
        ...     )
        ... )
        shape: (4, 3)
        ┌──────────────┬────────┬─────────┐
        │ city         ┆ season ┆ weather │
        │ ---          ┆ ---    ┆ ---     │
        │ str          ┆ str    ┆ str     │
        ╞══════════════╪════════╪═════════╡
        │ Philadelphia ┆ Spring ┆ Sunny   │
        │ Philadelphia ┆ Summer ┆ Sunny   │
        │ Philadelphia ┆ Autumn ┆ Sunny   │
        │ Philadelphia ┆ Winter ┆ Sunny   │
        └──────────────┴────────┴─────────┘
        """
        pattern_pyexpr = parse_into_expression(pattern, str_as_lit=True)
        value_pyexpr = parse_into_expression(value, str_as_lit=True)
        return wrap_expr(
            self._pyexpr.str_replace_all(pattern_pyexpr, value_pyexpr, literal)
        )

    def reverse(self) -> Expr:
        """
        Returns string values in reversed order.

        Examples
        --------
        >>> df = pl.DataFrame({"text": ["foo", "bar", "man\u0303ana"]})
        >>> df.with_columns(pl.col("text").str.reverse().alias("reversed"))
        shape: (3, 2)
        ┌────────┬──────────┐
        │ text   ┆ reversed │
        │ ---    ┆ ---      │
        │ str    ┆ str      │
        ╞════════╪══════════╡
        │ foo    ┆ oof      │
        │ bar    ┆ rab      │
        │ mañana ┆ anañam   │
        └────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.str_reverse())

    def slice(
        self, offset: int | IntoExprColumn, length: int | IntoExprColumn | None = None
    ) -> Expr:
        """
        Extract a substring from each string value.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None` (default), the slice is taken to the
            end of the string.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Notes
        -----
        Both the `offset` and `length` inputs are defined in terms of the number
        of characters in the (UTF8) string. A character is defined as a
        `Unicode scalar value`_. A single character is represented by a single byte
        when working with ASCII text, and a maximum of 4 bytes otherwise.

        .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        Examples
        --------
        >>> df = pl.DataFrame({"s": ["pear", None, "papaya", "dragonfruit"]})
        >>> df.with_columns(pl.col("s").str.slice(-3).alias("slice"))
        shape: (4, 2)
        ┌─────────────┬───────┐
        │ s           ┆ slice │
        │ ---         ┆ ---   │
        │ str         ┆ str   │
        ╞═════════════╪═══════╡
        │ pear        ┆ ear   │
        │ null        ┆ null  │
        │ papaya      ┆ aya   │
        │ dragonfruit ┆ uit   │
        └─────────────┴───────┘

        Using the optional `length` parameter

        >>> df.with_columns(pl.col("s").str.slice(4, length=3).alias("slice"))
        shape: (4, 2)
        ┌─────────────┬───────┐
        │ s           ┆ slice │
        │ ---         ┆ ---   │
        │ str         ┆ str   │
        ╞═════════════╪═══════╡
        │ pear        ┆       │
        │ null        ┆ null  │
        │ papaya      ┆ ya    │
        │ dragonfruit ┆ onf   │
        └─────────────┴───────┘
        """
        offset_pyexpr = parse_into_expression(offset)
        length_pyexpr = parse_into_expression(length)
        return wrap_expr(self._pyexpr.str_slice(offset_pyexpr, length_pyexpr))

    def head(self, n: int | IntoExprColumn) -> Expr:
        """
        Return the first n characters of each string in a String Series.

        Parameters
        ----------
        n
            Length of the slice (integer or expression). Negative indexing is supported;
            see note (2) below.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Notes
        -----
        1) The `n` input is defined in terms of the number of characters in the (UTF8)
           string. A character is defined as a `Unicode scalar value`_. A single
           character is represented by a single byte when working with ASCII text, and a
           maximum of 4 bytes otherwise.

           .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        2) When the `n` input is negative, `head` returns characters up to the `n`th
           from the end of the string. For example, if `n = -3`, then all characters
           except the last three are returned.

        3) If the length of the string has fewer than `n` characters, the full string is
           returned.

        Examples
        --------
        Return up to the first 5 characters:

        >>> df = pl.DataFrame({"s": ["pear", None, "papaya", "dragonfruit"]})
        >>> df.with_columns(pl.col("s").str.head(5).alias("s_head_5"))
        shape: (4, 2)
        ┌─────────────┬──────────┐
        │ s           ┆ s_head_5 │
        │ ---         ┆ ---      │
        │ str         ┆ str      │
        ╞═════════════╪══════════╡
        │ pear        ┆ pear     │
        │ null        ┆ null     │
        │ papaya      ┆ papay    │
        │ dragonfruit ┆ drago    │
        └─────────────┴──────────┘

        Return characters determined by column `n`:

        >>> df = pl.DataFrame(
        ...     {
        ...         "s": ["pear", None, "papaya", "dragonfruit"],
        ...         "n": [3, 4, -2, -5],
        ...     }
        ... )
        >>> df.with_columns(pl.col("s").str.head("n").alias("s_head_n"))
        shape: (4, 3)
        ┌─────────────┬─────┬──────────┐
        │ s           ┆ n   ┆ s_head_n │
        │ ---         ┆ --- ┆ ---      │
        │ str         ┆ i64 ┆ str      │
        ╞═════════════╪═════╪══════════╡
        │ pear        ┆ 3   ┆ pea      │
        │ null        ┆ 4   ┆ null     │
        │ papaya      ┆ -2  ┆ papa     │
        │ dragonfruit ┆ -5  ┆ dragon   │
        └─────────────┴─────┴──────────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.str_head(n_pyexpr))

    def tail(self, n: int | IntoExprColumn) -> Expr:
        """
        Return the last n characters of each string in a String Series.

        Parameters
        ----------
        n
            Length of the slice (integer or expression). Negative indexing is supported;
            see note (2) below.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Notes
        -----
        1) The `n` input is defined in terms of the number of characters in the (UTF8)
           string. A character is defined as a `Unicode scalar value`_. A single
           character is represented by a single byte when working with ASCII text, and a
           maximum of 4 bytes otherwise.

           .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        2) When the `n` input is negative, `tail` returns characters starting from the
           `n`th from the beginning of the string. For example, if `n = -3`, then all
           characters except the first three are returned.

        3) If the length of the string has fewer than `n` characters, the full string is
           returned.

        Examples
        --------
        Return up to the last 5 characters:

        >>> df = pl.DataFrame({"s": ["pear", None, "papaya", "dragonfruit"]})
        >>> df.with_columns(pl.col("s").str.tail(5).alias("s_tail_5"))
        shape: (4, 2)
        ┌─────────────┬──────────┐
        │ s           ┆ s_tail_5 │
        │ ---         ┆ ---      │
        │ str         ┆ str      │
        ╞═════════════╪══════════╡
        │ pear        ┆ pear     │
        │ null        ┆ null     │
        │ papaya      ┆ apaya    │
        │ dragonfruit ┆ fruit    │
        └─────────────┴──────────┘

        Return characters determined by column `n`:

        >>> df = pl.DataFrame(
        ...     {
        ...         "s": ["pear", None, "papaya", "dragonfruit"],
        ...         "n": [3, 4, -2, -5],
        ...     }
        ... )
        >>> df.with_columns(pl.col("s").str.tail("n").alias("s_tail_n"))
        shape: (4, 3)
        ┌─────────────┬─────┬──────────┐
        │ s           ┆ n   ┆ s_tail_n │
        │ ---         ┆ --- ┆ ---      │
        │ str         ┆ i64 ┆ str      │
        ╞═════════════╪═════╪══════════╡
        │ pear        ┆ 3   ┆ ear      │
        │ null        ┆ 4   ┆ null     │
        │ papaya      ┆ -2  ┆ paya     │
        │ dragonfruit ┆ -5  ┆ nfruit   │
        └─────────────┴─────┴──────────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.str_tail(n_pyexpr))

    @deprecated(
        '`str.explode` is deprecated; use `str.split("").explode()` instead.'
        " Note that empty strings will result in null instead of being preserved."
        " To get the exact same behavior, split first and then use a `pl.when...then...otherwise`"
        " expression to handle the empty list before exploding."
    )
    def explode(self) -> Expr:
        """
        Returns a column with a separate row for every string character.

        .. deprecated:: 0.20.31
            Use the `.str.split("").explode()` method instead. Note that empty strings
            will result in null instead of being preserved. To get the exact same
            behavior, split first and then use a `pl.when...then...otherwise`
            expression to handle the empty list before exploding.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": ["foo", "bar"]})
        >>> df.select(pl.col("a").str.explode())  # doctest: +SKIP
        shape: (6, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ str │
        ╞═════╡
        │ f   │
        │ o   │
        │ o   │
        │ b   │
        │ a   │
        │ r   │
        └─────┘
        """
        split = self.split("")
        return F.when(split.ne_missing([])).then(split).otherwise([""]).explode()

    def to_integer(
        self,
        *,
        base: int | IntoExprColumn = 10,
        dtype: PolarsIntegerType = Int64,
        strict: bool = True,
    ) -> Expr:
        """
        Convert a String column into an Int64 column with base radix.

        Parameters
        ----------
        base
            Positive integer or expression which is the base of the string
            we are parsing.
            Default: 10.
        dtype
            Integer data type to cast the result to.
            Default: Int64.
        strict
            Bool, Default=True will raise any ParseError or overflow as ComputeError.
            False silently convert to Null.

        Returns
        -------
        Expr
            Expression of data type :class:`Int64`.

        Examples
        --------
        >>> df = pl.DataFrame({"bin": ["110", "101", "010", "invalid"]})
        >>> df.with_columns(
        ...     parsed=pl.col("bin").str.to_integer(
        ...         base=2, dtype=pl.Int32, strict=False
        ...     )
        ... )
        shape: (4, 2)
        ┌─────────┬────────┐
        │ bin     ┆ parsed │
        │ ---     ┆ ---    │
        │ str     ┆ i32    │
        ╞═════════╪════════╡
        │ 110     ┆ 6      │
        │ 101     ┆ 5      │
        │ 010     ┆ 2      │
        │ invalid ┆ null   │
        └─────────┴────────┘

        >>> df = pl.DataFrame({"hex": ["fa1e", "ff00", "cafe", None]})
        >>> df.with_columns(parsed=pl.col("hex").str.to_integer(base=16, strict=True))
        shape: (4, 2)
        ┌──────┬────────┐
        │ hex  ┆ parsed │
        │ ---  ┆ ---    │
        │ str  ┆ i64    │
        ╞══════╪════════╡
        │ fa1e ┆ 64030  │
        │ ff00 ┆ 65280  │
        │ cafe ┆ 51966  │
        │ null ┆ null   │
        └──────┴────────┘
        """
        base_pyexpr = parse_into_expression(base, str_as_lit=False)
        return wrap_expr(self._pyexpr.str_to_integer(base_pyexpr, dtype, strict))

    def contains_any(
        self, patterns: IntoExpr, *, ascii_case_insensitive: bool = False
    ) -> Expr:
        """
        Use the Aho-Corasick algorithm to find matches.

        Determines if any of the patterns are contained in the string.

        Parameters
        ----------
        patterns
            String patterns to search.
        ascii_case_insensitive
            Enable ASCII-aware case-insensitive matching.
            When this option is enabled, searching will be performed without respect
            to case for ASCII letters (a-z and A-Z) only.

        Notes
        -----
        This method supports matching on string literals only, and does not support
        regular expression matching.

        Examples
        --------
        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> df = pl.DataFrame(
        ...     {
        ...         "lyrics": [
        ...             "Everybody wants to rule the world",
        ...             "Tell me what you want, what you really really want",
        ...             "Can you feel the love tonight",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("lyrics").str.contains_any(["you", "me"]).alias("contains_any")
        ... )
        shape: (3, 2)
        ┌────────────────────────────────────────────────────┬──────────────┐
        │ lyrics                                             ┆ contains_any │
        │ ---                                                ┆ ---          │
        │ str                                                ┆ bool         │
        ╞════════════════════════════════════════════════════╪══════════════╡
        │ Everybody wants to rule the world                  ┆ false        │
        │ Tell me what you want, what you really really want ┆ true         │
        │ Can you feel the love tonight                      ┆ true         │
        └────────────────────────────────────────────────────┴──────────────┘
        """
        patterns_pyexpr = parse_into_expression(patterns, str_as_lit=False)
        return wrap_expr(
            self._pyexpr.str_contains_any(patterns_pyexpr, ascii_case_insensitive)
        )

    def replace_many(
        self,
        patterns: IntoExpr | Mapping[str, str],
        replace_with: IntoExpr | NoDefault = no_default,
        *,
        ascii_case_insensitive: bool = False,
    ) -> Expr:
        """
        Use the Aho-Corasick algorithm to replace many matches.

        Parameters
        ----------
        patterns
            String patterns to search and replace.
            Accepts expression input. Strings are parsed as column names, and other
            non-expression inputs are parsed as literals. Also accepts a mapping of
            patterns to their replacement as syntactic sugar for
            `replace_many(pl.Series(mapping.keys()), pl.Series(mapping.values()))`.
        replace_with
            Strings to replace where a pattern was a match.
            Accepts expression input. Non-expression inputs are parsed as literals.
            Length must match the length of `patterns` or have length 1. This can be
            broadcasted, so it supports many:one and many:many.
        ascii_case_insensitive
            Enable ASCII-aware case-insensitive matching.
            When this option is enabled, searching will be performed without respect
            to case for ASCII letters (a-z and A-Z) only.

        Notes
        -----
        This method supports matching on string literals only, and does not support
        regular expression matching.

        Examples
        --------
        Replace many patterns by passing sequences of equal length to the `patterns` and
        `replace_with` parameters.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> _ = pl.Config.set_tbl_width_chars(110)
        >>> df = pl.DataFrame(
        ...     {
        ...         "lyrics": [
        ...             "Everybody wants to rule the world",
        ...             "Tell me what you want, what you really really want",
        ...             "Can you feel the love tonight",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("lyrics")
        ...     .str.replace_many(
        ...         ["me", "you"],
        ...         ["you", "me"],
        ...     )
        ...     .alias("confusing")
        ... )
        shape: (3, 2)
        ┌────────────────────────────────────────────────────┬───────────────────────────────────────────────────┐
        │ lyrics                                             ┆ confusing                                         │
        │ ---                                                ┆ ---                                               │
        │ str                                                ┆ str                                               │
        ╞════════════════════════════════════════════════════╪═══════════════════════════════════════════════════╡
        │ Everybody wants to rule the world                  ┆ Everybody wants to rule the world                 │
        │ Tell me what you want, what you really really want ┆ Tell you what me want, what me really really want │
        │ Can you feel the love tonight                      ┆ Can me feel the love tonight                      │
        └────────────────────────────────────────────────────┴───────────────────────────────────────────────────┘

        Broadcast a replacement for many patterns by passing sequence of length 1 to the
        `replace_with` parameter.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> df = pl.DataFrame(
        ...     {
        ...         "lyrics": [
        ...             "Everybody wants to rule the world",
        ...             "Tell me what you want, what you really really want",
        ...             "Can you feel the love tonight",
        ...         ]
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.col("lyrics")
        ...     .str.replace_many(
        ...         ["me", "you", "they"],
        ...         [""],
        ...     )
        ...     .alias("removes_pronouns")
        ... )
        shape: (3, 2)
        ┌────────────────────────────────────────────────────┬────────────────────────────────────────────┐
        │ lyrics                                             ┆ removes_pronouns                           │
        │ ---                                                ┆ ---                                        │
        │ str                                                ┆ str                                        │
        ╞════════════════════════════════════════════════════╪════════════════════════════════════════════╡
        │ Everybody wants to rule the world                  ┆ Everybody wants to rule the world          │
        │ Tell me what you want, what you really really want ┆ Tell  what  want, what  really really want │
        │ Can you feel the love tonight                      ┆ Can  feel the love tonight                 │
        └────────────────────────────────────────────────────┴────────────────────────────────────────────┘

        Passing a mapping with patterns and replacements is also supported as syntactic
        sugar.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> _ = pl.Config.set_tbl_width_chars(110)
        >>> df = pl.DataFrame(
        ...     {
        ...         "lyrics": [
        ...             "Everybody wants to rule the world",
        ...             "Tell me what you want, what you really really want",
        ...             "Can you feel the love tonight",
        ...         ]
        ...     }
        ... )
        >>> mapping = {"me": "you", "you": "me", "want": "need"}
        >>> df.with_columns(
        ...     pl.col("lyrics").str.replace_many(mapping).alias("confusing")
        ... )
        shape: (3, 2)
        ┌────────────────────────────────────────────────────┬───────────────────────────────────────────────────┐
        │ lyrics                                             ┆ confusing                                         │
        │ ---                                                ┆ ---                                               │
        │ str                                                ┆ str                                               │
        ╞════════════════════════════════════════════════════╪═══════════════════════════════════════════════════╡
        │ Everybody wants to rule the world                  ┆ Everybody needs to rule the world                 │
        │ Tell me what you want, what you really really want ┆ Tell you what me need, what me really really need │
        │ Can you feel the love tonight                      ┆ Can me feel the love tonight                      │
        └────────────────────────────────────────────────────┴───────────────────────────────────────────────────┘
        """  # noqa: W505
        if replace_with is no_default:
            if not isinstance(patterns, Mapping):
                msg = "`replace_with` argument is required if `patterns` argument is not a Mapping type"
                raise TypeError(msg)
            # Early return in case of an empty mapping.
            if not patterns:
                return wrap_expr(self._pyexpr)
            replace_with = list(patterns.values())
            patterns = list(patterns.keys())

        patterns_pyexpr = parse_into_expression(
            patterns,  # type: ignore[arg-type]
            str_as_lit=False,
        )
        replace_with_pyexpr = parse_into_expression(replace_with, str_as_lit=True)
        return wrap_expr(
            self._pyexpr.str_replace_many(
                patterns_pyexpr, replace_with_pyexpr, ascii_case_insensitive
            )
        )

    @unstable()
    def extract_many(
        self,
        patterns: IntoExpr,
        *,
        ascii_case_insensitive: bool = False,
        overlapping: bool = False,
    ) -> Expr:
        """
        Use the Aho-Corasick algorithm to extract many matches.

        Parameters
        ----------
        patterns
            String patterns to search.
        ascii_case_insensitive
            Enable ASCII-aware case-insensitive matching.
            When this option is enabled, searching will be performed without respect
            to case for ASCII letters (a-z and A-Z) only.
        overlapping
            Whether matches may overlap.

        Notes
        -----
        This method supports matching on string literals only, and does not support
        regular expression matching.

        Examples
        --------
        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> df = pl.DataFrame({"values": ["discontent"]})
        >>> patterns = ["winter", "disco", "onte", "discontent"]
        >>> df.with_columns(
        ...     pl.col("values")
        ...     .str.extract_many(patterns, overlapping=False)
        ...     .alias("matches"),
        ...     pl.col("values")
        ...     .str.extract_many(patterns, overlapping=True)
        ...     .alias("matches_overlapping"),
        ... )
        shape: (1, 3)
        ┌────────────┬───────────┬─────────────────────────────────┐
        │ values     ┆ matches   ┆ matches_overlapping             │
        │ ---        ┆ ---       ┆ ---                             │
        │ str        ┆ list[str] ┆ list[str]                       │
        ╞════════════╪═══════════╪═════════════════════════════════╡
        │ discontent ┆ ["disco"] ┆ ["disco", "onte", "discontent"] │
        └────────────┴───────────┴─────────────────────────────────┘
        >>> df = pl.DataFrame(
        ...     {
        ...         "values": ["discontent", "rhapsody"],
        ...         "patterns": [
        ...             ["winter", "disco", "onte", "discontent"],
        ...             ["rhap", "ody", "coalesce"],
        ...         ],
        ...     }
        ... )
        >>> df.select(pl.col("values").str.extract_many("patterns"))
        shape: (2, 1)
        ┌─────────────────┐
        │ values          │
        │ ---             │
        │ list[str]       │
        ╞═════════════════╡
        │ ["disco"]       │
        │ ["rhap", "ody"] │
        └─────────────────┘
        """
        patterns_pyexpr = parse_into_expression(patterns, str_as_lit=False)
        return wrap_expr(
            self._pyexpr.str_extract_many(
                patterns_pyexpr, ascii_case_insensitive, overlapping
            )
        )

    @unstable()
    def find_many(
        self,
        patterns: IntoExpr,
        *,
        ascii_case_insensitive: bool = False,
        overlapping: bool = False,
    ) -> Expr:
        """
        Use the Aho-Corasick algorithm to find many matches.

        The function will return the bytes offset of the start of each match.
        The return type will be `List<UInt32>`

        Parameters
        ----------
        patterns
            String patterns to search.
        ascii_case_insensitive
            Enable ASCII-aware case-insensitive matching.
            When this option is enabled, searching will be performed without respect
            to case for ASCII letters (a-z and A-Z) only.
        overlapping
            Whether matches may overlap.

        Notes
        -----
        This method supports matching on string literals only, and does not support
        regular expression matching.

        Examples
        --------
        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> df = pl.DataFrame({"values": ["discontent"]})
        >>> patterns = ["winter", "disco", "onte", "discontent"]
        >>> df.with_columns(
        ...     pl.col("values")
        ...     .str.find_many(patterns, overlapping=False)
        ...     .alias("matches"),
        ...     pl.col("values")
        ...     .str.find_many(patterns, overlapping=True)
        ...     .alias("matches_overlapping"),
        ... )
        shape: (1, 3)
        ┌────────────┬───────────┬─────────────────────┐
        │ values     ┆ matches   ┆ matches_overlapping │
        │ ---        ┆ ---       ┆ ---                 │
        │ str        ┆ list[u32] ┆ list[u32]           │
        ╞════════════╪═══════════╪═════════════════════╡
        │ discontent ┆ [0]       ┆ [0, 4, 0]           │
        └────────────┴───────────┴─────────────────────┘
        >>> df = pl.DataFrame(
        ...     {
        ...         "values": ["discontent", "rhapsody"],
        ...         "patterns": [
        ...             ["winter", "disco", "onte", "discontent"],
        ...             ["rhap", "ody", "coalesce"],
        ...         ],
        ...     }
        ... )
        >>> df.select(pl.col("values").str.find_many("patterns"))
        shape: (2, 1)
        ┌───────────┐
        │ values    │
        │ ---       │
        │ list[u32] │
        ╞═══════════╡
        │ [0]       │
        │ [0, 5]    │
        └───────────┘
        """
        patterns_pyexpr = parse_into_expression(patterns, str_as_lit=False)
        return wrap_expr(
            self._pyexpr.str_find_many(
                patterns_pyexpr, ascii_case_insensitive, overlapping
            )
        )

    def join(self, delimiter: str = "", *, ignore_nulls: bool = True) -> Expr:
        """
        Vertically concatenate the string values in the column to a single string value.

        Parameters
        ----------
        delimiter
            The delimiter to insert between consecutive string values.
        ignore_nulls
            Ignore null values (default).
            If set to `False`, null values will be propagated. This means that
            if the column contains any null values, the output is null.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, None, 3]})
        >>> df.select(pl.col("foo").str.join("-"))
        shape: (1, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ str │
        ╞═════╡
        │ 1-3 │
        └─────┘
        >>> df.select(pl.col("foo").str.join(ignore_nulls=False))
        shape: (1, 1)
        ┌──────┐
        │ foo  │
        │ ---  │
        │ str  │
        ╞══════╡
        │ null │
        └──────┘
        """
        return wrap_expr(self._pyexpr.str_join(delimiter, ignore_nulls=ignore_nulls))

    @deprecated(
        "`str.concat` is deprecated; use `str.join` instead. Note also that the "
        "default `delimiter` for `str.join` is an empty string, not a hyphen."
    )
    def concat(
        self, delimiter: str | None = None, *, ignore_nulls: bool = True
    ) -> Expr:
        """
        Vertically concatenate the string values in the column to a single string value.

        .. deprecated:: 1.0.0
            Use :meth:`join` instead. Note that the default `delimiter` for :meth:`join`
            is an empty string instead of a hyphen.

        Parameters
        ----------
        delimiter
            The delimiter to insert between consecutive string values.
        ignore_nulls
            Ignore null values (default).
            If set to `False`, null values will be propagated. This means that
            if the column contains any null values, the output is null.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, None, 2]})
        >>> df.select(pl.col("foo").str.concat("-"))  # doctest: +SKIP
        shape: (1, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ str │
        ╞═════╡
        │ 1-2 │
        └─────┘
        >>> df.select(
        ...     pl.col("foo").str.concat("-", ignore_nulls=False)
        ... )  # doctest: +SKIP
        shape: (1, 1)
        ┌──────┐
        │ foo  │
        │ ---  │
        │ str  │
        ╞══════╡
        │ null │
        └──────┘
        """
        if delimiter is None:
            delimiter = "-"
        return self.join(delimiter, ignore_nulls=ignore_nulls)

    def escape_regex(self) -> Expr:
        r"""
        Returns string values with all regular expression meta characters escaped.

        Examples
        --------
        >>> df = pl.DataFrame({"text": ["abc", "def", None, "abc(\\w+)"]})
        >>> df.with_columns(pl.col("text").str.escape_regex().alias("escaped"))
         shape: (4, 2)
        ┌──────────┬──────────────┐
        │ text     ┆ escaped      │
        │ ---      ┆ ---          │
        │ str      ┆ str          │
        ╞══════════╪══════════════╡
        │ abc      ┆ abc          │
        │ def      ┆ def          │
        │ null     ┆ null         │
        │ abc(\w+) ┆ abc\(\\w\+\) │
        └──────────┴──────────────┘
        """
        return wrap_expr(self._pyexpr.str_escape_regex())

    def normalize(self, form: UnicodeForm = "NFC") -> Expr:
        """
        Returns the Unicode normal form of the string values.

        This uses the forms described in Unicode Standard Annex 15: <https://www.unicode.org/reports/tr15/>.

        Parameters
        ----------
        form : {'NFC', 'NFKC', 'NFD', 'NFKD'}
            Unicode form to use.

        Examples
        --------
        >>> df = pl.DataFrame({"text": ["01²", "ＫＡＤＯＫＡＷＡ"]})
        >>> new = df.with_columns(
        ...     nfc=pl.col("text").str.normalize("NFC"),
        ...     nfkc=pl.col("text").str.normalize("NFKC"),
        ... )
        >>> new
        shape: (2, 3)
        ┌──────────────────┬──────────────────┬──────────┐
        │ text             ┆ nfc              ┆ nfkc     │
        │ ---              ┆ ---              ┆ ---      │
        │ str              ┆ str              ┆ str      │
        ╞══════════════════╪══════════════════╪══════════╡
        │ 01²              ┆ 01²              ┆ 012      │
        │ ＫＡＤＯＫＡＷＡ    ┆ ＫＡＤＯＫＡＷＡ    ┆ KADOKAWA │
        └──────────────────┴──────────────────┴──────────┘
        >>> new.select(pl.all().str.len_bytes())
        shape: (2, 3)
        ┌──────┬─────┬──────┐
        │ text ┆ nfc ┆ nfkc │
        │ ---  ┆ --- ┆ ---  │
        │ u32  ┆ u32 ┆ u32  │
        ╞══════╪═════╪══════╡
        │ 4    ┆ 4   ┆ 3    │
        │ 24   ┆ 24  ┆ 8    │
        └──────┴─────┴──────┘
        """  # noqa: RUF002
        return wrap_expr(self._pyexpr.str_normalize(form))


def _validate_format_argument(format: str | None) -> None:
    if format is not None and ".%f" in format:
        message = (
            "Detected the pattern `.%f` in the chrono format string."
            " This pattern should not be used to parse values after a decimal point."
            " Use `%.f` instead."
            " See the full specification: https://docs.rs/chrono/latest/chrono/format/strftime"
        )
        warnings.warn(message, ChronoFormatWarning, stacklevel=find_stacklevel())
