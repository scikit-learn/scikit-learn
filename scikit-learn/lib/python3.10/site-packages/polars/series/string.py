from __future__ import annotations

from typing import TYPE_CHECKING

import polars._reexport as pl
import polars.functions as F
from polars._utils.deprecation import deprecate_nonkeyword_arguments, deprecated
from polars._utils.unstable import unstable
from polars._utils.various import no_default
from polars._utils.wrap import wrap_s
from polars.datatypes import Int64
from polars.datatypes.classes import Datetime
from polars.datatypes.constants import N_INFER_DEFAULT
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    import sys
    from collections.abc import Mapping

    from polars import Expr, Series
    from polars._plr import PySeries
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


@expr_dispatch
class StringNameSpace:
    """Series.str namespace."""

    _accessor = "str"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def to_date(
        self,
        format: str | None = None,
        *,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
    ) -> Series:
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

    def to_datetime(
        self,
        format: str | None = None,
        *,
        time_unit: TimeUnit | None = None,
        time_zone: str | None = None,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
        ambiguous: Ambiguous | pl.Series = "raise",
    ) -> pl.Series:
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
        if format is None and time_zone is None:
            if isinstance(ambiguous, str):
                ambiguous_s = pl.Series([ambiguous])
            else:
                ambiguous_s = ambiguous

            return wrap_s(
                self._s.str_to_datetime_infer(
                    time_unit,
                    strict,
                    exact,
                    ambiguous_s._s,
                )
            )
        else:
            ambiguous_expr = F.lit(ambiguous)
            s = wrap_s(self._s)
            return (
                s.to_frame()
                .select_seq(
                    F.col(s.name).str.to_datetime(
                        format,
                        time_unit=time_unit,
                        time_zone=time_zone,
                        strict=strict,
                        exact=exact,
                        cache=cache,
                        ambiguous=ambiguous_expr,
                    )
                )
                .to_series()
            )

    def to_time(
        self,
        format: str | None = None,
        *,
        strict: bool = True,
        cache: bool = True,
    ) -> Series:
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

    def strptime(
        self,
        dtype: PolarsTemporalType,
        format: str | None = None,
        *,
        strict: bool = True,
        exact: bool = True,
        cache: bool = True,
        ambiguous: Ambiguous | Series = "raise",
    ) -> Series:
        """
        Convert a String column into a Date/Datetime/Time column.

        Parameters
        ----------
        dtype
            The data type to convert to. Can be either Date, Datetime, or Time.
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
        if format is None and (
            dtype is Datetime
            or (isinstance(dtype, Datetime) and dtype.time_zone is None)
        ):
            time_unit = None
            if isinstance(dtype, Datetime):
                time_unit = dtype.time_unit

            return self.to_datetime(
                time_unit=time_unit,
                strict=strict,
                exact=exact,
                cache=cache,
                ambiguous=ambiguous,
            )
        else:
            ambiguous_expr = F.lit(ambiguous)
            s = wrap_s(self._s)
            return (
                s.to_frame()
                .select_seq(
                    F.col(s.name).str.strptime(
                        dtype,
                        format,
                        strict=strict,
                        exact=exact,
                        cache=cache,
                        ambiguous=ambiguous_expr,
                    )
                )
                .to_series()
            )

    @deprecate_nonkeyword_arguments(allowed_args=["self"], version="1.20.0")
    def to_decimal(
        self,
        inference_length: int = 100,
        *,
        scale: int | None = None,
    ) -> Series:
        """
        Convert a String column into a Decimal column.

        This method infers the needed parameters `precision` and `scale` if not
        given.

        .. versionchanged:: 1.20.0
            Parameter `inference_length` should now be passed as a keyword argument.

        Parameters
        ----------
        inference_length
            Number of elements to parse to determine the `precision` and `scale`
        scale
            Number of digits after the comma to use for the decimals.

        Examples
        --------
        >>> s = pl.Series(
        ...     ["40.12", "3420.13", "120134.19", "3212.98", "12.90", "143.09", "143.9"]
        ... )
        >>> s.str.to_decimal()
        shape: (7,)
        Series: '' [decimal[8,2]]
        [
            40.12
            3420.13
            120134.19
            3212.98
            12.90
            143.09
            143.90
        ]
        """
        if scale is not None:
            s = wrap_s(self._s)
            return (
                s.to_frame()
                .select_seq(F.col(s.name).str.to_decimal(scale=scale))
                .to_series()
            )
        else:
            return wrap_s(
                self._s.str_to_decimal_infer(inference_length=inference_length)
            )

    def len_bytes(self) -> Series:
        """
        Return the length of each string as the number of bytes.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

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
        >>> s = pl.Series(["Café", "345", "東京", None])
        >>> s.str.len_bytes()
        shape: (4,)
        Series: '' [u32]
        [
            5
            3
            6
            null
        ]
        """

    def len_chars(self) -> Series:
        """
        Return the length of each string as the number of characters.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

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
        >>> s = pl.Series(["Café", "345", "東京", None])
        >>> s.str.len_chars()
        shape: (4,)
        Series: '' [u32]
        [
            4
            3
            2
            null
        ]
        """

    def contains(
        self, pattern: str | Expr, *, literal: bool = False, strict: bool = True
    ) -> Series:
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

        Default (case-sensitive) match:

        >>> s = pl.Series("s", ["AAA", "aAa", "aaa"])
        >>> s.str.contains("AA").to_list()
        [True, False, False]

        Case-insensitive match, using an inline flag:

        >>> s = pl.Series("s", ["AAA", "aAa", "aaa"])
        >>> s.str.contains("(?i)AA").to_list()
        [True, True, True]

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series(["Crab", "cat and dog", "rab$bit", None])
        >>> s.str.contains("cat|bit")
        shape: (4,)
        Series: '' [bool]
        [
            false
            true
            true
            null
        ]
        >>> s.str.contains("rab$", literal=True)
        shape: (4,)
        Series: '' [bool]
        [
            false
            false
            true
            null
        ]
        """

    def find(
        self, pattern: str | Expr, *, literal: bool = False, strict: bool = True
    ) -> Series:
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

        >>> s = pl.Series("s", ["AAA", "aAa", "aaa"])

        Default (case-sensitive) match:

        >>> s.str.find("Aa").to_list()
        [None, 1, None]

        Case-insensitive match, using an inline flag:

        >>> s.str.find("(?i)Aa").to_list()
        [0, 0, 0]

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        See Also
        --------
        contains : Check if the string contains a substring that matches a pattern.

        Examples
        --------
        >>> s = pl.Series("txt", ["Crab", "Lobster", None, "Crustacean"])

        Find the index of the first substring matching a regex pattern:

        >>> s.str.find("a|e").rename("idx_rx")
        shape: (4,)
        Series: 'idx_rx' [u32]
        [
            2
            5
            null
            5
        ]

        Find the index of the first substring matching a literal pattern:

        >>> s.str.find("e", literal=True).rename("idx_lit")
        shape: (4,)
        Series: 'idx_lit' [u32]
        [
            null
            5
            null
            7
        ]

        Match against a pattern found in another column or (expression):

        >>> p = pl.Series("pat", ["a[bc]", "b.t", "[aeiuo]", "(?i)A[BC]"])
        >>> s.str.find(p).rename("idx")
        shape: (4,)
        Series: 'idx' [u32]
        [
            2
            2
            null
            5
        ]
        """

    def ends_with(self, suffix: str | Expr | None) -> Series:
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
        >>> s = pl.Series("fruits", ["apple", "mango", None])
        >>> s.str.ends_with("go")
        shape: (3,)
        Series: 'fruits' [bool]
        [
            false
            true
            null
        ]
        """

    def starts_with(self, prefix: str | Expr) -> Series:
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
        >>> s = pl.Series("fruits", ["apple", "mango", None])
        >>> s.str.starts_with("app")
        shape: (3,)
        Series: 'fruits' [bool]
        [
            true
            false
            null
        ]
        """

    def decode(self, encoding: TransferEncoding, *, strict: bool = True) -> Series:
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
        Series
            Series of data type :class:`Binary`.

        Examples
        --------
        >>> s = pl.Series("color", ["000000", "ffff00", "0000ff"])
        >>> s.str.decode("hex")
        shape: (3,)
        Series: 'color' [binary]
        [
                b"\x00\x00\x00"
                b"\xff\xff\x00"
                b"\x00\x00\xff"
        ]
        """

    def encode(self, encoding: TransferEncoding) -> Series:
        """
        Encode a value using the provided encoding.

        Parameters
        ----------
        encoding : {'hex', 'base64'}
            The encoding to use.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> s = pl.Series(["foo", "bar", None])
        >>> s.str.encode("hex")
        shape: (3,)
        Series: '' [str]
        [
            "666f6f"
            "626172"
            null
        ]
        """

    def json_decode(
        self,
        dtype: PolarsDataType | None = None,
        *,
        infer_schema_length: int | None = N_INFER_DEFAULT,
    ) -> Series:
        """
        Parse string values as JSON.

        Throws an error if invalid JSON strings are encountered.

        Parameters
        ----------
        dtype
            The dtype to cast the extracted value to. If None, the dtype will be
            inferred from the JSON value.
        infer_schema_length
            The maximum number of rows to scan for schema inference.
            If set to `None`, the full data may be scanned *(this is slow)*.

        See Also
        --------
        json_path_match : Extract the first match of json string with provided JSONPath
            expression.

        Examples
        --------
        >>> s = pl.Series("json", ['{"a":1, "b": true}', None, '{"a":2, "b": false}'])
        >>> s.str.json_decode()
        shape: (3,)
        Series: 'json' [struct[2]]
        [
                {1,true}
                null
                {2,false}
        ]
        """
        if dtype is not None:
            s = wrap_s(self._s)
            return (
                s.to_frame()
                .select_seq(F.col(s.name).str.json_decode(dtype))
                .to_series()
            )

        return wrap_s(self._s.str_json_decode(infer_schema_length))

    def json_path_match(self, json_path: IntoExprColumn) -> Series:
        """
        Extract the first match of JSON string with provided JSONPath expression.

        Throw errors if encounter invalid JSON strings.
        All return values will be cast to String regardless of the original value.

        Documentation on JSONPath standard can be found
        `here <https://goessner.net/articles/JsonPath/>`_.

        Parameters
        ----------
        json_path
            A valid JSON path query string.

        Returns
        -------
        Series
            Series of data type :class:`String`. Contains null values if the original
            value is null or the json_path returns nothing.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"json_val": ['{"a":"1"}', None, '{"a":2}', '{"a":2.1}', '{"a":true}']}
        ... )
        >>> df.select(pl.col("json_val").str.json_path_match("$.a"))[:, 0]
        shape: (5,)
        Series: 'json_val' [str]
        [
            "1"
            null
            "2"
            "2.1"
            "true"
        ]
        """

    def extract(self, pattern: IntoExprColumn, group_index: int = 1) -> Series:
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

        Returns
        -------
        Series
            Series of data type :class:`String`. Contains null values if the original
            value is null or regex captures nothing.

        Notes
        -----
        To modify regular expression behaviour (such as multi-line matching)
        with flags, use the inline `(?iLmsuxU)` syntax. For example:

        >>> s = pl.Series(
        ...     name="lines",
        ...     values=[
        ...         "I Like\nThose\nOdds",
        ...         "This is\nThe Way",
        ...     ],
        ... )
        >>> s.str.extract(r"(?m)^(T\w+)", 1).alias("matches")
        shape: (2,)
        Series: 'matches' [str]
        [
            "Those"
            "This"
        ]

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        Examples
        --------
        >>> s = pl.Series(
        ...     name="url",
        ...     values=[
        ...         "http://vote.com/ballon_dor?ref=polars&candidate=messi",
        ...         "http://vote.com/ballon_dor?candidate=ronaldo&ref=polars",
        ...         "http://vote.com/ballon_dor?error=404&ref=unknown",
        ...     ],
        ... )
        >>> s.str.extract(r"candidate=(\w+)", 1).alias("candidate")
        shape: (3,)
        Series: 'candidate' [str]
        [
            "messi"
            "ronaldo"
            null
        ]
        """

    def extract_all(self, pattern: str | Series) -> Series:
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

        >>> s = pl.Series(
        ...     name="email",
        ...     values=[
        ...         "real.email@spam.com",
        ...         "some_account@somewhere.net",
        ...         "abc.def.ghi.jkl@uvw.xyz.co.uk",
        ...     ],
        ... )
        >>> # extract name/domain parts from email, using verbose regex
        >>> s.str.extract_all(
        ...     r"""(?xi)   # activate 'verbose' and 'case-insensitive' flags
        ...       [         # (start character group)
        ...         A-Z     # letters
        ...         0-9     # digits
        ...         ._%+\-  # special chars
        ...       ]         # (end character group)
        ...       +         # 'one or more' quantifier
        ...     """
        ... ).alias("email_parts")
        shape: (3,)
        Series: 'email_parts' [list[str]]
        [
            ["real.email", "spam.com"]
            ["some_account", "somewhere.net"]
            ["abc.def.ghi.jkl", "uvw.xyz.co.uk"]
        ]

        See the regex crate's section on `grouping and flags
        <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_ for
        additional information about the use of inline expression modifiers.

        Returns
        -------
        Series
            Series of data type `List(String)`.

        Examples
        --------
        >>> s = pl.Series("foo", ["123 bla 45 asd", "xyz 678 910t", "bar", None])
        >>> s.str.extract_all(r"\d+")
        shape: (4,)
        Series: 'foo' [list[str]]
        [
            ["123", "45"]
            ["678", "910"]
            []
            null
        ]

        '''

    def extract_groups(self, pattern: str) -> Series:
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

        For example, we can access the first group via the string `"1"`::

            >>> (
            ...     pl.Series(["foo bar baz"])
            ...     .str.extract_groups(r"(\w+) (.+) (\w+)")
            ...     .struct["1"]
            ... )
            shape: (1,)
            Series: '1' [str]
            [
                "foo"
            ]

        Returns
        -------
        Series
            Series of data type :class:`Struct` with fields of data type
            :class:`String`.

        Examples
        --------
        >>> s = pl.Series(
        ...     name="url",
        ...     values=[
        ...         "http://vote.com/ballon_dor?candidate=messi&ref=python",
        ...         "http://vote.com/ballon_dor?candidate=weghorst&ref=polars",
        ...         "http://vote.com/ballon_dor?error=404&ref=rust",
        ...     ],
        ... )
        >>> s.str.extract_groups(r"candidate=(?<candidate>\w+)&ref=(?<ref>\w+)")
        shape: (3,)
        Series: 'url' [struct[2]]
        [
            {"messi","python"}
            {"weghorst","polars"}
            {null,null}
        ]
        """

    def count_matches(self, pattern: str | Series, *, literal: bool = False) -> Series:
        r"""
        Count all successive non-overlapping regex matches.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_. Can also be a :class:`Series` of
            regular expressions.
        literal
            Treat `pattern` as a literal string, not as a regular expression.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`. Returns null if the original
            value is null.

        Examples
        --------
        >>> s = pl.Series("foo", ["123 bla 45 asd", "xyz 678 910t", "bar", None])
        >>> # count digits
        >>> s.str.count_matches(r"\d")
        shape: (4,)
        Series: 'foo' [u32]
        [
            5
            6
            0
            null
        ]

        >>> s = pl.Series("bar", ["12 dbc 3xy", "cat\\w", "1zy3\\d\\d", None])
        >>> s.str.count_matches(r"\d", literal=True)
        shape: (4,)
        Series: 'bar' [u32]
        [
            0
            0
            2
            null
        ]
        """

    def split(self, by: IntoExpr, *, inclusive: bool = False) -> Series:
        """
        Split the string by a substring.

        Parameters
        ----------
        by
            Substring to split by.
        inclusive
            If True, include the split character/string in the results.

        Returns
        -------
        Series
            Series of data type `List(String)`.
        """

    def split_exact(self, by: IntoExpr, n: int, *, inclusive: bool = False) -> Series:
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

        Examples
        --------
        >>> df = pl.DataFrame({"x": ["a_1", None, "c", "d_4"]})
        >>> df["x"].str.split_exact("_", 1).alias("fields")
        shape: (4,)
        Series: 'fields' [struct[2]]
        [
                {"a","1"}
                {null,null}
                {"c",null}
                {"d","4"}
        ]

        Split string values in column x in exactly 2 parts and assign
        each part to a new column.

        >>> (
        ...     df["x"]
        ...     .str.split_exact("_", 1)
        ...     .struct.rename_fields(["first_part", "second_part"])
        ...     .alias("fields")
        ...     .to_frame()
        ...     .unnest("fields")
        ... )
        shape: (4, 2)
        ┌────────────┬─────────────┐
        │ first_part ┆ second_part │
        │ ---        ┆ ---         │
        │ str        ┆ str         │
        ╞════════════╪═════════════╡
        │ a          ┆ 1           │
        │ null       ┆ null        │
        │ c          ┆ null        │
        │ d          ┆ 4           │
        └────────────┴─────────────┘

        Returns
        -------
        Series
            Series of data type :class:`Struct` with fields of data type
            :class:`String`.
        """

    def splitn(self, by: IntoExpr, n: int) -> Series:
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

        Examples
        --------
        >>> df = pl.DataFrame({"s": ["foo bar", None, "foo-bar", "foo bar baz"]})
        >>> df["s"].str.splitn(" ", 2).alias("fields")
        shape: (4,)
        Series: 'fields' [struct[2]]
        [
                {"foo","bar"}
                {null,null}
                {"foo-bar",null}
                {"foo","bar baz"}
        ]

        Split string values in column s in exactly 2 parts and assign
        each part to a new column.

        >>> (
        ...     df["s"]
        ...     .str.splitn(" ", 2)
        ...     .struct.rename_fields(["first_part", "second_part"])
        ...     .alias("fields")
        ...     .to_frame()
        ...     .unnest("fields")
        ... )
        shape: (4, 2)
        ┌────────────┬─────────────┐
        │ first_part ┆ second_part │
        │ ---        ┆ ---         │
        │ str        ┆ str         │
        ╞════════════╪═════════════╡
        │ foo        ┆ bar         │
        │ null       ┆ null        │
        │ foo-bar    ┆ null        │
        │ foo        ┆ bar baz     │
        └────────────┴─────────────┘

        Returns
        -------
        Series
            Series of data type :class:`Struct` with fields of data type
            :class:`String`.
        """

    def replace(
        self, pattern: str, value: str, *, literal: bool = False, n: int = 1
    ) -> Series:
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
          use the inline `(?iLmsuxU)` syntax. (See the regex crate's section on
          `grouping and flags <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_
          for additional information about the use of inline expression modifiers).

        * The dollar sign (`$`) is a special character related to capture groups; if you
          want to replace some target pattern with characters that include a literal `$`
          you should escape it by doubling it up as `$$`, or set `literal=True` if you
          do not need a full regular expression pattern match. Otherwise, you will be
          referencing a (potentially non-existent) capture group.

          If not escaped, the `$0` in the replacement value (below) represents a capture
          group:

          .. code-block:: python

              >>> s = pl.Series("cents", ["000.25", "00.50", "0.75"])
              >>> s.str.replace(r"^(0+)\.", "$0.")
              shape: (3,)
              Series: 'cents' [str]
              [
                "000..25"
                "00..50"
                "0..75"
              ]

          To have `$` represent a literal value, it should be doubled-up as `$$`
          (or, for simpler find/replace operations, set `literal=True` if you do
          not require a full regular expression match):

          .. code-block:: python

              >>> s.str.replace(r"^(0+)\.", "$$0.")
              shape: (3,)
              Series: 'cents' [str]
              [
                "$0.25"
                "$0.50"
                "$0.75"
              ]

        Examples
        --------
        >>> s = pl.Series(["123abc", "abc456"])
        >>> s.str.replace(r"abc\b", "ABC")
        shape: (2,)
        Series: '' [str]
        [
            "123ABC"
            "abc456"
        ]

        Capture groups are supported. Use `$1` or `${1}` in the `value` string to refer
        to the first capture group in the `pattern`, `$2` or `${2}` to refer to the
        second capture group, and so on. You can also use *named* capture groups.

        >>> s = pl.Series(["hat", "hut"])
        >>> s.str.replace("h(.)t", "b${1}d")
        shape: (2,)
        Series: '' [str]
        [
            "bad"
            "bud"
        ]
        >>> s.str.replace("h(?<vowel>.)t", "b${vowel}d")
        shape: (2,)
        Series: '' [str]
        [
            "bad"
            "bud"
        ]

        Apply case-insensitive string replacement using the `(?i)` flag.

        >>> s = pl.Series("weather", ["Foggy", "Rainy", "Sunny"])
        >>> s.str.replace(r"(?i)foggy|rainy", "Sunny")
        shape: (3,)
        Series: 'weather' [str]
        [
            "Sunny"
            "Sunny"
            "Sunny"
        ]
        """

    def replace_all(self, pattern: str, value: str, *, literal: bool = False) -> Series:
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
          use the inline `(?iLmsuxU)` syntax. (See the regex crate's section on
          `grouping and flags <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_
          for additional information about the use of inline expression modifiers).

        * The dollar sign (`$`) is a special character related to capture groups; if you
          want to replace some target pattern with characters that include a literal `$`
          you should escape it by doubling it up as `$$`, or set `literal=True` if you
          do not need a full regular expression pattern match. Otherwise, you will be
          referencing a (potentially non-existent) capture group.

          In the example below we need to double up `$` (to represent a literal dollar
          sign, and then refer to the capture group using `$n` or `${n}`, hence the
          three consecutive `$` characters in the replacement value:

          .. code-block:: python

              >>> s = pl.Series("cost", ["#12.34", "#56.78"])
              >>> s.str.replace_all(r"#(\d+)", "$$${1}").alias("cost_usd")
              shape: (2,)
              Series: 'cost_usd' [str]
              [
                  "$12.34"
                  "$56.78"
              ]

        Examples
        --------
        >>> s = pl.Series(["123abc", "abc456"])
        >>> s.str.replace_all(r"abc\b", "ABC")
        shape: (2,)
        Series: '' [str]
        [
            "123ABC"
            "abc456"
        ]

        Capture groups are supported. Use `$1` or `${1}` in the `value` string to refer
        to the first capture group in the `pattern`, `$2` or `${2}` to refer to the
        second capture group, and so on. You can also use *named* capture groups.

        >>> s = pl.Series(["hat", "hut"])
        >>> s.str.replace_all("h(.)t", "b${1}d")
        shape: (2,)
        Series: '' [str]
        [
            "bad"
            "bud"
        ]
        >>> s.str.replace_all("h(?<vowel>.)t", "b${vowel}d")
        shape: (2,)
        Series: '' [str]
        [
            "bad"
            "bud"
        ]

        Apply case-insensitive string replacement using the `(?i)` flag.

        >>> s = pl.Series("weather", ["Foggy", "Rainy", "Sunny"])
        >>> s.str.replace_all(r"(?i)foggy|rainy", "Sunny")
        shape: (3,)
        Series: 'weather' [str]
        [
            "Sunny"
            "Sunny"
            "Sunny"
        ]
        """

    def strip_chars(self, characters: IntoExpr = None) -> Series:
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
        >>> s = pl.Series([" hello ", "\tworld"])
        >>> s.str.strip_chars()
        shape: (2,)
        Series: '' [str]
        [
                "hello"
                "world"
        ]

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so, unless that whitespace is
        also included in the string.

        >>> s.str.strip_chars("o ")
        shape: (2,)
        Series: '' [str]
        [
            "hell"
            "	world"
        ]
        """

    def strip_chars_start(self, characters: IntoExpr = None) -> Series:
        r"""
        Remove leading characters.

        Parameters
        ----------
        characters
            The set of characters to be removed. All combinations of this set of
            characters will be stripped from the start of the string. If set to None
            (default), all leading whitespace is removed instead.

        Examples
        --------
        >>> s = pl.Series([" hello ", "\tworld"])
        >>> s.str.strip_chars_start()
        shape: (2,)
        Series: '' [str]
        [
                "hello "
                "world"
        ]

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so.

        >>> s.str.strip_chars_start("wod\t")
        shape: (2,)
        Series: '' [str]
        [
                " hello "
                "rld"
        ]
        """

    def strip_chars_end(self, characters: IntoExpr = None) -> Series:
        r"""
        Remove trailing characters.

        Parameters
        ----------
        characters
            The set of characters to be removed. All combinations of this set of
            characters will be stripped from the end of the string. If set to None
            (default), all trailing whitespace is removed instead.

        Examples
        --------
        >>> s = pl.Series([" hello ", "world\t"])
        >>> s.str.strip_chars_end()
        shape: (2,)
        Series: '' [str]
        [
                " hello"
                "world"
        ]

        Characters can be stripped by passing a string as argument. Note that whitespace
        will not be stripped automatically when doing so.

        >>> s.str.strip_chars_end("orld\t")
        shape: (2,)
        Series: '' [str]
        [
            " hello "
            "w"
        ]
        """

    def strip_prefix(self, prefix: IntoExpr) -> Series:
        """
        Remove prefix.

        The prefix will be removed from the string exactly once, if found.

        Parameters
        ----------
        prefix
            The prefix to be removed.

        Examples
        --------
        >>> s = pl.Series(["foobar", "foofoobar", "foo", "bar"])
        >>> s.str.strip_prefix("foo")
        shape: (4,)
        Series: '' [str]
        [
                "bar"
                "foobar"
                ""
                "bar"
        ]
        """

    def strip_suffix(self, suffix: IntoExpr) -> Series:
        """
        Remove suffix.

        The suffix will be removed from the string exactly once, if found.

        Parameters
        ----------
        suffix
            The suffix to be removed.

        Examples
        --------
        >>> s = pl.Series(["foobar", "foobarbar", "foo", "bar"])
        >>> s.str.strip_suffix("bar")
        shape: (4,)
        Series: '' [str]
        [
                "foo"
                "foobar"
                "foo"
                ""
        ]
        """

    def pad_start(self, length: int | IntoExprColumn, fill_char: str = " ") -> Series:
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
        >>> s = pl.Series("a", ["cow", "monkey", "hippopotamus", None])
        >>> s.str.pad_start(8, "*")
        shape: (4,)
        Series: 'a' [str]
        [
            "*****cow"
            "**monkey"
            "hippopotamus"
            null
        ]
        """

    def pad_end(self, length: int | IntoExprColumn, fill_char: str = " ") -> Series:
        """
        Pad the end of the string until it reaches the given length.

        Parameters
        ----------
        length
            Pad the string until it reaches this length. Strings with length equal to or
            greater than this value are returned as-is.
        fill_char
            The character to pad the string with.

        See Also
        --------
        pad_start

        Examples
        --------
        >>> s = pl.Series(["cow", "monkey", "hippopotamus", None])
        >>> s.str.pad_end(8, "*")
        shape: (4,)
        Series: '' [str]
        [
            "cow*****"
            "monkey**"
            "hippopotamus"
            null
        ]
        """

    def zfill(self, length: int | IntoExprColumn) -> Series:
        """
        Pad the start of the string with zeros until it reaches the given length.

        A sign prefix (`-`) is handled by inserting the padding after the sign character
        rather than before.

        Parameters
        ----------
        length
            Pad the string until it reaches this length. Strings with length equal to or
            greater than this value are returned as-is.

        See Also
        --------
        pad_start

        Notes
        -----
        This method is intended for padding numeric strings. If your data contains
        non-ASCII characters, use :func:`pad_start` instead.

        Examples
        --------
        >>> s = pl.Series([-1, 123, 999999, None])
        >>> s.cast(pl.String).str.zfill(4)
        shape: (4,)
        Series: '' [str]
        [
                "-001"
                "0123"
                "999999"
                null
        ]
        """

    def to_lowercase(self) -> Series:
        """
        Modify strings to their lowercase equivalent.

        Examples
        --------
        >>> s = pl.Series("foo", ["CAT", "DOG"])
        >>> s.str.to_lowercase()
        shape: (2,)
        Series: 'foo' [str]
        [
            "cat"
            "dog"
        ]
        """

    def to_uppercase(self) -> Series:
        """
        Modify strings to their uppercase equivalent.

        Examples
        --------
        >>> s = pl.Series("foo", ["cat", "dog"])
        >>> s.str.to_uppercase()
        shape: (2,)
        Series: 'foo' [str]
        [
            "CAT"
            "DOG"
        ]
        """

    def to_titlecase(self) -> Series:
        """
        Modify strings to their titlecase equivalent.

        Notes
        -----
        This is a form of case transform where the first letter of each word is
        capitalized, with the rest of the word in lowercase. Non-alphanumeric
        characters define the word boundaries.

        Examples
        --------
        >>> s = pl.Series(
        ...     "quotes",
        ...     [
        ...         "'e.t. phone home'",
        ...         "you talkin' to me?",
        ...         "to infinity,and BEYOND!",
        ...     ],
        ... )
        >>> s.str.to_titlecase()
        shape: (3,)
        Series: 'quotes' [str]
        [
            "'E.T. Phone Home'"
            "You Talkin' To Me?"
            "To Infinity,And Beyond!"
        ]
        """

    def reverse(self) -> Series:
        """
        Returns string values in reversed order.

        Examples
        --------
        >>> s = pl.Series("text", ["foo", "bar", "man\u0303ana"])
        >>> s.str.reverse()
        shape: (3,)
        Series: 'text' [str]
        [
            "oof"
            "rab"
            "anañam"
        ]
        """

    def slice(
        self, offset: int | IntoExprColumn, length: int | IntoExprColumn | None = None
    ) -> Series:
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
        Series
            Series of data type :class:`String`.

        Notes
        -----
        Both the `offset` and `length` inputs are defined in terms of the number
        of characters in the (UTF8) string. A character is defined as a
        `Unicode scalar value`_. A single character is represented by a single byte
        when working with ASCII text, and a maximum of 4 bytes otherwise.

        .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        Examples
        --------
        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"])
        >>> s.str.slice(-3)
        shape: (4,)
        Series: '' [str]
        [
            "ear"
            null
            "aya"
            "uit"
        ]

        Using the optional `length` parameter

        >>> s.str.slice(4, length=3)
        shape: (4,)
        Series: '' [str]
        [
            ""
            null
            "ya"
            "onf"
        ]
        """

    def head(self, n: int | IntoExprColumn) -> Series:
        """
        Return the first n characters of each string in a String Series.

        Parameters
        ----------
        n
            Length of the slice (integer or expression). Negative indexing is supported;
            see note (2) below.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Notes
        -----
        1) The `n` input is defined in terms of the number of characters in the (UTF8)
           string. A character is defined as a `Unicode scalar value`_. A single
           character is represented by a single byte when working with ASCII text, and a
           maximum of 4 bytes otherwise.

           .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        2) When `n` is negative, `head` returns characters up to the `n`th from the end
           of the string. For example, if `n = -3`, then all characters except the last
           three are returned.

        3) If the length of the string has fewer than `n` characters, the full string is
           returned.

        Examples
        --------
        Return up to the first 5 characters.

        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"])
        >>> s.str.head(5)
        shape: (4,)
        Series: '' [str]
        [
            "pear"
            null
            "papay"
            "drago"
        ]

        Return up to the 3rd character from the end.

        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"])
        >>> s.str.head(-3)
        shape: (4,)
        Series: '' [str]
        [
            "p"
            null
            "pap"
            "dragonfr"
        ]
        """

    def tail(self, n: int | IntoExprColumn) -> Series:
        """
        Return the last n characters of each string in a String Series.

        Parameters
        ----------
        n
            Length of the slice (integer or expression). Negative indexing is supported;
            see note (2) below.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Notes
        -----
        1) The `n` input is defined in terms of the number of characters in the (UTF8)
           string. A character is defined as a `Unicode scalar value`_. A single
           character is represented by a single byte when working with ASCII text, and a
           maximum of 4 bytes otherwise.

           .. _Unicode scalar value: https://www.unicode.org/glossary/#unicode_scalar_value

        2) When `n` is negative, `tail` returns characters starting from the `n`th from
           the beginning of the string. For example, if `n = -3`, then all characters
           except the first three are returned.

        3) If the length of the string has fewer than `n` characters, the full string is
           returned.

        Examples
        --------
        Return up to the last 5 characters:

        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"])
        >>> s.str.tail(5)
        shape: (4,)
        Series: '' [str]
        [
            "pear"
            null
            "apaya"
            "fruit"
        ]

        Return from the 3rd character to the end:

        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"])
        >>> s.str.tail(-3)
        shape: (4,)
        Series: '' [str]
        [
            "r"
            null
            "aya"
            "gonfruit"
        ]
        """

    @deprecated(
        '`Series.str.explode` is deprecated; use `Series.str.split("").explode()` instead. '
        "Note that empty strings will result in null instead of being preserved. To get "
        "the exact same behavior, split first and then use a `pl.when...then...otherwise` "
        "expression to handle the empty list before exploding. "
    )
    def explode(self) -> Series:
        """
        Returns a column with a separate row for every string character.

        .. deprecated:: 0.20.31
            Use the `.str.split("").explode()` method instead. Note that empty strings
            will result in null instead of being preserved. To get the exact same
            behavior, split first and then use a `pl.when...then...otherwise`
            expression to handle the empty list before exploding.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> s = pl.Series("a", ["foo", "bar"])
        >>> s.str.explode()  # doctest: +SKIP
        shape: (6,)
        Series: 'a' [str]
        [
                "f"
                "o"
                "o"
                "b"
                "a"
                "r"
        ]
        """

    def to_integer(
        self,
        *,
        base: int | IntoExprColumn = 10,
        dtype: PolarsIntegerType = Int64,
        strict: bool = True,
    ) -> Series:
        """
        Convert an String column into a column of dtype with base radix.

        Parameters
        ----------
        base
            Positive integer or expression which is the base of the string
            we are parsing.
            Default: 10.
        dtype
            Polars integer type to cast to.
            Default: :class:`Int64`.
        strict
            Bool, Default=True will raise any ParseError or overflow as ComputeError.
            False silently convert to Null.

        Returns
        -------
        Series
            Series of data.

        Examples
        --------
        >>> s = pl.Series("bin", ["110", "101", "010", "invalid"])
        >>> s.str.to_integer(base=2, dtype=pl.Int32, strict=False)
        shape: (4,)
        Series: 'bin' [i32]
        [
                6
                5
                2
                null
        ]

        >>> s = pl.Series("hex", ["fa1e", "ff00", "cafe", None])
        >>> s.str.to_integer(base=16)
        shape: (4,)
        Series: 'hex' [i64]
        [
                64030
                65280
                51966
                null
        ]
        """

    def contains_any(
        self, patterns: Series | list[str], *, ascii_case_insensitive: bool = False
    ) -> Series:
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
        >>> s = pl.Series(
        ...     "lyrics",
        ...     [
        ...         "Everybody wants to rule the world",
        ...         "Tell me what you want, what you really really want",
        ...         "Can you feel the love tonight",
        ...     ],
        ... )
        >>> s.str.contains_any(["you", "me"])
        shape: (3,)
        Series: 'lyrics' [bool]
        [
            false
            true
            true
        ]
        """

    def replace_many(
        self,
        patterns: Series | list[str] | Mapping[str, str],
        replace_with: Series | list[str] | str | NoDefault = no_default,
        *,
        ascii_case_insensitive: bool = False,
    ) -> Series:
        """
        Use the Aho-Corasick algorithm to replace many matches.

        Parameters
        ----------
        patterns
            String patterns to search and replace.
            Also accepts a mapping of patterns to their replacement as syntactic sugar
            for `replace_many(pl.Series(mapping.keys()), pl.Series(mapping.values()))`.
        replace_with
            Strings to replace where a pattern was a match.
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
        Replace many patterns by passing lists of equal length to the `patterns` and
        `replace_with` parameters.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> s = pl.Series(
        ...     "lyrics",
        ...     [
        ...         "Everybody wants to rule the world",
        ...         "Tell me what you want, what you really really want",
        ...         "Can you feel the love tonight",
        ...     ],
        ... )
        >>> s.str.replace_many(["you", "me"], ["me", "you"])
        shape: (3,)
        Series: 'lyrics' [str]
        [
            "Everybody wants to rule the world"
            "Tell you what me want, what me really really want"
            "Can me feel the love tonight"
        ]

        Broadcast a replacement for many patterns by passing a sequence of length 1 to
        the `replace_with` parameter.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> s = pl.Series(
        ...     "lyrics",
        ...     [
        ...         "Everybody wants to rule the world",
        ...         "Tell me what you want, what you really really want",
        ...         "Can you feel the love tonight",
        ...     ],
        ... )
        >>> s.str.replace_many(["me", "you", "they"], [""])
        shape: (3,)
        Series: 'lyrics' [str]
        [
            "Everybody wants to rule the world"
            "Tell  what  want, what  really really want"
            "Can  feel the love tonight"
        ]

        Passing a mapping with patterns and replacements is also supported as syntactic
        sugar.

        >>> _ = pl.Config.set_fmt_str_lengths(100)
        >>> s = pl.Series(
        ...     "lyrics",
        ...     [
        ...         "Everybody wants to rule the world",
        ...         "Tell me what you want, what you really really want",
        ...         "Can you feel the love tonight",
        ...     ],
        ... )
        >>> mapping = {"me": "you", "you": "me", "want": "need"}
        >>> s.str.replace_many(mapping)
        shape: (3,)
        Series: 'lyrics' [str]
        [
            "Everybody needs to rule the world"
            "Tell you what me need, what me really really need"
            "Can me feel the love tonight"
        ]
        """

    @unstable()
    def extract_many(
        self,
        patterns: Series | list[str],
        *,
        ascii_case_insensitive: bool = False,
        overlapping: bool = False,
    ) -> Series:
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
        >>> s = pl.Series("values", ["discontent"])
        >>> patterns = ["winter", "disco", "onte", "discontent"]
        >>> s.str.extract_many(patterns, overlapping=True)
        shape: (1,)
        Series: 'values' [list[str]]
        [
            ["disco", "onte", "discontent"]
        ]

        """

    @unstable()
    def find_many(
        self,
        patterns: IntoExpr,
        *,
        ascii_case_insensitive: bool = False,
        overlapping: bool = False,
    ) -> Series:
        """
        Use the Aho-Corasick algorithm to find all matches.

        The function returns the byte offset of the start of each match.
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

    def join(self, delimiter: str = "", *, ignore_nulls: bool = True) -> Series:
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
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> s = pl.Series([1, None, 3])
        >>> s.str.join("-")
        shape: (1,)
        Series: '' [str]
        [
            "1-3"
        ]
        >>> s.str.join(ignore_nulls=False)
        shape: (1,)
        Series: '' [str]
        [
            null
        ]
        """

    @deprecated(
        "`Series.str.concat` is deprecated; use `Series.str.join` instead. Note also "
        "that the default `delimiter` for `str.join` is an empty string, not a hyphen."
    )
    def concat(
        self, delimiter: str | None = None, *, ignore_nulls: bool = True
    ) -> Series:
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
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> pl.Series([1, None, 2]).str.concat("-")  # doctest: +SKIP
        shape: (1,)
        Series: '' [str]
        [
            "1-2"
        ]
        >>> pl.Series([1, None, 2]).str.concat(ignore_nulls=False)  # doctest: +SKIP
        shape: (1,)
        Series: '' [str]
        [
            null
        ]
        """

    def escape_regex(self) -> Series:
        r"""
        Returns string values with all regular expression meta characters escaped.

        Returns
        -------
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> pl.Series(["abc", "def", None, "abc(\\w+)"]).str.escape_regex()
        shape: (4,)
        Series: '' [str]
        [
            "abc"
            "def"
            null
            "abc\(\\w\+\)"
        ]
        """

    def normalize(self, form: UnicodeForm = "NFC") -> Series:
        """
        Returns the Unicode normal form of the string values.

        This uses the forms described in Unicode Standard Annex 15: <https://www.unicode.org/reports/tr15/>.

        Parameters
        ----------
        form : {'NFC', 'NFKC', 'NFD', 'NFKD'}
            Unicode form to use.

        Examples
        --------
        >>> s = pl.Series(["01²", "ＫＡＤＯＫＡＷＡ"])
        >>> s.str.normalize("NFC")
        shape: (2,)
        Series: '' [str]
        [
                "01²"
                "ＫＡＤＯＫＡＷＡ"
        ]
        >>> s.str.normalize("NFKC")
        shape: (2,)
        Series: '' [str]
        [
                "012"
                "KADOKAWA"
        ]
        """  # noqa: RUF002
