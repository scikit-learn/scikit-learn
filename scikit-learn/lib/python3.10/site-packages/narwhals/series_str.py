from __future__ import annotations

from typing import Any, Generic

from narwhals.dependencies import is_narwhals_series
from narwhals.typing import SeriesT


class SeriesStringNamespace(Generic[SeriesT]):
    def __init__(self, series: SeriesT) -> None:
        self._narwhals_series = series

    def len_chars(self) -> SeriesT:
        r"""Return the length of each string as the number of characters.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(["foo", "345", None])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.len_chars().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [u32]
            [
                    3
                    3
                    null
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.len_chars()
        )

    def _extract_compliant(self, arg: Any) -> Any:
        return arg._compliant_series if is_narwhals_series(arg) else arg

    def replace(
        self, pattern: str, value: str | SeriesT, *, literal: bool = False, n: int = 1
    ) -> SeriesT:
        r"""Replace first matching regex/literal substring with a new string value.

        Arguments:
            pattern: A valid regular expression pattern.
            value: String that will replace the matched substring.
            literal: Treat `pattern` as a literal string.
            n: Number of matches to replace.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["123abc", "abc abc123"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.replace("abc", "").to_native()
            0        123
            1     abc123
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.replace(
                pattern, self._extract_compliant(value), literal=literal, n=n
            )
        )

    def replace_all(
        self, pattern: str, value: str | SeriesT, *, literal: bool = False
    ) -> SeriesT:
        r"""Replace all matching regex/literal substring with a new string value.

        Arguments:
            pattern: A valid regular expression pattern.
            value: String that will replace the matched substring.
            literal: Treat `pattern` as a literal string.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["123abc", "abc abc123"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.replace_all("abc", "").to_native()
            0     123
            1     123
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.replace_all(
                pattern, self._extract_compliant(value), literal=literal
            )
        )

    def strip_chars(self, characters: str | None = None) -> SeriesT:
        r"""Remove leading and trailing characters.

        Arguments:
            characters: The set of characters to be removed. All combinations of this set of characters will be stripped from the start and end of the string. If set to None (default), all leading and trailing whitespace is removed instead.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(["apple", "\nmango"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.strip_chars().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [str]
            [
                    "apple"
                    "mango"
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.strip_chars(characters)
        )

    def starts_with(self, prefix: str) -> SeriesT:
        r"""Check if string values start with a substring.

        Arguments:
            prefix: prefix substring

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["apple", "mango", None])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.starts_with("app").to_native()
            0     True
            1    False
            2     None
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.starts_with(prefix)
        )

    def ends_with(self, suffix: str) -> SeriesT:
        r"""Check if string values end with a substring.

        Arguments:
            suffix: suffix substring

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["apple", "mango", None])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.ends_with("ngo").to_native()
            0    False
            1     True
            2     None
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.ends_with(suffix)
        )

    def contains(self, pattern: str, *, literal: bool = False) -> SeriesT:
        r"""Check if string contains a substring that matches a pattern.

        Arguments:
            pattern: A Character sequence or valid regular expression pattern.
            literal: If True, treats the pattern as a literal string.
                     If False, assumes the pattern is a regular expression.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([["cat", "dog", "rabbit and parrot"]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.contains("cat|parrot").to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                true,
                false,
                true
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.contains(pattern, literal=literal)
        )

    def slice(self, offset: int, length: int | None = None) -> SeriesT:
        r"""Create subslices of the string values of a Series.

        Arguments:
            offset: Start index. Negative indexing is supported.
            length: Length of the slice. If set to `None` (default), the slice is taken to the
                end of the string.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["pear", None, "papaya"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.slice(4, 3).to_native()  # doctest: +NORMALIZE_WHITESPACE
            0
            1    None
            2      ya
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.slice(
                offset=offset, length=length
            )
        )

    def split(self, by: str) -> SeriesT:
        r"""Split the string values of a Series by a substring.

        Arguments:
            by: Substring to split by.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(["foo bar", "foo_bar"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.split("_").to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [list[str]]
            [
                    ["foo bar"]
                    ["foo", "bar"]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.split(by=by)
        )

    def head(self, n: int = 5) -> SeriesT:
        r"""Take the first n elements of each string.

        Arguments:
            n: Number of elements to take. Negative indexing is supported (see note (1.))

        Notes:
            1. When the `n` input is negative, `head` returns characters up to the n-th from the end of the string.
            For example, if `n = -3`, then all characters except the last three are returned.
            2. If the length of the string has fewer than `n` characters, the full string is returned.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([["taata", "taatatata", "zukkyun"]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.head().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                "taata",
                "taata",
                "zukky"
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.slice(offset=0, length=n)
        )

    def tail(self, n: int = 5) -> SeriesT:
        r"""Take the last n elements of each string.

        Arguments:
            n: Number of elements to take. Negative indexing is supported (see note (1.))

        Notes:
            1. When the `n` input is negative, `tail` returns characters starting from the n-th from the beginning of
            the string. For example, if `n = -3`, then all characters except the first three are returned.
            2. If the length of the string has fewer than `n` characters, the full string is returned.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([["taata", "taatatata", "zukkyun"]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.tail().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                "taata",
                "atata",
                "kkyun"
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.slice(offset=-n, length=None)
        )

    def to_uppercase(self) -> SeriesT:
        r"""Transform string to uppercase variant.

        Notes:
            The PyArrow backend will convert 'ß' to 'ẞ' instead of 'SS'.
            For more info see: https://github.com/apache/arrow/issues/34599
            There may be other unicode-edge-case-related variations across implementations.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["apple", None])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.to_uppercase().to_native()
            0    APPLE
            1     None
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.to_uppercase()
        )

    def to_lowercase(self) -> SeriesT:
        r"""Transform string to lowercase variant.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["APPLE", None])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.to_lowercase().to_native()
            0    apple
            1     None
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.to_lowercase()
        )

    def to_datetime(self, format: str | None = None) -> SeriesT:
        """Parse Series with strings to a Series with Datetime dtype.

        Notes:
            - pandas defaults to nanosecond time unit, Polars to microsecond.
              Prior to pandas 2.0, nanoseconds were the only time unit supported
              in pandas, with no ability to set any other one. The ability to
              set the time unit in pandas, if the version permits, will arrive.
            - timezone-aware strings are all converted to and parsed as UTC.

        Warning:
            As different backends auto-infer format in different ways, if `format=None`
            there is no guarantee that the result will be equal.

        Arguments:
            format: Format to use for conversion. If set to None (default), the format is
                inferred from the data.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(["2020-01-01", "2020-01-02"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.to_datetime(
            ...     format="%Y-%m-%d"
            ... ).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (2,)
            Series: '' [datetime[μs]]
            [
                    2020-01-01 00:00:00
                    2020-01-02 00:00:00
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.to_datetime(format=format)
        )

    def to_date(self, format: str | None = None) -> SeriesT:
        """Convert to date dtype.

        Warning:
            As different backends auto-infer format in different ways, if `format=None`
            there is no guarantee that the result will be equal.

        Arguments:
            format: Format to use for conversion. If set to None (default), the format is
                inferred from the data.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([["2020-01-01", "2020-01-02"]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.to_date(format="%Y-%m-%d").to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2020-01-01,
                2020-01-02
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.to_date(format=format)
        )

    def to_titlecase(self) -> SeriesT:
        """Modify strings to their titlecase equivalent.

        Notes:
            This is a form of case transform where the first letter of each word is
            capitalized, with the rest of the word in lowercase.

        Warning:
            Different backends might follow different rules to determine what a "word" is:

            - polars uses **non-alphanumeric** characters to define the word boundaries.
            - pandas-like and pyarrow use **non-alphabetic** characters to define
                the word boundaries, matching the behavior of
                [`str.title`](https://docs.python.org/3/library/stdtypes.html#str.title).

            As an example of such difference, in the former case the string `"with123numbers"`
            is mapped to `"With123numbers"` (notice lowercase **n** after the digits), while
            in the latter to `"With123Numbers"` (notice uppercase **N** after the digits).

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array(
            ...     [
            ...         [
            ...             "'e.t. phone home'",
            ...             "you talkin' to me?",
            ...             "to infinity,and BEYOND!",
            ...         ]
            ...     ]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.to_titlecase().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                "'E.T. Phone Home'",
                "You Talkin' To Me?",
                "To Infinity,And Beyond!"
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.to_titlecase()
        )

    def zfill(self, width: int) -> SeriesT:
        r"""Pad strings with zeros on the left.

        Arguments:
            width: The target width of the string. If the string is shorter than this width, it will be padded with zeros on the left.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["+1", "-23", "456", "123456"])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.str.zfill(5).to_native()
            0     +0001
            1     -0023
            2     00456
            3    123456
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.str.zfill(width)
        )
