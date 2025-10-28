from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.unstable import unstable
from polars._utils.wrap import wrap_s
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    from polars import Series
    from polars._plr import PySeries


@expr_dispatch
class CatNameSpace:
    """Namespace for categorical related series."""

    _accessor = "cat"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def get_categories(self) -> Series:
        """
        Get the categories stored in this data type.

        Examples
        --------
        >>> s = pl.Series(["foo", "bar", "foo", "foo", "ham"], dtype=pl.Categorical)
        >>> s.cat.get_categories()  # doctest: +SKIP
        shape: (3,)
        Series: '' [str]
        [
            "foo"
            "bar"
            "ham"
        ]
        """

    def is_local(self) -> bool:
        """
        Return whether or not the column is a local categorical.

        Always returns false.
        """
        return self._s.cat_is_local()

    def to_local(self) -> Series:
        """Simply returns the column as-is, local representations are deprecated."""
        return wrap_s(self._s.cat_to_local())

    @unstable()
    def uses_lexical_ordering(self) -> bool:
        """
        Indicate whether the Series uses lexical ordering.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Examples
        --------
        >>> s = pl.Series(["b", "a", "b"]).cast(pl.Categorical)
        >>> s.cat.uses_lexical_ordering()
        True
        """
        return self._s.cat_uses_lexical_ordering()

    def len_bytes(self) -> Series:
        """
        Return the byte-length of the string representation of each value.

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
        >>> s = pl.Series(["Café", "345", "東京", None], dtype=pl.Categorical)
        >>> s.cat.len_bytes()
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
        Return the number of characters of the string representation of each value.

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
        >>> s = pl.Series(["Café", "345", "東京", None], dtype=pl.Categorical)
        >>> s.cat.len_chars()
        shape: (4,)
        Series: '' [u32]
        [
            4
            3
            2
            null
        ]
        """

    def starts_with(self, prefix: str) -> Series:
        """
        Check if string representations of values start with a substring.

        Parameters
        ----------
        prefix
            Prefix substring.

        See Also
        --------
        contains : Check if the string repr contains a substring that matches a pattern.
        ends_with : Check if string repr ends with a substring.

        Examples
        --------
        >>> s = pl.Series("fruits", ["apple", "mango", None], dtype=pl.Categorical)
        >>> s.cat.starts_with("app")
        shape: (3,)
        Series: 'fruits' [bool]
        [
            true
            false
            null
        ]
        """

    def ends_with(self, suffix: str) -> Series:
        """
        Check if string representations of values end with a substring.

        Parameters
        ----------
        suffix
            Suffix substring.

        See Also
        --------
        contains : Check if the string repr contains a substring that matches a pattern.
        starts_with : Check if string repr starts with a substring.

        Examples
        --------
        >>> s = pl.Series("fruits", ["apple", "mango", None], dtype=pl.Categorical)
        >>> s.cat.ends_with("go")
        shape: (3,)
        Series: 'fruits' [bool]
        [
            false
            true
            null
        ]
        """

    def slice(self, offset: int, length: int | None = None) -> Series:
        """
        Extract a substring from the string representation of each string value.

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
        >>> s = pl.Series(["pear", None, "papaya", "dragonfruit"], dtype=pl.Categorical)
        >>> s.cat.slice(-3)
        shape: (4,)
        Series: '' [str]
        [
            "ear"
            null
            "aya"
            "uit"
        ]

        Using the optional `length` parameter

        >>> s.cat.slice(4, length=3)
        shape: (4,)
        Series: '' [str]
        [
            ""
            null
            "ya"
            "onf"
        ]
        """
