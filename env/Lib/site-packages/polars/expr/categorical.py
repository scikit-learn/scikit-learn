from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.wrap import wrap_expr

if TYPE_CHECKING:
    from polars import Expr


class ExprCatNameSpace:
    """Namespace for categorical related expressions."""

    _accessor = "cat"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def get_categories(self) -> Expr:
        """
        Get the categories stored in this data type.

        Examples
        --------
        >>> df = pl.Series(
        ...     "cats", ["foo", "bar", "foo", "foo", "ham"], dtype=pl.Categorical
        ... ).to_frame()
        >>> df.select(pl.col("cats").cat.get_categories())
        shape: (3, 1)
        ┌──────┐
        │ cats │
        │ ---  │
        │ str  │
        ╞══════╡
        │ foo  │
        │ bar  │
        │ ham  │
        └──────┘
        """
        return wrap_expr(self._pyexpr.cat_get_categories())

    def len_bytes(self) -> Expr:
        """
        Return the byte-length of the string representation of each value.

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
        >>> df = pl.DataFrame(
        ...     {"a": pl.Series(["Café", "345", "東京", None], dtype=pl.Categorical)}
        ... )
        >>> df.with_columns(
        ...     pl.col("a").cat.len_bytes().alias("n_bytes"),
        ...     pl.col("a").cat.len_chars().alias("n_chars"),
        ... )
        shape: (4, 3)
        ┌──────┬─────────┬─────────┐
        │ a    ┆ n_bytes ┆ n_chars │
        │ ---  ┆ ---     ┆ ---     │
        │ cat  ┆ u32     ┆ u32     │
        ╞══════╪═════════╪═════════╡
        │ Café ┆ 5       ┆ 4       │
        │ 345  ┆ 3       ┆ 3       │
        │ 東京 ┆ 6       ┆ 2       │
        │ null ┆ null    ┆ null    │
        └──────┴─────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.cat_len_bytes())

    def len_chars(self) -> Expr:
        """
        Return the number of characters of the string representation of each value.

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
        >>> df = pl.DataFrame(
        ...     {"a": pl.Series(["Café", "345", "東京", None], dtype=pl.Categorical)}
        ... )
        >>> df.with_columns(
        ...     pl.col("a").cat.len_chars().alias("n_chars"),
        ...     pl.col("a").cat.len_bytes().alias("n_bytes"),
        ... )
        shape: (4, 3)
        ┌──────┬─────────┬─────────┐
        │ a    ┆ n_chars ┆ n_bytes │
        │ ---  ┆ ---     ┆ ---     │
        │ cat  ┆ u32     ┆ u32     │
        ╞══════╪═════════╪═════════╡
        │ Café ┆ 4       ┆ 5       │
        │ 345  ┆ 3       ┆ 3       │
        │ 東京 ┆ 2       ┆ 6       │
        │ null ┆ null    ┆ null    │
        └──────┴─────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.cat_len_chars())
