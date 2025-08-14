from __future__ import annotations

from typing import TYPE_CHECKING, Generic, TypeVar

if TYPE_CHECKING:
    from narwhals.expr import Expr
    from narwhals.typing import NonNestedLiteral

ExprT = TypeVar("ExprT", bound="Expr")


class ExprListNamespace(Generic[ExprT]):
    def __init__(self, expr: ExprT) -> None:
        self._expr = expr

    def len(self) -> ExprT:
        """Return the number of elements in each list.

        Null values count towards the total.

        Returns:
            A new expression.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 2], [3, 4, None], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_len=nw.col("a").list.len())
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |shape: (4, 2)           |
            |┌──────────────┬───────┐|
            |│ a            ┆ a_len │|
            |│ ---          ┆ ---   │|
            |│ list[i64]    ┆ u32   │|
            |╞══════════════╪═══════╡|
            |│ [1, 2]       ┆ 2     │|
            |│ [3, 4, null] ┆ 3     │|
            |│ null         ┆ null  │|
            |│ []           ┆ 0     │|
            |└──────────────┴───────┘|
            └────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).list.len()
        )

    def unique(self) -> ExprT:
        """Get the unique/distinct values in the list.

        Null values are included in the result. The order of unique values is not guaranteed.

        Returns:
            A new expression.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 1, 2], [3, 3, None], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_unique=nw.col("a").list.unique())
            ┌────────────────────────────┐
            |     Narwhals DataFrame     |
            |----------------------------|
            |shape: (4, 2)               |
            |┌──────────────┬───────────┐|
            |│ a            ┆ a_unique  │|
            |│ ---          ┆ ---       │|
            |│ list[i64]    ┆ list[i64] │|
            |╞══════════════╪═══════════╡|
            |│ [1, 1, 2]    ┆ [1, 2]    │|
            |│ [3, 3, null] ┆ [null, 3] │|
            |│ null         ┆ null      │|
            |│ []           ┆ []        │|
            |└──────────────┴───────────┘|
            └────────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).list.unique()
        )

    def contains(self, item: NonNestedLiteral) -> ExprT:
        """Check if sublists contain the given item.

        Arguments:
            item: Item that will be checked for membership.

        Returns:
            A new expression.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 2], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_contains_1=nw.col("a").list.contains(1))
            ┌────────────────────────────┐
            |     Narwhals DataFrame     |
            |----------------------------|
            |shape: (3, 2)               |
            |┌───────────┬──────────────┐|
            |│ a         ┆ a_contains_1 │|
            |│ ---       ┆ ---          │|
            |│ list[i64] ┆ bool         │|
            |╞═══════════╪══════════════╡|
            |│ [1, 2]    ┆ true         │|
            |│ null      ┆ null         │|
            |│ []        ┆ false        │|
            |└───────────┴──────────────┘|
            └────────────────────────────┘
        """
        return self._expr._with_elementwise(
            lambda plx: self._expr._to_compliant_expr(plx).list.contains(item)
        )
