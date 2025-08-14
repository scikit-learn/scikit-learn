from __future__ import annotations

from typing import TYPE_CHECKING, Any, Generic, TypeVar

from narwhals._expression_parsing import all_exprs_are_scalar_like
from narwhals._utils import flatten, tupleify
from narwhals.exceptions import InvalidOperationError
from narwhals.typing import DataFrameT

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

    from narwhals._compliant.typing import CompliantExprAny
    from narwhals.dataframe import LazyFrame
    from narwhals.expr import Expr

LazyFrameT = TypeVar("LazyFrameT", bound="LazyFrame[Any]")


class GroupBy(Generic[DataFrameT]):
    def __init__(
        self,
        df: DataFrameT,
        keys: Sequence[str] | Sequence[CompliantExprAny],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._df: DataFrameT = df
        self._keys = keys
        self._grouped = self._df._compliant_frame.group_by(
            self._keys, drop_null_keys=drop_null_keys
        )

    def agg(self, *aggs: Expr | Iterable[Expr], **named_aggs: Expr) -> DataFrameT:
        """Compute aggregations for each group of a group by operation.

        Arguments:
            aggs: Aggregations to compute for each group of the group by operation,
                specified as positional arguments.
            named_aggs: Additional aggregations, specified as keyword arguments.

        Returns:
            A new Dataframe.

        Examples:
            Group by one column or by multiple columns and call `agg` to compute
            the grouped sum of another column.

            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {
            ...         "a": ["a", "b", "a", "b", "c"],
            ...         "b": [1, 2, 1, 3, 3],
            ...         "c": [5, 4, 3, 2, 1],
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>>
            >>> df.group_by("a").agg(nw.col("b").sum()).sort("a")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  a  2      |
            |     1  b  5      |
            |     2  c  3      |
            └──────────────────┘
            >>>
            >>> df.group_by("a", "b").agg(nw.col("c").sum()).sort("a", "b").to_native()
               a  b  c
            0  a  1  8
            1  b  2  4
            2  b  3  2
            3  c  3  1
        """
        flat_aggs = tuple(flatten(aggs))
        if not all_exprs_are_scalar_like(*flat_aggs, **named_aggs):
            msg = (
                "Found expression which does not aggregate.\n\n"
                "All expressions passed to GroupBy.agg must aggregate.\n"
                "For example, `df.group_by('a').agg(nw.col('b').sum())` is valid,\n"
                "but `df.group_by('a').agg(nw.col('b'))` is not."
            )
            raise InvalidOperationError(msg)
        plx = self._df.__narwhals_namespace__()
        compliant_aggs = (
            *(x._to_compliant_expr(plx) for x in flat_aggs),
            *(
                value.alias(key)._to_compliant_expr(plx)
                for key, value in named_aggs.items()
            ),
        )
        return self._df._with_compliant(self._grouped.agg(*compliant_aggs))

    def __iter__(self) -> Iterator[tuple[Any, DataFrameT]]:
        yield from (
            (tupleify(key), self._df._with_compliant(df))
            for (key, df) in self._grouped.__iter__()
        )


class LazyGroupBy(Generic[LazyFrameT]):
    def __init__(
        self,
        df: LazyFrameT,
        keys: Sequence[str] | Sequence[CompliantExprAny],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._df: LazyFrameT = df
        self._keys = keys
        self._grouped = self._df._compliant_frame.group_by(
            self._keys, drop_null_keys=drop_null_keys
        )

    def agg(self, *aggs: Expr | Iterable[Expr], **named_aggs: Expr) -> LazyFrameT:
        """Compute aggregations for each group of a group by operation.

        Arguments:
            aggs: Aggregations to compute for each group of the group by operation,
                specified as positional arguments.
            named_aggs: Additional aggregations, specified as keyword arguments.

        Returns:
            A new LazyFrame.

        Examples:
            Group by one column or by multiple columns and call `agg` to compute
            the grouped sum of another column.

            >>> import polars as pl
            >>> import narwhals as nw
            >>> from narwhals.typing import IntoFrameT
            >>> lf_native = pl.LazyFrame(
            ...     {
            ...         "a": ["a", "b", "a", "b", "c"],
            ...         "b": [1, 2, 1, 3, 3],
            ...         "c": [5, 4, 3, 2, 1],
            ...     }
            ... )
            >>> lf = nw.from_native(lf_native)
            >>>
            >>> nw.to_native(lf.group_by("a").agg(nw.col("b").sum()).sort("a")).collect()
            shape: (3, 2)
            ┌─────┬─────┐
            │ a   ┆ b   │
            │ --- ┆ --- │
            │ str ┆ i64 │
            ╞═════╪═════╡
            │ a   ┆ 2   │
            │ b   ┆ 5   │
            │ c   ┆ 3   │
            └─────┴─────┘
            >>>
            >>> lf.group_by("a", "b").agg(nw.sum("c")).sort("a", "b").collect()
            ┌───────────────────┐
            |Narwhals DataFrame |
            |-------------------|
            |shape: (4, 3)      |
            |┌─────┬─────┬─────┐|
            |│ a   ┆ b   ┆ c   │|
            |│ --- ┆ --- ┆ --- │|
            |│ str ┆ i64 ┆ i64 │|
            |╞═════╪═════╪═════╡|
            |│ a   ┆ 1   ┆ 8   │|
            |│ b   ┆ 2   ┆ 4   │|
            |│ b   ┆ 3   ┆ 2   │|
            |│ c   ┆ 3   ┆ 1   │|
            |└─────┴─────┴─────┘|
            └───────────────────┘
        """
        flat_aggs = tuple(flatten(aggs))
        if not all_exprs_are_scalar_like(*flat_aggs, **named_aggs):
            msg = (
                "Found expression which does not aggregate.\n\n"
                "All expressions passed to GroupBy.agg must aggregate.\n"
                "For example, `df.group_by('a').agg(nw.col('b').sum())` is valid,\n"
                "but `df.group_by('a').agg(nw.col('b'))` is not."
            )
            raise InvalidOperationError(msg)
        plx = self._df.__narwhals_namespace__()
        compliant_aggs = (
            *(x._to_compliant_expr(plx) for x in flat_aggs),
            *(
                value.alias(key)._to_compliant_expr(plx)
                for key, value in named_aggs.items()
            ),
        )
        return self._df._with_compliant(self._grouped.agg(*compliant_aggs))
