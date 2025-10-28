from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Callable

from polars._utils.parse import parse_into_expression
from polars._utils.wrap import wrap_expr

if TYPE_CHECKING:
    from polars import Expr
    from polars._typing import IntoExpr, IntoExprColumn


class ExprArrayNameSpace:
    """Namespace for array related expressions."""

    _accessor = "arr"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def len(self) -> Expr:
        """
        Return the number of elements in each array.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.len())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ u32 │
        ╞═════╡
        │ 2   │
        │ 2   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_len())

    def slice(
        self,
        offset: int | str | Expr,
        length: int | str | Expr | None = None,
        *,
        as_array: bool = False,
    ) -> Expr:
        """
        Slice every subarray.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None` (default), the slice is taken to the
            end of the list.
        as_array
            Return result as a fixed-length `Array`, otherwise as a `List`.
            If true `length` and `offset` must be constant values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.slice(0, 1))
        shape: (2, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ list[i64] │
        ╞═══════════╡
        │ [1]       │
        │ [4]       │
        └───────────┘
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.slice(0, 1, as_array=True))
        shape: (2, 1)
        ┌───────────────┐
        │ a             │
        │ ---           │
        │ array[i64, 1] │
        ╞═══════════════╡
        │ [1]           │
        │ [4]           │
        └───────────────┘
        """
        offset_pyexpr = parse_into_expression(offset)
        length_pyexpr = parse_into_expression(length) if length is not None else None
        return wrap_expr(self._pyexpr.arr_slice(offset_pyexpr, length_pyexpr, as_array))

    def head(self, n: int | str | Expr = 5, *, as_array: bool = False) -> Expr:
        """
        Get the first `n` elements of the sub-arrays.

        Parameters
        ----------
        n
            Number of values to return for each sublist.
        as_array
            Return result as a fixed-length `Array`, otherwise as a `List`.
            If true `n` must be a constant value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.head(1))
        shape: (2, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ list[i64] │
        ╞═══════════╡
        │ [1]       │
        │ [4]       │
        └───────────┘
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.head(1, as_array=True))
        shape: (2, 1)
        ┌───────────────┐
        │ a             │
        │ ---           │
        │ array[i64, 1] │
        ╞═══════════════╡
        │ [1]           │
        │ [4]           │
        └───────────────┘
        """
        return self.slice(0, n, as_array=as_array)

    def tail(self, n: int | str | Expr = 5, *, as_array: bool = False) -> Expr:
        """
        Slice the last `n` values of every sublist.

        Parameters
        ----------
        n
            Number of values to return for each sublist.
        as_array
            Return result as a fixed-length `Array`, otherwise as a `List`.
            If true `n` must be a constant value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.tail(1))
        shape: (2, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ list[i64] │
        ╞═══════════╡
        │ [2]       │
        │ [3]       │
        └───────────┘
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.tail(1, as_array=True))
        shape: (2, 1)
        ┌───────────────┐
        │ a             │
        │ ---           │
        │ array[i64, 1] │
        ╞═══════════════╡
        │ [2]           │
        │ [3]           │
        └───────────────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.arr_tail(n_pyexpr, as_array))

    def min(self) -> Expr:
        """
        Compute the min values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.min())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 3   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_min())

    def max(self) -> Expr:
        """
        Compute the max values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.max())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 4   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_max())

    def sum(self) -> Expr:
        """
        Compute the sum values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.sum())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 3   │
        │ 7   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_sum())

    def std(self, ddof: int = 1) -> Expr:
        """
        Compute the std of the values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.std())
        shape: (2, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ f64      │
        ╞══════════╡
        │ 0.707107 │
        │ 0.707107 │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arr_std(ddof))

    def var(self, ddof: int = 1) -> Expr:
        """
        Compute the var of the values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.var())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 0.5 │
        │ 0.5 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_var(ddof))

    def mean(self) -> Expr:
        """
        Compute the mean of the values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2, 3], [1, 1, 16]]},
        ...     schema={"a": pl.Array(pl.Int64, 3)},
        ... )
        >>> df.select(pl.col("a").arr.mean())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 2.0 │
        │ 6.0 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_mean())

    def median(self) -> Expr:
        """
        Compute the median of the values of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [4, 3]]},
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.select(pl.col("a").arr.median())
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 1.5 │
        │ 3.5 │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_median())

    def unique(self, *, maintain_order: bool = False) -> Expr:
        """
        Get the unique/distinct values in the array.

        Parameters
        ----------
        maintain_order
            Maintain order of data. This requires more work.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[1, 1, 2]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 3)},
        ... )
        >>> df.select(pl.col("a").arr.unique())
        shape: (1, 1)
        ┌───────────┐
        │ a         │
        │ ---       │
        │ list[i64] │
        ╞═══════════╡
        │ [1, 2]    │
        └───────────┘
        """
        return wrap_expr(self._pyexpr.arr_unique(maintain_order))

    def n_unique(self) -> Expr:
        """
        Count the number of unique values in every sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[1, 1, 2], [2, 3, 4]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 3)},
        ... )
        >>> df.with_columns(n_unique=pl.col("a").arr.n_unique())
        shape: (2, 2)
        ┌───────────────┬──────────┐
        │ a             ┆ n_unique │
        │ ---           ┆ ---      │
        │ array[i64, 3] ┆ u32      │
        ╞═══════════════╪══════════╡
        │ [1, 1, 2]     ┆ 2        │
        │ [2, 3, 4]     ┆ 3        │
        └───────────────┴──────────┘
        """
        return wrap_expr(self._pyexpr.arr_n_unique())

    def to_list(self) -> Expr:
        """
        Convert an Array column into a List column with the same inner data type.

        Returns
        -------
        Expr
            Expression of data type :class:`List`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={"a": [[1, 2], [3, 4]]},
        ...     schema={"a": pl.Array(pl.Int8, 2)},
        ... )
        >>> df.select(pl.col("a").arr.to_list())
        shape: (2, 1)
        ┌──────────┐
        │ a        │
        │ ---      │
        │ list[i8] │
        ╞══════════╡
        │ [1, 2]   │
        │ [3, 4]   │
        └──────────┘
        """
        return wrap_expr(self._pyexpr.arr_to_list())

    def any(self) -> Expr:
        """
        Evaluate whether any boolean value is true for every subarray.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "a": [
        ...             [True, True],
        ...             [False, True],
        ...             [False, False],
        ...             [None, None],
        ...             None,
        ...         ]
        ...     },
        ...     schema={"a": pl.Array(pl.Boolean, 2)},
        ... )
        >>> df.with_columns(any=pl.col("a").arr.any())
        shape: (5, 2)
        ┌────────────────┬───────┐
        │ a              ┆ any   │
        │ ---            ┆ ---   │
        │ array[bool, 2] ┆ bool  │
        ╞════════════════╪═══════╡
        │ [true, true]   ┆ true  │
        │ [false, true]  ┆ true  │
        │ [false, false] ┆ false │
        │ [null, null]   ┆ false │
        │ null           ┆ null  │
        └────────────────┴───────┘
        """
        return wrap_expr(self._pyexpr.arr_any())

    def all(self) -> Expr:
        """
        Evaluate whether all boolean values are true for every subarray.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     data={
        ...         "a": [
        ...             [True, True],
        ...             [False, True],
        ...             [False, False],
        ...             [None, None],
        ...             None,
        ...         ]
        ...     },
        ...     schema={"a": pl.Array(pl.Boolean, 2)},
        ... )
        >>> df.with_columns(all=pl.col("a").arr.all())
        shape: (5, 2)
        ┌────────────────┬───────┐
        │ a              ┆ all   │
        │ ---            ┆ ---   │
        │ array[bool, 2] ┆ bool  │
        ╞════════════════╪═══════╡
        │ [true, true]   ┆ true  │
        │ [false, true]  ┆ false │
        │ [false, false] ┆ false │
        │ [null, null]   ┆ true  │
        │ null           ┆ null  │
        └────────────────┴───────┘
        """
        return wrap_expr(self._pyexpr.arr_all())

    def sort(self, *, descending: bool = False, nulls_last: bool = False) -> Expr:
        """
        Sort the arrays in this column.

        Parameters
        ----------
        descending
            Sort in descending order.
        nulls_last
            Place null values last.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[3, 2, 1], [9, 1, 2]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 3)},
        ... )
        >>> df.with_columns(sort=pl.col("a").arr.sort())
        shape: (2, 2)
        ┌───────────────┬───────────────┐
        │ a             ┆ sort          │
        │ ---           ┆ ---           │
        │ array[i64, 3] ┆ array[i64, 3] │
        ╞═══════════════╪═══════════════╡
        │ [3, 2, 1]     ┆ [1, 2, 3]     │
        │ [9, 1, 2]     ┆ [1, 2, 9]     │
        └───────────────┴───────────────┘
        >>> df.with_columns(sort=pl.col("a").arr.sort(descending=True))
        shape: (2, 2)
        ┌───────────────┬───────────────┐
        │ a             ┆ sort          │
        │ ---           ┆ ---           │
        │ array[i64, 3] ┆ array[i64, 3] │
        ╞═══════════════╪═══════════════╡
        │ [3, 2, 1]     ┆ [3, 2, 1]     │
        │ [9, 1, 2]     ┆ [9, 2, 1]     │
        └───────────────┴───────────────┘
        """
        return wrap_expr(self._pyexpr.arr_sort(descending, nulls_last))

    def reverse(self) -> Expr:
        """
        Reverse the arrays in this column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[3, 2, 1], [9, 1, 2]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 3)},
        ... )
        >>> df.with_columns(reverse=pl.col("a").arr.reverse())
        shape: (2, 2)
        ┌───────────────┬───────────────┐
        │ a             ┆ reverse       │
        │ ---           ┆ ---           │
        │ array[i64, 3] ┆ array[i64, 3] │
        ╞═══════════════╪═══════════════╡
        │ [3, 2, 1]     ┆ [1, 2, 3]     │
        │ [9, 1, 2]     ┆ [2, 1, 9]     │
        └───────────────┴───────────────┘
        """
        return wrap_expr(self._pyexpr.arr_reverse())

    def arg_min(self) -> Expr:
        """
        Retrieve the index of the minimal value in every sub-array.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[1, 2], [2, 1]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.with_columns(arg_min=pl.col("a").arr.arg_min())
        shape: (2, 2)
        ┌───────────────┬─────────┐
        │ a             ┆ arg_min │
        │ ---           ┆ ---     │
        │ array[i64, 2] ┆ u32     │
        ╞═══════════════╪═════════╡
        │ [1, 2]        ┆ 0       │
        │ [2, 1]        ┆ 1       │
        └───────────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.arr_arg_min())

    def arg_max(self) -> Expr:
        """
        Retrieve the index of the maximum value in every sub-array.

        Returns
        -------
        Expr
            Expression of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [[1, 2], [2, 1]],
        ...     },
        ...     schema={"a": pl.Array(pl.Int64, 2)},
        ... )
        >>> df.with_columns(arg_max=pl.col("a").arr.arg_max())
        shape: (2, 2)
        ┌───────────────┬─────────┐
        │ a             ┆ arg_max │
        │ ---           ┆ ---     │
        │ array[i64, 2] ┆ u32     │
        ╞═══════════════╪═════════╡
        │ [1, 2]        ┆ 1       │
        │ [2, 1]        ┆ 0       │
        └───────────────┴─────────┘
        """
        return wrap_expr(self._pyexpr.arr_arg_max())

    def get(self, index: int | IntoExprColumn, *, null_on_oob: bool = False) -> Expr:
        """
        Get the value by index in the sub-arrays.

        So index `0` would return the first item of every sublist
        and index `-1` would return the last item of every sublist
        if an index is out of bounds, it will return a `None`.

        Parameters
        ----------
        index
            Index to return per sub-array
        null_on_oob
            Behavior if an index is out of bounds:
            True -> set as null
            False -> raise an error

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"arr": [[1, 2, 3], [4, 5, 6], [7, 8, 9]], "idx": [1, -2, 0]},
        ...     schema={"arr": pl.Array(pl.Int32, 3), "idx": pl.Int32},
        ... )
        >>> df.with_columns(get=pl.col("arr").arr.get("idx", null_on_oob=True))
        shape: (3, 3)
        ┌───────────────┬─────┬─────┐
        │ arr           ┆ idx ┆ get │
        │ ---           ┆ --- ┆ --- │
        │ array[i32, 3] ┆ i32 ┆ i32 │
        ╞═══════════════╪═════╪═════╡
        │ [1, 2, 3]     ┆ 1   ┆ 2   │
        │ [4, 5, 6]     ┆ -2  ┆ 5   │
        │ [7, 8, 9]     ┆ 0   ┆ 7   │
        └───────────────┴─────┴─────┘
        """
        index_pyexpr = parse_into_expression(index)
        return wrap_expr(self._pyexpr.arr_get(index_pyexpr, null_on_oob))

    def first(self) -> Expr:
        """
        Get the first value of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [[1, 2, 3], [4, 5, 6], [7, 8, 9]]},
        ...     schema={"a": pl.Array(pl.Int32, 3)},
        ... )
        >>> df.with_columns(first=pl.col("a").arr.first())
        shape: (3, 2)
        ┌───────────────┬───────┐
        │ a             ┆ first │
        │ ---           ┆ ---   │
        │ array[i32, 3] ┆ i32   │
        ╞═══════════════╪═══════╡
        │ [1, 2, 3]     ┆ 1     │
        │ [4, 5, 6]     ┆ 4     │
        │ [7, 8, 9]     ┆ 7     │
        └───────────────┴───────┘
        """
        return self.get(0, null_on_oob=True)

    def last(self) -> Expr:
        """
        Get the last value of the sub-arrays.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [[1, 2, 3], [4, 5, 6], [7, 9, 8]]},
        ...     schema={"a": pl.Array(pl.Int32, 3)},
        ... )
        >>> df.with_columns(last=pl.col("a").arr.last())
        shape: (3, 2)
        ┌───────────────┬──────┐
        │ a             ┆ last │
        │ ---           ┆ ---  │
        │ array[i32, 3] ┆ i32  │
        ╞═══════════════╪══════╡
        │ [1, 2, 3]     ┆ 3    │
        │ [4, 5, 6]     ┆ 6    │
        │ [7, 9, 8]     ┆ 8    │
        └───────────────┴──────┘
        """
        return self.get(-1, null_on_oob=True)

    def join(self, separator: IntoExprColumn, *, ignore_nulls: bool = True) -> Expr:
        """
        Join all string items in a sub-array and place a separator between them.

        This errors if inner type of array `!= String`.

        Parameters
        ----------
        separator
            string to separate the items with
        ignore_nulls
            Ignore null values (default).

            If set to ``False``, null values will be propagated.
            If the sub-list contains any null values, the output is ``None``.

        Returns
        -------
        Expr
            Expression of data type :class:`String`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"s": [["a", "b"], ["x", "y"]], "separator": ["*", "_"]},
        ...     schema={
        ...         "s": pl.Array(pl.String, 2),
        ...         "separator": pl.String,
        ...     },
        ... )
        >>> df.with_columns(join=pl.col("s").arr.join(pl.col("separator")))
        shape: (2, 3)
        ┌───────────────┬───────────┬──────┐
        │ s             ┆ separator ┆ join │
        │ ---           ┆ ---       ┆ ---  │
        │ array[str, 2] ┆ str       ┆ str  │
        ╞═══════════════╪═══════════╪══════╡
        │ ["a", "b"]    ┆ *         ┆ a*b  │
        │ ["x", "y"]    ┆ _         ┆ x_y  │
        └───────────────┴───────────┴──────┘
        """
        separator_pyexpr = parse_into_expression(separator, str_as_lit=True)
        return wrap_expr(self._pyexpr.arr_join(separator_pyexpr, ignore_nulls))

    def explode(self) -> Expr:
        """
        Returns a column with a separate row for every array element.

        Returns
        -------
        Expr
            Expression with the data type of the array elements.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [[1, 2, 3], [4, 5, 6]]}, schema={"a": pl.Array(pl.Int64, 3)}
        ... )
        >>> df.select(pl.col("a").arr.explode())
        shape: (6, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        │ 4   │
        │ 5   │
        │ 6   │
        └─────┘
        """
        return wrap_expr(self._pyexpr.arr_explode())

    def contains(self, item: IntoExpr, *, nulls_equal: bool = True) -> Expr:
        """
        Check if sub-arrays contain the given item.

        Parameters
        ----------
        item
            Item that will be checked for membership
        nulls_equal : bool, default True
            If True, treat null as a distinct value. Null values will not propagate.

        Returns
        -------
        Expr
            Expression of data type :class:`Boolean`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [["a", "b"], ["x", "y"], ["a", "c"]]},
        ...     schema={"a": pl.Array(pl.String, 2)},
        ... )
        >>> df.with_columns(contains=pl.col("a").arr.contains("a"))
        shape: (3, 2)
        ┌───────────────┬──────────┐
        │ a             ┆ contains │
        │ ---           ┆ ---      │
        │ array[str, 2] ┆ bool     │
        ╞═══════════════╪══════════╡
        │ ["a", "b"]    ┆ true     │
        │ ["x", "y"]    ┆ false    │
        │ ["a", "c"]    ┆ true     │
        └───────────────┴──────────┘
        """
        item_pyexpr = parse_into_expression(item, str_as_lit=True)
        return wrap_expr(self._pyexpr.arr_contains(item_pyexpr, nulls_equal))

    def count_matches(self, element: IntoExpr) -> Expr:
        """
        Count how often the value produced by `element` occurs.

        Parameters
        ----------
        element
            An expression that produces a single value

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [[1, 2], [1, 1], [2, 2]]}, schema={"a": pl.Array(pl.Int64, 2)}
        ... )
        >>> df.with_columns(number_of_twos=pl.col("a").arr.count_matches(2))
        shape: (3, 2)
        ┌───────────────┬────────────────┐
        │ a             ┆ number_of_twos │
        │ ---           ┆ ---            │
        │ array[i64, 2] ┆ u32            │
        ╞═══════════════╪════════════════╡
        │ [1, 2]        ┆ 1              │
        │ [1, 1]        ┆ 0              │
        │ [2, 2]        ┆ 2              │
        └───────────────┴────────────────┘
        """
        element_pyexpr = parse_into_expression(element, str_as_lit=True)
        return wrap_expr(self._pyexpr.arr_count_matches(element_pyexpr))

    def to_struct(
        self, fields: Sequence[str] | Callable[[int], str] | None = None
    ) -> Expr:
        """
        Convert the Series of type `Array` to a Series of type `Struct`.

        Parameters
        ----------
        fields
            If the name and number of the desired fields is known in advance
            a list of field names can be given, which will be assigned by index.
            Otherwise, to dynamically assign field names, a custom function can be
            used; if neither are set, fields will be `field_0, field_1 .. field_n`.

        Examples
        --------
        Convert array to struct with default field name assignment:

        >>> df = pl.DataFrame(
        ...     {"n": [[0, 1, 2], [3, 4, 5]]}, schema={"n": pl.Array(pl.Int8, 3)}
        ... )
        >>> df.with_columns(struct=pl.col("n").arr.to_struct())
        shape: (2, 2)
        ┌──────────────┬───────────┐
        │ n            ┆ struct    │
        │ ---          ┆ ---       │
        │ array[i8, 3] ┆ struct[3] │
        ╞══════════════╪═══════════╡
        │ [0, 1, 2]    ┆ {0,1,2}   │
        │ [3, 4, 5]    ┆ {3,4,5}   │
        └──────────────┴───────────┘

        Convert array to struct with field name assignment by function/index:

        >>> df = pl.DataFrame(
        ...     {"n": [[0, 1, 2], [3, 4, 5]]}, schema={"n": pl.Array(pl.Int8, 3)}
        ... )
        >>> df.select(pl.col("n").arr.to_struct(fields=lambda idx: f"n{idx}")).rows(
        ...     named=True
        ... )
        [{'n': {'n0': 0, 'n1': 1, 'n2': 2}}, {'n': {'n0': 3, 'n1': 4, 'n2': 5}}]

        Convert array to struct with field name assignment by
        index from a list of names:

        >>> df.select(pl.col("n").arr.to_struct(fields=["c1", "c2", "c3"])).rows(
        ...     named=True
        ... )
        [{'n': {'c1': 0, 'c2': 1, 'c3': 2}}, {'n': {'c1': 3, 'c2': 4, 'c3': 5}}]
        """
        if isinstance(fields, Sequence):
            field_names = list(fields)
            pyexpr = self._pyexpr.arr_to_struct(None)
            return wrap_expr(pyexpr).struct.rename_fields(field_names)
        else:
            pyexpr = self._pyexpr.arr_to_struct(fields)
            return wrap_expr(pyexpr)

    def shift(self, n: int | IntoExprColumn = 1) -> Expr:
        """
        Shift array values by the given number of indices.

        Parameters
        ----------
        n
            Number of indices to shift forward. If a negative value is passed, values
            are shifted in the opposite direction instead.

        Notes
        -----
        This method is similar to the `LAG` operation in SQL when the value for `n`
        is positive. With a negative value for `n`, it is similar to `LEAD`.

        Examples
        --------
        By default, array values are shifted forward by one index.

        >>> df = pl.DataFrame(
        ...     {"a": [[1, 2, 3], [4, 5, 6]]}, schema={"a": pl.Array(pl.Int64, 3)}
        ... )
        >>> df.with_columns(shift=pl.col("a").arr.shift())
        shape: (2, 2)
        ┌───────────────┬───────────────┐
        │ a             ┆ shift         │
        │ ---           ┆ ---           │
        │ array[i64, 3] ┆ array[i64, 3] │
        ╞═══════════════╪═══════════════╡
        │ [1, 2, 3]     ┆ [null, 1, 2]  │
        │ [4, 5, 6]     ┆ [null, 4, 5]  │
        └───────────────┴───────────────┘

        Pass a negative value to shift in the opposite direction instead.

        >>> df.with_columns(shift=pl.col("a").arr.shift(-2))
        shape: (2, 2)
        ┌───────────────┬─────────────────┐
        │ a             ┆ shift           │
        │ ---           ┆ ---             │
        │ array[i64, 3] ┆ array[i64, 3]   │
        ╞═══════════════╪═════════════════╡
        │ [1, 2, 3]     ┆ [3, null, null] │
        │ [4, 5, 6]     ┆ [6, null, null] │
        └───────────────┴─────────────────┘
        """
        n_pyexpr = parse_into_expression(n)
        return wrap_expr(self._pyexpr.arr_shift(n_pyexpr))

    def eval(self, expr: Expr, *, as_list: bool = False) -> Expr:
        """
        Run any polars expression against the arrays' elements.

        Parameters
        ----------
        expr
            Expression to run. Note that you can select an element with `pl.element()`
        as_list
            Collect the resulting data as a list. This allows for expressions which
            output a variable amount of data.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 8, 3], "b": [4, 5, 2]})
        >>> df.with_columns(rank=pl.concat_arr("a", "b").arr.eval(pl.element().rank()))
        shape: (3, 3)
        ┌─────┬─────┬───────────────┐
        │ a   ┆ b   ┆ rank          │
        │ --- ┆ --- ┆ ---           │
        │ i64 ┆ i64 ┆ array[f64, 2] │
        ╞═════╪═════╪═══════════════╡
        │ 1   ┆ 4   ┆ [1.0, 2.0]    │
        │ 8   ┆ 5   ┆ [2.0, 1.0]    │
        │ 3   ┆ 2   ┆ [2.0, 1.0]    │
        └─────┴─────┴───────────────┘

        See Also
        --------
        polars.Expr.arr.agg: Evaluate any expression and automatically explode.
        polars.Expr.list.eval: Same for the List datatype.
        """
        return wrap_expr(self._pyexpr.arr_eval(expr._pyexpr, as_list=as_list))

    def agg(self, expr: Expr) -> Expr:
        """
        Run any polars aggregation expression against the arrays' elements.

        Parameters
        ----------
        expr
            Expression to run. Note that you can select an element with `pl.element()`.

        Examples
        --------
        >>> df = pl.Series(
        ...     "a", [[1, None], [42, 13], [None, None]], pl.Array(pl.Int64, 2)
        ... ).to_frame()
        >>> df.with_columns(null_count=pl.col.a.arr.agg(pl.element().null_count()))
        shape: (3, 2)
        ┌───────────────┬────────────┐
        │ a             ┆ null_count │
        │ ---           ┆ ---        │
        │ array[i64, 2] ┆ u32        │
        ╞═══════════════╪════════════╡
        │ [1, null]     ┆ 1          │
        │ [42, 13]      ┆ 0          │
        │ [null, null]  ┆ 2          │
        └───────────────┴────────────┘
        >>> df.with_columns(no_nulls=pl.col.a.arr.agg(pl.element().drop_nulls()))
        shape: (3, 2)
        ┌───────────────┬───────────┐
        │ a             ┆ no_nulls  │
        │ ---           ┆ ---       │
        │ array[i64, 2] ┆ list[i64] │
        ╞═══════════════╪═══════════╡
        │ [1, null]     ┆ [1]       │
        │ [42, 13]      ┆ [42, 13]  │
        │ [null, null]  ┆ []        │
        └───────────────┴───────────┘

        See Also
        --------
        polars.Expr.arr.eval: Evaluate any expression without automatic explode.
        polars.Expr.list.agg: Same for the List datatype.
        """
        return wrap_expr(self._pyexpr.arr_agg(expr._pyexpr))
