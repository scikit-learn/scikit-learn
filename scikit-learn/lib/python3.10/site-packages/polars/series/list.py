from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable

from polars import functions as F
from polars._utils.unstable import unstable
from polars._utils.wrap import wrap_s
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    from collections.abc import Collection

    from polars import Expr, Series
    from polars._plr import PySeries
    from polars._typing import (
        IntoExpr,
        IntoExprColumn,
        ListToStructWidthStrategy,
        NullBehavior,
    )


@expr_dispatch
class ListNameSpace:
    """Namespace for list related methods."""

    _accessor = "list"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def all(self) -> Series:
        """
        Evaluate whether all boolean values in a list are true.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Notes
        -----
        If there are no non-null elements in a row, the output is `True`.

        Examples
        --------
        >>> s = pl.Series(
        ...     [[True, True], [False, True], [False, False], [None], [], None],
        ...     dtype=pl.List(pl.Boolean),
        ... )
        >>> s.list.all()
        shape: (6,)
        Series: '' [bool]
        [
            true
            false
            false
            true
            true
            null
        ]
        """

    def any(self) -> Series:
        """
        Evaluate whether any boolean value in a list is true.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Notes
        -----
        If there are no non-null elements in a row, the output is `False`.

        Examples
        --------
        >>> s = pl.Series(
        ...     [[True, True], [False, True], [False, False], [None], [], None],
        ...     dtype=pl.List(pl.Boolean),
        ... )
        >>> s.list.any()
        shape: (6,)
        Series: '' [bool]
        [
            true
            true
            false
            false
            false
            null
        ]
        """

    def len(self) -> Series:
        """
        Return the number of elements in each list.

        Null values count towards the total.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

        Examples
        --------
        >>> s = pl.Series([[1, 2, None], [5]])
        >>> s.list.len()
        shape: (2,)
        Series: '' [u32]
        [
            3
            1
        ]
        """

    def drop_nulls(self) -> Series:
        """
        Drop all null values in the list.

        The original order of the remaining elements is preserved.

        Examples
        --------
        >>> s = pl.Series("values", [[None, 1, None, 2], [None], [3, 4]])
        >>> s.list.drop_nulls()
        shape: (3,)
        Series: 'values' [list[i64]]
        [
            [1, 2]
            []
            [3, 4]
        ]
        """

    def sample(
        self,
        n: int | IntoExprColumn | None = None,
        *,
        fraction: float | IntoExprColumn | None = None,
        with_replacement: bool = False,
        shuffle: bool = False,
        seed: int | None = None,
    ) -> Series:
        """
        Sample from this list.

        Parameters
        ----------
        n
            Number of items to return. Cannot be used with `fraction`. Defaults to 1 if
            `fraction` is None.
        fraction
            Fraction of items to return. Cannot be used with `n`.
        with_replacement
            Allow values to be sampled more than once.
        shuffle
            Shuffle the order of sampled data points.
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated for each sample operation.

        Examples
        --------
        >>> s = pl.Series("values", [[1, 2, 3], [4, 5]])
        >>> s.list.sample(n=pl.Series("n", [2, 1]), seed=1)
        shape: (2,)
        Series: 'values' [list[i64]]
        [
            [2, 3]
            [5]
        ]
        """

    def sum(self) -> Series:
        """
        Sum all the arrays in the list.

        Notes
        -----
        If there are no non-null elements in a row, the output is `0`.

        Examples
        --------
        >>> s = pl.Series("values", [[1], [2, 3]])
        >>> s.list.sum()
        shape: (2,)
        Series: 'values' [i64]
        [
            1
            5
        ]
        """

    def max(self) -> Series:
        """
        Compute the max value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[4, 1], [2, 3]])
        >>> s.list.max()
        shape: (2,)
        Series: 'values' [i64]
        [
            4
            3
        ]
        """

    def min(self) -> Series:
        """
        Compute the min value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[4, 1], [2, 3]])
        >>> s.list.min()
        shape: (2,)
        Series: 'values' [i64]
        [
            1
            2
        ]
        """

    def mean(self) -> Series:
        """
        Compute the mean value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[3, 1], [3, 3]])
        >>> s.list.mean()
        shape: (2,)
        Series: 'values' [f64]
        [
            2.0
            3.0
        ]
        """

    def median(self) -> Series:
        """
        Compute the median value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[-1, 0, 1], [1, 10]])
        >>> s.list.median()
        shape: (2,)
        Series: 'values' [f64]
        [
                0.0
                5.5
        ]
        """

    def std(self, ddof: int = 1) -> Series:
        """
        Compute the std value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[-1, 0, 1], [1, 10]])
        >>> s.list.std()
        shape: (2,)
        Series: 'values' [f64]
        [
                1.0
                6.363961
        ]
        """

    def var(self, ddof: int = 1) -> Series:
        """
        Compute the var value of the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("values", [[-1, 0, 1], [1, 10]])
        >>> s.list.var()
        shape: (2,)
        Series: 'values' [f64]
        [
                1.0
                40.5
        ]
        """

    def sort(
        self,
        *,
        descending: bool = False,
        nulls_last: bool = False,
        multithreaded: bool = True,
    ) -> Series:
        """
        Sort the arrays in this column.

        Parameters
        ----------
        descending
            Sort in descending order.
        nulls_last
            Place null values last.
        multithreaded
            Sort using multiple threads.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [9, 1, 2]])
        >>> s.list.sort()
        shape: (2,)
        Series: 'a' [list[i64]]
        [
                [1, 2, 3]
                [1, 2, 9]
        ]
        >>> s.list.sort(descending=True)
        shape: (2,)
        Series: 'a' [list[i64]]
        [
                [3, 2, 1]
                [9, 2, 1]
        ]
        """

    def reverse(self) -> Series:
        """
        Reverse the arrays in the list.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [9, 1, 2]])
        >>> s.list.reverse()
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [1, 2, 3]
            [2, 1, 9]
        ]
        """

    def unique(self, *, maintain_order: bool = False) -> Series:
        """
        Get the unique/distinct values in the list.

        Parameters
        ----------
        maintain_order
            Maintain order of data. This requires more work.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 1, 2], [2, 3, 3]])
        >>> s.list.unique()
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [1, 2]
            [2, 3]
        ]
        """

    def n_unique(self) -> Series:
        """
        Count the number of unique values in every sub-lists.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 1, 2], [2, 3, 4]])
        >>> s.list.n_unique()
        shape: (2,)
        Series: 'a' [u32]
        [
            2
            3
        ]
        """

    def concat(self, other: list[Series] | Series | list[Any]) -> Series:
        """
        Concat the arrays in a Series dtype List in linear time.

        Parameters
        ----------
        other
            Columns to concat into a List Series

        Examples
        --------
        >>> s1 = pl.Series("a", [["a", "b"], ["c"]])
        >>> s2 = pl.Series("b", [["c"], ["d", None]])
        >>> s1.list.concat(s2)
        shape: (2,)
        Series: 'a' [list[str]]
        [
            ["a", "b", "c"]
            ["c", "d", null]
        ]
        """

    def get(
        self,
        index: int | Series | list[int],
        *,
        null_on_oob: bool = False,
    ) -> Series:
        """
        Get the value by index in the sublists.

        So index `0` would return the first item of every sublist
        and index `-1` would return the last item of every sublist
        if an index is out of bounds, it will return a `None`.

        Parameters
        ----------
        index
            Index to return per sublist
        null_on_oob
            Behavior if an index is out of bounds:

            * True -> set as null
            * False -> raise an error

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [], [1, 2]])
        >>> s.list.get(0, null_on_oob=True)
        shape: (3,)
        Series: 'a' [i64]
        [
            3
            null
            1
        ]
        """

    def gather(
        self,
        indices: Series | list[int] | list[list[int]],
        *,
        null_on_oob: bool = False,
    ) -> Series:
        """
        Take sublists by multiple indices.

        The indices may be defined in a single column, or by sublists in another
        column of dtype `List`.

        Parameters
        ----------
        indices
            Indices to return per sublist
        null_on_oob
            Behavior if an index is out of bounds:
            True -> set as null
            False -> raise an error
            Note that defaulting to raising an error is much cheaper

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [], [1, 2]])
        >>> s.list.gather([0, 2], null_on_oob=True)
        shape: (3,)
        Series: 'a' [list[i64]]
        [
            [3, 1]
            [null, null]
            [1, null]
        ]
        """

    def gather_every(
        self, n: int | IntoExprColumn, offset: int | IntoExprColumn = 0
    ) -> Series:
        """
        Take every n-th value start from offset in sublists.

        Parameters
        ----------
        n
            Gather every n-th element.
        offset
            Starting index.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3], [], [6, 7, 8, 9]])
        >>> s.list.gather_every(2, offset=1)
        shape: (3,)
        Series: 'a' [list[i64]]
        [
            [2]
            []
            [7, 9]
        ]
        """

    def __getitem__(self, item: int) -> Series:
        return self.get(item)

    def join(self, separator: IntoExprColumn, *, ignore_nulls: bool = True) -> Series:
        """
        Join all string items in a sublist and place a separator between them.

        This errors if inner type of list `!= String`.

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
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> s = pl.Series([["foo", "bar"], ["hello", "world"]])
        >>> s.list.join(separator="-")
        shape: (2,)
        Series: '' [str]
        [
            "foo-bar"
            "hello-world"
        ]
        """

    def first(self) -> Series:
        """
        Get the first value of the sublists.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [], [1, 2]])
        >>> s.list.first()
        shape: (3,)
        Series: 'a' [i64]
        [
            3
            null
            1
        ]
        """

    def last(self) -> Series:
        """
        Get the last value of the sublists.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [], [1, 2]])
        >>> s.list.last()
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            null
            2
        ]
        """

    @unstable()
    def item(self) -> Series:
        """
        Get the single value of the sublists.

        This errors if the sublist length is not exactly one.

        See Also
        --------
        :meth:`Series.list.get` : Get the value by index in the sublists.

        Examples
        --------
        >>> s = pl.Series("a", [[1], [4], [6]])
        >>> s.list.item()
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            4
            6
        ]
        >>> df = pl.Series("a", [[3, 2, 1], [1], [2]])
        >>> df.list.item()
        Traceback (most recent call last):
        ...
        polars.exceptions.ComputeError: aggregation 'item' expected a single value, got 3 values
        """  # noqa: W505

    def contains(self, item: IntoExpr, *, nulls_equal: bool = True) -> Series:
        """
        Check if sublists contain the given item.

        Parameters
        ----------
        item
            Item that will be checked for membership
        nulls_equal : bool, default True
            If True, treat null as a distinct value. Null values will not propagate.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [], [1, 2]])
        >>> s.list.contains(1)
        shape: (3,)
        Series: 'a' [bool]
        [
            true
            false
            true
        ]
        """

    def arg_min(self) -> Series:
        """
        Retrieve the index of the minimal value in every sublist.

        Returns
        -------
        Series
            Series of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [2, 1]])
        >>> s.list.arg_min()
        shape: (2,)
        Series: 'a' [u32]
        [
            0
            1
        ]
        """

    def arg_max(self) -> Series:
        """
        Retrieve the index of the maximum value in every sublist.

        Returns
        -------
        Series
            Series of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [2, 1]])
        >>> s.list.arg_max()
        shape: (2,)
        Series: 'a' [u32]
        [
            1
            0
        ]
        """

    def diff(self, n: int = 1, null_behavior: NullBehavior = "ignore") -> Series:
        """
        Calculate the first discrete difference between shifted items of every sublist.

        Parameters
        ----------
        n
            Number of slots to shift.
        null_behavior : {'ignore', 'drop'}
            How to handle null values.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3, 4], [10, 2, 1]])
        >>> s.list.diff()
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [null, 1, … 1]
            [null, -8, -1]
        ]

        >>> s.list.diff(n=2)
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [null, null, … 2]
            [null, null, -9]
        ]

        >>> s.list.diff(n=2, null_behavior="drop")
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [2, 2]
            [-9]
        ]
        """

    def shift(self, n: int | IntoExprColumn = 1) -> Series:
        """
        Shift list values by the given number of indices.

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
        By default, list values are shifted forward by one index.

        >>> s = pl.Series([[1, 2, 3], [4, 5]])
        >>> s.list.shift()
        shape: (2,)
        Series: '' [list[i64]]
        [
                [null, 1, 2]
                [null, 4]
        ]

        Pass a negative value to shift in the opposite direction instead.

        >>> s.list.shift(-2)
        shape: (2,)
        Series: '' [list[i64]]
        [
                [3, null, null]
                [null, null]
        ]
        """

    def slice(self, offset: int | Expr, length: int | Expr | None = None) -> Series:
        """
        Slice every sublist.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None` (default), the slice is taken to the
            end of the list.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3, 4], [10, 2, 1]])
        >>> s.list.slice(1, 2)
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [2, 3]
            [2, 1]
        ]
        """

    def head(self, n: int | Expr = 5) -> Series:
        """
        Slice the first `n` values of every sublist.

        Parameters
        ----------
        n
            Number of values to return for each sublist.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3, 4], [10, 2, 1]])
        >>> s.list.head(2)
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [1, 2]
            [10, 2]
        ]
        """

    def tail(self, n: int | Expr = 5) -> Series:
        """
        Slice the last `n` values of every sublist.

        Parameters
        ----------
        n
            Number of values to return for each sublist.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3, 4], [10, 2, 1]])
        >>> s.list.tail(2)
        shape: (2,)
        Series: 'a' [list[i64]]
        [
            [3, 4]
            [2, 1]
        ]
        """

    def explode(self) -> Series:
        """
        Returns a column with a separate row for every list element.

        Returns
        -------
        Series
            Series with the data type of the list elements.

        See Also
        --------
        Series.reshape : Reshape this Series to a flat Series or a Series of Lists.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3], [4, 5, 6]])
        >>> s.list.explode()
        shape: (6,)
        Series: 'a' [i64]
        [
            1
            2
            3
            4
            5
            6
        ]
        """

    def count_matches(self, element: IntoExpr) -> Series:
        """
        Count how often the value produced by `element` occurs.

        Parameters
        ----------
        element
            An expression that produces a single value

        Examples
        --------
        >>> s = pl.Series("a", [[0], [1], [1, 2, 3, 2], [1, 2, 1], [4, 4]])
        >>> s.list.count_matches(1)
        shape: (5,)
        Series: 'a' [u32]
        [
            0
            1
            1
            2
            0
        ]
        """

    def to_array(self, width: int) -> Series:
        """
        Convert a List column into an Array column with the same inner data type.

        Parameters
        ----------
        width
            Width of the resulting Array column.

        Returns
        -------
        Series
            Series of data type :class:`Array`.

        Examples
        --------
        >>> s = pl.Series([[1, 2], [3, 4]], dtype=pl.List(pl.Int8))
        >>> s.list.to_array(2)
        shape: (2,)
        Series: '' [array[i8, 2]]
        [
                [1, 2]
                [3, 4]
        ]
        """

    def to_struct(
        self,
        n_field_strategy: ListToStructWidthStrategy = "first_non_null",
        fields: Callable[[int], str] | Sequence[str] | None = None,
    ) -> Series:
        """
        Convert the series of type `List` to a series of type `Struct`.

        Parameters
        ----------
        n_field_strategy : {'first_non_null', 'max_width'}
            Strategy to determine the number of fields of the struct.

            * "first_non_null": set number of fields equal to the length of the
              first non zero-length sublist.
            * "max_width": set number of fields as max length of all sublists.
        fields
            If the name and number of the desired fields is known in advance
            a list of field names can be given, which will be assigned by index.
            Otherwise, to dynamically assign field names, a custom function can be
            used; if neither are set, fields will be `field_0, field_1 .. field_n`.

        Examples
        --------
        Convert list to struct with default field name assignment:

        >>> s1 = pl.Series("n", [[0, 1, 2], [0, 1]])
        >>> s2 = s1.list.to_struct()
        >>> s2
        shape: (2,)
        Series: 'n' [struct[3]]
        [
            {0,1,2}
            {0,1,null}
        ]
        >>> s2.struct.fields
        ['field_0', 'field_1', 'field_2']

        Convert list to struct with field name assignment by function/index:

        >>> s3 = s1.list.to_struct(fields=lambda idx: f"n{idx:02}")
        >>> s3.struct.fields
        ['n00', 'n01', 'n02']

        Convert list to struct with field name assignment by index from a list of names:

        >>> s1.list.to_struct(fields=["one", "two", "three"]).struct.unnest()
        shape: (2, 3)
        ┌─────┬─────┬───────┐
        │ one ┆ two ┆ three │
        │ --- ┆ --- ┆ ---   │
        │ i64 ┆ i64 ┆ i64   │
        ╞═════╪═════╪═══════╡
        │ 0   ┆ 1   ┆ 2     │
        │ 0   ┆ 1   ┆ null  │
        └─────┴─────┴───────┘
        """
        if isinstance(fields, Sequence):
            s = wrap_s(self._s)
            return (
                s.to_frame()
                .select_seq(F.col(s.name).list.to_struct(fields=fields))
                .to_series()
            )

        return wrap_s(self._s.list_to_struct(n_field_strategy, fields))

    def eval(self, expr: Expr, *, parallel: bool = False) -> Series:
        """
        Run any polars expression against the lists' elements.

        Parameters
        ----------
        expr
            Expression to run. Note that you can select an element with `pl.first()`, or
            `pl.col()`
        parallel
            Run all expression parallel. Don't activate this blindly.
            Parallelism is worth it if there is enough work to do per thread.

            This likely should not be use in the group by context, because we already
            parallel execution per group

        Examples
        --------
        >>> s = pl.Series("a", [[1, 4], [8, 5], [3, 2]])
        >>> s.list.eval(pl.element().rank())
        shape: (3,)
        Series: 'a' [list[f64]]
        [
            [1.0, 2.0]
            [2.0, 1.0]
            [2.0, 1.0]
        ]
        """

    def agg(self, expr: Expr) -> Series:
        """

        Run any polars aggregation expression against the list' elements.

        Parameters
        ----------
        expr
            Expression to run. Note that you can select an element with `pl.element()`.

        Examples
        --------
        >>> s = pl.Series("a", [[1, None], [42, 13], [None, None]])
        >>> s.list.agg(pl.element().null_count())
        shape: (3,)
        Series: 'a' [u32]
        [
            1
            0
            2
        ]
        >>> s.list.agg(pl.element().drop_nulls())
        shape: (3,)
        Series: 'a' [list[i64]]
        [
            [1]
            [42, 13]
            []
        ]
        """

    def filter(self, predicate: Expr) -> Series:
        """
        Filter elements in each list by a boolean expression, returning a new Series of lists.

        Parameters
        ----------
        predicate
            A boolean expression evaluated on each list element.
            Use `pl.element()` to refer to the current element.

        Examples
        --------
        >>> import polars as pl
        >>> s = pl.Series("a", [[1, 4], [8, 5], [3, 2]])
        >>> s.list.filter(pl.element() % 2 == 0)
        shape: (3,)
        Series: 'a' [list[i64]]
        [
            [4]
            [8]
            [2]
        ]
        """  # noqa: W505

    def set_union(self, other: Series | Collection[Any]) -> Series:
        """
        Compute the SET UNION between the elements in this list and the elements of `other`.

        Parameters
        ----------
        other
            Right hand side of the set operation.

        Examples
        --------
        >>> a = pl.Series([[1, 2, 3], [], [None, 3], [5, 6, 7]])
        >>> b = pl.Series([[2, 3, 4], [3], [3, 4, None], [6, 8]])
        >>> a.list.set_union(b)  # doctest: +IGNORE_RESULT
        shape: (4,)
        Series: '' [list[i64]]
        [
                [1, 2, 3, 4]
                [3]
                [null, 3, 4]
                [5, 6, 7, 8]
        ]
        """  # noqa: W505

    def set_difference(self, other: Series | Collection[Any]) -> Series:
        """
        Compute the SET DIFFERENCE between the elements in this list and the elements of `other`.

        Parameters
        ----------
        other
            Right hand side of the set operation.

        See Also
        --------
        polars.Series.list.diff: Calculates the n-th discrete difference of every sublist.

        Examples
        --------
        >>> a = pl.Series([[1, 2, 3], [], [None, 3], [5, 6, 7]])
        >>> b = pl.Series([[2, 3, 4], [3], [3, 4, None], [6, 8]])
        >>> a.list.set_difference(b)
        shape: (4,)
        Series: '' [list[i64]]
        [
                [1]
                []
                []
                [5, 7]
        ]
        """  # noqa: W505

    def set_intersection(self, other: Series | Collection[Any]) -> Series:
        """
        Compute the SET INTERSECTION between the elements in this list and the elements of `other`.

        Parameters
        ----------
        other
            Right hand side of the set operation.

        Examples
        --------
        >>> a = pl.Series([[1, 2, 3], [], [None, 3], [5, 6, 7]])
        >>> b = pl.Series([[2, 3, 4], [3], [3, 4, None], [6, 8]])
        >>> a.list.set_intersection(b)
        shape: (4,)
        Series: '' [list[i64]]
        [
                [2, 3]
                []
                [null, 3]
                [6]
        ]
        """  # noqa: W505

    def set_symmetric_difference(self, other: Series | Collection[Any]) -> Series:
        """
        Compute the SET SYMMETRIC DIFFERENCE between the elements in this list and the elements of `other`.

        Parameters
        ----------
        other
            Right hand side of the set operation.

        Examples
        --------
        >>> a = pl.Series([[1, 2, 3], [], [None, 3], [5, 6, 7]])
        >>> b = pl.Series([[2, 3, 4], [3], [3, 4, None], [6, 8]])
        >>> a.list.set_symmetric_difference(b)
        shape: (4,)
        Series: '' [list[i64]]
        [
            [1, 4]
            [3]
            [4]
            [5, 7, 8]
        ]
        """  # noqa: W505
