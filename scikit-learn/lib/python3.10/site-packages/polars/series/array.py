from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from polars import functions as F
from polars._utils.wrap import wrap_s
from polars.series.utils import expr_dispatch

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars import Series
    from polars._plr import PySeries
    from polars._typing import IntoExpr, IntoExprColumn
    from polars.expr.expr import Expr


@expr_dispatch
class ArrayNameSpace:
    """Namespace for array related methods."""

    _accessor = "arr"

    def __init__(self, series: Series) -> None:
        self._s: PySeries = series._s

    def min(self) -> Series:
        """
        Compute the min values of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.min()
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            3
        ]
        """

    def max(self) -> Series:
        """
        Compute the max values of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.max()
        shape: (2,)
        Series: 'a' [i64]
        [
            2
            4
        ]
        """

    def sum(self) -> Series:
        """
        Compute the sum values of the sub-arrays.

        Notes
        -----
        If there are no non-null elements in a row, the output is `0`.

        Examples
        --------
        >>> s = pl.Series([[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.sum()
        shape: (2,)
        Series: '' [i64]
        [
            3
            7
        ]
        """

    def std(self, ddof: int = 1) -> Series:
        """
        Compute the std of the values of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.std()
        shape: (2,)
        Series: 'a' [f64]
        [
            0.707107
            0.707107
        ]
        """

    def var(self, ddof: int = 1) -> Series:
        """
        Compute the var of the values of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.var()
        shape: (2,)
        Series: 'a' [f64]
        [
                0.5
                0.5
        ]
        """

    def median(self) -> Series:
        """
        Compute the median of the values of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.median()
        shape: (2,)
        Series: 'a' [f64]
        [
            1.5
            3.5
        ]
        """

    def unique(self, *, maintain_order: bool = False) -> Series:
        """
        Get the unique/distinct values in the array.

        Parameters
        ----------
        maintain_order
            Maintain order of data. This requires more work.

        Returns
        -------
        Series
            Series of data type :class:`List`.

        Examples
        --------
        >>> s = pl.Series([[1, 1, 2], [3, 4, 5]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.unique()
        shape: (2,)
        Series: '' [list[i64]]
        [
            [1, 2]
            [3, 4, 5]
        ]
        """

    def n_unique(self) -> Series:
        """
        Count the number of unique values in every sub-arrays.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 4]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.n_unique()
        shape: (2,)
        Series: 'a' [u32]
        [
            2
            1
        ]
        """

    def to_list(self) -> Series:
        """
        Convert an Array column into a List column with the same inner data type.

        Returns
        -------
        Series
            Series of data type :class:`List`.

        Examples
        --------
        >>> s = pl.Series([[1, 2], [3, 4]], dtype=pl.Array(pl.Int8, 2))
        >>> s.arr.to_list()
        shape: (2,)
        Series: '' [list[i8]]
        [
                [1, 2]
                [3, 4]
        ]
        """

    def any(self) -> Series:
        """
        Evaluate whether any boolean value is true for every subarray.

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
        ...     [[True, True], [False, True], [False, False], [None, None], None],
        ...     dtype=pl.Array(pl.Boolean, 2),
        ... )
        >>> s.arr.any()
        shape: (5,)
        Series: '' [bool]
        [
            true
            true
            false
            false
            null
        ]
        """

    def len(self) -> Series:
        """
        Return the number of elements in each array.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
        >>> s.arr.len()
        shape: (2,)
        Series: 'a' [u32]
        [
            2
            2
        ]
        """

    def slice(
        self,
        offset: int | Expr,
        length: int | Expr | None = None,
        *,
        as_array: bool = False,
    ) -> Series:
        """
        Slice the sub-arrays.

        Parameters
        ----------
        offset
            The starting index of the slice.
        length
            The length of the slice.
        as_array
            Return the result as a Series of data type :class:`.Array`.

        Returns
        -------
        Series
            Series of data type :class:`.List` or :class:`.Array` if `as_array=True`.

        Examples
        --------
        >>> s = pl.Series(
        ...     [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]],
        ...     dtype=pl.Array(pl.Int64, 6),
        ... )
        >>> s.arr.slice(1)
        shape: (2,)
        Series: '' [list[i64]]
        [
            [2, 3, … 6]
            [8, 9, … 12]
        ]
        >>> s.arr.slice(1, 3, as_array=True)
        shape: (2,)
        Series: '' [array[i64, 3]]
        [
            [2, 3, 4]
            [8, 9, 10]
        ]
        >>> s.arr.slice(-2)
        shape: (2,)
        Series: '' [list[i64]]
        [
            [5, 6]
            [11, 12]
        ]
        """

    def head(self, n: int | Expr = 5, *, as_array: bool = False) -> Series:
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
        >>> s = pl.Series(
        ...     [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]],
        ...     dtype=pl.Array(pl.Int64, 6),
        ... )
        >>> s.arr.head()
        shape: (2,)
        Series: '' [list[i64]]
        [
            [1, 2, … 5]
            [7, 8, … 11]
        ]
        >>> s.arr.head(3, as_array=True)
        shape: (2,)
        Series: '' [array[i64, 3]]
        [
            [1, 2, 3]
            [7, 8, 9]
        ]
        """

    def tail(self, n: int | Expr = 5, *, as_array: bool = False) -> Series:
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
        >>> s = pl.Series(
        ...     [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]],
        ...     dtype=pl.Array(pl.Int64, 6),
        ... )
        >>> s.arr.tail()
        shape: (2,)
        Series: '' [list[i64]]
        [
            [2, 3, … 6]
            [8, 9, … 12]
        ]
        >>> s.arr.tail(3, as_array=True)
        shape: (2,)
        Series: '' [array[i64, 3]]
        [
            [4, 5, 6]
            [10, 11, 12]
        ]
        """

    def all(self) -> Series:
        """
        Evaluate whether all boolean values are true for every subarray.

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
        ...     [[True, True], [False, True], [False, False], [None, None], None],
        ...     dtype=pl.Array(pl.Boolean, 2),
        ... )
        >>> s.arr.all()
        shape: (5,)
        Series: '' [bool]
        [
            true
            false
            false
            true
            null
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
        >>> s = pl.Series("a", [[3, 2, 1], [9, 1, 2]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.sort()
        shape: (2,)
        Series: 'a' [array[i64, 3]]
        [
            [1, 2, 3]
            [1, 2, 9]
        ]
        >>> s.arr.sort(descending=True)
        shape: (2,)
        Series: 'a' [array[i64, 3]]
        [
            [3, 2, 1]
            [9, 2, 1]
        ]

        """

    def reverse(self) -> Series:
        """
        Reverse the arrays in this column.

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [9, 1, 2]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.reverse()
        shape: (2,)
        Series: 'a' [array[i64, 3]]
        [
            [1, 2, 3]
            [2, 1, 9]
        ]

        """

    def arg_min(self) -> Series:
        """
        Retrieve the index of the minimal value in every sub-array.

        Returns
        -------
        Series
            Series of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> s = pl.Series("a", [[3, 2, 1], [9, 1, 2]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.arg_min()
        shape: (2,)
        Series: 'a' [u32]
        [
            2
            1
        ]

        """

    def arg_max(self) -> Series:
        """
        Retrieve the index of the maximum value in every sub-array.

        Returns
        -------
        Series
            Series of data type :class:`UInt32` or :class:`UInt64`
            (depending on compilation).

        Examples
        --------
        >>> s = pl.Series("a", [[0, 9, 3], [9, 1, 2]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.arg_max()
        shape: (2,)
        Series: 'a' [u32]
        [
            1
            0
        ]

        """

    def get(self, index: int | IntoExprColumn, *, null_on_oob: bool = False) -> Series:
        """
        Get the value by index in the sub-arrays.

        So index `0` would return the first item of every sublist
        and index `-1` would return the last item of every sublist
        if an index is out of bounds, it will return a `None`.

        Parameters
        ----------
        index
            Index to return per sublist
        null_on_oob
            Behavior if an index is out of bounds:
            True -> set as null
            False -> raise an error

        Returns
        -------
        Series
            Series of innter data type.

        Examples
        --------
        >>> s = pl.Series(
        ...     "a", [[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=pl.Array(pl.Int32, 3)
        ... )
        >>> s.arr.get(pl.Series([1, -2, 0]), null_on_oob=True)
        shape: (3,)
        Series: 'a' [i32]
        [
            2
            5
            7
        ]

        """

    def first(self) -> Series:
        """
        Get the first value of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series(
        ...     "a", [[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=pl.Array(pl.Int32, 3)
        ... )
        >>> s.arr.first()
        shape: (3,)
        Series: 'a' [i32]
        [
            1
            4
            7
        ]

        """

    def last(self) -> Series:
        """
        Get the last value of the sub-arrays.

        Examples
        --------
        >>> s = pl.Series(
        ...     "a", [[1, 2, 3], [4, 5, 6], [7, 9, 8]], dtype=pl.Array(pl.Int32, 3)
        ... )
        >>> s.arr.last()
        shape: (3,)
        Series: 'a' [i32]
        [
            3
            6
            8
        ]

        """

    def join(self, separator: IntoExprColumn, *, ignore_nulls: bool = True) -> Series:
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
        Series
            Series of data type :class:`String`.

        Examples
        --------
        >>> s = pl.Series([["x", "y"], ["a", "b"]], dtype=pl.Array(pl.String, 2))
        >>> s.arr.join(separator="-")
        shape: (2,)
        Series: '' [str]
        [
            "x-y"
            "a-b"
        ]

        """

    def explode(self) -> Series:
        """
        Returns a column with a separate row for every array element.

        Returns
        -------
        Series
            Series with the data type of the array elements.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3], [4, 5, 6]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.explode()
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

    def contains(self, item: IntoExpr, *, nulls_equal: bool = True) -> Series:
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
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series(
        ...     "a", [[3, 2, 1], [1, 2, 3], [4, 5, 6]], dtype=pl.Array(pl.Int32, 3)
        ... )
        >>> s.arr.contains(1)
        shape: (3,)
        Series: 'a' [bool]
        [
            true
            true
            false
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
        >>> s = pl.Series("a", [[1, 2, 3], [2, 2, 2]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.count_matches(2)
        shape: (2,)
        Series: 'a' [u32]
        [
            1
            3
        ]

        """

    def to_struct(
        self,
        fields: Callable[[int], str] | Sequence[str] | None = None,
    ) -> Series:
        """
        Convert the series of type `Array` to a series of type `Struct`.

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

        >>> s1 = pl.Series("n", [[0, 1, 2], [3, 4, 5]], dtype=pl.Array(pl.Int8, 3))
        >>> s2 = s1.arr.to_struct()
        >>> s2
        shape: (2,)
        Series: 'n' [struct[3]]
        [
            {0,1,2}
            {3,4,5}
        ]
        >>> s2.struct.fields
        ['field_0', 'field_1', 'field_2']

        Convert array to struct with field name assignment by function/index:

        >>> s3 = s1.arr.to_struct(fields=lambda idx: f"n{idx:02}")
        >>> s3.struct.fields
        ['n00', 'n01', 'n02']

        Convert array to struct with field name assignment by
        index from a list of names:

        >>> s1.arr.to_struct(fields=["one", "two", "three"]).struct.unnest()
        shape: (2, 3)
        ┌─────┬─────┬───────┐
        │ one ┆ two ┆ three │
        │ --- ┆ --- ┆ ---   │
        │ i8  ┆ i8  ┆ i8    │
        ╞═════╪═════╪═══════╡
        │ 0   ┆ 1   ┆ 2     │
        │ 3   ┆ 4   ┆ 5     │
        └─────┴─────┴───────┘
        """
        s = wrap_s(self._s)
        return s.to_frame().select(F.col(s.name).arr.to_struct(fields)).to_series()

    def shift(self, n: int | IntoExprColumn = 1) -> Series:
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

        >>> s = pl.Series([[1, 2, 3], [4, 5, 6]], dtype=pl.Array(pl.Int64, 3))
        >>> s.arr.shift()
        shape: (2,)
        Series: '' [array[i64, 3]]
        [
            [null, 1, 2]
            [null, 4, 5]
        ]

        Pass a negative value to shift in the opposite direction instead.

        >>> s.arr.shift(-2)
        shape: (2,)
        Series: '' [array[i64, 3]]
        [
            [3, null, null]
            [6, null, null]
        ]
        """

    def eval(self, expr: Expr, *, as_list: bool = False) -> Series:
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
        >>> s = pl.Series("a", [[1, 4], [8, 5], [3, 2]], pl.Array(pl.Int64, 2))
        >>> s.arr.eval(pl.element().rank())
        shape: (3,)
        Series: 'a' [array[f64, 2]]
        [
            [1.0, 2.0]
            [2.0, 1.0]
            [2.0, 1.0]
        ]
        """

    def agg(self, expr: Expr) -> Series:
        """
        Run any polars aggregation expression against the arrays' elements.

        Parameters
        ----------
        expr
            Expression to run. Note that you can select an element with `pl.element()`.

        Examples
        --------
        >>> s = pl.Series(
        ...     "a", [[1, None], [42, 13], [None, None]], pl.Array(pl.Int64, 2)
        ... )
        >>> s.arr.agg(pl.element().null_count())
        shape: (3,)
        Series: 'a' [u32]
        [
            1
            0
            2
        ]
        >>> s.arr.agg(pl.element().drop_nulls())
        shape: (3,)
        Series: 'a' [list[i64]]
        [
            [1]
            [42, 13]
            []
        ]
        """
