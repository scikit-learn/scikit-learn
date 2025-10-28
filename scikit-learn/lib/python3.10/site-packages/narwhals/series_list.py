from __future__ import annotations

from typing import TYPE_CHECKING, Generic

from narwhals.typing import SeriesT

if TYPE_CHECKING:
    from narwhals.typing import NonNestedLiteral


class SeriesListNamespace(Generic[SeriesT]):
    def __init__(self, series: SeriesT) -> None:
        self._narwhals_series = series

    def len(self) -> SeriesT:
        """Return the number of elements in each list.

        Null values count towards the total.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> s_native = pa.chunked_array([[[1, 2], [3, 4, None], None, []]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.list.len().to_native()  # doctest: +ELLIPSIS
            <pyarrow.lib.ChunkedArray object at ...>
            [
              [
                2,
                3,
                null,
                0
              ]
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.list.len()
        )

    def unique(self) -> SeriesT:
        """Get the unique/distinct values in the list.

        Null values are included in the result. The order of unique values is not guaranteed.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([[1, 1, 2], [3, 3, None], None, []])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.list.unique().to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (4,)
            Series: '' [list[i64]]
            [
               [1, 2]
               [null, 3]
               null
               []
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.list.unique()
        )

    def contains(self, item: NonNestedLiteral) -> SeriesT:
        """Check if sublists contain the given item.

        Arguments:
            item: Item that will be checked for membership.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([[1, 2], None, []])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.list.contains(1).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [bool]
            [
                    true
                    null
                    false
            ]
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.list.contains(item)
        )

    def get(self, index: int) -> SeriesT:
        """Return the value by index in each list.

        Negative indices are not accepted.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series([[1, 2], [3, 4, None], [None, 5]])
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.list.get(1).to_native()  # doctest: +NORMALIZE_WHITESPACE
            shape: (3,)
            Series: '' [i64]
            [
                    2
                    4
                    5
            ]
        """
        if not isinstance(index, int):
            msg = (
                f"Index must be of type 'int'. Got type '{type(index).__name__}' instead."
            )
            raise TypeError(msg)

        if index < 0:
            msg = f"Index {index} is out of bounds: should be greater than or equal to 0."
            raise ValueError(msg)

        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.list.get(index)
        )
