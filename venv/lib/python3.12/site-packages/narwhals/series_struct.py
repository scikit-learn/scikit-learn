from __future__ import annotations

from typing import Generic

from narwhals.typing import SeriesT


class SeriesStructNamespace(Generic[SeriesT]):
    def __init__(self, series: SeriesT) -> None:
        self._narwhals_series = series

    def field(self, name: str) -> SeriesT:
        r"""Retrieve a Struct field as a new expression.

        Arguments:
            name: Name of the struct field to retrieve.

        Returns:
            A new Series.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> s_native = pl.Series(
            ...     [{"id": "0", "name": "john"}, {"id": "1", "name": "jane"}]
            ... )
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.struct.field("name").to_list()
            ['john', 'jane']
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.struct.field(name)
        )
