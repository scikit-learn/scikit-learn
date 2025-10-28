from __future__ import annotations

from typing import Generic

from narwhals.typing import SeriesT


class SeriesCatNamespace(Generic[SeriesT]):
    def __init__(self, series: SeriesT) -> None:
        self._narwhals_series = series

    def get_categories(self) -> SeriesT:
        """Get unique categories from column.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> s_native = pd.Series(["apple", "mango", "mango"], dtype="category")
            >>> s = nw.from_native(s_native, series_only=True)
            >>> s.cat.get_categories().to_native()
            0    apple
            1    mango
            dtype: object
        """
        return self._narwhals_series._with_compliant(
            self._narwhals_series._compliant_series.cat.get_categories()
        )
