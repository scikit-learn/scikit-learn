from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.wrap import wrap_df

if TYPE_CHECKING:
    from polars import DataFrame
    from polars._plr import PyInProcessQuery


class InProcessQuery:
    """
    A placeholder for an in process query.

    This can be used to do something else while a query is running.
    The queries can be cancelled. You can peek if the query is finished,
    or you can await the result.
    """

    def __init__(self, ipq: PyInProcessQuery) -> None:
        self._inner = ipq

    def cancel(self) -> None:
        """Cancel the query at earliest convenience."""
        self._inner.cancel()

    def fetch(self) -> DataFrame | None:
        """
        Fetch the result.

        If it is ready, a materialized DataFrame is returned.
        If it is not ready it will return `None`.
        """
        if (out := self._inner.fetch()) is not None:
            return wrap_df(out)
        else:
            return None

    def fetch_blocking(self) -> DataFrame:
        """Await the result synchronously."""
        return wrap_df(self._inner.fetch_blocking())
