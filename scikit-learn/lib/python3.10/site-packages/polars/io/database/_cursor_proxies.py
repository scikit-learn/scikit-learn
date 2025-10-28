from __future__ import annotations

from typing import TYPE_CHECKING, Any

from polars._dependencies import import_optional
from polars.io.database._utils import _run_async

if TYPE_CHECKING:
    import sys
    from collections.abc import Coroutine, Iterable

    import pyarrow as pa

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self


class ODBCCursorProxy:
    """Cursor proxy for ODBC connections (requires `arrow-odbc`)."""

    def __init__(self, connection_string: str) -> None:
        self.connection_string = connection_string
        self.execute_options: dict[str, Any] = {}
        self.query: str | None = None

    def close(self) -> None:
        """Close the cursor."""
        # n/a: nothing to close

    def execute(self, query: str, **execute_options: Any) -> None:
        """Execute a query (n/a: just store query for the fetch* methods)."""
        self.execute_options = execute_options
        self.query = query

    def fetch_arrow_table(
        self, batch_size: int = 10_000, *, fetch_all: bool = False
    ) -> pa.Table:
        """Fetch all results as a pyarrow Table."""
        from pyarrow import Table

        return Table.from_batches(
            self.fetch_record_batches(batch_size=batch_size, fetch_all=True)
        )

    def fetch_record_batches(
        self, batch_size: int = 10_000, *, fetch_all: bool = False
    ) -> Iterable[pa.RecordBatch]:
        """Fetch results as an iterable of RecordBatches."""
        from arrow_odbc import read_arrow_batches_from_odbc
        from pyarrow import RecordBatch

        n_batches = 0
        batch_reader = read_arrow_batches_from_odbc(
            query=self.query,
            batch_size=batch_size,
            connection_string=self.connection_string,
            **self.execute_options,
        )
        for batch in batch_reader:
            yield batch
            n_batches += 1

        if n_batches == 0 and fetch_all:
            # empty result set; return empty batch with accurate schema
            yield RecordBatch.from_pylist([], schema=batch_reader.schema)

    # note: internally arrow-odbc always reads batches
    fetchall = fetch_arrow_table
    fetchmany = fetch_record_batches


class SurrealDBCursorProxy:
    """Cursor proxy for both SurrealDB and AsyncSurrealDB connections."""

    _cached_result: list[dict[str, Any]] | None = None

    def __init__(self, client: Any) -> None:
        surrealdb = import_optional("surrealdb")
        self.is_async = isinstance(client, surrealdb.AsyncSurrealDB)
        self.execute_options: dict[str, Any] = {}
        self.client = client
        self.query: str = None  # type: ignore[assignment]

    @staticmethod
    async def _unpack_result_async(
        result: Coroutine[Any, Any, list[dict[str, Any]]],
    ) -> Coroutine[Any, Any, list[dict[str, Any]]]:
        """Unpack the async query result."""
        response = (await result)[0]
        if response["status"] != "OK":
            raise RuntimeError(response["result"])
        return response["result"]

    @staticmethod
    def _unpack_result(
        result: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        """Unpack the query result."""
        response = result[0]
        if response["status"] != "OK":
            raise RuntimeError(response["result"])
        return response["result"]

    def close(self) -> None:
        """Close the cursor."""
        # no-op; never close a user's Surreal session

    def execute(self, query: str, **execute_options: Any) -> Self:
        """Execute a query (n/a: just store query for the fetch* methods)."""
        self._cached_result = None
        self.execute_options = execute_options
        self.query = query
        return self

    def fetchall(self) -> list[dict[str, Any]]:
        """Fetch all results (as a list of dictionaries)."""
        return (
            _run_async(
                self._unpack_result_async(
                    result=self.client.query(
                        query=self.query,
                        variables=(self.execute_options or None),
                    ),
                )
            )
            if self.is_async
            else self._unpack_result(
                result=self.client.query(
                    query=self.query,
                    variables=(self.execute_options or None),
                ),
            )
        )

    def fetchmany(self, size: int) -> list[dict[str, Any]]:
        """Fetch results in batches (simulated)."""
        # first 'fetchmany' call acquires/caches the result object
        if self._cached_result is None:
            self._cached_result = self.fetchall()

        # return batches from the result, actively removing from the cache
        # as we go, so as not to hold on to additional copies when done
        result = self._cached_result[:size]
        del self._cached_result[:size]
        return result
