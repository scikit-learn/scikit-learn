from __future__ import annotations

import re
from collections.abc import Coroutine, Sequence
from contextlib import suppress
from inspect import Parameter, signature
from typing import TYPE_CHECKING, Any

from polars import functions as F
from polars._utils.various import parse_version, qualified_type_name
from polars.convert import from_arrow
from polars.datatypes import N_INFER_DEFAULT
from polars.exceptions import (
    DuplicateError,
    ModuleUpgradeRequiredError,
    UnsuitableSQLError,
)
from polars.io.database._arrow_registry import ARROW_DRIVER_REGISTRY
from polars.io.database._cursor_proxies import ODBCCursorProxy, SurrealDBCursorProxy
from polars.io.database._inference import dtype_from_cursor_description
from polars.io.database._utils import _run_async

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterable, Iterator
    from types import TracebackType

    import pyarrow as pa

    from polars.io.database._arrow_registry import ArrowDriverProperties

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

    from sqlalchemy.sql.elements import TextClause
    from sqlalchemy.sql.expression import Selectable

    from polars import DataFrame
    from polars._typing import ConnectionOrCursor, Cursor, SchemaDict

_INVALID_QUERY_TYPES = {
    "ALTER",
    "ANALYZE",
    "CREATE",
    "DELETE",
    "DROP",
    "GRANT",
    "INSERT",
    "REPLACE",
    "REVOKE",
    "UPDATE",
    "UPSERT",
    "USE",
    "VACUUM",
}


class CloseAfterFrameIter:
    """Allows cursor close to be deferred until the last batch is returned."""

    def __init__(self, frames: Any, *, cursor: Cursor) -> None:
        self._iter_frames = frames
        self._cursor = cursor

    def __iter__(self) -> Iterator[DataFrame]:
        yield from self._iter_frames

        if hasattr(self._cursor, "close"):
            self._cursor.close()


class ConnectionExecutor:
    """Abstraction for querying databases with user-supplied connection objects."""

    # indicate if we can/should close the cursor on scope exit. note that we
    # should never close the underlying connection, or a user-supplied cursor.
    can_close_cursor: bool = False

    def __init__(self, connection: ConnectionOrCursor) -> None:
        self.driver_name = (
            "arrow_odbc_proxy"
            if isinstance(connection, ODBCCursorProxy)
            else type(connection).__module__.split(".", 1)[0].lower()
        )
        if self.driver_name == "surrealdb":
            connection = SurrealDBCursorProxy(client=connection)

        self.cursor = self._normalise_cursor(connection)
        self.result: Any = None

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        # if we created it and are finished with it, we can
        # close the cursor (but NOT the connection)
        if self._is_alchemy_async(self.cursor):
            from sqlalchemy.ext.asyncio import AsyncConnection

            if isinstance(self.cursor, AsyncConnection):
                _run_async(self._close_async_cursor())
        elif self.can_close_cursor and hasattr(self.cursor, "close"):
            self.cursor.close()

    def __repr__(self) -> str:
        return f"<{type(self).__name__} module={self.driver_name!r}>"

    @staticmethod
    def _apply_overrides(df: DataFrame, schema_overrides: SchemaDict) -> DataFrame:
        """Apply schema overrides to a DataFrame."""
        existing_schema = df.schema
        if cast_cols := [
            F.col(col).cast(dtype)
            for col, dtype in schema_overrides.items()
            if col in existing_schema and dtype != existing_schema[col]
        ]:
            df = df.with_columns(cast_cols)
        return df

    async def _close_async_cursor(self) -> None:
        if self.can_close_cursor and hasattr(self.cursor, "close"):
            from sqlalchemy.ext.asyncio.exc import AsyncContextNotStarted

            with suppress(AsyncContextNotStarted):
                await self.cursor.close()

    @staticmethod
    def _check_module_version(module_name: str, minimum_version: str) -> None:
        """Check the module version against a minimum required version."""
        mod = __import__(module_name)
        with suppress(AttributeError):
            module_version: tuple[int, ...] | None = None
            for version_attr in ("__version__", "version"):
                if isinstance(ver := getattr(mod, version_attr, None), str):
                    module_version = parse_version(ver)
                    break
            if module_version and module_version < parse_version(minimum_version):
                msg = f"`read_database` queries require at least {module_name} version {minimum_version}"
                raise ModuleUpgradeRequiredError(msg)

    def _fetch_arrow(
        self,
        driver_properties: ArrowDriverProperties,
        *,
        batch_size: int | None,
        iter_batches: bool,
    ) -> Iterable[pa.RecordBatch]:
        """Yield Arrow data as a generator of one or more RecordBatches or Tables."""
        fetch_batches = driver_properties["fetch_batches"]
        if not iter_batches or fetch_batches is None:
            fetch_method = driver_properties["fetch_all"]
            yield getattr(self.result, fetch_method)()
        else:
            size = [batch_size] if driver_properties["exact_batch_size"] else []
            repeat_batch_calls = driver_properties["repeat_batch_calls"]
            fetchmany_arrow = getattr(self.result, fetch_batches)
            if not repeat_batch_calls:
                yield from fetchmany_arrow(*size)
            else:
                while True:
                    arrow = fetchmany_arrow(*size)
                    if not arrow:
                        break
                    yield arrow

    @staticmethod
    def _fetchall_rows(result: Cursor, *, is_alchemy: bool) -> Iterable[Sequence[Any]]:
        """Fetch row data in a single call, returning the complete result set."""
        rows = result.fetchall()
        return (
            rows
            if rows and (is_alchemy or isinstance(rows[0], (list, tuple, dict)))
            else [tuple(row) for row in rows]
        )

    def _fetchmany_rows(
        self, result: Cursor, *, batch_size: int | None, is_alchemy: bool
    ) -> Iterable[Sequence[Any]]:
        """Fetch row data incrementally, yielding over the complete result set."""
        while True:
            rows = result.fetchmany(batch_size)
            if not rows:
                break
            elif is_alchemy or isinstance(rows[0], (list, tuple, dict)):
                yield rows
            else:
                yield [tuple(row) for row in rows]

    def _from_arrow(
        self,
        *,
        batch_size: int | None,
        iter_batches: bool,
        schema_overrides: SchemaDict | None,
        infer_schema_length: int | None,
    ) -> DataFrame | Iterator[DataFrame] | None:
        """Return resultset data in Arrow format for frame init."""
        from polars import DataFrame

        try:
            # all ADBC drivers have the same method names
            driver = (
                "adbc" if self.driver_name.startswith("adbc_") else self.driver_name
            )
            driver_properties_list = ARROW_DRIVER_REGISTRY.get(driver, [])
            for i, driver_properties in enumerate(driver_properties_list, start=1):
                if ver := driver_properties["minimum_version"]:
                    # for ADBC drivers, the minimum version constraint is on the driver
                    # manager rather than the driver itself
                    driver_to_check = (
                        "adbc_driver_manager" if driver == "adbc" else self.driver_name
                    )
                    # if the minimum version constraint is not met, try additional
                    # driver properties with lower constraints
                    try:
                        self._check_module_version(driver_to_check, ver)
                    except ModuleUpgradeRequiredError:
                        if i < len(driver_properties_list):
                            continue
                        raise

                if iter_batches and (
                    driver_properties["exact_batch_size"] and not batch_size
                ):
                    msg = (
                        f"Cannot set `iter_batches` for {self.driver_name} "
                        "without also setting a non-zero `batch_size`"
                    )
                    raise ValueError(msg)  # noqa: TRY301

                frames = (
                    self._apply_overrides(batch, (schema_overrides or {}))
                    if isinstance(batch, DataFrame)
                    else from_arrow(batch, schema_overrides=schema_overrides)
                    for batch in self._fetch_arrow(
                        driver_properties,
                        iter_batches=iter_batches,
                        batch_size=batch_size,
                    )
                )
                return frames if iter_batches else next(frames)  # type: ignore[arg-type,return-value]
        except Exception as err:
            # eg: valid turbodbc/snowflake connection, but no arrow support
            # compiled in to the underlying driver (or on this connection)
            arrow_not_supported = (
                "does not support Apache Arrow",
                "Apache Arrow format is not supported",
            )
            if not any(e in str(err) for e in arrow_not_supported):
                raise

        return None

    def _from_rows(
        self,
        *,
        batch_size: int | None,
        iter_batches: bool,
        schema_overrides: SchemaDict | None,
        infer_schema_length: int | None,
    ) -> DataFrame | Iterator[DataFrame] | None:
        """Return resultset data row-wise for frame init."""
        from polars import DataFrame

        if iter_batches and not batch_size:
            msg = (
                "Cannot set `iter_batches` without also setting a non-zero `batch_size`"
            )
            raise ValueError(msg)

        if is_async := isinstance(original_result := self.result, Coroutine):
            self.result = _run_async(self.result)
        try:
            if hasattr(self.result, "fetchall"):
                if is_alchemy := (self.driver_name == "sqlalchemy"):
                    if hasattr(self.result, "cursor"):
                        cursor_desc = [
                            (d[0], d[1:]) for d in self.result.cursor.description
                        ]
                    elif hasattr(self.result, "_metadata"):
                        cursor_desc = [(k, None) for k in self.result._metadata.keys]
                    else:
                        msg = f"Unable to determine metadata from query result; {self.result!r}"
                        raise ValueError(msg)

                elif hasattr(self.result, "description"):
                    cursor_desc = [(d[0], d[1:]) for d in self.result.description]
                else:
                    cursor_desc = []

                schema_overrides = self._inject_type_overrides(
                    description=cursor_desc,
                    schema_overrides=(schema_overrides or {}),
                )
                result_columns = [nm for nm, _ in cursor_desc]
                frames = (
                    DataFrame(
                        data=rows,
                        schema=result_columns or None,
                        schema_overrides=schema_overrides,
                        infer_schema_length=infer_schema_length,
                        orient="row",
                    )
                    for rows in (
                        self._fetchmany_rows(
                            self.result,
                            batch_size=batch_size,
                            is_alchemy=is_alchemy,
                        )
                        if iter_batches
                        else [self._fetchall_rows(self.result, is_alchemy=is_alchemy)]  # type: ignore[list-item]
                    )
                )
                return frames if iter_batches else next(frames)  # type: ignore[arg-type]
            return None
        finally:
            if is_async:
                original_result.close()

    def _inject_type_overrides(
        self,
        description: list[tuple[str, Any]],
        schema_overrides: SchemaDict,
    ) -> SchemaDict:
        """
        Attempt basic dtype inference from a cursor description.

        Notes
        -----
        This is limited; the `type_code` description attr may contain almost anything,
        from strings or python types to driver-specific codes, classes, enums, etc.
        We currently only do the additional inference from string/python type values.
        (Further refinement will require per-driver module knowledge and lookups).
        """
        dupe_check = set()
        for nm, desc in description:
            if nm in dupe_check:
                msg = f"column {nm!r} appears more than once in the query/result cursor"
                raise DuplicateError(msg)
            elif desc is not None and nm not in schema_overrides:
                dtype = dtype_from_cursor_description(self.cursor, desc)
                if dtype is not None:
                    schema_overrides[nm] = dtype  # type: ignore[index]
            dupe_check.add(nm)

        return schema_overrides

    @staticmethod
    def _is_alchemy_async(conn: Any) -> bool:
        """Check if the given connection is SQLALchemy async."""
        try:
            from sqlalchemy.ext.asyncio import (
                AsyncConnection,
                AsyncSession,
                async_sessionmaker,
            )

            return isinstance(conn, (AsyncConnection, AsyncSession, async_sessionmaker))
        except ImportError:
            return False

    @staticmethod
    def _is_alchemy_engine(conn: Any) -> bool:
        """Check if the given connection is a SQLAlchemy Engine."""
        from sqlalchemy.engine import Engine

        if isinstance(conn, Engine):
            return True
        try:
            from sqlalchemy.ext.asyncio import AsyncEngine

            return isinstance(conn, AsyncEngine)
        except ImportError:
            return False

    @staticmethod
    def _is_alchemy_object(conn: Any) -> bool:
        """Check if the given connection is a SQLAlchemy object (of any kind)."""
        return type(conn).__module__.split(".", 1)[0] == "sqlalchemy"

    @staticmethod
    def _is_alchemy_session(conn: Any) -> bool:
        """Check if the given connection is a SQLAlchemy Session object."""
        from sqlalchemy.ext.asyncio import AsyncSession
        from sqlalchemy.orm import Session, sessionmaker

        if isinstance(conn, (AsyncSession, Session, sessionmaker)):
            return True

        try:
            from sqlalchemy.ext.asyncio import async_sessionmaker

            return isinstance(conn, async_sessionmaker)
        except ImportError:
            return False

    @staticmethod
    def _is_alchemy_result(result: Any) -> bool:
        """Check if the given result is a SQLAlchemy Result object."""
        try:
            from sqlalchemy.engine import CursorResult

            if isinstance(result, CursorResult):
                return True

            from sqlalchemy.ext.asyncio import AsyncResult

            return isinstance(result, AsyncResult)
        except ImportError:
            return False

    def _normalise_cursor(self, conn: Any) -> Cursor:
        """Normalise a connection object such that we have the query executor."""
        if self.driver_name == "sqlalchemy":
            if self._is_alchemy_session(conn):
                return conn
            else:
                # where possible, use the raw connection to access arrow integration
                if conn.engine.driver == "databricks-sql-python":
                    self.driver_name = "databricks"
                    return conn.engine.raw_connection().cursor()
                elif conn.engine.driver == "duckdb_engine":
                    self.driver_name = "duckdb"
                    return conn
                elif self._is_alchemy_engine(conn):
                    # note: if we create it, we can close it
                    self.can_close_cursor = True
                    return conn.connect()
                else:
                    return conn

        elif hasattr(conn, "cursor"):
            # connection has a dedicated cursor; prefer over direct execute
            cursor = cursor() if callable(cursor := conn.cursor) else cursor
            self.can_close_cursor = True
            return cursor

        elif hasattr(conn, "execute"):
            # can execute directly (given cursor, sqlalchemy connection, etc)
            return conn

        msg = (
            f"Unrecognised connection type {qualified_type_name(conn)!r}; no "
            "'execute' or 'cursor' method"
        )
        raise TypeError(msg)

    async def _sqlalchemy_async_execute(self, query: TextClause, **options: Any) -> Any:
        """Execute a query using an async SQLAlchemy connection."""
        is_session = self._is_alchemy_session(self.cursor)
        cursor = self.cursor.begin() if is_session else self.cursor  # type: ignore[attr-defined]
        async with cursor as conn:  # type: ignore[union-attr]
            if is_session and not hasattr(conn, "execute"):
                conn = conn.session
            result = await conn.execute(query, **options)
            return result

    def _sqlalchemy_setup(
        self, query: str | TextClause | Selectable, options: dict[str, Any]
    ) -> tuple[Any, dict[str, Any], str | TextClause | Selectable]:
        """Prepare a query for execution using a SQLAlchemy connection."""
        from sqlalchemy.orm import Session
        from sqlalchemy.sql import text
        from sqlalchemy.sql.elements import TextClause

        param_key = "parameters"
        cursor_execute = None
        if (
            isinstance(self.cursor, Session)
            and "parameters" in options
            and "params" not in options
        ):
            options = options.copy()
            options["params"] = options.pop("parameters")
            param_key = "params"

        params = options.get(param_key)
        is_async = self._is_alchemy_async(self.cursor)
        if (
            not is_async
            and isinstance(params, Sequence)
            and hasattr(self.cursor, "exec_driver_sql")
        ):
            cursor_execute = self.cursor.exec_driver_sql
            if isinstance(query, TextClause):
                query = str(query)
            if isinstance(params, list) and not all(
                isinstance(p, (dict, tuple)) for p in params
            ):
                options[param_key] = tuple(params)

        elif isinstance(query, str):
            query = text(query)

        if cursor_execute is None:
            cursor_execute = (
                self._sqlalchemy_async_execute if is_async else self.cursor.execute
            )
        return cursor_execute, options, query

    def execute(
        self,
        query: str | TextClause | Selectable,
        *,
        options: dict[str, Any] | None = None,
        select_queries_only: bool = True,
    ) -> Self:
        """Execute a query and reference the result set."""
        if select_queries_only and isinstance(query, str):
            q = re.search(r"\w{3,}", re.sub(r"/\*(.|[\r\n])*?\*/", "", query))
            if (query_type := "" if not q else q.group(0)) in _INVALID_QUERY_TYPES:
                msg = f"{query_type} statements are not valid 'read' queries"
                raise UnsuitableSQLError(msg)

        options = options or {}

        if self._is_alchemy_object(self.cursor):
            cursor_execute, options, query = self._sqlalchemy_setup(query, options)
        else:
            cursor_execute = self.cursor.execute

        # note: some cursor execute methods (eg: sqlite3) only take positional
        # params, hence the slightly convoluted resolution of the 'options' dict
        try:
            params = signature(cursor_execute).parameters
        except ValueError:
            params = {}  # type: ignore[assignment]

        if not options or any(
            p.kind in (Parameter.KEYWORD_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
            for p in params.values()
        ):
            result = cursor_execute(query, **options)
        else:
            positional_options = (
                options[o] for o in (params or options) if (not options or o in options)
            )
            result = cursor_execute(query, *positional_options)

        # note: some cursors execute in-place, some access results via a property
        result = self.cursor if (result is None or result is True) else result
        if self.driver_name == "duckdb" and self._is_alchemy_result(result):
            result = result.cursor

        self.result = result
        return self

    def to_polars(
        self,
        *,
        iter_batches: bool = False,
        batch_size: int | None = None,
        schema_overrides: SchemaDict | None = None,
        infer_schema_length: int | None = N_INFER_DEFAULT,
    ) -> DataFrame | Iterator[DataFrame]:
        """
        Convert the result set to a DataFrame.

        Wherever possible we try to return arrow-native data directly; only
        fall back to initialising with row-level data if no other option.
        """
        if self.result is None:
            msg = "cannot return a frame before executing a query"
            raise RuntimeError(msg)

        can_close = self.can_close_cursor

        if defer_cursor_close := (iter_batches and can_close):
            self.can_close_cursor = False

        for frame_init in (
            self._from_arrow,  # init from arrow-native data (where support exists)
            self._from_rows,  # row-wise fallback (sqlalchemy, dbapi2, pyodbc, etc)
        ):
            frame = frame_init(
                batch_size=batch_size,
                iter_batches=iter_batches,
                schema_overrides=schema_overrides,
                infer_schema_length=infer_schema_length,
            )
            if frame is not None:
                if defer_cursor_close:
                    frame = (
                        df
                        for df in CloseAfterFrameIter(
                            frame,
                            cursor=self.result,
                        )
                    )
                return frame

        msg = (
            f"Currently no support for {self.driver_name!r} connection {self.cursor!r}"
        )
        raise NotImplementedError(msg)
