from __future__ import annotations

import re
from typing import TYPE_CHECKING, Any, Literal, overload

from polars.datatypes import N_INFER_DEFAULT
from polars.dependencies import import_optional
from polars.io.database._cursor_proxies import ODBCCursorProxy
from polars.io.database._executor import ConnectionExecutor

if TYPE_CHECKING:
    from collections.abc import Iterator

    from sqlalchemy.sql.elements import TextClause
    from sqlalchemy.sql.expression import Selectable

    from polars import DataFrame
    from polars._typing import ConnectionOrCursor, DbReadEngine, SchemaDict


@overload
def read_database(
    query: str | TextClause | Selectable,
    connection: ConnectionOrCursor | str,
    *,
    iter_batches: Literal[False] = ...,
    batch_size: int | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    execute_options: dict[str, Any] | None = ...,
) -> DataFrame: ...


@overload
def read_database(
    query: str | TextClause | Selectable,
    connection: ConnectionOrCursor | str,
    *,
    iter_batches: Literal[True],
    batch_size: int | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    execute_options: dict[str, Any] | None = ...,
) -> Iterator[DataFrame]: ...


@overload
def read_database(
    query: str | TextClause | Selectable,
    connection: ConnectionOrCursor | str,
    *,
    iter_batches: bool,
    batch_size: int | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    execute_options: dict[str, Any] | None = ...,
) -> DataFrame | Iterator[DataFrame]: ...


def read_database(
    query: str | TextClause | Selectable,
    connection: ConnectionOrCursor | str,
    *,
    iter_batches: bool = False,
    batch_size: int | None = None,
    schema_overrides: SchemaDict | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    execute_options: dict[str, Any] | None = None,
) -> DataFrame | Iterator[DataFrame]:
    """
    Read the results of a SQL query into a DataFrame, given a connection object.

    Parameters
    ----------
    query
        SQL query to execute (if using a SQLAlchemy connection object this can
        be a suitable "Selectable", otherwise it is expected to be a string).
    connection
        An instantiated connection (or cursor/client object) that the query can be
        executed against. Can also pass a valid ODBC connection string, identified as
        such if it contains the string "Driver={...}", in which case the `arrow-odbc`
        package will be used to establish the connection and return Arrow-native data
        to Polars. Async driver connections are also supported, though this is currently
        considered unstable.

        .. warning::
            Use of asynchronous connections is currently considered **unstable**, and
            unexpected issues may arise; if this happens, please report them.
    iter_batches
        Return an iterator of DataFrames, where each DataFrame represents a batch of
        data returned by the query; this can be useful for processing large resultsets
        in a memory-efficient manner. If supported by the backend, this value is passed
        to the underlying query execution method (note that very low values will
        typically result in poor performance as it will cause many round-trips to the
        database as the data is returned). If the backend does not support changing
        the batch size then a single DataFrame is yielded from the iterator.
    batch_size
        Indicate the size of each batch when `iter_batches` is True (note that you can
        still set this when `iter_batches` is False, in which case the resulting
        DataFrame is constructed internally using batched return before being returned
        to you. Note that some backends (such as Snowflake) may support batch operation
        but not allow for an explicit size to be set; in this case you will still
        receive batches but their size is determined by the backend (in which case any
        value set here will be ignored).
    schema_overrides
        A dictionary mapping column names to dtypes, used to override the schema
        inferred from the query cursor or given by the incoming Arrow data (depending
        on driver/backend). This can be useful if the given types can be more precisely
        defined (for example, if you know that a given column can be declared as `u32`
        instead of `i64`).
    infer_schema_length
        The maximum number of rows to scan for schema inference. If set to `None`, the
        full data may be scanned *(this can be slow)*. This parameter only applies if
        the data is read as a sequence of rows and the `schema_overrides` parameter
        is not set for the given column; Arrow-aware drivers also ignore this value.
    execute_options
        These options will be passed through into the underlying query execution method
        as kwargs. In the case of connections made using an ODBC string (which use
        `arrow-odbc`) these options are passed to the `read_arrow_batches_from_odbc`
        method.

    Notes
    -----
    * This function supports a wide range of native database drivers (ranging from local
      databases such as SQLite to large cloud databases such as Snowflake), as well as
      generic libraries such as ADBC, SQLAlchemy and various flavours of ODBC. If the
      backend supports returning Arrow data directly then this facility will be used to
      efficiently instantiate the DataFrame; otherwise, the DataFrame is initialised
      from row-wise data.

    * Support for Arrow Flight SQL data is available via the `adbc-driver-flightsql`
      package; see https://arrow.apache.org/adbc/current/driver/flight_sql.html for
      more details about using this driver (notable databases implementing Flight SQL
      include Dremio and InfluxDB).

    * The `read_database_uri` function can be noticeably faster than `read_database`
      if you are using a SQLAlchemy or DBAPI2 connection, as `connectorx` and `adbc`
      optimise translation of the result set into Arrow format. Note that you can
      determine a connection's URI from a SQLAlchemy engine object by calling
      `conn.engine.url.render_as_string(hide_password=False)`.

    * If Polars has to create a cursor from your connection in order to execute the
      query then that cursor will be automatically closed when the query completes;
      however, Polars will *never* close any other open connection or cursor.

    * Polars is able to support more than just relational databases and SQL queries
      through this function. For example, you can load local graph database results
      from a `KùzuDB` connection in conjunction with a Cypher query, or use SurrealQL
      with SurrealDB.

    See Also
    --------
    read_database_uri : Create a DataFrame from a SQL query using a URI string.

    Examples
    --------
    Instantiate a DataFrame from a SQL query against a user-supplied connection:

    >>> df = pl.read_database(
    ...     query="SELECT * FROM test_data",
    ...     connection=user_conn,
    ...     schema_overrides={"normalised_score": pl.UInt8},
    ... )  # doctest: +SKIP

    Use a parameterised SQLAlchemy query, passing named values via `execute_options`:

    >>> df = pl.read_database(
    ...     query="SELECT * FROM test_data WHERE metric > :value",
    ...     connection=alchemy_conn,
    ...     execute_options={"parameters": {"value": 0}},
    ... )  # doctest: +SKIP

    Use 'qmark' style parameterisation; values are still passed via `execute_options`,
    but in this case the "parameters" value is a sequence of literals, not a dict:

    >>> df = pl.read_database(
    ...     query="SELECT * FROM test_data WHERE metric > ?",
    ...     connection=alchemy_conn,
    ...     execute_options={"parameters": [0]},
    ... )  # doctest: +SKIP

    Instantiate a DataFrame using an ODBC connection string (requires the `arrow-odbc`
    package) setting upper limits on the buffer size of variadic text/binary columns,
    returning the result as an iterator over DataFrames that each contain 1000 rows:

    >>> for df in pl.read_database(
    ...     query="SELECT * FROM test_data",
    ...     connection="Driver={PostgreSQL};Server=localhost;Port=5432;Database=test;Uid=usr;Pwd=",
    ...     execute_options={"max_text_size": 512, "max_binary_size": 1024},
    ...     iter_batches=True,
    ...     batch_size=1000,
    ... ):
    ...     do_something(df)  # doctest: +SKIP

    Load graph database results from a `KùzuDB` connection and a Cypher query:

    >>> df = pl.read_database(
    ...     query="MATCH (a:User)-[f:Follows]->(b:User) RETURN a.name, f.since, b.name",
    ...     connection=kuzu_db_conn,
    ... )  # doctest: +SKIP

    Load data from an asynchronous SQLAlchemy driver/engine; note that asynchronous
    connections and sessions are also supported here:

    >>> from sqlalchemy.ext.asyncio import create_async_engine
    >>> async_engine = create_async_engine("sqlite+aiosqlite:///test.db")
    >>> df = pl.read_database(
    ...     query="SELECT * FROM test_data",
    ...     connection=async_engine,
    ... )  # doctest: +SKIP

    Load data from an asynchronous SurrealDB client connection object; note that
    both the WS (`Surreal`) and HTTP (`SurrealHTTP`) clients are supported:

    >>> import asyncio
    >>> async def surreal_query_to_frame(query: str, url: str):
    ...     async with Surreal(url) as client:
    ...         await client.use(namespace="test", database="test")
    ...         return pl.read_database(query=query, connection=client)
    >>> df = asyncio.run(
    ...     surreal_query_to_frame(
    ...         query="SELECT * FROM test_data",
    ...         url="ws://localhost:8000/rpc",
    ...     )
    ... )  # doctest: +SKIP

    """  # noqa: W505
    if isinstance(connection, str):
        # check for odbc connection string
        if re.search(r"\bdriver\s*=\s*{[^}]+?}", connection, re.IGNORECASE):
            _ = import_optional(
                module_name="arrow_odbc",
                err_prefix="use of ODBC connection string requires the",
                err_suffix="package",
            )
            connection = ODBCCursorProxy(connection)
        elif "://" in connection:
            # otherwise looks like a mistaken call to read_database_uri
            msg = "use of string URI is invalid here; call `read_database_uri` instead"
            raise ValueError(msg)
        else:
            msg = "unable to identify string connection as valid ODBC (no driver)"
            raise ValueError(msg)

    # return frame from arbitrary connections using the executor abstraction
    with ConnectionExecutor(connection) as cx:
        return cx.execute(
            query=query,
            options=execute_options,
        ).to_polars(
            batch_size=batch_size,
            iter_batches=iter_batches,
            schema_overrides=schema_overrides,
            infer_schema_length=infer_schema_length,
        )


@overload
def read_database_uri(
    query: str,
    uri: str,
    *,
    partition_on: str | None = None,
    partition_range: tuple[int, int] | None = None,
    partition_num: int | None = None,
    protocol: str | None = None,
    engine: Literal["adbc"],
    schema_overrides: SchemaDict | None = None,
    execute_options: dict[str, Any] | None = None,
) -> DataFrame: ...


@overload
def read_database_uri(
    query: list[str] | str,
    uri: str,
    *,
    partition_on: str | None = None,
    partition_range: tuple[int, int] | None = None,
    partition_num: int | None = None,
    protocol: str | None = None,
    engine: Literal["connectorx"] | None = None,
    schema_overrides: SchemaDict | None = None,
    execute_options: None = None,
) -> DataFrame: ...


@overload
def read_database_uri(
    query: str,
    uri: str,
    *,
    partition_on: str | None = None,
    partition_range: tuple[int, int] | None = None,
    partition_num: int | None = None,
    protocol: str | None = None,
    engine: DbReadEngine | None = None,
    schema_overrides: None = None,
    execute_options: dict[str, Any] | None = None,
) -> DataFrame: ...


def read_database_uri(
    query: list[str] | str,
    uri: str,
    *,
    partition_on: str | None = None,
    partition_range: tuple[int, int] | None = None,
    partition_num: int | None = None,
    protocol: str | None = None,
    engine: DbReadEngine | None = None,
    schema_overrides: SchemaDict | None = None,
    execute_options: dict[str, Any] | None = None,
) -> DataFrame:
    """
    Read the results of a SQL query into a DataFrame, given a URI.

    Parameters
    ----------
    query
        Raw SQL query (or queries).
    uri
        A connectorx or ADBC connection URI string that starts with the backend's
        driver name, for example:

        * "postgresql://user:pass@server:port/database"
        * "snowflake://user:pass@account/database/schema?warehouse=warehouse&role=role"

        The caller is responsible for escaping any special characters in the string,
        which will be passed "as-is" to the underlying engine (this is most often
        required when coming across special characters in the password).
    partition_on
        The column on which to partition the result (connectorx).
    partition_range
        The value range of the partition column (connectorx).
    partition_num
        How many partitions to generate (connectorx).
    protocol
        Backend-specific transfer protocol directive (connectorx); see connectorx
        documentation for more details.
    engine : {'connectorx', 'adbc'}
        Selects the engine used for reading the database (defaulting to connectorx):

        * `'connectorx'`
          Supports a range of databases, such as PostgreSQL, Redshift, MySQL, MariaDB,
          Clickhouse, Oracle, BigQuery, SQL Server, and so on. For an up-to-date list
          please see the connectorx docs:
          https://github.com/sfu-db/connector-x#supported-sources--destinations
        * `'adbc'`
          Currently there is limited support for this engine, with a relatively small
          number of drivers available, most of which are still in development. For
          an up-to-date list of drivers please see the ADBC docs:
          https://arrow.apache.org/adbc/
    schema_overrides
        A dictionary mapping column names to dtypes, used to override the schema
        given in the data returned by the query.
    execute_options
        These options will be passed to the underlying query execution method as
        kwargs. Note that connectorx does not support this parameter.

    Notes
    -----
    For `connectorx`, ensure that you have `connectorx>=0.3.2`. The documentation
    is available `here <https://sfu-db.github.io/connector-x/intro.html>`_.

    For `adbc` you will need to have installed `pyarrow` and the ADBC driver associated
    with the backend you are connecting to, eg: `adbc-driver-postgresql`.

    If your password contains special characters, you will need to escape them.
    This will usually require the use of a URL-escaping function, for example:

    >>> from urllib.parse import quote, quote_plus
    >>> quote_plus("pass word?")
    'pass+word%3F'
    >>> quote("pass word?")
    'pass%20word%3F'

    See Also
    --------
    read_database : Create a DataFrame from a SQL query using a connection object.

    Examples
    --------
    Create a DataFrame from a SQL query using a single thread:

    >>> uri = "postgresql://username:password@server:port/database"
    >>> query = "SELECT * FROM lineitem"
    >>> pl.read_database_uri(query, uri)  # doctest: +SKIP

    Create a DataFrame in parallel using 10 threads by automatically partitioning
    the provided SQL on the partition column:

    >>> uri = "postgresql://username:password@server:port/database"
    >>> query = "SELECT * FROM lineitem"
    >>> pl.read_database_uri(
    ...     query,
    ...     uri,
    ...     partition_on="partition_col",
    ...     partition_num=10,
    ...     engine="connectorx",
    ... )  # doctest: +SKIP

    Create a DataFrame in parallel using 2 threads by explicitly providing two
    SQL queries:

    >>> uri = "postgresql://username:password@server:port/database"
    >>> queries = [
    ...     "SELECT * FROM lineitem WHERE partition_col <= 10",
    ...     "SELECT * FROM lineitem WHERE partition_col > 10",
    ... ]
    >>> pl.read_database_uri(queries, uri, engine="connectorx")  # doctest: +SKIP

    Read data from Snowflake using the ADBC driver:

    >>> df = pl.read_database_uri(
    ...     "SELECT * FROM test_table",
    ...     "snowflake://user:pass@company-org/testdb/public?warehouse=test&role=myrole",
    ...     engine="adbc",
    ... )  # doctest: +SKIP
    """
    from polars.io.database._utils import _read_sql_adbc, _read_sql_connectorx

    if not isinstance(uri, str):
        msg = f"expected connection to be a URI string; found {type(uri).__name__!r}"
        raise TypeError(msg)
    elif engine is None:
        engine = "connectorx"

    if engine == "connectorx":
        if execute_options:
            msg = "the 'connectorx' engine does not support use of `execute_options`"
            raise ValueError(msg)
        return _read_sql_connectorx(
            query,
            connection_uri=uri,
            partition_on=partition_on,
            partition_range=partition_range,
            partition_num=partition_num,
            protocol=protocol,
            schema_overrides=schema_overrides,
        )
    elif engine == "adbc":
        if not isinstance(query, str):
            msg = "only a single SQL query string is accepted for adbc"
            raise ValueError(msg)
        return _read_sql_adbc(
            query,
            connection_uri=uri,
            schema_overrides=schema_overrides,
            execute_options=execute_options,
        )
    else:
        msg = f"engine must be one of {{'connectorx', 'adbc'}}, got {engine!r}"
        raise ValueError(msg)
