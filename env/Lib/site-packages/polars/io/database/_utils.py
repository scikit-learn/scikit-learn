from __future__ import annotations

import re
from typing import TYPE_CHECKING, Any

from polars.convert import from_arrow
from polars.dependencies import import_optional

if TYPE_CHECKING:
    from collections.abc import Coroutine

    from polars import DataFrame
    from polars._typing import SchemaDict


def _run_async(co: Coroutine[Any, Any, Any]) -> Any:
    """Run asynchronous code as if it was synchronous."""
    import asyncio

    from polars._utils.unstable import issue_unstable_warning
    from polars.dependencies import import_optional

    issue_unstable_warning(
        "Use of asynchronous connections is currently considered unstable "
        "and unexpected issues may arise; if this happens, please report them."
    )
    nest_asyncio = import_optional("nest_asyncio")
    nest_asyncio.apply()
    return asyncio.run(co)


def _read_sql_connectorx(
    query: str | list[str],
    connection_uri: str,
    partition_on: str | None = None,
    partition_range: tuple[int, int] | None = None,
    partition_num: int | None = None,
    protocol: str | None = None,
    schema_overrides: SchemaDict | None = None,
) -> DataFrame:
    cx = import_optional("connectorx")
    try:
        tbl = cx.read_sql(
            conn=connection_uri,
            query=query,
            return_type="arrow2",
            partition_on=partition_on,
            partition_range=partition_range,
            partition_num=partition_num,
            protocol=protocol,
        )
    except BaseException as err:
        # basic sanitisation of /user:pass/ credentials exposed in connectorx errs
        errmsg = re.sub("://[^:]+:[^:]+@", "://***:***@", str(err))
        raise type(err)(errmsg) from err

    return from_arrow(tbl, schema_overrides=schema_overrides)  # type: ignore[return-value]


def _read_sql_adbc(
    query: str,
    connection_uri: str,
    schema_overrides: SchemaDict | None,
    execute_options: dict[str, Any] | None = None,
) -> DataFrame:
    with _open_adbc_connection(connection_uri) as conn, conn.cursor() as cursor:
        cursor.execute(query, **(execute_options or {}))
        tbl = cursor.fetch_arrow_table()
    return from_arrow(tbl, schema_overrides=schema_overrides)  # type: ignore[return-value]


def _open_adbc_connection(connection_uri: str) -> Any:
    driver_name = connection_uri.split(":", 1)[0].lower()

    # map uri prefix to module when not 1:1
    module_suffix_map: dict[str, str] = {
        "postgres": "postgresql",
    }
    module_suffix = module_suffix_map.get(driver_name, driver_name)
    module_name = f"adbc_driver_{module_suffix}.dbapi"

    adbc_driver = import_optional(
        module_name,
        err_prefix="ADBC",
        err_suffix="driver not detected",
        install_message=f"If ADBC supports this database, please run: pip install adbc-driver-{driver_name} pyarrow",
    )

    # some backends require the driver name to be stripped from the URI
    if driver_name in ("sqlite", "snowflake"):
        connection_uri = re.sub(f"^{driver_name}:/{{,3}}", "", connection_uri)

    return adbc_driver.connect(connection_uri)
