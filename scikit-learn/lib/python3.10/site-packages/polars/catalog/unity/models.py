from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from polars._utils.unstable import issue_unstable_warning
from polars.exceptions import DuplicateError
from polars.schema import Schema

if TYPE_CHECKING:
    from datetime import datetime

    from polars.datatypes.classes import DataType


@dataclass
class CatalogInfo:
    """Information for a catalog within a metastore."""

    name: str
    comment: str | None
    properties: dict[str, str]
    options: dict[str, str]
    storage_location: str | None
    created_at: datetime | None
    created_by: str | None
    updated_at: datetime | None
    updated_by: str | None


@dataclass
class NamespaceInfo:
    """
    Information for a namespace within a catalog.

    This is also known by the name "schema" in unity catalog terminology.
    """

    name: str
    comment: str | None
    properties: dict[str, str]
    storage_location: str | None
    created_at: datetime | None
    created_by: str | None
    updated_at: datetime | None
    updated_by: str | None


@dataclass
class TableInfo:
    """Information for a catalog table."""

    name: str
    comment: str | None
    table_id: str
    table_type: TableType
    storage_location: str | None
    data_source_format: DataSourceFormat | None
    columns: list[ColumnInfo] | None
    properties: dict[str, str]
    created_at: datetime | None
    created_by: str | None
    updated_at: datetime | None
    updated_by: str | None

    def get_polars_schema(self) -> Schema | None:
        """
        Get the native polars schema of this table.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
        """
        issue_unstable_warning(
            "`get_polars_schema` functionality is considered unstable."
        )
        if self.columns is None:
            return None

        schema = Schema()

        for column_info in self.columns:
            if column_info.name in schema:
                msg = f"duplicate column name: {column_info.name}"
                raise DuplicateError(msg)
            schema[column_info.name] = column_info.get_polars_dtype()

        return schema


@dataclass
class ColumnInfo:
    """Information for a column within a catalog table."""

    name: str
    type_name: str
    type_text: str
    type_json: str
    position: int | None
    comment: str | None
    partition_index: int | None

    def get_polars_dtype(self) -> DataType:
        """
        Get the native polars datatype of this column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
        """
        issue_unstable_warning(
            "`get_polars_dtype` functionality is considered unstable."
        )

        from polars._plr import PyCatalogClient

        return PyCatalogClient.type_json_to_polars_type(self.type_json)


TableType = Literal[
    "MANAGED",
    "EXTERNAL",
    "VIEW",
    "MATERIALIZED_VIEW",
    "STREAMING_TABLE",
    "MANAGED_SHALLOW_CLONE",
    "FOREIGN",
    "EXTERNAL_SHALLOW_CLONE",
]

DataSourceFormat = Literal[
    "DELTA",
    "CSV",
    "JSON",
    "AVRO",
    "PARQUET",
    "ORC",
    "TEXT",
    "UNITY_CATALOG",
    "DELTASHARING",
    "DATABRICKS_FORMAT",
    "REDSHIFT_FORMAT",
    "SNOWFLAKE_FORMAT",
    "SQLDW_FORMAT",
    "SALESFORCE_FORMAT",
    "BIGQUERY_FORMAT",
    "NETSUITE_FORMAT",
    "WORKDAY_RAAS_FORMAT",
    "HIVE_SERDE",
    "HIVE_CUSTOM",
    "VECTOR_INDEX_FORMAT",
]
