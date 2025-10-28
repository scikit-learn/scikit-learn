"""Functions for reading data."""

from polars.io.avro import read_avro
from polars.io.clipboard import read_clipboard
from polars.io.csv import read_csv, read_csv_batched, scan_csv
from polars.io.database import read_database, read_database_uri
from polars.io.delta import read_delta, scan_delta
from polars.io.iceberg import scan_iceberg
from polars.io.ipc import read_ipc, read_ipc_schema, read_ipc_stream, scan_ipc
from polars.io.json import read_json
from polars.io.ndjson import read_ndjson, scan_ndjson
from polars.io.parquet import (
    read_parquet,
    read_parquet_metadata,
    read_parquet_schema,
    scan_parquet,
)
from polars.io.partition import (
    BasePartitionContext,
    KeyedPartition,
    KeyedPartitionContext,
    PartitionByKey,
    PartitionMaxSize,
    PartitionParted,
)
from polars.io.plugins import _defer as defer
from polars.io.pyarrow_dataset import scan_pyarrow_dataset
from polars.io.scan_options import ScanCastOptions
from polars.io.spreadsheet import read_excel, read_ods

__all__ = [
    "defer",
    "PartitionByKey",
    "PartitionMaxSize",
    "PartitionParted",
    "KeyedPartition",
    "BasePartitionContext",
    "KeyedPartitionContext",
    "read_avro",
    "read_clipboard",
    "read_csv",
    "read_csv_batched",
    "read_database",
    "read_database_uri",
    "read_delta",
    "read_excel",
    "read_ipc",
    "read_ipc_schema",
    "read_ipc_stream",
    "read_json",
    "read_ndjson",
    "read_ods",
    "read_parquet",
    "read_parquet_metadata",
    "read_parquet_schema",
    "scan_csv",
    "scan_delta",
    "scan_iceberg",
    "scan_ipc",
    "scan_ndjson",
    "scan_parquet",
    "scan_pyarrow_dataset",
    "ScanCastOptions",
]
