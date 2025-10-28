"""
Polars: Blazingly fast DataFrames
=================================

Polars is a fast, open-source library for data manipulation with an expressive, typed API.

Basic usage:

   >>> import polars as pl
   >>> df = pl.DataFrame(
   ...     {
   ...         "name": ["Alice", "Bob", "Charlie"],
   ...         "age": [25, 30, 35],
   ...         "city": ["New York", "London", "Tokyo"],
   ...     }
   ... )
   >>> df.filter(pl.col("age") > 28)
   shape: (2, 3)
   ┌─────────┬─────┬────────┐
   │ name    ┆ age ┆ city   │
   │ ---     ┆ --- ┆ ---    │
   │ str     ┆ i64 ┆ str    │
   ╞═════════╪═════╪════════╡
   │ Bob     ┆ 30  ┆ London │
   │ Charlie ┆ 35  ┆ Tokyo  │
   └─────────┴─────┴────────┘

User Guide: https://docs.pola.rs/
Python API Documentation: https://docs.pola.rs/api/python/stable/
Source Code: https://github.com/pola-rs/polars
"""  # noqa: D400, W505, D205

import contextlib

with contextlib.suppress(ImportError):  # Module not available when building docs
    # We also configure the allocator before importing the Polars Rust bindings.
    # See https://github.com/pola-rs/polars/issues/18088,
    # https://github.com/pola-rs/polars/pull/21829.
    import os

    jemalloc_conf = "dirty_decay_ms:500,muzzy_decay_ms:-1"
    if os.environ.get("POLARS_THP") == "1":
        jemalloc_conf += ",thp:always,metadata_thp:always"
    if override := os.environ.get("_RJEM_MALLOC_CONF"):
        jemalloc_conf += "," + override
    os.environ["_RJEM_MALLOC_CONF"] = jemalloc_conf

    # Initialize polars on the rust side. This function is highly
    # unsafe and should only be called once.
    from polars._plr import __register_startup_deps

    __register_startup_deps()

from typing import TYPE_CHECKING, Any

from polars import api, exceptions, plugins, selectors
from polars._utils.polars_version import get_polars_version as _get_polars_version

# TODO: remove need for importing wrap utils at top level
from polars._utils.wrap import wrap_df, wrap_s  # noqa: F401
from polars.catalog.unity import Catalog
from polars.config import Config
from polars.convert import (
    from_arrow,
    from_dataframe,
    from_dict,
    from_dicts,
    from_numpy,
    from_pandas,
    from_records,
    from_repr,
    from_torch,
    json_normalize,
)
from polars.dataframe import DataFrame
from polars.datatype_expr import DataTypeExpr
from polars.datatypes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    Categories,
    DataType,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Field,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Null,
    Object,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Unknown,
    Utf8,
)
from polars.expr import Expr
from polars.functions import (
    align_frames,
    all,
    all_horizontal,
    any,
    any_horizontal,
    approx_n_unique,
    arange,
    arctan2,
    arctan2d,
    arg_sort_by,
    arg_where,
    business_day_count,
    coalesce,
    col,
    collect_all,
    collect_all_async,
    concat,
    concat_arr,
    concat_list,
    concat_str,
    corr,
    count,
    cov,
    cum_count,
    cum_fold,
    cum_reduce,
    cum_sum,
    cum_sum_horizontal,
    date,
    date_range,
    date_ranges,
    datetime,
    datetime_range,
    datetime_ranges,
    dtype_of,
    duration,
    element,
    escape_regex,
    exclude,
    explain_all,
    field,
    first,
    fold,
    format,
    from_epoch,
    groups,
    head,
    implode,
    int_range,
    int_ranges,
    last,
    len,
    linear_space,
    linear_spaces,
    lit,
    map_batches,
    map_groups,
    max,
    max_horizontal,
    mean,
    mean_horizontal,
    median,
    min,
    min_horizontal,
    n_unique,
    nth,
    ones,
    quantile,
    reduce,
    repeat,
    rolling_corr,
    rolling_cov,
    row_index,
    select,
    self_dtype,
    set_random_seed,
    sql_expr,
    std,
    struct,
    struct_with_fields,
    sum,
    sum_horizontal,
    tail,
    time,
    time_range,
    time_ranges,
    union,
    var,
    when,
    zeros,
)
from polars.interchange import CompatLevel
from polars.io import (
    BasePartitionContext,
    KeyedPartition,
    KeyedPartitionContext,
    PartitionByKey,
    PartitionMaxSize,
    PartitionParted,
    ScanCastOptions,
    defer,
    read_avro,
    read_clipboard,
    read_csv,
    read_csv_batched,
    read_database,
    read_database_uri,
    read_delta,
    read_excel,
    read_ipc,
    read_ipc_schema,
    read_ipc_stream,
    read_json,
    read_ndjson,
    read_ods,
    read_parquet,
    read_parquet_metadata,
    read_parquet_schema,
    scan_csv,
    scan_delta,
    scan_iceberg,
    scan_ipc,
    scan_ndjson,
    scan_parquet,
    scan_pyarrow_dataset,
)
from polars.io.cloud import (
    CredentialProvider,
    CredentialProviderAWS,
    CredentialProviderAzure,
    CredentialProviderFunction,
    CredentialProviderFunctionReturn,
    CredentialProviderGCP,
)
from polars.lazyframe import GPUEngine, LazyFrame, QueryOptFlags
from polars.meta import (
    build_info,
    get_index_type,
    show_versions,
    thread_pool_size,
    threadpool_size,
)
from polars.schema import Schema
from polars.series import Series
from polars.sql import SQLContext, sql
from polars.string_cache import (
    StringCache,
    disable_string_cache,
    enable_string_cache,
    using_string_cache,
)

__version__: str = _get_polars_version()
del _get_polars_version

__all__ = [
    # modules
    "api",
    "exceptions",
    "plugins",
    "selectors",
    # core classes
    "DataFrame",
    "Expr",
    "LazyFrame",
    "Series",
    # Engine configuration
    "GPUEngine",
    # schema
    "Schema",
    # datatype_expr
    "DataTypeExpr",
    # datatypes
    "Array",
    "Binary",
    "Boolean",
    "Categorical",
    "Categories",
    "DataType",
    "Date",
    "Datetime",
    "Decimal",
    "Duration",
    "Enum",
    "Field",
    "Float32",
    "Float64",
    "Int8",
    "Int16",
    "Int32",
    "Int64",
    "Int128",
    "List",
    "Null",
    "Object",
    "String",
    "Struct",
    "Time",
    "UInt8",
    "UInt16",
    "UInt32",
    "UInt64",
    "UInt128",
    "Unknown",
    "Utf8",
    # polars.io
    "defer",
    "KeyedPartition",
    "BasePartitionContext",
    "KeyedPartitionContext",
    "PartitionByKey",
    "PartitionMaxSize",
    "PartitionParted",
    "ScanCastOptions",
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
    "Catalog",
    # polars.io.cloud
    "CredentialProvider",
    "CredentialProviderAWS",
    "CredentialProviderAzure",
    "CredentialProviderFunction",
    "CredentialProviderFunctionReturn",
    "CredentialProviderGCP",
    # polars.stringcache
    "StringCache",
    "disable_string_cache",
    "enable_string_cache",
    "using_string_cache",
    # polars.config
    "Config",
    # polars.functions.whenthen
    "when",
    # polars.functions
    "align_frames",
    "arg_where",
    "business_day_count",
    "concat",
    "union",
    "dtype_of",
    "struct_with_fields",
    "date_range",
    "date_ranges",
    "datetime_range",
    "datetime_ranges",
    "element",
    "ones",
    "repeat",
    "self_dtype",
    "time_range",
    "time_ranges",
    "zeros",
    "escape_regex",
    # polars.functions.aggregation
    "all",
    "all_horizontal",
    "any",
    "any_horizontal",
    "cum_sum",
    "cum_sum_horizontal",
    "max",
    "max_horizontal",
    "mean_horizontal",
    "min",
    "min_horizontal",
    "sum",
    "sum_horizontal",
    # polars.functions.lazy
    "approx_n_unique",
    "arange",
    "arctan2",
    "arctan2d",
    "arg_sort_by",
    "coalesce",
    "col",
    "collect_all",
    "collect_all_async",
    "concat_arr",
    "concat_list",
    "concat_str",
    "corr",
    "count",
    "cov",
    "cum_count",
    "cum_fold",
    "cum_reduce",
    "date",
    "datetime",
    "duration",
    "exclude",
    "explain_all",
    "field",
    "first",
    "fold",
    "format",
    "from_epoch",
    "groups",
    "head",
    "implode",
    "int_range",
    "int_ranges",
    "last",
    "linear_space",
    "linear_spaces",
    "lit",
    "map_batches",
    "map_groups",
    "mean",
    "median",
    "n_unique",
    "nth",
    "quantile",
    "reduce",
    "rolling_corr",
    "rolling_cov",
    "row_index",
    "select",
    "std",
    "struct",
    "tail",
    "time",
    "var",
    # polars.functions.len
    "len",
    # polars.functions.random
    "set_random_seed",
    # polars.convert
    "from_arrow",
    "from_dataframe",
    "from_dict",
    "from_dicts",
    "from_numpy",
    "from_pandas",
    "from_records",
    "from_repr",
    "from_torch",
    "json_normalize",
    # polars.meta
    "build_info",
    "get_index_type",
    "show_versions",
    "thread_pool_size",
    "threadpool_size",
    # polars.sql
    "SQLContext",
    "sql",
    "sql_expr",
    "CompatLevel",
    # optimization
    "QueryOptFlags",
]


if not TYPE_CHECKING:
    with contextlib.suppress(ImportError):  # Module not available when building docs
        import polars._plr as plr

    # This causes typechecking to resolve any Polars module attribute
    # as Any regardless of existence so we check for TYPE_CHECKING, see #24334.
    def __getattr__(name: str) -> Any:
        # Backwards compatibility for plugins. This used to be called `polars.polars`,
        # but is now `polars._plr`.
        if name == "polars":
            return plr
        elif name == "_allocator":
            return plr._allocator

        # Deprecate re-export of exceptions at top-level
        if name in dir(exceptions):
            from polars._utils.deprecation import issue_deprecation_warning

            issue_deprecation_warning(
                message=(
                    f"accessing `{name}` from the top-level `polars` module was deprecated "
                    "in version 1.0.0. Import it directly from the `polars.exceptions` module "
                    f"instead, e.g.: `from polars.exceptions import {name}`"
                ),
            )
            return getattr(exceptions, name)

        # Deprecate data type groups at top-level
        import polars.datatypes.group as dtgroup

        if name in dir(dtgroup):
            from polars._utils.deprecation import issue_deprecation_warning

            issue_deprecation_warning(
                message=(
                    f"`{name}` was deprecated in version 1.0.0. Define your own data type groups or "
                    "use the `polars.selectors` module for selecting columns of a certain data type."
                ),
            )
            return getattr(dtgroup, name)

        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)
