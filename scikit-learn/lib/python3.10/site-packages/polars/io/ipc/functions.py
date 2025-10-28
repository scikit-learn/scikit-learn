from __future__ import annotations

import contextlib
import os
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Literal

import polars._reexport as pl
import polars.functions as F
from polars._dependencies import import_optional
from polars._utils.deprecation import deprecate_renamed_parameter
from polars._utils.various import (
    is_str_sequence,
    normalize_filepath,
)
from polars._utils.wrap import wrap_df, wrap_ldf
from polars.io._utils import (
    get_sources,
    is_glob_pattern,
    is_local_file,
    parse_columns_arg,
    parse_row_index_args,
    prepare_file_arg,
)
from polars.io.cloud.credential_provider._builder import (
    _init_credential_provider_builder,
)
from polars.io.scan_options._options import ScanOptions

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyDataFrame, PyLazyFrame
    from polars._plr import read_ipc_schema as _read_ipc_schema

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars import DataFrame, DataType, LazyFrame
    from polars._typing import SchemaDict
    from polars.io.cloud import CredentialProviderFunction


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def read_ipc(
    source: str | Path | IO[bytes] | bytes,
    *,
    columns: list[int] | list[str] | None = None,
    n_rows: int | None = None,
    use_pyarrow: bool = False,
    memory_map: bool = True,
    storage_options: dict[str, Any] | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    rechunk: bool = True,
) -> DataFrame:
    """
    Read into a DataFrame from Arrow IPC (Feather v2) file.

    See "File or Random Access format" on https://arrow.apache.org/docs/python/ipc.html.
    Arrow IPC files are also known as Feather (v2) files.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). If `fsspec` is installed, it might be used
        to open remote files. For file-like objects, the stream position may not be
        updated accordingly after reading.
    columns
        Columns to select. Accepts a list of column indices (starting at zero) or a list
        of column names.
    n_rows
        Stop reading from IPC file after reading `n_rows`.
        Only valid when `use_pyarrow=False`.
    use_pyarrow
        Use pyarrow or the native Rust reader.
    memory_map
        Try to memory map the file. This can greatly improve performance on repeated
        queries as the OS may cache pages.
        Only uncompressed IPC files can be memory mapped.
    storage_options
        Extra options that make sense for `fsspec.open()` or a particular storage
        connection, e.g. host, port, username, password, etc.
    row_index_name
        Insert a row index column with the given name into the DataFrame as the first
        column. If set to `None` (default), no row index column is created.
    row_index_offset
        Start the row index at this offset. Cannot be negative.
        Only used if `row_index_name` is set.
    rechunk
        Make sure that all data is contiguous.

    Returns
    -------
    DataFrame

    See Also
    --------
    scan_ipc : Lazily read from an IPC file or multiple files via glob patterns.

    Warnings
    --------
    Calling `read_ipc().lazy()` is an antipattern as this forces Polars to materialize
    a full csv file and therefore cannot push any optimizations into the reader.
    Therefore always prefer `scan_ipc` if you want to work with `LazyFrame` s.

    If `memory_map` is set, the bytes on disk are mapped 1:1 to memory.
    That means that you cannot write to the same filename.
    E.g. `pl.read_ipc("my_file.arrow").write_ipc("my_file.arrow")` will fail.
    """
    if (
        # Check that it is not a BytesIO object
        isinstance(v := source, (str, Path))
    ) and (
        # HuggingFace only for now ⊂( ◜◒◝ )⊃
        (is_hf := str(v).startswith("hf://"))
        # Also dispatch on FORCE_ASYNC, so that this codepath gets run
        # through by our test suite during CI.
        or os.getenv("POLARS_FORCE_ASYNC") == "1"
        # TODO: Dispatch all paths to `scan_ipc` - this will need a breaking
        # change to the `storage_options` parameter.
    ):
        if is_hf and use_pyarrow:
            msg = "`use_pyarrow=True` is not supported for Hugging Face"
            raise ValueError(msg)

        lf = scan_ipc(
            source,
            n_rows=n_rows,
            storage_options=storage_options,
            row_index_name=row_index_name,
            row_index_offset=row_index_offset,
            rechunk=rechunk,
        )

        if columns:
            if isinstance(columns[0], int):
                lf = lf.select(F.nth(columns))  # type: ignore[arg-type]
            else:
                lf = lf.select(columns)

        df = lf.collect()

        return df

    if use_pyarrow and n_rows and not memory_map:
        msg = "`n_rows` cannot be used with `use_pyarrow=True` and `memory_map=False`"
        raise ValueError(msg)

    with prepare_file_arg(
        source, use_pyarrow=use_pyarrow, storage_options=storage_options
    ) as data:
        if use_pyarrow:
            pyarrow_feather = import_optional(
                "pyarrow.feather",
                err_prefix="",
                err_suffix="is required when using 'read_ipc(..., use_pyarrow=True)'",
            )
            tbl = pyarrow_feather.read_table(
                data,
                memory_map=memory_map,
                columns=columns,
            )
            df = pl.DataFrame._from_arrow(tbl, rechunk=rechunk)
            if row_index_name is not None:
                df = df.with_row_index(row_index_name, row_index_offset)
            if n_rows is not None:
                df = df.slice(0, n_rows)
            return df

        return _read_ipc_impl(
            data,
            columns=columns,
            n_rows=n_rows,
            row_index_name=row_index_name,
            row_index_offset=row_index_offset,
            rechunk=rechunk,
            memory_map=memory_map,
        )


def _read_ipc_impl(
    source: str | Path | IO[bytes] | bytes,
    *,
    columns: Sequence[int] | Sequence[str] | None = None,
    n_rows: int | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    rechunk: bool = True,
    memory_map: bool = True,
) -> DataFrame:
    if isinstance(source, (str, Path)):
        source = normalize_filepath(source, check_not_directory=False)
    if isinstance(columns, str):
        columns = [columns]

    if isinstance(source, str) and is_glob_pattern(source) and is_local_file(source):
        scan = scan_ipc(
            source,
            n_rows=n_rows,
            rechunk=rechunk,
            row_index_name=row_index_name,
            row_index_offset=row_index_offset,
        )
        if columns is None:
            df = scan.collect()
        elif is_str_sequence(columns, allow_str=False):
            df = scan.select(columns).collect()
        else:
            msg = (
                "cannot use glob patterns and integer based projection as `columns` argument"
                "\n\nUse columns: List[str]"
            )
            raise TypeError(msg)
        return df

    projection, columns = parse_columns_arg(columns)
    pydf = PyDataFrame.read_ipc(
        source,
        columns,
        projection,
        n_rows,
        parse_row_index_args(row_index_name, row_index_offset),
        memory_map=memory_map,
    )
    return wrap_df(pydf)


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def read_ipc_stream(
    source: str | Path | IO[bytes] | bytes,
    *,
    columns: list[int] | list[str] | None = None,
    n_rows: int | None = None,
    use_pyarrow: bool = False,
    storage_options: dict[str, Any] | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    rechunk: bool = True,
) -> DataFrame:
    """
    Read into a DataFrame from Arrow IPC record batch stream.

    See "Streaming format" on https://arrow.apache.org/docs/python/ipc.html.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). If `fsspec` is installed, it might be used
        to open remote files. For file-like objects, the stream position may not be
        updated accordingly after reading.
    columns
        Columns to select. Accepts a list of column indices (starting at zero) or a list
        of column names.
    n_rows
        Stop reading from IPC stream after reading `n_rows`.
        Only valid when `use_pyarrow=False`.
    use_pyarrow
        Use pyarrow or the native Rust reader.
    storage_options
        Extra options that make sense for `fsspec.open()` or a particular storage
        connection, e.g. host, port, username, password, etc.
    row_index_name
        Insert a row index column with the given name into the DataFrame as the first
        column. If set to `None` (default), no row index column is created.
    row_index_offset
        Start the row index at this offset. Cannot be negative.
        Only used if `row_index_name` is set.
    rechunk
        Make sure that all data is contiguous.

    Returns
    -------
    DataFrame
    """
    with prepare_file_arg(
        source, use_pyarrow=use_pyarrow, storage_options=storage_options
    ) as data:
        if use_pyarrow:
            pyarrow_ipc = import_optional(
                "pyarrow.ipc",
                err_prefix="",
                err_suffix="is required when using 'read_ipc_stream(..., use_pyarrow=True)'",
            )
            with pyarrow_ipc.RecordBatchStreamReader(data) as reader:
                tbl = reader.read_all()
                df = pl.DataFrame._from_arrow(tbl, rechunk=rechunk)
                if row_index_name is not None:
                    df = df.with_row_index(row_index_name, row_index_offset)
                if n_rows is not None:
                    df = df.slice(0, n_rows)
                return df

        return _read_ipc_stream_impl(
            data,
            columns=columns,
            n_rows=n_rows,
            row_index_name=row_index_name,
            row_index_offset=row_index_offset,
            rechunk=rechunk,
        )


def _read_ipc_stream_impl(
    source: str | Path | IO[bytes] | bytes,
    *,
    columns: Sequence[int] | Sequence[str] | None = None,
    n_rows: int | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    rechunk: bool = True,
) -> DataFrame:
    if isinstance(source, (str, Path)):
        source = normalize_filepath(source, check_not_directory=False)
    if isinstance(columns, str):
        columns = [columns]

    projection, columns = parse_columns_arg(columns)
    pydf = PyDataFrame.read_ipc_stream(
        source,
        columns,
        projection,
        n_rows,
        parse_row_index_args(row_index_name, row_index_offset),
        rechunk,
    )
    return wrap_df(pydf)


def read_ipc_schema(source: str | Path | IO[bytes] | bytes) -> dict[str, DataType]:
    """
    Get the schema of an IPC file without reading data.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). For file-like objects, the stream position
        may not be updated accordingly after reading.

    Returns
    -------
    dict
        Dictionary mapping column names to datatypes
    """
    if isinstance(source, (str, Path)):
        source = normalize_filepath(source, check_not_directory=False)

    return _read_ipc_schema(source)


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def scan_ipc(
    source: (
        str
        | Path
        | IO[bytes]
        | bytes
        | list[str]
        | list[Path]
        | list[IO[bytes]]
        | list[bytes]
    ),
    *,
    n_rows: int | None = None,
    cache: bool = True,
    rechunk: bool = False,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    glob: bool = True,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    memory_map: bool = True,
    retries: int = 2,
    file_cache_ttl: int | None = None,
    hive_partitioning: bool | None = None,
    hive_schema: SchemaDict | None = None,
    try_parse_hive_dates: bool = True,
    include_file_paths: str | None = None,
) -> LazyFrame:
    """
    Lazily read from an Arrow IPC (Feather v2) file or multiple files via glob patterns.

    This allows the query optimizer to push down predicates and projections to the scan
    level, thereby potentially reducing memory overhead.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    Parameters
    ----------
    source
        Path(s) to a file or directory
        When needing to authenticate for scanning cloud locations, see the
        `storage_options` parameter.
    n_rows
        Stop reading from IPC file after reading `n_rows`.
    cache
        Cache the result after reading.
    rechunk
        Reallocate to contiguous memory when all chunks/ files are parsed.
    row_index_name
        If not None, this will insert a row index column with give name into the
        DataFrame
    row_index_offset
        Offset to start the row index column (only use if the name is set)
    glob
        Expand path given via globbing rules.
    storage_options
        Options that indicate how to connect to a cloud provider.

        The cloud providers currently supported are AWS, GCP, and Azure.
        See supported keys here:

        * `aws <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html>`_
        * `gcp <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html>`_
        * `azure <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html>`_
        * Hugging Face (`hf://`): Accepts an API key under the `token` parameter: \
          `{'token': '...'}`, or by setting the `HF_TOKEN` environment variable.

        If `storage_options` is not provided, Polars will try to infer the information
        from environment variables.
    credential_provider
        Provide a function that can be called to provide cloud storage
        credentials. The function is expected to return a dictionary of
        credential keys along with an optional credential expiry time.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

    memory_map
        Try to memory map the file. This can greatly improve performance on repeated
        queries as the OS may cache pages.
        Only uncompressed IPC files can be memory mapped.
    retries
        Number of retries if accessing a cloud instance fails.
    file_cache_ttl
        Amount of time to keep downloaded cloud files since their last access time,
        in seconds. Uses the `POLARS_FILE_CACHE_TTL` environment variable
        (which defaults to 1 hour) if not given.
    hive_partitioning
        Infer statistics and schema from Hive partitioned URL and use them
        to prune reads. This is unset by default (i.e. `None`), meaning it is
        automatically enabled when a single directory is passed, and otherwise
        disabled.
    hive_schema
        The column names and data types of the columns by which the data is partitioned.
        If set to `None` (default), the schema of the Hive partitions is inferred.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    try_parse_hive_dates
        Whether to try parsing hive values as date/datetime types.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
    """
    # Memory Mapping is now a no-op
    _ = memory_map

    sources = get_sources(source)

    credential_provider_builder = _init_credential_provider_builder(
        credential_provider, sources, storage_options, "scan_parquet"
    )
    del credential_provider

    pylf = PyLazyFrame.new_from_ipc(
        sources=sources,
        scan_options=ScanOptions(
            row_index=(
                (row_index_name, row_index_offset)
                if row_index_name is not None
                else None
            ),
            pre_slice=(0, n_rows) if n_rows is not None else None,
            include_file_paths=include_file_paths,
            glob=glob,
            hive_partitioning=hive_partitioning,
            hive_schema=hive_schema,
            try_parse_hive_dates=try_parse_hive_dates,
            rechunk=rechunk,
            cache=cache,
            storage_options=(
                list(storage_options.items()) if storage_options is not None else None
            ),
            credential_provider=credential_provider_builder,
            retries=retries,
        ),
        file_cache_ttl=file_cache_ttl,
    )

    return wrap_ldf(pylf)
