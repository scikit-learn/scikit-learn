from __future__ import annotations

import contextlib
import io
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any

import polars.functions as F
from polars import concat as plconcat
from polars._dependencies import import_optional
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    issue_deprecation_warning,
)
from polars._utils.unstable import issue_unstable_warning
from polars._utils.various import (
    is_int_sequence,
    normalize_filepath,
)
from polars._utils.wrap import wrap_ldf
from polars.convert import from_arrow
from polars.io._utils import (
    get_sources,
    prepare_file_arg,
)
from polars.io.cloud.credential_provider._builder import (
    _init_credential_provider_builder,
)
from polars.io.scan_options._options import ScanOptions

with contextlib.suppress(ImportError):
    from polars._plr import PyLazyFrame
    from polars._plr import read_parquet_metadata as _read_parquet_metadata

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from polars import DataFrame, DataType, LazyFrame
    from polars._typing import (
        ColumnMapping,
        DefaultFieldValues,
        DeletionFiles,
        FileSource,
        ParallelStrategy,
        SchemaDict,
    )
    from polars.io.cloud import CredentialProviderFunction
    from polars.io.scan_options import ScanCastOptions


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def read_parquet(
    source: FileSource,
    *,
    columns: list[int] | list[str] | None = None,
    n_rows: int | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    parallel: ParallelStrategy = "auto",
    use_statistics: bool = True,
    hive_partitioning: bool | None = None,
    glob: bool = True,
    schema: SchemaDict | None = None,
    hive_schema: SchemaDict | None = None,
    try_parse_hive_dates: bool = True,
    rechunk: bool = False,
    low_memory: bool = False,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    retries: int = 2,
    use_pyarrow: bool = False,
    pyarrow_options: dict[str, Any] | None = None,
    memory_map: bool = True,
    include_file_paths: str | None = None,
    missing_columns: Literal["insert", "raise"] = "raise",
    allow_missing_columns: bool | None = None,
) -> DataFrame:
    """
    Read into a DataFrame from a parquet file.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    Parameters
    ----------
    source
        Path(s) to a file or directory
        When needing to authenticate for scanning cloud locations, see the
        `storage_options` parameter.

        File-like objects are supported (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). For file-like objects, the stream position
        may not be updated accordingly after reading.
    columns
        Columns to select. Accepts a list of column indices (starting at zero) or a list
        of column names.
    n_rows
        Stop reading from parquet file after reading `n_rows`.
        Only valid when `use_pyarrow=False`.
    row_index_name
        Insert a row index column with the given name into the DataFrame as the first
        column. If set to `None` (default), no row index column is created.
    row_index_offset
        Start the row index at this offset. Cannot be negative.
        Only used if `row_index_name` is set.
    parallel : {'auto', 'columns', 'row_groups', 'none'}
        This determines the direction of parallelism. 'auto' will try to determine the
        optimal direction.
    use_statistics
        Use statistics in the parquet to determine if pages
        can be skipped from reading.
    hive_partitioning
        Infer statistics and schema from Hive partitioned URL and use them
        to prune reads. This is unset by default (i.e. `None`), meaning it is
        automatically enabled when a single directory is passed, and otherwise
        disabled.
    glob
        Expand path given via globbing rules.
    schema
        Specify the datatypes of the columns. The datatypes must match the
        datatypes in the file(s). If there are extra columns that are not in the
        file(s), consider also passing `missing_columns='insert'`.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    hive_schema
        The column names and data types of the columns by which the data is partitioned.
        If set to `None` (default), the schema of the Hive partitions is inferred.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    try_parse_hive_dates
        Whether to try parsing hive values as date/datetime types.
    rechunk
        Make sure that all columns are contiguous in memory by
        aggregating the chunks into a single array.
    low_memory
        Reduce memory pressure at the expense of performance.
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
    retries
        Number of retries if accessing a cloud instance fails.
    use_pyarrow
        Use PyArrow instead of the Rust-native Parquet reader. The PyArrow reader is
        more stable.
    pyarrow_options
        Keyword arguments for `pyarrow.parquet.read_table
        <https://arrow.apache.org/docs/python/generated/pyarrow.parquet.read_table.html>`_.
    memory_map
        Memory map underlying file. This will likely increase performance.
        Only used when `use_pyarrow=True`.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
        Only valid when `use_pyarrow=False`.
    missing_columns
        Configuration for behavior when columns defined in the schema
        are missing from the data:

        * `insert`: Inserts the missing columns using NULLs as the row values.
        * `raise`: Raises an error.

    allow_missing_columns
        When reading a list of parquet files, if a column existing in the first
        file cannot be found in subsequent files, the default behavior is to
        raise an error. However, if `allow_missing_columns` is set to
        `True`, a full-NULL column is returned instead of erroring for the files
        that do not contain the column.

        .. deprecated:: 1.30.0
            Use the parameter `missing_columns` instead and pass one of
            `('insert', 'raise')`.

    Returns
    -------
    DataFrame

    See Also
    --------
    scan_parquet: Lazily read from a parquet file or multiple files via glob patterns.
    scan_pyarrow_dataset

    Warnings
    --------
    Calling `read_parquet().lazy()` is an antipattern as this forces Polars to
    materialize a full parquet file and therefore cannot push any optimizations
    into the reader. Therefore always prefer `scan_parquet` if you want to work
    with `LazyFrame` s.

    """
    if schema is not None:
        msg = "the `schema` parameter of `read_parquet` is considered unstable."
        issue_unstable_warning(msg)

    if hive_schema is not None:
        msg = "the `hive_schema` parameter of `read_parquet` is considered unstable."
        issue_unstable_warning(msg)

    # Dispatch to pyarrow if requested
    if use_pyarrow:
        if n_rows is not None:
            msg = "`n_rows` cannot be used with `use_pyarrow=True`"
            raise ValueError(msg)
        if include_file_paths is not None:
            msg = "`include_file_paths` cannot be used with `use_pyarrow=True`"
            raise ValueError(msg)
        if schema is not None:
            msg = "`schema` cannot be used with `use_pyarrow=True`"
            raise ValueError(msg)
        if hive_schema is not None:
            msg = (
                "cannot use `hive_partitions` with `use_pyarrow=True`"
                "\n\nHint: Pass `pyarrow_options` instead with a 'partitioning' entry."
            )
            raise TypeError(msg)
        return _read_parquet_with_pyarrow(
            source,
            columns=columns,
            storage_options=storage_options,
            pyarrow_options=pyarrow_options,
            memory_map=memory_map,
            rechunk=rechunk,
        )

    if allow_missing_columns is not None:
        issue_deprecation_warning(
            "the parameter `allow_missing_columns` for `read_parquet` is deprecated. "
            "Use the parameter `missing_columns` instead and pass one of "
            "`('insert', 'raise')`.",
            version="1.30.0",
        )

        missing_columns = "insert" if allow_missing_columns else "raise"

    # For other inputs, defer to `scan_parquet`
    lf = scan_parquet(
        source,
        n_rows=n_rows,
        row_index_name=row_index_name,
        row_index_offset=row_index_offset,
        parallel=parallel,
        use_statistics=use_statistics,
        hive_partitioning=hive_partitioning,
        schema=schema,
        hive_schema=hive_schema,
        try_parse_hive_dates=try_parse_hive_dates,
        rechunk=rechunk,
        low_memory=low_memory,
        cache=False,
        storage_options=storage_options,
        credential_provider=credential_provider,
        retries=retries,
        glob=glob,
        include_file_paths=include_file_paths,
        missing_columns=missing_columns,
    )

    if columns is not None:
        if is_int_sequence(columns):
            lf = lf.select(F.nth(columns))
        else:
            lf = lf.select(columns)

    return lf.collect()


def _read_parquet_with_pyarrow(
    source: str
    | Path
    | IO[bytes]
    | bytes
    | list[str]
    | list[Path]
    | list[IO[bytes]]
    | list[bytes],
    *,
    columns: list[int] | list[str] | None = None,
    storage_options: dict[str, Any] | None = None,
    pyarrow_options: dict[str, Any] | None = None,
    memory_map: bool = True,
    rechunk: bool = True,
) -> DataFrame:
    pyarrow_parquet = import_optional(
        "pyarrow.parquet",
        err_prefix="",
        err_suffix="is required when using `read_parquet(..., use_pyarrow=True)`",
    )
    pyarrow_options = pyarrow_options or {}

    sources: list[str | Path | IO[bytes] | bytes | list[str] | list[Path]] = []
    if isinstance(source, list):
        if len(source) > 0 and isinstance(source[0], (bytes, io.IOBase)):
            sources = source  # type: ignore[assignment]
        else:
            sources = [source]  # type: ignore[list-item]
    else:
        sources = [source]

    results: list[DataFrame] = []
    for source in sources:
        with prepare_file_arg(
            source,  # type: ignore[arg-type]
            use_pyarrow=True,
            storage_options=storage_options,
        ) as source_prep:
            pa_table = pyarrow_parquet.read_table(
                source_prep,
                memory_map=memory_map,
                columns=columns,
                **pyarrow_options,
            )
        result = from_arrow(pa_table, rechunk=rechunk)
        results.append(result)  # type: ignore[arg-type]

    if len(results) == 1:
        return results[0]
    else:
        return plconcat(results)


def read_parquet_schema(source: str | Path | IO[bytes] | bytes) -> dict[str, DataType]:
    """
    Get the schema of a Parquet file without reading data.

    If you would like to read the schema of a cloud file with authentication
    configuration, it is recommended use `scan_parquet` - e.g.
    `scan_parquet(..., storage_options=...).collect_schema()`.

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

    See Also
    --------
    scan_parquet
    """
    return scan_parquet(source).collect_schema()


def read_parquet_metadata(
    source: str | Path | IO[bytes] | bytes,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    retries: int = 2,
) -> dict[str, str]:
    """
    Get file-level custom metadata of a Parquet file without reading data.

    .. warning::
        This functionality is considered **experimental**. It may be removed or
        changed at any point without it being considered a breaking change.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). For file-like objects, the stream position
        may not be updated accordingly after reading.
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
    retries
        Number of retries if accessing a cloud instance fails.

    Returns
    -------
    dict
        Dictionary with the metadata. Empty if no custom metadata is available.
    """
    if isinstance(source, (str, Path)):
        source = normalize_filepath(source, check_not_directory=False)

    credential_provider_builder = _init_credential_provider_builder(
        credential_provider, source, storage_options, "scan_parquet"
    )
    del credential_provider

    return _read_parquet_metadata(
        source,
        storage_options=(
            list(storage_options.items()) if storage_options is not None else None
        ),
        credential_provider=credential_provider_builder,
        retries=retries,
    )


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def scan_parquet(
    source: FileSource,
    *,
    n_rows: int | None = None,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    parallel: ParallelStrategy = "auto",
    use_statistics: bool = True,
    hive_partitioning: bool | None = None,
    glob: bool = True,
    hidden_file_prefix: str | Sequence[str] | None = None,
    schema: SchemaDict | None = None,
    hive_schema: SchemaDict | None = None,
    try_parse_hive_dates: bool = True,
    rechunk: bool = False,
    low_memory: bool = False,
    cache: bool = True,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    retries: int = 2,
    include_file_paths: str | None = None,
    missing_columns: Literal["insert", "raise"] = "raise",
    allow_missing_columns: bool | None = None,
    extra_columns: Literal["ignore", "raise"] = "raise",
    cast_options: ScanCastOptions | None = None,
    _column_mapping: ColumnMapping | None = None,
    _default_values: DefaultFieldValues | None = None,
    _deletion_files: DeletionFiles | None = None,
    _table_statistics: DataFrame | None = None,
    _row_count: tuple[int, int] | None = None,
) -> LazyFrame:
    """
    Lazily read from a local or cloud-hosted parquet file (or files).

    This function allows the query optimizer to push down predicates and projections to
    the scan level, typically increasing performance and reducing memory overhead.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    .. versionchanged:: 1.30.0
        * The `allow_missing_columns` is deprecated in favor of `missing_columns`.

    Parameters
    ----------
    source
        Path(s) to a file or directory
        When needing to authenticate for scanning cloud locations, see the
        `storage_options` parameter.
    n_rows
        Stop reading from parquet file after reading `n_rows`.
    row_index_name
        If not None, this will insert a row index column with the given name into the
        DataFrame
    row_index_offset
        Offset to start the row index column (only used if the name is set)
    parallel : {'auto', 'columns', 'row_groups', 'prefiltered', 'none'}
        This determines the direction and strategy of parallelism. 'auto' will
        try to determine the optimal direction.

        The `prefiltered` strategy first evaluates the pushed-down predicates in
        parallel and determines a mask of which rows to read. Then, it
        parallelizes over both the columns and the row groups while filtering
        out rows that do not need to be read. This can provide significant
        speedups for large files (i.e. many row-groups) with a predicate that
        filters clustered rows or filters heavily. In other cases,
        `prefiltered` may slow down the scan compared other strategies.

        The `prefiltered` settings falls back to `auto` if no predicate is
        given.

        .. warning::
            The `prefiltered` strategy is considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

    use_statistics
        Use statistics in the parquet to determine if pages
        can be skipped from reading.
    hive_partitioning
        Infer statistics and schema from hive partitioned URL and use them
        to prune reads.
    glob
        Expand path given via globbing rules.
    hidden_file_prefix
        Skip reading files whose names begin with the specified prefixes.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    schema
        Specify the datatypes of the columns. The datatypes must match the
        datatypes in the file(s). If there are extra columns that are not in the
        file(s), consider also passing `missing_columns='insert'`.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    hive_schema
        The column names and data types of the columns by which the data is partitioned.
        If set to `None` (default), the schema of the Hive partitions is inferred.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    try_parse_hive_dates
        Whether to try parsing hive values as date/datetime types.
    rechunk
        In case of reading multiple files via a glob pattern rechunk the final DataFrame
        into contiguous memory chunks.
    low_memory
        Reduce memory pressure at the expense of performance.
    cache
        Cache the result after reading.
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
    retries
        Number of retries if accessing a cloud instance fails.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
    missing_columns
        Configuration for behavior when columns defined in the schema
        are missing from the data:

        * `insert`: Inserts the missing columns using NULLs as the row values.
        * `raise`: Raises an error.

    allow_missing_columns
        When reading a list of parquet files, if a column existing in the first
        file cannot be found in subsequent files, the default behavior is to
        raise an error. However, if `allow_missing_columns` is set to
        `True`, a full-NULL column is returned instead of erroring for the files
        that do not contain the column.

        .. deprecated:: 1.30.0
            Use the parameter `missing_columns` instead and pass one of
            `('insert', 'raise')`.
    extra_columns
        Configuration for behavior when extra columns outside of the
        defined schema are encountered in the data:

        * `ignore`: Silently ignores.
        * `raise`: Raises an error.

    cast_options
        Configuration for column type-casting during scans. Useful for datasets
        containing files that have differing schemas.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

    See Also
    --------
    read_parquet
    scan_pyarrow_dataset

    Examples
    --------
    Scan a local Parquet file.

    >>> pl.scan_parquet("path/to/file.parquet")  # doctest: +SKIP

    Scan a file on AWS S3.

    >>> source = "s3://bucket/*.parquet"
    >>> pl.scan_parquet(source)  # doctest: +SKIP
    >>> storage_options = {
    ...     "aws_access_key_id": "<secret>",
    ...     "aws_secret_access_key": "<secret>",
    ...     "aws_region": "us-east-1",
    ... }
    >>> pl.scan_parquet(source, storage_options=storage_options)  # doctest: +SKIP
    """
    if schema is not None:
        msg = "the `schema` parameter of `scan_parquet` is considered unstable."
        issue_unstable_warning(msg)

    if hive_schema is not None:
        msg = "the `hive_schema` parameter of `scan_parquet` is considered unstable."
        issue_unstable_warning(msg)

    if cast_options is not None:
        msg = "The `cast_options` parameter of `scan_parquet` is considered unstable."
        issue_unstable_warning(msg)

    if hidden_file_prefix is not None:
        msg = "The `hidden_file_prefix` parameter of `scan_parquet` is considered unstable."
        issue_unstable_warning(msg)

    if allow_missing_columns is not None:
        issue_deprecation_warning(
            "the parameter `allow_missing_columns` for `scan_parquet` is deprecated. "
            "Use the parameter `missing_columns` instead and pass one of "
            "`('insert', 'raise')`.",
            version="1.30.0",
        )

        missing_columns = "insert" if allow_missing_columns else "raise"

    sources = get_sources(source)

    credential_provider_builder = _init_credential_provider_builder(
        credential_provider, sources, storage_options, "scan_parquet"
    )

    del credential_provider

    pylf = PyLazyFrame.new_from_parquet(
        sources=sources,
        schema=schema,
        parallel=parallel,
        low_memory=low_memory,
        use_statistics=use_statistics,
        scan_options=ScanOptions(
            row_index=(
                (row_index_name, row_index_offset)
                if row_index_name is not None
                else None
            ),
            pre_slice=(0, n_rows) if n_rows is not None else None,
            cast_options=cast_options,
            extra_columns=extra_columns,
            missing_columns=missing_columns,
            include_file_paths=include_file_paths,
            glob=glob,
            hidden_file_prefix=(
                [hidden_file_prefix]
                if isinstance(hidden_file_prefix, str)
                else hidden_file_prefix
            ),
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
            column_mapping=_column_mapping,
            default_values=_default_values,
            deletion_files=_deletion_files,
            table_statistics=_table_statistics,
            row_count=_row_count,
        ),
    )

    return wrap_ldf(pylf)
