from __future__ import annotations

import contextlib
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Literal

from polars._utils.deprecation import deprecate_renamed_parameter
from polars._utils.various import is_path_or_str_sequence, normalize_filepath
from polars._utils.wrap import wrap_ldf
from polars.datatypes import N_INFER_DEFAULT
from polars.io._utils import parse_row_index_args
from polars.io.cloud.credential_provider._builder import (
    _init_credential_provider_builder,
)

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyLazyFrame

if TYPE_CHECKING:
    from polars import DataFrame, LazyFrame
    from polars._typing import SchemaDefinition
    from polars.io.cloud import CredentialProviderFunction


def read_ndjson(
    source: str
    | Path
    | IO[str]
    | IO[bytes]
    | bytes
    | list[str]
    | list[Path]
    | list[IO[str]]
    | list[IO[bytes]],
    *,
    schema: SchemaDefinition | None = None,
    schema_overrides: SchemaDefinition | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    batch_size: int | None = 1024,
    n_rows: int | None = None,
    low_memory: bool = False,
    rechunk: bool = False,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    ignore_errors: bool = False,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    retries: int = 2,
    file_cache_ttl: int | None = None,
    include_file_paths: str | None = None,
) -> DataFrame:
    r"""
    Read into a DataFrame from a newline delimited JSON file.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). For file-like objects, the stream position
        may not be updated accordingly after reading.
    schema : Sequence of str, (str,DataType) pairs, or a {str:DataType,} dict
        The DataFrame schema may be declared in several ways:

        * As a dict of {name:type} pairs; if type is None, it will be auto-inferred.
        * As a list of column names; in this case types are automatically inferred.
        * As a list of (name,type) pairs; this is equivalent to the dictionary form.

        If you supply a list of column names that does not match the names in the
        underlying data, the names given here will overwrite them. The number
        of names given in the schema should match the underlying data dimensions.
    schema_overrides : dict, default None
        Support type specification or override of one or more columns; note that
        any dtypes inferred from the schema param will be overridden.
    infer_schema_length
        The maximum number of rows to scan for schema inference.
        If set to `None`, the full data may be scanned *(this is slow)*.
    batch_size
        Number of rows to read in each batch.
    n_rows
        Stop reading from JSON file after reading `n_rows`.
    low_memory
        Reduce memory pressure at the expense of performance.
    rechunk
        Reallocate to contiguous memory when all chunks/ files are parsed.
    row_index_name
        If not None, this will insert a row index column with give name into the
        DataFrame
    row_index_offset
        Offset to start the row index column (only use if the name is set)
    ignore_errors
        Return `Null` if parsing fails because of schema mismatches.
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
    file_cache_ttl
        Amount of time to keep downloaded cloud files since their last access time,
        in seconds. Uses the `POLARS_FILE_CACHE_TTL` environment variable
        (which defaults to 1 hour) if not given.
    include_file_paths
        Include the path of the source file(s) as a column with this name.

    See Also
    --------
    scan_ndjson : Lazily read from an NDJSON file or multiple files via glob patterns.

    Warnings
    --------
    Calling `read_ndjson().lazy()` is an antipattern as this forces Polars to
    materialize a full ndjson file and therefore cannot push any optimizations into
    the reader. Therefore always prefer `scan_ndjson` if you want to work with
    `LazyFrame` s.

    Examples
    --------
    >>> from io import StringIO
    >>> json_str = '{"foo":1,"bar":6}\n{"foo":2,"bar":7}\n{"foo":3,"bar":8}\n'
    >>> pl.read_ndjson(StringIO(json_str))
    shape: (3, 2)
    ┌─────┬─────┐
    │ foo ┆ bar │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 6   │
    │ 2   ┆ 7   │
    │ 3   ┆ 8   │
    └─────┴─────┘
    """
    credential_provider_builder = _init_credential_provider_builder(
        credential_provider, source, storage_options, "read_ndjson"
    )

    del credential_provider

    return scan_ndjson(
        source,
        schema=schema,
        schema_overrides=schema_overrides,
        infer_schema_length=infer_schema_length,
        batch_size=batch_size,
        n_rows=n_rows,
        low_memory=low_memory,
        rechunk=rechunk,
        row_index_name=row_index_name,
        row_index_offset=row_index_offset,
        ignore_errors=ignore_errors,
        include_file_paths=include_file_paths,
        retries=retries,
        storage_options=storage_options,
        credential_provider=credential_provider_builder,  # type: ignore[arg-type]
        file_cache_ttl=file_cache_ttl,
    ).collect()


@deprecate_renamed_parameter("row_count_name", "row_index_name", version="0.20.4")
@deprecate_renamed_parameter("row_count_offset", "row_index_offset", version="0.20.4")
def scan_ndjson(
    source: (
        str
        | Path
        | IO[str]
        | IO[bytes]
        | bytes
        | list[str]
        | list[Path]
        | list[IO[str]]
        | list[IO[bytes]]
    ),
    *,
    schema: SchemaDefinition | None = None,
    schema_overrides: SchemaDefinition | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    batch_size: int | None = 1024,
    n_rows: int | None = None,
    low_memory: bool = False,
    rechunk: bool = False,
    row_index_name: str | None = None,
    row_index_offset: int = 0,
    ignore_errors: bool = False,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    retries: int = 2,
    file_cache_ttl: int | None = None,
    include_file_paths: str | None = None,
) -> LazyFrame:
    """
    Lazily read from a newline delimited JSON file or multiple files via glob patterns.

    This allows the query optimizer to push down predicates and projections to the scan
    level, thereby potentially reducing memory overhead.

    .. versionchanged:: 0.20.4
        * The `row_count_name` parameter was renamed `row_index_name`.
        * The `row_count_offset` parameter was renamed `row_index_offset`.

    Parameters
    ----------
    source
        Path to a file.
    schema : Sequence of str, (str,DataType) pairs, or a {str:DataType,} dict
        The DataFrame schema may be declared in several ways:

        * As a dict of {name:type} pairs; if type is None, it will be auto-inferred.
        * As a list of column names; in this case types are automatically inferred.
        * As a list of (name,type) pairs; this is equivalent to the dictionary form.

        If you supply a list of column names that does not match the names in the
        underlying data, the names given here will overwrite them. The number
        of names given in the schema should match the underlying data dimensions.
    schema_overrides : dict, default None
        Support type specification or override of one or more columns; note that
        any dtypes inferred from the schema param will be overridden.
    infer_schema_length
        The maximum number of rows to scan for schema inference.
        If set to `None`, the full data may be scanned *(this is slow)*.
    batch_size
        Number of rows to read in each batch.
    n_rows
        Stop reading from JSON file after reading `n_rows`.
    low_memory
        Reduce memory pressure at the expense of performance.
    rechunk
        Reallocate to contiguous memory when all chunks/ files are parsed.
    row_index_name
        If not None, this will insert a row index column with give name into the
        DataFrame
    row_index_offset
        Offset to start the row index column (only use if the name is set)
    ignore_errors
        Return `Null` if parsing fails because of schema mismatches.
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
    file_cache_ttl
        Amount of time to keep downloaded cloud files since their last access time,
        in seconds. Uses the `POLARS_FILE_CACHE_TTL` environment variable
        (which defaults to 1 hour) if not given.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
    """
    sources: list[str] | list[Path] | list[IO[str]] | list[IO[bytes]] = []
    if isinstance(source, (str, Path)):
        source = normalize_filepath(source, check_not_directory=False)
    elif isinstance(source, list):
        if is_path_or_str_sequence(source):
            sources = [
                normalize_filepath(source, check_not_directory=False)
                for source in source
            ]
        else:
            sources = source

        source = None  # type: ignore[assignment]

    if infer_schema_length == 0:
        msg = "'infer_schema_length' should be positive"
        raise ValueError(msg)

    credential_provider_builder = _init_credential_provider_builder(
        credential_provider, source, storage_options, "scan_ndjson"
    )

    del credential_provider

    if storage_options:
        storage_options = list(storage_options.items())  # type: ignore[assignment]
    else:
        # Handle empty dict input
        storage_options = None

    pylf = PyLazyFrame.new_from_ndjson(
        source,
        sources,
        infer_schema_length=infer_schema_length,
        schema=schema,
        schema_overrides=schema_overrides,
        batch_size=batch_size,
        n_rows=n_rows,
        low_memory=low_memory,
        rechunk=rechunk,
        row_index=parse_row_index_args(row_index_name, row_index_offset),
        ignore_errors=ignore_errors,
        include_file_paths=include_file_paths,
        retries=retries,
        cloud_options=storage_options,
        credential_provider=credential_provider_builder,
        file_cache_ttl=file_cache_ttl,
    )
    return wrap_ldf(pylf)
