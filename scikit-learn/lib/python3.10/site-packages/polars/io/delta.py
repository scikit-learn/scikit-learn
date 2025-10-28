from __future__ import annotations

import warnings
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from polars._dependencies import _DELTALAKE_AVAILABLE, deltalake
from polars.datatypes import Null, Time
from polars.datatypes.convert import unpack_dtypes
from polars.io.cloud._utils import _get_path_scheme
from polars.io.parquet import scan_parquet
from polars.io.pyarrow_dataset.functions import scan_pyarrow_dataset
from polars.io.scan_options.cast_options import ScanCastOptions
from polars.schema import Schema

if TYPE_CHECKING:
    from typing import Literal

    from deltalake import DeltaTable

    from polars import DataFrame, DataType, LazyFrame
    from polars.io.cloud import CredentialProviderFunction


def read_delta(
    source: str | Path | DeltaTable,
    *,
    version: int | str | datetime | None = None,
    columns: list[str] | None = None,
    rechunk: bool | None = None,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    delta_table_options: dict[str, Any] | None = None,
    use_pyarrow: bool = False,
    pyarrow_options: dict[str, Any] | None = None,
) -> DataFrame:
    """
    Reads into a DataFrame from a Delta lake table.

    Parameters
    ----------
    source
        DeltaTable or a Path or URI to the root of the Delta lake table.

        Note: For Local filesystem, absolute and relative paths are supported but
        for the supported object storages - GCS, Azure and S3 full URI must be provided.
    version
        Numerical version or timestamp version of the Delta lake table.

        Note: If `version` is not provided, the latest version of delta lake
        table is read.
    columns
        Columns to select. Accepts a list of column names.
    rechunk
        Make sure that all columns are contiguous in memory by
        aggregating the chunks into a single array.
    storage_options
        Extra options for the storage backends supported by `deltalake`.
        For cloud storages, this may include configurations for authentication etc.

        More info is available `here
        <https://delta-io.github.io/delta-rs/usage/loading-table/>`__.
    credential_provider
        Provide a function that can be called to provide cloud storage
        credentials. The function is expected to return a dictionary of
        credential keys along with an optional credential expiry time.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    delta_table_options
        Additional keyword arguments while reading a Delta lake Table.
    use_pyarrow
        Flag to enable pyarrow dataset reads.
    pyarrow_options
        Keyword arguments while converting a Delta lake Table to pyarrow table.

    Returns
    -------
    DataFrame

    Examples
    --------
    Reads a Delta table from local filesystem.
    Note: Since version is not provided, the latest version of the delta table is read.

    >>> table_path = "/path/to/delta-table/"
    >>> pl.read_delta(table_path)  # doctest: +SKIP

    Reads a specific version of the Delta table from local filesystem.
    Note: This will fail if the provided version of the delta table does not exist.

    >>> pl.read_delta(table_path, version=1)  # doctest: +SKIP

    Time travel a delta table from local filesystem using a timestamp version.

    >>> pl.read_delta(
    ...     table_path, version=datetime(2020, 1, 1, tzinfo=timezone.utc)
    ... )  # doctest: +SKIP

    Reads a Delta table from AWS S3.
    See a list of supported storage options for S3 `here
    <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html#variants>`__.

    >>> table_path = "s3://bucket/path/to/delta-table/"
    >>> storage_options = {
    ...     "AWS_ACCESS_KEY_ID": "THE_AWS_ACCESS_KEY_ID",
    ...     "AWS_SECRET_ACCESS_KEY": "THE_AWS_SECRET_ACCESS_KEY",
    ... }
    >>> pl.read_delta(table_path, storage_options=storage_options)  # doctest: +SKIP

    Reads a Delta table from Google Cloud storage (GCS).
    See a list of supported storage options for GCS `here
    <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html#variants>`__.

    >>> table_path = "gs://bucket/path/to/delta-table/"
    >>> storage_options = {"SERVICE_ACCOUNT": "SERVICE_ACCOUNT_JSON_ABSOLUTE_PATH"}
    >>> pl.read_delta(table_path, storage_options=storage_options)  # doctest: +SKIP

    Reads a Delta table from Azure.

    Following type of table paths are supported,

    * az://<container>/<path>
    * adl://<container>/<path>
    * abfs://<container>/<path>

    See a list of supported storage options for Azure `here
    <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html#variants>`__.

    >>> table_path = "az://container/path/to/delta-table/"
    >>> storage_options = {
    ...     "AZURE_STORAGE_ACCOUNT_NAME": "AZURE_STORAGE_ACCOUNT_NAME",
    ...     "AZURE_STORAGE_ACCOUNT_KEY": "AZURE_STORAGE_ACCOUNT_KEY",
    ... }
    >>> pl.read_delta(table_path, storage_options=storage_options)  # doctest: +SKIP

    Reads a Delta table with additional delta specific options. In the below example,
    `without_files` option is used which loads the table without file tracking
    information.

    >>> table_path = "/path/to/delta-table/"
    >>> delta_table_options = {"without_files": True}
    >>> pl.read_delta(
    ...     table_path, delta_table_options=delta_table_options
    ... )  # doctest: +SKIP
    """
    df = scan_delta(
        source=source,
        version=version,
        storage_options=storage_options,
        credential_provider=credential_provider,
        delta_table_options=delta_table_options,
        use_pyarrow=use_pyarrow,
        pyarrow_options=pyarrow_options,
        rechunk=rechunk,
    )

    if columns is not None:
        df = df.select(columns)
    return df.collect()


def scan_delta(
    source: str | Path | DeltaTable,
    *,
    version: int | str | datetime | None = None,
    storage_options: dict[str, Any] | None = None,
    credential_provider: CredentialProviderFunction | Literal["auto"] | None = "auto",
    delta_table_options: dict[str, Any] | None = None,
    use_pyarrow: bool = False,
    pyarrow_options: dict[str, Any] | None = None,
    rechunk: bool | None = None,
) -> LazyFrame:
    """
    Lazily read from a Delta lake table.

    Parameters
    ----------
    source
        DeltaTable or a Path or URI to the root of the Delta lake table.

        Note: For Local filesystem, absolute and relative paths are supported but
        for the supported object storages - GCS, Azure and S3 full URI must be provided.
    version
        Numerical version or timestamp version of the Delta lake table.

        Note: If `version` is not provided, the latest version of delta lake
        table is read.
    storage_options
        Extra options for the storage backends supported by `deltalake`.
        For cloud storages, this may include configurations for authentication etc.

        More info is available `here
        <https://delta-io.github.io/delta-rs/usage/loading-table/>`__.
    credential_provider
        Provide a function that can be called to provide cloud storage
        credentials. The function is expected to return a dictionary of
        credential keys along with an optional credential expiry time.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    delta_table_options
        Additional keyword arguments while reading a Delta lake Table.
    use_pyarrow
        Flag to enable pyarrow dataset reads.
    pyarrow_options
        Keyword arguments while converting a Delta lake Table to pyarrow table.
        Use this parameter when filtering on partitioned columns or to read
        from a 'fsspec' supported filesystem.
    rechunk
        Make sure that all columns are contiguous in memory by
        aggregating the chunks into a single array.

    Returns
    -------
    LazyFrame

    Examples
    --------
    Creates a scan for a Delta table from local filesystem.
    Note: Since version is not provided, the latest version of the delta table is read.

    >>> table_path = "/path/to/delta-table/"
    >>> pl.scan_delta(table_path).collect()  # doctest: +SKIP

    Creates a scan for a specific version of the Delta table from local filesystem.
    Note: This will fail if the provided version of the delta table does not exist.

    >>> pl.scan_delta(table_path, version=1).collect()  # doctest: +SKIP

    Time travel a delta table from local filesystem using a timestamp version.

    >>> pl.scan_delta(
    ...     table_path, version=datetime(2020, 1, 1, tzinfo=timezone.utc)
    ... ).collect()  # doctest: +SKIP

    Creates a scan for a Delta table from AWS S3.
    See a list of supported storage options for S3 `here
    <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html#variants>`__.

    >>> table_path = "s3://bucket/path/to/delta-table/"
    >>> storage_options = {
    ...     "AWS_REGION": "eu-central-1",
    ...     "AWS_ACCESS_KEY_ID": "THE_AWS_ACCESS_KEY_ID",
    ...     "AWS_SECRET_ACCESS_KEY": "THE_AWS_SECRET_ACCESS_KEY",
    ... }
    >>> pl.scan_delta(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for a Delta table from Google Cloud storage (GCS).
    See a list of supported storage options for GCS `here
    <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html#variants>`__.

    >>> table_path = "gs://bucket/path/to/delta-table/"
    >>> storage_options = {"SERVICE_ACCOUNT": "SERVICE_ACCOUNT_JSON_ABSOLUTE_PATH"}
    >>> pl.scan_delta(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for a Delta table from Azure.
    Supported options for Azure are available `here
    <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html#variants>`__.

    Following type of table paths are supported,

    * az://<container>/<path>
    * adl://<container>/<path>
    * abfs[s]://<container>/<path>

    >>> table_path = "az://container/path/to/delta-table/"
    >>> storage_options = {
    ...     "AZURE_STORAGE_ACCOUNT_NAME": "AZURE_STORAGE_ACCOUNT_NAME",
    ...     "AZURE_STORAGE_ACCOUNT_KEY": "AZURE_STORAGE_ACCOUNT_KEY",
    ... }
    >>> pl.scan_delta(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for a Delta table with additional delta specific options.
    In the below example, `without_files` option is used which loads the table without
    file tracking information.

    >>> table_path = "/path/to/delta-table/"
    >>> delta_table_options = {"without_files": True}
    >>> pl.scan_delta(
    ...     table_path, delta_table_options=delta_table_options
    ... ).collect()  # doctest: +SKIP
    """
    _check_if_delta_available()

    credential_provider_creds = {}

    from deltalake import DeltaTable

    from polars.io.cloud.credential_provider._builder import (
        _init_credential_provider_builder,
    )
    from polars.io.cloud.credential_provider._providers import (
        _get_credentials_from_provider_expiry_aware,
    )

    if not isinstance(source, DeltaTable):
        credential_provider_builder = _init_credential_provider_builder(
            credential_provider, source, storage_options, "scan_delta"
        )
    elif credential_provider is not None and credential_provider != "auto":
        msg = "cannot use credential_provider when passing a DeltaTable object"
        raise ValueError(msg)
    else:
        credential_provider_builder = None

    del credential_provider

    if credential_provider_builder and (
        provider := credential_provider_builder.build_credential_provider()
    ):
        credential_provider_creds = (
            _get_credentials_from_provider_expiry_aware(provider) or {}
        )

    dl_tbl = _get_delta_lake_table(
        table_path=source,
        version=version,
        storage_options=(
            {**(storage_options or {}), **credential_provider_creds}
            if storage_options is not None or credential_provider_builder is not None
            else None
        ),
        delta_table_options=delta_table_options,
    )

    if isinstance(source, DeltaTable) and (
        source._storage_options is not None or storage_options is not None
    ):
        storage_options = {**(source._storage_options or {}), **(storage_options or {})}

    if use_pyarrow:
        pyarrow_options = pyarrow_options or {}
        pa_ds = dl_tbl.to_pyarrow_dataset(**pyarrow_options)
        return scan_pyarrow_dataset(pa_ds)

    if pyarrow_options is not None:
        msg = "To make use of pyarrow_options, set use_pyarrow to True"
        raise ValueError(msg)

    from deltalake.exceptions import DeltaProtocolError
    from deltalake.table import (
        MAX_SUPPORTED_READER_VERSION,
        NOT_SUPPORTED_READER_VERSION,
        SUPPORTED_READER_FEATURES,
    )

    table_protocol = dl_tbl.protocol()
    if (
        table_protocol.min_reader_version > MAX_SUPPORTED_READER_VERSION
        or table_protocol.min_reader_version == NOT_SUPPORTED_READER_VERSION
    ):
        msg = (
            f"The table's minimum reader version is {table_protocol.min_reader_version} "
            f"but polars delta scanner only supports version 1 or {MAX_SUPPORTED_READER_VERSION} with these reader features: {SUPPORTED_READER_FEATURES}"
        )
        raise DeltaProtocolError(msg)
    if (
        table_protocol.min_reader_version >= 3
        and table_protocol.reader_features is not None
    ):
        missing_features = {*table_protocol.reader_features}.difference(
            SUPPORTED_READER_FEATURES
        )
        if len(missing_features) > 0:
            msg = f"The table has set these reader features: {missing_features} but these are not yet supported by the polars delta scanner."
            raise DeltaProtocolError(msg)

    delta_schema = dl_tbl.schema()
    polars_schema = Schema(delta_schema)
    partition_columns = dl_tbl.metadata().partition_columns

    def _split_schema(
        schema: Schema, partition_columns: list[str]
    ) -> tuple[Schema, Schema]:
        if len(partition_columns) == 0:
            return schema, Schema([])
        main_schema = []
        hive_schema = []

        for name, dtype in schema.items():
            if name in partition_columns:
                hive_schema.append((name, dtype))
            else:
                main_schema.append((name, dtype))

        return Schema(main_schema), Schema(hive_schema)

    # Required because main_schema cannot contain hive columns currently
    main_schema, hive_schema = _split_schema(polars_schema, partition_columns)

    file_uris = dl_tbl.file_uris()

    # LakeFS has an S3 compatible API, for reading therefore it's safe to do this.
    # Deltalake internally has an integration for writing commits
    if dl_tbl.table_uri.startswith("lakefs://"):
        file_uris = [file_uri.replace("lakefs://", "s3://") for file_uri in file_uris]

    return scan_parquet(
        file_uris,
        schema=main_schema,
        hive_schema=hive_schema if len(partition_columns) > 0 else None,
        cast_options=ScanCastOptions._default_iceberg(),
        missing_columns="insert",
        extra_columns="ignore",
        hive_partitioning=len(partition_columns) > 0,
        storage_options=storage_options,
        credential_provider=credential_provider_builder,  # type: ignore[arg-type]
        rechunk=rechunk or False,
    )


def _resolve_delta_lake_uri(table_uri: str | Path, *, strict: bool = True) -> str:
    resolved_uri = str(
        Path(table_uri).expanduser().resolve(strict)
        if _get_path_scheme(table_uri) is None
        else table_uri
    )

    return resolved_uri


def _get_delta_lake_table(
    table_path: str | Path | DeltaTable,
    version: int | str | datetime | None = None,
    storage_options: dict[str, Any] | None = None,
    delta_table_options: dict[str, Any] | None = None,
) -> deltalake.DeltaTable:
    """
    Initialize a Delta lake table for use in read and scan operations.

    Notes
    -----
    Make sure to install deltalake>=0.8.0. Read the documentation
    `here <https://delta-io.github.io/delta-rs/usage/installation/>`_.
    """
    _check_if_delta_available()

    if isinstance(table_path, deltalake.DeltaTable):
        if any(
            [
                version is not None,
                storage_options is not None,
                delta_table_options is not None,
            ]
        ):
            warnings.warn(
                """When supplying a DeltaTable directly, `version`, `storage_options`, and `delta_table_options` are ignored.
                To silence this warning, don't supply those parameters.""",
                RuntimeWarning,
                stacklevel=1,
            )
        return table_path
    if delta_table_options is None:
        delta_table_options = {}
    resolved_uri = _resolve_delta_lake_uri(table_path)
    if not isinstance(version, (str, datetime)):
        dl_tbl = deltalake.DeltaTable(
            resolved_uri,
            version=version,
            storage_options=storage_options,
            **delta_table_options,
        )
    else:
        dl_tbl = deltalake.DeltaTable(
            table_path,
            storage_options=storage_options,
            **delta_table_options,
        )
        dl_tbl.load_as_version(version)

    return dl_tbl


def _check_if_delta_available() -> None:
    if not _DELTALAKE_AVAILABLE:
        msg = "deltalake is not installed\n\nPlease run: pip install deltalake"
        raise ModuleNotFoundError(msg)


def _check_for_unsupported_types(dtypes: list[DataType]) -> None:
    schema_dtypes = unpack_dtypes(*dtypes)
    unsupported_types = {Time, Null}
    # Note that this overlap check does NOT work correctly for Categorical, so
    # if Categorical is added back to unsupported_types a different check will
    # need to be used.

    if overlap := schema_dtypes & unsupported_types:
        msg = f"dataframe contains unsupported data types: {overlap!r}"
        raise TypeError(msg)
