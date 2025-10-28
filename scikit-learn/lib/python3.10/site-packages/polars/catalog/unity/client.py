from __future__ import annotations

import contextlib
import importlib
import os
import sys
from typing import TYPE_CHECKING, Any, Literal

from polars._utils.unstable import issue_unstable_warning
from polars._utils.wrap import wrap_ldf
from polars.catalog.unity.models import (
    CatalogInfo,
    ColumnInfo,
    NamespaceInfo,
    TableInfo,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from datetime import datetime

    import deltalake

    from polars._typing import SchemaDict
    from polars.catalog.unity.models import DataSourceFormat, TableType
    from polars.dataframe.frame import DataFrame
    from polars.io.cloud import (
        CredentialProviderFunction,
        CredentialProviderFunctionReturn,
    )
    from polars.io.cloud.credential_provider._builder import CredentialProviderBuilder
    from polars.lazyframe import LazyFrame

with contextlib.suppress(ImportError):
    from polars._plr import PyCatalogClient

    PyCatalogClient.init_classes(
        catalog_info_cls=CatalogInfo,
        namespace_info_cls=NamespaceInfo,
        table_info_cls=TableInfo,
        column_info_cls=ColumnInfo,
    )


class Catalog:
    """
    Unity catalog client.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    """

    def __init__(
        self,
        workspace_url: str,
        *,
        bearer_token: str | None = "auto",
        require_https: bool = True,
    ) -> None:
        """
        Initialize a catalog client.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        workspace_url
            URL of the workspace, or alternatively the URL of the Unity catalog
            API endpoint.
        bearer_token
            Bearer token to authenticate with. This can also be set to:

            * "auto": Automatically retrieve bearer tokens from the environment.
            * "databricks-sdk": Use the Databricks SDK to retrieve and use the
              bearer token from the environment.
        require_https
            Require the `workspace_url` to use HTTPS.
        """
        issue_unstable_warning("`Catalog` functionality is considered unstable.")

        if require_https and not workspace_url.startswith("https://"):
            msg = (
                f"a non-HTTPS workspace_url was given ({workspace_url}). To "
                "allow non-HTTPS URLs, pass require_https=False."
            )
            raise ValueError(msg)

        if bearer_token == "databricks-sdk" or (
            bearer_token == "auto"
            # For security, in "auto" mode, only retrieve/use the token if:
            # * We are running inside a Databricks environment
            # * The `workspace_url` is pointing to Databricks and uses HTTPS
            and "DATABRICKS_RUNTIME_VERSION" in os.environ
            and workspace_url.startswith("https://")
            and (
                workspace_url.removeprefix("https://")
                .split("/", 1)[0]
                .endswith(".cloud.databricks.com")
            )
        ):
            bearer_token = self._get_databricks_token()

        if bearer_token == "auto":
            bearer_token = None

        self._client = PyCatalogClient.new(workspace_url, bearer_token)

    def list_catalogs(self) -> list[CatalogInfo]:
        """
        List the available catalogs.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
        """
        return self._client.list_catalogs()

    def list_namespaces(self, catalog_name: str) -> list[NamespaceInfo]:
        """
        List the available namespaces (unity schema) under the specified catalog.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        """
        return self._client.list_namespaces(catalog_name)

    def list_tables(self, catalog_name: str, namespace: str) -> list[TableInfo]:
        """
        List the available tables under the specified schema.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        """
        return self._client.list_tables(catalog_name, namespace)

    def get_table_info(
        self, catalog_name: str, namespace: str, table_name: str
    ) -> TableInfo:
        """
        Retrieve the metadata of the specified table.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        table_name
            Name of the table.
        """
        return self._client.get_table_info(catalog_name, namespace, table_name)

    def _get_table_credentials(
        self, table_id: str, *, write: bool
    ) -> tuple[dict[str, str] | None, dict[str, str], int]:
        return self._client.get_table_credentials(table_id=table_id, write=write)

    def scan_table(
        self,
        catalog_name: str,
        namespace: str,
        table_name: str,
        *,
        delta_table_version: int | str | datetime | None = None,
        delta_table_options: dict[str, Any] | None = None,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
    ) -> LazyFrame:
        """
        Retrieve the metadata of the specified table.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        table_name
            Name of the table.
        delta_table_version
            Version of the table to scan (Deltalake only).
        delta_table_options
            Additional keyword arguments while reading a Deltalake table.
        storage_options
            Options that indicate how to connect to a cloud provider.

            The cloud providers currently supported are AWS, GCP, and Azure.
            See supported keys here:

            * `aws <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html>`_
            * `gcp <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html>`_
            * `azure <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html>`_
            * Hugging Face (`hf://`): Accepts an API key under the `token` parameter: \
            `{'token': '...'}`, or by setting the `HF_TOKEN` environment variable.

            If `storage_options` is not provided, Polars will try to infer the
            information from environment variables.
        credential_provider
            Provide a function that can be called to provide cloud storage
            credentials. The function is expected to return a dictionary of
            credential keys along with an optional credential expiry time.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        retries
            Number of retries if accessing a cloud instance fails.

        """
        table_info = self.get_table_info(catalog_name, namespace, table_name)
        storage_location, data_source_format = _extract_location_and_data_format(
            table_info, "scan table"
        )

        credential_provider, storage_options = self._init_credentials(  # type: ignore[assignment]
            credential_provider,
            storage_options,
            table_info,
            write=False,
            caller_name="Catalog.scan_table",
        )

        if data_source_format in ["DELTA", "DELTASHARING"]:
            from polars.io.delta import scan_delta

            return scan_delta(
                storage_location,
                version=delta_table_version,
                delta_table_options=delta_table_options,
                storage_options=storage_options,
                credential_provider=credential_provider,
            )

        if delta_table_version is not None:
            msg = (
                "cannot apply delta_table_version for table of type "
                f"{data_source_format}"
            )
            raise ValueError(msg)

        if delta_table_options is not None:
            msg = (
                "cannot apply delta_table_options for table of type "
                f"{data_source_format}"
            )
            raise ValueError(msg)

        if storage_options:
            storage_options = list(storage_options.items())  # type: ignore[assignment]
        else:
            # Handle empty dict input
            storage_options = None

        return wrap_ldf(
            self._client.scan_table(
                catalog_name,
                namespace,
                table_name,
                credential_provider=credential_provider,
                cloud_options=storage_options,
                retries=retries,
            )
        )

    def write_table(
        self,
        df: DataFrame,
        catalog_name: str,
        namespace: str,
        table_name: str,
        *,
        delta_mode: Literal[
            "error", "append", "overwrite", "ignore", "merge"
        ] = "error",
        delta_write_options: dict[str, Any] | None = None,
        delta_merge_options: dict[str, Any] | None = None,
        storage_options: dict[str, str] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
    ) -> None | deltalake.table.TableMerger:
        """
        Write a DataFrame to a catalog table.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        df
            DataFrame to write.
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        table_name
            Name of the table.
        delta_mode : {'error', 'append', 'overwrite', 'ignore', 'merge'}
            (For delta tables) How to handle existing data.

            - If 'error', throw an error if the table already exists (default).
            - If 'append', will add new data.
            - If 'overwrite', will replace table with new data.
            - If 'ignore', will not write anything if table already exists.
            - If 'merge', return a `TableMerger` object to merge data from the DataFrame
              with the existing data.
        delta_write_options
            (For delta tables) Additional keyword arguments while writing a
            Delta lake Table.
            See a list of supported write options `here <https://delta-io.github.io/delta-rs/api/delta_writer/#deltalake.write_deltalake>`__.
        delta_merge_options
            (For delta tables) Keyword arguments which are required to `MERGE` a
            Delta lake Table.
            See a list of supported merge options `here <https://delta-io.github.io/delta-rs/api/delta_table/#deltalake.DeltaTable.merge>`__.
        storage_options
            Options that indicate how to connect to a cloud provider.

            The cloud providers currently supported are AWS, GCP, and Azure.
            See supported keys here:

            * `aws <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html>`_
            * `gcp <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html>`_
            * `azure <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html>`_
            * Hugging Face (`hf://`): Accepts an API key under the `token` parameter: \
            `{'token': '...'}`, or by setting the `HF_TOKEN` environment variable.

            If `storage_options` is not provided, Polars will try to infer the
            information from environment variables.
        credential_provider
            Provide a function that can be called to provide cloud storage
            credentials. The function is expected to return a dictionary of
            credential keys along with an optional credential expiry time.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        """
        table_info = self.get_table_info(catalog_name, namespace, table_name)
        storage_location, data_source_format = _extract_location_and_data_format(
            table_info, "scan table"
        )

        credential_provider, storage_options = self._init_credentials(  # type: ignore[assignment]
            credential_provider,
            storage_options,
            table_info,
            write=True,
            caller_name="Catalog.write_table",
        )

        if data_source_format in ["DELTA", "DELTASHARING"]:
            return df.write_delta(  # type: ignore[misc]
                storage_location,
                storage_options=storage_options,
                credential_provider=credential_provider,
                mode=delta_mode,
                delta_write_options=delta_write_options,
                delta_merge_options=delta_merge_options,
            )  # type: ignore[call-overload]

        else:
            msg = (
                "write_table: table format of "
                f"{catalog_name}.{namespace}.{table_name} "
                f"({data_source_format}) is unsupported."
            )
            raise NotImplementedError(msg)

    def create_catalog(
        self,
        catalog_name: str,
        *,
        comment: str | None = None,
        storage_root: str | None = None,
    ) -> CatalogInfo:
        """
        Create a catalog.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        comment
            Leaves a comment about the catalog.
        storage_root
            Base location at which to store the catalog.
        """
        return self._client.create_catalog(
            catalog_name=catalog_name, comment=comment, storage_root=storage_root
        )

    def delete_catalog(
        self,
        catalog_name: str,
        *,
        force: bool = False,
    ) -> None:
        """
        Delete a catalog.

        Note that depending on the table type and catalog server, this may not
        delete the actual data files from storage. For more details, please
        consult the documentation of the catalog provider you are using.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        force
            Forcibly delete the catalog even if it is not empty.
        """
        self._client.delete_catalog(catalog_name=catalog_name, force=force)

    def create_namespace(
        self,
        catalog_name: str,
        namespace: str,
        *,
        comment: str | None = None,
        storage_root: str | None = None,
    ) -> NamespaceInfo:
        """
        Create a namespace (unity schema) in the catalog.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        comment
            Leaves a comment about the table.
        storage_root
            Base location at which to store the namespace.
        """
        return self._client.create_namespace(
            catalog_name=catalog_name,
            namespace=namespace,
            comment=comment,
            storage_root=storage_root,
        )

    def delete_namespace(
        self,
        catalog_name: str,
        namespace: str,
        *,
        force: bool = False,
    ) -> None:
        """
        Delete a namespace (unity schema) in the catalog.

        Note that depending on the table type and catalog server, this may not
        delete the actual data files from storage. For more details, please
        consult the documentation of the catalog provider you are using.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        force
            Forcibly delete the namespace even if it is not empty.
        """
        self._client.delete_namespace(
            catalog_name=catalog_name, namespace=namespace, force=force
        )

    def create_table(
        self,
        catalog_name: str,
        namespace: str,
        table_name: str,
        *,
        schema: SchemaDict | None,
        table_type: TableType,
        data_source_format: DataSourceFormat | None = None,
        comment: str | None = None,
        storage_root: str | None = None,
        properties: dict[str, str] | None = None,
    ) -> TableInfo:
        """
        Create a table in the catalog.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        table_name
            Name of the table.
        schema
            Schema of the table.
        table_type
            Type of the table
        data_source_format
            Storage format of the table.
        comment
            Leaves a comment about the table.
        storage_root
            Base location at which to store the table.
        properties
            Extra key-value metadata to store.
        """
        return self._client.create_table(
            catalog_name=catalog_name,
            namespace=namespace,
            table_name=table_name,
            schema=schema,
            table_type=table_type,
            data_source_format=data_source_format,
            comment=comment,
            storage_root=storage_root,
            properties=list((properties or {}).items()),
        )

    def delete_table(
        self,
        catalog_name: str,
        namespace: str,
        table_name: str,
    ) -> None:
        """
        Delete the table stored at this location.

        Note that depending on the table type and catalog server, this may not
        delete the actual data files from storage. For more details, please
        consult the documentation of the catalog provider you are using.

        If you would like to perform manual deletions, the storage location of
        the files can be found using `get_table_info`.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        catalog_name
            Name of the catalog.
        namespace
            Name of the namespace (unity schema).
        table_name
            Name of the table.
        """
        self._client.delete_table(
            catalog_name=catalog_name,
            namespace=namespace,
            table_name=table_name,
        )

    def _init_credentials(
        self,
        credential_provider: CredentialProviderFunction | Literal["auto"] | None,
        storage_options: dict[str, Any] | None,
        table_info: TableInfo,
        *,
        write: bool,
        caller_name: str,
    ) -> tuple[
        CredentialProviderBuilder | None,
        dict[str, Any] | None,
    ]:
        from polars.io.cloud.credential_provider._builder import (
            CredentialProviderBuilder,
        )

        if credential_provider != "auto":
            if credential_provider:
                return CredentialProviderBuilder.from_initialized_provider(
                    credential_provider
                ), storage_options
            else:
                return None, storage_options

        verbose = os.getenv("POLARS_VERBOSE") == "1"

        catalog_credential_provider = CatalogCredentialProvider(
            self, table_info.table_id, write=write
        )

        try:
            v = catalog_credential_provider._credentials_iter()
            storage_update_options = next(v)

            if storage_update_options:
                storage_options = {**(storage_options or {}), **storage_update_options}

            for _ in v:
                pass

        except Exception as e:
            if verbose:
                table_name = table_info.name
                table_id = table_info.table_id
                msg = (
                    f"error auto-initializing CatalogCredentialProvider: {e!r} "
                    f"{table_name = } ({table_id = }) ({write = })"
                )
                print(msg, file=sys.stderr)
        else:
            if verbose:
                table_name = table_info.name
                table_id = table_info.table_id
                msg = (
                    "auto-selected CatalogCredentialProvider for "
                    f"{table_name = } ({table_id = })"
                )
                print(msg, file=sys.stderr)

            return CredentialProviderBuilder.from_initialized_provider(
                catalog_credential_provider
            ), storage_options

        # This should generally not happen, but if using the temporary
        # credentials API fails for whatever reason, we fallback to our built-in
        # credential provider resolution.

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )

        return _init_credential_provider_builder(
            "auto", table_info.storage_location, storage_options, caller_name
        ), storage_options

    @classmethod
    def _get_databricks_token(cls) -> str:
        if importlib.util.find_spec("databricks.sdk") is None:
            msg = "could not get Databricks token: databricks-sdk is not installed"
            raise ImportError(msg)

        # We code like this to bypass linting
        m = importlib.import_module("databricks.sdk.core").__dict__

        return m["DefaultCredentials"]()(m["Config"]())()["Authorization"][7:]


class CatalogCredentialProvider:
    """Retrieves credentials from the Unity catalog temporary credentials API."""

    def __init__(self, catalog: Catalog, table_id: str, *, write: bool) -> None:
        self.catalog = catalog
        self.table_id = table_id
        self.write = write

    def __call__(self) -> CredentialProviderFunctionReturn:  # noqa: D102
        _, (creds, expiry) = self._credentials_iter()
        return creds, expiry

    def _credentials_iter(
        self,
    ) -> Generator[Any]:
        creds, storage_update_options, expiry = self.catalog._get_table_credentials(
            self.table_id, write=self.write
        )

        yield storage_update_options

        if not creds:
            table_id = self.table_id
            msg = (
                "did not receive credentials from temporary credentials API for "
                f"{table_id = }"
            )
            raise Exception(msg)  # noqa: TRY002

        yield creds, expiry


def _extract_location_and_data_format(
    table_info: TableInfo, operation: str
) -> tuple[str, DataSourceFormat]:
    if table_info.storage_location is None:
        msg = f"cannot {operation}: no storage_location found"
        raise ValueError(msg)

    if table_info.data_source_format is None:
        msg = f"cannot {operation}: no data_source_format found"
        raise ValueError(msg)

    return table_info.storage_location, table_info.data_source_format
