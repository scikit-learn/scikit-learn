from __future__ import annotations

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import partial
from time import perf_counter
from typing import TYPE_CHECKING, Any, Literal

import polars._reexport as pl
from polars._utils.logging import eprint, verbose, verbose_print_sensitive
from polars.exceptions import ComputeError
from polars.io.iceberg._utils import (
    IcebergStatisticsLoader,
    IdentityTransformedPartitionValuesBuilder,
    _scan_pyarrow_dataset_impl,
    try_convert_pyarrow_predicate,
)
from polars.io.scan_options.cast_options import ScanCastOptions

if TYPE_CHECKING:
    import pyarrow as pa
    import pyiceberg.schema
    from pyiceberg.table import Table

    from polars.lazyframe.frame import LazyFrame


class IcebergDataset:
    """Dataset interface for PyIceberg."""

    def __init__(
        self,
        source: str | Table,
        *,
        snapshot_id: int | None = None,
        iceberg_storage_properties: dict[str, Any] | None = None,
        reader_override: Literal["native", "pyiceberg"] | None = None,
        use_metadata_statistics: bool = True,
        fast_deletion_count: bool = True,
        use_pyiceberg_filter: bool = True,
    ) -> None:
        self._metadata_path = None
        self._table = None
        self._snapshot_id = snapshot_id
        self._iceberg_storage_properties = iceberg_storage_properties
        self._reader_override: Literal["native", "pyiceberg"] | None = reader_override
        self._use_metadata_statistics = use_metadata_statistics
        self._fast_deletion_count = fast_deletion_count
        self._use_pyiceberg_filter = use_pyiceberg_filter

        # Accept either a path or a table object. The one we don't have is
        # lazily initialized when needed.

        if isinstance(source, str):
            self._metadata_path = source
        else:
            self._table = source

    #
    # PythonDatasetProvider interface functions
    #

    def schema(self) -> pa.schema:
        """Fetch the schema of the table."""
        return self.arrow_schema()

    def arrow_schema(self) -> pa.schema:
        """Fetch the arrow schema of the table."""
        from pyiceberg.io.pyarrow import schema_to_pyarrow

        return schema_to_pyarrow(self.table().schema())

    def to_dataset_scan(
        self,
        *,
        existing_resolved_version_key: str | None = None,
        limit: int | None = None,
        projection: list[str] | None = None,
        filter_columns: list[str] | None = None,
        pyarrow_predicate: str | None = None,
    ) -> tuple[LazyFrame, str] | None:
        """Construct a LazyFrame scan."""
        if (
            scan_data := self._to_dataset_scan_impl(
                existing_resolved_version_key=existing_resolved_version_key,
                limit=limit,
                projection=projection,
                filter_columns=filter_columns,
                pyarrow_predicate=pyarrow_predicate,
            )
        ) is None:
            return None

        return scan_data.to_lazyframe(), scan_data.snapshot_id_key

    def _to_dataset_scan_impl(
        self,
        *,
        existing_resolved_version_key: str | None = None,
        limit: int | None = None,
        projection: list[str] | None = None,
        filter_columns: list[str] | None = None,
        pyarrow_predicate: str | None = None,
    ) -> _NativeIcebergScanData | _PyIcebergScanData | None:
        from pyiceberg.io.pyarrow import schema_to_pyarrow

        import polars._utils.logging

        verbose = polars._utils.logging.verbose()

        iceberg_table_filter = None

        if (
            pyarrow_predicate is not None
            and self._use_metadata_statistics
            and self._use_pyiceberg_filter
        ):
            iceberg_table_filter = try_convert_pyarrow_predicate(pyarrow_predicate)

        if verbose:
            pyarrow_predicate_display = (
                "Some(<redacted>)" if pyarrow_predicate is not None else "None"
            )
            iceberg_table_filter_display = (
                "Some(<redacted>)" if iceberg_table_filter is not None else "None"
            )

            eprint(
                "IcebergDataset: to_dataset_scan(): "
                f"snapshot ID: {self._snapshot_id}, "
                f"limit: {limit}, "
                f"projection: {projection}, "
                f"filter_columns: {filter_columns}, "
                f"pyarrow_predicate: {pyarrow_predicate_display}, "
                f"iceberg_table_filter: {iceberg_table_filter_display}, "
                f"self._use_metadata_statistics: {self._use_metadata_statistics}"
            )

        verbose_print_sensitive(
            lambda: f"IcebergDataset: to_dataset_scan(): {pyarrow_predicate = }, {iceberg_table_filter = }"
        )

        tbl = self.table()

        if verbose:
            eprint(
                "IcebergDataset: to_dataset_scan(): "
                f"tbl.metadata.current_snapshot_id: {tbl.metadata.current_snapshot_id}"
            )

        snapshot_id = self._snapshot_id
        schema_id = None

        if snapshot_id is not None:
            snapshot = tbl.snapshot_by_id(snapshot_id)

            if snapshot is None:
                msg = f"iceberg snapshot ID not found: {snapshot_id}"
                raise ValueError(msg)

            schema_id = snapshot.schema_id

            if schema_id is None:
                msg = (
                    f"IcebergDataset: requested snapshot {snapshot_id} "
                    "did not contain a schema ID"
                )
                raise ValueError(msg)

            iceberg_schema = tbl.schemas()[schema_id]
            snapshot_id_key = f"{snapshot.snapshot_id}"
        else:
            iceberg_schema = tbl.schema()
            schema_id = tbl.metadata.current_schema_id

            snapshot_id_key = (
                f"{v.snapshot_id}" if (v := tbl.current_snapshot()) is not None else ""
            )

        if (
            existing_resolved_version_key is not None
            and existing_resolved_version_key == snapshot_id_key
        ):
            if verbose:
                eprint(
                    "IcebergDataset: to_dataset_scan(): early return "
                    f"({snapshot_id_key = })"
                )

            return None

        # Take from parameter first then envvar
        reader_override = self._reader_override or os.getenv(
            "POLARS_ICEBERG_READER_OVERRIDE"
        )

        if reader_override and reader_override not in ["native", "pyiceberg"]:
            msg = (
                "iceberg: unknown value for reader_override: "
                f"'{reader_override}', expected one of ('native', 'pyiceberg')"
            )
            raise ValueError(msg)

        fallback_reason = (
            "forced reader_override='pyiceberg'"
            if reader_override == "pyiceberg"
            else f"unsupported table format version: {tbl.format_version}"
            if not tbl.format_version <= 2
            else None
        )

        selected_fields = ("*",) if projection is None else tuple(projection)

        projected_iceberg_schema = (
            iceberg_schema
            if selected_fields == ("*",)
            else iceberg_schema.select(*selected_fields)
        )

        sources = []
        missing_field_defaults = IdentityTransformedPartitionValuesBuilder(
            tbl,
            projected_iceberg_schema,
        )
        statistics_loader: IcebergStatisticsLoader | None = (
            IcebergStatisticsLoader(tbl, iceberg_schema.select(*filter_columns))
            if self._use_metadata_statistics and filter_columns is not None
            else None
        )
        deletion_files: dict[int, list[str]] = {}
        total_physical_rows: int = 0
        total_deleted_rows: int = 0

        if reader_override != "pyiceberg" and not fallback_reason:
            from pyiceberg.manifest import DataFileContent, FileFormat

            if verbose:
                eprint("IcebergDataset: to_dataset_scan(): begin path expansion")

            start_time = perf_counter()

            scan = tbl.scan(
                snapshot_id=snapshot_id,
                limit=limit,
                selected_fields=selected_fields,
            )

            if iceberg_table_filter is not None:
                scan = scan.filter(iceberg_table_filter)

            total_deletion_files = 0

            for i, file_info in enumerate(scan.plan_files()):
                if file_info.file.file_format != FileFormat.PARQUET:
                    fallback_reason = (
                        f"non-parquet format: {file_info.file.file_format}"
                    )
                    break

                if file_info.delete_files:
                    deletion_files[i] = []

                    for deletion_file in file_info.delete_files:
                        if deletion_file.content != DataFileContent.POSITION_DELETES:
                            fallback_reason = (
                                "unsupported deletion file type: "
                                f"{deletion_file.content}"
                            )
                            break

                        if deletion_file.file_format != FileFormat.PARQUET:
                            fallback_reason = (
                                "unsupported deletion file format: "
                                f"{deletion_file.file_format}"
                            )
                            break

                        deletion_files[i].append(deletion_file.file_path)
                        total_deletion_files += 1
                        total_deleted_rows += deletion_file.record_count

                if fallback_reason:
                    break

                missing_field_defaults.push_partition_values(
                    current_index=i,
                    partition_spec_id=file_info.file.spec_id,
                    partition_values=file_info.file.partition,
                )

                if statistics_loader is not None:
                    statistics_loader.push_file_statistics(file_info.file)

                total_physical_rows += file_info.file.record_count

                sources.append(file_info.file.file_path)

            if verbose:
                elapsed = perf_counter() - start_time
                eprint(
                    "IcebergDataset: to_dataset_scan(): "
                    f"finish path expansion ({elapsed:.3f}s)"
                )

        if not fallback_reason:
            if verbose:
                s = "" if len(sources) == 1 else "s"
                s2 = "" if total_deletion_files == 1 else "s"

                eprint(
                    "IcebergDataset: to_dataset_scan(): "
                    f"native scan_parquet(): "
                    f"{len(sources)} source{s}, "
                    f"snapshot ID: {snapshot_id}, "
                    f"schema ID: {schema_id}, "
                    f"{total_deletion_files} deletion file{s2}"
                )

            # The arrow schema returned by `schema_to_pyarrow` will contain
            # 'PARQUET:field_id'
            column_mapping = schema_to_pyarrow(iceberg_schema)

            identity_transformed_values = missing_field_defaults.finish()

            min_max_statistics = (
                statistics_loader.finish(len(sources), identity_transformed_values)
                if statistics_loader is not None
                else None
            )

            storage_options = (
                _convert_iceberg_to_object_store_storage_options(
                    self._iceberg_storage_properties
                )
                if self._iceberg_storage_properties is not None
                else None
            )

            return _NativeIcebergScanData(
                sources=sources,
                projected_iceberg_schema=projected_iceberg_schema,
                column_mapping=column_mapping,
                default_values=identity_transformed_values,
                deletion_files=deletion_files,
                min_max_statistics=min_max_statistics,
                statistics_loader=statistics_loader,
                storage_options=storage_options,
                row_count=(
                    (total_physical_rows, total_deleted_rows)
                    if (
                        self._use_metadata_statistics
                        and (self._fast_deletion_count or total_deleted_rows == 0)
                    )
                    else None
                ),
                _snapshot_id_key=snapshot_id_key,
            )

        elif reader_override == "native":
            msg = f"iceberg reader_override='native' failed: {fallback_reason}"
            raise ComputeError(msg)

        if verbose:
            eprint(
                "IcebergDataset: to_dataset_scan(): "
                f"fallback to python[pyiceberg] scan: {fallback_reason}"
            )

        func = partial(
            _scan_pyarrow_dataset_impl,
            tbl,
            snapshot_id=snapshot_id,
            n_rows=limit,
            with_columns=projection,
            iceberg_table_filter=iceberg_table_filter,
        )

        arrow_schema = schema_to_pyarrow(tbl.schema())

        lf = pl.LazyFrame._scan_python_function(
            arrow_schema,
            func,
            pyarrow=True,
            is_pure=True,
        )

        return _PyIcebergScanData(lf=lf, _snapshot_id_key=snapshot_id_key)

    #
    # Accessors
    #

    def metadata_path(self) -> str:
        """Fetch the metadata path."""
        if self._metadata_path is None:
            if self._table is None:
                msg = "impl error: both metadata_path and table are None"
                raise ValueError(msg)

            self._metadata_path = self.table().metadata_location

        return self._metadata_path

    def table(self) -> Table:
        """Fetch the PyIceberg Table object."""
        if self._table is None:
            if self._metadata_path is None:
                msg = "impl error: both metadata_path and table are None"
                raise ValueError(msg)

            if verbose():
                eprint(f"IcebergDataset: construct table from {self._metadata_path = }")

            from pyiceberg.table import StaticTable

            self._table = StaticTable.from_metadata(
                metadata_location=self._metadata_path,
                properties=self._iceberg_storage_properties or {},
            )

        return self._table

    #
    # Serialization functions
    #
    # We don't serialize the iceberg table object - the remote machine should
    # use their own permissions to reconstruct the table object from the path.
    #

    def __getstate__(self) -> dict[str, Any]:
        state = {
            "metadata_path": self.metadata_path(),
            "snapshot_id": self._snapshot_id,
            "iceberg_storage_properties": self._iceberg_storage_properties,
            "reader_override": self._reader_override,
            "use_metadata_statistics": self._use_metadata_statistics,
            "fast_deletion_count": self._fast_deletion_count,
            "use_pyiceberg_filter": self._use_pyiceberg_filter,
        }

        if verbose():
            path_repr = state["metadata_path"]
            snapshot_id = f"'{v}'" if (v := state["snapshot_id"]) is not None else None
            keys_repr = _redact_dict_values(state["iceberg_storage_properties"])
            reader_override = state["reader_override"]
            use_metadata_statistics = state["use_metadata_statistics"]
            fast_deletion_count = state["fast_deletion_count"]
            use_pyiceberg_filter = state["use_pyiceberg_filter"]

            eprint(
                "IcebergDataset: getstate(): "
                f"path: '{path_repr}', "
                f"snapshot_id: {snapshot_id}, "
                f"iceberg_storage_properties: {keys_repr}, "
                f"reader_override: {reader_override}, "
                f"use_metadata_statistics: {use_metadata_statistics}, "
                f"fast_deletion_count: {fast_deletion_count}, "
                f"use_pyiceberg_filter: {use_pyiceberg_filter}"
            )

        return state

    def __setstate__(self, state: dict[str, Any]) -> None:
        if verbose():
            path_repr = state["metadata_path"]
            snapshot_id = f"'{v}'" if (v := state["snapshot_id"]) is not None else None
            keys_repr = _redact_dict_values(state["iceberg_storage_properties"])
            reader_override = state["reader_override"]
            use_metadata_statistics = state["use_metadata_statistics"]
            fast_deletion_count = state["fast_deletion_count"]
            use_pyiceberg_filter = state["use_pyiceberg_filter"]

            eprint(
                "IcebergDataset: getstate(): "
                f"path: '{path_repr}', "
                f"snapshot_id: '{snapshot_id}', "
                f"iceberg_storage_properties: {keys_repr}, "
                f"reader_override: {reader_override}, "
                f"use_metadata_statistics: {use_metadata_statistics}, "
                f"fast_deletion_count: {fast_deletion_count}, "
                f"use_pyiceberg_filter: {use_pyiceberg_filter}"
            )

        IcebergDataset.__init__(
            self,
            state["metadata_path"],
            snapshot_id=state["snapshot_id"],
            iceberg_storage_properties=state["iceberg_storage_properties"],
            reader_override=state["reader_override"],
            use_metadata_statistics=state["use_metadata_statistics"],
            fast_deletion_count=state["fast_deletion_count"],
            use_pyiceberg_filter=state["use_pyiceberg_filter"],
        )


class _ResolvedScanDataBase(ABC):
    @abstractmethod
    def to_lazyframe(self) -> pl.LazyFrame: ...

    @property
    @abstractmethod
    def snapshot_id_key(self) -> str: ...


@dataclass
class _NativeIcebergScanData(_ResolvedScanDataBase):
    """Resolved parameters for a native Iceberg scan."""

    sources: list[str]
    projected_iceberg_schema: pyiceberg.schema.Schema
    column_mapping: pa.Schema
    default_values: dict[int, pl.Series | str]
    deletion_files: dict[int, list[str]]
    min_max_statistics: pl.DataFrame | None
    # This is here for test purposes, as the `min_max_statistics` on this
    # dataclass contain coalesced values from `default_values`, a test may
    # access the statistics loader directly to inspect the values before
    # coalescing.
    statistics_loader: IcebergStatisticsLoader | None
    storage_options: dict[str, str] | None
    # (physical, deleted)
    row_count: tuple[int, int] | None
    _snapshot_id_key: str

    def to_lazyframe(self) -> pl.LazyFrame:
        from polars.io.parquet.functions import scan_parquet

        return scan_parquet(
            self.sources,
            cast_options=ScanCastOptions._default_iceberg(),
            missing_columns="insert",
            extra_columns="ignore",
            storage_options=self.storage_options,
            _column_mapping=("iceberg-column-mapping", self.column_mapping),
            _default_values=("iceberg", self.default_values),
            _deletion_files=("iceberg-position-delete", self.deletion_files),
            _table_statistics=self.min_max_statistics,
            _row_count=self.row_count,
        )

    @property
    def snapshot_id_key(self) -> str:
        return self._snapshot_id_key


@dataclass
class _PyIcebergScanData(_ResolvedScanDataBase):
    """Resolved parameters for reading via PyIceberg."""

    # We're not interested in inspecting anything for the pyiceberg scan, so
    # this class is just a wrapper.
    lf: pl.LazyFrame
    _snapshot_id_key: str

    def to_lazyframe(self) -> pl.LazyFrame:
        return self.lf

    @property
    def snapshot_id_key(self) -> str:
        return self._snapshot_id_key


def _redact_dict_values(obj: Any) -> Any:
    return (
        {k: "REDACTED" for k in obj.keys()}  # noqa: SIM118
        if isinstance(obj, dict)
        else f"<{type(obj).__name__} object>"
        if obj is not None
        else "None"
    )


def _convert_iceberg_to_object_store_storage_options(
    iceberg_storage_properties: dict[str, str],
) -> dict[str, str]:
    storage_options = {}

    for k, v in iceberg_storage_properties.items():
        if (
            translated_key := ICEBERG_TO_OBJECT_STORE_CONFIG_KEY_MAP.get(k)
        ) is not None:
            storage_options[translated_key] = v
        elif "." not in k:
            # Pass-through non-Iceberg config keys, as they may be native config
            # keys. We identify Iceberg keys by checking for a dot - from
            # observation nearly all Iceberg config keys contain dots, whereas
            # native config keys do not contain them.
            storage_options[k] = v

        # Otherwise, unknown keys are ignored / not passed. This is to avoid
        # interfering with credential provider auto-init, which bails on
        # unknown keys.

    return storage_options


# https://py.iceberg.apache.org/configuration/#fileio
# This does not contain all keys - some have no object-store equivalent.
ICEBERG_TO_OBJECT_STORE_CONFIG_KEY_MAP: dict[str, str] = {
    # S3
    "s3.endpoint": "aws_endpoint_url",
    "s3.access-key-id": "aws_access_key_id",
    "s3.secret-access-key": "aws_secret_access_key",
    "s3.session-token": "aws_session_token",
    "s3.region": "aws_region",
    "s3.proxy-uri": "proxy_url",
    "s3.connect-timeout": "connect_timeout",
    "s3.request-timeout": "timeout",
    "s3.force-virtual-addressing": "aws_virtual_hosted_style_request",
    # Azure
    "adls.account-name": "azure_storage_account_name",
    "adls.account-key": "azure_storage_account_key",
    "adls.sas-token": "azure_storage_sas_key",
    "adls.tenant-id": "azure_storage_tenant_id",
    "adls.client-id": "azure_storage_client_id",
    "adls.client-secret": "azure_storage_client_secret",
    "adls.account-host": "azure_storage_authority_host",
    "adls.token": "azure_storage_token",
    # Google storage
    "gcs.oauth2.token": "bearer_token",
    # HuggingFace
    "hf.token": "token",
}
