from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars._typing import (
        ColumnMapping,
        DefaultFieldValues,
        DeletionFiles,
        SchemaDict,
    )
    from polars.dataframe.frame import DataFrame
    from polars.io.cloud.credential_provider._builder import CredentialProviderBuilder
    from polars.io.scan_options.cast_options import ScanCastOptions

from dataclasses import dataclass


# TODO: Add `kw_only=True` after 3.9 support dropped
@dataclass
class ScanOptions:
    """
    Holds scan options that are generic over scan type.

    For internal use. Most of the options will parse into `UnifiedScanArgs`.
    """

    row_index: tuple[str, int] | None = None
    # (i64, usize)
    pre_slice: tuple[int, int] | None = None
    cast_options: ScanCastOptions | None = None
    extra_columns: Literal["ignore", "raise"] = "raise"
    missing_columns: Literal["insert", "raise"] = "raise"
    include_file_paths: str | None = None

    # For path expansion
    glob: bool = True
    hidden_file_prefix: Sequence[str] | None = None

    # Hive
    # Note: `None` means auto.
    hive_partitioning: bool | None = None
    hive_schema: SchemaDict | None = None
    try_parse_hive_dates: bool = True

    rechunk: bool = False
    cache: bool = True

    # Cloud
    storage_options: list[tuple[str, str]] | None = None
    credential_provider: CredentialProviderBuilder | None = None
    retries: int = 2

    column_mapping: ColumnMapping | None = None
    default_values: DefaultFieldValues | None = None
    deletion_files: DeletionFiles | None = None
    table_statistics: DataFrame | None = None
    # (physical, deleted)
    row_count: tuple[int, int] | None = None
