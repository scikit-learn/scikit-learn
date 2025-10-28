from __future__ import annotations

import contextlib
from collections.abc import Sequence
from typing import TYPE_CHECKING

from polars._utils.various import (
    _process_null_values,
    normalize_filepath,
)
from polars._utils.wrap import wrap_df
from polars.datatypes import N_INFER_DEFAULT, parse_into_dtype
from polars.io._utils import parse_columns_arg, parse_row_index_args
from polars.io.csv._utils import _update_columns

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyBatchedCsv

if TYPE_CHECKING:
    from pathlib import Path

    from polars import DataFrame
    from polars._typing import CsvEncoding, PolarsDataType, SchemaDict


class BatchedCsvReader:
    """Read a CSV file in batches."""

    def __init__(
        self,
        source: str | Path,
        *,
        has_header: bool = True,
        columns: Sequence[int] | Sequence[str] | None = None,
        separator: str = ",",
        comment_prefix: str | None = None,
        quote_char: str | None = '"',
        skip_rows: int = 0,
        skip_lines: int = 0,
        schema_overrides: SchemaDict | Sequence[PolarsDataType] | None = None,
        null_values: str | Sequence[str] | dict[str, str] | None = None,
        missing_utf8_is_empty_string: bool = False,
        ignore_errors: bool = False,
        try_parse_dates: bool = False,
        n_threads: int | None = None,
        infer_schema_length: int | None = N_INFER_DEFAULT,
        batch_size: int = 50_000,
        n_rows: int | None = None,
        encoding: CsvEncoding = "utf8",
        low_memory: bool = False,
        rechunk: bool = True,
        skip_rows_after_header: int = 0,
        row_index_name: str | None = None,
        row_index_offset: int = 0,
        eol_char: str = "\n",
        new_columns: Sequence[str] | None = None,
        raise_if_empty: bool = True,
        truncate_ragged_lines: bool = False,
        decimal_comma: bool = False,
    ) -> None:
        path = normalize_filepath(source, check_not_directory=False)

        dtype_list: Sequence[tuple[str, PolarsDataType]] | None = None
        dtype_slice: Sequence[PolarsDataType] | None = None
        if schema_overrides is not None:
            if isinstance(schema_overrides, dict):
                dtype_list = []
                for k, v in schema_overrides.items():
                    dtype_list.append((k, parse_into_dtype(v)))
            elif isinstance(schema_overrides, Sequence):
                dtype_slice = schema_overrides
            else:
                msg = "`schema_overrides` arg should be list or dict"
                raise TypeError(msg)

        processed_null_values = _process_null_values(null_values)
        projection, columns = parse_columns_arg(columns)

        self._reader = PyBatchedCsv.new(
            infer_schema_length=infer_schema_length,
            chunk_size=batch_size,
            has_header=has_header,
            ignore_errors=ignore_errors,
            n_rows=n_rows,
            skip_rows=skip_rows,
            skip_lines=skip_lines,
            projection=projection,
            separator=separator,
            rechunk=rechunk,
            columns=columns,
            encoding=encoding,
            n_threads=n_threads,
            path=path,
            schema_overrides=dtype_list,
            overwrite_dtype_slice=dtype_slice,
            low_memory=low_memory,
            comment_prefix=comment_prefix,
            quote_char=quote_char,
            null_values=processed_null_values,
            missing_utf8_is_empty_string=missing_utf8_is_empty_string,
            try_parse_dates=try_parse_dates,
            skip_rows_after_header=skip_rows_after_header,
            row_index=parse_row_index_args(row_index_name, row_index_offset),
            eol_char=eol_char,
            raise_if_empty=raise_if_empty,
            truncate_ragged_lines=truncate_ragged_lines,
            decimal_comma=decimal_comma,
        )
        self.new_columns = new_columns

    def next_batches(self, n: int) -> list[DataFrame] | None:
        """
        Read `n` batches from the reader.

        These batches will be parallelized over the available threads.

        Parameters
        ----------
        n
            Number of chunks to fetch; ideally this is >= number of threads.

        Examples
        --------
        >>> reader = pl.read_csv_batched(
        ...     "./pdsh/tables_scale_100/lineitem.tbl",
        ...     separator="|",
        ...     try_parse_dates=True,
        ... )  # doctest: +SKIP
        >>> reader.next_batches(5)  # doctest: +SKIP

        Returns
        -------
        list of DataFrames
        """
        if (batches := self._reader.next_batches(n)) is not None:
            if self.new_columns:
                return [
                    _update_columns(wrap_df(df), self.new_columns) for df in batches
                ]
            else:
                return [wrap_df(df) for df in batches]
        return None
