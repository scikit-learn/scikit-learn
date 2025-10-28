from __future__ import annotations

import os
import re
import warnings
from collections import defaultdict
from collections.abc import Sequence
from datetime import time
from glob import glob
from io import BufferedReader, BytesIO, StringIO, TextIOWrapper
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Callable, NoReturn, overload

import polars._reexport as pl
from polars import from_arrow
from polars import functions as F
from polars._dependencies import _PYARROW_AVAILABLE, import_optional
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    issue_deprecation_warning,
)
from polars._utils.various import deduplicate_names, normalize_filepath, parse_version
from polars.datatypes import (
    N_INFER_DEFAULT,
    Boolean,
    Date,
    Datetime,
    Duration,
    Int64,
    Null,
    String,
    Time,
    UInt8,
)
from polars.datatypes.group import FLOAT_DTYPES, INTEGER_DTYPES, NUMERIC_DTYPES
from polars.exceptions import (
    ModuleUpgradeRequiredError,
    NoDataError,
    ParameterCollisionError,
)
from polars.functions import concat
from polars.io._utils import looks_like_url, process_file_url
from polars.io.csv.functions import read_csv

if TYPE_CHECKING:
    from typing import Literal

    from polars._typing import ExcelSpreadsheetEngine, FileSource, SchemaDict


def _sources(source: FileSource) -> tuple[Any, bool]:
    """Unpack any glob patterns, standardise file paths."""
    read_multiple_workbooks = True
    sources: list[Any] = []

    if isinstance(source, memoryview):
        source = source.tobytes()
    if not isinstance(source, Sequence) or isinstance(source, (bytes, str)):
        read_multiple_workbooks = False
        source = [source]  # type: ignore[assignment]

    for src in source:  # type: ignore[union-attr]
        if isinstance(src, (str, os.PathLike)) and not Path(src).exists():
            src = os.path.expanduser(str(src))  # noqa: PTH111
            if looks_like_url(src):
                sources.append(src)
                continue
            sources.extend(files := glob(src, recursive=True))  # noqa: PTH207
            if not files:
                msg = f"no workbook found at path {src!r}"
                raise FileNotFoundError(msg)
            read_multiple_workbooks = True
        else:
            if isinstance(src, os.PathLike):
                src = str(src)
            sources.append(src)

    return sources, read_multiple_workbooks


def _standardize_duplicates(s: str) -> str:
    """Standardize columns with '_duplicated_n' names."""
    return re.sub(r"_duplicated_(\d+)", repl=r"\1", string=s)


def _unpack_read_results(
    frames: list[pl.DataFrame] | list[dict[str, pl.DataFrame]],
    *,
    read_multiple_workbooks: bool,
) -> Any:
    if not frames:
        msg = "no data found in the given workbook(s) and sheet(s)"
        raise NoDataError(msg)

    if not read_multiple_workbooks:
        # one sheet from one workbook
        return frames[0]

    if isinstance(frames[0], pl.DataFrame):
        # one sheet from multiple workbooks
        return concat(frames, how="vertical_relaxed")  # type: ignore[type-var]
    else:
        # multiple sheets from multiple workbooks
        sheet_frames = defaultdict(list)
        for res in frames:
            for sheet, df in res.items():  # type: ignore[union-attr]
                sheet_frames[sheet].append(df)
        return {k: concat(v, how="vertical_relaxed") for k, v in sheet_frames.items()}


@overload
def read_excel(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: str,
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    read_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_excel(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: None = ...,
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    read_options: dict[str, Any] | None = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_excel(
    source: FileSource,
    *,
    sheet_id: int,
    sheet_name: str,
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    read_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> NoReturn: ...


# note: 'ignore' required as mypy thinks that the return value for
# Literal[0] overlaps with the return value for other integers
@overload  # type: ignore[overload-overlap]
def read_excel(
    source: FileSource,
    *,
    sheet_id: Literal[0] | Sequence[int],
    sheet_name: None = ...,
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    read_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> dict[str, pl.DataFrame]: ...


@overload
def read_excel(
    source: FileSource,
    *,
    sheet_id: int,
    sheet_name: None = ...,
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    read_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_excel(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: list[str] | tuple[str],
    table_name: str | None = ...,
    engine: ExcelSpreadsheetEngine = ...,
    engine_options: dict[str, Any] | None = ...,
    read_options: dict[str, Any] | None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | str | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> dict[str, pl.DataFrame]: ...


@deprecate_renamed_parameter("xlsx2csv_options", "engine_options", version="0.20.6")
@deprecate_renamed_parameter("read_csv_options", "read_options", version="0.20.7")
def read_excel(
    source: FileSource,
    *,
    sheet_id: int | Sequence[int] | None = None,
    sheet_name: str | list[str] | tuple[str] | None = None,
    table_name: str | None = None,
    engine: ExcelSpreadsheetEngine = "calamine",
    engine_options: dict[str, Any] | None = None,
    read_options: dict[str, Any] | None = None,
    has_header: bool = True,
    columns: Sequence[int] | Sequence[str] | str | None = None,
    schema_overrides: SchemaDict | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    include_file_paths: str | None = None,
    drop_empty_rows: bool = True,
    drop_empty_cols: bool = True,
    raise_if_empty: bool = True,
) -> pl.DataFrame | dict[str, pl.DataFrame]:
    """
    Read Excel spreadsheet data into a DataFrame.

    .. versionadded:: 1.20
        Support loading data from named table objects with `table_name` parameter.
    .. versionadded:: 1.18
        Support loading data from a list (or glob pattern) of multiple workbooks.
    .. versionchanged:: 1.0
        Default engine is now "calamine" (was "xlsx2csv").
    .. versionchanged:: 0.20.7
        The `read_csv_options` parameter was renamed `read_options`.
    .. versionchanged:: 0.20.6
        The `xlsx2csv_options` parameter was renamed `engine_options`.

    Parameters
    ----------
    source
        Path(s) to a file or a file-like object (by "file-like object" we refer to
        objects that have a `read()` method, such as a file handler like the builtin
        `open` function, or a `BytesIO` instance). For file-like objects, the stream
        position may not be updated after reading.
    sheet_id
        Sheet number(s) to convert (set `0` to load all sheets as DataFrames) and
        return a `{sheetname:frame,}` dict. (Defaults to `1` if neither this nor
        `sheet_name` are specified). Can also take a sequence of sheet numbers.
    sheet_name
        Sheet name(s) to convert; cannot be used in conjunction with `sheet_id`. If
        more than one is given then a `{sheetname:frame,}` dict is returned.
    table_name
        Name of a specific table to read; note that table names are unique across
        the workbook, so additionally specifying a sheet id or name is optional;
        if one of those parameters *is* specified, an error will be raised if
        the named table is not found in that particular sheet.
    engine : {'calamine', 'openpyxl', 'xlsx2csv'}
        Library used to parse the spreadsheet file; defaults to "calamine".

        * "calamine": this engine can be used for reading all major types of Excel
          Workbook (`.xlsx`, `.xlsb`, `.xls`) and is dramatically faster than the
          other options, using the `fastexcel` module to bind the Rust-based Calamine
          parser.
        * "openpyxl": this engine is significantly slower than both `calamine` and
          `xlsx2csv`, but can provide a useful fallback if you are otherwise unable
          to read data from your workbook.
        * "xlsx2csv": converts the data to an in-memory CSV before using the native
          polars `read_csv` method to parse the result.
    engine_options
        Additional options passed to the underlying engine's primary parsing
        constructor (given below), if supported:

        * "calamine": n/a (can only provide `read_options`)
        * "openpyxl": `load_workbook <https://openpyxl.readthedocs.io/en/stable/api/openpyxl.reader.excel.html#openpyxl.reader.excel.load_workbook>`_
        * "xlsx2csv": `Xlsx2csv <https://github.com/dilshod/xlsx2csv/blob/f35734aa453d65102198a77e7b8cd04928e6b3a2/xlsx2csv.py#L157>`_
    read_options
        Options passed to the underlying engine method that reads the sheet data.
        Where supported, this allows for additional control over parsing. The
        specific read methods associated with each engine are:

        * "calamine": `load_sheet_by_name <https://fastexcel.toucantoco.dev/fastexcel.html#ExcelReader.load_sheet_by_name>`_
          (or `load_table <https://fastexcel.toucantoco.dev/fastexcel.html#ExcelReader.load_table>`_
          if using the `table_name` parameter).
        * "openpyxl": n/a (can only provide `engine_options`)
        * "xlsx2csv": see :meth:`read_csv`
    has_header
        Indicate if the first row of the table data is a header or not. If False,
        column names will be autogenerated in the following format: `column_x`, with
        `x` being an enumeration over every column in the dataset, starting at 1.
    columns
        Columns to read from the sheet; if not specified, all columns are read. Can
        be given as a sequence of column names or indices, or a single column name.
    schema_overrides
        Support type specification or override of one or more columns.
    infer_schema_length
        The maximum number of rows to scan for schema inference. If set to `None`, the
        entire dataset is scanned to determine the dtypes, which can slow parsing for
        large workbooks. Note that only the "calamine" and "xlsx2csv" engines support
        this parameter.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
    drop_empty_rows
        Indicate whether to omit empty rows when reading data into the DataFrame.
    drop_empty_cols
        Indicate whether to omit empty columns (with no headers) when reading data into
        the DataFrame (note that empty column identification may vary depending on the
        underlying engine being used).
    raise_if_empty
        When there is no data in the sheet,`NoDataError` is raised. If this parameter
        is set to False, an empty DataFrame (with no columns) is returned instead.

    Returns
    -------
    DataFrame
        If reading a single sheet.
    dict
        If reading multiple sheets, a "{sheetname: DataFrame, ...}" dict is returned.

    See Also
    --------
    read_ods

    Notes
    -----
    * Where possible, prefer the default "calamine" engine for reading Excel Workbooks,
      as it is significantly faster than the other options.
    * When using the `xlsx2csv` engine the target Excel sheet is first converted
      to CSV using `xlsx2csv.Xlsx2csv(source).convert()` and then parsed with Polars'
      :func:`read_csv` function. You can pass additional options to `read_options`
      to influence this part of the parsing pipeline.
    * If you want to read multiple sheets and set *different* options (`read_options`,
      `schema_overrides`, etc), you should make separate calls as the options are set
      globally, not on a per-sheet basis.

    Examples
    --------
    Read the "data" worksheet from an Excel file into a DataFrame.

    >>> pl.read_excel(
    ...     source="test.xlsx",
    ...     sheet_name="data",
    ... )  # doctest: +SKIP

    If the correct dtypes can't be determined, use the `schema_overrides` parameter
    to specify them, or increase the inference length with `infer_schema_length`.

    >>> pl.read_excel(
    ...     source="test.xlsx",
    ...     schema_overrides={"dt": pl.Date},
    ...     infer_schema_length=None,
    ... )  # doctest: +SKIP

    Using the `xlsx2csv` engine, read table data from sheet 3 in an Excel workbook as a
    DataFrame while skipping empty lines in the sheet. As sheet 3 does not have a header
    row, you can pass the necessary additional settings for this to the `read_options`
    parameter; these will be passed to :func:`read_csv`.

    >>> pl.read_excel(
    ...     source="test.xlsx",
    ...     sheet_id=3,
    ...     engine="xlsx2csv",
    ...     engine_options={"skip_empty_lines": True},
    ...     read_options={"has_header": False, "new_columns": ["a", "b", "c"]},
    ... )  # doctest: +SKIP
    """
    sources, read_multiple_workbooks = _sources(source)
    frames: list[pl.DataFrame] | list[dict[str, pl.DataFrame]] = [  # type: ignore[assignment]
        _read_spreadsheet(
            src,
            sheet_id=sheet_id,
            sheet_name=sheet_name,
            table_name=table_name,
            engine=engine,
            engine_options=engine_options,
            read_options=read_options,
            schema_overrides=schema_overrides,
            infer_schema_length=infer_schema_length,
            include_file_paths=include_file_paths,
            raise_if_empty=raise_if_empty,
            has_header=has_header,
            columns=columns,
            drop_empty_rows=drop_empty_rows,
            drop_empty_cols=drop_empty_cols,
        )
        for src in sources
    ]
    return _unpack_read_results(
        frames=frames,
        read_multiple_workbooks=read_multiple_workbooks,
    )


@overload
def read_ods(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: str,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_ods(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_ods(
    source: FileSource,
    *,
    sheet_id: int,
    sheet_name: str,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> NoReturn: ...


@overload  # type: ignore[overload-overlap]
def read_ods(
    source: FileSource,
    *,
    sheet_id: Literal[0] | Sequence[int],
    sheet_name: None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> dict[str, pl.DataFrame]: ...


@overload
def read_ods(
    source: FileSource,
    *,
    sheet_id: int,
    sheet_name: None = ...,
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> pl.DataFrame: ...


@overload
def read_ods(
    source: FileSource,
    *,
    sheet_id: None = ...,
    sheet_name: list[str] | tuple[str],
    has_header: bool = ...,
    columns: Sequence[int] | Sequence[str] | None = ...,
    schema_overrides: SchemaDict | None = ...,
    infer_schema_length: int | None = ...,
    include_file_paths: str | None = ...,
    drop_empty_rows: bool = ...,
    drop_empty_cols: bool = ...,
    raise_if_empty: bool = ...,
) -> dict[str, pl.DataFrame]: ...


def read_ods(
    source: FileSource,
    *,
    sheet_id: int | Sequence[int] | None = None,
    sheet_name: str | list[str] | tuple[str] | None = None,
    has_header: bool = True,
    columns: Sequence[int] | Sequence[str] | None = None,
    schema_overrides: SchemaDict | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    include_file_paths: str | None = None,
    drop_empty_rows: bool = True,
    drop_empty_cols: bool = True,
    raise_if_empty: bool = True,
) -> pl.DataFrame | dict[str, pl.DataFrame]:
    """
    Read OpenOffice (ODS) spreadsheet data into a DataFrame.

    Parameters
    ----------
    source
        Path to a file or a file-like object (by "file-like object" we refer to objects
        that have a `read()` method, such as a file handler like the builtin `open`
        function, or a `BytesIO` instance). For file-like objects, the stream position
        may not be updated accordingly after reading.
    sheet_id
        Sheet number(s) to convert, starting from 1 (set `0` to load *all* worksheets
        as DataFrames) and return a `{sheetname:frame,}` dict. (Defaults to `1` if
        neither this nor `sheet_name` are specified). Can also take a sequence of sheet
        numbers.
    sheet_name
        Sheet name(s) to convert; cannot be used in conjunction with `sheet_id`. If
        more than one is given then a `{sheetname:frame,}` dict is returned.
    has_header
        Indicate if the first row of the table data is a header or not. If False,
        column names will be autogenerated in the following format: `column_x`, with
        `x` being an enumeration over every column in the dataset, starting at 1.
    columns
        Columns to read from the sheet; if not specified, all columns are read. Can
        be given as a sequence of column names or indices.
    schema_overrides
        Support type specification or override of one or more columns.
    infer_schema_length
        The maximum number of rows to scan for schema inference. If set to `None`, the
        entire dataset is scanned to determine the dtypes, which can slow parsing for
        large workbooks.
    include_file_paths
        Include the path of the source file(s) as a column with this name.
    drop_empty_rows
        Indicate whether to omit empty rows when reading data into the DataFrame.
    drop_empty_cols
        Indicate whether to omit empty columns (with no headers) when reading data into
        the DataFrame (note that empty column identification may vary depending on the
        underlying engine being used).
    raise_if_empty
        When there is no data in the sheet,`NoDataError` is raised. If this parameter
        is set to False, an empty DataFrame (with no columns) is returned instead.

    Returns
    -------
    DataFrame, or a `{sheetname: DataFrame, ...}` dict if reading multiple sheets.

    See Also
    --------
    read_excel

    Examples
    --------
    Read the "data" worksheet from an OpenOffice spreadsheet file into a DataFrame.

    >>> pl.read_ods(
    ...     source="test.ods",
    ...     sheet_name="data",
    ... )  # doctest: +SKIP

    If the correct dtypes can't be determined, use the `schema_overrides` parameter
    to specify them, or increase the inference length with `infer_schema_length`.

    >>> pl.read_ods(
    ...     source="test.ods",
    ...     sheet_id=3,
    ...     schema_overrides={"dt": pl.Date},
    ...     raise_if_empty=False,
    ... )  # doctest: +SKIP
    """
    sources, read_multiple_workbooks = _sources(source)
    frames: list[pl.DataFrame] | list[dict[str, pl.DataFrame]] = [  # type: ignore[assignment]
        _read_spreadsheet(
            src,
            sheet_id=sheet_id,
            sheet_name=sheet_name,
            table_name=None,
            engine="calamine",
            engine_options={},
            read_options=None,
            schema_overrides=schema_overrides,
            infer_schema_length=infer_schema_length,
            include_file_paths=include_file_paths,
            raise_if_empty=raise_if_empty,
            drop_empty_rows=drop_empty_rows,
            drop_empty_cols=drop_empty_cols,
            has_header=has_header,
            columns=columns,
        )
        for src in sources
    ]
    return _unpack_read_results(
        frames=frames,
        read_multiple_workbooks=read_multiple_workbooks,
    )


def _read_spreadsheet(
    source: str | IO[bytes] | bytes,
    *,
    sheet_id: int | Sequence[int] | None,
    sheet_name: str | Sequence[str] | None,
    table_name: str | None,
    engine: ExcelSpreadsheetEngine,
    engine_options: dict[str, Any] | None = None,
    read_options: dict[str, Any] | None = None,
    schema_overrides: SchemaDict | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    include_file_paths: str | None = None,
    columns: Sequence[int] | Sequence[str] | str | None = None,
    has_header: bool = True,
    raise_if_empty: bool = True,
    drop_empty_rows: bool = True,
    drop_empty_cols: bool = True,
) -> pl.DataFrame | dict[str, pl.DataFrame]:
    if isinstance(source, str):
        source = normalize_filepath(source)
        if looks_like_url(source):
            source = process_file_url(source)

    if isinstance(columns, str):
        columns = [columns]

    read_options = _get_read_options(
        read_options,
        engine=engine,
        columns=columns,
        has_header=has_header,
        infer_schema_length=infer_schema_length,
    )
    engine_options = (engine_options or {}).copy()
    schema_overrides = dict(schema_overrides or {})

    # establish the reading function, parser, and available worksheets
    reader_fn, parser, worksheets = _initialise_spreadsheet_parser(
        engine, source, engine_options
    )
    try:
        # parse data from the indicated sheet(s)
        sheet_names, return_multiple_sheets = _get_sheet_names(
            sheet_id, sheet_name, table_name, worksheets
        )
        parsed_sheets = {
            name: reader_fn(
                parser=parser,
                sheet_name=name,
                schema_overrides=schema_overrides,
                read_options=read_options,
                raise_if_empty=raise_if_empty,
                columns=columns,
                table_name=table_name,
                drop_empty_rows=drop_empty_rows,
                drop_empty_cols=drop_empty_cols,
            )
            for name in sheet_names
        }
    finally:
        if hasattr(parser, "close"):
            parser.close()

    if not parsed_sheets:
        param, value = ("id", sheet_id) if sheet_name is None else ("name", sheet_name)
        msg = f"no matching sheets found when `sheet_{param}` is {value!r}"
        raise ValueError(msg)

    if include_file_paths:
        workbook = source if isinstance(source, str) else "in-mem"
        parsed_sheets = {
            name: frame.with_columns(F.lit(workbook).alias(include_file_paths))
            for name, frame in parsed_sheets.items()
        }
    if return_multiple_sheets:
        return parsed_sheets
    return next(iter(parsed_sheets.values()))


def _get_read_options(
    read_options: dict[str, Any] | None,
    *,
    engine: ExcelSpreadsheetEngine,
    columns: Sequence[int] | Sequence[str] | None,
    infer_schema_length: int | None,
    has_header: bool,
) -> dict[str, Any]:
    """Normalise top-level parameters to engine-specific 'read_options' dict."""
    read_options = (read_options or {}).copy()

    if engine == "calamine":
        if ("use_columns" in read_options) and columns:
            msg = 'cannot specify both `columns` and `read_options["use_columns"]`'
            raise ParameterCollisionError(msg)
        elif read_options.get("header_row") is not None and has_header is False:
            msg = 'the values of `has_header` and `read_options["header_row"]` are not compatible'
            raise ParameterCollisionError(msg)
        elif ("schema_sample_rows" in read_options) and (
            infer_schema_length != N_INFER_DEFAULT
        ):
            msg = 'cannot specify both `infer_schema_length` and `read_options["schema_sample_rows"]`'
            raise ParameterCollisionError(msg)

        read_options["schema_sample_rows"] = infer_schema_length
        if has_header is False and "header_row" not in read_options:
            read_options["header_row"] = None

    elif engine == "xlsx2csv":
        if ("columns" in read_options) and columns:
            msg = 'cannot specify both `columns` and `read_options["columns"]`'
            raise ParameterCollisionError(msg)
        elif (
            "has_header" in read_options
            and read_options["has_header"] is not has_header
        ):
            msg = 'the values of `has_header` and `read_options["has_header"]` are not compatible'
            raise ParameterCollisionError(msg)
        elif ("infer_schema_length" in read_options) and (
            infer_schema_length != N_INFER_DEFAULT
        ):
            msg = 'cannot specify both `infer_schema_length` and `read_options["infer_schema_length"]`'
            raise ParameterCollisionError(msg)

        read_options["infer_schema_length"] = infer_schema_length
        if "has_header" not in read_options:
            read_options["has_header"] = has_header
    else:
        read_options["infer_schema_length"] = infer_schema_length
        read_options["has_header"] = has_header

    return read_options


def _get_sheet_names(
    sheet_id: int | Sequence[int] | None,
    sheet_name: str | Sequence[str] | None,
    table_name: str | None,
    worksheets: list[dict[str, Any]],
) -> tuple[list[str], bool]:
    """Establish sheets to read; indicate if we are returning a dict frames."""
    if sheet_id is not None and sheet_name is not None:
        msg = f"cannot specify both `sheet_name` ({sheet_name!r}) and `sheet_id` ({sheet_id!r})"
        raise ValueError(msg)

    sheet_names = []
    if sheet_id is None and sheet_name is None:
        name = None if table_name else worksheets[0]["name"]
        sheet_names.append(name)
        return_multiple_sheets = False
    elif sheet_id == 0:
        sheet_names.extend(ws["name"] for ws in worksheets)
        return_multiple_sheets = True
    else:
        return_multiple_sheets = (
            (isinstance(sheet_name, Sequence) and not isinstance(sheet_name, str))
            or isinstance(sheet_id, Sequence)
            or sheet_id == 0
        )
        if names := (
            (sheet_name,) if isinstance(sheet_name, str) else sheet_name or ()
        ):
            known_sheet_names = {ws["name"] for ws in worksheets}
            for name in names:
                if name not in known_sheet_names:
                    msg = f"no matching sheet found when `sheet_name` is {name!r}"
                    raise ValueError(msg)
                sheet_names.append(name)
        else:
            ids = (sheet_id,) if isinstance(sheet_id, int) else sheet_id or ()
            sheet_names_by_idx = {
                idx: ws["name"]
                for idx, ws in enumerate(worksheets, start=1)
                if (sheet_id == 0 or ws["index"] in ids or ws["name"] in names)
            }
            for idx in ids:
                if (name := sheet_names_by_idx.get(idx)) is None:
                    msg = f"no matching sheet found when `sheet_id` is {idx}"
                    raise ValueError(msg)
                sheet_names.append(name)

    return sheet_names, return_multiple_sheets  # type: ignore[return-value]


def _initialise_spreadsheet_parser(
    engine: str | None,
    source: str | IO[bytes] | bytes,
    engine_options: dict[str, Any],
) -> tuple[Callable[..., pl.DataFrame], Any, list[dict[str, Any]]]:
    """Instantiate the indicated spreadsheet parser and establish related properties."""
    if isinstance(source, str) and not Path(source).exists():
        raise FileNotFoundError(source)

    if engine == "xlsx2csv":  # default
        xlsx2csv = import_optional("xlsx2csv")

        # establish sensible defaults for unset options
        for option, value in {
            "exclude_hidden_sheets": False,
            "skip_empty_lines": False,
            "skip_hidden_rows": False,
            "floatformat": "%f",
        }.items():
            engine_options.setdefault(option, value)

        if isinstance(source, bytes):
            source = BytesIO(source)

        parser = xlsx2csv.Xlsx2csv(source, **engine_options)
        sheets = parser.workbook.sheets
        return _read_spreadsheet_xlsx2csv, parser, sheets

    elif engine == "openpyxl":
        openpyxl = import_optional("openpyxl")
        if isinstance(source, bytes):
            source = BytesIO(source)

        parser = openpyxl.load_workbook(source, data_only=True, **engine_options)
        sheets = [{"index": i + 1, "name": ws.title} for i, ws in enumerate(parser)]
        return _read_spreadsheet_openpyxl, parser, sheets

    elif engine == "calamine":
        fastexcel = import_optional("fastexcel", min_version="0.7.0")
        reading_bytesio, reading_bytes = (
            isinstance(source, BytesIO),
            isinstance(source, bytes),
        )
        if (reading_bytesio or reading_bytes) and parse_version(
            module_version := fastexcel.__version__
        ) < (0, 10):
            msg = f"`fastexcel` >= 0.10 is required to read bytes; found {module_version})"
            raise ModuleUpgradeRequiredError(msg)

        if reading_bytesio:
            source = source.getvalue()  # type: ignore[union-attr]
        elif isinstance(source, (BufferedReader, TextIOWrapper)):
            if "b" not in source.mode:
                msg = f"file {source.name!r} must be opened in binary mode"
                raise OSError(msg)
            elif (filename := source.name) and Path(filename).exists():
                source = filename
            else:
                source = source.read()

        parser = fastexcel.read_excel(source, **engine_options)
        sheets = [
            {"index": i + 1, "name": nm} for i, nm in enumerate(parser.sheet_names)
        ]
        return _read_spreadsheet_calamine, parser, sheets

    msg = f"unrecognized engine: {engine!r}"
    raise NotImplementedError(msg)


def _csv_buffer_to_frame(
    csv: StringIO,
    *,
    separator: str,
    read_options: dict[str, Any],
    schema_overrides: SchemaDict | None,
    drop_empty_rows: bool,
    drop_empty_cols: bool,
    raise_if_empty: bool,
) -> pl.DataFrame:
    """Translate StringIO buffer containing delimited data as a DataFrame."""
    # handle (completely) empty sheet data
    if csv.tell() == 0:
        return _empty_frame(raise_if_empty)

    # otherwise rewind the buffer and parse as csv
    csv.seek(0)

    if read_options is None:
        read_options = {}

    date_cols = []
    if schema_overrides:
        if csv_dtypes := read_options.get("dtypes", {}):
            issue_deprecation_warning(
                "the `dtypes` parameter for `read_csv` is deprecated. It has been renamed to `schema_overrides`.",
                version="0.20.31",
            )

        csv_schema_overrides = read_options.get("schema_overrides", csv_dtypes)
        if set(csv_schema_overrides).intersection(schema_overrides):
            msg = "cannot specify columns in both `schema_overrides` and `read_options['dtypes']`"
            raise ParameterCollisionError(msg)

        overrides, schema_overrides = {**csv_schema_overrides, **schema_overrides}, {}
        for nm, dtype in overrides.items():
            if dtype != Date:
                schema_overrides[nm] = dtype
            else:
                date_cols.append(nm)

        read_options = read_options.copy()
        read_options["schema_overrides"] = schema_overrides

    df = _drop_null_data(
        df=read_csv(
            csv,
            separator=separator,
            **read_options,
        ),
        raise_if_empty=raise_if_empty,
        drop_empty_rows=drop_empty_rows,
        drop_empty_cols=drop_empty_cols,
    )
    if date_cols:
        date_casts, schema = {}, df.schema
        for nm in date_cols:
            if schema[nm] == String:
                date_casts[nm] = (
                    F.col(nm)
                    .str.replace(r"(?:[ T]00:00:00(?:\.0+)?)$", "")
                    .str.to_date()
                )
        if date_casts:
            df = df.with_columns(**date_casts)
    return df


def _drop_null_data(
    df: pl.DataFrame,
    *,
    raise_if_empty: bool,
    drop_empty_rows: bool = True,
    drop_empty_cols: bool = True,
) -> pl.DataFrame:
    """If DataFrame contains columns/rows that contain only nulls, drop them."""
    null_cols: list[str] = []
    if drop_empty_cols:
        for col_name in df.columns:
            # note that if multiple unnamed columns are found then all but the first one
            # will be named as "_duplicated_{n}" (or "__UNNAMED__{n}" from calamine)
            if col_name == "" or re.match(r"(_duplicated_|__UNNAMED__)\d+$", col_name):
                col = df[col_name]
                if (
                    col.dtype == Null
                    or col.null_count() == df.height
                    or (
                        col.dtype in NUMERIC_DTYPES
                        and col.replace(0, None).null_count() == df.height
                    )
                ):
                    null_cols.append(col_name)
        if null_cols:
            df = df.drop(*null_cols)

    if df.height == df.width == 0:
        return _empty_frame(raise_if_empty)
    if drop_empty_rows:
        return df.filter(~F.all_horizontal(F.all().is_null()))
    return df


def _empty_frame(raise_if_empty: bool) -> pl.DataFrame:  # noqa: FBT001
    if raise_if_empty:
        msg = (
            "empty Excel sheet"
            "\n\nIf you want to read this as an empty DataFrame, set `raise_if_empty=False`."
        )
        raise NoDataError(msg)
    return pl.DataFrame()


def _reorder_columns(
    df: pl.DataFrame, columns: Sequence[int] | Sequence[str] | None
) -> pl.DataFrame:
    if columns:
        from polars.selectors import by_index, by_name

        cols = by_index(*columns) if isinstance(columns[0], int) else by_name(*columns)
        df = df.select(cols)
    return df


def _read_spreadsheet_calamine(
    parser: Any,
    *,
    sheet_name: str | None,
    read_options: dict[str, Any],
    schema_overrides: SchemaDict | None,
    columns: Sequence[int] | Sequence[str] | None,
    table_name: str | None = None,
    drop_empty_rows: bool,
    drop_empty_cols: bool,
    raise_if_empty: bool,
) -> pl.DataFrame:
    # if we have 'schema_overrides' and a more recent version of `fastexcel`
    # we can pass translated dtypes to the engine to refine the initial parse
    fastexcel = import_optional("fastexcel")
    fastexcel_version = parse_version(original_version := fastexcel.__version__)

    if fastexcel_version < (0, 9) and "schema_sample_rows" in read_options:
        msg = f"a more recent version of `fastexcel` is required for 'schema_sample_rows' (>= 0.9; found {original_version})"
        raise ModuleUpgradeRequiredError(msg)
    if fastexcel_version < (0, 10, 2) and "use_columns" in read_options:
        msg = f"a more recent version of `fastexcel` is required for 'use_columns' (>= 0.10.2; found {original_version})"
        raise ModuleUpgradeRequiredError(msg)
    if table_name and fastexcel_version < (0, 12):
        msg = f"a more recent version of `fastexcel` is required for 'table_name' (>= 0.12.0; found {original_version})"
        raise ValueError(msg)

    if columns:
        if not isinstance(columns, list):
            columns = list(columns)  # type: ignore[assignment]
        read_options["use_columns"] = columns

    schema_overrides = schema_overrides or {}
    if read_options.get("schema_sample_rows") == 0:
        # ref: https://github.com/ToucanToco/fastexcel/issues/236
        del read_options["schema_sample_rows"]
        read_options["dtypes"] = (
            "string"
            if fastexcel_version >= (0, 12, 1)
            else dict.fromkeys(range(16384), "string")
        )
    elif schema_overrides and fastexcel_version >= (0, 10):
        parser_dtypes = read_options.get("dtypes", {})
        for name, dtype in schema_overrides.items():
            if name not in parser_dtypes:
                if (base_dtype := dtype.base_type()) in INTEGER_DTYPES:
                    parser_dtypes[name] = "int"
                elif base_dtype in FLOAT_DTYPES:
                    parser_dtypes[name] = "float"
                elif base_dtype == String:
                    parser_dtypes[name] = "string"
                elif base_dtype == Duration:
                    parser_dtypes[name] = "duration"
                elif base_dtype == Boolean:
                    parser_dtypes[name] = "boolean"

        read_options["dtypes"] = parser_dtypes

    if fastexcel_version < (0, 11, 2):
        ws = parser.load_sheet_by_name(name=sheet_name, **read_options)
        df = ws.to_polars()
    else:
        if table_name:
            if col_names := read_options.get("use_columns"):
                selected_col_names = set(col_names)
                read_options["use_columns"] = lambda col: col.name in selected_col_names

            xl_table = parser.load_table(table_name, **read_options)

            if sheet_name and sheet_name != xl_table.sheet_name:
                msg = f"table named {table_name!r} not found in sheet {sheet_name!r}"
                raise RuntimeError(msg)
            df = xl_table.to_polars()

        elif _PYARROW_AVAILABLE:
            # eager loading is faster / more memory-efficient, but requires pyarrow
            ws_arrow = parser.load_sheet_eager(sheet_name, **read_options)
            df = from_arrow(ws_arrow)
        else:
            ws_arrow = parser.load_sheet(sheet_name, **read_options)
            df = from_arrow(ws_arrow)

        if read_options.get("header_row", False) is None and not read_options.get(
            "column_names"
        ):
            df.columns = [f"column_{i}" for i in range(1, df.width + 1)]

    df = _drop_null_data(
        df,
        raise_if_empty=raise_if_empty,
        drop_empty_rows=drop_empty_rows,
        drop_empty_cols=drop_empty_cols,
    )

    # note: even if we applied parser dtypes we still re-apply schema_overrides
    # natively as we can refine integer/float types, temporal precision, etc.
    if schema_overrides:
        lf, schema = df.lazy(), df.schema
        str_to_temporal, updated_overrides = [], {}
        for nm, tp in schema_overrides.items():
            if schema[nm] != String:
                updated_overrides[nm] = tp
            elif tp == Datetime:
                str_to_temporal.append(
                    F.col(nm).str.to_datetime(
                        time_unit=getattr(tp, "time_unit", None),
                        time_zone=getattr(tp, "time_zone", None),
                    )
                )
            elif tp == Date:
                dt_str = F.col(nm).str.replace(r"(?:[ T]00:00:00(?:\.0+)?)$", "")
                str_to_temporal.append(dt_str.str.to_date())
            elif tp == Time:
                str_to_temporal.append(F.col(nm).str.to_time())
            else:
                updated_overrides[nm] = tp

        if str_to_temporal:
            lf = lf.with_columns(*str_to_temporal)
        if updated_overrides:
            lf = lf.cast(dtypes=updated_overrides)
        df = lf.collect()

    # standardise on string dtype for null columns in empty frame
    if df.is_empty():
        df = df.cast({Null: String})

    # further refine dtypes
    type_checks = []
    for c, dtype in df.schema.items():
        if c not in schema_overrides:
            # may read integer data as float; cast back to int where possible.
            if dtype in FLOAT_DTYPES:
                check_cast = [
                    F.col(c).floor().eq_missing(F.col(c)) & F.col(c).is_not_nan(),
                    F.col(c).cast(Int64),
                ]
                type_checks.append(check_cast)
            # do a similar check for datetime columns that have only 00:00:00 times.
            elif dtype == Datetime:
                check_cast = [
                    F.col(c).dt.time().eq(time(0, 0, 0)),
                    F.col(c).cast(Date),
                ]
                type_checks.append(check_cast)

    if type_checks:
        apply_cast = df.select(d[0].all(ignore_nulls=True) for d in type_checks).row(0)
        if downcast := [
            cast for apply, (_, cast) in zip(apply_cast, type_checks) if apply
        ]:
            df = df.with_columns(*downcast)

    return df


def _read_spreadsheet_openpyxl(
    parser: Any,
    *,
    sheet_name: str | None,
    read_options: dict[str, Any],
    schema_overrides: SchemaDict | None,
    columns: Sequence[int] | Sequence[str] | None,
    table_name: str | None = None,
    drop_empty_rows: bool,
    drop_empty_cols: bool,
    raise_if_empty: bool,
) -> pl.DataFrame:
    """Use the 'openpyxl' library to read data from the given worksheet."""
    infer_schema_length = read_options.pop("infer_schema_length", None)
    has_header = read_options.pop("has_header", True)
    schema_overrides = schema_overrides or {}
    no_inference = infer_schema_length == 0
    header: list[str | None] = []

    if table_name and not sheet_name:
        sheet_name, n_tables = None, 0
        for sheet in parser.worksheets:
            n_tables += 1
            if table_name in sheet.tables:
                ws, sheet_name = sheet, sheet.title
                break
        if sheet_name is None:
            msg = (
                f"table named {table_name!r} not found in sheet {sheet_name!r}"
                if n_tables
                else f"no named tables found in sheet {sheet_name!r} (looking for {table_name!r})"
            )
            raise RuntimeError(msg)
    else:
        ws = parser[sheet_name]

    # prefer detection of actual table objects; otherwise read
    # data in the used worksheet range, dropping null columns
    if tables := getattr(ws, "tables", None):
        table = tables[table_name] if table_name else next(iter(tables.values()))
        rows = list(ws[table.ref])
        if not rows:
            return _empty_frame(raise_if_empty)
        if has_header:
            header.extend(cell.value for cell in rows.pop(0))
        else:
            header.extend(f"column_{n}" for n in range(1, len(rows[0]) + 1))
        if table.totalsRowCount:
            rows = rows[: -table.totalsRowCount]
        rows_iter = rows
    elif table_name:
        msg = f"no named tables found in sheet {sheet_name!r} (looking for {table_name!r})"
        raise RuntimeError(msg)
    else:
        if not has_header:
            if not (rows_iter := list(ws.iter_rows())):
                return _empty_frame(raise_if_empty)
            n_cols = len(rows_iter[0])
            header = [f"column_{n}" for n in range(1, n_cols + 1)]
        else:
            rows_iter = ws.iter_rows()
            for row in rows_iter:
                row_values = [cell.value for cell in row]
                if any(v is not None for v in row_values):
                    header.extend(row_values)
                    break

    dtype = String if no_inference else None
    series_data = []
    for name, column_data in zip(header, zip(*rows_iter)):
        if name or not drop_empty_cols:
            values = [cell.value for cell in column_data]
            if no_inference or (dtype := schema_overrides.get(name)) == String:  # type: ignore[assignment,arg-type]
                # note: if we initialise the series with mixed-type data (eg: str/int)
                # then the non-strings will become null, so we handle the cast here
                values = [str(v) if (v is not None) else v for v in values]

            if (tp := schema_overrides.get(name)) in (Date, Datetime, Time):  # type: ignore[operator,arg-type]
                s = pl.Series(name, values, strict=False)
                if s.dtype == String:
                    if tp == Datetime:
                        s = s.str.to_datetime(
                            time_unit=getattr(tp, "time_unit", None),
                            time_zone=getattr(tp, "time_zone", None),
                        )
                    elif tp == Date:
                        s = s.str.replace(
                            r"(?:[ T]00:00:00(?:\.0+)?)$", ""
                        ).str.to_date()
                    elif tp == Time:
                        s = s.str.to_time()
            else:
                s = pl.Series(name, values, dtype=dtype, strict=False)
            series_data.append(s)

    names = deduplicate_names(s.name for s in series_data)
    df = pl.DataFrame(
        dict(zip(names, series_data)),
        schema_overrides=schema_overrides,
        infer_schema_length=infer_schema_length,
        strict=False,
    )
    df = _drop_null_data(
        df,
        raise_if_empty=raise_if_empty,
        drop_empty_rows=drop_empty_rows,
        drop_empty_cols=drop_empty_cols,
    )
    df = _reorder_columns(df, columns)
    return df


def _read_spreadsheet_xlsx2csv(
    parser: Any,
    *,
    sheet_name: str | None,
    read_options: dict[str, Any],
    schema_overrides: SchemaDict | None,
    columns: Sequence[int] | Sequence[str] | None,
    table_name: str | None = None,
    drop_empty_rows: bool,
    drop_empty_cols: bool,
    raise_if_empty: bool,
) -> pl.DataFrame:
    """Use the 'xlsx2csv' library to read data from the given worksheet."""
    if table_name:
        msg = "the `table_name` parameter is not supported by the 'xlsx2csv' engine"
        raise ValueError(msg)

    csv_buffer = StringIO()
    with warnings.catch_warnings():
        # xlsx2csv version 0.8.4 throws a DeprecationWarning in Python 3.13
        # https://github.com/dilshod/xlsx2csv/pull/287
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        parser.convert(outfile=csv_buffer, sheetname=sheet_name)

    read_options.setdefault("truncate_ragged_lines", True)
    if columns:
        read_options["columns"] = columns

    cast_to_boolean = []
    if schema_overrides:
        for col, dtype in schema_overrides.items():
            if dtype == Boolean:
                schema_overrides[col] = UInt8  # type: ignore[index]
                cast_to_boolean.append(F.col(col).cast(Boolean))

    df = _csv_buffer_to_frame(
        csv_buffer,
        separator=",",
        read_options=read_options,
        schema_overrides=schema_overrides,
        raise_if_empty=raise_if_empty,
        drop_empty_rows=drop_empty_rows,
        drop_empty_cols=drop_empty_cols,
    )
    if cast_to_boolean:
        df = df.with_columns(*cast_to_boolean)

    df = df.rename(_standardize_duplicates)
    return _reorder_columns(df, columns)
