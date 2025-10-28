from __future__ import annotations

import contextlib
from io import BytesIO, StringIO
from pathlib import Path
from typing import TYPE_CHECKING

from polars._utils.various import normalize_filepath
from polars._utils.wrap import wrap_df
from polars.datatypes import N_INFER_DEFAULT

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyDataFrame

if TYPE_CHECKING:
    from io import IOBase

    from polars import DataFrame
    from polars._typing import SchemaDefinition


def read_json(
    source: str | Path | IOBase | bytes,
    *,
    schema: SchemaDefinition | None = None,
    schema_overrides: SchemaDefinition | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
) -> DataFrame:
    """
    Read into a DataFrame from a JSON file.

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

    See Also
    --------
    read_ndjson

    Examples
    --------
    >>> from io import StringIO
    >>> json_str = '[{"foo":1,"bar":6},{"foo":2,"bar":7},{"foo":3,"bar":8}]'
    >>> pl.read_json(StringIO(json_str))
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

    With the schema defined.

    >>> pl.read_json(StringIO(json_str), schema={"foo": pl.Int64, "bar": pl.Float64})
    shape: (3, 2)
    ┌─────┬─────┐
    │ foo ┆ bar │
    │ --- ┆ --- │
    │ i64 ┆ f64 │
    ╞═════╪═════╡
    │ 1   ┆ 6.0 │
    │ 2   ┆ 7.0 │
    │ 3   ┆ 8.0 │
    └─────┴─────┘
    """
    if isinstance(source, StringIO):
        source = BytesIO(source.getvalue().encode())
    elif isinstance(source, (str, Path)):
        source = normalize_filepath(source)

    pydf = PyDataFrame.read_json(
        source,
        infer_schema_length=infer_schema_length,
        schema=schema,
        schema_overrides=schema_overrides,
    )
    return wrap_df(pydf)
