from __future__ import annotations

import io
import itertools
import re
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Any, Literal, overload

import polars._reexport as pl
from polars import functions as F
from polars._dependencies import _check_for_pyarrow
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._utils.construction.dataframe import (
    arrow_to_pydf,
    dict_to_pydf,
    numpy_to_pydf,
    pandas_to_pydf,
    sequence_to_pydf,
)
from polars._utils.construction.series import arrow_to_pyseries, pandas_to_pyseries
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    issue_deprecation_warning,
)
from polars._utils.pycapsule import is_pycapsule, pycapsule_to_frame
from polars._utils.various import (
    _cast_repr_strings_with_schema,
    issue_warning,
    qualified_type_name,
)
from polars._utils.wrap import wrap_df, wrap_s
from polars.datatypes import N_INFER_DEFAULT, Categorical, String
from polars.exceptions import NoDataError

if TYPE_CHECKING:
    from collections.abc import Mapping

    from polars import DataFrame, Series
    from polars._dependencies import numpy as np
    from polars._dependencies import torch
    from polars._typing import (
        ArrowArrayExportable,
        ArrowStreamExportable,
        Orientation,
        PolarsDataType,
        SchemaDefinition,
        SchemaDict,
    )
    from polars.interchange.protocol import SupportsInterchange


def from_dict(
    data: Mapping[str, Sequence[object] | Mapping[str, Sequence[object]] | Series],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
) -> DataFrame:
    """
    Construct a DataFrame from a dictionary of sequences.

    This operation clones data, unless you pass a `{str: pl.Series,}` dict.

    Parameters
    ----------
    data : dict of sequences
        Two-dimensional data represented as a dictionary. dict must contain
        Sequences.
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
        any dtypes inferred from the columns param will be overridden.
    strict : bool, default True
        Throw an error if any `data` value does not exactly match the given or inferred
        data type for that column. If set to `False`, values that do not match the data
        type are cast to that data type or, if casting is not possible, set to null
        instead.

    Returns
    -------
    DataFrame

    Examples
    --------
    >>> df = pl.from_dict({"a": [1, 2], "b": [3, 4]})
    >>> df
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    │ 2   ┆ 4   │
    └─────┴─────┘
    """
    return wrap_df(
        dict_to_pydf(
            data,
            schema=schema,
            schema_overrides=schema_overrides,
            strict=strict,
        )
    )


def from_dicts(
    data: Iterable[Mapping[str, Any]],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    infer_schema_length: int | None = N_INFER_DEFAULT,
) -> DataFrame:
    """
    Construct a DataFrame from a sequence of dictionaries. This operation clones data.

    Parameters
    ----------
    data
        Sequence with dictionaries mapping column name to value
    schema : Sequence of str, (str,DataType) pairs, or a {str:DataType,} dict
        The DataFrame schema may be declared in several ways:

        * As a dict of {name:type} pairs; if type is None, it will be auto-inferred.
        * As a list of column names; in this case types are automatically inferred.
        * As a list of (name,type) pairs; this is equivalent to the dictionary form.

        If a list of column names is supplied that does NOT match the names in the
        underlying data, the names given here will overwrite the actual fields in
        the order that they appear - however, in this case it is typically clearer
        to rename after loading the frame.

        If you want to drop some of the fields found in the input dictionaries, a
        *partial* schema can be declared, in which case omitted fields will not be
        loaded. Similarly, you can extend the loaded frame with empty columns by
        adding them to the schema.
    schema_overrides : dict, default None
        Support override of inferred types for one or more columns.
    strict : bool, default True
        Throw an error if any `data` value does not exactly match the given or inferred
        data type for that column. If set to `False`, values that do not match the data
        type are cast to that data type or, if casting is not possible, set to null
        instead.
    infer_schema_length
        The maximum number of rows to scan for schema inference.
        If set to `None`, the full data may be scanned *(this is slow)*.

    Returns
    -------
    DataFrame

    Examples
    --------
    >>> data = [{"a": 1, "b": 4}, {"a": 2, "b": 5}, {"a": 3, "b": 6}]
    >>> df = pl.from_dicts(data)
    >>> df
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 4   │
    │ 2   ┆ 5   │
    │ 3   ┆ 6   │
    └─────┴─────┘

    Declaring a partial `schema` will drop the omitted columns.

    >>> df = pl.from_dicts(data, schema={"a": pl.Int32})
    >>> df
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i32 │
    ╞═════╡
    │ 1   │
    │ 2   │
    │ 3   │
    └─────┘

    Can also use the `schema` param to extend the loaded columns with one
    or more additional (empty) columns that are not present in the input dicts:

    >>> pl.from_dicts(
    ...     data,
    ...     schema=["a", "b", "c", "d"],
    ...     schema_overrides={"c": pl.Float64, "d": pl.String},
    ... )
    shape: (3, 4)
    ┌─────┬─────┬──────┬──────┐
    │ a   ┆ b   ┆ c    ┆ d    │
    │ --- ┆ --- ┆ ---  ┆ ---  │
    │ i64 ┆ i64 ┆ f64  ┆ str  │
    ╞═════╪═════╪══════╪══════╡
    │ 1   ┆ 4   ┆ null ┆ null │
    │ 2   ┆ 5   ┆ null ┆ null │
    │ 3   ┆ 6   ┆ null ┆ null │
    └─────┴─────┴──────┴──────┘
    """
    if not data and not (schema or schema_overrides):
        msg = "no data, cannot infer schema"
        raise NoDataError(msg)

    return pl.DataFrame(
        data,
        schema=schema,
        schema_overrides=schema_overrides,
        strict=strict,
        infer_schema_length=infer_schema_length,
    )


def from_records(
    data: Sequence[Any],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    orient: Orientation | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
) -> DataFrame:
    """
    Construct a DataFrame from a sequence of sequences. This operation clones data.

    Note that this is slower than creating from columnar memory.

    Parameters
    ----------
    data : Sequence of sequences
        Two-dimensional data represented as a sequence of sequences.
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
        any dtypes inferred from the columns param will be overridden.
    strict : bool, default True
        Throw an error if any `data` value does not exactly match the given or inferred
        data type for that column. If set to `False`, values that do not match the data
        type are cast to that data type or, if casting is not possible, set to null
        instead.
    orient : {None, 'col', 'row'}
        Whether to interpret two-dimensional data as columns or as rows. If None,
        the orientation is inferred by matching the columns and data dimensions. If
        this does not yield conclusive results, column orientation is used.
    infer_schema_length
        The maximum number of rows to scan for schema inference.
        If set to `None`, the full data may be scanned *(this is slow)*.

    Returns
    -------
    DataFrame

    Examples
    --------
    >>> data = [[1, 2, 3], [4, 5, 6]]
    >>> df = pl.from_records(data, schema=["a", "b"])
    >>> df
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 4   │
    │ 2   ┆ 5   │
    │ 3   ┆ 6   │
    └─────┴─────┘
    """
    if not isinstance(data, Sequence):
        msg = (
            f"expected data of type Sequence, got {type(data).__name__!r}"
            "\n\nHint: Try passing your data to the DataFrame constructor instead,"
            " e.g. `pl.DataFrame(data)`."
        )
        raise TypeError(msg)

    return wrap_df(
        sequence_to_pydf(
            data,
            schema=schema,
            schema_overrides=schema_overrides,
            strict=strict,
            orient=orient,
            infer_schema_length=infer_schema_length,
        )
    )


def from_numpy(
    data: np.ndarray[Any, Any],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    orient: Orientation | None = None,
) -> DataFrame:
    """
    Construct a DataFrame from a NumPy ndarray. This operation clones data.

    Note that this is slower than creating from columnar memory.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Two-dimensional data represented as a NumPy ndarray.
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
        any dtypes inferred from the columns param will be overridden.
    orient : {None, 'col', 'row'}
        Whether to interpret two-dimensional data as columns or as rows. If None,
        the orientation is inferred by matching the columns and data dimensions. If
        this does not yield conclusive results, column orientation is used.

    Returns
    -------
    DataFrame

    Examples
    --------
    >>> import numpy as np
    >>> data = np.array([[1, 2, 3], [4, 5, 6]])
    >>> df = pl.from_numpy(data, schema=["a", "b"], orient="col")
    >>> df
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 4   │
    │ 2   ┆ 5   │
    │ 3   ┆ 6   │
    └─────┴─────┘
    """
    return wrap_df(
        numpy_to_pydf(
            data=data,
            schema=schema,
            schema_overrides=schema_overrides,
            orient=orient,
        )
    )


def from_torch(
    tensor: torch.Tensor,
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    orient: Orientation | None = None,
    force: bool = False,
) -> DataFrame:
    """
    Construct a DataFrame from a PyTorch Tensor.

    Parameters
    ----------
    tensor : :class:`torch.Tensor`
        A PyTorch `Tensor` object of one or more dimensions.
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
        any dtypes inferred from the columns param will be overridden.
    orient : {None, 'col', 'row'}
        Whether to interpret two-dimensional data as columns or as rows. If None,
        the orientation is inferred by matching the columns and data dimensions. If
        this does not yield conclusive results, column orientation is used.
    force : bool
        If False, the conversion is performed only if the Tensor is on CPU, does not
        require grad, does not have its conjugate bit set, and is of a dtype (and
        layout) that NumPy supports; this will typically be zero-copy. If True, it
        is equivalent to calling `.detach().cpu().resolve_conj().resolve_neg()`
        before passing the Tensor to Polars.

    Returns
    -------
    DataFrame

    Examples
    --------
    >>> import torch
    >>> data = torch.tensor(
    ...     [
    ...         [1234.5, 200.0, 3000.5],
    ...         [8000.0, 500.5, 6000.0],
    ...     ]
    ... )
    >>> df = pl.from_torch(
    ...     data,
    ...     schema=["colx", "coly", "colz"],
    ...     schema_overrides={"colz": pl.Float64},
    ... )
    >>> df
    shape: (2, 3)
    ┌────────┬───────┬────────┐
    │ colx   ┆ coly  ┆ colz   │
    │ ---    ┆ ---   ┆ ---    │
    │ f32    ┆ f32   ┆ f64    │
    ╞════════╪═══════╪════════╡
    │ 1234.5 ┆ 200.0 ┆ 3000.5 │
    │ 8000.0 ┆ 500.5 ┆ 6000.0 │
    └────────┴───────┴────────┘
    """
    return wrap_df(
        numpy_to_pydf(
            data=tensor.numpy(force=force),
            schema=schema,
            schema_overrides=schema_overrides,
            orient=orient,
        )
    )


# Note: we cannot @overload the typing (Series vs DataFrame) here, as pyarrow
# does not (yet?) implement any support for type hints; attempts to hint here
# will simply result in mypy inferring "Any", which isn't at all useful...


def from_arrow(
    data: (
        pa.Table
        | pa.Array
        | pa.ChunkedArray
        | pa.RecordBatch
        | Iterable[pa.RecordBatch | pa.Table]
        | ArrowArrayExportable
        | ArrowStreamExportable
    ),
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    rechunk: bool = True,
) -> DataFrame | Series:
    """
    Create a DataFrame or Series from an Arrow Table or Array.

    This operation will be zero copy for the most part. Types that are not
    supported by Polars may be cast to the closest supported type.

    Hint: You can also directly pass arrow tables to `pl.DataFrame()` / arrow
    arrays to `pl.Series()` if the output type is known to avoid typing issues.

    Parameters
    ----------
    data : :class:`pyarrow.Table`, :class:`pyarrow.Array`, one or more :class:`pyarrow.RecordBatch`
        Data representing an Arrow Table, Array, sequence of RecordBatches or Tables, or other
        object that supports the Arrow PyCapsule interface.
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
    rechunk : bool, default True
        Make sure that all data is in contiguous memory.

    Returns
    -------
    DataFrame or Series

    Examples
    --------
    Constructing a DataFrame from an Arrow Table:

    >>> import pyarrow as pa
    >>> data = pa.table({"a": [1, 2, 3], "b": [4, 5, 6]})
    >>> pl.from_arrow(data)
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 4   │
    │ 2   ┆ 5   │
    │ 3   ┆ 6   │
    └─────┴─────┘

    Constructing a Series from an Arrow Array:

    >>> import pyarrow as pa
    >>> data = pa.array([1, 2, 3])
    >>> pl.from_arrow(data, schema={"s": pl.Int32})
    shape: (3,)
    Series: 's' [i32]
    [
        1
        2
        3
    ]
    """  # noqa: W505
    if is_pycapsule(data) and not _check_for_pyarrow(data):
        return pycapsule_to_frame(
            data,
            schema=schema,
            schema_overrides=schema_overrides,
            rechunk=rechunk,
        )

    elif isinstance(data, (pa.Table, pa.RecordBatch)):
        return wrap_df(
            arrow_to_pydf(
                data=data,
                rechunk=rechunk,
                schema=schema,
                schema_overrides=schema_overrides,
            )
        )
    elif isinstance(data, (pa.Array, pa.ChunkedArray)):
        name = getattr(data, "_name", "") or ""
        s = wrap_s(arrow_to_pyseries(name, data, rechunk=rechunk))
        s = pl.DataFrame(
            data=s,
            schema=schema,
            schema_overrides=schema_overrides,
        ).to_series()
        return s if (name or schema or schema_overrides) else s.alias("")

    elif not data:
        return pl.DataFrame(
            schema=schema,
            schema_overrides=schema_overrides,
        )

    if isinstance(data, Iterable):
        pa_table = pa.Table.from_batches(
            itertools.chain.from_iterable(
                (b.to_batches() if isinstance(b, pa.Table) else [b]) for b in data
            )
        )
        return wrap_df(
            arrow_to_pydf(
                data=pa_table,
                rechunk=rechunk,
                schema=schema,
                schema_overrides=schema_overrides,
            )
        )

    msg = f"expected PyArrow Table, Array, or one or more RecordBatches; got {qualified_type_name(data)!r}"
    raise TypeError(msg)


@overload
def from_pandas(
    data: pd.DataFrame,
    *,
    schema_overrides: SchemaDict | None = ...,
    rechunk: bool = ...,
    nan_to_null: bool = ...,
    include_index: bool = ...,
) -> DataFrame: ...


@overload
def from_pandas(
    data: pd.Series[Any] | pd.Index[Any] | pd.DatetimeIndex,
    *,
    schema_overrides: SchemaDict | None = ...,
    rechunk: bool = ...,
    nan_to_null: bool = ...,
    include_index: Literal[False] = ...,
) -> Series: ...


@overload
def from_pandas(
    data: pd.Series[Any],
    *,
    schema_overrides: SchemaDict | None = ...,
    rechunk: bool = ...,
    nan_to_null: bool = ...,
    include_index: Literal[True],
) -> DataFrame: ...


def from_pandas(
    data: pd.DataFrame | pd.Series[Any] | pd.Index[Any] | pd.DatetimeIndex,
    *,
    schema_overrides: SchemaDict | None = None,
    rechunk: bool = True,
    nan_to_null: bool = True,
    include_index: bool = False,
) -> DataFrame | Series:
    """
    Construct a Polars DataFrame or Series from a pandas DataFrame, Series, or Index.

    This operation may clone data. If you want to ensure that in-place modifications
    of the output don't affect the input, you may want to consider one of the following:

    - Enable `Copy-On-Write <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`_
      in pandas.
    - Call :meth:`DataFrame.clone` on the output of `from_pandas`.

    This requires that :mod:`pandas` and :mod:`pyarrow` are installed.

    Parameters
    ----------
    data : :class:`pandas.DataFrame` or :class:`pandas.Series` or :class:`pandas.Index`
        Data represented as a pandas DataFrame, Series, or Index.
    schema_overrides : dict, default None
        Support override of inferred types for one or more columns.
    rechunk : bool, default True
        Make sure that all data is in contiguous memory.
    nan_to_null : bool, default True
        If data contains `NaN` values PyArrow will convert the `NaN` to `None`
    include_index : bool, default False
        Load any non-default pandas indexes as columns.

        .. note::
            If the input is a pandas ``DataFrame`` and has a nameless index
            which just enumerates the rows, then it will not be included in the
            result, regardless of this parameter. If you want to be sure to include it,
            please call ``.reset_index()`` prior to calling this function.

    Returns
    -------
    DataFrame

    Examples
    --------
    Constructing a :class:`DataFrame` from a :class:`pandas.DataFrame`:

    >>> import pandas as pd
    >>> pd_df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["a", "b", "c"])
    >>> df = pl.from_pandas(pd_df)
    >>> df
        shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 2   ┆ 3   │
    │ 4   ┆ 5   ┆ 6   │
    └─────┴─────┴─────┘

    Constructing a Series from a :class:`pandas.Series`:

    >>> import pandas as pd
    >>> pd_series = pd.Series([1, 2, 3], name="pd")
    >>> df = pl.from_pandas(pd_series)
    >>> df
    shape: (3,)
    Series: 'pd' [i64]
    [
        1
        2
        3
    ]
    """
    if include_index and isinstance(data, pd.Series):
        data = data.reset_index()

    if isinstance(data, (pd.Series, pd.Index, pd.DatetimeIndex)):
        return wrap_s(pandas_to_pyseries("", data, nan_to_null=nan_to_null))
    elif isinstance(data, pd.DataFrame):
        return wrap_df(
            pandas_to_pydf(
                data,
                schema_overrides=schema_overrides,
                rechunk=rechunk,
                nan_to_null=nan_to_null,
                include_index=include_index,
            )
        )
    else:
        msg = f"expected pandas DataFrame or Series, got {qualified_type_name(data)!r}"
        raise TypeError(msg)


@deprecate_renamed_parameter("tbl", "data", version="0.20.17")
def from_repr(data: str) -> DataFrame | Series:
    """
    Construct a Polars DataFrame or Series from its string representation.

    .. versionchanged:: 0.20.17
        The `tbl` parameter was renamed to `data`.

    Parameters
    ----------
    data
        A string containing a polars DataFrame or Series repr; does not need
        to be trimmed of whitespace (or leading prompts) as the repr will be
        found/extracted automatically.

    Notes
    -----
    This function handles the default UTF8_FULL (and UTF8_FULL_CONDENSED) DataFrame
    tables, with or without rounded corners. Truncated columns/rows are omitted,
    wrapped headers are accounted for, and dtypes are automatically identified.

    Currently compound/nested dtypes such as List and Struct are not supported;
    neither are Object dtypes. The DuckDB table/relation repr is also compatible
    with this function.

    See Also
    --------
    polars.DataFrame.to_init_repr
    polars.Series.to_init_repr

    Examples
    --------
    From DataFrame table repr:

    >>> df = pl.from_repr(
    ...     '''
    ...     Out[3]:
    ...     shape: (1, 5)
    ...     ┌───────────┬────────────┬───┬───────┬────────────────────────────────┐
    ...     │ source_ac ┆ source_cha ┆ … ┆ ident ┆ timestamp                      │
    ...     │ tor_id    ┆ nnel_id    ┆   ┆ ---   ┆ ---                            │
    ...     │ ---       ┆ ---        ┆   ┆ str   ┆ datetime[μs, Asia/Tokyo]       │
    ...     │ i32       ┆ i64        ┆   ┆       ┆                                │
    ...     ╞═══════════╪════════════╪═══╪═══════╪════════════════════════════════╡
    ...     │ 123456780 ┆ 9876543210 ┆ … ┆ a:b:c ┆ 2023-03-25 10:56:59.663053 JST │
    ...     │ …         ┆ …          ┆ … ┆ …     ┆ …                              │
    ...     │ 803065983 ┆ 2055938745 ┆ … ┆ x:y:z ┆ 2023-03-25 12:38:18.050545 JST │
    ...     └───────────┴────────────┴───┴───────┴────────────────────────────────┘
    ... '''
    ... )
    >>> df
    shape: (2, 4)
    ┌─────────────────┬───────────────────┬───────┬────────────────────────────────┐
    │ source_actor_id ┆ source_channel_id ┆ ident ┆ timestamp                      │
    │ ---             ┆ ---               ┆ ---   ┆ ---                            │
    │ i32             ┆ i64               ┆ str   ┆ datetime[μs, Asia/Tokyo]       │
    ╞═════════════════╪═══════════════════╪═══════╪════════════════════════════════╡
    │ 123456780       ┆ 9876543210        ┆ a:b:c ┆ 2023-03-25 10:56:59.663053 JST │
    │ 803065983       ┆ 2055938745        ┆ x:y:z ┆ 2023-03-25 12:38:18.050545 JST │
    └─────────────────┴───────────────────┴───────┴────────────────────────────────┘

    From Series repr:

    >>> s = pl.from_repr(
    ...     '''
    ...     shape: (3,)
    ...     Series: 's' [bool]
    ...     [
    ...        true
    ...        false
    ...        true
    ...     ]
    ...     '''
    ... )
    >>> s.to_list()
    [True, False, True]
    """
    # find DataFrame table...
    m = re.search(r"([┌╭].*?[┘╯])", data, re.DOTALL)
    if m is not None:
        return _from_dataframe_repr(m)

    # ...or Series in the given string
    m = re.search(
        pattern=r"(?:shape: (\(\d+,\))\n.*?)?Series:\s+([^\n]+)\s+\[([^\n]+)](.*)",
        string=data,
        flags=re.DOTALL,
    )
    if m is not None:
        return _from_series_repr(m)

    msg = "input string does not contain DataFrame or Series"
    raise ValueError(msg)


def _from_dataframe_repr(m: re.Match[str]) -> DataFrame:
    """Reconstruct a DataFrame from a regex-matched table repr."""
    from polars.datatypes.convert import dtype_short_repr_to_dtype
    from polars.io.database._inference import dtype_from_database_typename

    def _dtype_from_name(tp: str | None) -> PolarsDataType | None:
        return (
            None
            if tp is None
            else (
                dtype_short_repr_to_dtype(tp)
                or dtype_from_database_typename(tp, raise_unmatched=False)
            )
        )

    # extract elements from table structure
    lines = m.group().split("\n")[1:-1]
    rows = [
        [re.sub(r"^[\W+]*│", "", elem).strip() for elem in row]
        for row in [re.split("[│┆|]", row.lstrip("#. ").rstrip("│ ")) for row in lines]
        if len(row) > 1 or not re.search("├[╌┼]+┤", row[0])
    ]

    # determine beginning/end of the header block
    table_body_start = 2
    found_header_divider = False
    for idx, (elem, *_) in enumerate(rows):
        if re.match(r"^\W*[╞]", elem):
            found_header_divider = True
            table_body_start = idx
            break

    # handle headers with wrapped column names and determine headers/dtypes
    header_rows = rows[:table_body_start]
    header_block: list[Sequence[str]]
    if (
        not found_header_divider
        and len(header_rows) == 2
        and not any("---" in h for h in header_rows)
    ):
        header_block = list(zip(*header_rows))
    else:
        header_block = ["".join(h).split("---") for h in zip(*header_rows)]

    dtypes: list[str | None]
    if all(len(h) == 1 for h in header_block):
        headers = [h[0] for h in header_block]
        dtypes = [None] * len(headers)
    else:
        headers, dtypes = (list(h) for h in itertools.zip_longest(*header_block))

    body = rows[table_body_start + 1 :]
    if not headers[0] and not dtypes[0]:
        body = [row[1:] for row in body]
        headers = headers[1:]
        dtypes = dtypes[1:]

    no_dtypes = all(d is None for d in dtypes)

    # transpose rows into columns, detect/omit truncated columns
    coldata = list(zip(*(row for row in body if not all((e == "…") for e in row))))
    for el in ("…", "..."):
        if el in headers:
            idx = headers.index(el)
            for table_elem in (headers, dtypes):
                table_elem.pop(idx)
            if coldata:
                coldata.pop(idx)

    # init cols as String Series, handle "null" -> None, create schema from repr dtype
    data = [
        pl.Series([(None if v in ("null", "NULL") else v) for v in cd], dtype=String)
        for cd in coldata
    ]
    schema = dict(zip(headers, (_dtype_from_name(d) for d in dtypes)))
    if schema and data and (n_extend_cols := (len(schema) - len(data))) > 0:
        empty_data = [None] * len(data[0])
        data.extend((pl.Series(empty_data, dtype=String)) for _ in range(n_extend_cols))

    for dtype in set(schema.values()):
        if dtype is not None and (dtype.is_nested() or dtype.is_object()):
            msg = (
                f"`from_repr` does not support data type {dtype.base_type().__name__!r}"
            )
            raise NotImplementedError(msg)

    # construct DataFrame from string series and cast from repr to native dtype
    df = pl.DataFrame(data=data, orient="col", schema=list(schema))
    if no_dtypes:
        if df.is_empty():
            # if no dtypes *and* empty, default to string
            return df.with_columns(F.all().cast(String))
        else:
            # otherwise, take a trip through our CSV inference logic
            if all(tp == String for tp in df.schema.values()):
                from polars.io import read_csv

                buf = io.BytesIO()
                df.write_csv(file=buf)
                buf.seek(0)
                df = read_csv(
                    buf,
                    new_columns=df.columns,
                    try_parse_dates=True,
                    infer_schema_length=None,
                )
            return df
    elif schema and not data:
        return df.cast(schema)  # type: ignore[arg-type]
    else:
        return _cast_repr_strings_with_schema(df, schema)


def _from_series_repr(m: re.Match[str]) -> Series:
    """Reconstruct a Series from a regex-matched series repr."""
    from polars.datatypes.convert import dtype_short_repr_to_dtype

    shape = m.groups()[0]
    name = m.groups()[1][1:-1]
    length = int(shape[1:-2] if shape else -1)
    dtype = dtype_short_repr_to_dtype(m.groups()[2])

    if length == 0:
        string_values = []
    else:
        string_values = [
            v.strip()
            for v in re.findall(r"[\s>#]*(?:\t|\s{2,})([^\n]*)\n", m.groups()[-1])
        ]
        if string_values == ["[", "]"]:
            string_values = []
        else:
            start: int | None = None
            end: int | None = None
            for idx, v in enumerate(string_values):
                if start is None and v.lstrip("#> ") == "[":
                    start = idx
                if v.lstrip("#> ") == "]":
                    end = idx
            if start is not None and end is not None:
                string_values = string_values[start + 1 : end]

    values = string_values[:length] if length > 0 else string_values
    values = [(None if v == "null" else v) for v in values if v not in ("…", "...")]

    if not values:
        return pl.Series(name=name, values=values, dtype=dtype)
    else:
        srs = pl.Series(name=name, values=values, dtype=String)
        if dtype is None:
            return srs
        elif dtype in (Categorical, String):
            return srs.str.replace('^"(.*)"$', r"$1").cast(dtype)

        return _cast_repr_strings_with_schema(
            srs.to_frame(), schema={srs.name: dtype}
        ).to_series()


def from_dataframe(
    df: SupportsInterchange | ArrowArrayExportable | ArrowStreamExportable,
    *,
    allow_copy: bool | None = None,
    rechunk: bool = True,
) -> DataFrame:
    """
    Build a Polars DataFrame from any dataframe supporting the PyCapsule Interface.

    .. versionchanged:: 1.23.0

       `from_dataframe` uses the PyCapsule Interface instead of the Dataframe
       Interchange Protocol for conversion, only using the latter as a fallback.

    Parameters
    ----------
    df
        Object supporting the dataframe PyCapsule Interface.
    allow_copy
        Allow memory to be copied to perform the conversion. If set to False, may cause
        conversions that are not zero-copy to fail.

        .. deprecated: 1.23.0
            `allow_copy` is deprecated and will be removed in a future version.
    rechunk : bool, default True
        Make sure that all data is in contiguous memory.

    Notes
    -----
    - Details on the PyCapsule Interface:
      https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html.
    - Details on the Python dataframe interchange protocol:
      https://data-apis.org/dataframe-protocol/latest/index.html.
      Using a dedicated function like :func:`from_pandas` or :func:`from_arrow` is
      a more efficient method of conversion.

    Examples
    --------
    Convert a pandas dataframe to Polars.

    >>> import pandas as pd
    >>> df_pd = pd.DataFrame({"a": [1, 2], "b": [3.0, 4.0], "c": ["x", "y"]})
    >>> pl.from_dataframe(df_pd)
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ f64 ┆ str │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3.0 ┆ x   │
    │ 2   ┆ 4.0 ┆ y   │
    └─────┴─────┴─────┘
    """
    if allow_copy is not None:
        issue_deprecation_warning(
            "`allow_copy` is deprecated and will be removed in a future version.",
            version="1.23",
        )
    else:
        allow_copy = True
    if is_pycapsule(df):
        try:
            return pycapsule_to_frame(df, rechunk=rechunk)
        except Exception as exc:
            issue_warning(
                f"Failed to convert dataframe using PyCapsule Interface with exception: {exc!r}.\n"
                "Falling back to Dataframe Interchange Protocol, which is known to be less robust.",
                UserWarning,
            )
    from polars.interchange.from_dataframe import from_dataframe

    result = from_dataframe(df, allow_copy=allow_copy)  # type: ignore[arg-type]
    if rechunk:
        return result.rechunk()
    return result
