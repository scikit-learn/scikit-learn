"""Module containing logic related to eager DataFrames."""

from __future__ import annotations

import contextlib
import os
import random
from collections import defaultdict
from collections.abc import (
    Generator,
    Iterable,
    Mapping,
    Sequence,
    Sized,
)
from io import BytesIO, StringIO
from pathlib import Path
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    NoReturn,
    TypeVar,
    cast,
    get_args,
    overload,
)

import polars._reexport as pl
from polars import functions as F
from polars._dependencies import (
    _ALTAIR_AVAILABLE,
    _GREAT_TABLES_AVAILABLE,
    _PANDAS_AVAILABLE,
    _PYARROW_AVAILABLE,
    _check_for_numpy,
    _check_for_pandas,
    _check_for_pyarrow,
    _check_for_torch,
    altair,
    great_tables,
    import_optional,
    torch,
)
from polars._dependencies import numpy as np
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._typing import DbWriteMode, JaxExportType, TorchExportType
from polars._utils.construction import (
    arrow_to_pydf,
    dataframe_to_pydf,
    dict_to_pydf,
    iterable_to_pydf,
    numpy_to_pydf,
    pandas_to_pydf,
    sequence_to_pydf,
    series_to_pydf,
)
from polars._utils.convert import parse_as_duration_string
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    deprecated,
    issue_deprecation_warning,
)
from polars._utils.getitem import get_df_item_by_key
from polars._utils.parse import parse_into_expression
from polars._utils.pycapsule import is_pycapsule, pycapsule_to_frame
from polars._utils.serde import serialize_polars_object
from polars._utils.unstable import issue_unstable_warning, unstable
from polars._utils.various import (
    is_bool_sequence,
    no_default,
    normalize_filepath,
    parse_version,
    qualified_type_name,
    require_same_type,
    scale_bytes,
    warn_null_comparison,
)
from polars._utils.wrap import wrap_expr, wrap_ldf, wrap_s
from polars.dataframe._html import NotebookFormatter
from polars.dataframe.group_by import DynamicGroupBy, GroupBy, RollingGroupBy
from polars.dataframe.plotting import DataFramePlot
from polars.datatypes import (
    N_INFER_DEFAULT,
    Boolean,
    Float32,
    Float64,
    Int32,
    Int64,
    Null,
    Object,
    String,
    Struct,
    UInt16,
    UInt32,
    UInt64,
)
from polars.datatypes.group import INTEGER_DTYPES
from polars.exceptions import (
    ColumnNotFoundError,
    InvalidOperationError,
    ModuleUpgradeRequiredError,
    NoRowsReturnedError,
    TooManyRowsReturnedError,
)
from polars.functions import col, lit
from polars.interchange.protocol import CompatLevel
from polars.schema import Schema
from polars.selectors import _expand_selector_dicts, _expand_selectors

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyDataFrame
    from polars._plr import dtype_str_repr as _dtype_str_repr
    from polars._plr import write_clipboard_string as _write_clipboard_string

if TYPE_CHECKING:
    import sys
    from collections.abc import Collection, Iterator, Mapping
    from datetime import timedelta
    from io import IOBase
    from typing import Literal

    import deltalake
    import jax
    import numpy.typing as npt
    import pyiceberg
    from great_tables import GT
    from xlsxwriter import Workbook
    from xlsxwriter.worksheet import Worksheet

    from polars import DataType, Expr, LazyFrame, Series
    from polars._typing import (
        AsofJoinStrategy,
        AvroCompression,
        ClosedInterval,
        ColumnFormatDict,
        ColumnNameOrSelector,
        ColumnTotalsDefinition,
        ColumnWidthsDefinition,
        ComparisonOperator,
        ConditionalFormatDict,
        ConnectionOrCursor,
        CsvQuoteStyle,
        DbWriteEngine,
        EngineType,
        FillNullStrategy,
        FrameInitTypes,
        IndexOrder,
        IntoExpr,
        IntoExprColumn,
        IpcCompression,
        JoinStrategy,
        JoinValidation,
        Label,
        MaintainOrderJoin,
        MultiColSelector,
        MultiIndexSelector,
        OneOrMoreDataTypes,
        Orientation,
        ParquetCompression,
        ParquetMetadata,
        PartitioningScheme,
        PivotAgg,
        PolarsDataType,
        PythonDataType,
        QuantileMethod,
        RowTotalsDefinition,
        SchemaDefinition,
        SchemaDict,
        SelectorType,
        SerializationFormat,
        SingleColSelector,
        SingleIndexSelector,
        SizeUnit,
        StartBy,
        UniqueKeepStrategy,
        UnstackDirection,
    )
    from polars._utils.various import NoDefault
    from polars.interchange.dataframe import PolarsDataFrame
    from polars.io.cloud import CredentialProviderFunction
    from polars.ml.torch import PolarsDataset

    if sys.version_info >= (3, 10):
        from typing import Concatenate, ParamSpec
    else:
        from typing_extensions import Concatenate, ParamSpec

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004

    T = TypeVar("T")
    P = ParamSpec("P")


class DataFrame:
    """
    Two-dimensional data structure representing data as a table with rows and columns.

    Parameters
    ----------
    data : dict, Sequence, ndarray, Series, or pandas.DataFrame
        Two-dimensional data in various forms; dict input must contain Sequences,
        Generators, or a `range`. Sequence may contain Series or other Sequences.
    schema : Sequence of str, (str,DataType) pairs, or a {str:DataType,} dict
        The schema of the resulting DataFrame. The schema may be declared in several
        ways:

        * As a dict of {name:type} pairs; if type is None, it will be auto-inferred.
        * As a list of column names; in this case types are automatically inferred.
        * As a list of (name,type) pairs; this is equivalent to the dictionary form.

        If you supply a list of column names that does not match the names in the
        underlying data, the names given here will overwrite them. The number
        of names given in the schema should match the underlying data dimensions.

        If set to `None` (default), the schema is inferred from the data.
    schema_overrides : dict, default None
        Support type specification or override of one or more columns; note that
        any dtypes inferred from the schema param will be overridden.

        The number of entries in the schema should match the underlying data
        dimensions, unless a sequence of dictionaries is being passed, in which case
        a *partial* schema can be declared to prevent specific fields from being loaded.
    strict : bool, default True
        Throw an error if any `data` value does not exactly match the given or inferred
        data type for that column. If set to `False`, values that do not match the data
        type are cast to that data type or, if casting is not possible, set to null
        instead.
    orient : {'col', 'row'}, default None
        Whether to interpret two-dimensional data as columns or as rows. If None,
        the orientation is inferred by matching the columns and data dimensions. If
        this does not yield conclusive results, column orientation is used.
    infer_schema_length : int or None
        The maximum number of rows to scan for schema inference. If set to `None`, the
        full data may be scanned *(this can be slow)*. This parameter only applies if
        the input data is a sequence or generator of rows; other input is read as-is.
    nan_to_null : bool, default False
        If the data comes from one or more numpy arrays, can optionally convert input
        data np.nan values to null instead. This is a no-op for all other input data.

    Notes
    -----
    Polars explicitly does not support subclassing of its core data types. See
    the following GitHub issue for possible workarounds:
    https://github.com/pola-rs/polars/issues/2846#issuecomment-1711799869

    Examples
    --------
    Constructing a DataFrame from a dictionary:

    >>> data = {"a": [1, 2], "b": [3, 4]}
    >>> df = pl.DataFrame(data)
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

    Notice that the dtypes are automatically inferred as polars Int64:

    >>> df.dtypes
    [Int64, Int64]

    To specify a more detailed/specific frame schema you can supply the `schema`
    parameter with a dictionary of (name,dtype) pairs...

    >>> data = {"col1": [0, 2], "col2": [3, 7]}
    >>> df2 = pl.DataFrame(data, schema={"col1": pl.Float32, "col2": pl.Int64})
    >>> df2
    shape: (2, 2)
    ┌──────┬──────┐
    │ col1 ┆ col2 │
    │ ---  ┆ ---  │
    │ f32  ┆ i64  │
    ╞══════╪══════╡
    │ 0.0  ┆ 3    │
    │ 2.0  ┆ 7    │
    └──────┴──────┘

    ...a sequence of (name,dtype) pairs...

    >>> data = {"col1": [1, 2], "col2": [3, 4]}
    >>> df3 = pl.DataFrame(data, schema=[("col1", pl.Float32), ("col2", pl.Int64)])
    >>> df3
    shape: (2, 2)
    ┌──────┬──────┐
    │ col1 ┆ col2 │
    │ ---  ┆ ---  │
    │ f32  ┆ i64  │
    ╞══════╪══════╡
    │ 1.0  ┆ 3    │
    │ 2.0  ┆ 4    │
    └──────┴──────┘

    ...or a list of typed Series.

    >>> data = [
    ...     pl.Series("col1", [1, 2], dtype=pl.Float32),
    ...     pl.Series("col2", [3, 4], dtype=pl.Int64),
    ... ]
    >>> df4 = pl.DataFrame(data)
    >>> df4
    shape: (2, 2)
    ┌──────┬──────┐
    │ col1 ┆ col2 │
    │ ---  ┆ ---  │
    │ f32  ┆ i64  │
    ╞══════╪══════╡
    │ 1.0  ┆ 3    │
    │ 2.0  ┆ 4    │
    └──────┴──────┘

    Constructing a DataFrame from a numpy ndarray, specifying column names:

    >>> import numpy as np
    >>> data = np.array([(1, 2), (3, 4)], dtype=np.int64)
    >>> df5 = pl.DataFrame(data, schema=["a", "b"], orient="col")
    >>> df5
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    │ 2   ┆ 4   │
    └─────┴─────┘

    Constructing a DataFrame from a list of lists, row orientation specified:

    >>> data = [[1, 2, 3], [4, 5, 6]]
    >>> df6 = pl.DataFrame(data, schema=["a", "b", "c"], orient="row")
    >>> df6
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 2   ┆ 3   │
    │ 4   ┆ 5   ┆ 6   │
    └─────┴─────┴─────┘
    """

    _df: PyDataFrame
    _accessors: ClassVar[set[str]] = {"plot", "style"}

    def __init__(
        self,
        data: FrameInitTypes | None = None,
        schema: SchemaDefinition | None = None,
        *,
        schema_overrides: SchemaDict | None = None,
        strict: bool = True,
        orient: Orientation | None = None,
        infer_schema_length: int | None = N_INFER_DEFAULT,
        nan_to_null: bool = False,
    ) -> None:
        if data is None:
            self._df = dict_to_pydf(
                {}, schema=schema, schema_overrides=schema_overrides
            )

        elif isinstance(data, dict):
            self._df = dict_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                nan_to_null=nan_to_null,
            )

        elif isinstance(data, (list, tuple, Sequence)):
            self._df = sequence_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                orient=orient,
                infer_schema_length=infer_schema_length,
                nan_to_null=nan_to_null,
            )

        elif isinstance(data, pl.Series):
            self._df = series_to_pydf(
                data, schema=schema, schema_overrides=schema_overrides, strict=strict
            )

        elif _check_for_numpy(data) and isinstance(data, np.ndarray):
            self._df = numpy_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                orient=orient,
                nan_to_null=nan_to_null,
            )

        elif _check_for_pyarrow(data) and isinstance(data, pa.Table):
            self._df = arrow_to_pydf(
                data, schema=schema, schema_overrides=schema_overrides, strict=strict
            )

        elif _check_for_pandas(data) and isinstance(data, pd.DataFrame):
            self._df = pandas_to_pydf(
                data, schema=schema, schema_overrides=schema_overrides, strict=strict
            )

        elif _check_for_torch(data) and isinstance(data, torch.Tensor):
            self._df = numpy_to_pydf(
                data.numpy(force=False),
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                orient=orient,
                nan_to_null=nan_to_null,
            )

        elif (
            not hasattr(data, "__arrow_c_stream__")
            and not isinstance(data, Sized)
            and isinstance(data, (Generator, Iterable))
        ):
            self._df = iterable_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                orient=orient,
                infer_schema_length=infer_schema_length,
            )

        elif isinstance(data, pl.DataFrame):
            self._df = dataframe_to_pydf(
                data, schema=schema, schema_overrides=schema_overrides, strict=strict
            )

        elif is_pycapsule(data):
            self._df = pycapsule_to_frame(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
            )._df
        else:
            msg = (
                f"DataFrame constructor called with unsupported type {type(data).__name__!r}"
                " for the `data` parameter"
            )
            raise TypeError(msg)

    @classmethod
    def deserialize(
        cls, source: str | Path | IOBase, *, format: SerializationFormat = "binary"
    ) -> DataFrame:
        """
        Read a serialized DataFrame from a file.

        Parameters
        ----------
        source
            Path to a file or a file-like object (by file-like object, we refer to
            objects that have a `read()` method, such as a file handler (e.g.
            via builtin `open` function) or `BytesIO`).
        format
            The format with which the DataFrame was serialized. Options:

            - `"binary"`: Deserialize from binary format (bytes). This is the default.
            - `"json"`: Deserialize from JSON format (string).

        See Also
        --------
        DataFrame.serialize

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        >>> import io
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})
        >>> bytes = df.serialize()
        >>> pl.DataFrame.deserialize(io.BytesIO(bytes))
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ f64 │
        ╞═════╪═════╡
        │ 1   ┆ 4.0 │
        │ 2   ┆ 5.0 │
        │ 3   ┆ 6.0 │
        └─────┴─────┘
        """
        if isinstance(source, StringIO):
            source = BytesIO(source.getvalue().encode())
        elif isinstance(source, (str, Path)):
            source = normalize_filepath(source)

        if format == "binary":
            deserializer = PyDataFrame.deserialize_binary
        elif format == "json":
            deserializer = PyDataFrame.deserialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return cls._from_pydf(deserializer(source))

    @classmethod
    def _from_pydf(cls, py_df: PyDataFrame) -> DataFrame:
        """Construct Polars DataFrame from FFI PyDataFrame object."""
        df = cls.__new__(cls)
        df._df = py_df
        return df

    @classmethod
    def _from_arrow(
        cls,
        data: pa.Table | pa.RecordBatch,
        schema: SchemaDefinition | None = None,
        *,
        schema_overrides: SchemaDict | None = None,
        rechunk: bool = True,
    ) -> DataFrame:
        """
        Construct a DataFrame from an Arrow table.

        This operation will be zero copy for the most part. Types that are not
        supported by Polars may be cast to the closest supported type.

        Parameters
        ----------
        data : arrow Table, RecordBatch, or sequence of sequences
            Data representing an Arrow Table or RecordBatch.
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
        rechunk : bool, default True
            Make sure that all data is in contiguous memory.
        """
        return cls._from_pydf(
            arrow_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                rechunk=rechunk,
            )
        )

    @classmethod
    def _from_pandas(
        cls,
        data: pd.DataFrame,
        schema: SchemaDefinition | None = None,
        *,
        schema_overrides: SchemaDict | None = None,
        rechunk: bool = True,
        nan_to_null: bool = True,
        include_index: bool = False,
    ) -> DataFrame:
        """
        Construct a Polars DataFrame from a pandas DataFrame.

        Parameters
        ----------
        data : pandas DataFrame
            Two-dimensional data represented as a pandas DataFrame.
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
        rechunk : bool, default True
            Make sure that all data is in contiguous memory.
        nan_to_null : bool, default True
            If the data contains NaN values they will be converted to null/None.
        include_index : bool, default False
            Load any non-default pandas indexes as columns.
        """
        return cls._from_pydf(
            pandas_to_pydf(
                data,
                schema=schema,
                schema_overrides=schema_overrides,
                rechunk=rechunk,
                nan_to_null=nan_to_null,
                include_index=include_index,
            )
        )

    def _replace(self, column: str, new_column: Series) -> DataFrame:
        """Replace a column by a new Series (in place)."""
        self._df.replace(column, new_column._s)
        return self

    @classmethod
    def _import_columns(cls, pointer: int, width: int) -> DataFrame:
        return cls._from_pydf(PyDataFrame._import_columns(pointer, width))

    @property
    @unstable()
    def plot(self) -> DataFramePlot:
        """
        Create a plot namespace.

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        .. versionchanged:: 1.6.0
            In prior versions of Polars, HvPlot was the plotting backend. If you would
            like to restore the previous plotting functionality, all you need to do
            is add `import hvplot.polars` at the top of your script and replace
            `df.plot` with `df.hvplot`.

        Polars does not implement plotting logic itself, but instead defers to
        `Altair <https://altair-viz.github.io/>`_:

        - `df.plot.line(**kwargs)`
          is shorthand for
          `alt.Chart(df).mark_line(tooltip=True).encode(**kwargs).interactive()`
        - `df.plot.point(**kwargs)`
          is shorthand for
          `alt.Chart(df).mark_point(tooltip=True).encode(**kwargs).interactive()` (and
          `plot.scatter` is provided as an alias)
        - `df.plot.bar(**kwargs)`
          is shorthand for
          `alt.Chart(df).mark_bar(tooltip=True).encode(**kwargs).interactive()`
        - for any other attribute `attr`, `df.plot.attr(**kwargs)`
          is shorthand for
          `alt.Chart(df).mark_attr(tooltip=True).encode(**kwargs).interactive()`

        For configuration, we suggest reading
        `Chart Configuration <https://altair-viz.github.io/altair-tutorial/notebooks/08-Configuration.html>`_.
        For example, you can:

        - Change the width/height/title with
          ``.properties(width=500, height=350, title="My amazing plot")``.
        - Change the x-axis label rotation with ``.configure_axisX(labelAngle=30)``.
        - Change the opacity of the points in your scatter plot with
          ``.configure_point(opacity=.5)``.

        Examples
        --------
        Scatter plot:

        >>> df = pl.DataFrame(
        ...     {
        ...         "length": [1, 4, 6],
        ...         "width": [4, 5, 6],
        ...         "species": ["setosa", "setosa", "versicolor"],
        ...     }
        ... )
        >>> df.plot.point(x="length", y="width", color="species")  # doctest: +SKIP

        Set the x-axis title by using ``altair.X``:

        >>> import altair as alt
        >>> df.plot.point(
        ...     x=alt.X("length", title="Length"), y="width", color="species"
        ... )  # doctest: +SKIP

        Line plot:

        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "date": [date(2020, 1, 2), date(2020, 1, 3), date(2020, 1, 4)] * 2,
        ...         "price": [1, 4, 6, 1, 5, 2],
        ...         "stock": ["a", "a", "a", "b", "b", "b"],
        ...     }
        ... )
        >>> df.plot.line(x="date", y="price", color="stock")  # doctest: +SKIP

        Bar plot:

        >>> df = pl.DataFrame(
        ...     {
        ...         "day": ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"] * 2,
        ...         "group": ["a"] * 7 + ["b"] * 7,
        ...         "value": [1, 3, 2, 4, 5, 6, 1, 1, 3, 2, 4, 5, 1, 2],
        ...     }
        ... )
        >>> df.plot.bar(
        ...     x="day", y="value", color="day", column="group"
        ... )  # doctest: +SKIP

        Or, to make a stacked version of the plot above:

        >>> df.plot.bar(x="day", y="value", color="group")  # doctest: +SKIP
        """
        if not _ALTAIR_AVAILABLE or parse_version(altair.__version__) < (5, 4, 0):
            msg = "altair>=5.4.0 is required for `.plot`"
            raise ModuleUpgradeRequiredError(msg)
        return DataFramePlot(self)

    @property
    @unstable()
    def style(self) -> GT:
        """
        Create a Great Table for styling.

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Polars does not implement styling logic itself, but instead defers to
        the Great Tables package. Please see the `Great Tables reference <https://posit-dev.github.io/great-tables/reference/>`_
        for more information and documentation.

        Examples
        --------
        Import some styling helpers, and create example data:

        >>> import polars.selectors as cs
        >>> from great_tables import loc, style
        >>> df = pl.DataFrame(
        ...     {
        ...         "site_id": [0, 1, 2],
        ...         "measure_a": [5, 4, 6],
        ...         "measure_b": [7, 3, 3],
        ...     }
        ... )

        Emphasize the site_id as row names:

        >>> df.style.tab_stub(rowname_col="site_id")  # doctest: +SKIP

        Fill the background for the highest measure_a value row:

        >>> df.style.tab_style(
        ...     style.fill("yellow"),
        ...     loc.body(rows=pl.col("measure_a") == pl.col("measure_a").max()),
        ... )  # doctest: +SKIP

        Put a spanner (high-level label) over measure columns:

        >>> df.style.tab_spanner(
        ...     "Measures", cs.starts_with("measure")
        ... )  # doctest: +SKIP

        Format measure_b values to two decimal places:

        >>> df.style.fmt_number("measure_b", decimals=2)  # doctest: +SKIP
        """
        if not _GREAT_TABLES_AVAILABLE:
            msg = "great_tables is required for `.style`"
            raise ModuleNotFoundError(msg)

        return great_tables.GT(self)

    @property
    def shape(self) -> tuple[int, int]:
        """
        Get the shape of the DataFrame.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5]})
        >>> df.shape
        (5, 1)
        """
        return self._df.shape()

    @property
    def height(self) -> int:
        """
        Get the number of rows.

        Returns
        -------
        int

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3, 4, 5]})
        >>> df.height
        5
        """
        return self._df.height()

    @property
    def width(self) -> int:
        """
        Get the number of columns.

        Returns
        -------
        int

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4, 5, 6],
        ...     }
        ... )
        >>> df.width
        2
        """
        return self._df.width()

    @property
    def columns(self) -> list[str]:
        """
        Get or set column names.

        Returns
        -------
        list of str
            A list containing the name of each column in order.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.columns
        ['foo', 'bar', 'ham']

        Set column names:

        >>> df.columns = ["apple", "banana", "orange"]
        >>> df
        shape: (3, 3)
        ┌───────┬────────┬────────┐
        │ apple ┆ banana ┆ orange │
        │ ---   ┆ ---    ┆ ---    │
        │ i64   ┆ i64    ┆ str    │
        ╞═══════╪════════╪════════╡
        │ 1     ┆ 6      ┆ a      │
        │ 2     ┆ 7      ┆ b      │
        │ 3     ┆ 8      ┆ c      │
        └───────┴────────┴────────┘
        """
        return self._df.columns()

    @columns.setter
    def columns(self, names: Sequence[str]) -> None:
        """
        Change the column names of the `DataFrame`.

        Parameters
        ----------
        names
            A list with new names for the `DataFrame`.
            The length of the list should be equal to the width of the `DataFrame`.
        """
        self._df.set_column_names(names)

    @property
    def dtypes(self) -> list[DataType]:
        """
        Get the column data types.

        The data types can also be found in column headers when printing the DataFrame.

        Returns
        -------
        list of DataType
            A list containing the data type of each column in order.

        See Also
        --------
        schema

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.dtypes
        [Int64, Float64, String]
        >>> df
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6.0 ┆ a   │
        │ 2   ┆ 7.0 ┆ b   │
        │ 3   ┆ 8.0 ┆ c   │
        └─────┴─────┴─────┘
        """
        return self._df.dtypes()

    @property
    def flags(self) -> dict[str, dict[str, bool]]:
        """
        Get flags that are set on the columns of this DataFrame.

        Returns
        -------
        dict
            Mapping from column names to column flags.
        """
        return {name: self[name].flags for name in self.columns}

    @property
    def schema(self) -> Schema:
        """
        Get an ordered mapping of column names to their data type.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.schema
        Schema({'foo': Int64, 'bar': Float64, 'ham': String})
        """
        return Schema(zip(self.columns, self.dtypes), check_dtypes=False)

    def __array__(
        self, dtype: npt.DTypeLike | None = None, copy: bool | None = None
    ) -> np.ndarray[Any, Any]:
        """
        Return a NumPy ndarray with the given data type.

        This method ensures a Polars DataFrame can be treated as a NumPy ndarray.
        It enables `np.asarray` and NumPy universal functions.

        See the NumPy documentation for more information:
        https://numpy.org/doc/stable/user/basics.interoperability.html#the-array-method
        """
        if copy is None:
            writable, allow_copy = False, True
        elif copy is True:
            writable, allow_copy = True, True
        elif copy is False:
            writable, allow_copy = False, False
        else:
            msg = f"invalid input for `copy`: {copy!r}"
            raise TypeError(msg)

        arr = self.to_numpy(writable=writable, allow_copy=allow_copy)

        if dtype is not None and dtype != arr.dtype:
            if copy is False:
                # TODO: Only raise when data must be copied
                msg = f"copy not allowed: cast from {arr.dtype} to {dtype} prohibited"
                raise RuntimeError(msg)

            arr = arr.__array__(dtype)

        return arr

    def __dataframe__(
        self,
        nan_as_null: bool = False,  # noqa: FBT001
        allow_copy: bool = True,  # noqa: FBT001
    ) -> PolarsDataFrame:
        """
        Convert to a dataframe object implementing the dataframe interchange protocol.

        Parameters
        ----------
        nan_as_null
            Overwrite null values in the data with `NaN`.

            .. warning::
                This functionality has not been implemented and the parameter will be
                removed in a future version.
                Setting this to `True` will raise a `NotImplementedError`.
        allow_copy
            Allow memory to be copied to perform the conversion. If set to `False`,
            causes conversions that are not zero-copy to fail.

        Notes
        -----
        Details on the Python dataframe interchange protocol:
        https://data-apis.org/dataframe-protocol/latest/index.html

        Examples
        --------
        Convert a Polars DataFrame to a generic dataframe object and access some
        properties.

        >>> df = pl.DataFrame({"a": [1, 2], "b": [3.0, 4.0], "c": ["x", "y"]})
        >>> dfi = df.__dataframe__()
        >>> dfi.num_rows()
        2
        >>> dfi.get_column(1).dtype
        (<DtypeKind.FLOAT: 2>, 64, 'g', '=')
        """
        if nan_as_null:
            msg = (
                "functionality for `nan_as_null` has not been implemented and the"
                " parameter will be removed in a future version"
                "\n\nUse the default `nan_as_null=False`."
            )
            raise NotImplementedError(msg)

        from polars.interchange.dataframe import PolarsDataFrame

        return PolarsDataFrame(self, allow_copy=allow_copy)

    def _comp(self, other: Any, op: ComparisonOperator) -> DataFrame:
        """Compare a DataFrame with another object."""
        if isinstance(other, DataFrame):
            return self._compare_to_other_df(other, op)
        else:
            return self._compare_to_non_df(other, op)

    def _compare_to_other_df(
        self,
        other: DataFrame,
        op: ComparisonOperator,
    ) -> DataFrame:
        """Compare a DataFrame with another DataFrame."""
        if self.columns != other.columns:
            msg = "DataFrame columns do not match"
            raise ValueError(msg)
        if self.shape != other.shape:
            msg = "DataFrame dimensions do not match"
            raise ValueError(msg)

        suffix = "__POLARS_CMP_OTHER"
        other_renamed = other.select(F.all().name.suffix(suffix))
        combined = F.concat([self, other_renamed], how="horizontal")

        if op == "eq":
            expr = [F.col(n) == F.col(f"{n}{suffix}") for n in self.columns]
        elif op == "neq":
            expr = [F.col(n) != F.col(f"{n}{suffix}") for n in self.columns]
        elif op == "gt":
            expr = [F.col(n) > F.col(f"{n}{suffix}") for n in self.columns]
        elif op == "lt":
            expr = [F.col(n) < F.col(f"{n}{suffix}") for n in self.columns]
        elif op == "gt_eq":
            expr = [F.col(n) >= F.col(f"{n}{suffix}") for n in self.columns]
        elif op == "lt_eq":
            expr = [F.col(n) <= F.col(f"{n}{suffix}") for n in self.columns]
        else:
            msg = f"unexpected comparison operator {op!r}"
            raise ValueError(msg)

        return combined.select(expr)

    def _compare_to_non_df(
        self,
        other: Any,
        op: ComparisonOperator,
    ) -> DataFrame:
        """Compare a DataFrame with a non-DataFrame object."""
        warn_null_comparison(other)
        if op == "eq":
            return self.select(F.all() == other)
        elif op == "neq":
            return self.select(F.all() != other)
        elif op == "gt":
            return self.select(F.all() > other)
        elif op == "lt":
            return self.select(F.all() < other)
        elif op == "gt_eq":
            return self.select(F.all() >= other)
        elif op == "lt_eq":
            return self.select(F.all() <= other)
        else:
            msg = f"unexpected comparison operator {op!r}"
            raise ValueError(msg)

    def _div(self, other: Any, *, floordiv: bool) -> DataFrame:
        if isinstance(other, pl.Series):
            if floordiv:
                return self.select(F.all() // lit(other))
            return self.select(F.all() / lit(other))

        elif not isinstance(other, DataFrame):
            s = _prepare_other_arg(other, length=self.height)
            other = DataFrame([s.alias(f"n{i}") for i in range(self.width)])

        orig_dtypes = other.dtypes
        # TODO: Dispatch to a native floordiv
        other = self._cast_all_from_to(other, INTEGER_DTYPES, Float64)
        df = self._from_pydf(self._df.div_df(other._df))

        df = (
            df
            if not floordiv
            else df.with_columns([s.floor() for s in df if s.dtype.is_float()])
        )
        if floordiv:
            int_casts = [
                col(column).cast(tp)
                for i, (column, tp) in enumerate(self.schema.items())
                if tp.is_integer()
                and (orig_dtypes[i].is_integer() or orig_dtypes[i] == Null)
            ]
            if int_casts:
                return df.with_columns(int_casts)
        return df

    def _cast_all_from_to(
        self, df: DataFrame, from_: frozenset[PolarsDataType], to: PolarsDataType
    ) -> DataFrame:
        casts = [s.cast(to).alias(s.name) for s in df if s.dtype in from_]
        return df.with_columns(casts) if casts else df

    def __floordiv__(self, other: DataFrame | Series | int | float) -> DataFrame:
        return self._div(other, floordiv=True)

    def __truediv__(self, other: DataFrame | Series | int | float) -> DataFrame:
        return self._div(other, floordiv=False)

    def __bool__(self) -> NoReturn:
        msg = (
            "the truth value of a DataFrame is ambiguous"
            "\n\nHint: to check if a DataFrame contains any values, use `is_empty()`."
        )
        raise TypeError(msg)

    def __eq__(self, other: object) -> DataFrame:  # type: ignore[override]
        return self._comp(other, "eq")

    def __ne__(self, other: object) -> DataFrame:  # type: ignore[override]
        return self._comp(other, "neq")

    def __gt__(self, other: Any) -> DataFrame:
        return self._comp(other, "gt")

    def __lt__(self, other: Any) -> DataFrame:
        return self._comp(other, "lt")

    def __ge__(self, other: Any) -> DataFrame:
        return self._comp(other, "gt_eq")

    def __le__(self, other: Any) -> DataFrame:
        return self._comp(other, "lt_eq")

    def __getstate__(self) -> bytes:
        return self.serialize()

    def __setstate__(self, state: bytes) -> None:
        self._df = self.deserialize(BytesIO(state))._df

    def __mul__(self, other: DataFrame | Series | int | float) -> DataFrame:
        if isinstance(other, DataFrame):
            return self._from_pydf(self._df.mul_df(other._df))

        other = _prepare_other_arg(other)
        return self._from_pydf(self._df.mul(other._s))

    def __rmul__(self, other: int | float) -> DataFrame:
        return self * other

    def __add__(
        self, other: DataFrame | Series | int | float | bool | str
    ) -> DataFrame:
        if isinstance(other, DataFrame):
            return self._from_pydf(self._df.add_df(other._df))
        other = _prepare_other_arg(other)
        return self._from_pydf(self._df.add(other._s))

    def __radd__(
        self, other: DataFrame | Series | int | float | bool | str
    ) -> DataFrame:
        if isinstance(other, str):
            return self.select((lit(other) + F.col("*")).name.keep())
        return self + other

    def __sub__(self, other: DataFrame | Series | int | float) -> DataFrame:
        if isinstance(other, DataFrame):
            return self._from_pydf(self._df.sub_df(other._df))
        other = _prepare_other_arg(other)
        return self._from_pydf(self._df.sub(other._s))

    def __mod__(self, other: DataFrame | Series | int | float) -> DataFrame:
        if isinstance(other, DataFrame):
            return self._from_pydf(self._df.rem_df(other._df))
        other = _prepare_other_arg(other)
        return self._from_pydf(self._df.rem(other._s))

    def __str__(self) -> str:
        return self._df.as_str()

    def __repr__(self) -> str:
        return self.__str__()

    def __contains__(self, key: str) -> bool:
        return key in self.columns

    def __iter__(self) -> Iterator[Series]:
        return self.iter_columns()

    def __reversed__(self) -> Iterator[Series]:
        return reversed(self.get_columns())

    # `str` overlaps with `Sequence[str]`
    # We can ignore this but we must keep this overload ordering
    @overload
    def __getitem__(
        self, key: tuple[SingleIndexSelector, SingleColSelector]
    ) -> Any: ...

    @overload
    def __getitem__(  # type: ignore[overload-overlap]
        self, key: str | tuple[MultiIndexSelector, SingleColSelector]
    ) -> Series: ...

    @overload
    def __getitem__(
        self,
        key: (
            SingleIndexSelector
            | MultiIndexSelector
            | MultiColSelector
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> DataFrame: ...

    def __getitem__(
        self,
        key: (
            SingleIndexSelector
            | SingleColSelector
            | MultiColSelector
            | MultiIndexSelector
            | tuple[SingleIndexSelector, SingleColSelector]
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, SingleColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> DataFrame | Series | Any:
        """
        Get part of the DataFrame as a new DataFrame, Series, or scalar.

        Parameters
        ----------
        key
            Rows / columns to select. This is easiest to explain via example. Suppose
            we have a DataFrame with columns `'a'`, `'d'`, `'c'`, `'d'`. Here is what
            various types of `key` would do:

            - `df[0, 'a']` extracts the first element of column `'a'` and returns a
              scalar.
            - `df[0]` extracts the first row and returns a Dataframe.
            - `df['a']` extracts column `'a'` and returns a Series.
            - `df[0:2]` extracts the first two rows and returns a Dataframe.
            - `df[0:2, 'a']` extracts the first two rows from column `'a'` and returns
              a Series.
            - `df[0:2, 0]` extracts the first two rows from the first column and returns
              a Series.
            - `df[[0, 1], [0, 1, 2]]` extracts the first two rows and the first three
              columns and returns a Dataframe.
            - `df[0: 2, ['a', 'c']]` extracts the first two rows from columns `'a'` and
              `'c'` and returns a Dataframe.
            - `df[:, 0: 2]` extracts all rows from the first two columns and returns a
              Dataframe.
            - `df[:, 'a': 'c']` extracts all rows and all columns positioned between
              `'a'` and `'c'` *inclusive* and returns a Dataframe. In our example,
              that would extract columns `'a'`, `'d'`, and `'c'`.

        Returns
        -------
        DataFrame, Series, or scalar, depending on `key`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [1, 2, 3], "d": [4, 5, 6], "c": [1, 3, 2], "b": [7, 8, 9]}
        ... )
        >>> df[0]
        shape: (1, 4)
        ┌─────┬─────┬─────┬─────┐
        │ a   ┆ d   ┆ c   ┆ b   │
        │ --- ┆ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╪═════╡
        │ 1   ┆ 4   ┆ 1   ┆ 7   │
        └─────┴─────┴─────┴─────┘
        >>> df[0, "a"]
        1
        >>> df["a"]
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            3
        ]
        >>> df[0:2]
        shape: (2, 4)
        ┌─────┬─────┬─────┬─────┐
        │ a   ┆ d   ┆ c   ┆ b   │
        │ --- ┆ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╪═════╡
        │ 1   ┆ 4   ┆ 1   ┆ 7   │
        │ 2   ┆ 5   ┆ 3   ┆ 8   │
        └─────┴─────┴─────┴─────┘
        >>> df[0:2, "a"]
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            2
        ]
        >>> df[0:2, 0]
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            2
        ]
        >>> df[[0, 1], [0, 1, 2]]
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ d   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 4   ┆ 1   │
        │ 2   ┆ 5   ┆ 3   │
        └─────┴─────┴─────┘
        >>> df[0:2, ["a", "c"]]
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ c   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 1   │
        │ 2   ┆ 3   │
        └─────┴─────┘
        >>> df[:, 0:2]
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ d   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 4   │
        │ 2   ┆ 5   │
        │ 3   ┆ 6   │
        └─────┴─────┘
        >>> df[:, "a":"c"]
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ d   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 4   ┆ 1   │
        │ 2   ┆ 5   ┆ 3   │
        │ 3   ┆ 6   ┆ 2   │
        └─────┴─────┴─────┘
        """
        return get_df_item_by_key(self, key)

    def __setitem__(
        self,
        key: str | Sequence[int] | Sequence[str] | tuple[Any, str | int],
        value: Any,
    ) -> None:  # pragma: no cover
        """
        Modify DataFrame elements in place, using assignment syntax.

        Parameters
        ----------
        key : str | Sequence[int] | Sequence[str] | tuple[Any, str | int]
            Specifies the location(s) within the DataFrame to assign new values.
            The behavior varies based on the type of `key`:

            - Str: `df["a"] = value`:
                Not supported. Raises a `TypeError`. Use `df.with_columns(...)`
                to add or modify columns.

            - Sequence[str]: `df[["a", "b"]] = value`:
                Assigns multiple columns at once. `value` must be a 2D array-like
                structure with the same number of columns as the list
                of column names provided.

            - tuple[Any, str | int]: `df[row_idx, "a"] = value`:
                Assigns a new value to a specific element in the DataFrame, where
                `row_idx` specifies the row and `"a"` specifies the column.

            - `df[row_idx, col_idx] = value`:
                Similar to the above, but `col_idx` is the integer index of the column.

        value : Any
            The new value(s) to assign. The expected structure of `value` depends on the
            form of `key`:

            - For multiple column assignment (`df[["a", "b"]] = value`), `value` should
              be a 2D array-like object with shape (n_rows, n_columns).

            - For single element assignment (`df[row_idx, "a"] = value`), `value` should
              be a scalar.

        Raises
        ------
        TypeError
            If an unsupported assignment is attempted, such as assigning a Series
            directly to a column using `df["a"] = series`.

        ValueError
            If the shape of `value` does not match the expected shape based on `key`.

        Examples
        --------
        Sequence[str] :  `df[["a", "b"]] = value`:

        >>> import numpy as np
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        >>> df[["a", "b"]] = np.array([[10, 40], [20, 50], [30, 60]])
        >>> df
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 10  ┆ 40  │
        │ 20  ┆ 50  │
        │ 30  ┆ 60  │
        └─────┴─────┘

        tuple[Any, str | int] : `df[row_idx, "a"] = value`:

        >>> df[1, "a"] = 100
        >>> df
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 10  ┆ 40  │
        │ 100 ┆ 50  │
        │ 30  ┆ 60  │
        └─────┴─────┘

        `df[row_idx, col_idx] = value`:

        >>> df[0, 1] = 30
        >>> df
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 10  ┆ 30  │
        │ 100 ┆ 50  │
        │ 30  ┆ 60  │
        └─────┴─────┘
        """
        # df["foo"] = series
        if isinstance(key, str):
            msg = (
                "DataFrame object does not support `Series` assignment by index"
                "\n\nUse `DataFrame.with_columns`."
            )
            raise TypeError(msg)

        # df[["C", "D"]]
        elif isinstance(key, list):
            # TODO: Use python sequence constructors
            value = np.array(value)
            if value.ndim != 2:
                msg = "can only set multiple columns with 2D matrix"
                raise ValueError(msg)
            if value.shape[1] != len(key):
                msg = "matrix columns should be equal to list used to determine column names"
                raise ValueError(msg)

            # TODO: we can parallelize this by calling from_numpy
            columns = []
            for i, name in enumerate(key):
                columns.append(pl.Series(name, value[:, i]))
            self._df = self.with_columns(columns)._df

        # df[a, b]
        elif isinstance(key, tuple):
            row_selection, col_selection = key

            if (
                isinstance(row_selection, pl.Series) and row_selection.dtype == Boolean
            ) or is_bool_sequence(row_selection):
                msg = (
                    "not allowed to set DataFrame by boolean mask in the row position"
                    "\n\nConsider using `DataFrame.with_columns`."
                )
                raise TypeError(msg)

            # get series column selection
            if isinstance(col_selection, str):
                s = self.__getitem__(col_selection)
            elif isinstance(col_selection, int):
                s = self[:, col_selection]
            else:
                msg = f"unexpected column selection {col_selection!r}"
                raise TypeError(msg)

            # dispatch to __setitem__ of Series to do modification
            s[row_selection] = value

            # now find the location to place series
            # df[idx]
            if isinstance(col_selection, int):
                self.replace_column(col_selection, s)
            # df["foo"]
            elif isinstance(col_selection, str):
                self._replace(col_selection, s)
        else:
            msg = (
                f"cannot use `__setitem__` on DataFrame"
                f" with key {key!r} of type {type(key).__name__!r}"
                f" and value {value!r} of type {type(value).__name__!r}"
            )
            raise TypeError(msg)

    def __len__(self) -> int:
        return self.height

    def __copy__(self) -> DataFrame:
        return self.clone()

    def __deepcopy__(self, memo: None = None) -> DataFrame:
        return self.clone()

    def _ipython_key_completions_(self) -> list[str]:
        return self.columns

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object:
        """
        Export a DataFrame via the Arrow PyCapsule Interface.

        https://arrow.apache.org/docs/dev/format/CDataInterface/PyCapsuleInterface.html
        """
        return self._df.__arrow_c_stream__(requested_schema)

    def _repr_html_(self, *, _from_series: bool = False) -> str:
        """
        Format output data in HTML for display in Jupyter Notebooks.

        Output rows and columns can be modified by setting the following ENVIRONMENT
        variables:

        * POLARS_FMT_MAX_COLS: set the number of columns
        * POLARS_FMT_MAX_ROWS: set the number of rows
        """
        max_cols = int(os.environ.get("POLARS_FMT_MAX_COLS", default=75))
        if max_cols < 0:
            max_cols = self.width

        max_rows = int(os.environ.get("POLARS_FMT_MAX_ROWS", default=10))
        if max_rows < 0:
            max_rows = self.height

        return "".join(
            NotebookFormatter(
                self,
                max_cols=max_cols,
                max_rows=max_rows,
                from_series=_from_series,
            ).render()
        )

    def collect_schema(self) -> Schema:
        """
        Get an ordered mapping of column names to their data type.

        This is an alias for the :attr:`schema` property.

        See Also
        --------
        schema

        Notes
        -----
        This method is included to facilitate writing code that is generic for both
        DataFrame and LazyFrame.

        Examples
        --------
        Determine the schema.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.collect_schema()
        Schema({'foo': Int64, 'bar': Float64, 'ham': String})

        Access various properties of the schema using the :class:`Schema` object.

        >>> schema = df.collect_schema()
        >>> schema["bar"]
        Float64
        >>> schema.names()
        ['foo', 'bar', 'ham']
        >>> schema.dtypes()
        [Int64, Float64, String]
        >>> schema.len()
        3
        """
        return self.schema

    def item(self, row: int | None = None, column: int | str | None = None) -> Any:
        """
        Return the DataFrame as a scalar, or return the element at the given row/column.

        Parameters
        ----------
        row
            Optional row index.
        column
            Optional column index or name.

        See Also
        --------
        row : Get the values of a single row, either by index or by predicate.

        Notes
        -----
        If row/col not provided, this is equivalent to `df[0,0]`, with a check that
        the shape is (1,1). With row/col, this is equivalent to `df[row,col]`.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        >>> df.select((pl.col("a") * pl.col("b")).sum()).item()
        32
        >>> df.item(1, 1)
        5
        >>> df.item(2, "b")
        6
        """
        if row is None and column is None:
            if self.shape != (1, 1):
                msg = (
                    "can only call `.item()` if the dataframe is of shape (1, 1),"
                    " or if explicit row/col values are provided;"
                    f" frame has shape {self.shape!r}"
                )
                raise ValueError(msg)
            return self._df.to_series(0).get_index(0)

        elif row is None or column is None:
            msg = "cannot call `.item()` with only one of `row` or `column`"
            raise ValueError(msg)

        s = (
            self._df.to_series(column)
            if isinstance(column, int)
            else self._df.get_column(column)
        )
        return s.get_index_signed(row)

    @deprecate_renamed_parameter("future", "compat_level", version="1.1")
    def to_arrow(self, *, compat_level: CompatLevel | None = None) -> pa.Table:
        """
        Collect the underlying arrow arrays in an Arrow Table.

        This operation is mostly zero copy.

        Data types that do copy:
            - CategoricalType

        .. versionchanged:: 1.1
            The `future` parameter was renamed `compat_level`.

        Parameters
        ----------
        compat_level
            Use a specific compatibility level
            when exporting Polars' internal data structures.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"foo": [1, 2, 3, 4, 5, 6], "bar": ["a", "b", "c", "d", "e", "f"]}
        ... )
        >>> df.to_arrow()
        pyarrow.Table
        foo: int64
        bar: large_string
        ----
        foo: [[1,2,3,4,5,6]]
        bar: [["a","b","c","d","e","f"]]
        """
        if not self.width:  # 0x0 dataframe, cannot infer schema from batches
            return pa.table({})

        compat_level_py: int | bool
        if compat_level is None:
            compat_level_py = False
        elif isinstance(compat_level, CompatLevel):
            compat_level_py = compat_level._version

        record_batches = self._df.to_arrow(compat_level_py)
        return pa.Table.from_batches(record_batches)

    @overload
    def to_dict(self, *, as_series: Literal[True] = ...) -> dict[str, Series]: ...

    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...

    @overload
    def to_dict(
        self, *, as_series: bool
    ) -> dict[str, Series] | dict[str, list[Any]]: ...

    def to_dict(
        self, *, as_series: bool = True
    ) -> dict[str, Series] | dict[str, list[Any]]:
        """
        Convert DataFrame to a dictionary mapping column name to values.

        Parameters
        ----------
        as_series
            True -> Values are Series
            False -> Values are List[Any]

        See Also
        --------
        rows_by_key
        to_dicts

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [1, 2, 3, 4, 5],
        ...         "fruits": ["banana", "banana", "apple", "apple", "banana"],
        ...         "B": [5, 4, 3, 2, 1],
        ...         "cars": ["beetle", "audi", "beetle", "beetle", "beetle"],
        ...         "optional": [28, 300, None, 2, -30],
        ...     }
        ... )
        >>> df
        shape: (5, 5)
        ┌─────┬────────┬─────┬────────┬──────────┐
        │ A   ┆ fruits ┆ B   ┆ cars   ┆ optional │
        │ --- ┆ ---    ┆ --- ┆ ---    ┆ ---      │
        │ i64 ┆ str    ┆ i64 ┆ str    ┆ i64      │
        ╞═════╪════════╪═════╪════════╪══════════╡
        │ 1   ┆ banana ┆ 5   ┆ beetle ┆ 28       │
        │ 2   ┆ banana ┆ 4   ┆ audi   ┆ 300      │
        │ 3   ┆ apple  ┆ 3   ┆ beetle ┆ null     │
        │ 4   ┆ apple  ┆ 2   ┆ beetle ┆ 2        │
        │ 5   ┆ banana ┆ 1   ┆ beetle ┆ -30      │
        └─────┴────────┴─────┴────────┴──────────┘
        >>> df.to_dict(as_series=False)
        {'A': [1, 2, 3, 4, 5],
        'fruits': ['banana', 'banana', 'apple', 'apple', 'banana'],
        'B': [5, 4, 3, 2, 1],
        'cars': ['beetle', 'audi', 'beetle', 'beetle', 'beetle'],
        'optional': [28, 300, None, 2, -30]}
        >>> df.to_dict(as_series=True)
        {'A': shape: (5,)
        Series: 'A' [i64]
        [
            1
            2
            3
            4
            5
        ], 'fruits': shape: (5,)
        Series: 'fruits' [str]
        [
            "banana"
            "banana"
            "apple"
            "apple"
            "banana"
        ], 'B': shape: (5,)
        Series: 'B' [i64]
        [
            5
            4
            3
            2
            1
        ], 'cars': shape: (5,)
        Series: 'cars' [str]
        [
            "beetle"
            "audi"
            "beetle"
            "beetle"
            "beetle"
        ], 'optional': shape: (5,)
        Series: 'optional' [i64]
        [
            28
            300
            null
            2
            -30
        ]}
        """
        if as_series:
            return {s.name: s for s in self}
        else:
            return {s.name: s.to_list() for s in self}

    def to_dicts(self) -> list[dict[str, Any]]:
        """
        Convert every row to a dictionary of Python-native values.

        Notes
        -----
        If you have `ns`-precision temporal values you should be aware that Python
        natively only supports up to `μs`-precision; `ns`-precision values will be
        truncated to microseconds on conversion to Python. If this matters to your
        use-case you should export to a different format (such as Arrow or NumPy).

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> df.to_dicts()
        [{'foo': 1, 'bar': 4}, {'foo': 2, 'bar': 5}, {'foo': 3, 'bar': 6}]
        """
        return self.rows(named=True)

    def to_numpy(
        self,
        *,
        order: IndexOrder = "fortran",
        writable: bool = False,
        allow_copy: bool = True,
        structured: bool = False,
        use_pyarrow: bool | None = None,
    ) -> np.ndarray[Any, Any]:
        """
        Convert this DataFrame to a NumPy ndarray.

        This operation copies data only when necessary. The conversion is zero copy when
        all of the following hold:

        - The DataFrame is fully contiguous in memory, with all Series back-to-back and
          all Series consisting of a single chunk.
        - The data type is an integer or float.
        - The DataFrame contains no null values.
        - The `order` parameter is set to `fortran` (default).
        - The `writable` parameter is set to `False` (default).

        Parameters
        ----------
        order
            The index order of the returned NumPy array, either C-like or
            Fortran-like. In general, using the Fortran-like index order is faster.
            However, the C-like order might be more appropriate to use for downstream
            applications to prevent cloning data, e.g. when reshaping into a
            one-dimensional array.
        writable
            Ensure the resulting array is writable. This will force a copy of the data
            if the array was created without copy, as the underlying Arrow data is
            immutable.
        allow_copy
            Allow memory to be copied to perform the conversion. If set to `False`,
            causes conversions that are not zero-copy to fail.
        structured
            Return a `structured array`_ with a data type that corresponds to the
            DataFrame schema. If set to `False` (default), a 2D ndarray is
            returned instead.

            .. _structured array: https://numpy.org/doc/stable/user/basics.rec.html

        use_pyarrow
            Use `pyarrow.Array.to_numpy
            <https://arrow.apache.org/docs/python/generated/pyarrow.Array.html#pyarrow.Array.to_numpy>`_

            function for the conversion to NumPy if necessary.

            .. deprecated:: 0.20.28
                Polars now uses its native engine by default for conversion to NumPy.

        Examples
        --------
        Numeric data without nulls can be converted without copying data in some cases.
        The resulting array will not be writable.

        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> arr = df.to_numpy()
        >>> arr
        array([[1],
               [2],
               [3]])
        >>> arr.flags.writeable
        False

        Set `writable=True` to force data copy to make the array writable.

        >>> df.to_numpy(writable=True).flags.writeable
        True

        If the DataFrame contains different numeric data types, the resulting data type
        will be the supertype. This requires data to be copied. Integer types with
        nulls are cast to a float type with `nan` representing a null value.

        >>> df = pl.DataFrame({"a": [1, 2, None], "b": [4.0, 5.0, 6.0]})
        >>> df.to_numpy()
        array([[ 1.,  4.],
               [ 2.,  5.],
               [nan,  6.]])

        Set `allow_copy=False` to raise an error if data would be copied.

        >>> s.to_numpy(allow_copy=False)  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        RuntimeError: copy not allowed: cannot convert to a NumPy array without copying data

        Polars defaults to F-contiguous order. Use `order="c"` to force the resulting
        array to be C-contiguous.

        >>> df.to_numpy(order="c").flags.c_contiguous
        True

        DataFrames with mixed types will result in an array with an object dtype.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.5, 7.0, 8.5],
        ...         "ham": ["a", "b", "c"],
        ...     },
        ...     schema_overrides={"foo": pl.UInt8, "bar": pl.Float32},
        ... )
        >>> df.to_numpy()
        array([[1, 6.5, 'a'],
               [2, 7.0, 'b'],
               [3, 8.5, 'c']], dtype=object)

        Set `structured=True` to convert to a structured array, which can better
        preserve individual column data such as name and data type.

        >>> df.to_numpy(structured=True)
        array([(1, 6.5, 'a'), (2, 7. , 'b'), (3, 8.5, 'c')],
              dtype=[('foo', 'u1'), ('bar', '<f4'), ('ham', '<U1')])
        """  # noqa: W505
        if use_pyarrow is not None:
            issue_deprecation_warning(
                "the `use_pyarrow` parameter for `DataFrame.to_numpy` is deprecated."
                " Polars now uses its native engine by default for conversion to NumPy.",
                version="0.20.28",
            )

        if structured:
            if not allow_copy and not self.is_empty():
                msg = "copy not allowed: cannot create structured array without copying data"
                raise RuntimeError(msg)

            arrays = []
            struct_dtype = []
            for s in self.iter_columns():
                if s.dtype == Struct:
                    arr = s.struct.unnest().to_numpy(
                        structured=True,
                        allow_copy=True,
                        use_pyarrow=use_pyarrow,
                    )
                else:
                    arr = s.to_numpy(use_pyarrow=use_pyarrow)

                if s.dtype == String and not s.has_nulls():
                    arr = arr.astype(str, copy=False)
                arrays.append(arr)
                struct_dtype.append((s.name, arr.dtype, arr.shape[1:]))

            out = np.empty(self.height, dtype=struct_dtype)
            for idx, c in enumerate(self.columns):
                out[c] = arrays[idx]
            return out

        return self._df.to_numpy(order, writable=writable, allow_copy=allow_copy)

    @overload
    def to_jax(
        self,
        return_type: Literal["array"] = ...,
        *,
        device: jax.Device | str | None = ...,
        label: str | Expr | Sequence[str | Expr] | None = ...,
        features: str | Expr | Sequence[str | Expr] | None = ...,
        dtype: PolarsDataType | None = ...,
        order: IndexOrder = ...,
    ) -> jax.Array: ...

    @overload
    def to_jax(
        self,
        return_type: Literal["dict"],
        *,
        device: jax.Device | str | None = ...,
        label: str | Expr | Sequence[str | Expr] | None = ...,
        features: str | Expr | Sequence[str | Expr] | None = ...,
        dtype: PolarsDataType | None = ...,
        order: IndexOrder = ...,
    ) -> dict[str, jax.Array]: ...

    @unstable()
    def to_jax(
        self,
        return_type: JaxExportType = "array",
        *,
        device: jax.Device | str | None = None,
        label: str | Expr | Sequence[str | Expr] | None = None,
        features: str | Expr | Sequence[str | Expr] | None = None,
        dtype: PolarsDataType | None = None,
        order: IndexOrder = "fortran",
    ) -> jax.Array | dict[str, jax.Array]:
        """
        Convert DataFrame to a Jax Array, or dict of Jax Arrays.

        .. versionadded:: 0.20.27

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Parameters
        ----------
        return_type : {"array", "dict"}
            Set return type; a Jax Array, or dict of Jax Arrays.
        device
            Specify the jax `Device` on which the array will be created; can provide
            a string (such as "cpu", "gpu", or "tpu") in which case the device is
            retrieved as `jax.devices(string)[0]`. For more specific control you
            can supply the instantiated `Device` directly. If None, arrays are
            created on the default device.
        label
            One or more column names, expressions, or selectors that label the feature
            data; results in a `{"label": ..., "features": ...}` dict being returned
            when `return_type` is "dict" instead of a `{"col": array, }` dict.
        features
            One or more column names, expressions, or selectors that contain the feature
            data; if omitted, all columns that are not designated as part of the label
            are used. Only applies when `return_type` is "dict".
        dtype
            Unify the dtype of all returned arrays; this casts any column that is
            not already of the required dtype before converting to Array. Note that
            export will be single-precision (32bit) unless the Jax config/environment
            directs otherwise (eg: "jax_enable_x64" was set True in the config object
            at startup, or "JAX_ENABLE_X64" is set to "1" in the environment).
        order : {"c", "fortran"}
            The index order of the returned Jax array, either C-like (row-major) or
            Fortran-like (column-major).

        See Also
        --------
        to_dummies
        to_numpy
        to_torch

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "lbl": [0, 1, 2, 3],
        ...         "feat1": [1, 0, 0, 1],
        ...         "feat2": [1.5, -0.5, 0.0, -2.25],
        ...     }
        ... )

        Standard return type (2D Array), on the standard device:

        >>> df.to_jax()
        Array([[ 0.  ,  1.  ,  1.5 ],
               [ 1.  ,  0.  , -0.5 ],
               [ 2.  ,  0.  ,  0.  ],
               [ 3.  ,  1.  , -2.25]], dtype=float32)

        Create the Array on the default GPU device:

        >>> a = df.to_jax(device="gpu")  # doctest: +SKIP
        >>> a.device()  # doctest: +SKIP
        GpuDevice(id=0, process_index=0)

        Create the Array on a specific GPU device:

        >>> gpu_device = jax.devices("gpu")[1]  # doctest: +SKIP
        >>> a = df.to_jax(device=gpu_device)  # doctest: +SKIP
        >>> a.device()  # doctest: +SKIP
        GpuDevice(id=1, process_index=0)

        As a dictionary of individual Arrays:

        >>> df.to_jax("dict")
        {'lbl': Array([0, 1, 2, 3], dtype=int32),
         'feat1': Array([1, 0, 0, 1], dtype=int32),
         'feat2': Array([ 1.5 , -0.5 ,  0.  , -2.25], dtype=float32)}

        As a "label" and "features" dictionary; note that as "features" is not
        declared, it defaults to all the columns that are not in "label":

        >>> df.to_jax("dict", label="lbl")
        {'label': Array([[0],
                [1],
                [2],
                [3]], dtype=int32),
         'features': Array([[ 1.  ,  1.5 ],
                [ 0.  , -0.5 ],
                [ 0.  ,  0.  ],
                [ 1.  , -2.25]], dtype=float32)}

        As a "label" and "features" dictionary where each is designated using
        a col or selector expression (which can also be used to cast the data
        if the label and features are better-represented with different dtypes):

        >>> import polars.selectors as cs
        >>> df.to_jax(
        ...     return_type="dict",
        ...     features=cs.float(),
        ...     label=pl.col("lbl").cast(pl.UInt8),
        ... )
        {'label': Array([[0],
                [1],
                [2],
                [3]], dtype=uint8),
         'features': Array([[ 1.5 ],
                [-0.5 ],
                [ 0.  ],
                [-2.25]], dtype=float32)}
        """
        if return_type != "dict" and (label is not None or features is not None):
            msg = "`label` and `features` only apply when `return_type` is 'dict'"
            raise ValueError(msg)
        elif return_type == "dict" and label is None and features is not None:
            msg = "`label` is required if setting `features` when `return_type='dict'"
            raise ValueError(msg)

        jx = import_optional(
            "jax",
            install_message="Please see `https://jax.readthedocs.io/en/latest/installation.html` "
            "for specific installation recommendations for the Jax package",
        )
        enabled_double_precision = jx.config.jax_enable_x64 or bool(
            int(os.environ.get("JAX_ENABLE_X64", "0"))
        )
        if dtype:
            frame = self.cast(dtype)
        elif not enabled_double_precision:
            # enforce single-precision unless environment/config directs otherwise
            frame = self.cast({Float64: Float32, Int64: Int32, UInt64: UInt32})
        else:
            frame = self

        if isinstance(device, str):
            device = jx.devices(device)[0]

        with contextlib.nullcontext() if device is None else jx.default_device(device):
            if return_type == "array":
                # note: jax arrays are immutable, so can avoid a copy (vs torch)
                from polars.ml.utilities import frame_to_numpy

                arr = frame_to_numpy(
                    df=frame,
                    order=order,
                    writable=False,
                    target="Jax Array",
                )
                return jx.numpy.asarray(a=arr, order="K")

            elif return_type == "dict":
                if label is not None:
                    # return a {"label": array(s), "features": array(s)} dict
                    label_frame = frame.select(label)
                    features_frame = (
                        frame.select(features)
                        if features is not None
                        else frame.drop(*label_frame.columns)
                    )
                    return {
                        "label": label_frame.to_jax(),
                        "features": features_frame.to_jax(),
                    }
                else:
                    # return a {"col": array} dict
                    return {srs.name: srs.to_jax() for srs in frame}
            else:
                valid_jax_types = ", ".join(get_args(JaxExportType))
                msg = f"invalid `return_type`: {return_type!r}\nExpected one of: {valid_jax_types}"
                raise ValueError(msg)

    @overload
    def to_torch(
        self,
        return_type: Literal["tensor"] = ...,
        *,
        label: str | Expr | Sequence[str | Expr] | None = ...,
        features: str | Expr | Sequence[str | Expr] | None = ...,
        dtype: PolarsDataType | None = ...,
    ) -> torch.Tensor: ...

    @overload
    def to_torch(
        self,
        return_type: Literal["dataset"],
        *,
        label: str | Expr | Sequence[str | Expr] | None = ...,
        features: str | Expr | Sequence[str | Expr] | None = ...,
        dtype: PolarsDataType | None = ...,
    ) -> PolarsDataset: ...

    @overload
    def to_torch(
        self,
        return_type: Literal["dict"],
        *,
        label: str | Expr | Sequence[str | Expr] | None = ...,
        features: str | Expr | Sequence[str | Expr] | None = ...,
        dtype: PolarsDataType | None = ...,
    ) -> dict[str, torch.Tensor]: ...

    @unstable()
    def to_torch(
        self,
        return_type: TorchExportType = "tensor",
        *,
        label: str | Expr | Sequence[str | Expr] | None = None,
        features: str | Expr | Sequence[str | Expr] | None = None,
        dtype: PolarsDataType | None = None,
    ) -> torch.Tensor | dict[str, torch.Tensor] | PolarsDataset:
        """
        Convert DataFrame to a PyTorch Tensor, Dataset, or dict of Tensors.

        .. versionadded:: 0.20.23

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Parameters
        ----------
        return_type : {"tensor", "dataset", "dict"}
            Set return type; a PyTorch Tensor, PolarsDataset (a frame-specialized
            TensorDataset), or dict of Tensors.
        label
            One or more column names, expressions, or selectors that label the feature
            data; when `return_type` is "dataset", the PolarsDataset will return
            `(features, label)` tensor tuples for each row. Otherwise, it returns
            `(features,)` tensor tuples where the feature contains all the row data.
        features
            One or more column names, expressions, or selectors that contain the feature
            data; if omitted, all columns that are not designated as part of the label
            are used.
        dtype
            Unify the dtype of all returned tensors; this casts any column that is
            not of the required dtype before converting to Tensor. This includes
            the label column *unless* the label is an expression (such as
            `pl.col("label_column").cast(pl.Int16)`).

        See Also
        --------
        to_dummies
        to_jax
        to_numpy

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "lbl": [0, 1, 2, 3],
        ...         "feat1": [1, 0, 0, 1],
        ...         "feat2": [1.5, -0.5, 0.0, -2.25],
        ...     }
        ... )

        Standard return type (Tensor), with f32 supertype:

        >>> df.to_torch(dtype=pl.Float32)
        tensor([[ 0.0000,  1.0000,  1.5000],
                [ 1.0000,  0.0000, -0.5000],
                [ 2.0000,  0.0000,  0.0000],
                [ 3.0000,  1.0000, -2.2500]])

        As a dictionary of individual Tensors:

        >>> df.to_torch("dict")
        {'lbl': tensor([0, 1, 2, 3]),
         'feat1': tensor([1, 0, 0, 1]),
         'feat2': tensor([ 1.5000, -0.5000,  0.0000, -2.2500], dtype=torch.float64)}

        As a "label" and "features" dictionary; note that as "features" is not
        declared, it defaults to all the columns that are not in "label":

        >>> df.to_torch("dict", label="lbl", dtype=pl.Float32)
        {'label': tensor([[0.],
                 [1.],
                 [2.],
                 [3.]]),
         'features': tensor([[ 1.0000,  1.5000],
                 [ 0.0000, -0.5000],
                 [ 0.0000,  0.0000],
                 [ 1.0000, -2.2500]])}

        As a PolarsDataset, with f64 supertype:

        >>> ds = df.to_torch("dataset", dtype=pl.Float64)
        >>> ds[3]
        (tensor([ 3.0000,  1.0000, -2.2500], dtype=torch.float64),)
        >>> ds[:2]
        (tensor([[ 0.0000,  1.0000,  1.5000],
                 [ 1.0000,  0.0000, -0.5000]], dtype=torch.float64),)
        >>> ds[[0, 3]]
        (tensor([[ 0.0000,  1.0000,  1.5000],
                 [ 3.0000,  1.0000, -2.2500]], dtype=torch.float64),)

        As a convenience the PolarsDataset can opt in to half-precision data
        for experimentation (usually this would be set on the model/pipeline):

        >>> list(ds.half())
        [(tensor([0.0000, 1.0000, 1.5000], dtype=torch.float16),),
         (tensor([ 1.0000,  0.0000, -0.5000], dtype=torch.float16),),
         (tensor([2., 0., 0.], dtype=torch.float16),),
         (tensor([ 3.0000,  1.0000, -2.2500], dtype=torch.float16),)]

        Pass PolarsDataset to a DataLoader, designating the label:

        >>> from torch.utils.data import DataLoader
        >>> ds = df.to_torch("dataset", label="lbl")
        >>> dl = DataLoader(ds, batch_size=2)
        >>> batches = list(dl)
        >>> batches[0]
        [tensor([[ 1.0000,  1.5000],
                 [ 0.0000, -0.5000]], dtype=torch.float64), tensor([0, 1])]

        Note that labels can be given as expressions, allowing them to have
        a dtype independent of the feature columns (multi-column labels are
        supported).

        >>> ds = df.to_torch(
        ...     return_type="dataset",
        ...     dtype=pl.Float32,
        ...     label=pl.col("lbl").cast(pl.Int16),
        ... )
        >>> ds[:2]
        (tensor([[ 1.0000,  1.5000],
                 [ 0.0000, -0.5000]]), tensor([0, 1], dtype=torch.int16))

        Easily integrate with (for example) scikit-learn and other datasets:

        >>> from sklearn.datasets import fetch_california_housing  # doctest: +SKIP
        >>> housing = fetch_california_housing()  # doctest: +SKIP
        >>> df = pl.DataFrame(
        ...     data=housing.data,
        ...     schema=housing.feature_names,
        ... ).with_columns(
        ...     Target=housing.target,
        ... )  # doctest: +SKIP
        >>> train = df.to_torch("dataset", label="Target")  # doctest: +SKIP
        >>> loader = DataLoader(
        ...     train,
        ...     shuffle=True,
        ...     batch_size=64,
        ... )  # doctest: +SKIP
        """
        if return_type not in ("dataset", "dict") and (
            label is not None or features is not None
        ):
            msg = "`label` and `features` only apply when `return_type` is 'dataset' or 'dict'"
            raise ValueError(msg)
        elif return_type == "dict" and label is None and features is not None:
            msg = "`label` is required if setting `features` when `return_type='dict'"
            raise ValueError(msg)

        torch = import_optional("torch")

        # Cast columns.
        if dtype in (UInt16, UInt32, UInt64):
            msg = f"PyTorch does not support u16, u32, or u64 dtypes; given {dtype}"
            raise ValueError(msg)

        to_dtype = dtype or {UInt16: Int32, UInt32: Int64, UInt64: Int64}

        if label is not None:
            label_frame = self.select(label)
            # Avoid casting the label if it's an expression.
            if not isinstance(label, pl.Expr):
                label_frame = label_frame.cast(to_dtype)  # type: ignore[arg-type]
            features_frame = (
                self.select(features)
                if features is not None
                else self.drop(*label_frame.columns)
            ).cast(to_dtype)  # type: ignore[arg-type]
            frame = F.concat([label_frame, features_frame], how="horizontal")
        else:
            frame = (self.select(features) if features is not None else self).cast(
                to_dtype  # type: ignore[arg-type]
            )

        if return_type == "tensor":
            # note: torch tensors are not immutable, so we must consider them writable
            from polars.ml.utilities import frame_to_numpy

            arr = frame_to_numpy(frame, writable=True, target="Tensor")
            return torch.from_numpy(arr)

        elif return_type == "dict":
            if label is not None:
                # return a {"label": tensor(s), "features": tensor(s)} dict
                return {
                    "label": label_frame.to_torch(),
                    "features": features_frame.to_torch(),
                }
            else:
                # return a {"col": tensor} dict
                return {srs.name: srs.to_torch() for srs in frame}

        elif return_type == "dataset":
            # return a torch Dataset object
            from polars.ml.torch import PolarsDataset

            pds_label = None if label is None else label_frame.columns
            return PolarsDataset(frame, label=pds_label, features=features)
        else:
            valid_torch_types = ", ".join(get_args(TorchExportType))
            msg = f"invalid `return_type`: {return_type!r}\nExpected one of: {valid_torch_types}"
            raise ValueError(msg)

    def to_pandas(
        self,
        *,
        use_pyarrow_extension_array: bool = False,
        **kwargs: Any,
    ) -> pd.DataFrame:
        """
        Convert this DataFrame to a pandas DataFrame.

        This operation copies data if `use_pyarrow_extension_array` is not enabled.

        Parameters
        ----------
        use_pyarrow_extension_array
            Use PyArrow-backed extension arrays instead of NumPy arrays for the columns
            of the pandas DataFrame. This allows zero copy operations and preservation
            of null values. Subsequent operations on the resulting pandas DataFrame may
            trigger conversion to NumPy if those operations are not supported by PyArrow
            compute functions.
        **kwargs
            Additional keyword arguments to be passed to
            :meth:`pyarrow.Table.to_pandas`.

        Returns
        -------
        :class:`pandas.DataFrame`

        Notes
        -----
        This operation requires that both :mod:`pandas` and :mod:`pyarrow` are
        installed.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.to_pandas()
           foo  bar ham
        0    1  6.0   a
        1    2  7.0   b
        2    3  8.0   c

        Null values in numeric columns are converted to `NaN`.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, None],
        ...         "bar": [6.0, None, 8.0],
        ...         "ham": [None, "b", "c"],
        ...     }
        ... )
        >>> df.to_pandas()
           foo  bar   ham
        0  1.0  6.0  None
        1  2.0  NaN     b
        2  NaN  8.0     c

        Pass `use_pyarrow_extension_array=True` to get a pandas DataFrame with columns
        backed by PyArrow extension arrays. This will preserve null values.

        >>> df.to_pandas(use_pyarrow_extension_array=True)
            foo   bar   ham
        0     1   6.0  <NA>
        1     2  <NA>     b
        2  <NA>   8.0     c
        >>> _.dtypes
        foo           int64[pyarrow]
        bar          double[pyarrow]
        ham    large_string[pyarrow]
        dtype: object
        """
        if use_pyarrow_extension_array:
            if parse_version(pd.__version__) < (1, 5):
                msg = f'pandas>=1.5.0 is required for `to_pandas("use_pyarrow_extension_array=True")`, found Pandas {pd.__version__!r}'
                raise ModuleUpgradeRequiredError(msg)
            if not _PYARROW_AVAILABLE or parse_version(pa.__version__) < (8, 0):
                msg = "pyarrow>=8.0.0 is required for `to_pandas(use_pyarrow_extension_array=True)`"
                if _PYARROW_AVAILABLE:
                    msg += f", found pyarrow {pa.__version__!r}."
                    raise ModuleUpgradeRequiredError(msg)
                else:
                    raise ModuleNotFoundError(msg)

        # handle Object columns separately (Arrow does not convert them correctly)
        if Object in self.dtypes:
            return self._to_pandas_with_object_columns(
                use_pyarrow_extension_array=use_pyarrow_extension_array, **kwargs
            )

        return self._to_pandas_without_object_columns(
            self, use_pyarrow_extension_array=use_pyarrow_extension_array, **kwargs
        )

    def _to_pandas_with_object_columns(
        self,
        *,
        use_pyarrow_extension_array: bool,
        **kwargs: Any,
    ) -> pd.DataFrame:
        # Find which columns are of type pl.Object, and which aren't:
        object_columns = []
        not_object_columns = []
        for i, dtype in enumerate(self.dtypes):
            if dtype.is_object():
                object_columns.append(i)
            else:
                not_object_columns.append(i)

        # Export columns that aren't pl.Object, in the same order:
        if not_object_columns:
            df_without_objects = self[:, not_object_columns]
            pandas_df = self._to_pandas_without_object_columns(
                df_without_objects,
                use_pyarrow_extension_array=use_pyarrow_extension_array,
                **kwargs,
            )
        else:
            pandas_df = pd.DataFrame()

        # Add columns that are pl.Object, using Series' custom to_pandas()
        # logic for this case. We do this in order, so the original index for
        # the next column in this dataframe is correct for the partially
        # constructed Pandas dataframe, since there are no additional or
        # missing columns to the inserted column's left.
        for i in object_columns:
            name = self.columns[i]
            pandas_df.insert(i, name, self.to_series(i).to_pandas())

        return pandas_df

    def _to_pandas_without_object_columns(
        self,
        df: DataFrame,
        *,
        use_pyarrow_extension_array: bool,
        **kwargs: Any,
    ) -> pd.DataFrame:
        if not df.width:  # Empty dataframe, cannot infer schema from batches
            return pd.DataFrame()

        record_batches = df._df.to_pandas()
        tbl = pa.Table.from_batches(record_batches)
        if use_pyarrow_extension_array:
            return tbl.to_pandas(
                self_destruct=True,
                split_blocks=True,
                types_mapper=lambda pa_dtype: pd.ArrowDtype(pa_dtype),
                **kwargs,
            )

        date_as_object = kwargs.pop("date_as_object", False)
        return tbl.to_pandas(date_as_object=date_as_object, **kwargs)

    def to_series(self, index: int = 0) -> Series:
        """
        Select column as Series at index location.

        Parameters
        ----------
        index
            Location of selection.

        See Also
        --------
        get_column

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.to_series(1)
        shape: (3,)
        Series: 'bar' [i64]
        [
                6
                7
                8
        ]
        """
        return wrap_s(self._df.to_series(index))

    def to_init_repr(self, n: int = 1000) -> str:
        """
        Convert DataFrame to instantiable string representation.

        Parameters
        ----------
        n
            Only use first n rows.

        See Also
        --------
        polars.Series.to_init_repr
        polars.from_repr

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     [
        ...         pl.Series("foo", [1, 2, 3], dtype=pl.UInt8),
        ...         pl.Series("bar", [6.0, 7.0, 8.0], dtype=pl.Float32),
        ...         pl.Series("ham", ["a", "b", "c"], dtype=pl.String),
        ...     ]
        ... )
        >>> print(df.to_init_repr())
        pl.DataFrame(
            [
                pl.Series('foo', [1, 2, 3], dtype=pl.UInt8),
                pl.Series('bar', [6.0, 7.0, 8.0], dtype=pl.Float32),
                pl.Series('ham', ['a', 'b', 'c'], dtype=pl.String),
            ]
        )

        >>> df_from_str_repr = eval(df.to_init_repr())
        >>> df_from_str_repr
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ u8  ┆ f32 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6.0 ┆ a   │
        │ 2   ┆ 7.0 ┆ b   │
        │ 3   ┆ 8.0 ┆ c   │
        └─────┴─────┴─────┘
        """
        output = StringIO()
        output.write("pl.DataFrame(\n    [\n")

        for i in range(self.width):
            output.write("        ")
            output.write(self.to_series(i).to_init_repr(n))
            output.write(",\n")

        output.write("    ]\n)\n")

        return output.getvalue()

    @overload
    def serialize(
        self, file: None = ..., *, format: Literal["binary"] = ...
    ) -> bytes: ...

    @overload
    def serialize(self, file: None = ..., *, format: Literal["json"]) -> str: ...

    @overload
    def serialize(
        self, file: IOBase | str | Path, *, format: SerializationFormat = ...
    ) -> None: ...

    def serialize(
        self,
        file: IOBase | str | Path | None = None,
        *,
        format: SerializationFormat = "binary",
    ) -> bytes | str | None:
        r"""
        Serialize this DataFrame to a file or string in JSON format.

        Parameters
        ----------
        file
            File path or writable file-like object to which the result will be written.
            If set to `None` (default), the output is returned as a string instead.
        format
            The format in which to serialize. Options:

            - `"binary"`: Serialize to binary format (bytes). This is the default.
            - `"json"`: Serialize to JSON format (string).

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        Serialize the DataFrame into a binary representation.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...     }
        ... )
        >>> bytes = df.serialize()
        >>> type(bytes)
        <class 'bytes'>

        The bytes can later be deserialized back into a DataFrame.

        >>> import io
        >>> pl.DataFrame.deserialize(io.BytesIO(bytes))
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
        if format == "binary":
            serializer = self._df.serialize_binary
        elif format == "json":
            serializer = self._df.serialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return serialize_polars_object(serializer, file, format)

    @overload
    def write_json(self, file: None = ...) -> str: ...

    @overload
    def write_json(self, file: IOBase | str | Path) -> None: ...

    def write_json(self, file: IOBase | str | Path | None = None) -> str | None:
        """
        Serialize to JSON representation.

        Parameters
        ----------
        file
            File path or writable file-like object to which the result will be written.
            If set to `None` (default), the output is returned as a string instead.

        See Also
        --------
        DataFrame.write_ndjson

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...     }
        ... )
        >>> df.write_json()
        '[{"foo":1,"bar":6},{"foo":2,"bar":7},{"foo":3,"bar":8}]'
        """

        def write_json_to_string() -> str:
            with BytesIO() as buf:
                self._df.write_json(buf)
                json_bytes = buf.getvalue()
            return json_bytes.decode("utf8")

        if file is None:
            return write_json_to_string()
        elif isinstance(file, StringIO):
            json_str = write_json_to_string()
            file.write(json_str)
            return None
        elif isinstance(file, (str, Path)):
            file = normalize_filepath(file)
            self._df.write_json(file)
            return None
        else:
            self._df.write_json(file)
            return None

    @overload
    def write_ndjson(self, file: None = None) -> str: ...

    @overload
    def write_ndjson(self, file: str | Path | IO[bytes] | IO[str]) -> None: ...

    def write_ndjson(
        self, file: str | Path | IO[bytes] | IO[str] | None = None
    ) -> str | None:
        r"""
        Serialize to newline delimited JSON representation.

        Parameters
        ----------
        file
            File path or writable file-like object to which the result will be written.
            If set to `None` (default), the output is returned as a string instead.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...     }
        ... )
        >>> df.write_ndjson()
        '{"foo":1,"bar":6}\n{"foo":2,"bar":7}\n{"foo":3,"bar":8}\n'
        """
        should_return_buffer = False
        target: str | Path | IO[bytes] | IO[str]
        if file is None:
            target = cast("IO[bytes]", BytesIO())
            should_return_buffer = True
        elif isinstance(file, (str, os.PathLike)):
            target = normalize_filepath(file)
        else:
            target = file

        engine: EngineType = "in-memory"

        from polars.lazyframe.opt_flags import QueryOptFlags

        self.lazy().sink_ndjson(
            target,
            optimizations=QueryOptFlags._eager(),
            engine=engine,
        )

        if should_return_buffer:
            return str(target.getvalue(), encoding="utf-8")  # type: ignore[union-attr]

        return None

    @overload
    def write_csv(
        self,
        file: None = None,
        *,
        include_bom: bool = ...,
        include_header: bool = ...,
        separator: str = ...,
        line_terminator: str = ...,
        quote_char: str = ...,
        batch_size: int = ...,
        datetime_format: str | None = ...,
        date_format: str | None = ...,
        time_format: str | None = ...,
        float_scientific: bool | None = ...,
        float_precision: int | None = ...,
        decimal_comma: bool = ...,
        null_value: str | None = ...,
        quote_style: CsvQuoteStyle | None = ...,
        storage_options: dict[str, Any] | None = ...,
        credential_provider: CredentialProviderFunction | Literal["auto"] | None = ...,
        retries: int = ...,
    ) -> str: ...

    @overload
    def write_csv(
        self,
        file: str | Path | IO[str] | IO[bytes],
        *,
        include_bom: bool = ...,
        include_header: bool = ...,
        separator: str = ...,
        line_terminator: str = ...,
        quote_char: str = ...,
        batch_size: int = ...,
        datetime_format: str | None = ...,
        date_format: str | None = ...,
        time_format: str | None = ...,
        float_scientific: bool | None = ...,
        float_precision: int | None = ...,
        decimal_comma: bool = ...,
        null_value: str | None = ...,
        quote_style: CsvQuoteStyle | None = ...,
        storage_options: dict[str, Any] | None = ...,
        credential_provider: CredentialProviderFunction | Literal["auto"] | None = ...,
        retries: int = ...,
    ) -> None: ...

    def write_csv(
        self,
        file: str | Path | IO[str] | IO[bytes] | None = None,
        *,
        include_bom: bool = False,
        include_header: bool = True,
        separator: str = ",",
        line_terminator: str = "\n",
        quote_char: str = '"',
        batch_size: int = 1024,
        datetime_format: str | None = None,
        date_format: str | None = None,
        time_format: str | None = None,
        float_scientific: bool | None = None,
        float_precision: int | None = None,
        decimal_comma: bool = False,
        null_value: str | None = None,
        quote_style: CsvQuoteStyle | None = None,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
    ) -> str | None:
        """
        Write to comma-separated values (CSV) file.

        Parameters
        ----------
        file
            File path or writable file-like object to which the result will be written.
            If set to `None` (default), the output is returned as a string instead.
        include_bom
            Whether to include UTF-8 BOM in the CSV output.
        include_header
            Whether to include header in the CSV output.
        separator
            Separate CSV fields with this symbol.
        line_terminator
            String used to end each row.
        quote_char
            Byte to use as quoting character.
        batch_size
            Number of rows that will be processed per thread.
        datetime_format
            A format string, with the specifiers defined by the
            `chrono <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            Rust crate. If no format specified, the default fractional-second
            precision is inferred from the maximum timeunit found in the frame's
            Datetime cols (if any).
        date_format
            A format string, with the specifiers defined by the
            `chrono <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            Rust crate.
        time_format
            A format string, with the specifiers defined by the
            `chrono <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_
            Rust crate.
        float_scientific
            Whether to use scientific form always (true), never (false), or
            automatically (None) for `Float32` and `Float64` datatypes.
        float_precision
            Number of decimal places to write, applied to both `Float32` and
            `Float64` datatypes.
        decimal_comma
            Use a comma as the decimal separator instead of a point in standard
            notation. Floats will be encapsulated in quotes if necessary; set the
            field separator to override.
        null_value
            A string representing null values (defaulting to the empty string).
        quote_style : {'necessary', 'always', 'non_numeric', 'never'}
            Determines the quoting strategy used.

            - necessary (default): This puts quotes around fields only when necessary.
              They are necessary when fields contain a quote,
              separator or record terminator.
              Quotes are also necessary when writing an empty record
              (which is indistinguishable from a record with one empty field).
              This is the default.
            - always: This puts quotes around every field. Always.
            - never: This never puts quotes around fields, even if that results in
              invalid CSV data (e.g.: by not quoting strings containing the separator).
            - non_numeric: This puts quotes around all fields that are non-numeric.
              Namely, when writing a field that does not parse as a valid float
              or integer, then quotes will be used even if they aren`t strictly
              necessary.
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

        Examples
        --------
        >>> import pathlib
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> path: pathlib.Path = dirpath / "new_file.csv"
        >>> df.write_csv(path, separator=",")
        """
        from polars.io.csv._utils import _check_arg_is_1byte

        _check_arg_is_1byte("separator", separator, can_be_empty=False)
        _check_arg_is_1byte("quote_char", quote_char, can_be_empty=True)
        if not null_value:
            null_value = None

        should_return_buffer = False
        target: str | Path | IO[bytes] | IO[str]
        if file is None:
            target = cast("IO[bytes]", BytesIO())
            should_return_buffer = True
        elif isinstance(file, (str, os.PathLike)):
            target = normalize_filepath(file)
        else:
            target = file

        engine: EngineType = "in-memory"

        from polars.lazyframe.opt_flags import QueryOptFlags

        self.lazy().sink_csv(
            target,
            include_bom=include_bom,
            include_header=include_header,
            separator=separator,
            line_terminator=line_terminator,
            quote_char=quote_char,
            batch_size=batch_size,
            datetime_format=datetime_format,
            date_format=date_format,
            time_format=time_format,
            float_scientific=float_scientific,
            float_precision=float_precision,
            decimal_comma=decimal_comma,
            null_value=null_value,
            quote_style=quote_style,
            storage_options=storage_options,
            credential_provider=credential_provider,
            retries=retries,
            optimizations=QueryOptFlags._eager(),
            engine=engine,
        )

        if should_return_buffer:
            return str(target.getvalue(), encoding="utf-8")  # type: ignore[union-attr]

        return None

    def write_clipboard(self, *, separator: str = "\t", **kwargs: Any) -> None:
        """
        Copy `DataFrame` in csv format to the system clipboard with `write_csv`.

        Useful for pasting into Excel or other similar spreadsheet software.

        Parameters
        ----------
        separator
            Separate CSV fields with this symbol.
        kwargs
            Additional arguments to pass to `write_csv`.

        See Also
        --------
        polars.read_clipboard: Read a DataFrame from the clipboard.
        write_csv: Write to comma-separated values (CSV) file.
        """
        result: str = self.write_csv(file=None, separator=separator, **kwargs)
        _write_clipboard_string(result)

    def write_avro(
        self,
        file: str | Path | IO[bytes],
        compression: AvroCompression = "uncompressed",
        name: str = "",
    ) -> None:
        """
        Write to Apache Avro file.

        Parameters
        ----------
        file
            File path or writable file-like object to which the data will be written.
        compression : {'uncompressed', 'snappy', 'deflate'}
            Compression method. Defaults to "uncompressed".
        name
            Schema name. Defaults to empty string.

        Examples
        --------
        >>> import pathlib
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> path: pathlib.Path = dirpath / "new_file.avro"
        >>> df.write_avro(path)
        """
        if compression is None:
            compression = "uncompressed"
        if isinstance(file, (str, Path)):
            file = normalize_filepath(file)
        if name is None:
            name = ""

        self._df.write_avro(file, compression, name)

    def write_excel(
        self,
        workbook: str | Workbook | IO[bytes] | Path | None = None,
        worksheet: str | Worksheet | None = None,
        *,
        position: tuple[int, int] | str = "A1",
        table_style: str | dict[str, Any] | None = None,
        table_name: str | None = None,
        column_formats: ColumnFormatDict | None = None,
        dtype_formats: dict[OneOrMoreDataTypes, str] | None = None,
        conditional_formats: ConditionalFormatDict | None = None,
        header_format: dict[str, Any] | None = None,
        column_totals: ColumnTotalsDefinition | None = None,
        column_widths: ColumnWidthsDefinition | None = None,
        row_totals: RowTotalsDefinition | None = None,
        row_heights: dict[int | tuple[int, ...], int] | int | None = None,
        sparklines: dict[str, Sequence[str] | dict[str, Any]] | None = None,
        formulas: dict[str, str | dict[str, str]] | None = None,
        float_precision: int = 3,
        include_header: bool = True,
        autofilter: bool = True,
        autofit: bool = False,
        hidden_columns: Sequence[str] | SelectorType | None = None,
        hide_gridlines: bool = False,
        sheet_zoom: int | None = None,
        freeze_panes: (
            str
            | tuple[int, int]
            | tuple[str, int, int]
            | tuple[int, int, int, int]
            | None
        ) = None,
    ) -> Workbook:
        """
        Write frame data to a table in an Excel workbook/worksheet.

        Parameters
        ----------
        workbook : {str, Workbook}
            String name or path of the workbook to create, BytesIO object, file opened
            in binary-mode, or an `xlsxwriter.Workbook` object that has not been closed.
            If None, writes to `dataframe.xlsx` in the working directory.
        worksheet : {str, Worksheet}
            Name of target worksheet or an `xlsxwriter.Worksheet` object (in which
            case `workbook` must be the parent `xlsxwriter.Workbook` object); if None,
            writes to "Sheet1" when creating a new workbook (note that writing to an
            existing workbook requires a valid existing -or new- worksheet name).
        position : {str, tuple}
            Table position in Excel notation (eg: "A1"), or a (row,col) integer tuple.
        table_style : {str, dict}
            A named Excel table style, such as "Table Style Medium 4", or a dictionary
            of `{"key":value,}` options containing one or more of the following keys:
            "style", "first_column", "last_column", "banded_columns, "banded_rows".
        table_name : str
            Name of the output table object in the worksheet; can then be referred to
            in the sheet by formulae/charts, or by subsequent `xlsxwriter` operations.
        column_formats : dict
            A `{colname(s):str,}` or `{selector:str,}` dictionary for applying an
            Excel format string to the given columns. Formats defined here (such as
            "dd/mm/yyyy", "0.00%", etc) will override any defined in `dtype_formats`.
        dtype_formats : dict
            A `{dtype:str,}` dictionary that sets the default Excel format for the
            given dtype. (This can be overridden on a per-column basis by the
            `column_formats` param).
        conditional_formats : dict
            A dictionary of colname (or selector) keys to a format str, dict, or list
            that defines conditional formatting options for the specified columns.

            * If supplying a string typename, should be one of the valid `xlsxwriter`
              types such as "3_color_scale", "data_bar", etc.
            * If supplying a dictionary you can make use of any/all `xlsxwriter`
              supported options, including icon sets, formulae, etc.
            * Supplying multiple columns as a tuple/key will apply a single format
              across all columns - this is effective in creating a heatmap, as the
              min/max values will be determined across the entire range, not per-column.
            * Finally, you can also supply a list made up from the above options
              in order to apply *more* than one conditional format to the same range.
        header_format : dict
            A `{key:value,}` dictionary of `xlsxwriter` format options to apply
            to the table header row, such as `{"bold":True, "font_color":"#702963"}`.
        column_totals : {bool, list, dict}
            Add a column-total row to the exported table.

            * If True, all numeric columns will have an associated total using "sum".
            * If passing a string, it must be one of the valid total function names
              and all numeric columns will have an associated total using that function.
            * If passing a list of colnames, only those given will have a total.
            * For more control, pass a `{colname:funcname,}` dict.

            Valid column-total function names are "average", "count_nums", "count",
            "max", "min", "std_dev", "sum", and "var".
        column_widths : {dict, int}
            A `{colname:int,}` or `{selector:int,}` dict or a single integer that
            sets (or overrides if autofitting) table column widths, in integer pixel
            units. If given as an integer the same value is used for all table columns.
        row_totals : {dict, list, bool}
            Add a row-total column to the right-hand side of the exported table.

            * If True, a column called "total" will be added at the end of the table
              that applies a "sum" function row-wise across all numeric columns.
            * If passing a list/sequence of column names, only the matching columns
              will participate in the sum.
            * Can also pass a `{colname:columns,}` dictionary to create one or
              more total columns with distinct names, referencing different columns.
        row_heights : {dict, int}
            An int or `{row_index:int,}` dictionary that sets the height of the given
            rows (if providing a dictionary) or all rows (if providing an integer) that
            intersect with the table body (including any header and total row) in
            integer pixel units. Note that `row_index` starts at zero and will be
            the header row (unless `include_header` is False).
        sparklines : dict
            A `{colname:list,}` or `{colname:dict,}` dictionary defining one or more
            sparklines to be written into a new column in the table.

            * If passing a list of colnames (used as the source of the sparkline data)
              the default sparkline settings are used (eg: line chart with no markers).
            * For more control an `xlsxwriter`-compliant options dict can be supplied,
              in which case three additional polars-specific keys are available:
              "columns", "insert_before", and "insert_after". These allow you to define
              the source columns and position the sparkline(s) with respect to other
              table columns. If no position directive is given, sparklines are added to
              the end of the table (eg: to the far right) in the order they are given.
        formulas : dict
            A `{colname:formula,}` or `{colname:dict,}` dictionary defining one or
            more formulas to be written into a new column in the table. Note that you
            are strongly advised to use structured references in your formulae wherever
            possible to make it simple to reference columns by name.

            * If providing a string formula (such as "=[@colx]*[@coly]") the column will
              be added to the end of the table (eg: to the far right), after any default
              sparklines and before any row_totals.
            * For the most control supply an options dictionary with the following keys:
              "formula" (mandatory), one of "insert_before" or "insert_after", and
              optionally "return_dtype". The latter is used to appropriately format the
              output of the formula and allow it to participate in row/column totals.
        float_precision : int
            Default number of decimals displayed for floating point columns (note that
            this is purely a formatting directive; the actual values are not rounded).
        include_header : bool
            Indicate if the table should be created with a header row.
        autofilter : bool
            If the table has headers, provide autofilter capability.
        autofit : bool
            Calculate individual column widths from the data.
        hidden_columns : str | list
             A column name, list of column names, or a selector representing table
             columns to mark as hidden in the output worksheet.
        hide_gridlines : bool
            Do not display any gridlines on the output worksheet.
        sheet_zoom : int
            Set the default zoom level of the output worksheet.
        freeze_panes : str | (str, int, int) | (int, int) | (int, int, int, int)
            Freeze workbook panes.

            * If (row, col) is supplied, panes are split at the top-left corner of the
              specified cell, which are 0-indexed. Thus, to freeze only the top row,
              supply (1, 0).
            * Alternatively, cell notation can be used to supply the cell. For example,
              "A2" indicates the split occurs at the top-left of cell A2, which is the
              equivalent of (1, 0).
            * If (row, col, top_row, top_col) are supplied, the panes are split based on
              the `row` and `col`, and the scrolling region is initialized to begin at
              the `top_row` and `top_col`. Thus, to freeze only the top row and have the
              scrolling region begin at row 10, column D (5th col), supply (1, 0, 9, 4).
              Using cell notation for (row, col), supplying ("A2", 9, 4) is equivalent.

        Notes
        -----
        * A list of compatible `xlsxwriter` format property names can be found here:
          https://xlsxwriter.readthedocs.io/format.html#format-methods-and-format-properties

        * Conditional formatting dictionaries should provide xlsxwriter-compatible
          definitions; polars will take care of how they are applied on the worksheet
          with respect to the relative sheet/column position. For supported options,
          see: https://xlsxwriter.readthedocs.io/working_with_conditional_formats.html

        * Similarly, sparkline option dictionaries should contain xlsxwriter-compatible
          key/values, as well as a mandatory polars "columns" key that defines the
          sparkline source data; these source columns should all be adjacent. Two other
          polars-specific keys are available to help define where the sparkline appears
          in the table: "insert_after", and "insert_before". The value associated with
          these keys should be the name of a column in the exported table.
          https://xlsxwriter.readthedocs.io/working_with_sparklines.html

        * Formula dictionaries *must* contain a key called "formula", and then optional
          "insert_after", "insert_before", and/or "return_dtype" keys. These additional
          keys allow the column to be injected into the table at a specific location,
          and/or to define the return type of the formula (eg: "Int64", "Float64", etc).
          Formulas that refer to table columns should use Excel's structured references
          syntax to ensure the formula is applied correctly and is table-relative.
          https://support.microsoft.com/en-us/office/using-structured-references-with-excel-tables-f5ed2452-2337-4f71-bed3-c8ae6d2b276e

        * If you want unformatted output, you can use a selector to apply the "General"
          format to all columns (or all *non-temporal* columns to preserve formatting
          of date/datetime columns), eg: `column_formats={~cs.temporal(): "General"}`.

        Examples
        --------
        Instantiate a basic DataFrame:

        >>> from random import uniform
        >>> from datetime import date
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "dtm": [date(2023, 1, 1), date(2023, 1, 2), date(2023, 1, 3)],
        ...         "num": [uniform(-500, 500), uniform(-500, 500), uniform(-500, 500)],
        ...         "val": [10_000, 20_000, 30_000],
        ...     }
        ... )

        Export to "dataframe.xlsx" (the default workbook name, if not specified) in the
        working directory, add column totals on all numeric columns ("sum" by default),
        then autofit:

        >>> df.write_excel(column_totals=True, autofit=True)  # doctest: +SKIP

        Write frame to a specific location on the sheet, set a named table style,
        apply US-style date formatting, increase floating point formatting precision,
        apply a non-default column total function to a specific column, autofit:

        >>> df.write_excel(  # doctest: +SKIP
        ...     position="B4",
        ...     table_style="Table Style Light 16",
        ...     dtype_formats={pl.Date: "mm/dd/yyyy"},
        ...     column_totals={"num": "average"},
        ...     float_precision=6,
        ...     autofit=True,
        ... )

        Write the same frame to a named worksheet twice, applying different styles
        and conditional formatting to each table, adding custom-formatted table
        titles using explicit `xlsxwriter` integration:

        >>> from xlsxwriter import Workbook
        >>> with Workbook("multi_frame.xlsx") as wb:  # doctest: +SKIP
        ...     # basic/default conditional formatting
        ...     df.write_excel(
        ...         workbook=wb,
        ...         worksheet="data",
        ...         position=(3, 1),  # specify position as (row,col) coordinates
        ...         conditional_formats={"num": "3_color_scale", "val": "data_bar"},
        ...         table_style="Table Style Medium 4",
        ...     )
        ...
        ...     # advanced conditional formatting, custom styles
        ...     df.write_excel(
        ...         workbook=wb,
        ...         worksheet="data",
        ...         position=(df.height + 7, 1),
        ...         table_style={
        ...             "style": "Table Style Light 4",
        ...             "first_column": True,
        ...         },
        ...         conditional_formats={
        ...             "num": {
        ...                 "type": "3_color_scale",
        ...                 "min_color": "#76933c",
        ...                 "mid_color": "#c4d79b",
        ...                 "max_color": "#ebf1de",
        ...             },
        ...             "val": {
        ...                 "type": "data_bar",
        ...                 "data_bar_2010": True,
        ...                 "bar_color": "#9bbb59",
        ...                 "bar_negative_color_same": True,
        ...                 "bar_negative_border_color_same": True,
        ...             },
        ...         },
        ...         column_formats={"num": "#,##0.000;[White]-#,##0.000"},
        ...         column_widths={"val": 125},
        ...         autofit=True,
        ...     )
        ...
        ...     # add some table titles (with a custom format)
        ...     ws = wb.get_worksheet_by_name("data")
        ...     fmt_title = wb.add_format(
        ...         {
        ...             "font_color": "#4f6228",
        ...             "font_size": 12,
        ...             "italic": True,
        ...             "bold": True,
        ...         }
        ...     )
        ...     ws.write(2, 1, "Basic/default conditional formatting", fmt_title)
        ...     ws.write(df.height + 6, 1, "Custom conditional formatting", fmt_title)

        Export a table containing two different types of sparklines. Use default
        options for the "trend" sparkline and customized options (and positioning)
        for the "+/-" `win_loss` sparkline, with non-default integer formatting,
        column totals, a subtle two-tone heatmap and hidden worksheet gridlines:

        >>> df = pl.DataFrame(
        ...     {
        ...         "id": ["aaa", "bbb", "ccc", "ddd", "eee"],
        ...         "q1": [100, 55, -20, 0, 35],
        ...         "q2": [30, -10, 15, 60, 20],
        ...         "q3": [-50, 0, 40, 80, 80],
        ...         "q4": [75, 55, 25, -10, -55],
        ...     }
        ... )
        >>> df.write_excel(  # doctest: +SKIP
        ...     table_style="Table Style Light 2",
        ...     # apply accounting format to all flavours of integer
        ...     dtype_formats={dt: "#,##0_);(#,##0)" for dt in [pl.Int32, pl.Int64]},
        ...     sparklines={
        ...         # default options; just provide source cols
        ...         "trend": ["q1", "q2", "q3", "q4"],
        ...         # customized sparkline type, with positioning directive
        ...         "+/-": {
        ...             "columns": ["q1", "q2", "q3", "q4"],
        ...             "insert_after": "id",
        ...             "type": "win_loss",
        ...         },
        ...     },
        ...     conditional_formats={
        ...         # create a unified multi-column heatmap
        ...         ("q1", "q2", "q3", "q4"): {
        ...             "type": "2_color_scale",
        ...             "min_color": "#95b3d7",
        ...             "max_color": "#ffffff",
        ...         },
        ...     },
        ...     column_totals=["q1", "q2", "q3", "q4"],
        ...     row_totals=True,
        ...     hide_gridlines=True,
        ... )

        Export a table containing an Excel formula-based column that calculates a
        standardised Z-score, showing use of structured references in conjunction
        with positioning directives, column totals, and custom formatting.

        >>> df = pl.DataFrame(
        ...     {
        ...         "id": ["a123", "b345", "c567", "d789", "e101"],
        ...         "points": [99, 45, 50, 85, 35],
        ...     }
        ... )
        >>> df.write_excel(  # doctest: +SKIP
        ...     table_style={
        ...         "style": "Table Style Medium 15",
        ...         "first_column": True,
        ...     },
        ...     column_formats={
        ...         "id": {"font": "Consolas"},
        ...         "points": {"align": "center"},
        ...         "z-score": {"align": "center"},
        ...     },
        ...     column_totals="average",
        ...     formulas={
        ...         "z-score": {
        ...             # use structured references to refer to the table columns and 'totals' row
        ...             "formula": "=STANDARDIZE([@points], [[#Totals],[points]], STDEV([points]))",
        ...             "insert_after": "points",
        ...             "return_dtype": pl.Float64,
        ...         }
        ...     },
        ...     hide_gridlines=True,
        ...     sheet_zoom=125,
        ... )

        Create and reference a Worksheet object directly, adding a basic chart.
        Taking advantage of structured references to set chart series values and
        categories is *strongly* recommended so you do not have to calculate
        cell positions with respect to the frame data and worksheet:

        >>> with Workbook("basic_chart.xlsx") as wb:  # doctest: +SKIP
        ...     # create worksheet object and write frame data to it
        ...     ws = wb.add_worksheet("demo")
        ...     df.write_excel(
        ...         workbook=wb,
        ...         worksheet=ws,
        ...         table_name="DataTable",
        ...         table_style="Table Style Medium 26",
        ...         hide_gridlines=True,
        ...     )
        ...     # create chart object, point to the written table
        ...     # data using structured references, and style it
        ...     chart = wb.add_chart({"type": "column"})
        ...     chart.set_title({"name": "Example Chart"})
        ...     chart.set_legend({"none": True})
        ...     chart.set_style(38)
        ...     chart.add_series(
        ...         {  # note the use of structured references
        ...             "values": "=DataTable[points]",
        ...             "categories": "=DataTable[id]",
        ...             "data_labels": {"value": True},
        ...         }
        ...     )
        ...     # add chart to the worksheet
        ...     ws.insert_chart("D1", chart)

        Export almost entirely unformatted data (no numeric styling or standardised
        floating point precision), omit autofilter, but keep date/datetime formatting:

        >>> import polars.selectors as cs
        >>> df = pl.DataFrame(
        ...     {
        ...         "n1": [-100, None, 200, 555],
        ...         "n2": [987.4321, -200, 44.444, 555.5],
        ...     }
        ... )
        >>> df.write_excel(  # doctest: +SKIP
        ...     column_formats={~cs.temporal(): "General"},
        ...     autofilter=False,
        ... )
        """  # noqa: W505
        from polars.io.spreadsheet._write_utils import (
            _unpack_multi_column_dict,
            _xl_apply_conditional_formats,
            _xl_inject_sparklines,
            _xl_setup_table_columns,
            _xl_setup_table_options,
            _xl_setup_workbook,
            _xl_unique_table_name,
            _XLFormatCache,
        )

        xlsxwriter = import_optional("xlsxwriter", err_prefix="Excel export requires")
        from xlsxwriter.utility import xl_cell_to_rowcol

        # setup workbook/worksheet
        wb, ws, can_close = _xl_setup_workbook(workbook, worksheet)
        df, is_empty = self, self.is_empty()

        # note: `_xl_setup_table_columns` converts nested data (List, Struct, etc.) to
        # string, so we keep a reference to the original so that column selection with
        # selectors that target such types remains correct
        df_original = df

        # setup table format/columns
        fmt_cache = _XLFormatCache(wb)
        column_formats = column_formats or {}
        table_style, table_options = _xl_setup_table_options(table_style)
        table_name = table_name or _xl_unique_table_name(wb)
        table_columns, column_formats, df = _xl_setup_table_columns(  # type: ignore[assignment]
            df=df,
            format_cache=fmt_cache,
            column_formats=column_formats,
            column_totals=column_totals,
            dtype_formats=dtype_formats,
            header_format=header_format,
            float_precision=float_precision,
            table_style=table_style,
            row_totals=row_totals,
            sparklines=sparklines,
            formulas=formulas,
        )

        # normalise cell refs (eg: "B3" => (2,1)) and establish table start/finish,
        # accounting for potential presence/absence of headers and a totals row.
        table_start = (
            xl_cell_to_rowcol(position) if isinstance(position, str) else position
        )
        table_finish = (
            table_start[0]
            + df.height
            + int(is_empty)
            - int(not include_header)
            + int(bool(column_totals)),
            table_start[1] + df.width - 1,
        )

        excel_max_valid_rows = 1048575
        excel_max_valid_cols = 16384

        if (
            table_finish[0] > excel_max_valid_rows
            or table_finish[1] > excel_max_valid_cols
        ):
            msg = f"writing {df.height}x{df.width} frame at {position!r} does not fit worksheet dimensions of {excel_max_valid_rows} rows and {excel_max_valid_cols} columns"
            raise InvalidOperationError(msg)

        # write table structure and formats into the target sheet
        if not is_empty or include_header:
            ws.add_table(
                *table_start,
                *table_finish,
                {
                    "data": df.rows(),
                    "style": table_style,
                    "columns": table_columns,
                    "header_row": include_header,
                    "autofilter": autofilter,
                    "total_row": bool(column_totals) and not is_empty,
                    "name": table_name,
                    **table_options,
                },
            )

            # apply conditional formats
            if conditional_formats:
                _xl_apply_conditional_formats(
                    df=df,
                    ws=ws,
                    conditional_formats=conditional_formats,
                    table_start=table_start,
                    include_header=include_header,
                    format_cache=fmt_cache,
                )

        # additional column-level properties
        if hidden_columns is None:
            hidden = set()
        elif isinstance(hidden_columns, str):
            hidden = {hidden_columns}
        else:
            hidden = set(_expand_selectors(df_original, hidden_columns))

        # Autofit section needs to be present above column_widths section
        # to ensure that parameters provided in the column_widths section
        # are not overwritten by autofit
        #
        # table/rows all written; apply (optional) autofit
        if autofit and not is_empty:
            xlv = xlsxwriter.__version__
            if parse_version(xlv) < (3, 0, 8):
                msg = f"`autofit=True` requires xlsxwriter 3.0.8 or higher, found {xlv}"
                raise ModuleUpgradeRequiredError(msg)
            ws.autofit()

        if isinstance(column_widths, int):
            column_widths = dict.fromkeys(df.columns, column_widths)
        else:
            column_widths = _expand_selector_dicts(  # type: ignore[assignment]
                df_original, column_widths, expand_keys=True, expand_values=False
            )
        column_widths = _unpack_multi_column_dict(column_widths or {})  # type: ignore[assignment]

        for column in df.columns:
            options = {"hidden": True} if column in hidden else {}
            col_idx = table_start[1] + df.get_column_index(column)
            if column in column_widths:  # type: ignore[operator]
                ws.set_column_pixels(
                    col_idx,
                    col_idx,
                    column_widths[column],  # type: ignore[index]
                    None,
                    options,
                )
            elif options:
                ws.set_column(col_idx, col_idx, None, None, options)

        # finally, inject any sparklines into the table
        for column, params in (sparklines or {}).items():
            _xl_inject_sparklines(
                ws,
                df,
                table_start,
                column,
                include_header=include_header,
                params=params,
            )

        # worksheet options
        if hide_gridlines:
            ws.hide_gridlines(2)
        if sheet_zoom:
            ws.set_zoom(sheet_zoom)
        if row_heights:
            if isinstance(row_heights, int):
                for idx in range(table_start[0], table_finish[0] + 1):
                    ws.set_row_pixels(idx, row_heights)
            elif isinstance(row_heights, dict):
                for idx, height in _unpack_multi_column_dict(row_heights).items():  # type: ignore[assignment]
                    ws.set_row_pixels(idx, height)

        if freeze_panes:
            if isinstance(freeze_panes, str):
                ws.freeze_panes(freeze_panes)
            else:
                ws.freeze_panes(*freeze_panes)

        if can_close:
            wb.close()
        return wb

    @overload
    def write_ipc(
        self,
        file: None,
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
    ) -> BytesIO: ...

    @overload
    def write_ipc(
        self,
        file: str | Path | IO[bytes],
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
    ) -> None: ...

    @deprecate_renamed_parameter("future", "compat_level", version="1.1")
    def write_ipc(
        self,
        file: str | Path | IO[bytes] | None,
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
    ) -> BytesIO | None:
        """
        Write to Arrow IPC binary stream or Feather file.

        See "File or Random Access format" in https://arrow.apache.org/docs/python/ipc.html.

        .. versionchanged:: 1.1
            The `future` parameter was renamed `compat_level`.

        Parameters
        ----------
        file
            Path or writable file-like object to which the IPC data will be
            written. If set to `None`, the output is returned as a BytesIO object.
        compression : {'uncompressed', 'lz4', 'zstd'}
            Compression method. Defaults to "uncompressed".
        compat_level
            Use a specific compatibility level
            when exporting Polars' internal data structures.
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

        Examples
        --------
        >>> import pathlib
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> path: pathlib.Path = dirpath / "new_file.arrow"
        >>> df.write_ipc(path)
        """
        return_bytes = file is None
        target: str | Path | IO[bytes]
        if file is None:
            target = BytesIO()
        else:
            target = file

        from polars.lazyframe.opt_flags import QueryOptFlags

        self.lazy().sink_ipc(
            target,
            compression=compression,
            compat_level=compat_level,
            storage_options=storage_options,
            credential_provider=credential_provider,
            retries=retries,
            optimizations=QueryOptFlags._eager(),
            engine="in-memory",
        )
        return target if return_bytes else None  # type: ignore[return-value]

    @overload
    def write_ipc_stream(
        self,
        file: None,
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
    ) -> BytesIO: ...

    @overload
    def write_ipc_stream(
        self,
        file: str | Path | IO[bytes],
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
    ) -> None: ...

    @deprecate_renamed_parameter("future", "compat_level", version="1.1")
    def write_ipc_stream(
        self,
        file: str | Path | IO[bytes] | None,
        *,
        compression: IpcCompression = "uncompressed",
        compat_level: CompatLevel | None = None,
    ) -> BytesIO | None:
        """
        Write to Arrow IPC record batch stream.

        See "Streaming format" in https://arrow.apache.org/docs/python/ipc.html.

        .. versionchanged:: 1.1
            The `future` parameter was renamed `compat_level`.

        Parameters
        ----------
        file
            Path or writable file-like object to which the IPC record batch data will
            be written. If set to `None`, the output is returned as a BytesIO object.
        compression : {'uncompressed', 'lz4', 'zstd'}
            Compression method. Defaults to "uncompressed".
        compat_level
            Use a specific compatibility level
            when exporting Polars' internal data structures.

        Examples
        --------
        >>> import pathlib
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> path: pathlib.Path = dirpath / "new_file.arrow"
        >>> df.write_ipc_stream(path)
        """
        return_bytes = file is None
        if return_bytes:
            file = BytesIO()
        elif isinstance(file, (str, Path)):
            file = normalize_filepath(file)

        compat_level_py: int | bool
        if compat_level is None:
            compat_level_py = True
        elif isinstance(compat_level, CompatLevel):
            compat_level_py = compat_level._version

        if compression is None:
            compression = "uncompressed"

        self._df.write_ipc_stream(file, compression, compat_level_py)
        return file if return_bytes else None  # type: ignore[return-value]

    def write_parquet(
        self,
        file: str | Path | IO[bytes],
        *,
        compression: ParquetCompression = "zstd",
        compression_level: int | None = None,
        statistics: bool | str | dict[str, bool] = True,
        row_group_size: int | None = None,
        data_page_size: int | None = None,
        use_pyarrow: bool = False,
        pyarrow_options: dict[str, Any] | None = None,
        partition_by: str | Sequence[str] | None = None,
        partition_chunk_size_bytes: int = 4_294_967_296,
        storage_options: dict[str, Any] | None = None,
        credential_provider: (
            CredentialProviderFunction | Literal["auto"] | None
        ) = "auto",
        retries: int = 2,
        metadata: ParquetMetadata | None = None,
        mkdir: bool = False,
    ) -> None:
        """
        Write to Apache Parquet file.

        Parameters
        ----------
        file
            File path or writable file-like object to which the result will be written.
            This should be a path to a directory if writing a partitioned dataset.
        compression : {'lz4', 'uncompressed', 'snappy', 'gzip', 'lzo', 'brotli', 'zstd'}
            Choose "zstd" for good compression performance.
            Choose "lz4" for fast compression/decompression.
            Choose "snappy" for more backwards compatibility guarantees
            when you deal with older parquet readers.
        compression_level
            The level of compression to use. Higher compression means smaller files on
            disk.

            - "gzip" : min-level: 0, max-level: 9, default: 6.
            - "brotli" : min-level: 0, max-level: 11, default: 1.
            - "zstd" : min-level: 1, max-level: 22, default: 3.

        statistics
            Write statistics to the parquet headers. This is the default behavior.

            Possible values:

            - `True`: enable default set of statistics (default). Some
              statistics may be disabled.
            - `False`: disable all statistics
            - "full": calculate and write all available statistics. Cannot be
              combined with `use_pyarrow`.
            - `{ "statistic-key": True / False, ... }`. Cannot be combined with
              `use_pyarrow`. Available keys:

              - "min": column minimum value (default: `True`)
              - "max": column maximum value (default: `True`)
              - "distinct_count": number of unique column values (default: `False`)
              - "null_count": number of null values in column (default: `True`)
        row_group_size
            Size of the row groups in number of rows. Defaults to 512^2 rows.
        data_page_size
            Size of the data page in bytes. Defaults to 1024^2 bytes.
        use_pyarrow
            Use C++ parquet implementation vs Rust parquet implementation.
            At the moment C++ supports more features.
        pyarrow_options
            Arguments passed to `pyarrow.parquet.write_table`.

            If you pass `partition_cols` here, the dataset will be written
            using `pyarrow.parquet.write_to_dataset`.
            The `partition_cols` parameter leads to write the dataset to a directory.
            Similar to Spark's partitioned datasets.
        partition_by
            Column(s) to partition by. A partitioned dataset will be written if this is
            specified. This parameter is considered unstable and is subject to change.
        partition_chunk_size_bytes
            Approximate size to split DataFrames within a single partition when
            writing. Note this is calculated using the size of the DataFrame in
            memory - the size of the output file may differ depending on the
            file format / compression.
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
        metadata
            A dictionary or callback to add key-values to the file-level Parquet
            metadata.

            .. warning::
                This functionality is considered **experimental**. It may be removed or
                changed at any point without it being considered a breaking change.
        mkdir: bool
            Recursively create all the directories in the path.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.

        Examples
        --------
        >>> import pathlib
        >>>
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> path: pathlib.Path = dirpath / "new_file.parquet"
        >>> df.write_parquet(path)

        We can use pyarrow with use_pyarrow_write_to_dataset=True
        to write partitioned datasets. The following example will
        write the first row to ../watermark=1/*.parquet and the
        other rows to ../watermark=2/*.parquet.

        >>> df = pl.DataFrame({"a": [1, 2, 3], "watermark": [1, 2, 2]})
        >>> path: pathlib.Path = dirpath / "partitioned_object"
        >>> df.write_parquet(
        ...     path,
        ...     use_pyarrow=True,
        ...     pyarrow_options={"partition_cols": ["watermark"]},
        ... )
        """
        if compression is None:
            compression = "uncompressed"
        if isinstance(file, (str, Path)):
            if partition_by is not None or (
                pyarrow_options is not None and pyarrow_options.get("partition_cols")
            ):
                file = normalize_filepath(file, check_not_directory=False)
            else:
                file = normalize_filepath(file)

        if use_pyarrow:
            if statistics == "full" or isinstance(statistics, dict):
                msg = "write_parquet with `use_pyarrow=True` allows only boolean values for `statistics`"
                raise ValueError(msg)
            if metadata is not None:
                msg = "write_parquet with `use_pyarrow=True` cannot be combined with `metadata`"
                raise ValueError(msg)
            if mkdir:
                msg = "write_parquet with `use_pyarrow=True` cannot be combined with `mkdir`"
                raise ValueError(msg)

            tbl = self.to_arrow()
            data = {}

            for i, column in enumerate(tbl):
                # extract the name before casting
                name = f"column_{i}" if column._name is None else column._name

                data[name] = column

            tbl = pa.table(data)

            # do not remove this import!
            # needed below
            import pyarrow.parquet  # noqa: F401

            if pyarrow_options is None:
                pyarrow_options = {}
            pyarrow_options["compression"] = (
                None if compression == "uncompressed" else compression
            )
            pyarrow_options["compression_level"] = compression_level
            pyarrow_options["write_statistics"] = statistics
            pyarrow_options["row_group_size"] = row_group_size
            pyarrow_options["data_page_size"] = data_page_size

            if pyarrow_options.get("partition_cols"):
                pa.parquet.write_to_dataset(
                    table=tbl,
                    root_path=file,
                    **(pyarrow_options or {}),
                )
            else:
                pa.parquet.write_table(
                    table=tbl,
                    where=file,
                    **(pyarrow_options or {}),
                )

            return

        target: str | Path | IO[bytes] | PartitioningScheme = file
        engine: EngineType = "in-memory"
        if partition_by is not None:
            if not isinstance(file, str):
                msg = "expected file to be a `str` since partition-by is set"
                raise TypeError(msg)

            from polars.io import PartitionByKey

            target = PartitionByKey(file, by=partition_by)
            mkdir = True
            engine = "streaming"

        from polars.lazyframe.opt_flags import QueryOptFlags

        self.lazy().sink_parquet(
            target,
            compression=compression,
            compression_level=compression_level,
            statistics=statistics,
            row_group_size=row_group_size,
            data_page_size=data_page_size,
            storage_options=storage_options,
            credential_provider=credential_provider,
            retries=retries,
            metadata=metadata,
            engine=engine,
            mkdir=mkdir,
            optimizations=QueryOptFlags._eager(),
        )

    def write_database(
        self,
        table_name: str,
        connection: ConnectionOrCursor | str,
        *,
        if_table_exists: DbWriteMode = "fail",
        engine: DbWriteEngine | None = None,
        engine_options: dict[str, Any] | None = None,
    ) -> int:
        """
        Write the data in a Polars DataFrame to a database.

        .. versionadded:: 0.20.26
            Support for instantiated connection objects in addition to URI strings, and
            a new `engine_options` parameter.

        Parameters
        ----------
        table_name
            Schema-qualified name of the table to create or append to in the target
            SQL database. If your table name contains special characters, it should
            be quoted.
        connection
            An existing SQLAlchemy or ADBC connection against the target database, or
            a URI string that will be used to instantiate such a connection, such as:

            * "postgresql://user:pass@server:port/database"
            * "sqlite:////path/to/database.db"
        if_table_exists : {'append', 'replace', 'fail'}
            The insert mode:

            * 'replace' will create a new database table, overwriting an existing one.
            * 'append' will append to an existing table.
            * 'fail' will fail if table already exists.
        engine : {'sqlalchemy', 'adbc'}
            Select the engine to use for writing frame data; only necessary when
            supplying a URI string (defaults to 'sqlalchemy' if unset)
        engine_options
            Additional options to pass to the insert method associated with the engine
            specified by the option `engine`.

            * Setting `engine` to "sqlalchemy" currently inserts using Pandas' `to_sql`
              method (though this will eventually be phased out in favor of a native
              solution).
            * Setting `engine` to "adbc" inserts using the ADBC cursor's `adbc_ingest`
              method. Note that when passing an instantiated connection object, PyArrow
              is required for SQLite and Snowflake drivers.

        Examples
        --------
        Insert into a temporary table using a PostgreSQL URI and the ADBC engine:

        >>> df.write_database(
        ...     table_name="target_table",
        ...     connection="postgresql://user:pass@server:port/database",
        ...     engine="adbc",
        ...     engine_options={"temporary": True},
        ... )  # doctest: +SKIP

        Insert into a table using a `pyodbc` SQLAlchemy connection to SQL Server
        that was instantiated with "fast_executemany=True" to improve performance:

        >>> pyodbc_uri = (
        ...     "mssql+pyodbc://user:pass@server:1433/test?"
        ...     "driver=ODBC+Driver+18+for+SQL+Server"
        ... )
        >>> engine = create_engine(pyodbc_uri, fast_executemany=True)  # doctest: +SKIP
        >>> df.write_database(
        ...     table_name="target_table",
        ...     connection=engine,
        ... )  # doctest: +SKIP

        Returns
        -------
        int
            The number of rows affected, if the driver provides this information.
            Otherwise, returns -1.
        """
        if if_table_exists not in (valid_write_modes := get_args(DbWriteMode)):
            allowed = ", ".join(repr(m) for m in valid_write_modes)
            msg = f"write_database `if_table_exists` must be one of {{{allowed}}}, got {if_table_exists!r}"
            raise ValueError(msg)

        connection_module_root = type(connection).__module__.split(".", 1)[0]

        if engine is None:
            if isinstance(connection, str) or connection_module_root == "sqlalchemy":
                engine = "sqlalchemy"
            elif connection_module_root.startswith("adbc"):
                engine = "adbc"

        def unpack_table_name(name: str) -> tuple[str | None, str | None, str]:
            """Unpack optionally qualified table name to catalog/schema/table tuple."""
            from csv import reader as delimited_read

            components: list[str | None] = next(delimited_read([name], delimiter="."))  # type: ignore[arg-type]
            if len(components) > 3:
                msg = f"`table_name` appears to be invalid: '{name}'"
                raise ValueError(msg)
            catalog, schema, tbl = ([None] * (3 - len(components))) + components
            return catalog, schema, tbl  # type: ignore[return-value]

        if engine == "adbc":
            from polars.io.database._utils import (
                _get_adbc_module_name_from_uri,
                _import_optional_adbc_driver,
                _is_adbc_snowflake_conn,
                _open_adbc_connection,
            )

            conn, can_close_conn = (
                (_open_adbc_connection(connection), True)
                if isinstance(connection, str)
                else (connection, False)
            )

            driver_manager = import_optional("adbc_driver_manager")

            # base class for ADBC connections
            if not isinstance(conn, driver_manager.dbapi.Connection):
                msg = (
                    f"unrecognised connection type {qualified_type_name(connection)!r}"
                )
                raise TypeError(msg)

            driver_manager_str_version = getattr(driver_manager, "__version__", "0.0")
            driver_manager_version = parse_version(driver_manager_str_version)

            if if_table_exists == "fail":
                # if the table exists, 'create' will raise an error,
                # resulting in behaviour equivalent to 'fail'
                mode = "create"
            elif if_table_exists == "replace":
                if driver_manager_version < (0, 7):
                    msg = (
                        "`if_table_exists = 'replace'` requires ADBC version >= 0.7, "
                        f"found {driver_manager_str_version}"
                    )
                    raise ModuleUpgradeRequiredError(msg)
                mode = "replace"
            elif if_table_exists == "append":
                mode = "append"
            else:
                msg = (
                    f"unexpected value for `if_table_exists`: {if_table_exists!r}"
                    f"\n\nChoose one of {{'fail', 'replace', 'append'}}"
                )
                raise ValueError(msg)

            with (
                conn if can_close_conn else contextlib.nullcontext(),
                conn.cursor() as cursor,
            ):
                catalog, db_schema, unpacked_table_name = unpack_table_name(table_name)
                n_rows: int

                # We can reliably introspect the underlying driver from a URI
                # We can also introspect instantiated connections when PyArrow is
                # installed. Otherwise, the underlying driver is unknown
                # Ref: https://github.com/apache/arrow-adbc/issues/2828
                if isinstance(connection, str):
                    adbc_module_name = _get_adbc_module_name_from_uri(connection)
                elif _PYARROW_AVAILABLE:
                    adbc_module_name = (
                        f"adbc_driver_{conn.adbc_get_info()['vendor_name'].lower()}"
                    )
                else:
                    adbc_module_name = "Unknown"

                if adbc_module_name != "Unknown":
                    adbc_driver = _import_optional_adbc_driver(
                        adbc_module_name, dbapi_submodule=False
                    )
                    adbc_driver_str_version = getattr(adbc_driver, "__version__", "0.0")
                else:
                    adbc_driver = "Unknown"
                    # If we can't introspect the driver, guess that it has the same
                    # version as the driver manager. This is what happens by default
                    # when installed
                    adbc_driver_str_version = driver_manager_str_version

                adbc_driver_version = parse_version(adbc_driver_str_version)

                if adbc_module_name.split("_")[-1] == "sqlite":
                    catalog, db_schema = db_schema, None

                    # note: ADBC didnt't support 'replace' until adbc-driver-sqlite
                    # version 0.11 (it was released for other drivers in version 0.7)
                    if (
                        driver_manager_version >= (0, 7)
                        and adbc_driver_version < (0, 11)
                        and if_table_exists == "replace"
                    ):
                        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
                        mode = "create"

                # For Snowflake, we convert to PyArrow until string_view columns can be
                # written. Ref: https://github.com/apache/arrow-adbc/issues/3420
                is_snowflake_driver = (
                    "snowflake" in adbc_module_name
                    if _PYARROW_AVAILABLE
                    else _is_adbc_snowflake_conn(conn)
                )
                if is_snowflake_driver and not _PYARROW_AVAILABLE:
                    msg = (
                        "write_database with Snowflake driver requires 'pyarrow'.\n"
                        "Please install using the command `pip install pyarrow`."
                    )
                    raise ModuleNotFoundError(msg)

                # As of adbc_driver_manager 1.6.0, adbc_ingest can take a Polars
                # DataFrame via the PyCapsule interface
                data = (
                    self
                    if (driver_manager_version >= (1, 6)) and not is_snowflake_driver
                    else self.to_arrow()
                )

                # use of schema-qualified table names was released in
                # adbc-driver-manager 0.7.0 and is working without bugs from driver
                # version (e.g., adbc-driver-postgresql) version 0.8.0
                if driver_manager_version >= (0, 7) and adbc_driver_version >= (0, 8):
                    n_rows = cursor.adbc_ingest(
                        unpacked_table_name,
                        data=data,
                        mode=mode,
                        catalog_name=catalog,
                        db_schema_name=db_schema,
                        **(engine_options or {}),
                    )
                elif db_schema is not None:
                    adbc_driver_pypi_name = (
                        adbc_module_name.replace("_", "-")
                        if adbc_module_name != "Unknown"
                        else "adbc-driver-<driver>"
                    )
                    msg = (
                        "use of schema-qualified table names requires "
                        "adbc-driver-manager version >= 0.7.0, found "
                        f"{driver_manager_str_version} and {adbc_driver_pypi_name} "
                        f"version >= 0.8.0, found {adbc_driver_str_version}"
                    )
                    raise ModuleUpgradeRequiredError(
                        # https://github.com/apache/arrow-adbc/issues/1000
                        # https://github.com/apache/arrow-adbc/issues/1109
                        msg
                    )
                else:
                    n_rows = cursor.adbc_ingest(
                        table_name=unpacked_table_name,
                        data=data,
                        mode=mode,
                        **(engine_options or {}),
                    )
                conn.commit()
            return n_rows

        elif engine == "sqlalchemy":
            if not _PANDAS_AVAILABLE:
                msg = "writing with 'sqlalchemy' engine currently requires pandas.\n\nInstall with: pip install pandas"
                raise ModuleNotFoundError(msg)
            elif (pd_version := parse_version(pd.__version__)) < (1, 5):
                msg = f"writing with 'sqlalchemy' engine requires pandas >= 1.5; found {pd.__version__!r}"
                raise ModuleUpgradeRequiredError(msg)

            import_optional(
                module_name="sqlalchemy",
                min_version=("2.0" if pd_version >= (2, 2) else "1.4"),
                min_err_prefix="pandas >= 2.2 requires",
            )
            # note: the catalog (database) should be a part of the connection string
            from sqlalchemy.engine import Connectable, create_engine
            from sqlalchemy.orm import Session

            sa_object: Connectable
            if isinstance(connection, str):
                sa_object = create_engine(connection)
            elif isinstance(connection, Session):
                sa_object = connection.connection()
            elif isinstance(connection, Connectable):
                sa_object = connection
            else:
                msg = (
                    f"unrecognised connection type {qualified_type_name(connection)!r}"
                )
                raise TypeError(msg)

            catalog, db_schema, unpacked_table_name = unpack_table_name(table_name)
            if catalog:
                msg = f"Unexpected three-part table name; provide the database/catalog ({catalog!r}) on the connection URI"
                raise ValueError(msg)

            # ensure conversion to pandas uses the pyarrow extension array option
            # so that we can make use of the sql/db export *without* copying data
            res: int | None = self.to_pandas(
                use_pyarrow_extension_array=True,
            ).to_sql(
                name=unpacked_table_name,
                schema=db_schema,
                con=sa_object,
                if_exists=if_table_exists,
                index=False,
                **(engine_options or {}),
            )
            return -1 if res is None else res

        elif isinstance(engine, str):
            msg = f"engine {engine!r} is not supported"
            raise ValueError(msg)
        else:
            msg = f"unrecognised connection type {qualified_type_name(connection)!r}"
            raise TypeError(msg)

    @unstable()
    def write_iceberg(
        self,
        target: str | pyiceberg.table.Table,
        mode: Literal["append", "overwrite"],
    ) -> None:
        """
        Write DataFrame to an Iceberg table.

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Parameters
        ----------
        target
            Name of the table or the Table object representing an Iceberg table.
        mode : {'append', 'overwrite'}
            How to handle existing data.

            - If 'append', will add new data.
            - If 'overwrite', will replace table with new data.

        """
        from pyiceberg.catalog import load_catalog

        if isinstance(target, str):
            catalog = load_catalog()
            table = catalog.load_table(target)
        else:
            table = target

        data = self.to_arrow(compat_level=CompatLevel.oldest())

        if mode == "append":
            table.append(data)
        else:
            table.overwrite(data)

    @overload
    def write_delta(
        self,
        target: str | Path | deltalake.DeltaTable,
        *,
        mode: Literal["error", "append", "overwrite", "ignore"] = ...,
        overwrite_schema: bool | None = ...,
        storage_options: dict[str, str] | None = ...,
        credential_provider: CredentialProviderFunction | Literal["auto"] | None = ...,
        delta_write_options: dict[str, Any] | None = ...,
    ) -> None: ...

    @overload
    def write_delta(
        self,
        target: str | Path | deltalake.DeltaTable,
        *,
        mode: Literal["merge"],
        overwrite_schema: bool | None = ...,
        storage_options: dict[str, str] | None = ...,
        credential_provider: CredentialProviderFunction | Literal["auto"] | None = ...,
        delta_merge_options: dict[str, Any],
    ) -> deltalake.table.TableMerger: ...

    def write_delta(
        self,
        target: str | Path | deltalake.DeltaTable,
        *,
        mode: Literal["error", "append", "overwrite", "ignore", "merge"] = "error",
        overwrite_schema: bool | None = None,
        storage_options: dict[str, str] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        delta_write_options: dict[str, Any] | None = None,
        delta_merge_options: dict[str, Any] | None = None,
    ) -> deltalake.table.TableMerger | None:
        """
        Write DataFrame as delta table.

        Parameters
        ----------
        target
            URI of a table or a DeltaTable object.
        mode : {'error', 'append', 'overwrite', 'ignore', 'merge'}
            How to handle existing data.

            - If 'error', throw an error if the table already exists (default).
            - If 'append', will add new data.
            - If 'overwrite', will replace table with new data.
            - If 'ignore', will not write anything if table already exists.
            - If 'merge', return a `TableMerger` object to merge data from the DataFrame
              with the existing data.
        overwrite_schema
            If True, allows updating the schema of the table.

            .. deprecated:: 0.20.14
                Use the parameter `delta_write_options` instead and pass
                `{"schema_mode": "overwrite"}`.
        storage_options
            Extra options for the storage backends supported by `deltalake`.
            For cloud storages, this may include configurations for authentication etc.

            - See a list of supported storage options for S3 `here <https://docs.rs/object_store/latest/object_store/aws/enum.AmazonS3ConfigKey.html#variants>`__.
            - See a list of supported storage options for GCS `here <https://docs.rs/object_store/latest/object_store/gcp/enum.GoogleConfigKey.html#variants>`__.
            - See a list of supported storage options for Azure `here <https://docs.rs/object_store/latest/object_store/azure/enum.AzureConfigKey.html#variants>`__.
        credential_provider
            Provide a function that can be called to provide cloud storage
            credentials. The function is expected to return a dictionary of
            credential keys along with an optional credential expiry time.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        delta_write_options
            Additional keyword arguments while writing a Delta lake Table.
            See a list of supported write options `here <https://delta-io.github.io/delta-rs/api/delta_writer/#deltalake.write_deltalake>`__.
        delta_merge_options
            Keyword arguments which are required to `MERGE` a Delta lake Table.
            See a list of supported merge options `here <https://delta-io.github.io/delta-rs/api/delta_table/#deltalake.DeltaTable.merge>`__.

        Raises
        ------
        TypeError
            If the DataFrame contains unsupported data types.
        ArrowInvalidError
            If the DataFrame contains data types that could not be cast to their
            primitive type.
        TableNotFoundError
            If the delta table doesn't exist and MERGE action is triggered

        Notes
        -----
        The Polars data types :class:`Null` and :class:`Time` are not supported
        by the delta protocol specification and will raise a TypeError. Columns
        using The :class:`Categorical` data type will be converted to
        normal (non-categorical) strings when written.

        Polars columns are always nullable. To write data to a delta table with
        non-nullable columns, a custom pyarrow schema has to be passed to the
        `delta_write_options`. See the last example below.

        Examples
        --------
        Write a dataframe to the local filesystem as a Delta Lake table.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> table_path = "/path/to/delta-table/"
        >>> df.write_delta(table_path)  # doctest: +SKIP

        Append data to an existing Delta Lake table on the local filesystem.
        Note that this will fail if the schema of the new data does not match the
        schema of the existing table.

        >>> df.write_delta(table_path, mode="append")  # doctest: +SKIP

        Overwrite a Delta Lake table as a new version.
        If the schemas of the new and old data are the same, specifying the
        `schema_mode` is not required.

        >>> existing_table_path = "/path/to/delta-table/"
        >>> df.write_delta(
        ...     existing_table_path,
        ...     mode="overwrite",
        ...     delta_write_options={"schema_mode": "overwrite"},
        ... )  # doctest: +SKIP

        Write a DataFrame as a Delta Lake table to a cloud object store like S3.

        >>> table_path = "s3://bucket/prefix/to/delta-table/"
        >>> df.write_delta(
        ...     table_path,
        ...     storage_options={
        ...         "AWS_REGION": "THE_AWS_REGION",
        ...         "AWS_ACCESS_KEY_ID": "THE_AWS_ACCESS_KEY_ID",
        ...         "AWS_SECRET_ACCESS_KEY": "THE_AWS_SECRET_ACCESS_KEY",
        ...     },
        ... )  # doctest: +SKIP

        Write DataFrame as a Delta Lake table with non-nullable columns.

        >>> import pyarrow as pa
        >>> existing_table_path = "/path/to/delta-table/"
        >>> df.write_delta(
        ...     existing_table_path,
        ...     delta_write_options={
        ...         "schema": pa.schema([pa.field("foo", pa.int64(), nullable=False)])
        ...     },
        ... )  # doctest: +SKIP

        Write DataFrame as a Delta Lake table with zstd compression.
        For all `delta_write_options` keyword arguments, check the deltalake docs
        `here
        <https://delta-io.github.io/delta-rs/api/delta_writer/#deltalake.write_deltalake>`__,
        and for Writer Properties in particular `here
        <https://delta-io.github.io/delta-rs/api/delta_writer/#deltalake.WriterProperties>`__.

        >>> import deltalake
        >>> df.write_delta(
        ...     table_path,
        ...     delta_write_options={
        ...         "writer_properties": deltalake.WriterProperties(compression="zstd"),
        ...     },
        ... )  # doctest: +SKIP

        Merge the DataFrame with an existing Delta Lake table.
        For all `TableMerger` methods, check the deltalake docs
        `here <https://delta-io.github.io/delta-rs/api/delta_table/delta_table_merger/>`__.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> table_path = "/path/to/delta-table/"
        >>> (
        ...     df.write_delta(
        ...         "table_path",
        ...         mode="merge",
        ...         delta_merge_options={
        ...             "predicate": "s.foo = t.foo",
        ...             "source_alias": "s",
        ...             "target_alias": "t",
        ...         },
        ...     )
        ...     .when_matched_update_all()
        ...     .when_not_matched_insert_all()
        ...     .execute()
        ... )  # doctest: +SKIP
        """
        if overwrite_schema is not None:
            issue_deprecation_warning(
                "the parameter `overwrite_schema` for `write_delta` is deprecated."
                ' Use the parameter `delta_write_options` instead and pass `{"schema_mode": "overwrite"}`.',
                version="0.20.14",
            )

        from polars.io.delta import (
            _check_for_unsupported_types,
            _check_if_delta_available,
            _resolve_delta_lake_uri,
        )

        _check_if_delta_available()

        from deltalake import DeltaTable, write_deltalake

        _check_for_unsupported_types(self.dtypes)

        if isinstance(target, (str, Path)):
            target = _resolve_delta_lake_uri(str(target), strict=False)

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )
        from polars.io.cloud.credential_provider._providers import (
            _get_credentials_from_provider_expiry_aware,
        )

        if not isinstance(target, DeltaTable):
            credential_provider_builder = _init_credential_provider_builder(
                credential_provider, target, storage_options, "write_delta"
            )
        elif credential_provider is not None and credential_provider != "auto":
            msg = "cannot use credential_provider when passing a DeltaTable object"
            raise ValueError(msg)
        else:
            credential_provider_builder = None

        del credential_provider

        credential_provider_creds = {}

        if credential_provider_builder and (
            provider := credential_provider_builder.build_credential_provider()
        ):
            credential_provider_creds = (
                _get_credentials_from_provider_expiry_aware(provider) or {}
            )

        # We aren't calling into polars-native write functions so we just update
        # the storage_options here.
        storage_options = (
            {**(storage_options or {}), **credential_provider_creds}
            if storage_options is not None or credential_provider_builder is not None
            else None
        )

        if mode == "merge":
            if delta_merge_options is None:
                msg = "you need to pass delta_merge_options with at least a given predicate for `MERGE` to work."
                raise ValueError(msg)
            if isinstance(target, str):
                dt = DeltaTable(table_uri=target, storage_options=storage_options)
            else:
                dt = target

            return dt.merge(self, **delta_merge_options)

        else:
            if delta_write_options is None:
                delta_write_options = {}

            if overwrite_schema:
                delta_write_options["schema_mode"] = "overwrite"

            write_deltalake(
                table_or_uri=target,
                data=self,
                mode=mode,
                storage_options=storage_options,
                **delta_write_options,
            )
            return None

    def estimated_size(self, unit: SizeUnit = "b") -> int | float:
        """
        Return an estimation of the total (heap) allocated size of the `DataFrame`.

        Estimated size is given in the specified unit (bytes by default).

        This estimation is the sum of the size of its buffers, validity, including
        nested arrays. Multiple arrays may share buffers and bitmaps. Therefore, the
        size of 2 arrays is not the sum of the sizes computed from this function. In
        particular, [`StructArray`]'s size is an upper bound.

        When an array is sliced, its allocated size remains constant because the buffer
        unchanged. However, this function will yield a smaller number. This is because
        this function returns the visible size of the buffer, not its total capacity.

        FFI buffers are included in this estimation.

        Notes
        -----
        For data with Object dtype, the estimated size only reports the pointer
        size, which is a huge underestimation.

        Parameters
        ----------
        unit : {'b', 'kb', 'mb', 'gb', 'tb'}
            Scale the returned size to the given unit.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "x": list(reversed(range(1_000_000))),
        ...         "y": [v / 1000 for v in range(1_000_000)],
        ...         "z": [str(v) for v in range(1_000_000)],
        ...     },
        ...     schema=[("x", pl.UInt32), ("y", pl.Float64), ("z", pl.String)],
        ... )
        >>> df.estimated_size()
        17888890
        >>> df.estimated_size("mb")
        17.0601749420166
        """
        sz = self._df.estimated_size()
        return scale_bytes(sz, unit)

    def transpose(
        self,
        *,
        include_header: bool = False,
        header_name: str = "column",
        column_names: str | Iterable[str] | None = None,
    ) -> DataFrame:
        """
        Transpose a DataFrame over the diagonal.

        Parameters
        ----------
        include_header
            If set, the column names will be added as first column.
        header_name
            If `include_header` is set, this determines the name of the column that will
            be inserted.
        column_names
            Optional iterable yielding strings or a string naming an existing column.
            These will name the value (non-header) columns in the transposed data.

        Notes
        -----
        This is a very expensive operation. Perhaps you can do it differently.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        >>> df.transpose(include_header=True)
        shape: (2, 4)
        ┌────────┬──────────┬──────────┬──────────┐
        │ column ┆ column_0 ┆ column_1 ┆ column_2 │
        │ ---    ┆ ---      ┆ ---      ┆ ---      │
        │ str    ┆ i64      ┆ i64      ┆ i64      │
        ╞════════╪══════════╪══════════╪══════════╡
        │ a      ┆ 1        ┆ 2        ┆ 3        │
        │ b      ┆ 4        ┆ 5        ┆ 6        │
        └────────┴──────────┴──────────┴──────────┘

        Replace the auto-generated column names with a list

        >>> df.transpose(include_header=False, column_names=["x", "y", "z"])
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ x   ┆ y   ┆ z   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 2   ┆ 3   │
        │ 4   ┆ 5   ┆ 6   │
        └─────┴─────┴─────┘

        Include the header as a separate column

        >>> df.transpose(
        ...     include_header=True, header_name="foo", column_names=["x", "y", "z"]
        ... )
        shape: (2, 4)
        ┌─────┬─────┬─────┬─────┐
        │ foo ┆ x   ┆ y   ┆ z   │
        │ --- ┆ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 2   ┆ 3   │
        │ b   ┆ 4   ┆ 5   ┆ 6   │
        └─────┴─────┴─────┴─────┘

        Replace the auto-generated column with column names from a generator function

        >>> def name_generator():
        ...     base_name = "my_column_"
        ...     count = 0
        ...     while True:
        ...         yield f"{base_name}{count}"
        ...         count += 1
        >>> df.transpose(include_header=False, column_names=name_generator())
        shape: (2, 3)
        ┌─────────────┬─────────────┬─────────────┐
        │ my_column_0 ┆ my_column_1 ┆ my_column_2 │
        │ ---         ┆ ---         ┆ ---         │
        │ i64         ┆ i64         ┆ i64         │
        ╞═════════════╪═════════════╪═════════════╡
        │ 1           ┆ 2           ┆ 3           │
        │ 4           ┆ 5           ┆ 6           │
        └─────────────┴─────────────┴─────────────┘

        Use an existing column as the new column names

        >>> df = pl.DataFrame(dict(id=["i", "j", "k"], a=[1, 2, 3], b=[4, 5, 6]))
        >>> df.transpose(column_names="id")
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ i   ┆ j   ┆ k   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 2   ┆ 3   │
        │ 4   ┆ 5   ┆ 6   │
        └─────┴─────┴─────┘
        >>> df.transpose(include_header=True, header_name="new_id", column_names="id")
        shape: (2, 4)
        ┌────────┬─────┬─────┬─────┐
        │ new_id ┆ i   ┆ j   ┆ k   │
        │ ---    ┆ --- ┆ --- ┆ --- │
        │ str    ┆ i64 ┆ i64 ┆ i64 │
        ╞════════╪═════╪═════╪═════╡
        │ a      ┆ 1   ┆ 2   ┆ 3   │
        │ b      ┆ 4   ┆ 5   ┆ 6   │
        └────────┴─────┴─────┴─────┘
        """
        keep_names_as = header_name if include_header else None
        column_names_: Sequence[str] | None
        if isinstance(column_names, Generator):
            column_names_ = [next(column_names) for _ in range(self.height)]
        else:
            column_names_ = column_names  # type: ignore[assignment]
        return self._from_pydf(self._df.transpose(keep_names_as, column_names_))

    def reverse(self) -> DataFrame:
        """
        Reverse the DataFrame.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "key": ["a", "b", "c"],
        ...         "val": [1, 2, 3],
        ...     }
        ... )
        >>> df.reverse()
        shape: (3, 2)
        ┌─────┬─────┐
        │ key ┆ val │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ c   ┆ 3   │
        │ b   ┆ 2   │
        │ a   ┆ 1   │
        └─────┴─────┘
        """
        return self.select(F.col("*").reverse())

    def rename(
        self, mapping: Mapping[str, str] | Callable[[str], str], *, strict: bool = True
    ) -> DataFrame:
        """
        Rename column names.

        Parameters
        ----------
        mapping
            Key value pairs that map from old name to new name, or a function
            that takes the old name as input and returns the new name.
        strict
            Validate that all column names exist in the current schema,
            and throw an exception if any do not. (Note that this parameter
            is a no-op when passing a function to `mapping`).

        See Also
        --------
        Expr.name.replace

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"foo": [1, 2, 3], "bar": [6, 7, 8], "ham": ["a", "b", "c"]}
        ... )
        >>> df.rename({"foo": "apple"})
        shape: (3, 3)
        ┌───────┬─────┬─────┐
        │ apple ┆ bar ┆ ham │
        │ ---   ┆ --- ┆ --- │
        │ i64   ┆ i64 ┆ str │
        ╞═══════╪═════╪═════╡
        │ 1     ┆ 6   ┆ a   │
        │ 2     ┆ 7   ┆ b   │
        │ 3     ┆ 8   ┆ c   │
        └───────┴─────┴─────┘
        >>> df.rename(lambda column_name: "c" + column_name[1:])
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ coo ┆ car ┆ cam │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .rename(mapping, strict=strict)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def insert_column(self, index: int, column: IntoExprColumn) -> DataFrame:
        """
        Insert a Series (or expression) at a certain column index.

        This operation is in place.

        Parameters
        ----------
        index
            Index at which to insert the new column.
        column
            `Series` or expression to insert.

        Examples
        --------
        Insert a new Series column at the given index:

        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> s = pl.Series("baz", [97, 98, 99])
        >>> df.insert_column(1, s)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ baz ┆ bar │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 97  ┆ 4   │
        │ 2   ┆ 98  ┆ 5   │
        │ 3   ┆ 99  ┆ 6   │
        └─────┴─────┴─────┘

        Insert a new expression column at the given index:

        >>> df = pl.DataFrame(
        ...     {"a": [2, 4, 2], "b": [0.5, 4, 10], "c": ["xx", "yy", "zz"]}
        ... )
        >>> expr = (pl.col("b") / pl.col("a")).alias("b_div_a")
        >>> df.insert_column(2, expr)
        shape: (3, 4)
        ┌─────┬──────┬─────────┬─────┐
        │ a   ┆ b    ┆ b_div_a ┆ c   │
        │ --- ┆ ---  ┆ ---     ┆ --- │
        │ i64 ┆ f64  ┆ f64     ┆ str │
        ╞═════╪══════╪═════════╪═════╡
        │ 2   ┆ 0.5  ┆ 0.25    ┆ xx  │
        │ 4   ┆ 4.0  ┆ 1.0     ┆ yy  │
        │ 2   ┆ 10.0 ┆ 5.0     ┆ zz  │
        └─────┴──────┴─────────┴─────┘
        """
        if (original_index := index) < 0:
            index = self.width + index
            if index < 0:
                msg = f"column index {original_index} is out of range (frame has {self.width} columns)"
                raise IndexError(msg)
        elif index > self.width:
            msg = f"column index {original_index} is out of range (frame has {self.width} columns)"
            raise IndexError(msg)

        if isinstance(column, pl.Series):
            self._df.insert_column(index, column._s)
        else:
            if isinstance(column, str):
                column = F.col(column)
            if isinstance(column, pl.Expr):
                cols = self.columns
                cols.insert(index, column)  # type: ignore[arg-type]
                self._df = self.select(cols)._df
            else:
                msg = f"column must be a Series or Expr, got {column!r} (type={qualified_type_name(column)})"
                raise TypeError(msg)
        return self

    def filter(
        self,
        *predicates: (
            IntoExprColumn
            | Iterable[IntoExprColumn]
            | bool
            | list[bool]
            | np.ndarray[Any, Any]
        ),
        **constraints: Any,
    ) -> DataFrame:
        """
        Filter rows, retaining those that match the given predicate expression(s).

        The original order of the remaining rows is preserved.

        Only rows where the predicate resolves as True are retained; when the
        predicate result is False (or null), the row is discarded.

        Parameters
        ----------
        predicates
            Expression(s) that evaluate to a boolean Series.
        constraints
            Column filters; use `name = value` to filter columns by the supplied value.
            Each constraint will behave the same as `pl.col(name).eq(value)`, and
            be implicitly joined with the other filter conditions using `&`.

        Notes
        -----
        If you are transitioning from Pandas, and performing filter operations based on
        the comparison of two or more columns, please note that in Polars any comparison
        involving `null` values will result in a `null` result, *not* boolean True or
        False. As a result, these rows will not be retained. Ensure that null values
        are handled appropriately to avoid unexpected behaviour (see examples below).

        See Also
        --------
        remove

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, None, 4, None, 0],
        ...         "bar": [6, 7, 8, None, None, 9, 0],
        ...         "ham": ["a", "b", "c", None, "d", "e", "f"],
        ...     }
        ... )

        Filter rows matching a condition:

        >>> df.filter(pl.col("foo") > 1)
        shape: (3, 3)
        ┌─────┬──────┬─────┐
        │ foo ┆ bar  ┆ ham │
        │ --- ┆ ---  ┆ --- │
        │ i64 ┆ i64  ┆ str │
        ╞═════╪══════╪═════╡
        │ 2   ┆ 7    ┆ b   │
        │ 3   ┆ 8    ┆ c   │
        │ 4   ┆ null ┆ d   │
        └─────┴──────┴─────┘

        Filter on multiple conditions, combined with and/or operators:

        >>> df.filter(
        ...     (pl.col("foo") < 3) & (pl.col("ham") == "a"),
        ... )
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘

        >>> df.filter(
        ...     (pl.col("foo") == 1) | (pl.col("ham") == "c"),
        ... )
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘

        Provide multiple filters using `*args` syntax:

        >>> df.filter(
        ...     pl.col("foo") <= 2,
        ...     ~pl.col("ham").is_in(["b", "c"]),
        ... )
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 0   ┆ 0   ┆ f   │
        └─────┴─────┴─────┘

        Provide multiple filters using `**kwargs` syntax:

        >>> df.filter(foo=2, ham="b")
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 2   ┆ 7   ┆ b   │
        └─────┴─────┴─────┘

        Filter by comparing two columns against each other:

        >>> df.filter(
        ...     pl.col("foo") == pl.col("bar"),
        ... )
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 0   ┆ 0   ┆ f   │
        └─────┴─────┴─────┘

        >>> df.filter(
        ...     pl.col("foo") != pl.col("bar"),
        ... )
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘

        Notice how the row with `None` values is filtered out. In order to keep the
        same behavior as pandas, use:

        >>> df.filter(
        ...     pl.col("foo").ne_missing(pl.col("bar")),
        ... )
        shape: (5, 3)
        ┌──────┬──────┬─────┐
        │ foo  ┆ bar  ┆ ham │
        │ ---  ┆ ---  ┆ --- │
        │ i64  ┆ i64  ┆ str │
        ╞══════╪══════╪═════╡
        │ 1    ┆ 6    ┆ a   │
        │ 2    ┆ 7    ┆ b   │
        │ 3    ┆ 8    ┆ c   │
        │ 4    ┆ null ┆ d   │
        │ null ┆ 9    ┆ e   │
        └──────┴──────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .filter(*predicates, **constraints)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def remove(
        self,
        *predicates: (
            IntoExprColumn
            | Iterable[IntoExprColumn]
            | bool
            | list[bool]
            | np.ndarray[Any, Any]
        ),
        **constraints: Any,
    ) -> DataFrame:
        """
        Remove rows, dropping those that match the given predicate expression(s).

        The original order of the remaining rows is preserved.

        Rows where the filter predicate does not evaluate to True are retained
        (this includes rows where the predicate evaluates as `null`).

        Parameters
        ----------
        predicates
            Expression that evaluates to a boolean Series.
        constraints
            Column filters; use `name = value` to filter columns using the supplied
            value. Each constraint behaves the same as `pl.col(name).eq(value)`,
            and is implicitly joined with the other filter conditions using `&`.

        Notes
        -----
        If you are transitioning from Pandas, and performing filter operations based on
        the comparison of two or more columns, please note that in Polars any comparison
        involving `null` values will result in a `null` result, *not* boolean True or
        False. As a result, these rows will not be removed. Ensure that null values
        are handled appropriately to avoid unexpected behaviour (see examples below).

        See Also
        --------
        filter

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [2, 3, None, 4, 0],
        ...         "bar": [5, 6, None, None, 0],
        ...         "ham": ["a", "b", None, "c", "d"],
        ...     }
        ... )

        Remove rows matching a condition:

        >>> df.remove(pl.col("bar") >= 5)
        shape: (3, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ 4    ┆ null ┆ c    │
        │ 0    ┆ 0    ┆ d    │
        └──────┴──────┴──────┘

        Discard rows based on multiple conditions, combined with and/or operators:

        >>> df.remove(
        ...     (pl.col("foo") >= 0) & (pl.col("bar") >= 0),
        ... )
        shape: (2, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ 4    ┆ null ┆ c    │
        └──────┴──────┴──────┘

        >>> df.remove(
        ...     (pl.col("foo") >= 0) | (pl.col("bar") >= 0),
        ... )
        shape: (1, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        └──────┴──────┴──────┘

        Provide multiple constraints using `*args` syntax:

        >>> df.remove(
        ...     pl.col("ham").is_not_null(),
        ...     pl.col("bar") >= 0,
        ... )
        shape: (2, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ 4    ┆ null ┆ c    │
        └──────┴──────┴──────┘

        Provide constraints(s) using `**kwargs` syntax:

        >>> df.remove(foo=0, bar=0)
        shape: (4, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ 2    ┆ 5    ┆ a    │
        │ 3    ┆ 6    ┆ b    │
        │ null ┆ null ┆ null │
        │ 4    ┆ null ┆ c    │
        └──────┴──────┴──────┘

        Remove rows by comparing two columns against each other:

        >>> df.remove(
        ...     pl.col("foo").ne_missing(pl.col("bar")),
        ... )
        shape: (2, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ 0    ┆ 0    ┆ d    │
        └──────┴──────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .remove(*predicates, **constraints)
            .collect(optimizations=QueryOptFlags._eager())
        )

    @overload
    def glimpse(
        self,
        *,
        max_items_per_column: int = ...,
        max_colname_length: int = ...,
        return_as_string: Literal[False] = ...,
    ) -> None: ...

    @overload
    def glimpse(
        self,
        *,
        max_items_per_column: int = ...,
        max_colname_length: int = ...,
        return_as_string: Literal[True],
    ) -> str: ...

    @overload
    def glimpse(
        self,
        *,
        max_items_per_column: int = ...,
        max_colname_length: int = ...,
        return_as_string: bool,
    ) -> str | None: ...

    def glimpse(
        self,
        *,
        max_items_per_column: int = 10,
        max_colname_length: int = 50,
        return_as_string: bool = False,
    ) -> str | None:
        """
        Return a dense preview of the DataFrame.

        The formatting shows one line per column so that wide dataframes display
        cleanly. Each line shows the column name, the data type, and the first
        few values.

        Parameters
        ----------
        max_items_per_column
            Maximum number of items to show per column.
        max_colname_length
            Maximum length of the displayed column names; values that exceed this
            value are truncated with a trailing ellipsis.
        return_as_string
            If True, return the preview as a string instead of printing to stdout.

        See Also
        --------
        describe, head, tail

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1.0, 2.8, 3.0],
        ...         "b": [4, 5, None],
        ...         "c": [True, False, True],
        ...         "d": [None, "b", "c"],
        ...         "e": ["usd", "eur", None],
        ...         "f": [date(2020, 1, 1), date(2021, 1, 2), date(2022, 1, 1)],
        ...     }
        ... )
        >>> df.glimpse()
        Rows: 3
        Columns: 6
        $ a  <f64> 1.0, 2.8, 3.0
        $ b  <i64> 4, 5, None
        $ c <bool> True, False, True
        $ d  <str> None, 'b', 'c'
        $ e  <str> 'usd', 'eur', None
        $ f <date> 2020-01-01, 2021-01-02, 2022-01-01
        """
        # always print at most this number of values (mainly ensures that
        # we do not cast long arrays to strings, which would be slow)
        max_n_values = min(max_items_per_column, self.height)
        schema = self.schema

        def _parse_column(col_name: str, dtype: PolarsDataType) -> tuple[str, str, str]:
            fn = repr if schema[col_name] == String else str
            values = self[:max_n_values, col_name].to_list()
            val_str = ", ".join(fn(v) for v in values)
            if len(col_name) > max_colname_length:
                col_name = col_name[: (max_colname_length - 1)] + "…"
            return col_name, f"<{_dtype_str_repr(dtype)}>", val_str

        data = [_parse_column(s, dtype) for s, dtype in self.schema.items()]

        # determine column layout widths
        max_col_name = max((len(col_name) for col_name, _, _ in data))
        max_col_dtype = max((len(dtype_str) for _, dtype_str, _ in data))

        # print header
        output = StringIO()
        output.write(f"Rows: {self.height}\nColumns: {self.width}\n")

        # print individual columns: one row per column
        for col_name, dtype_str, val_str in data:
            output.write(
                f"$ {col_name:<{max_col_name}} {dtype_str:>{max_col_dtype}} {val_str}\n"
            )

        s = output.getvalue()
        if return_as_string:
            return s

        print(s, end=None)
        return None

    def describe(
        self,
        percentiles: Sequence[float] | float | None = (0.25, 0.50, 0.75),
        *,
        interpolation: QuantileMethod = "nearest",
    ) -> DataFrame:
        """
        Summary statistics for a DataFrame.

        Parameters
        ----------
        percentiles
            One or more percentiles to include in the summary statistics.
            All values must be in the range `[0, 1]`.

        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method used when calculating percentiles.

        Notes
        -----
        The median is included by default as the 50% percentile.

        Warnings
        --------
        We do not guarantee the output of `describe` to be stable. It will show
        statistics that we deem informative, and may be updated in the future.
        Using `describe` programmatically (versus interactive exploration) is
        not recommended for this reason.

        See Also
        --------
        glimpse

        Examples
        --------
        >>> from datetime import date, time
        >>> df = pl.DataFrame(
        ...     {
        ...         "float": [1.0, 2.8, 3.0],
        ...         "int": [40, 50, None],
        ...         "bool": [True, False, True],
        ...         "str": ["zz", "xx", "yy"],
        ...         "date": [date(2020, 1, 1), date(2021, 7, 5), date(2022, 12, 31)],
        ...         "time": [time(10, 20, 30), time(14, 45, 50), time(23, 15, 10)],
        ...     }
        ... )

        Show default frame statistics:

        >>> df.describe()
        shape: (9, 7)
        ┌────────────┬──────────┬──────────┬──────────┬──────┬─────────────────────┬──────────┐
        │ statistic  ┆ float    ┆ int      ┆ bool     ┆ str  ┆ date                ┆ time     │
        │ ---        ┆ ---      ┆ ---      ┆ ---      ┆ ---  ┆ ---                 ┆ ---      │
        │ str        ┆ f64      ┆ f64      ┆ f64      ┆ str  ┆ str                 ┆ str      │
        ╞════════════╪══════════╪══════════╪══════════╪══════╪═════════════════════╪══════════╡
        │ count      ┆ 3.0      ┆ 2.0      ┆ 3.0      ┆ 3    ┆ 3                   ┆ 3        │
        │ null_count ┆ 0.0      ┆ 1.0      ┆ 0.0      ┆ 0    ┆ 0                   ┆ 0        │
        │ mean       ┆ 2.266667 ┆ 45.0     ┆ 0.666667 ┆ null ┆ 2021-07-02 16:00:00 ┆ 16:07:10 │
        │ std        ┆ 1.101514 ┆ 7.071068 ┆ null     ┆ null ┆ null                ┆ null     │
        │ min        ┆ 1.0      ┆ 40.0     ┆ 0.0      ┆ xx   ┆ 2020-01-01          ┆ 10:20:30 │
        │ 25%        ┆ 2.8      ┆ 40.0     ┆ null     ┆ null ┆ 2021-07-05          ┆ 14:45:50 │
        │ 50%        ┆ 2.8      ┆ 50.0     ┆ null     ┆ null ┆ 2021-07-05          ┆ 14:45:50 │
        │ 75%        ┆ 3.0      ┆ 50.0     ┆ null     ┆ null ┆ 2022-12-31          ┆ 23:15:10 │
        │ max        ┆ 3.0      ┆ 50.0     ┆ 1.0      ┆ zz   ┆ 2022-12-31          ┆ 23:15:10 │
        └────────────┴──────────┴──────────┴──────────┴──────┴─────────────────────┴──────────┘

        Customize which percentiles are displayed, applying linear interpolation:

        >>> with pl.Config(tbl_rows=12):
        ...     df.describe(
        ...         percentiles=[0.1, 0.3, 0.5, 0.7, 0.9],
        ...         interpolation="linear",
        ...     )
        shape: (11, 7)
        ┌────────────┬──────────┬──────────┬──────────┬──────┬─────────────────────┬──────────┐
        │ statistic  ┆ float    ┆ int      ┆ bool     ┆ str  ┆ date                ┆ time     │
        │ ---        ┆ ---      ┆ ---      ┆ ---      ┆ ---  ┆ ---                 ┆ ---      │
        │ str        ┆ f64      ┆ f64      ┆ f64      ┆ str  ┆ str                 ┆ str      │
        ╞════════════╪══════════╪══════════╪══════════╪══════╪═════════════════════╪══════════╡
        │ count      ┆ 3.0      ┆ 2.0      ┆ 3.0      ┆ 3    ┆ 3                   ┆ 3        │
        │ null_count ┆ 0.0      ┆ 1.0      ┆ 0.0      ┆ 0    ┆ 0                   ┆ 0        │
        │ mean       ┆ 2.266667 ┆ 45.0     ┆ 0.666667 ┆ null ┆ 2021-07-02 16:00:00 ┆ 16:07:10 │
        │ std        ┆ 1.101514 ┆ 7.071068 ┆ null     ┆ null ┆ null                ┆ null     │
        │ min        ┆ 1.0      ┆ 40.0     ┆ 0.0      ┆ xx   ┆ 2020-01-01          ┆ 10:20:30 │
        │ 10%        ┆ 1.36     ┆ 41.0     ┆ null     ┆ null ┆ 2020-04-20          ┆ 11:13:34 │
        │ 30%        ┆ 2.08     ┆ 43.0     ┆ null     ┆ null ┆ 2020-11-26          ┆ 12:59:42 │
        │ 50%        ┆ 2.8      ┆ 45.0     ┆ null     ┆ null ┆ 2021-07-05          ┆ 14:45:50 │
        │ 70%        ┆ 2.88     ┆ 47.0     ┆ null     ┆ null ┆ 2022-02-07          ┆ 18:09:34 │
        │ 90%        ┆ 2.96     ┆ 49.0     ┆ null     ┆ null ┆ 2022-09-13          ┆ 21:33:18 │
        │ max        ┆ 3.0      ┆ 50.0     ┆ 1.0      ┆ zz   ┆ 2022-12-31          ┆ 23:15:10 │
        └────────────┴──────────┴──────────┴──────────┴──────┴─────────────────────┴──────────┘
        """  # noqa: W505
        if not self.columns:
            msg = "cannot describe a DataFrame that has no columns"
            raise TypeError(msg)

        return self.lazy().describe(
            percentiles=percentiles, interpolation=interpolation
        )

    def get_column_index(self, name: str) -> int:
        """
        Find the index of a column by name.

        Parameters
        ----------
        name
            Name of the column to find.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"foo": [1, 2, 3], "bar": [6, 7, 8], "ham": ["a", "b", "c"]}
        ... )
        >>> df.get_column_index("ham")
        2
        >>> df.get_column_index("sandwich")  # doctest: +SKIP
        ColumnNotFoundError: sandwich
        """
        return self._df.get_column_index(name)

    def replace_column(self, index: int, column: Series) -> DataFrame:
        """
        Replace a column at an index location.

        This operation is in place.

        Parameters
        ----------
        index
            Column index.
        column
            Series that will replace the column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> s = pl.Series("apple", [10, 20, 30])
        >>> df.replace_column(0, s)
        shape: (3, 3)
        ┌───────┬─────┬─────┐
        │ apple ┆ bar ┆ ham │
        │ ---   ┆ --- ┆ --- │
        │ i64   ┆ i64 ┆ str │
        ╞═══════╪═════╪═════╡
        │ 10    ┆ 6   ┆ a   │
        │ 20    ┆ 7   ┆ b   │
        │ 30    ┆ 8   ┆ c   │
        └───────┴─────┴─────┘
        """
        if index < 0:
            index = self.width + index
        self._df.replace_column(index, column._s)
        return self

    def sort(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        *more_by: IntoExpr,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool | Sequence[bool] = False,
        multithreaded: bool = True,
        maintain_order: bool = False,
    ) -> DataFrame:
        """
        Sort the dataframe by the given columns.

        Parameters
        ----------
        by
            Column(s) to sort by. Accepts expression input, including selectors. Strings
            are parsed as column names.
        *more_by
            Additional columns to sort by, specified as positional arguments.
        descending
            Sort in descending order. When sorting by multiple columns, can be specified
            per column by passing a sequence of booleans.
        nulls_last
            Place null values last; can specify a single boolean applying to all columns
            or a sequence of booleans for per-column control.
        multithreaded
            Sort using multiple threads.
        maintain_order
            Whether the order should be maintained if elements are equal.

        Examples
        --------
        Pass a single column name to sort by that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None],
        ...         "b": [6.0, 5.0, 4.0],
        ...         "c": ["a", "c", "b"],
        ...     }
        ... )
        >>> df.sort("a")
        shape: (3, 3)
        ┌──────┬─────┬─────┐
        │ a    ┆ b   ┆ c   │
        │ ---  ┆ --- ┆ --- │
        │ i64  ┆ f64 ┆ str │
        ╞══════╪═════╪═════╡
        │ null ┆ 4.0 ┆ b   │
        │ 1    ┆ 6.0 ┆ a   │
        │ 2    ┆ 5.0 ┆ c   │
        └──────┴─────┴─────┘

        Sorting by expressions is also supported.

        >>> df.sort(pl.col("a") + pl.col("b") * 2, nulls_last=True)
        shape: (3, 3)
        ┌──────┬─────┬─────┐
        │ a    ┆ b   ┆ c   │
        │ ---  ┆ --- ┆ --- │
        │ i64  ┆ f64 ┆ str │
        ╞══════╪═════╪═════╡
        │ 2    ┆ 5.0 ┆ c   │
        │ 1    ┆ 6.0 ┆ a   │
        │ null ┆ 4.0 ┆ b   │
        └──────┴─────┴─────┘

        Sort by multiple columns by passing a list of columns.

        >>> df.sort(["c", "a"], descending=True)
        shape: (3, 3)
        ┌──────┬─────┬─────┐
        │ a    ┆ b   ┆ c   │
        │ ---  ┆ --- ┆ --- │
        │ i64  ┆ f64 ┆ str │
        ╞══════╪═════╪═════╡
        │ 2    ┆ 5.0 ┆ c   │
        │ null ┆ 4.0 ┆ b   │
        │ 1    ┆ 6.0 ┆ a   │
        └──────┴─────┴─────┘

        Or use positional arguments to sort by multiple columns in the same way.

        >>> df.sort("c", "a", descending=[False, True])
        shape: (3, 3)
        ┌──────┬─────┬─────┐
        │ a    ┆ b   ┆ c   │
        │ ---  ┆ --- ┆ --- │
        │ i64  ┆ f64 ┆ str │
        ╞══════╪═════╪═════╡
        │ 1    ┆ 6.0 ┆ a   │
        │ null ┆ 4.0 ┆ b   │
        │ 2    ┆ 5.0 ┆ c   │
        └──────┴─────┴─────┘
        """
        from polars.lazyframe import QueryOptFlags

        return (
            self.lazy()
            .sort(
                by,
                *more_by,
                descending=descending,
                nulls_last=nulls_last,
                multithreaded=multithreaded,
                maintain_order=maintain_order,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    def sql(self, query: str, *, table_name: str = "self") -> DataFrame:
        """
        Execute a SQL query against the DataFrame.

        .. versionadded:: 0.20.24

        .. warning::
            This functionality is considered **unstable**, although it is close to
            being considered stable. It may be changed at any point without it being
            considered a breaking change.

        Parameters
        ----------
        query
            SQL query to execute.
        table_name
            Optionally provide an explicit name for the table that represents the
            calling frame (defaults to "self").

        Notes
        -----
        * The calling frame is automatically registered as a table in the SQL context
          under the name "self". If you want access to the DataFrames and LazyFrames
          found in the current globals, use the top-level :meth:`pl.sql <polars.sql>`.
        * More control over registration and execution behaviour is available by
          using the :class:`SQLContext` object.
        * The SQL query executes in lazy mode before being collected and returned
          as a DataFrame.

        See Also
        --------
        SQLContext

        Examples
        --------
        >>> from datetime import date
        >>> df1 = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["zz", "yy", "xx"],
        ...         "c": [date(1999, 12, 31), date(2010, 10, 10), date(2077, 8, 8)],
        ...     }
        ... )

        Query the DataFrame using SQL:

        >>> df1.sql("SELECT c, b FROM self WHERE a > 1")
        shape: (2, 2)
        ┌────────────┬─────┐
        │ c          ┆ b   │
        │ ---        ┆ --- │
        │ date       ┆ str │
        ╞════════════╪═════╡
        │ 2010-10-10 ┆ yy  │
        │ 2077-08-08 ┆ xx  │
        └────────────┴─────┘

        Apply transformations to a DataFrame using SQL, aliasing "self" to "frame".

        >>> df1.sql(
        ...     query='''
        ...         SELECT
        ...             a,
        ...             (a % 2 == 0) AS a_is_even,
        ...             CONCAT_WS(':', b, b) AS b_b,
        ...             EXTRACT(year FROM c) AS year,
        ...             0::float4 AS "zero",
        ...         FROM frame
        ...     ''',
        ...     table_name="frame",
        ... )
        shape: (3, 5)
        ┌─────┬───────────┬───────┬──────┬──────┐
        │ a   ┆ a_is_even ┆ b_b   ┆ year ┆ zero │
        │ --- ┆ ---       ┆ ---   ┆ ---  ┆ ---  │
        │ i64 ┆ bool      ┆ str   ┆ i32  ┆ f32  │
        ╞═════╪═══════════╪═══════╪══════╪══════╡
        │ 1   ┆ false     ┆ zz:zz ┆ 1999 ┆ 0.0  │
        │ 2   ┆ true      ┆ yy:yy ┆ 2010 ┆ 0.0  │
        │ 3   ┆ false     ┆ xx:xx ┆ 2077 ┆ 0.0  │
        └─────┴───────────┴───────┴──────┴──────┘
        """
        from polars.sql import SQLContext

        issue_unstable_warning(
            "`sql` is considered **unstable** (although it is close to being considered stable)."
        )
        with SQLContext(register_globals=False, eager=True) as ctx:
            name = table_name if table_name else "self"
            ctx.register(name=name, frame=self)
            return ctx.execute(query)

    @deprecate_renamed_parameter("descending", "reverse", version="1.0.0")
    def top_k(
        self,
        k: int,
        *,
        by: IntoExpr | Iterable[IntoExpr],
        reverse: bool | Sequence[bool] = False,
    ) -> DataFrame:
        """
        Return the `k` largest rows.

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        .. versionchanged:: 1.0.0
            The `descending` parameter was renamed `reverse`.

        Parameters
        ----------
        k
            Number of rows to return.
        by
            Column(s) used to determine the top rows.
            Accepts expression input. Strings are parsed as column names.
        reverse
            Consider the `k` smallest elements of the `by` column(s) (instead of the `k`
            largest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        bottom_k

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [2, 1, 1, 3, 2, 1],
        ...     }
        ... )

        Get the rows which contain the 4 largest values in column b.

        >>> df.top_k(4, by="b")
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ b   ┆ 3   │
        │ a   ┆ 2   │
        │ b   ┆ 2   │
        │ b   ┆ 1   │
        └─────┴─────┘

        Get the rows which contain the 4 largest values when sorting on column b and a.

        >>> df.top_k(4, by=["b", "a"])
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ b   ┆ 3   │
        │ b   ┆ 2   │
        │ a   ┆ 2   │
        │ c   ┆ 1   │
        └─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .top_k(k, by=by, reverse=reverse)
            .collect(
                optimizations=QueryOptFlags(
                    projection_pushdown=False,
                    predicate_pushdown=False,
                    comm_subplan_elim=False,
                    slice_pushdown=True,
                )
            )
        )

    @deprecate_renamed_parameter("descending", "reverse", version="1.0.0")
    def bottom_k(
        self,
        k: int,
        *,
        by: IntoExpr | Iterable[IntoExpr],
        reverse: bool | Sequence[bool] = False,
    ) -> DataFrame:
        """
        Return the `k` smallest rows.

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        .. versionchanged:: 1.0.0
            The `descending` parameter was renamed `reverse`.

        Parameters
        ----------
        k
            Number of rows to return.
        by
            Column(s) used to determine the bottom rows.
            Accepts expression input. Strings are parsed as column names.
        reverse
            Consider the `k` largest elements of the `by` column(s) (instead of the `k`
            smallest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        top_k

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [2, 1, 1, 3, 2, 1],
        ...     }
        ... )

        Get the rows which contain the 4 smallest values in column b.

        >>> df.bottom_k(4, by="b")
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ b   ┆ 1   │
        │ a   ┆ 1   │
        │ c   ┆ 1   │
        │ a   ┆ 2   │
        └─────┴─────┘

        Get the rows which contain the 4 smallest values when sorting on column a and b.

        >>> df.bottom_k(4, by=["a", "b"])
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ a   ┆ 1   │
        │ a   ┆ 2   │
        │ b   ┆ 1   │
        │ b   ┆ 2   │
        └─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .bottom_k(k, by=by, reverse=reverse)
            .collect(
                optimizations=QueryOptFlags(
                    projection_pushdown=False,
                    predicate_pushdown=False,
                    comm_subplan_elim=False,
                    slice_pushdown=True,
                )
            )
        )

    def equals(self, other: DataFrame, *, null_equal: bool = True) -> bool:
        """
        Check whether the DataFrame is equal to another DataFrame.

        Parameters
        ----------
        other
            DataFrame to compare with.
        null_equal
            Consider null values as equal.

        See Also
        --------
        polars.testing.assert_frame_equal

        Examples
        --------
        >>> df1 = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df2 = pl.DataFrame(
        ...     {
        ...         "foo": [3, 2, 1],
        ...         "bar": [8.0, 7.0, 6.0],
        ...         "ham": ["c", "b", "a"],
        ...     }
        ... )
        >>> df1.equals(df1)
        True
        >>> df1.equals(df2)
        False
        """
        require_same_type(self, other)
        return self._df.equals(other._df, null_equal=null_equal)

    def slice(self, offset: int, length: int | None = None) -> DataFrame:
        """
        Get a slice of this DataFrame.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None`, all rows starting at the offset
            will be selected.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.slice(1, 2)
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 2   ┆ 7.0 ┆ b   │
        │ 3   ┆ 8.0 ┆ c   │
        └─────┴─────┴─────┘
        """
        if (length is not None) and length < 0:
            length = self.height - offset + length
        return self._from_pydf(self._df.slice(offset, length))

    def head(self, n: int = 5) -> DataFrame:
        """
        Get the first `n` rows.

        Parameters
        ----------
        n
            Number of rows to return. If a negative value is passed, return all rows
            except the last `abs(n)`.

        See Also
        --------
        tail, glimpse, slice

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> df.head(3)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘

        Pass a negative value to get all rows `except` the last `abs(n)`.

        >>> df.head(-3)
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        └─────┴─────┴─────┘
        """
        if n < 0:
            n = max(0, self.height + n)
        return self._from_pydf(self._df.head(n))

    def tail(self, n: int = 5) -> DataFrame:
        """
        Get the last `n` rows.

        Parameters
        ----------
        n
            Number of rows to return. If a negative value is passed, return all rows
            except the first `abs(n)`.

        See Also
        --------
        head, slice

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> df.tail(3)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 3   ┆ 8   ┆ c   │
        │ 4   ┆ 9   ┆ d   │
        │ 5   ┆ 10  ┆ e   │
        └─────┴─────┴─────┘

        Pass a negative value to get all rows `except` the first `abs(n)`.

        >>> df.tail(-3)
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 4   ┆ 9   ┆ d   │
        │ 5   ┆ 10  ┆ e   │
        └─────┴─────┴─────┘
        """
        if n < 0:
            n = max(0, self.height + n)
        return self._from_pydf(self._df.tail(n))

    def limit(self, n: int = 5) -> DataFrame:
        """
        Get the first `n` rows.

        Alias for :func:`DataFrame.head`.

        Parameters
        ----------
        n
            Number of rows to return. If a negative value is passed, return all rows
            except the last `abs(n)`.

        See Also
        --------
        head

        Examples
        --------
        Get the first 3 rows of a DataFrame.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 4, 5],
        ...         "bar": [6, 7, 8, 9, 10],
        ...         "ham": ["a", "b", "c", "d", "e"],
        ...     }
        ... )
        >>> df.limit(3)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘
        """
        return self.head(n)

    def drop_nans(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
    ) -> DataFrame:
        """
        Drop all rows that contain one or more NaN values.

        The original order of the remaining rows is preserved.

        Parameters
        ----------
        subset
            Column name(s) for which NaN values are considered; if set to `None`
            (default), use all columns (note that only floating-point columns
            can contain NaNs).

        See Also
        --------
        drop_nulls

        Notes
        -----
        A NaN value is not the same as a null value.
        To drop null values, use :func:`drop_nulls`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [-20.5, float("nan"), 80.0],
        ...         "bar": [float("nan"), 110.0, 25.5],
        ...         "ham": ["xxx", "yyy", None],
        ...     }
        ... )

        The default behavior of this method is to drop rows where any single
        value in the row is NaN:

        >>> df.drop_nans()
        shape: (1, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ f64  ┆ f64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ 80.0 ┆ 25.5 ┆ null │
        └──────┴──────┴──────┘

        This behaviour can be constrained to consider only a subset of columns, as
        defined by name, or with a selector. For example, dropping rows only if
        there is a NaN in the "bar" column:

        >>> df.drop_nans(subset=["bar"])
        shape: (2, 3)
        ┌──────┬───────┬──────┐
        │ foo  ┆ bar   ┆ ham  │
        │ ---  ┆ ---   ┆ ---  │
        │ f64  ┆ f64   ┆ str  │
        ╞══════╪═══════╪══════╡
        │ NaN  ┆ 110.0 ┆ yyy  │
        │ 80.0 ┆ 25.5  ┆ null │
        └──────┴───────┴──────┘

        Dropping a row only if *all* values are NaN requires a different formulation:

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [float("nan"), float("nan"), float("nan"), float("nan")],
        ...         "b": [10.0, 2.5, float("nan"), 5.25],
        ...         "c": [65.75, float("nan"), float("nan"), 10.5],
        ...     }
        ... )
        >>> df.filter(~pl.all_horizontal(pl.all().is_nan()))
        shape: (3, 3)
        ┌─────┬──────┬───────┐
        │ a   ┆ b    ┆ c     │
        │ --- ┆ ---  ┆ ---   │
        │ f64 ┆ f64  ┆ f64   │
        ╞═════╪══════╪═══════╡
        │ NaN ┆ 10.0 ┆ 65.75 │
        │ NaN ┆ 2.5  ┆ NaN   │
        │ NaN ┆ 5.25 ┆ 10.5  │
        └─────┴──────┴───────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy().drop_nans(subset).collect(optimizations=QueryOptFlags._eager())
        )

    def drop_nulls(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
    ) -> DataFrame:
        """
        Drop all rows that contain one or more null values.

        The original order of the remaining rows is preserved.

        Parameters
        ----------
        subset
            Column name(s) for which null values are considered.
            If set to `None` (default), use all columns.

        See Also
        --------
        drop_nans

        Notes
        -----
        A null value is not the same as a NaN value.
        To drop NaN values, use :func:`drop_nans`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, None, 8],
        ...         "ham": ["a", "b", None],
        ...     }
        ... )

        The default behavior of this method is to drop rows where any single
        value of the row is null.

        >>> df.drop_nulls()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘

        This behaviour can be constrained to consider only a subset of columns, as
        defined by name or with a selector. For example, dropping rows if there is
        a null in any of the integer columns:

        >>> import polars.selectors as cs
        >>> df.drop_nulls(subset=cs.integer())
        shape: (2, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ i64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 1   ┆ 6   ┆ a    │
        │ 3   ┆ 8   ┆ null │
        └─────┴─────┴──────┘

        Below are some additional examples that show how to drop null
        values based on other conditions.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [None, None, None, None],
        ...         "b": [1, 2, None, 1],
        ...         "c": [1, None, None, 1],
        ...     }
        ... )
        >>> df
        shape: (4, 3)
        ┌──────┬──────┬──────┐
        │ a    ┆ b    ┆ c    │
        │ ---  ┆ ---  ┆ ---  │
        │ null ┆ i64  ┆ i64  │
        ╞══════╪══════╪══════╡
        │ null ┆ 1    ┆ 1    │
        │ null ┆ 2    ┆ null │
        │ null ┆ null ┆ null │
        │ null ┆ 1    ┆ 1    │
        └──────┴──────┴──────┘

        Drop a row only if all values are null:

        >>> df.filter(~pl.all_horizontal(pl.all().is_null()))
        shape: (3, 3)
        ┌──────┬─────┬──────┐
        │ a    ┆ b   ┆ c    │
        │ ---  ┆ --- ┆ ---  │
        │ null ┆ i64 ┆ i64  │
        ╞══════╪═════╪══════╡
        │ null ┆ 1   ┆ 1    │
        │ null ┆ 2   ┆ null │
        │ null ┆ 1   ┆ 1    │
        └──────┴─────┴──────┘

        Drop a column if all values are null:

        >>> df[[s.name for s in df if not (s.null_count() == df.height)]]
        shape: (4, 2)
        ┌──────┬──────┐
        │ b    ┆ c    │
        │ ---  ┆ ---  │
        │ i64  ┆ i64  │
        ╞══════╪══════╡
        │ 1    ┆ 1    │
        │ 2    ┆ null │
        │ null ┆ null │
        │ 1    ┆ 1    │
        └──────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy().drop_nulls(subset).collect(optimizations=QueryOptFlags._eager())
        )

    def pipe(
        self,
        function: Callable[Concatenate[DataFrame, P], T],
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> T:
        """
        Offers a structured way to apply a sequence of user-defined functions (UDFs).

        Parameters
        ----------
        function
            Callable; will receive the frame as the first parameter,
            followed by any given args/kwargs.
        *args
            Arguments to pass to the UDF.
        **kwargs
            Keyword arguments to pass to the UDF.

        Notes
        -----
        It is recommended to use LazyFrame when piping operations, in order
        to fully take advantage of query optimization and parallelization.
        See :meth:`df.lazy() <polars.DataFrame.lazy>`.

        Examples
        --------
        >>> def cast_str_to_int(data, col_name):
        ...     return data.with_columns(pl.col(col_name).cast(pl.Int64))
        >>> df = pl.DataFrame({"a": [1, 2, 3, 4], "b": ["10", "20", "30", "40"]})
        >>> df.pipe(cast_str_to_int, col_name="b")
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 10  │
        │ 2   ┆ 20  │
        │ 3   ┆ 30  │
        │ 4   ┆ 40  │
        └─────┴─────┘

        >>> df = pl.DataFrame({"b": [1, 2], "a": [3, 4]})
        >>> df
        shape: (2, 2)
        ┌─────┬─────┐
        │ b   ┆ a   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 3   │
        │ 2   ┆ 4   │
        └─────┴─────┘
        >>> df.pipe(lambda tdf: tdf.select(sorted(tdf.columns)))
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 3   ┆ 1   │
        │ 4   ┆ 2   │
        └─────┴─────┘
        """
        return function(self, *args, **kwargs)

    def map_columns(
        self,
        column_names: str | Sequence[str] | pl.Selector,
        function: Callable[[Series], Series],
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> DataFrame:
        """
        Apply eager functions to columns of a DataFrame.

        Users should always prefer :meth:`with_columns` unless they are using
        expressions that are only possible on `Series` and not on `Expr`. This is almost
        never the case, except for a very select few functions that cannot know the
        output datatype without looking at the data.

        Parameters
        ----------
        column_names
            The columns to apply the UDF to.
        function
            Callable; will receive a column series as the first parameter,
            followed by any given args/kwargs.
        *args
            Arguments to pass to the UDF.
        **kwargs
            Keyword arguments to pass to the UDF.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3, 4], "b": ["10", "20", "30", "40"]})
        >>> df.map_columns("a", lambda s: s.shrink_dtype())
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i8  ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ 10  │
        │ 2   ┆ 20  │
        │ 3   ┆ 30  │
        │ 4   ┆ 40  │
        └─────┴─────┘

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ['{"x":"a"}', None, '{"x":"b"}', None],
        ...         "b": ['{"a":1, "b": true}', None, '{"a":2, "b": false}', None],
        ...     }
        ... )
        >>> df.map_columns(["a", "b"], lambda s: s.str.json_decode())
        shape: (4, 2)
        ┌───────────┬───────────┐
        │ a         ┆ b         │
        │ ---       ┆ ---       │
        │ struct[1] ┆ struct[2] │
        ╞═══════════╪═══════════╡
        │ {"a"}     ┆ {1,true}  │
        │ null      ┆ null      │
        │ {"b"}     ┆ {2,false} │
        │ null      ┆ null      │
        └───────────┴───────────┘
        >>> import polars.selectors as cs
        >>> df.map_columns(cs.all(), lambda s: s.str.json_decode())
        shape: (4, 2)
        ┌───────────┬───────────┐
        │ a         ┆ b         │
        │ ---       ┆ ---       │
        │ struct[1] ┆ struct[2] │
        ╞═══════════╪═══════════╡
        │ {"a"}     ┆ {1,true}  │
        │ null      ┆ null      │
        │ {"b"}     ┆ {2,false} │
        │ null      ┆ null      │
        └───────────┴───────────┘

        See Also
        --------
        with_columns
        """
        c_names: list[str]
        if isinstance(column_names, (pl.Selector, pl.Expr)):
            from polars.selectors import expand_selector

            c_names = list(expand_selector(self, column_names))
        elif isinstance(column_names, str):
            c_names = [column_names]
        else:
            c_names = list(column_names)

        return self.with_columns(
            **{c: function(self[c], *args, **kwargs) for c in c_names}
        )

    def with_row_index(self, name: str = "index", offset: int = 0) -> DataFrame:
        """
        Add a row index as the first column in the DataFrame.

        Parameters
        ----------
        name
            Name of the index column.
        offset
            Start the index at this offset. Cannot be negative.

        Notes
        -----
        The resulting column does not have any special properties. It is a regular
        column of type `UInt32` (or `UInt64` in `polars[rt64]`).

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> df.with_row_index()
        shape: (3, 3)
        ┌───────┬─────┬─────┐
        │ index ┆ a   ┆ b   │
        │ ---   ┆ --- ┆ --- │
        │ u32   ┆ i64 ┆ i64 │
        ╞═══════╪═════╪═════╡
        │ 0     ┆ 1   ┆ 2   │
        │ 1     ┆ 3   ┆ 4   │
        │ 2     ┆ 5   ┆ 6   │
        └───────┴─────┴─────┘
        >>> df.with_row_index("id", offset=1000)
        shape: (3, 3)
        ┌──────┬─────┬─────┐
        │ id   ┆ a   ┆ b   │
        │ ---  ┆ --- ┆ --- │
        │ u32  ┆ i64 ┆ i64 │
        ╞══════╪═════╪═════╡
        │ 1000 ┆ 1   ┆ 2   │
        │ 1001 ┆ 3   ┆ 4   │
        │ 1002 ┆ 5   ┆ 6   │
        └──────┴─────┴─────┘

        An index column can also be created using the expressions :func:`int_range`
        and :func:`len`.

        >>> df.select(
        ...     pl.int_range(pl.len(), dtype=pl.UInt32).alias("index"),
        ...     pl.all(),
        ... )
        shape: (3, 3)
        ┌───────┬─────┬─────┐
        │ index ┆ a   ┆ b   │
        │ ---   ┆ --- ┆ --- │
        │ u32   ┆ i64 ┆ i64 │
        ╞═══════╪═════╪═════╡
        │ 0     ┆ 1   ┆ 2   │
        │ 1     ┆ 3   ┆ 4   │
        │ 2     ┆ 5   ┆ 6   │
        └───────┴─────┴─────┘
        """
        try:
            return self._from_pydf(self._df.with_row_index(name, offset))
        except OverflowError:
            issue = "negative" if offset < 0 else "greater than the maximum index value"
            msg = f"`offset` input for `with_row_index` cannot be {issue}, got {offset}"
            raise ValueError(msg) from None

    @deprecated(
        "`DataFrame.with_row_count` is deprecated; use `with_row_index` instead."
        " Note that the default column name has changed from 'row_nr' to 'index'."
    )
    def with_row_count(self, name: str = "row_nr", offset: int = 0) -> DataFrame:
        """
        Add a column at index 0 that counts the rows.

        .. deprecated:: 0.20.4
            Use the :meth:`with_row_index` method instead.
            Note that the default column name has changed from 'row_nr' to 'index'.

        Parameters
        ----------
        name
            Name of the column to add.
        offset
            Start the row count at this offset. Default = 0

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> df.with_row_count()  # doctest: +SKIP
        shape: (3, 3)
        ┌────────┬─────┬─────┐
        │ row_nr ┆ a   ┆ b   │
        │ ---    ┆ --- ┆ --- │
        │ u32    ┆ i64 ┆ i64 │
        ╞════════╪═════╪═════╡
        │ 0      ┆ 1   ┆ 2   │
        │ 1      ┆ 3   ┆ 4   │
        │ 2      ┆ 5   ┆ 6   │
        └────────┴─────┴─────┘
        """
        return self.with_row_index(name, offset)

    def group_by(
        self,
        *by: IntoExpr | Iterable[IntoExpr],
        maintain_order: bool = False,
        **named_by: IntoExpr,
    ) -> GroupBy:
        """
        Start a group by operation.

        Parameters
        ----------
        *by
            Column(s) to group by. Accepts expression input. Strings are parsed as
            column names.
        maintain_order
            Ensure that the order of the groups is consistent with the input data.
            This is slower than a default group by.
            Settings this to `True` blocks the possibility
            to run on the streaming engine.

            .. note::
                Within each group, the order of rows is always preserved, regardless
                of this argument.
        **named_by
            Additional columns to group by, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        Returns
        -------
        GroupBy
            Object which can be used to perform aggregations.

        Examples
        --------
        Group by one column and call `agg` to compute the grouped sum of another
        column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "c"],
        ...         "b": [1, 2, 1, 3, 3],
        ...         "c": [5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> df.group_by("a").agg(pl.col("b").sum())  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ a   ┆ 2   │
        │ b   ┆ 5   │
        │ c   ┆ 3   │
        └─────┴─────┘

        Set `maintain_order=True` to ensure the order of the groups is consistent with
        the input.

        >>> df.group_by("a", maintain_order=True).agg(pl.col("c"))
        shape: (3, 2)
        ┌─────┬───────────┐
        │ a   ┆ c         │
        │ --- ┆ ---       │
        │ str ┆ list[i64] │
        ╞═════╪═══════════╡
        │ a   ┆ [5, 3]    │
        │ b   ┆ [4, 2]    │
        │ c   ┆ [1]       │
        └─────┴───────────┘

        Group by multiple columns by passing a list of column names.

        >>> df.group_by(["a", "b"]).agg(pl.max("c"))  # doctest: +IGNORE_RESULT
        shape: (4, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 5   │
        │ b   ┆ 2   ┆ 4   │
        │ b   ┆ 3   ┆ 2   │
        │ c   ┆ 3   ┆ 1   │
        └─────┴─────┴─────┘

        Or use positional arguments to group by multiple columns in the same way.
        Expressions are also accepted.

        >>> df.group_by("a", pl.col("b") // 2).agg(pl.col("c").mean())  # doctest: +SKIP
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ f64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 0   ┆ 4.0 │
        │ b   ┆ 1   ┆ 3.0 │
        │ c   ┆ 1   ┆ 1.0 │
        └─────┴─────┴─────┘

        The `GroupBy` object returned by this method is iterable, returning the name
        and data of each group.

        >>> for name, data in df.group_by("a"):  # doctest: +SKIP
        ...     print(name)
        ...     print(data)
        ('a',)
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 5   │
        │ a   ┆ 1   ┆ 3   │
        └─────┴─────┴─────┘
        ('b',)
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 2   ┆ 4   │
        │ b   ┆ 3   ┆ 2   │
        └─────┴─────┴─────┘
        ('c',)
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ c   ┆ 3   ┆ 1   │
        └─────┴─────┴─────┘
        """
        for value in named_by.values():
            if not isinstance(value, (str, pl.Expr, pl.Series)):
                msg = (
                    f"Expected Polars expression or object convertible to one, got {type(value)}.\n\n"
                    "Hint: if you tried\n"
                    f"    group_by(by={value!r})\n"
                    "then you probably want to use this instead:\n"
                    f"    group_by({value!r})"
                )
                raise TypeError(msg)
        return GroupBy(self, *by, **named_by, maintain_order=maintain_order)

    @deprecate_renamed_parameter("by", "group_by", version="0.20.14")
    def rolling(
        self,
        index_column: IntoExpr,
        *,
        period: str | timedelta,
        offset: str | timedelta | None = None,
        closed: ClosedInterval = "right",
        group_by: IntoExpr | Iterable[IntoExpr] | None = None,
    ) -> RollingGroupBy:
        """
        Create rolling groups based on a temporal or integer column.

        Different from a `group_by_dynamic` the windows are now determined by the
        individual values and are not of constant intervals. For constant intervals use
        :func:`DataFrame.group_by_dynamic`.

        If you have a time series `<t_0, t_1, ..., t_n>`, then by default the
        windows created will be

            * (t_0 - period, t_0]
            * (t_1 - period, t_1]
            * ...
            * (t_n - period, t_n]

        whereas if you pass a non-default `offset`, then the windows will be

            * (t_0 + offset, t_0 + offset + period]
            * (t_1 + offset, t_1 + offset + period]
            * ...
            * (t_n + offset, t_n + offset + period]

        The `period` and `offset` arguments are created either from a timedelta, or
        by using the following string language:

        - 1ns   (1 nanosecond)
        - 1us   (1 microsecond)
        - 1ms   (1 millisecond)
        - 1s    (1 second)
        - 1m    (1 minute)
        - 1h    (1 hour)
        - 1d    (1 calendar day)
        - 1w    (1 calendar week)
        - 1mo   (1 calendar month)
        - 1q    (1 calendar quarter)
        - 1y    (1 calendar year)
        - 1i    (1 index count)

        Or combine them:
        "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

        By "calendar day", we mean the corresponding time on the next day (which may
        not be 24 hours, due to daylight savings). Similarly for "calendar week",
        "calendar month", "calendar quarter", and "calendar year".

        .. versionchanged:: 0.20.14
            The `by` parameter was renamed `group_by`.

        Parameters
        ----------
        index_column
            Column used to group based on the time window.
            Often of type Date/Datetime.
            This column must be sorted in ascending order (or, if `group_by` is
            specified, then it must be sorted in ascending order within each group).

            In case of a rolling operation on indices, dtype needs to be one of
            {UInt32, UInt64, Int32, Int64}. Note that the first three get temporarily
            cast to Int64, so if performance matters use an Int64 column.
        period
            Length of the window - must be non-negative.
        offset
            Offset of the window. Default is `-period`.
        closed : {'right', 'left', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive).
        group_by
            Also group by this column/these columns

        Returns
        -------
        RollingGroupBy
            Object you can call `.agg` on to aggregate by groups, the result
            of which will be sorted by `index_column` (but note that if `group_by`
            columns are passed, it will only be sorted within each group).

        See Also
        --------
        group_by_dynamic

        Examples
        --------
        >>> dates = [
        ...     "2020-01-01 13:45:48",
        ...     "2020-01-01 16:42:13",
        ...     "2020-01-01 16:45:09",
        ...     "2020-01-02 18:12:48",
        ...     "2020-01-03 19:45:32",
        ...     "2020-01-08 23:16:43",
        ... ]
        >>> df = pl.DataFrame({"dt": dates, "a": [3, 7, 5, 9, 2, 1]}).with_columns(
        ...     pl.col("dt").str.strptime(pl.Datetime).set_sorted()
        ... )
        >>> out = df.rolling(index_column="dt", period="2d").agg(
        ...     [
        ...         pl.sum("a").alias("sum_a"),
        ...         pl.min("a").alias("min_a"),
        ...         pl.max("a").alias("max_a"),
        ...     ]
        ... )
        >>> assert out["sum_a"].to_list() == [3, 10, 15, 24, 11, 1]
        >>> assert out["max_a"].to_list() == [3, 7, 7, 9, 9, 1]
        >>> assert out["min_a"].to_list() == [3, 3, 3, 3, 2, 1]
        >>> out
        shape: (6, 4)
        ┌─────────────────────┬───────┬───────┬───────┐
        │ dt                  ┆ sum_a ┆ min_a ┆ max_a │
        │ ---                 ┆ ---   ┆ ---   ┆ ---   │
        │ datetime[μs]        ┆ i64   ┆ i64   ┆ i64   │
        ╞═════════════════════╪═══════╪═══════╪═══════╡
        │ 2020-01-01 13:45:48 ┆ 3     ┆ 3     ┆ 3     │
        │ 2020-01-01 16:42:13 ┆ 10    ┆ 3     ┆ 7     │
        │ 2020-01-01 16:45:09 ┆ 15    ┆ 3     ┆ 7     │
        │ 2020-01-02 18:12:48 ┆ 24    ┆ 3     ┆ 9     │
        │ 2020-01-03 19:45:32 ┆ 11    ┆ 2     ┆ 9     │
        │ 2020-01-08 23:16:43 ┆ 1     ┆ 1     ┆ 1     │
        └─────────────────────┴───────┴───────┴───────┘

        If you use an index count in `period` or `offset`, then it's based on the
        values in `index_column`:

        >>> df = pl.DataFrame({"int": [0, 4, 5, 6, 8], "value": [1, 4, 2, 4, 1]})
        >>> df.rolling("int", period="3i").agg(pl.col("int").alias("aggregated"))
        shape: (5, 2)
        ┌─────┬────────────┐
        │ int ┆ aggregated │
        │ --- ┆ ---        │
        │ i64 ┆ list[i64]  │
        ╞═════╪════════════╡
        │ 0   ┆ [0]        │
        │ 4   ┆ [4]        │
        │ 5   ┆ [4, 5]     │
        │ 6   ┆ [4, 5, 6]  │
        │ 8   ┆ [6, 8]     │
        └─────┴────────────┘

        If you want the index count to be based on row number, then you may want to
        combine `rolling` with :meth:`.with_row_index`.
        """
        return RollingGroupBy(
            self,
            index_column=index_column,
            period=period,
            offset=offset,
            closed=closed,
            group_by=group_by,
        )

    @deprecate_renamed_parameter("by", "group_by", version="0.20.14")
    def group_by_dynamic(
        self,
        index_column: IntoExpr,
        *,
        every: str | timedelta,
        period: str | timedelta | None = None,
        offset: str | timedelta | None = None,
        include_boundaries: bool = False,
        closed: ClosedInterval = "left",
        label: Label = "left",
        group_by: IntoExpr | Iterable[IntoExpr] | None = None,
        start_by: StartBy = "window",
    ) -> DynamicGroupBy:
        """
        Group based on a time value (or index value of type Int32, Int64).

        Time windows are calculated and rows are assigned to windows. Different from a
        normal group by is that a row can be member of multiple groups.
        By default, the windows look like:

        - [start, start + period)
        - [start + every, start + every + period)
        - [start + 2*every, start + 2*every + period)
        - ...

        where `start` is determined by `start_by`, `offset`, `every`, and the earliest
        datapoint. See the `start_by` argument description for details.

        .. warning::
            The index column must be sorted in ascending order. If `group_by` is passed, then
            the index column must be sorted in ascending order within each group.

        .. versionchanged:: 0.20.14
            The `by` parameter was renamed `group_by`.

        Parameters
        ----------
        index_column
            Column used to group based on the time window.
            Often of type Date/Datetime.
            This column must be sorted in ascending order (or, if `group_by` is specified,
            then it must be sorted in ascending order within each group).

            In case of a dynamic group by on indices, dtype needs to be one of
            {Int32, Int64}. Note that Int32 gets temporarily cast to Int64, so if
            performance matters use an Int64 column.
        every
            interval of the window
        period
            length of the window, if None it will equal 'every'
        offset
            offset of the window, does not take effect if `start_by` is 'datapoint'.
            Defaults to zero.
        include_boundaries
            Add the lower and upper bound of the window to the "_lower_boundary" and
            "_upper_boundary" columns. This will impact performance because it's harder to
            parallelize
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive).
        label : {'left', 'right', 'datapoint'}
            Define which label to use for the window:

            - 'left': lower boundary of the window
            - 'right': upper boundary of the window
            - 'datapoint': the first value of the index column in the given window.
              If you don't need the label to be at one of the boundaries, choose this
              option for maximum performance
        group_by
            Also group by this column/these columns
        start_by : {'window', 'datapoint', 'monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday'}
            The strategy to determine the start of the first window by.

            * 'window': Start by taking the earliest timestamp, truncating it with
              `every`, and then adding `offset`.
              Note that weekly windows start on Monday.
            * 'datapoint': Start from the first encountered data point.
            * a day of the week (only takes effect if `every` contains `'w'`):

              * 'monday': Start the window on the Monday before the first data point.
              * 'tuesday': Start the window on the Tuesday before the first data point.
              * ...
              * 'sunday': Start the window on the Sunday before the first data point.

              The resulting window is then shifted back until the earliest datapoint
              is in or in front of it.

        Returns
        -------
        DynamicGroupBy
            Object you can call `.agg` on to aggregate by groups, the result
            of which will be sorted by `index_column` (but note that if `group_by` columns are
            passed, it will only be sorted within each group).

        See Also
        --------
        rolling

        Notes
        -----
        1) If you're coming from pandas, then

           .. code-block:: python

               # polars
               df.group_by_dynamic("ts", every="1d").agg(pl.col("value").sum())

           is equivalent to

           .. code-block:: python

               # pandas
               df.set_index("ts").resample("D")["value"].sum().reset_index()

           though note that, unlike pandas, polars doesn't add extra rows for empty
           windows. If you need `index_column` to be evenly spaced, then please combine
           with :func:`DataFrame.upsample`.

        2) The `every`, `period` and `offset` arguments are created with
           the following string language:

           - 1ns   (1 nanosecond)
           - 1us   (1 microsecond)
           - 1ms   (1 millisecond)
           - 1s    (1 second)
           - 1m    (1 minute)
           - 1h    (1 hour)
           - 1d    (1 calendar day)
           - 1w    (1 calendar week)
           - 1mo   (1 calendar month)
           - 1q    (1 calendar quarter)
           - 1y    (1 calendar year)
           - 1i    (1 index count)

           Or combine them (except in `every`):
           "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

           By "calendar day", we mean the corresponding time on the next day (which may
           not be 24 hours, due to daylight savings). Similarly for "calendar week",
           "calendar month", "calendar quarter", and "calendar year".

           In case of a group_by_dynamic on an integer column, the windows are defined by:

           - "1i"      # length 1
           - "10i"     # length 10

        Examples
        --------
        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "time": pl.datetime_range(
        ...             start=datetime(2021, 12, 16),
        ...             end=datetime(2021, 12, 16, 3),
        ...             interval="30m",
        ...             eager=True,
        ...         ),
        ...         "n": range(7),
        ...     }
        ... )
        >>> df
        shape: (7, 2)
        ┌─────────────────────┬─────┐
        │ time                ┆ n   │
        │ ---                 ┆ --- │
        │ datetime[μs]        ┆ i64 │
        ╞═════════════════════╪═════╡
        │ 2021-12-16 00:00:00 ┆ 0   │
        │ 2021-12-16 00:30:00 ┆ 1   │
        │ 2021-12-16 01:00:00 ┆ 2   │
        │ 2021-12-16 01:30:00 ┆ 3   │
        │ 2021-12-16 02:00:00 ┆ 4   │
        │ 2021-12-16 02:30:00 ┆ 5   │
        │ 2021-12-16 03:00:00 ┆ 6   │
        └─────────────────────┴─────┘

        Group by windows of 1 hour.

        >>> df.group_by_dynamic("time", every="1h", closed="right").agg(pl.col("n"))
        shape: (4, 2)
        ┌─────────────────────┬───────────┐
        │ time                ┆ n         │
        │ ---                 ┆ ---       │
        │ datetime[μs]        ┆ list[i64] │
        ╞═════════════════════╪═══════════╡
        │ 2021-12-15 23:00:00 ┆ [0]       │
        │ 2021-12-16 00:00:00 ┆ [1, 2]    │
        │ 2021-12-16 01:00:00 ┆ [3, 4]    │
        │ 2021-12-16 02:00:00 ┆ [5, 6]    │
        └─────────────────────┴───────────┘

        The window boundaries can also be added to the aggregation result

        >>> df.group_by_dynamic(
        ...     "time", every="1h", include_boundaries=True, closed="right"
        ... ).agg(pl.col("n").mean())
        shape: (4, 4)
        ┌─────────────────────┬─────────────────────┬─────────────────────┬─────┐
        │ _lower_boundary     ┆ _upper_boundary     ┆ time                ┆ n   │
        │ ---                 ┆ ---                 ┆ ---                 ┆ --- │
        │ datetime[μs]        ┆ datetime[μs]        ┆ datetime[μs]        ┆ f64 │
        ╞═════════════════════╪═════════════════════╪═════════════════════╪═════╡
        │ 2021-12-15 23:00:00 ┆ 2021-12-16 00:00:00 ┆ 2021-12-15 23:00:00 ┆ 0.0 │
        │ 2021-12-16 00:00:00 ┆ 2021-12-16 01:00:00 ┆ 2021-12-16 00:00:00 ┆ 1.5 │
        │ 2021-12-16 01:00:00 ┆ 2021-12-16 02:00:00 ┆ 2021-12-16 01:00:00 ┆ 3.5 │
        │ 2021-12-16 02:00:00 ┆ 2021-12-16 03:00:00 ┆ 2021-12-16 02:00:00 ┆ 5.5 │
        └─────────────────────┴─────────────────────┴─────────────────────┴─────┘

        When closed="left", the window excludes the right end of interval:
        [lower_bound, upper_bound)

        >>> df.group_by_dynamic("time", every="1h", closed="left").agg(pl.col("n"))
        shape: (4, 2)
        ┌─────────────────────┬───────────┐
        │ time                ┆ n         │
        │ ---                 ┆ ---       │
        │ datetime[μs]        ┆ list[i64] │
        ╞═════════════════════╪═══════════╡
        │ 2021-12-16 00:00:00 ┆ [0, 1]    │
        │ 2021-12-16 01:00:00 ┆ [2, 3]    │
        │ 2021-12-16 02:00:00 ┆ [4, 5]    │
        │ 2021-12-16 03:00:00 ┆ [6]       │
        └─────────────────────┴───────────┘

        When closed="both" the time values at the window boundaries belong to 2 groups.

        >>> df.group_by_dynamic("time", every="1h", closed="both").agg(pl.col("n"))
        shape: (4, 2)
        ┌─────────────────────┬───────────┐
        │ time                ┆ n         │
        │ ---                 ┆ ---       │
        │ datetime[μs]        ┆ list[i64] │
        ╞═════════════════════╪═══════════╡
        │ 2021-12-16 00:00:00 ┆ [0, 1, 2] │
        │ 2021-12-16 01:00:00 ┆ [2, 3, 4] │
        │ 2021-12-16 02:00:00 ┆ [4, 5, 6] │
        │ 2021-12-16 03:00:00 ┆ [6]       │
        └─────────────────────┴───────────┘

        Dynamic group bys can also be combined with grouping on normal keys

        >>> df = df.with_columns(groups=pl.Series(["a", "a", "a", "b", "b", "a", "a"]))
        >>> df
        shape: (7, 3)
        ┌─────────────────────┬─────┬────────┐
        │ time                ┆ n   ┆ groups │
        │ ---                 ┆ --- ┆ ---    │
        │ datetime[μs]        ┆ i64 ┆ str    │
        ╞═════════════════════╪═════╪════════╡
        │ 2021-12-16 00:00:00 ┆ 0   ┆ a      │
        │ 2021-12-16 00:30:00 ┆ 1   ┆ a      │
        │ 2021-12-16 01:00:00 ┆ 2   ┆ a      │
        │ 2021-12-16 01:30:00 ┆ 3   ┆ b      │
        │ 2021-12-16 02:00:00 ┆ 4   ┆ b      │
        │ 2021-12-16 02:30:00 ┆ 5   ┆ a      │
        │ 2021-12-16 03:00:00 ┆ 6   ┆ a      │
        └─────────────────────┴─────┴────────┘
        >>> df.group_by_dynamic(
        ...     "time",
        ...     every="1h",
        ...     closed="both",
        ...     group_by="groups",
        ...     include_boundaries=True,
        ... ).agg(pl.col("n"))
        shape: (6, 5)
        ┌────────┬─────────────────────┬─────────────────────┬─────────────────────┬───────────┐
        │ groups ┆ _lower_boundary     ┆ _upper_boundary     ┆ time                ┆ n         │
        │ ---    ┆ ---                 ┆ ---                 ┆ ---                 ┆ ---       │
        │ str    ┆ datetime[μs]        ┆ datetime[μs]        ┆ datetime[μs]        ┆ list[i64] │
        ╞════════╪═════════════════════╪═════════════════════╪═════════════════════╪═══════════╡
        │ a      ┆ 2021-12-16 00:00:00 ┆ 2021-12-16 01:00:00 ┆ 2021-12-16 00:00:00 ┆ [0, 1, 2] │
        │ a      ┆ 2021-12-16 01:00:00 ┆ 2021-12-16 02:00:00 ┆ 2021-12-16 01:00:00 ┆ [2]       │
        │ a      ┆ 2021-12-16 02:00:00 ┆ 2021-12-16 03:00:00 ┆ 2021-12-16 02:00:00 ┆ [5, 6]    │
        │ a      ┆ 2021-12-16 03:00:00 ┆ 2021-12-16 04:00:00 ┆ 2021-12-16 03:00:00 ┆ [6]       │
        │ b      ┆ 2021-12-16 01:00:00 ┆ 2021-12-16 02:00:00 ┆ 2021-12-16 01:00:00 ┆ [3, 4]    │
        │ b      ┆ 2021-12-16 02:00:00 ┆ 2021-12-16 03:00:00 ┆ 2021-12-16 02:00:00 ┆ [4]       │
        └────────┴─────────────────────┴─────────────────────┴─────────────────────┴───────────┘

        Dynamic group by on an index column

        >>> df = pl.DataFrame(
        ...     {
        ...         "idx": pl.int_range(0, 6, eager=True),
        ...         "A": ["A", "A", "B", "B", "B", "C"],
        ...     }
        ... )
        >>> (
        ...     df.group_by_dynamic(
        ...         "idx",
        ...         every="2i",
        ...         period="3i",
        ...         include_boundaries=True,
        ...         closed="right",
        ...     ).agg(pl.col("A").alias("A_agg_list"))
        ... )
        shape: (4, 4)
        ┌─────────────────┬─────────────────┬─────┬─────────────────┐
        │ _lower_boundary ┆ _upper_boundary ┆ idx ┆ A_agg_list      │
        │ ---             ┆ ---             ┆ --- ┆ ---             │
        │ i64             ┆ i64             ┆ i64 ┆ list[str]       │
        ╞═════════════════╪═════════════════╪═════╪═════════════════╡
        │ -2              ┆ 1               ┆ -2  ┆ ["A", "A"]      │
        │ 0               ┆ 3               ┆ 0   ┆ ["A", "B", "B"] │
        │ 2               ┆ 5               ┆ 2   ┆ ["B", "B", "C"] │
        │ 4               ┆ 7               ┆ 4   ┆ ["C"]           │
        └─────────────────┴─────────────────┴─────┴─────────────────┘
        """  # noqa: W505
        return DynamicGroupBy(
            self,
            index_column=index_column,
            every=every,
            period=period,
            offset=offset,
            label=label,
            include_boundaries=include_boundaries,
            closed=closed,
            group_by=group_by,
            start_by=start_by,
        )

    @deprecate_renamed_parameter("by", "group_by", version="0.20.14")
    def upsample(
        self,
        time_column: str,
        *,
        every: str | timedelta,
        group_by: str | Sequence[str] | None = None,
        maintain_order: bool = False,
    ) -> DataFrame:
        """
        Upsample a DataFrame at a regular frequency.

        The `every` argument is created with the following string language:

        - 1ns   (1 nanosecond)
        - 1us   (1 microsecond)
        - 1ms   (1 millisecond)
        - 1s    (1 second)
        - 1m    (1 minute)
        - 1h    (1 hour)
        - 1d    (1 calendar day)
        - 1w    (1 calendar week)
        - 1mo   (1 calendar month)
        - 1q    (1 calendar quarter)
        - 1y    (1 calendar year)
        - 1i    (1 index count)

        Or combine them:

        - "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

        By "calendar day", we mean the corresponding time on the next day (which may
        not be 24 hours, due to daylight savings). Similarly for "calendar week",
        "calendar month", "calendar quarter", and "calendar year".

        .. versionchanged:: 0.20.14
            The `by` parameter was renamed `group_by`.

        Parameters
        ----------
        time_column
            Time column will be used to determine a date_range.
            Note that this column has to be sorted for the output to make sense.
        every
            Interval will start 'every' duration.
        group_by
            First group by these columns and then upsample for every group.
        maintain_order
            Keep the ordering predictable. This is slower.

        Returns
        -------
        DataFrame
            Result will be sorted by `time_column` (but note that if `group_by` columns
            are passed, it will only be sorted within each group).

        Examples
        --------
        Upsample a DataFrame by a certain interval.

        >>> from datetime import datetime
        >>> df = pl.DataFrame(
        ...     {
        ...         "time": [
        ...             datetime(2021, 2, 1),
        ...             datetime(2021, 4, 1),
        ...             datetime(2021, 5, 1),
        ...             datetime(2021, 6, 1),
        ...         ],
        ...         "groups": ["A", "B", "A", "B"],
        ...         "values": [0, 1, 2, 3],
        ...     }
        ... ).set_sorted("time")
        >>> df.upsample(
        ...     time_column="time", every="1mo", group_by="groups", maintain_order=True
        ... ).select(pl.all().fill_null(strategy="forward"))
        shape: (7, 3)
        ┌─────────────────────┬────────┬────────┐
        │ time                ┆ groups ┆ values │
        │ ---                 ┆ ---    ┆ ---    │
        │ datetime[μs]        ┆ str    ┆ i64    │
        ╞═════════════════════╪════════╪════════╡
        │ 2021-02-01 00:00:00 ┆ A      ┆ 0      │
        │ 2021-03-01 00:00:00 ┆ A      ┆ 0      │
        │ 2021-04-01 00:00:00 ┆ A      ┆ 0      │
        │ 2021-05-01 00:00:00 ┆ A      ┆ 2      │
        │ 2021-04-01 00:00:00 ┆ B      ┆ 1      │
        │ 2021-05-01 00:00:00 ┆ B      ┆ 1      │
        │ 2021-06-01 00:00:00 ┆ B      ┆ 3      │
        └─────────────────────┴────────┴────────┘
        """
        if group_by is None:
            group_by = []
        if isinstance(group_by, str):
            group_by = [group_by]

        every = parse_as_duration_string(every)

        return self._from_pydf(
            self._df.upsample(group_by, time_column, every, maintain_order)
        )

    def join_asof(
        self,
        other: DataFrame,
        *,
        left_on: str | None | Expr = None,
        right_on: str | None | Expr = None,
        on: str | None | Expr = None,
        by_left: str | Sequence[str] | None = None,
        by_right: str | Sequence[str] | None = None,
        by: str | Sequence[str] | None = None,
        strategy: AsofJoinStrategy = "backward",
        suffix: str = "_right",
        tolerance: str | int | float | timedelta | None = None,
        allow_parallel: bool = True,
        force_parallel: bool = False,
        coalesce: bool = True,
        allow_exact_matches: bool = True,
        check_sortedness: bool = True,
    ) -> DataFrame:
        """
        Perform an asof join.

        This is similar to a left-join except that we match on nearest key rather than
        equal keys.

        Both DataFrames must be sorted by the `on` key (within each `by` group, if
        specified).

        For each row in the left DataFrame:

          - A "backward" search selects the last row in the right DataFrame whose
            'on' key is less than or equal to the left's key.

          - A "forward" search selects the first row in the right DataFrame whose
            'on' key is greater than or equal to the left's key.

          - A "nearest" search selects the last row in the right DataFrame whose value
            is nearest to the left's key. String keys are not currently supported for a
            nearest search.

        The default is "backward".

        Parameters
        ----------
        other
            Lazy DataFrame to join with.
        left_on
            Join column of the left DataFrame.
        right_on
            Join column of the right DataFrame.
        on
            Join column of both DataFrames. If set, `left_on` and `right_on` should be
            None.
        by
            Join on these columns before doing asof join
        by_left
            Join on these columns before doing asof join
        by_right
            Join on these columns before doing asof join
        strategy : {'backward', 'forward', 'nearest'}
            Join strategy.
        suffix
            Suffix to append to columns with a duplicate name.
        tolerance
            Numeric tolerance. By setting this the join will only be done if the near
            keys are within this distance. If an asof join is done on columns of dtype
            "Date", "Datetime", "Duration" or "Time", use either a datetime.timedelta
            object or the following string language:

                - 1ns   (1 nanosecond)
                - 1us   (1 microsecond)
                - 1ms   (1 millisecond)
                - 1s    (1 second)
                - 1m    (1 minute)
                - 1h    (1 hour)
                - 1d    (1 calendar day)
                - 1w    (1 calendar week)
                - 1mo   (1 calendar month)
                - 1q    (1 calendar quarter)
                - 1y    (1 calendar year)

                Or combine them:
                "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

                By "calendar day", we mean the corresponding time on the next day
                (which may not be 24 hours, due to daylight savings). Similarly for
                "calendar week", "calendar month", "calendar quarter", and
                "calendar year".

        allow_parallel
            Allow the physical plan to optionally evaluate the computation of both
            DataFrames up to the join in parallel.
        force_parallel
            Force the physical plan to evaluate the computation of both DataFrames up to
            the join in parallel.
        coalesce
            Coalescing behavior (merging of `on` / `left_on` / `right_on` columns):

            - *True*: Always coalesce join columns.
            - *False*: Never coalesce join columns.

            Note that joining on any other expressions than `col`
            will turn off coalescing.
        allow_exact_matches
            Whether exact matches are valid join predicates.

            - If True, allow matching with the same ``on`` value
                (i.e. less-than-or-equal-to / greater-than-or-equal-to)
            - If False, don't match the same ``on`` value
                (i.e., strictly less-than / strictly greater-than).
        check_sortedness
            Check the sortedness of the asof keys. If the keys are not sorted Polars
            will error. Currently, sortedness cannot be checked if 'by' groups are
            provided.

        Examples
        --------
        >>> from datetime import date
        >>> gdp = pl.DataFrame(
        ...     {
        ...         "date": pl.date_range(
        ...             date(2016, 1, 1),
        ...             date(2020, 1, 1),
        ...             "1y",
        ...             eager=True,
        ...         ),
        ...         "gdp": [4164, 4411, 4566, 4696, 4827],
        ...     }
        ... )
        >>> gdp
        shape: (5, 2)
        ┌────────────┬──────┐
        │ date       ┆ gdp  │
        │ ---        ┆ ---  │
        │ date       ┆ i64  │
        ╞════════════╪══════╡
        │ 2016-01-01 ┆ 4164 │
        │ 2017-01-01 ┆ 4411 │
        │ 2018-01-01 ┆ 4566 │
        │ 2019-01-01 ┆ 4696 │
        │ 2020-01-01 ┆ 4827 │
        └────────────┴──────┘

        >>> population = pl.DataFrame(
        ...     {
        ...         "date": [date(2016, 3, 1), date(2018, 8, 1), date(2019, 1, 1)],
        ...         "population": [82.19, 82.66, 83.12],
        ...     }
        ... ).sort("date")
        >>> population
        shape: (3, 2)
        ┌────────────┬────────────┐
        │ date       ┆ population │
        │ ---        ┆ ---        │
        │ date       ┆ f64        │
        ╞════════════╪════════════╡
        │ 2016-03-01 ┆ 82.19      │
        │ 2018-08-01 ┆ 82.66      │
        │ 2019-01-01 ┆ 83.12      │
        └────────────┴────────────┘

        Note how the dates don't quite match. If we join them using `join_asof` and
        `strategy='backward'`, then each date from `population` which doesn't have an
        exact match is matched with the closest earlier date from `gdp`:

        >>> population.join_asof(gdp, on="date", strategy="backward")
        shape: (3, 3)
        ┌────────────┬────────────┬──────┐
        │ date       ┆ population ┆ gdp  │
        │ ---        ┆ ---        ┆ ---  │
        │ date       ┆ f64        ┆ i64  │
        ╞════════════╪════════════╪══════╡
        │ 2016-03-01 ┆ 82.19      ┆ 4164 │
        │ 2018-08-01 ┆ 82.66      ┆ 4566 │
        │ 2019-01-01 ┆ 83.12      ┆ 4696 │
        └────────────┴────────────┴──────┘

        Note how:

        - date `2016-03-01` from `population` is matched with `2016-01-01` from `gdp`;
        - date `2018-08-01` from `population` is matched with `2018-01-01` from `gdp`.

        You can verify this by passing `coalesce=False`:

        >>> population.join_asof(gdp, on="date", strategy="backward", coalesce=False)
        shape: (3, 4)
        ┌────────────┬────────────┬────────────┬──────┐
        │ date       ┆ population ┆ date_right ┆ gdp  │
        │ ---        ┆ ---        ┆ ---        ┆ ---  │
        │ date       ┆ f64        ┆ date       ┆ i64  │
        ╞════════════╪════════════╪════════════╪══════╡
        │ 2016-03-01 ┆ 82.19      ┆ 2016-01-01 ┆ 4164 │
        │ 2018-08-01 ┆ 82.66      ┆ 2018-01-01 ┆ 4566 │
        │ 2019-01-01 ┆ 83.12      ┆ 2019-01-01 ┆ 4696 │
        └────────────┴────────────┴────────────┴──────┘

        If we instead use `strategy='forward'`, then each date from `population` which
        doesn't have an exact match is matched with the closest later date from `gdp`:

        >>> population.join_asof(gdp, on="date", strategy="forward")
        shape: (3, 3)
        ┌────────────┬────────────┬──────┐
        │ date       ┆ population ┆ gdp  │
        │ ---        ┆ ---        ┆ ---  │
        │ date       ┆ f64        ┆ i64  │
        ╞════════════╪════════════╪══════╡
        │ 2016-03-01 ┆ 82.19      ┆ 4411 │
        │ 2018-08-01 ┆ 82.66      ┆ 4696 │
        │ 2019-01-01 ┆ 83.12      ┆ 4696 │
        └────────────┴────────────┴──────┘

        Note how:

        - date `2016-03-01` from `population` is matched with `2017-01-01` from `gdp`;
        - date `2018-08-01` from `population` is matched with `2019-01-01` from `gdp`.

        Finally, `strategy='nearest'` gives us a mix of the two results above, as each
        date from `population` which doesn't have an exact match is matched with the
        closest date from `gdp`, regardless of whether it's earlier or later:

        >>> population.join_asof(gdp, on="date", strategy="nearest")
        shape: (3, 3)
        ┌────────────┬────────────┬──────┐
        │ date       ┆ population ┆ gdp  │
        │ ---        ┆ ---        ┆ ---  │
        │ date       ┆ f64        ┆ i64  │
        ╞════════════╪════════════╪══════╡
        │ 2016-03-01 ┆ 82.19      ┆ 4164 │
        │ 2018-08-01 ┆ 82.66      ┆ 4696 │
        │ 2019-01-01 ┆ 83.12      ┆ 4696 │
        └────────────┴────────────┴──────┘

        Note how:

        - date `2016-03-01` from `population` is matched with `2016-01-01` from `gdp`;
        - date `2018-08-01` from `population` is matched with `2019-01-01` from `gdp`.

        They `by` argument allows joining on another column first, before the asof join.
        In this example we join by `country` first, then asof join by date, as above.

        >>> gdp_dates = pl.date_range(  # fmt: skip
        ...     date(2016, 1, 1), date(2020, 1, 1), "1y", eager=True
        ... )
        >>> gdp2 = pl.DataFrame(
        ...     {
        ...         "country": ["Germany"] * 5 + ["Netherlands"] * 5,
        ...         "date": pl.concat([gdp_dates, gdp_dates]),
        ...         "gdp": [4164, 4411, 4566, 4696, 4827, 784, 833, 914, 910, 909],
        ...     }
        ... ).sort("country", "date")
        >>>
        >>> gdp2
        shape: (10, 3)
        ┌─────────────┬────────────┬──────┐
        │ country     ┆ date       ┆ gdp  │
        │ ---         ┆ ---        ┆ ---  │
        │ str         ┆ date       ┆ i64  │
        ╞═════════════╪════════════╪══════╡
        │ Germany     ┆ 2016-01-01 ┆ 4164 │
        │ Germany     ┆ 2017-01-01 ┆ 4411 │
        │ Germany     ┆ 2018-01-01 ┆ 4566 │
        │ Germany     ┆ 2019-01-01 ┆ 4696 │
        │ Germany     ┆ 2020-01-01 ┆ 4827 │
        │ Netherlands ┆ 2016-01-01 ┆ 784  │
        │ Netherlands ┆ 2017-01-01 ┆ 833  │
        │ Netherlands ┆ 2018-01-01 ┆ 914  │
        │ Netherlands ┆ 2019-01-01 ┆ 910  │
        │ Netherlands ┆ 2020-01-01 ┆ 909  │
        └─────────────┴────────────┴──────┘
        >>> pop2 = pl.DataFrame(
        ...     {
        ...         "country": ["Germany"] * 3 + ["Netherlands"] * 3,
        ...         "date": [
        ...             date(2016, 3, 1),
        ...             date(2018, 8, 1),
        ...             date(2019, 1, 1),
        ...             date(2016, 3, 1),
        ...             date(2018, 8, 1),
        ...             date(2019, 1, 1),
        ...         ],
        ...         "population": [82.19, 82.66, 83.12, 17.11, 17.32, 17.40],
        ...     }
        ... ).sort("country", "date")
        >>>
        >>> pop2
        shape: (6, 3)
        ┌─────────────┬────────────┬────────────┐
        │ country     ┆ date       ┆ population │
        │ ---         ┆ ---        ┆ ---        │
        │ str         ┆ date       ┆ f64        │
        ╞═════════════╪════════════╪════════════╡
        │ Germany     ┆ 2016-03-01 ┆ 82.19      │
        │ Germany     ┆ 2018-08-01 ┆ 82.66      │
        │ Germany     ┆ 2019-01-01 ┆ 83.12      │
        │ Netherlands ┆ 2016-03-01 ┆ 17.11      │
        │ Netherlands ┆ 2018-08-01 ┆ 17.32      │
        │ Netherlands ┆ 2019-01-01 ┆ 17.4       │
        └─────────────┴────────────┴────────────┘
        >>> pop2.join_asof(gdp2, by="country", on="date", strategy="nearest")
        shape: (6, 4)
        ┌─────────────┬────────────┬────────────┬──────┐
        │ country     ┆ date       ┆ population ┆ gdp  │
        │ ---         ┆ ---        ┆ ---        ┆ ---  │
        │ str         ┆ date       ┆ f64        ┆ i64  │
        ╞═════════════╪════════════╪════════════╪══════╡
        │ Germany     ┆ 2016-03-01 ┆ 82.19      ┆ 4164 │
        │ Germany     ┆ 2018-08-01 ┆ 82.66      ┆ 4696 │
        │ Germany     ┆ 2019-01-01 ┆ 83.12      ┆ 4696 │
        │ Netherlands ┆ 2016-03-01 ┆ 17.11      ┆ 784  │
        │ Netherlands ┆ 2018-08-01 ┆ 17.32      ┆ 910  │
        │ Netherlands ┆ 2019-01-01 ┆ 17.4       ┆ 910  │
        └─────────────┴────────────┴────────────┴──────┘
        """
        require_same_type(self, other)

        if on is not None:
            if not isinstance(on, (str, pl.Expr)):
                msg = (
                    f"expected `on` to be str or Expr, got {qualified_type_name(on)!r}"
                )
                raise TypeError(msg)
        else:
            if not isinstance(left_on, (str, pl.Expr)):
                msg = f"expected `left_on` to be str or Expr, got {qualified_type_name(left_on)!r}"
                raise TypeError(msg)
            elif not isinstance(right_on, (str, pl.Expr)):
                msg = f"expected `right_on` to be str or Expr, got {qualified_type_name(right_on)!r}"
                raise TypeError(msg)

        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .join_asof(
                other.lazy(),
                left_on=left_on,
                right_on=right_on,
                on=on,
                by_left=by_left,
                by_right=by_right,
                by=by,
                strategy=strategy,
                suffix=suffix,
                tolerance=tolerance,
                allow_parallel=allow_parallel,
                force_parallel=force_parallel,
                coalesce=coalesce,
                allow_exact_matches=allow_exact_matches,
                check_sortedness=check_sortedness,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    @deprecate_renamed_parameter("join_nulls", "nulls_equal", version="1.24")
    def join(
        self,
        other: DataFrame,
        on: str | Expr | Sequence[str | Expr] | None = None,
        how: JoinStrategy = "inner",
        *,
        left_on: str | Expr | Sequence[str | Expr] | None = None,
        right_on: str | Expr | Sequence[str | Expr] | None = None,
        suffix: str = "_right",
        validate: JoinValidation = "m:m",
        nulls_equal: bool = False,
        coalesce: bool | None = None,
        maintain_order: MaintainOrderJoin | None = None,
    ) -> DataFrame:
        """
        Join in SQL-like fashion.

        .. versionchanged:: 1.24
            The `join_nulls` parameter was renamed `nulls_equal`.

        Parameters
        ----------
        other
            DataFrame to join with.
        on
            Name(s) of the join columns in both DataFrames. If set, `left_on` and
            `right_on` should be None. This should not be specified if `how='cross'`.
        how : {'inner', 'left', 'right', 'full', 'semi', 'anti', 'cross'}
            Join strategy.

            .. list-table ::
               :header-rows: 0

               * - **inner**
                 - *(Default)* Returns rows that have matching values in both tables.
               * - **left**
                 - Returns all rows from the left table, and the matched rows from
                   the right table.
               * - **full**
                 - Returns all rows when there is a match in either left or right.
               * - **cross**
                 - Returns the Cartesian product of rows from both tables
               * - **semi**
                 - Returns rows from the left table that have a match in the right
                   table.
               * - **anti**
                 - Returns rows from the left table that have no match in the right
                   table.

        left_on
            Name(s) of the left join column(s).
        right_on
            Name(s) of the right join column(s).
        suffix
            Suffix to append to columns with a duplicate name.
        validate: {'m:m', 'm:1', '1:m', '1:1'}
            Checks if join is of specified type.

            .. list-table ::
               :header-rows: 0

               * - **m:m**
                 - *(Default)* Many-to-many (default). Does not result in checks.
               * - **1:1**
                 - One-to-one. Checks if join keys are unique in both left and
                   right datasets.
               * - **1:m**
                 - One-to-many. Checks if join keys are unique in left dataset.
               * - **m:1**
                 - Many-to-one. Check if join keys are unique in right dataset.

            .. note::
                This is currently not supported by the streaming engine.

        nulls_equal
            Join on null values. By default null values will never produce matches.
        coalesce
            Coalescing behavior (merging of join columns).

            .. list-table ::
               :header-rows: 0

               * - **None**
                 - *(Default)* Coalesce unless `how='full'` is specified.
               * - **True**
                 - Always coalesce join columns.
               * - **False**
                 - Never coalesce join columns.

            .. note::
                Joining on any other expressions than `col`
                will turn off coalescing.
        maintain_order : {'none', 'left', 'right', 'left_right', 'right_left'}
            Which DataFrame row order to preserve, if any.
            Do not rely on any observed ordering without explicitly setting this
            parameter, as your code may break in a future release.
            Not specifying any ordering can improve performance.

            .. list-table ::
               :header-rows: 0

               * - **none**
                 - *(Default)* No specific ordering is desired. The ordering might
                   differ across Polars versions or even between different runs.
               * - **left**
                 - Preserves the order of the left DataFrame.
               * - **right**
                 - Preserves the order of the right DataFrame.
               * - **left_right**
                 - First preserves the order of the left DataFrame, then the right.
               * - **right_left**
                 - First preserves the order of the right DataFrame, then the left.

        See Also
        --------
        join_asof

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> other_df = pl.DataFrame(
        ...     {
        ...         "apple": ["x", "y", "z"],
        ...         "ham": ["a", "b", "d"],
        ...     }
        ... )
        >>> df.join(other_df, on="ham")
        shape: (2, 4)
        ┌─────┬─────┬─────┬───────┐
        │ foo ┆ bar ┆ ham ┆ apple │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ i64 ┆ f64 ┆ str ┆ str   │
        ╞═════╪═════╪═════╪═══════╡
        │ 1   ┆ 6.0 ┆ a   ┆ x     │
        │ 2   ┆ 7.0 ┆ b   ┆ y     │
        └─────┴─────┴─────┴───────┘

        >>> df.join(other_df, on="ham", how="full")
        shape: (4, 5)
        ┌──────┬──────┬──────┬───────┬───────────┐
        │ foo  ┆ bar  ┆ ham  ┆ apple ┆ ham_right │
        │ ---  ┆ ---  ┆ ---  ┆ ---   ┆ ---       │
        │ i64  ┆ f64  ┆ str  ┆ str   ┆ str       │
        ╞══════╪══════╪══════╪═══════╪═══════════╡
        │ 1    ┆ 6.0  ┆ a    ┆ x     ┆ a         │
        │ 2    ┆ 7.0  ┆ b    ┆ y     ┆ b         │
        │ null ┆ null ┆ null ┆ z     ┆ d         │
        │ 3    ┆ 8.0  ┆ c    ┆ null  ┆ null      │
        └──────┴──────┴──────┴───────┴───────────┘

        >>> df.join(other_df, on="ham", how="full", coalesce=True)
        shape: (4, 4)
        ┌──────┬──────┬─────┬───────┐
        │ foo  ┆ bar  ┆ ham ┆ apple │
        │ ---  ┆ ---  ┆ --- ┆ ---   │
        │ i64  ┆ f64  ┆ str ┆ str   │
        ╞══════╪══════╪═════╪═══════╡
        │ 1    ┆ 6.0  ┆ a   ┆ x     │
        │ 2    ┆ 7.0  ┆ b   ┆ y     │
        │ null ┆ null ┆ d   ┆ z     │
        │ 3    ┆ 8.0  ┆ c   ┆ null  │
        └──────┴──────┴─────┴───────┘

        >>> df.join(other_df, on="ham", how="left")
        shape: (3, 4)
        ┌─────┬─────┬─────┬───────┐
        │ foo ┆ bar ┆ ham ┆ apple │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ i64 ┆ f64 ┆ str ┆ str   │
        ╞═════╪═════╪═════╪═══════╡
        │ 1   ┆ 6.0 ┆ a   ┆ x     │
        │ 2   ┆ 7.0 ┆ b   ┆ y     │
        │ 3   ┆ 8.0 ┆ c   ┆ null  │
        └─────┴─────┴─────┴───────┘

        >>> df.join(other_df, on="ham", how="semi")
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6.0 ┆ a   │
        │ 2   ┆ 7.0 ┆ b   │
        └─────┴─────┴─────┘

        >>> df.join(other_df, on="ham", how="anti")
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 3   ┆ 8.0 ┆ c   │
        └─────┴─────┴─────┘

        >>> df.join(other_df, how="cross")
        shape: (9, 5)
        ┌─────┬─────┬─────┬───────┬───────────┐
        │ foo ┆ bar ┆ ham ┆ apple ┆ ham_right │
        │ --- ┆ --- ┆ --- ┆ ---   ┆ ---       │
        │ i64 ┆ f64 ┆ str ┆ str   ┆ str       │
        ╞═════╪═════╪═════╪═══════╪═══════════╡
        │ 1   ┆ 6.0 ┆ a   ┆ x     ┆ a         │
        │ 1   ┆ 6.0 ┆ a   ┆ y     ┆ b         │
        │ 1   ┆ 6.0 ┆ a   ┆ z     ┆ d         │
        │ 2   ┆ 7.0 ┆ b   ┆ x     ┆ a         │
        │ 2   ┆ 7.0 ┆ b   ┆ y     ┆ b         │
        │ 2   ┆ 7.0 ┆ b   ┆ z     ┆ d         │
        │ 3   ┆ 8.0 ┆ c   ┆ x     ┆ a         │
        │ 3   ┆ 8.0 ┆ c   ┆ y     ┆ b         │
        │ 3   ┆ 8.0 ┆ c   ┆ z     ┆ d         │
        └─────┴─────┴─────┴───────┴───────────┘

        Notes
        -----
        For joining on columns with categorical data, see :class:`polars.StringCache`.
        """
        require_same_type(self, other)

        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .join(
                other=other.lazy(),
                left_on=left_on,
                right_on=right_on,
                on=on,
                how=how,
                suffix=suffix,
                validate=validate,
                nulls_equal=nulls_equal,
                coalesce=coalesce,
                maintain_order=maintain_order,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    @unstable()
    def join_where(
        self,
        other: DataFrame,
        *predicates: Expr | Iterable[Expr],
        suffix: str = "_right",
    ) -> DataFrame:
        """
        Perform a join based on one or multiple (in)equality predicates.

        This performs an inner join, so only rows where all predicates are true
        are included in the result, and a row from either DataFrame may be included
        multiple times in the result.

        .. note::
            The row order of the input DataFrames is not preserved.

        .. warning::
            This functionality is experimental. It may be
            changed at any point without it being considered a breaking change.

        Parameters
        ----------
        other
            DataFrame to join with.
        *predicates
            (In)Equality condition to join the two tables on.
            When a column name occurs in both tables, the proper suffix must
            be applied in the predicate.
        suffix
            Suffix to append to columns with a duplicate name.

        Examples
        --------
        Join two dataframes together based on two predicates which get AND-ed together.

        >>> east = pl.DataFrame(
        ...     {
        ...         "id": [100, 101, 102],
        ...         "dur": [120, 140, 160],
        ...         "rev": [12, 14, 16],
        ...         "cores": [2, 8, 4],
        ...     }
        ... )
        >>> west = pl.DataFrame(
        ...     {
        ...         "t_id": [404, 498, 676, 742],
        ...         "time": [90, 130, 150, 170],
        ...         "cost": [9, 13, 15, 16],
        ...         "cores": [4, 2, 1, 4],
        ...     }
        ... )
        >>> east.join_where(
        ...     west,
        ...     pl.col("dur") < pl.col("time"),
        ...     pl.col("rev") < pl.col("cost"),
        ... )
        shape: (5, 8)
        ┌─────┬─────┬─────┬───────┬──────┬──────┬──────┬─────────────┐
        │ id  ┆ dur ┆ rev ┆ cores ┆ t_id ┆ time ┆ cost ┆ cores_right │
        │ --- ┆ --- ┆ --- ┆ ---   ┆ ---  ┆ ---  ┆ ---  ┆ ---         │
        │ i64 ┆ i64 ┆ i64 ┆ i64   ┆ i64  ┆ i64  ┆ i64  ┆ i64         │
        ╞═════╪═════╪═════╪═══════╪══════╪══════╪══════╪═════════════╡
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 498  ┆ 130  ┆ 13   ┆ 2           │
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 676  ┆ 150  ┆ 15   ┆ 1           │
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 742  ┆ 170  ┆ 16   ┆ 4           │
        │ 101 ┆ 140 ┆ 14  ┆ 8     ┆ 676  ┆ 150  ┆ 15   ┆ 1           │
        │ 101 ┆ 140 ┆ 14  ┆ 8     ┆ 742  ┆ 170  ┆ 16   ┆ 4           │
        └─────┴─────┴─────┴───────┴──────┴──────┴──────┴─────────────┘

        To OR them together, use a single expression and the `|` operator.

        >>> east.join_where(
        ...     west,
        ...     (pl.col("dur") < pl.col("time")) | (pl.col("rev") < pl.col("cost")),
        ... )
        shape: (6, 8)
        ┌─────┬─────┬─────┬───────┬──────┬──────┬──────┬─────────────┐
        │ id  ┆ dur ┆ rev ┆ cores ┆ t_id ┆ time ┆ cost ┆ cores_right │
        │ --- ┆ --- ┆ --- ┆ ---   ┆ ---  ┆ ---  ┆ ---  ┆ ---         │
        │ i64 ┆ i64 ┆ i64 ┆ i64   ┆ i64  ┆ i64  ┆ i64  ┆ i64         │
        ╞═════╪═════╪═════╪═══════╪══════╪══════╪══════╪═════════════╡
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 498  ┆ 130  ┆ 13   ┆ 2           │
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 676  ┆ 150  ┆ 15   ┆ 1           │
        │ 100 ┆ 120 ┆ 12  ┆ 2     ┆ 742  ┆ 170  ┆ 16   ┆ 4           │
        │ 101 ┆ 140 ┆ 14  ┆ 8     ┆ 676  ┆ 150  ┆ 15   ┆ 1           │
        │ 101 ┆ 140 ┆ 14  ┆ 8     ┆ 742  ┆ 170  ┆ 16   ┆ 4           │
        │ 102 ┆ 160 ┆ 16  ┆ 4     ┆ 742  ┆ 170  ┆ 16   ┆ 4           │
        └─────┴─────┴─────┴───────┴──────┴──────┴──────┴─────────────┘
        """
        require_same_type(self, other)

        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .join_where(
                other.lazy(),
                *predicates,
                suffix=suffix,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    def map_rows(
        self,
        function: Callable[[tuple[Any, ...]], Any],
        return_dtype: PolarsDataType | None = None,
        *,
        inference_size: int = 256,
    ) -> DataFrame:
        """
        Apply a custom/user-defined function (UDF) over the rows of the DataFrame.

        .. warning::
            This method is much slower than the native expressions API.
            Only use it if you cannot implement your logic otherwise.

        The UDF will receive each row as a tuple of values: `udf(row)`.

        Implementing logic using a Python function is almost always *significantly*
        slower and more memory intensive than implementing the same logic using
        the native expression API because:

        - The native expression engine runs in Rust; UDFs run in Python.
        - Use of Python UDFs forces the DataFrame to be materialized in memory.
        - Polars-native expressions can be parallelised (UDFs typically cannot).
        - Polars-native expressions can be logically optimised (UDFs cannot).

        Wherever possible you should strongly prefer the native expression API
        to achieve the best performance.

        Parameters
        ----------
        function
            Custom function or lambda.
        return_dtype
            Output type of the operation. If none given, Polars tries to infer the type.
        inference_size
            Only used in the case when the custom function returns rows.
            This uses the first `n` rows to determine the output schema.

        Notes
        -----
        * The frame-level `map_rows` cannot track column names (as the UDF is a
          black-box that may arbitrarily drop, rearrange, transform, or add new
          columns); if you want to apply a UDF such that column names are preserved,
          you should use the expression-level `map_elements` syntax instead.

        * If your function is expensive and you don't want it to be called more than
          once for a given input, consider applying an `@lru_cache` decorator to it.
          If your data is suitable you may achieve *significant* speedups.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [-1, 5, 8]})

        Return a DataFrame by mapping each row to a tuple:

        >>> df.map_rows(lambda t: (t[0] * 2, t[1] * 3))
        shape: (3, 2)
        ┌──────────┬──────────┐
        │ column_0 ┆ column_1 │
        │ ---      ┆ ---      │
        │ i64      ┆ i64      │
        ╞══════════╪══════════╡
        │ 2        ┆ -3       │
        │ 4        ┆ 15       │
        │ 6        ┆ 24       │
        └──────────┴──────────┘

        However, it is much better to implement this with a native expression:

        >>> df.select(
        ...     pl.col("foo") * 2,
        ...     pl.col("bar") * 3,
        ... )  # doctest: +IGNORE_RESULT

        Return a DataFrame with a single column by mapping each row to a scalar:

        >>> df.map_rows(lambda t: (t[0] * 2 + t[1]))
        shape: (3, 1)
        ┌─────┐
        │ map │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 9   │
        │ 14  │
        └─────┘

        In this case it is better to use the following native expression:

        >>> df.select(pl.col("foo") * 2 + pl.col("bar"))  # doctest: +IGNORE_RESULT
        """
        # TODO: Enable warning for inefficient map
        # from polars._utils.udfs import warn_on_inefficient_map
        # warn_on_inefficient_map(function, columns=self.columns, map_target="frame)

        out, is_df = self._df.map_rows(function, return_dtype, inference_size)
        if is_df:
            return self._from_pydf(out)
        else:
            return wrap_s(out).to_frame()

    def hstack(
        self, columns: list[Series] | DataFrame, *, in_place: bool = False
    ) -> DataFrame:
        """
        Return a new DataFrame grown horizontally by stacking multiple Series to it.

        Parameters
        ----------
        columns
            Series to stack.
        in_place
            Modify in place.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> x = pl.Series("apple", [10, 20, 30])
        >>> df.hstack([x])
        shape: (3, 4)
        ┌─────┬─────┬─────┬───────┐
        │ foo ┆ bar ┆ ham ┆ apple │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ i64 ┆ i64 ┆ str ┆ i64   │
        ╞═════╪═════╪═════╪═══════╡
        │ 1   ┆ 6   ┆ a   ┆ 10    │
        │ 2   ┆ 7   ┆ b   ┆ 20    │
        │ 3   ┆ 8   ┆ c   ┆ 30    │
        └─────┴─────┴─────┴───────┘
        """
        if not isinstance(columns, list):
            columns = columns.get_columns()
        if in_place:
            self._df.hstack_mut([s._s for s in columns])
            return self
        else:
            return self._from_pydf(self._df.hstack([s._s for s in columns]))

    def vstack(self, other: DataFrame, *, in_place: bool = False) -> DataFrame:
        """
        Grow this DataFrame vertically by stacking a DataFrame to it.

        Parameters
        ----------
        other
            DataFrame to stack.
        in_place
            Modify in place.

        See Also
        --------
        extend

        Examples
        --------
        >>> df1 = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2],
        ...         "bar": [6, 7],
        ...         "ham": ["a", "b"],
        ...     }
        ... )
        >>> df2 = pl.DataFrame(
        ...     {
        ...         "foo": [3, 4],
        ...         "bar": [8, 9],
        ...         "ham": ["c", "d"],
        ...     }
        ... )
        >>> df1.vstack(df2)
        shape: (4, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 2   ┆ 7   ┆ b   │
        │ 3   ┆ 8   ┆ c   │
        │ 4   ┆ 9   ┆ d   │
        └─────┴─────┴─────┘
        """
        require_same_type(self, other)
        if in_place:
            self._df.vstack_mut(other._df)
            return self

        return self._from_pydf(self._df.vstack(other._df))

    def extend(self, other: DataFrame) -> DataFrame:
        """
        Extend the memory backed by this `DataFrame` with the values from `other`.

        Different from `vstack` which adds the chunks from `other` to the chunks of
        this `DataFrame`, `extend` appends the data from `other` to the underlying
        memory locations and thus may cause a reallocation.

        If this does not cause a reallocation, the resulting data structure will not
        have any extra chunks and thus will yield faster queries.

        Prefer `extend` over `vstack` when you want to do a query after a single
        append. For instance, during online operations where you add `n` rows and rerun
        a query.

        Prefer `vstack` over `extend` when you want to append many times before
        doing a query. For instance, when you read in multiple files and want to store
        them in a single `DataFrame`. In the latter case, finish the sequence of
        `vstack` operations with a `rechunk`.

        Parameters
        ----------
        other
            DataFrame to vertically add.

        Warnings
        --------
        This method modifies the dataframe in-place. The dataframe is returned for
        convenience only.

        See Also
        --------
        vstack

        Examples
        --------
        >>> df1 = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> df2 = pl.DataFrame({"foo": [10, 20, 30], "bar": [40, 50, 60]})
        >>> df1.extend(df2)
        shape: (6, 2)
        ┌─────┬─────┐
        │ foo ┆ bar │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 4   │
        │ 2   ┆ 5   │
        │ 3   ┆ 6   │
        │ 10  ┆ 40  │
        │ 20  ┆ 50  │
        │ 30  ┆ 60  │
        └─────┴─────┘
        """
        require_same_type(self, other)
        self._df.extend(other._df)
        return self

    def drop(
        self,
        *columns: ColumnNameOrSelector | Iterable[ColumnNameOrSelector],
        strict: bool = True,
    ) -> DataFrame:
        """
        Remove columns from the dataframe.

        Parameters
        ----------
        *columns
            Names of the columns that should be removed from the dataframe.
            Accepts column selector input.
        strict
            Validate that all column names exist in the current schema,
            and throw an exception if any do not.

        Examples
        --------
        Drop a single column by passing the name of that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.drop("ham")
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

        Drop multiple columns by passing a list of column names.

        >>> df.drop(["bar", "ham"])
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘

        Drop multiple columns by passing a selector.

        >>> import polars.selectors as cs
        >>> df.drop(cs.numeric())
        shape: (3, 1)
        ┌─────┐
        │ ham │
        │ --- │
        │ str │
        ╞═════╡
        │ a   │
        │ b   │
        │ c   │
        └─────┘

        Use positional arguments to drop multiple columns.

        >>> df.drop("foo", "ham")
        shape: (3, 1)
        ┌─────┐
        │ bar │
        │ --- │
        │ f64 │
        ╞═════╡
        │ 6.0 │
        │ 7.0 │
        │ 8.0 │
        └─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .drop(*columns, strict=strict)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def drop_in_place(self, name: str) -> Series:
        """
        Drop a single column in-place and return the dropped column.

        Parameters
        ----------
        name
            Name of the column to drop.

        Returns
        -------
        Series
            The dropped column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.drop_in_place("ham")
        shape: (3,)
        Series: 'ham' [str]
        [
            "a"
            "b"
            "c"
        ]
        """
        return wrap_s(self._df.drop_in_place(name))

    def cast(
        self,
        dtypes: (
            Mapping[
                ColumnNameOrSelector | PolarsDataType, PolarsDataType | PythonDataType
            ]
            | PolarsDataType
        ),
        *,
        strict: bool = True,
    ) -> DataFrame:
        """
        Cast DataFrame column(s) to the specified dtype(s).

        Parameters
        ----------
        dtypes
            Mapping of column names (or selector) to dtypes, or a single dtype
            to which all columns will be cast.
        strict
            Raise if cast is invalid on rows after predicates are pushed down.
            If `False`, invalid casts will produce null values.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": [date(2020, 1, 2), date(2021, 3, 4), date(2022, 5, 6)],
        ...     }
        ... )

        Cast specific frame columns to the specified dtypes:

        >>> df.cast({"foo": pl.Float32, "bar": pl.UInt8})
        shape: (3, 3)
        ┌─────┬─────┬────────────┐
        │ foo ┆ bar ┆ ham        │
        │ --- ┆ --- ┆ ---        │
        │ f32 ┆ u8  ┆ date       │
        ╞═════╪═════╪════════════╡
        │ 1.0 ┆ 6   ┆ 2020-01-02 │
        │ 2.0 ┆ 7   ┆ 2021-03-04 │
        │ 3.0 ┆ 8   ┆ 2022-05-06 │
        └─────┴─────┴────────────┘

        Cast all frame columns matching one dtype (or dtype group) to another dtype:

        >>> df.cast({pl.Date: pl.Datetime})
        shape: (3, 3)
        ┌─────┬─────┬─────────────────────┐
        │ foo ┆ bar ┆ ham                 │
        │ --- ┆ --- ┆ ---                 │
        │ i64 ┆ f64 ┆ datetime[μs]        │
        ╞═════╪═════╪═════════════════════╡
        │ 1   ┆ 6.0 ┆ 2020-01-02 00:00:00 │
        │ 2   ┆ 7.0 ┆ 2021-03-04 00:00:00 │
        │ 3   ┆ 8.0 ┆ 2022-05-06 00:00:00 │
        └─────┴─────┴─────────────────────┘

        Use selectors to define the columns being cast:

        >>> import polars.selectors as cs
        >>> df.cast({cs.numeric(): pl.UInt32, cs.temporal(): pl.String})
        shape: (3, 3)
        ┌─────┬─────┬────────────┐
        │ foo ┆ bar ┆ ham        │
        │ --- ┆ --- ┆ ---        │
        │ u32 ┆ u32 ┆ str        │
        ╞═════╪═════╪════════════╡
        │ 1   ┆ 6   ┆ 2020-01-02 │
        │ 2   ┆ 7   ┆ 2021-03-04 │
        │ 3   ┆ 8   ┆ 2022-05-06 │
        └─────┴─────┴────────────┘

        Cast all frame columns to the specified dtype:

        >>> df.cast(pl.String).to_dict(as_series=False)
        {'foo': ['1', '2', '3'],
         'bar': ['6.0', '7.0', '8.0'],
         'ham': ['2020-01-02', '2021-03-04', '2022-05-06']}
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .cast(dtypes, strict=strict)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def clear(self, n: int = 0) -> DataFrame:
        """
        Create an empty (n=0) or `n`-row null-filled (n>0) copy of the DataFrame.

        Returns a `n`-row null-filled DataFrame with an identical schema.
        `n` can be greater than the current number of rows in the DataFrame.

        Parameters
        ----------
        n
            Number of (null-filled) rows to return in the cleared frame.

        See Also
        --------
        clone : Cheap deepcopy/clone.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [None, 2, 3, 4],
        ...         "b": [0.5, None, 2.5, 13],
        ...         "c": [True, True, False, None],
        ...     }
        ... )
        >>> df.clear()
        shape: (0, 3)
        ┌─────┬─────┬──────┐
        │ a   ┆ b   ┆ c    │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ f64 ┆ bool │
        ╞═════╪═════╪══════╡
        └─────┴─────┴──────┘

        >>> df.clear(n=2)
        shape: (2, 3)
        ┌──────┬──────┬──────┐
        │ a    ┆ b    ┆ c    │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ f64  ┆ bool │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ null ┆ null ┆ null │
        └──────┴──────┴──────┘
        """
        if n < 0:
            msg = f"`n` should be greater than or equal to 0, got {n}"
            raise ValueError(msg)
        # faster path
        if n == 0:
            return self._from_pydf(self._df.clear())
        return self.__class__(
            {
                nm: pl.Series(name=nm, dtype=tp).extend_constant(None, n)
                for nm, tp in self.schema.items()
            }
        )

    def clone(self) -> DataFrame:
        """
        Create a copy of this DataFrame.

        This is a cheap operation that does not copy data.

        See Also
        --------
        clear : Create an empty copy of the current DataFrame, with identical
            schema but no data.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [0.5, 4, 10, 13],
        ...         "c": [True, True, False, True],
        ...     }
        ... )
        >>> df.clone()
        shape: (4, 3)
        ┌─────┬──────┬───────┐
        │ a   ┆ b    ┆ c     │
        │ --- ┆ ---  ┆ ---   │
        │ i64 ┆ f64  ┆ bool  │
        ╞═════╪══════╪═══════╡
        │ 1   ┆ 0.5  ┆ true  │
        │ 2   ┆ 4.0  ┆ true  │
        │ 3   ┆ 10.0 ┆ false │
        │ 4   ┆ 13.0 ┆ true  │
        └─────┴──────┴───────┘
        """
        return self._from_pydf(self._df.clone())

    def get_columns(self) -> list[Series]:
        """
        Get the DataFrame as a List of Series.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> df.get_columns()
        [shape: (3,)
        Series: 'foo' [i64]
        [
                1
                2
                3
        ], shape: (3,)
        Series: 'bar' [i64]
        [
                4
                5
                6
        ]]

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [0.5, 4, 10, 13],
        ...         "c": [True, True, False, True],
        ...     }
        ... )
        >>> df.get_columns()
        [shape: (4,)
        Series: 'a' [i64]
        [
            1
            2
            3
            4
        ], shape: (4,)
        Series: 'b' [f64]
        [
            0.5
            4.0
            10.0
            13.0
        ], shape: (4,)
        Series: 'c' [bool]
        [
            true
            true
            false
            true
        ]]
        """
        return [wrap_s(s) for s in self._df.get_columns()]

    @overload
    def get_column(self, name: str, *, default: Series | NoDefault = ...) -> Series: ...

    @overload
    def get_column(self, name: str, *, default: Any) -> Any: ...

    def get_column(
        self, name: str, *, default: Any | NoDefault = no_default
    ) -> Series | Any:
        """
        Get a single column by name.

        Parameters
        ----------
        name
            String name of the column to retrieve.
        default
            Value to return if the column does not exist; if not explicitly set and
            the column is not present a `ColumnNotFoundError` exception is raised.

        Returns
        -------
        Series (or arbitrary default value, if specified).

        See Also
        --------
        to_series

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> df.get_column("foo")
        shape: (3,)
        Series: 'foo' [i64]
        [
            1
            2
            3
        ]

        Missing column handling; can optionally provide an arbitrary default value
        to the method (otherwise a `ColumnNotFoundError` exception is raised).

        >>> df.get_column("baz", default=pl.Series("baz", ["?", "?", "?"]))
        shape: (3,)
        Series: 'baz' [str]
        [
            "?"
            "?"
            "?"
        ]
        >>> res = df.get_column("baz", default=None)
        >>> res is None
        True
        """
        try:
            return wrap_s(self._df.get_column(name))
        except ColumnNotFoundError:
            if default is no_default:
                raise
            return default

    def fill_null(
        self,
        value: Any | Expr | None = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
        *,
        matches_supertype: bool = True,
    ) -> DataFrame:
        """
        Fill null values using the specified value or strategy.

        Parameters
        ----------
        value
            Value used to fill null values.
        strategy : {None, 'forward', 'backward', 'min', 'max', 'mean', 'zero', 'one'}
            Strategy used to fill null values.
        limit
            Number of consecutive null values to fill when using the 'forward' or
            'backward' strategy.
        matches_supertype
            Fill all matching supertype of the fill `value`.

        Returns
        -------
        DataFrame
            DataFrame with None values replaced by the filling strategy.

        See Also
        --------
        fill_nan

        Notes
        -----
        A null value is not the same as a NaN value.
        To fill NaN values, use :func:`fill_nan`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, None, 4],
        ...         "b": [0.5, 4, None, 13],
        ...     }
        ... )
        >>> df.fill_null(99)
        shape: (4, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ f64  │
        ╞═════╪══════╡
        │ 1   ┆ 0.5  │
        │ 2   ┆ 4.0  │
        │ 99  ┆ 99.0 │
        │ 4   ┆ 13.0 │
        └─────┴──────┘
        >>> df.fill_null(strategy="forward")
        shape: (4, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ f64  │
        ╞═════╪══════╡
        │ 1   ┆ 0.5  │
        │ 2   ┆ 4.0  │
        │ 2   ┆ 4.0  │
        │ 4   ┆ 13.0 │
        └─────┴──────┘

        >>> df.fill_null(strategy="max")
        shape: (4, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ f64  │
        ╞═════╪══════╡
        │ 1   ┆ 0.5  │
        │ 2   ┆ 4.0  │
        │ 4   ┆ 13.0 │
        │ 4   ┆ 13.0 │
        └─────┴──────┘

        >>> df.fill_null(strategy="zero")
        shape: (4, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ f64  │
        ╞═════╪══════╡
        │ 1   ┆ 0.5  │
        │ 2   ┆ 4.0  │
        │ 0   ┆ 0.0  │
        │ 4   ┆ 13.0 │
        └─────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .fill_null(value, strategy, limit, matches_supertype=matches_supertype)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def fill_nan(self, value: Expr | int | float | None) -> DataFrame:
        """
        Fill floating point NaN values by an Expression evaluation.

        Parameters
        ----------
        value
            Value used to fill NaN values.

        Returns
        -------
        DataFrame
            DataFrame with NaN values replaced by the given value.

        See Also
        --------
        fill_null

        Notes
        -----
        A NaN value is not the same as a null value.
        To fill null values, use :func:`fill_null`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1.5, 2, float("nan"), 4],
        ...         "b": [0.5, 4, float("nan"), 13],
        ...     }
        ... )
        >>> df.fill_nan(99)
        shape: (4, 2)
        ┌──────┬──────┐
        │ a    ┆ b    │
        │ ---  ┆ ---  │
        │ f64  ┆ f64  │
        ╞══════╪══════╡
        │ 1.5  ┆ 0.5  │
        │ 2.0  ┆ 4.0  │
        │ 99.0 ┆ 99.0 │
        │ 4.0  ┆ 13.0 │
        └──────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().fill_nan(value).collect(optimizations=QueryOptFlags._eager())

    def explode(
        self,
        columns: ColumnNameOrSelector | Iterable[ColumnNameOrSelector],
        *more_columns: ColumnNameOrSelector,
    ) -> DataFrame:
        """
        Explode the dataframe to long format by exploding the given columns.

        Parameters
        ----------
        columns
            Column names, expressions, or a selector defining them. The underlying
            columns being exploded must be of the `List` or `Array` data type.
        *more_columns
            Additional names of columns to explode, specified as positional arguments.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "letters": ["a", "a", "b", "c"],
        ...         "numbers": [[1], [2, 3], [4, 5], [6, 7, 8]],
        ...     }
        ... )
        >>> df
        shape: (4, 2)
        ┌─────────┬───────────┐
        │ letters ┆ numbers   │
        │ ---     ┆ ---       │
        │ str     ┆ list[i64] │
        ╞═════════╪═══════════╡
        │ a       ┆ [1]       │
        │ a       ┆ [2, 3]    │
        │ b       ┆ [4, 5]    │
        │ c       ┆ [6, 7, 8] │
        └─────────┴───────────┘
        >>> df.explode("numbers")
        shape: (8, 2)
        ┌─────────┬─────────┐
        │ letters ┆ numbers │
        │ ---     ┆ ---     │
        │ str     ┆ i64     │
        ╞═════════╪═════════╡
        │ a       ┆ 1       │
        │ a       ┆ 2       │
        │ a       ┆ 3       │
        │ b       ┆ 4       │
        │ b       ┆ 5       │
        │ c       ┆ 6       │
        │ c       ┆ 7       │
        │ c       ┆ 8       │
        └─────────┴─────────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .explode(columns, *more_columns)
            .collect(optimizations=QueryOptFlags._eager())
        )

    @deprecate_renamed_parameter("columns", "on", version="1.0.0")
    def pivot(
        self,
        on: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        index: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        values: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        aggregate_function: PivotAgg | Expr | None = None,
        maintain_order: bool = True,
        sort_columns: bool = False,
        separator: str = "_",
    ) -> DataFrame:
        """
        Create a spreadsheet-style pivot table as a DataFrame.

        Only available in eager mode. See "Examples" section below for how to do a
        "lazy pivot" if you know the unique column values in advance.

        .. versionchanged:: 1.0.0
            The `columns` parameter was renamed `on`.

        Parameters
        ----------
        on
            The column(s) whose values will be used as the new columns of the output
            DataFrame.
        index
            The column(s) that remain from the input to the output. The output DataFrame will have one row
            for each unique combination of the `index`'s values.
            If None, all remaining columns not specified on `on` and `values` will be used. At least one
            of `index` and `values` must be specified.
        values
            The existing column(s) of values which will be moved under the new columns from index. If an
            aggregation is specified, these are the values on which the aggregation will be computed.
            If None, all remaining columns not specified on `on` and `index` will be used.
            At least one of `index` and `values` must be specified.
        aggregate_function
            Choose from:

            - None: no aggregation takes place, will raise error if multiple values are in group.
            - A predefined aggregate function string, one of
              {'min', 'max', 'first', 'last', 'sum', 'mean', 'median', 'len'}
            - An expression to do the aggregation. The expression can only access data from the respective
              'values' columns as generated by pivot, through `pl.element()`.
        maintain_order
            Ensure the values of `index` are sorted by discovery order.
        sort_columns
            Sort the transposed columns by name. Default is by order of discovery.
        separator
            Used as separator/delimiter in generated column names in case of multiple
            `values` columns.

        Returns
        -------
        DataFrame

        Notes
        -----
        In some other frameworks, you might know this operation as `pivot_wider`.

        Examples
        --------
        You can use `pivot` to reshape a dataframe from "long" to "wide" format.

        For example, suppose we have a dataframe of test scores achieved by some
        students, where each row represents a distinct test.

        >>> df = pl.DataFrame(
        ...     {
        ...         "name": ["Cady", "Cady", "Karen", "Karen"],
        ...         "subject": ["maths", "physics", "maths", "physics"],
        ...         "test_1": [98, 99, 61, 58],
        ...         "test_2": [100, 100, 60, 60],
        ...     }
        ... )
        >>> df
        shape: (4, 4)
        ┌───────┬─────────┬────────┬────────┐
        │ name  ┆ subject ┆ test_1 ┆ test_2 │
        │ ---   ┆ ---     ┆ ---    ┆ ---    │
        │ str   ┆ str     ┆ i64    ┆ i64    │
        ╞═══════╪═════════╪════════╪════════╡
        │ Cady  ┆ maths   ┆ 98     ┆ 100    │
        │ Cady  ┆ physics ┆ 99     ┆ 100    │
        │ Karen ┆ maths   ┆ 61     ┆ 60     │
        │ Karen ┆ physics ┆ 58     ┆ 60     │
        └───────┴─────────┴────────┴────────┘

        Using `pivot`, we can reshape so we have one row per student, with different
        subjects as columns, and their `test_1` scores as values:

        >>> df.pivot("subject", index="name", values="test_1")
        shape: (2, 3)
        ┌───────┬───────┬─────────┐
        │ name  ┆ maths ┆ physics │
        │ ---   ┆ ---   ┆ ---     │
        │ str   ┆ i64   ┆ i64     │
        ╞═══════╪═══════╪═════════╡
        │ Cady  ┆ 98    ┆ 99      │
        │ Karen ┆ 61    ┆ 58      │
        └───────┴───────┴─────────┘

        You can use selectors too - here we include all test scores in the pivoted table:

        >>> import polars.selectors as cs
        >>> df.pivot("subject", values=cs.starts_with("test"))
        shape: (2, 5)
        ┌───────┬──────────────┬────────────────┬──────────────┬────────────────┐
        │ name  ┆ test_1_maths ┆ test_1_physics ┆ test_2_maths ┆ test_2_physics │
        │ ---   ┆ ---          ┆ ---            ┆ ---          ┆ ---            │
        │ str   ┆ i64          ┆ i64            ┆ i64          ┆ i64            │
        ╞═══════╪══════════════╪════════════════╪══════════════╪════════════════╡
        │ Cady  ┆ 98           ┆ 99             ┆ 100          ┆ 100            │
        │ Karen ┆ 61           ┆ 58             ┆ 60           ┆ 60             │
        └───────┴──────────────┴────────────────┴──────────────┴────────────────┘

        If you end up with multiple values per cell, you can specify how to aggregate
        them with `aggregate_function`:

        >>> df = pl.DataFrame(
        ...     {
        ...         "ix": [1, 1, 2, 2, 1, 2],
        ...         "col": ["a", "a", "a", "a", "b", "b"],
        ...         "foo": [0, 1, 2, 2, 7, 1],
        ...         "bar": [0, 2, 0, 0, 9, 4],
        ...     }
        ... )
        >>> df.pivot("col", index="ix", aggregate_function="sum")
        shape: (2, 5)
        ┌─────┬───────┬───────┬───────┬───────┐
        │ ix  ┆ foo_a ┆ foo_b ┆ bar_a ┆ bar_b │
        │ --- ┆ ---   ┆ ---   ┆ ---   ┆ ---   │
        │ i64 ┆ i64   ┆ i64   ┆ i64   ┆ i64   │
        ╞═════╪═══════╪═══════╪═══════╪═══════╡
        │ 1   ┆ 1     ┆ 7     ┆ 2     ┆ 9     │
        │ 2   ┆ 4     ┆ 1     ┆ 0     ┆ 4     │
        └─────┴───────┴───────┴───────┴───────┘

        You can also pass a custom aggregation function using
        :meth:`polars.element`:

        >>> df = pl.DataFrame(
        ...     {
        ...         "col1": ["a", "a", "a", "b", "b", "b"],
        ...         "col2": ["x", "x", "x", "x", "y", "y"],
        ...         "col3": [6, 7, 3, 2, 5, 7],
        ...     }
        ... )
        >>> df.pivot(
        ...     "col2",
        ...     index="col1",
        ...     values="col3",
        ...     aggregate_function=pl.element().tanh().mean(),
        ... )
        shape: (2, 3)
        ┌──────┬──────────┬──────────┐
        │ col1 ┆ x        ┆ y        │
        │ ---  ┆ ---      ┆ ---      │
        │ str  ┆ f64      ┆ f64      │
        ╞══════╪══════════╪══════════╡
        │ a    ┆ 0.998347 ┆ null     │
        │ b    ┆ 0.964028 ┆ 0.999954 │
        └──────┴──────────┴──────────┘

        Note that `pivot` is only available in eager mode. If you know the unique
        column values in advance, you can use :meth:`polars.LazyFrame.group_by` to
        get the same result as above in lazy mode:

        >>> index = pl.col("col1")
        >>> on = pl.col("col2")
        >>> values = pl.col("col3")
        >>> unique_column_values = ["x", "y"]
        >>> aggregate_function = lambda col: col.tanh().mean()
        >>> df.lazy().group_by(index).agg(
        ...     aggregate_function(values.filter(on == value)).alias(value)
        ...     for value in unique_column_values
        ... ).collect()  # doctest: +IGNORE_RESULT
        shape: (2, 3)
        ┌──────┬──────────┬──────────┐
        │ col1 ┆ x        ┆ y        │
        │ ---  ┆ ---      ┆ ---      │
        │ str  ┆ f64      ┆ f64      │
        ╞══════╪══════════╪══════════╡
        │ a    ┆ 0.998347 ┆ null     │
        │ b    ┆ 0.964028 ┆ 0.999954 │
        └──────┴──────────┴──────────┘
        """  # noqa: W505
        on = _expand_selectors(self, on)
        if values is not None:
            values = _expand_selectors(self, values)
        if index is not None:
            index = _expand_selectors(self, index)

        if isinstance(aggregate_function, str):
            if aggregate_function == "first":
                aggregate_expr = F.element().first()._pyexpr
            elif aggregate_function == "sum":
                aggregate_expr = F.element().sum()._pyexpr
            elif aggregate_function == "max":
                aggregate_expr = F.element().max()._pyexpr
            elif aggregate_function == "min":
                aggregate_expr = F.element().min()._pyexpr
            elif aggregate_function == "mean":
                aggregate_expr = F.element().mean()._pyexpr
            elif aggregate_function == "median":
                aggregate_expr = F.element().median()._pyexpr
            elif aggregate_function == "last":
                aggregate_expr = F.element().last()._pyexpr
            elif aggregate_function == "len":
                aggregate_expr = F.len()._pyexpr
            elif aggregate_function == "count":
                issue_deprecation_warning(
                    "`aggregate_function='count'` input for `pivot` is deprecated."
                    " Please use `aggregate_function='len'`.",
                    version="0.20.5",
                )
                aggregate_expr = F.len()._pyexpr
            else:
                msg = f"invalid input for `aggregate_function` argument: {aggregate_function!r}"
                raise ValueError(msg)
        elif aggregate_function is None:
            aggregate_expr = None
        else:
            aggregate_expr = aggregate_function._pyexpr

        return self._from_pydf(
            self._df.pivot_expr(
                on,
                index,
                values,
                maintain_order,
                sort_columns,
                aggregate_expr,
                separator,
            )
        )

    def unpivot(
        self,
        on: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        *,
        index: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        variable_name: str | None = None,
        value_name: str | None = None,
    ) -> DataFrame:
        """
        Unpivot a DataFrame from wide to long format.

        Optionally leaves identifiers set.

        This function is useful to massage a DataFrame into a format where one or more
        columns are identifier variables (index) while all other columns, considered
        measured variables (on), are "unpivoted" to the row axis leaving just
        two non-identifier columns, 'variable' and 'value'.

        Parameters
        ----------
        on
            Column(s) or selector(s) to use as values variables; if `on`
            is empty all columns that are not in `index` will be used.
        index
            Column(s) or selector(s) to use as identifier variables.
        variable_name
            Name to give to the `variable` column. Defaults to "variable"
        value_name
            Name to give to the `value` column. Defaults to "value"

        Notes
        -----
        If you're coming from pandas, this is similar to `pandas.DataFrame.melt`,
        but with `index` replacing `id_vars` and `on` replacing `value_vars`.
        In other frameworks, you might know this operation as `pivot_longer`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["x", "y", "z"],
        ...         "b": [1, 3, 5],
        ...         "c": [2, 4, 6],
        ...     }
        ... )
        >>> import polars.selectors as cs
        >>> df.unpivot(cs.numeric(), index="a")
        shape: (6, 3)
        ┌─────┬──────────┬───────┐
        │ a   ┆ variable ┆ value │
        │ --- ┆ ---      ┆ ---   │
        │ str ┆ str      ┆ i64   │
        ╞═════╪══════════╪═══════╡
        │ x   ┆ b        ┆ 1     │
        │ y   ┆ b        ┆ 3     │
        │ z   ┆ b        ┆ 5     │
        │ x   ┆ c        ┆ 2     │
        │ y   ┆ c        ┆ 4     │
        │ z   ┆ c        ┆ 6     │
        └─────┴──────────┴───────┘
        """
        on = [] if on is None else _expand_selectors(self, on)
        index = [] if index is None else _expand_selectors(self, index)

        return self._from_pydf(self._df.unpivot(on, index, value_name, variable_name))

    def unstack(
        self,
        *,
        step: int,
        how: UnstackDirection = "vertical",
        columns: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        fill_values: list[Any] | None = None,
    ) -> DataFrame:
        """
        Unstack a long table to a wide form without doing an aggregation.

        This can be much faster than a pivot, because it can skip the grouping phase.

        Parameters
        ----------
        step
            Number of rows in the unstacked frame.
        how : { 'vertical', 'horizontal' }
            Direction of the unstack.
        columns
            Column name(s) or selector(s) to include in the operation.
            If set to `None` (default), use all columns.
        fill_values
            Fill values that don't fit the new size with this value.

        Examples
        --------
        >>> from string import ascii_uppercase
        >>> df = pl.DataFrame(
        ...     {
        ...         "x": list(ascii_uppercase[0:8]),
        ...         "y": pl.int_range(1, 9, eager=True),
        ...     }
        ... ).with_columns(
        ...     z=pl.int_ranges(pl.col("y"), pl.col("y") + 2, dtype=pl.UInt8),
        ... )
        >>> df
        shape: (8, 3)
        ┌─────┬─────┬──────────┐
        │ x   ┆ y   ┆ z        │
        │ --- ┆ --- ┆ ---      │
        │ str ┆ i64 ┆ list[u8] │
        ╞═════╪═════╪══════════╡
        │ A   ┆ 1   ┆ [1, 2]   │
        │ B   ┆ 2   ┆ [2, 3]   │
        │ C   ┆ 3   ┆ [3, 4]   │
        │ D   ┆ 4   ┆ [4, 5]   │
        │ E   ┆ 5   ┆ [5, 6]   │
        │ F   ┆ 6   ┆ [6, 7]   │
        │ G   ┆ 7   ┆ [7, 8]   │
        │ H   ┆ 8   ┆ [8, 9]   │
        └─────┴─────┴──────────┘
        >>> df.unstack(step=4, how="vertical")
        shape: (4, 6)
        ┌─────┬─────┬─────┬─────┬──────────┬──────────┐
        │ x_0 ┆ x_1 ┆ y_0 ┆ y_1 ┆ z_0      ┆ z_1      │
        │ --- ┆ --- ┆ --- ┆ --- ┆ ---      ┆ ---      │
        │ str ┆ str ┆ i64 ┆ i64 ┆ list[u8] ┆ list[u8] │
        ╞═════╪═════╪═════╪═════╪══════════╪══════════╡
        │ A   ┆ E   ┆ 1   ┆ 5   ┆ [1, 2]   ┆ [5, 6]   │
        │ B   ┆ F   ┆ 2   ┆ 6   ┆ [2, 3]   ┆ [6, 7]   │
        │ C   ┆ G   ┆ 3   ┆ 7   ┆ [3, 4]   ┆ [7, 8]   │
        │ D   ┆ H   ┆ 4   ┆ 8   ┆ [4, 5]   ┆ [8, 9]   │
        └─────┴─────┴─────┴─────┴──────────┴──────────┘
        >>> df.unstack(step=2, how="horizontal")
        shape: (4, 6)
        ┌─────┬─────┬─────┬─────┬──────────┬──────────┐
        │ x_0 ┆ x_1 ┆ y_0 ┆ y_1 ┆ z_0      ┆ z_1      │
        │ --- ┆ --- ┆ --- ┆ --- ┆ ---      ┆ ---      │
        │ str ┆ str ┆ i64 ┆ i64 ┆ list[u8] ┆ list[u8] │
        ╞═════╪═════╪═════╪═════╪══════════╪══════════╡
        │ A   ┆ B   ┆ 1   ┆ 2   ┆ [1, 2]   ┆ [2, 3]   │
        │ C   ┆ D   ┆ 3   ┆ 4   ┆ [3, 4]   ┆ [4, 5]   │
        │ E   ┆ F   ┆ 5   ┆ 6   ┆ [5, 6]   ┆ [6, 7]   │
        │ G   ┆ H   ┆ 7   ┆ 8   ┆ [7, 8]   ┆ [8, 9]   │
        └─────┴─────┴─────┴─────┴──────────┴──────────┘
        >>> import polars.selectors as cs
        >>> df.unstack(step=5, columns=cs.numeric(), fill_values=0)
        shape: (5, 2)
        ┌─────┬─────┐
        │ y_0 ┆ y_1 │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 6   │
        │ 2   ┆ 7   │
        │ 3   ┆ 8   │
        │ 4   ┆ 0   │
        │ 5   ┆ 0   │
        └─────┴─────┘
        """
        import math

        df = self.select(columns) if columns is not None else self

        height = df.height
        if how == "vertical":
            n_rows = step
            n_cols = math.ceil(height / n_rows)
        else:
            n_cols = step
            n_rows = math.ceil(height / n_cols)

        if n_fill := n_cols * n_rows - height:
            if not isinstance(fill_values, list):
                fill_values = [fill_values for _ in range(df.width)]

            df = df.select(
                s.extend_constant(next_fill, n_fill)
                for s, next_fill in zip(df, fill_values)
            )

        if how == "horizontal":
            df = (
                df.with_columns(
                    (F.int_range(0, n_cols * n_rows, eager=True) % n_cols).alias(
                        "__sort_order"
                    ),
                )
                .sort("__sort_order")
                .drop("__sort_order")
            )

        zfill_val = math.floor(math.log10(n_cols)) + 1
        slices = [
            s.slice(slice_nbr * n_rows, n_rows).alias(
                s.name + "_" + str(slice_nbr).zfill(zfill_val)
            )
            for s in df
            for slice_nbr in range(n_cols)
        ]

        return DataFrame(slices)

    @overload
    def partition_by(
        self,
        by: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *more_by: ColumnNameOrSelector,
        maintain_order: bool = ...,
        include_key: bool = ...,
        as_dict: Literal[False] = ...,
    ) -> list[DataFrame]: ...

    @overload
    def partition_by(
        self,
        by: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *more_by: ColumnNameOrSelector,
        maintain_order: bool = ...,
        include_key: bool = ...,
        as_dict: Literal[True],
    ) -> dict[tuple[Any, ...], DataFrame]: ...

    @overload
    def partition_by(
        self,
        by: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *more_by: ColumnNameOrSelector,
        maintain_order: bool = ...,
        include_key: bool = ...,
        as_dict: bool,
    ) -> list[DataFrame] | dict[tuple[Any, ...], DataFrame]: ...

    def partition_by(
        self,
        by: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *more_by: ColumnNameOrSelector,
        maintain_order: bool = True,
        include_key: bool = True,
        as_dict: bool = False,
    ) -> list[DataFrame] | dict[tuple[Any, ...], DataFrame]:
        """
        Group by the given columns and return the groups as separate dataframes.

        Parameters
        ----------
        by
            Column name(s) or selector(s) to group by.
        *more_by
            Additional names of columns to group by, specified as positional arguments.
        maintain_order
            Ensure that the order of the groups is consistent with the input data.
            This is slower than a default partition by operation.
        include_key
            Include the columns used to partition the DataFrame in the output.
        as_dict
            Return a dictionary instead of a list. The dictionary keys are tuples of
            the distinct group values that identify each group.

        Examples
        --------
        Pass a single column name to partition by that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "c"],
        ...         "b": [1, 2, 1, 3, 3],
        ...         "c": [5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> df.partition_by("a")  # doctest: +IGNORE_RESULT
        [shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 5   │
        │ a   ┆ 1   ┆ 3   │
        └─────┴─────┴─────┘,
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 2   ┆ 4   │
        │ b   ┆ 3   ┆ 2   │
        └─────┴─────┴─────┘,
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ c   ┆ 3   ┆ 1   │
        └─────┴─────┴─────┘]

        Partition by multiple columns by either passing a list of column names, or by
        specifying each column name as a positional argument.

        >>> df.partition_by("a", "b")  # doctest: +IGNORE_RESULT
        [shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 5   │
        │ a   ┆ 1   ┆ 3   │
        └─────┴─────┴─────┘,
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 2   ┆ 4   │
        └─────┴─────┴─────┘,
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 3   ┆ 2   │
        └─────┴─────┴─────┘,
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ c   ┆ 3   ┆ 1   │
        └─────┴─────┴─────┘]

        Return the partitions as a dictionary by specifying `as_dict=True`.

        >>> import polars.selectors as cs
        >>> df.partition_by(cs.string(), as_dict=True)  # doctest: +IGNORE_RESULT
        {('a',): shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 1   ┆ 5   │
        │ a   ┆ 1   ┆ 3   │
        └─────┴─────┴─────┘,
        ('b',): shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 2   ┆ 4   │
        │ b   ┆ 3   ┆ 2   │
        └─────┴─────┴─────┘,
        ('c',): shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ c   ┆ 3   ┆ 1   │
        └─────┴─────┴─────┘}
        """
        by_parsed = _expand_selectors(self, by, *more_by)

        partitions = [
            self._from_pydf(_df)
            for _df in self._df.partition_by(by_parsed, maintain_order, include_key)
        ]

        if as_dict:
            if include_key:
                names = [p.select(by_parsed).row(0) for p in partitions]
            else:
                if not maintain_order:  # Group keys cannot be matched to partitions
                    msg = "cannot use `partition_by` with `maintain_order=False, include_key=False, as_dict=True`"
                    raise ValueError(msg)
                names = self.select(by_parsed).unique(maintain_order=True).rows()

            return dict(zip(names, partitions))

        return partitions

    def shift(self, n: int = 1, *, fill_value: IntoExpr | None = None) -> DataFrame:
        """
        Shift values by the given number of indices.

        Parameters
        ----------
        n
            Number of indices to shift forward. If a negative value is passed, values
            are shifted in the opposite direction instead.
        fill_value
            Fill the resulting null values with this value. Accepts scalar expression
            input. Non-expression inputs are parsed as literals.

        Notes
        -----
        This method is similar to the `LAG` operation in SQL when the value for `n`
        is positive. With a negative value for `n`, it is similar to `LEAD`.

        Examples
        --------
        By default, values are shifted forward by one index.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [5, 6, 7, 8],
        ...     }
        ... )
        >>> df.shift()
        shape: (4, 2)
        ┌──────┬──────┐
        │ a    ┆ b    │
        │ ---  ┆ ---  │
        │ i64  ┆ i64  │
        ╞══════╪══════╡
        │ null ┆ null │
        │ 1    ┆ 5    │
        │ 2    ┆ 6    │
        │ 3    ┆ 7    │
        └──────┴──────┘

        Pass a negative value to shift in the opposite direction instead.

        >>> df.shift(-2)
        shape: (4, 2)
        ┌──────┬──────┐
        │ a    ┆ b    │
        │ ---  ┆ ---  │
        │ i64  ┆ i64  │
        ╞══════╪══════╡
        │ 3    ┆ 7    │
        │ 4    ┆ 8    │
        │ null ┆ null │
        │ null ┆ null │
        └──────┴──────┘

        Specify `fill_value` to fill the resulting null values.

        >>> df.shift(-2, fill_value=100)
        shape: (4, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 3   ┆ 7   │
        │ 4   ┆ 8   │
        │ 100 ┆ 100 │
        │ 100 ┆ 100 │
        └─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .shift(n, fill_value=fill_value)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def is_duplicated(self) -> Series:
        """
        Get a mask of all duplicated rows in this DataFrame.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 1],
        ...         "b": ["x", "y", "z", "x"],
        ...     }
        ... )
        >>> df.is_duplicated()
        shape: (4,)
        Series: '' [bool]
        [
                true
                false
                false
                true
        ]

        This mask can be used to visualize the duplicated lines like this:

        >>> df.filter(df.is_duplicated())
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ x   │
        │ 1   ┆ x   │
        └─────┴─────┘
        """
        return wrap_s(self._df.is_duplicated())

    def is_unique(self) -> Series:
        """
        Get a mask of all unique rows in this DataFrame.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 1],
        ...         "b": ["x", "y", "z", "x"],
        ...     }
        ... )
        >>> df.is_unique()
        shape: (4,)
        Series: '' [bool]
        [
                false
                true
                true
                false
        ]

        This mask can be used to visualize the unique lines like this:

        >>> df.filter(df.is_unique())
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 2   ┆ y   │
        │ 3   ┆ z   │
        └─────┴─────┘
        """
        return wrap_s(self._df.is_unique())

    def lazy(self) -> LazyFrame:
        """
        Start a lazy query from this point. This returns a `LazyFrame` object.

        Operations on a `LazyFrame` are not executed until this is triggered
        by calling one of:

        * :meth:`.collect() <polars.LazyFrame.collect>`
            (run on all data)
        * :meth:`.explain() <polars.LazyFrame.explain>`
            (print the query plan)
        * :meth:`.show_graph() <polars.LazyFrame.show_graph>`
            (show the query plan as graphviz graph)
        * :meth:`.collect_schema() <polars.LazyFrame.collect_schema>`
            (return the final frame schema)

        Lazy operations are recommended because they allow for query optimization and
        additional parallelism.

        Returns
        -------
        LazyFrame

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [None, 2, 3, 4],
        ...         "b": [0.5, None, 2.5, 13],
        ...         "c": [True, True, False, None],
        ...     }
        ... )
        >>> df.lazy()
        <LazyFrame at ...>
        """
        return wrap_ldf(self._df.lazy())

    def select(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> DataFrame:
        """
        Select columns from this DataFrame.

        Parameters
        ----------
        *exprs
            Column(s) to select, specified as positional arguments.
            Accepts expression input. Strings are parsed as column names,
            other non-expression inputs are parsed as literals.
        **named_exprs
            Additional columns to select, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        Examples
        --------
        Pass the name of a column to select that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.select("foo")
        shape: (3, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘

        Multiple columns can be selected by passing a list of column names.

        >>> df.select(["foo", "bar"])
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

        Multiple columns can also be selected using positional arguments instead of a
        list. Expressions are also accepted.

        >>> df.select(pl.col("foo"), pl.col("bar") + 1)
        shape: (3, 2)
        ┌─────┬─────┐
        │ foo ┆ bar │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 7   │
        │ 2   ┆ 8   │
        │ 3   ┆ 9   │
        └─────┴─────┘

        Use keyword arguments to easily name your expression inputs.

        >>> df.select(threshold=pl.when(pl.col("foo") > 2).then(10).otherwise(0))
        shape: (3, 1)
        ┌───────────┐
        │ threshold │
        │ ---       │
        │ i32       │
        ╞═══════════╡
        │ 0         │
        │ 0         │
        │ 10        │
        └───────────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .select(*exprs, **named_exprs)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def select_seq(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> DataFrame:
        """
        Select columns from this DataFrame.

        This will run all expression sequentially instead of in parallel.
        Use this when the work per expression is cheap.

        Parameters
        ----------
        *exprs
            Column(s) to select, specified as positional arguments.
            Accepts expression input. Strings are parsed as column names,
            other non-expression inputs are parsed as literals.
        **named_exprs
            Additional columns to select, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        See Also
        --------
        select
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .select_seq(*exprs, **named_exprs)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def with_columns(
        self,
        *exprs: IntoExpr | Iterable[IntoExpr],
        **named_exprs: IntoExpr,
    ) -> DataFrame:
        """
        Add columns to this DataFrame.

        Added columns will replace existing columns with the same name.

        Parameters
        ----------
        *exprs
            Column(s) to add, specified as positional arguments.
            Accepts expression input. Strings are parsed as column names, other
            non-expression inputs are parsed as literals.
        **named_exprs
            Additional columns to add, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        Returns
        -------
        DataFrame
            A new DataFrame with the columns added.

        Notes
        -----
        Creating a new DataFrame using this method does not create a new copy of
        existing data.

        Examples
        --------
        Pass an expression to add it as a new column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [0.5, 4, 10, 13],
        ...         "c": [True, True, False, True],
        ...     }
        ... )
        >>> df.with_columns((pl.col("a") ** 2).alias("a^2"))
        shape: (4, 4)
        ┌─────┬──────┬───────┬─────┐
        │ a   ┆ b    ┆ c     ┆ a^2 │
        │ --- ┆ ---  ┆ ---   ┆ --- │
        │ i64 ┆ f64  ┆ bool  ┆ i64 │
        ╞═════╪══════╪═══════╪═════╡
        │ 1   ┆ 0.5  ┆ true  ┆ 1   │
        │ 2   ┆ 4.0  ┆ true  ┆ 4   │
        │ 3   ┆ 10.0 ┆ false ┆ 9   │
        │ 4   ┆ 13.0 ┆ true  ┆ 16  │
        └─────┴──────┴───────┴─────┘

        Added columns will replace existing columns with the same name.

        >>> df.with_columns(pl.col("a").cast(pl.Float64))
        shape: (4, 3)
        ┌─────┬──────┬───────┐
        │ a   ┆ b    ┆ c     │
        │ --- ┆ ---  ┆ ---   │
        │ f64 ┆ f64  ┆ bool  │
        ╞═════╪══════╪═══════╡
        │ 1.0 ┆ 0.5  ┆ true  │
        │ 2.0 ┆ 4.0  ┆ true  │
        │ 3.0 ┆ 10.0 ┆ false │
        │ 4.0 ┆ 13.0 ┆ true  │
        └─────┴──────┴───────┘

        Multiple columns can be added using positional arguments.

        >>> df.with_columns(
        ...     (pl.col("a") ** 2).alias("a^2"),
        ...     (pl.col("b") / 2).alias("b/2"),
        ...     (pl.col("c").not_()).alias("not c"),
        ... )
        shape: (4, 6)
        ┌─────┬──────┬───────┬─────┬──────┬───────┐
        │ a   ┆ b    ┆ c     ┆ a^2 ┆ b/2  ┆ not c │
        │ --- ┆ ---  ┆ ---   ┆ --- ┆ ---  ┆ ---   │
        │ i64 ┆ f64  ┆ bool  ┆ i64 ┆ f64  ┆ bool  │
        ╞═════╪══════╪═══════╪═════╪══════╪═══════╡
        │ 1   ┆ 0.5  ┆ true  ┆ 1   ┆ 0.25 ┆ false │
        │ 2   ┆ 4.0  ┆ true  ┆ 4   ┆ 2.0  ┆ false │
        │ 3   ┆ 10.0 ┆ false ┆ 9   ┆ 5.0  ┆ true  │
        │ 4   ┆ 13.0 ┆ true  ┆ 16  ┆ 6.5  ┆ false │
        └─────┴──────┴───────┴─────┴──────┴───────┘

        Multiple columns can also be added by passing a list of expressions.

        >>> df.with_columns(
        ...     [
        ...         (pl.col("a") ** 2).alias("a^2"),
        ...         (pl.col("b") / 2).alias("b/2"),
        ...         (pl.col("c").not_()).alias("not c"),
        ...     ]
        ... )
        shape: (4, 6)
        ┌─────┬──────┬───────┬─────┬──────┬───────┐
        │ a   ┆ b    ┆ c     ┆ a^2 ┆ b/2  ┆ not c │
        │ --- ┆ ---  ┆ ---   ┆ --- ┆ ---  ┆ ---   │
        │ i64 ┆ f64  ┆ bool  ┆ i64 ┆ f64  ┆ bool  │
        ╞═════╪══════╪═══════╪═════╪══════╪═══════╡
        │ 1   ┆ 0.5  ┆ true  ┆ 1   ┆ 0.25 ┆ false │
        │ 2   ┆ 4.0  ┆ true  ┆ 4   ┆ 2.0  ┆ false │
        │ 3   ┆ 10.0 ┆ false ┆ 9   ┆ 5.0  ┆ true  │
        │ 4   ┆ 13.0 ┆ true  ┆ 16  ┆ 6.5  ┆ false │
        └─────┴──────┴───────┴─────┴──────┴───────┘

        Use keyword arguments to easily name your expression inputs.

        >>> df.with_columns(
        ...     ab=pl.col("a") * pl.col("b"),
        ...     not_c=pl.col("c").not_(),
        ... )
        shape: (4, 5)
        ┌─────┬──────┬───────┬──────┬───────┐
        │ a   ┆ b    ┆ c     ┆ ab   ┆ not_c │
        │ --- ┆ ---  ┆ ---   ┆ ---  ┆ ---   │
        │ i64 ┆ f64  ┆ bool  ┆ f64  ┆ bool  │
        ╞═════╪══════╪═══════╪══════╪═══════╡
        │ 1   ┆ 0.5  ┆ true  ┆ 0.5  ┆ false │
        │ 2   ┆ 4.0  ┆ true  ┆ 8.0  ┆ false │
        │ 3   ┆ 10.0 ┆ false ┆ 30.0 ┆ true  │
        │ 4   ┆ 13.0 ┆ true  ┆ 52.0 ┆ false │
        └─────┴──────┴───────┴──────┴───────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .with_columns(*exprs, **named_exprs)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def with_columns_seq(
        self,
        *exprs: IntoExpr | Iterable[IntoExpr],
        **named_exprs: IntoExpr,
    ) -> DataFrame:
        """
        Add columns to this DataFrame.

        Added columns will replace existing columns with the same name.

        This will run all expression sequentially instead of in parallel.
        Use this when the work per expression is cheap.

        Parameters
        ----------
        *exprs
            Column(s) to add, specified as positional arguments.
            Accepts expression input. Strings are parsed as column names, other
            non-expression inputs are parsed as literals.
        **named_exprs
            Additional columns to add, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        Returns
        -------
        DataFrame
            A new DataFrame with the columns added.

        See Also
        --------
        with_columns
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .with_columns_seq(*exprs, **named_exprs)
            .collect(optimizations=QueryOptFlags._eager())
        )

    @overload
    def n_chunks(self, strategy: Literal["first"] = ...) -> int: ...

    @overload
    def n_chunks(self, strategy: Literal["all"]) -> list[int]: ...

    def n_chunks(self, strategy: Literal["first", "all"] = "first") -> int | list[int]:
        """
        Get number of chunks used by the ChunkedArrays of this DataFrame.

        Parameters
        ----------
        strategy : {'first', 'all'}
            Return the number of chunks of the 'first' column,
            or 'all' columns in this DataFrame.


        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [0.5, 4, 10, 13],
        ...         "c": [True, True, False, True],
        ...     }
        ... )
        >>> df.n_chunks()
        1
        >>> df.n_chunks(strategy="all")
        [1, 1, 1]
        """
        if strategy == "first":
            return self._df.n_chunks()
        elif strategy == "all":
            return [s.n_chunks() for s in self.__iter__()]
        else:
            msg = (
                f"unexpected input for `strategy`: {strategy!r}"
                f"\n\nChoose one of {{'first', 'all'}}"
            )
            raise ValueError(msg)

    def max(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their maximum value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.max()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().max().collect(optimizations=QueryOptFlags._eager())

    def max_horizontal(self) -> Series:
        """
        Get the maximum value horizontally across columns.

        Returns
        -------
        Series
            A Series named `"max"`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4.0, 5.0, 6.0],
        ...     }
        ... )
        >>> df.max_horizontal()
        shape: (3,)
        Series: 'max' [f64]
        [
                4.0
                5.0
                6.0
        ]
        """
        return self.select(max=F.max_horizontal(F.all())).to_series()

    def min(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their minimum value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.min()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().min().collect(optimizations=QueryOptFlags._eager())

    def min_horizontal(self) -> Series:
        """
        Get the minimum value horizontally across columns.

        Returns
        -------
        Series
            A Series named `"min"`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4.0, 5.0, 6.0],
        ...     }
        ... )
        >>> df.min_horizontal()
        shape: (3,)
        Series: 'min' [f64]
        [
                1.0
                2.0
                3.0
        ]
        """
        return self.select(min=F.min_horizontal(F.all())).to_series()

    def sum(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their sum value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.sum()
        shape: (1, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ i64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 6   ┆ 21  ┆ null │
        └─────┴─────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().sum().collect(optimizations=QueryOptFlags._eager())

    def sum_horizontal(self, *, ignore_nulls: bool = True) -> Series:
        """
        Sum all values horizontally across columns.

        Parameters
        ----------
        ignore_nulls
            Ignore null values (default).
            If set to `False`, any null value in the input will lead to a null output.

        Returns
        -------
        Series
            A Series named `"sum"`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4.0, 5.0, 6.0],
        ...     }
        ... )
        >>> df.sum_horizontal()
        shape: (3,)
        Series: 'sum' [f64]
        [
                5.0
                7.0
                9.0
        ]
        """
        return self.select(
            sum=F.sum_horizontal(F.all(), ignore_nulls=ignore_nulls)
        ).to_series()

    def mean(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their mean value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...         "spam": [True, False, None],
        ...     }
        ... )
        >>> df.mean()
        shape: (1, 4)
        ┌─────┬─────┬──────┬──────┐
        │ foo ┆ bar ┆ ham  ┆ spam │
        │ --- ┆ --- ┆ ---  ┆ ---  │
        │ f64 ┆ f64 ┆ str  ┆ f64  │
        ╞═════╪═════╪══════╪══════╡
        │ 2.0 ┆ 7.0 ┆ null ┆ 0.5  │
        └─────┴─────┴──────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().mean().collect(optimizations=QueryOptFlags._eager())

    def mean_horizontal(self, *, ignore_nulls: bool = True) -> Series:
        """
        Take the mean of all values horizontally across columns.

        Parameters
        ----------
        ignore_nulls
            Ignore null values (default).
            If set to `False`, any null value in the input will lead to a null output.

        Returns
        -------
        Series
            A Series named `"mean"`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4.0, 5.0, 6.0],
        ...     }
        ... )
        >>> df.mean_horizontal()
        shape: (3,)
        Series: 'mean' [f64]
        [
                2.5
                3.5
                4.5
        ]
        """
        return self.select(
            mean=F.mean_horizontal(F.all(), ignore_nulls=ignore_nulls)
        ).to_series()

    def std(self, ddof: int = 1) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their standard deviation value.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.std()
        shape: (1, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ f64 ┆ f64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 1.0 ┆ 1.0 ┆ null │
        └─────┴─────┴──────┘
        >>> df.std(ddof=0)
        shape: (1, 3)
        ┌──────────┬──────────┬──────┐
        │ foo      ┆ bar      ┆ ham  │
        │ ---      ┆ ---      ┆ ---  │
        │ f64      ┆ f64      ┆ str  │
        ╞══════════╪══════════╪══════╡
        │ 0.816497 ┆ 0.816497 ┆ null │
        └──────────┴──────────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().std(ddof).collect(optimizations=QueryOptFlags._eager())

    def var(self, ddof: int = 1) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their variance value.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.var()
        shape: (1, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ f64 ┆ f64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 1.0 ┆ 1.0 ┆ null │
        └─────┴─────┴──────┘
        >>> df.var(ddof=0)
        shape: (1, 3)
        ┌──────────┬──────────┬──────┐
        │ foo      ┆ bar      ┆ ham  │
        │ ---      ┆ ---      ┆ ---  │
        │ f64      ┆ f64      ┆ str  │
        ╞══════════╪══════════╪══════╡
        │ 0.666667 ┆ 0.666667 ┆ null │
        └──────────┴──────────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().var(ddof).collect(optimizations=QueryOptFlags._eager())

    def median(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their median value.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.median()
        shape: (1, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ f64 ┆ f64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 2.0 ┆ 7.0 ┆ null │
        └─────┴─────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().median().collect(optimizations=QueryOptFlags._eager())

    def product(self) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their product values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": [0.5, 4, 10],
        ...         "c": [True, True, False],
        ...     }
        ... )

        >>> df.product()
        shape: (1, 3)
        ┌─────┬──────┬─────┐
        │ a   ┆ b    ┆ c   │
        │ --- ┆ ---  ┆ --- │
        │ i64 ┆ f64  ┆ i64 │
        ╞═════╪══════╪═════╡
        │ 6   ┆ 20.0 ┆ 0   │
        └─────┴──────┴─────┘
        """
        exprs = []
        for name, dt in self.schema.items():
            if dt.is_numeric() or isinstance(dt, Boolean):
                exprs.append(F.col(name).product())
            else:
                exprs.append(F.lit(None).alias(name))

        return self.select(exprs)

    def quantile(
        self, quantile: float, interpolation: QuantileMethod = "nearest"
    ) -> DataFrame:
        """
        Aggregate the columns of this DataFrame to their quantile value.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.quantile(0.5, "nearest")
        shape: (1, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ f64 ┆ f64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 2.0 ┆ 7.0 ┆ null │
        └─────┴─────┴──────┘
        """  # noqa: W505
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .quantile(quantile, interpolation)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def to_dummies(
        self,
        columns: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        *,
        separator: str = "_",
        drop_first: bool = False,
        drop_nulls: bool = False,
    ) -> DataFrame:
        """
        Convert categorical variables into dummy/indicator variables.

        Parameters
        ----------
        columns
            Column name(s) or selector(s) that should be converted to dummy
            variables. If set to `None` (default), convert all columns.
        separator
            Separator/delimiter used when generating column names.
        drop_first
            Remove the first category from the variables being encoded.
        drop_nulls
            If there are `None` values in the series, a `null` column is not generated

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2],
        ...         "bar": [3, 4],
        ...         "ham": ["a", "b"],
        ...     }
        ... )
        >>> df.to_dummies()
        shape: (2, 6)
        ┌───────┬───────┬───────┬───────┬───────┬───────┐
        │ foo_1 ┆ foo_2 ┆ bar_3 ┆ bar_4 ┆ ham_a ┆ ham_b │
        │ ---   ┆ ---   ┆ ---   ┆ ---   ┆ ---   ┆ ---   │
        │ u8    ┆ u8    ┆ u8    ┆ u8    ┆ u8    ┆ u8    │
        ╞═══════╪═══════╪═══════╪═══════╪═══════╪═══════╡
        │ 1     ┆ 0     ┆ 1     ┆ 0     ┆ 1     ┆ 0     │
        │ 0     ┆ 1     ┆ 0     ┆ 1     ┆ 0     ┆ 1     │
        └───────┴───────┴───────┴───────┴───────┴───────┘

        >>> df.to_dummies(drop_first=True)
        shape: (2, 3)
        ┌───────┬───────┬───────┐
        │ foo_2 ┆ bar_4 ┆ ham_b │
        │ ---   ┆ ---   ┆ ---   │
        │ u8    ┆ u8    ┆ u8    │
        ╞═══════╪═══════╪═══════╡
        │ 0     ┆ 0     ┆ 0     │
        │ 1     ┆ 1     ┆ 1     │
        └───────┴───────┴───────┘

        >>> import polars.selectors as cs
        >>> df.to_dummies(cs.integer(), separator=":")
        shape: (2, 5)
        ┌───────┬───────┬───────┬───────┬─────┐
        │ foo:1 ┆ foo:2 ┆ bar:3 ┆ bar:4 ┆ ham │
        │ ---   ┆ ---   ┆ ---   ┆ ---   ┆ --- │
        │ u8    ┆ u8    ┆ u8    ┆ u8    ┆ str │
        ╞═══════╪═══════╪═══════╪═══════╪═════╡
        │ 1     ┆ 0     ┆ 1     ┆ 0     ┆ a   │
        │ 0     ┆ 1     ┆ 0     ┆ 1     ┆ b   │
        └───────┴───────┴───────┴───────┴─────┘

        >>> df.to_dummies(cs.integer(), drop_first=True, separator=":")
        shape: (2, 3)
        ┌───────┬───────┬─────┐
        │ foo:2 ┆ bar:4 ┆ ham │
        │ ---   ┆ ---   ┆ --- │
        │ u8    ┆ u8    ┆ str │
        ╞═══════╪═══════╪═════╡
        │ 0     ┆ 0     ┆ a   │
        │ 1     ┆ 1     ┆ b   │
        └───────┴───────┴─────┘
        """
        if columns is not None:
            columns = _expand_selectors(self, columns)
        return self._from_pydf(
            self._df.to_dummies(columns, separator, drop_first, drop_nulls)
        )

    def unique(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
        *,
        keep: UniqueKeepStrategy = "any",
        maintain_order: bool = False,
    ) -> DataFrame:
        """
        Drop duplicate rows from this dataframe.

        Parameters
        ----------
        subset
            Column name(s) or selector(s), to consider when identifying
            duplicate rows. If set to `None` (default), use all columns.
        keep : {'first', 'last', 'any', 'none'}
            Which of the duplicate rows to keep.

            * 'any': Does not give any guarantee of which row is kept.
                     This allows more optimizations.
            * 'none': Don't keep duplicate rows.
            * 'first': Keep first unique row.
            * 'last': Keep last unique row.
        maintain_order
            Keep the same order as the original DataFrame. This is more expensive to
            compute.
            Settings this to `True` blocks the possibility
            to run on the streaming engine.

        Returns
        -------
        DataFrame
            DataFrame with unique rows.

        Warnings
        --------
        This method will fail if there is a column of type `List` in the DataFrame or
        subset.

        Notes
        -----
        If you're coming from pandas, this is similar to
        `pandas.DataFrame.drop_duplicates`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3, 1],
        ...         "bar": ["a", "a", "a", "a"],
        ...         "ham": ["b", "b", "b", "b"],
        ...     }
        ... )
        >>> df.unique(maintain_order=True)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ str ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ a   ┆ b   │
        │ 2   ┆ a   ┆ b   │
        │ 3   ┆ a   ┆ b   │
        └─────┴─────┴─────┘
        >>> df.unique(subset=["bar", "ham"], maintain_order=True)
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ str ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ a   ┆ b   │
        └─────┴─────┴─────┘
        >>> df.unique(keep="last", maintain_order=True)
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ str ┆ str │
        ╞═════╪═════╪═════╡
        │ 2   ┆ a   ┆ b   │
        │ 3   ┆ a   ┆ b   │
        │ 1   ┆ a   ┆ b   │
        └─────┴─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .unique(subset=subset, keep=keep, maintain_order=maintain_order)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def n_unique(self, subset: str | Expr | Sequence[str | Expr] | None = None) -> int:
        """
        Return the number of unique rows, or the number of unique row-subsets.

        Parameters
        ----------
        subset
            One or more columns/expressions that define what to count;
            omit to return the count of unique rows.

        Notes
        -----
        This method operates at the `DataFrame` level; to operate on subsets at the
        expression level you can make use of struct-packing instead, for example:

        >>> expr_unique_subset = pl.struct("a", "b").n_unique()

        If instead you want to count the number of unique values per-column, you can
        also use expression-level syntax to return a new frame containing that result:

        >>> df = pl.DataFrame(
        ...     [[1, 2, 3], [1, 2, 4]], schema=["a", "b", "c"], orient="row"
        ... )
        >>> df_nunique = df.select(pl.all().n_unique())

        In aggregate context there is also an equivalent method for returning the
        unique values per-group:

        >>> df_agg_nunique = df.group_by("a").n_unique()

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 1, 2, 3, 4, 5],
        ...         "b": [0.5, 0.5, 1.0, 2.0, 3.0, 3.0],
        ...         "c": [True, True, True, False, True, True],
        ...     }
        ... )
        >>> df.n_unique()
        5

        Simple columns subset.

        >>> df.n_unique(subset=["b", "c"])
        4

        Expression subset.

        >>> df.n_unique(
        ...     subset=[
        ...         (pl.col("a") // 2),
        ...         (pl.col("c") | (pl.col("b") >= 2)),
        ...     ],
        ... )
        3
        """
        if isinstance(subset, str):
            expr = F.col(subset)
        elif isinstance(subset, pl.Expr):
            expr = subset
        elif isinstance(subset, Sequence) and len(subset) == 1:
            expr = wrap_expr(parse_into_expression(subset[0]))
        else:
            struct_fields = F.all() if (subset is None) else subset
            expr = F.struct(struct_fields)

        from polars.lazyframe.opt_flags import QueryOptFlags

        df = (
            self.lazy()
            .select(expr.n_unique())
            .collect(optimizations=QueryOptFlags._eager())
        )
        return 0 if df.is_empty() else df.row(0)[0]

    @deprecated(
        "`DataFrame.approx_n_unique` is deprecated; "
        "use `select(pl.all().approx_n_unique())` instead."
    )
    def approx_n_unique(self) -> DataFrame:
        """
        Approximate count of unique values.

        .. deprecated:: 0.20.11
            Use the `select(pl.all().approx_n_unique())` method instead.

        This is done using the HyperLogLog++ algorithm for cardinality estimation.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> df.approx_n_unique()  # doctest: +SKIP
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ u32 ┆ u32 │
        ╞═════╪═════╡
        │ 4   ┆ 2   │
        └─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy().approx_n_unique().collect(optimizations=QueryOptFlags._eager())
        )

    def rechunk(self) -> DataFrame:
        """
        Rechunk the data in this DataFrame to a contiguous allocation.

        This will make sure all subsequent operations have optimal and predictable
        performance.
        """
        return self._from_pydf(self._df.rechunk())

    def null_count(self) -> DataFrame:
        """
        Create a new DataFrame that shows the null counts per column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, None, 3],
        ...         "bar": [6, 7, None],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.null_count()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ u32 ┆ u32 ┆ u32 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 1   ┆ 0   │
        └─────┴─────┴─────┘
        """
        return self._from_pydf(self._df.null_count())

    def sample(
        self,
        n: int | Series | None = None,
        *,
        fraction: float | Series | None = None,
        with_replacement: bool = False,
        shuffle: bool = False,
        seed: int | None = None,
    ) -> DataFrame:
        """
        Sample from this DataFrame.

        Parameters
        ----------
        n
            Number of items to return. Cannot be used with `fraction`. Defaults to 1 if
            `fraction` is None.
        fraction
            Fraction of items to return. Cannot be used with `n`.
        with_replacement
            Allow values to be sampled more than once.
        shuffle
            If set to True, the order of the sampled rows will be shuffled. If
            set to False (default), the order of the returned rows will be
            neither stable nor fully random.
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated for each sample operation.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.sample(n=2, seed=0)  # doctest: +IGNORE_RESULT
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 3   ┆ 8   ┆ c   │
        │ 2   ┆ 7   ┆ b   │
        └─────┴─────┴─────┘
        """
        if n is not None and fraction is not None:
            msg = "cannot specify both `n` and `fraction`"
            raise ValueError(msg)

        if seed is None:
            seed = random.randint(0, 10000)

        if n is None and fraction is not None:
            if not isinstance(fraction, pl.Series):
                fraction = pl.Series("frac", [fraction])

            return self._from_pydf(
                self._df.sample_frac(fraction._s, with_replacement, shuffle, seed)
            )

        if n is None:
            n = 1

        if not isinstance(n, pl.Series):
            n = pl.Series("", [n])

        return self._from_pydf(self._df.sample_n(n._s, with_replacement, shuffle, seed))

    def fold(self, operation: Callable[[Series, Series], Series]) -> Series:
        """
        Apply a horizontal reduction on a DataFrame.

        This can be used to effectively determine aggregations on a row level, and can
        be applied to any DataType that can be supercast (cast to a similar parent
        type).

        An example of the supercast rules when applying an arithmetic operation on two
        DataTypes are for instance:

        - Int8 + String = String
        - Float32 + Int64 = Float32
        - Float32 + Float64 = Float64

        Examples
        --------
        A horizontal sum operation:

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [2, 1, 3],
        ...         "b": [1, 2, 3],
        ...         "c": [1.0, 2.0, 3.0],
        ...     }
        ... )
        >>> df.fold(lambda s1, s2: s1 + s2)
        shape: (3,)
        Series: 'a' [f64]
        [
            4.0
            5.0
            9.0
        ]

        A horizontal minimum operation:

        >>> df = pl.DataFrame({"a": [2, 1, 3], "b": [1, 2, 3], "c": [1.0, 2.0, 3.0]})
        >>> df.fold(lambda s1, s2: s1.zip_with(s1 < s2, s2))
        shape: (3,)
        Series: 'a' [f64]
        [
            1.0
            1.0
            3.0
        ]

        A horizontal string concatenation:

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": ["foo", "bar", None],
        ...         "b": [1, 2, 3],
        ...         "c": [1.0, 2.0, 3.0],
        ...     }
        ... )
        >>> df.fold(lambda s1, s2: s1 + s2)
        shape: (3,)
        Series: 'a' [str]
        [
            "foo11.0"
            "bar22.0"
            null
        ]

        A horizontal boolean or, similar to a row-wise .any():

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [False, False, True],
        ...         "b": [False, True, False],
        ...     }
        ... )
        >>> df.fold(lambda s1, s2: s1 | s2)
        shape: (3,)
        Series: 'a' [bool]
        [
                false
                true
                true
        ]

        Parameters
        ----------
        operation
            function that takes two `Series` and returns a `Series`.
        """
        acc = self.to_series(0)

        for i in range(1, self.width):
            acc = operation(acc, self.to_series(i))
        return acc

    @overload
    def row(
        self,
        index: int | None = ...,
        *,
        by_predicate: Expr | None = ...,
        named: Literal[False] = ...,
    ) -> tuple[Any, ...]: ...

    @overload
    def row(
        self,
        index: int | None = ...,
        *,
        by_predicate: Expr | None = ...,
        named: Literal[True],
    ) -> dict[str, Any]: ...

    def row(
        self,
        index: int | None = None,
        *,
        by_predicate: Expr | None = None,
        named: bool = False,
    ) -> tuple[Any, ...] | dict[str, Any]:
        """
        Get the values of a single row, either by index or by predicate.

        Parameters
        ----------
        index
            Row index.
        by_predicate
            Select the row according to a given expression/predicate.
        named
            Return a dictionary instead of a tuple. The dictionary is a mapping of
            column name to row value. This is more expensive than returning a regular
            tuple, but allows for accessing values by column name.

        Returns
        -------
        tuple (default) or dictionary of row values

        Notes
        -----
        The `index` and `by_predicate` params are mutually exclusive. Additionally,
        to ensure clarity, the `by_predicate` parameter must be supplied by keyword.

        When using `by_predicate` it is an error condition if anything other than
        one row is returned; more than one row raises `TooManyRowsReturnedError`, and
        zero rows will raise `NoRowsReturnedError` (both inherit from `RowsError`).

        Warnings
        --------
        You should NEVER use this method to iterate over a DataFrame; if you require
        row-iteration you should strongly prefer use of `iter_rows()` instead.

        See Also
        --------
        iter_rows : Row iterator over frame data (does not materialise all rows).
        rows : Materialise all frame data as a list of rows (potentially expensive).
        item: Return dataframe element as a scalar.

        Examples
        --------
        Specify an index to return the row at the given index as a tuple.

        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> df.row(2)
        (3, 8, 'c')

        Specify `named=True` to get a dictionary instead with a mapping of column
        names to row values.

        >>> df.row(2, named=True)
        {'foo': 3, 'bar': 8, 'ham': 'c'}

        Use `by_predicate` to return the row that matches the given predicate.

        >>> df.row(by_predicate=(pl.col("ham") == "b"))
        (2, 7, 'b')
        """
        if index is not None and by_predicate is not None:
            msg = "cannot set both 'index' and 'by_predicate'; mutually exclusive"
            raise ValueError(msg)
        elif isinstance(index, pl.Expr):
            msg = "expressions should be passed to the `by_predicate` parameter"
            raise TypeError(msg)

        if index is not None:
            row = self._df.row_tuple(index)
            if named:
                return dict(zip(self.columns, row))
            else:
                return row

        elif by_predicate is not None:
            if not isinstance(by_predicate, pl.Expr):
                msg = f"expected `by_predicate` to be an expression, got {qualified_type_name(by_predicate)!r}"
                raise TypeError(msg)
            rows = self.filter(by_predicate).rows()
            n_rows = len(rows)
            if n_rows > 1:
                msg = f"predicate <{by_predicate!s}> returned {n_rows} rows"
                raise TooManyRowsReturnedError(msg)
            elif n_rows == 0:
                msg = f"predicate <{by_predicate!s}> returned no rows"
                raise NoRowsReturnedError(msg)

            row = rows[0]
            if named:
                return dict(zip(self.columns, row))
            else:
                return row
        else:
            msg = "one of `index` or `by_predicate` must be set"
            raise ValueError(msg)

    @overload
    def rows(self, *, named: Literal[False] = ...) -> list[tuple[Any, ...]]: ...

    @overload
    def rows(self, *, named: Literal[True]) -> list[dict[str, Any]]: ...

    def rows(
        self, *, named: bool = False
    ) -> list[tuple[Any, ...]] | list[dict[str, Any]]:
        """
        Returns all data in the DataFrame as a list of rows of python-native values.

        By default, each row is returned as a tuple of values given in the same order
        as the frame columns. Setting `named=True` will return rows of dictionaries
        instead.

        Parameters
        ----------
        named
            Return dictionaries instead of tuples. The dictionaries are a mapping of
            column name to row value. This is more expensive than returning a regular
            tuple, but allows for accessing values by column name.

        Notes
        -----
        If you have `ns`-precision temporal values you should be aware that Python
        natively only supports up to `μs`-precision; `ns`-precision values will be
        truncated to microseconds on conversion to Python. If this matters to your
        use-case you should export to a different format (such as Arrow or NumPy).

        Warnings
        --------
        Row-iteration is not optimal as the underlying data is stored in columnar form;
        where possible, prefer export via one of the dedicated export/output methods.
        You should also consider using `iter_rows` instead, to avoid materialising all
        the data at once; there is little performance difference between the two, but
        peak memory can be reduced if processing rows in batches.

        Returns
        -------
        list of row value tuples (default), or list of dictionaries (if `named=True`).

        See Also
        --------
        iter_rows : Row iterator over frame data (does not materialise all rows).
        rows_by_key : Materialises frame data as a key-indexed dictionary.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "x": ["a", "b", "b", "a"],
        ...         "y": [1, 2, 3, 4],
        ...         "z": [0, 3, 6, 9],
        ...     }
        ... )
        >>> df.rows()
        [('a', 1, 0), ('b', 2, 3), ('b', 3, 6), ('a', 4, 9)]
        >>> df.rows(named=True)
        [{'x': 'a', 'y': 1, 'z': 0},
         {'x': 'b', 'y': 2, 'z': 3},
         {'x': 'b', 'y': 3, 'z': 6},
         {'x': 'a', 'y': 4, 'z': 9}]
        """
        if named:
            # Load these into the local namespace for a minor performance boost
            dict_, zip_, columns = dict, zip, self.columns
            return [dict_(zip_(columns, row)) for row in self._df.row_tuples()]
        else:
            return self._df.row_tuples()

    @overload
    def rows_by_key(
        self,
        key: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        named: Literal[False] = ...,
        include_key: bool = ...,
        unique: Literal[False] = ...,
    ) -> dict[Any, list[Any]]: ...

    @overload
    def rows_by_key(
        self,
        key: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        named: Literal[False] = ...,
        include_key: bool = ...,
        unique: Literal[True],
    ) -> dict[Any, Any]: ...

    @overload
    def rows_by_key(
        self,
        key: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        named: Literal[True],
        include_key: bool = ...,
        unique: Literal[False] = ...,
    ) -> dict[Any, list[dict[str, Any]]]: ...

    @overload
    def rows_by_key(
        self,
        key: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        named: Literal[True],
        include_key: bool = ...,
        unique: Literal[True],
    ) -> dict[Any, dict[str, Any]]: ...

    def rows_by_key(
        self,
        key: ColumnNameOrSelector | Sequence[ColumnNameOrSelector],
        *,
        named: bool = False,
        include_key: bool = False,
        unique: bool = False,
    ) -> dict[Any, Any]:
        """
        Returns all data as a dictionary of python-native values keyed by some column.

        This method is like `rows`, but instead of returning rows in a flat list, rows
        are grouped by the values in the `key` column(s) and returned as a dictionary.

        Note that this method should not be used in place of native operations, due to
        the high cost of materializing all frame data out into a dictionary; it should
        be used only when you need to move the values out into a Python data structure
        or other object that cannot operate directly with Polars/Arrow.

        Parameters
        ----------
        key
            The column(s) to use as the key for the returned dictionary. If multiple
            columns are specified, the key will be a tuple of those values, otherwise
            it will be a string.
        named
            Return dictionary rows instead of tuples, mapping column name to row value.
        include_key
            Include key values inline with the associated data (by default the key
            values are omitted as a memory/performance optimisation, as they can be
            reoconstructed from the key).
        unique
            Indicate that the key is unique; this will result in a 1:1 mapping from
            key to a single associated row. Note that if the key is *not* actually
            unique the last row with the given key will be returned.

        Notes
        -----
        If you have `ns`-precision temporal values you should be aware that Python
        natively only supports up to `μs`-precision; `ns`-precision values will be
        truncated to microseconds on conversion to Python. If this matters to your
        use-case you should export to a different format (such as Arrow or NumPy).

        See Also
        --------
        rows : Materialize all frame data as a list of rows (potentially expensive).
        iter_rows : Row iterator over frame data (does not materialize all rows).
        to_dict : Convert DataFrame to a dictionary mapping column name to values.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "w": ["a", "b", "b", "a"],
        ...         "x": ["q", "q", "q", "k"],
        ...         "y": [1.0, 2.5, 3.0, 4.5],
        ...         "z": [9, 8, 7, 6],
        ...     }
        ... )

        Group rows by the given key column(s):

        >>> df.rows_by_key(key=["w"])
        defaultdict(<class 'list'>,
            {'a': [('q', 1.0, 9), ('k', 4.5, 6)],
             'b': [('q', 2.5, 8), ('q', 3.0, 7)]})

        Return the same row groupings as dictionaries:

        >>> df.rows_by_key(key=["w"], named=True)
        defaultdict(<class 'list'>,
            {'a': [{'x': 'q', 'y': 1.0, 'z': 9},
                   {'x': 'k', 'y': 4.5, 'z': 6}],
             'b': [{'x': 'q', 'y': 2.5, 'z': 8},
                   {'x': 'q', 'y': 3.0, 'z': 7}]})

        Return row groupings, assuming keys are unique:

        >>> df.rows_by_key(key=["z"], unique=True)
        {9: ('a', 'q', 1.0),
         8: ('b', 'q', 2.5),
         7: ('b', 'q', 3.0),
         6: ('a', 'k', 4.5)}

        Return row groupings as dictionaries, assuming keys are unique:

        >>> df.rows_by_key(key=["z"], named=True, unique=True)
        {9: {'w': 'a', 'x': 'q', 'y': 1.0},
         8: {'w': 'b', 'x': 'q', 'y': 2.5},
         7: {'w': 'b', 'x': 'q', 'y': 3.0},
         6: {'w': 'a', 'x': 'k', 'y': 4.5}}

        Return dictionary rows grouped by a compound key, including key values:

        >>> df.rows_by_key(key=["w", "x"], named=True, include_key=True)
        defaultdict(<class 'list'>,
            {('a', 'q'): [{'w': 'a', 'x': 'q', 'y': 1.0, 'z': 9}],
             ('b', 'q'): [{'w': 'b', 'x': 'q', 'y': 2.5, 'z': 8},
                          {'w': 'b', 'x': 'q', 'y': 3.0, 'z': 7}],
             ('a', 'k'): [{'w': 'a', 'x': 'k', 'y': 4.5, 'z': 6}]})
        """
        key = _expand_selectors(self, key)

        keys = (
            iter(self.get_column(key[0]))
            if len(key) == 1
            else self.select(key).iter_rows()
        )

        if include_key:
            values = self
        else:
            data_cols = [k for k in self.schema if k not in key]
            values = self.select(data_cols)

        zipped = zip(keys, values.iter_rows(named=named))  # type: ignore[call-overload]

        # if unique, we expect to write just one entry per key; otherwise, we're
        # returning a list of rows for each key, so append into a defaultdict.
        if unique:
            rows = dict(zipped)
        else:
            rows = defaultdict(list)
            for key, data in zipped:
                rows[key].append(data)

        return rows

    @overload
    def iter_rows(
        self, *, named: Literal[False] = ..., buffer_size: int = ...
    ) -> Iterator[tuple[Any, ...]]: ...

    @overload
    def iter_rows(
        self, *, named: Literal[True], buffer_size: int = ...
    ) -> Iterator[dict[str, Any]]: ...

    def iter_rows(
        self, *, named: bool = False, buffer_size: int = 512
    ) -> Iterator[tuple[Any, ...]] | Iterator[dict[str, Any]]:
        """
        Returns an iterator over the DataFrame of rows of python-native values.

        Parameters
        ----------
        named
            Return dictionaries instead of tuples. The dictionaries are a mapping of
            column name to row value. This is more expensive than returning a regular
            tuple, but allows for accessing values by column name.
        buffer_size
            Determines the number of rows that are buffered internally while iterating
            over the data; you should only modify this in very specific cases where the
            default value is determined not to be a good fit to your access pattern, as
            the speedup from using the buffer is significant (~2-4x). Setting this
            value to zero disables row buffering (not recommended).

        Notes
        -----
        If you have `ns`-precision temporal values you should be aware that Python
        natively only supports up to `μs`-precision; `ns`-precision values will be
        truncated to microseconds on conversion to Python. If this matters to your
        use-case you should export to a different format (such as Arrow or NumPy).

        Warnings
        --------
        Row iteration is not optimal as the underlying data is stored in columnar form;
        where possible, prefer export via one of the dedicated export/output methods
        that deals with columnar data.

        Returns
        -------
        iterator of tuples (default) or dictionaries (if named) of python row values

        See Also
        --------
        rows : Materialises all frame data as a list of rows (potentially expensive).
        rows_by_key : Materialises frame data as a key-indexed dictionary.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> [row[0] for row in df.iter_rows()]
        [1, 3, 5]
        >>> [row["b"] for row in df.iter_rows(named=True)]
        [2, 4, 6]
        """
        # load into the local namespace for a (minor) performance boost in the hot loops
        columns, get_row, dict_, zip_ = self.columns, self.row, dict, zip
        has_object = Object in self.dtypes

        # note: buffering rows results in a 2-4x speedup over individual calls
        # to ".row(i)", so it should only be disabled in extremely specific cases.
        if buffer_size and not has_object:
            for offset in range(0, self.height, buffer_size):
                zerocopy_slice = self.slice(offset, buffer_size)
                if named:
                    for row in zerocopy_slice.rows(named=False):
                        yield dict_(zip_(columns, row))
                else:
                    yield from zerocopy_slice.rows(named=False)
        elif named:
            for i in range(self.height):
                yield dict_(zip_(columns, get_row(i)))
        else:
            for i in range(self.height):
                yield get_row(i)

    def iter_columns(self) -> Iterator[Series]:
        """
        Returns an iterator over the columns of this DataFrame.

        Yields
        ------
        Series

        Notes
        -----
        Consider whether you can use :func:`all` instead.
        If you can, it will be more efficient.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> [s.name for s in df.iter_columns()]
        ['a', 'b']

        If you're using this to modify a dataframe's columns, e.g.

        >>> # Do NOT do this
        >>> pl.DataFrame(column * 2 for column in df.iter_columns())
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 2   ┆ 4   │
        │ 6   ┆ 8   │
        │ 10  ┆ 12  │
        └─────┴─────┘

        then consider whether you can use :func:`all` instead:

        >>> df.select(pl.all() * 2)
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 2   ┆ 4   │
        │ 6   ┆ 8   │
        │ 10  ┆ 12  │
        └─────┴─────┘
        """
        for s in self._df.get_columns():
            yield wrap_s(s)

    def iter_slices(self, n_rows: int = 10_000) -> Iterator[DataFrame]:
        r"""
        Returns a non-copying iterator of slices over the underlying DataFrame.

        Parameters
        ----------
        n_rows
            Determines the number of rows contained in each DataFrame slice.

        Examples
        --------
        >>> from datetime import date
        >>> df = pl.DataFrame(
        ...     data={
        ...         "a": range(17_500),
        ...         "b": date(2023, 1, 1),
        ...         "c": "klmnoopqrstuvwxyz",
        ...     },
        ...     schema_overrides={"a": pl.Int32},
        ... )
        >>> for idx, frame in enumerate(df.iter_slices()):
        ...     print(f"{type(frame).__name__}:[{idx}]:{len(frame)}")
        DataFrame:[0]:10000
        DataFrame:[1]:7500

        Using `iter_slices` is an efficient way to chunk-iterate over DataFrames and
        any supported frame export/conversion types; for example, as RecordBatches:

        >>> for frame in df.iter_slices(n_rows=15_000):
        ...     record_batch = frame.to_arrow().to_batches()[0]
        ...     print(f"{record_batch.schema}\n<< {len(record_batch)}")
        a: int32
        b: date32[day]
        c: large_string
        << 15000
        a: int32
        b: date32[day]
        c: large_string
        << 2500

        See Also
        --------
        iter_rows : Row iterator over frame data (does not materialise all rows).
        partition_by : Split into multiple DataFrames, partitioned by groups.
        """
        for offset in range(0, self.height, n_rows):
            yield self.slice(offset, n_rows)

    def shrink_to_fit(self, *, in_place: bool = False) -> DataFrame:
        """
        Shrink DataFrame memory usage.

        Shrinks to fit the exact capacity needed to hold the data.
        """
        if in_place:
            self._df.shrink_to_fit()
            return self
        else:
            df = self.clone()
            df._df.shrink_to_fit()
            return df

    def gather_every(self, n: int, offset: int = 0) -> DataFrame:
        """
        Take every nth row in the DataFrame and return as a new DataFrame.

        Parameters
        ----------
        n
            Gather every *n*-th row.
        offset
            Starting index.

        Examples
        --------
        >>> s = pl.DataFrame({"a": [1, 2, 3, 4], "b": [5, 6, 7, 8]})
        >>> s.gather_every(2)
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 5   │
        │ 3   ┆ 7   │
        └─────┴─────┘

        >>> s.gather_every(2, offset=1)
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 2   ┆ 6   │
        │ 4   ┆ 8   │
        └─────┴─────┘
        """
        return self.select(F.col("*").gather_every(n, offset))

    def hash_rows(
        self,
        seed: int = 0,
        seed_1: int | None = None,
        seed_2: int | None = None,
        seed_3: int | None = None,
    ) -> Series:
        """
        Hash and combine the rows in this DataFrame.

        The hash value is of type `UInt64`.

        Parameters
        ----------
        seed
            Random seed parameter. Defaults to 0.
        seed_1
            Random seed parameter. Defaults to `seed` if not set.
        seed_2
            Random seed parameter. Defaults to `seed` if not set.
        seed_3
            Random seed parameter. Defaults to `seed` if not set.

        Notes
        -----
        This implementation of `hash_rows` does not guarantee stable results
        across different Polars versions. Its stability is only guaranteed within a
        single version.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, None, 3, 4],
        ...         "ham": ["a", "b", None, "d"],
        ...     }
        ... )
        >>> df.hash_rows(seed=42)  # doctest: +IGNORE_RESULT
        shape: (4,)
        Series: '' [u64]
        [
            10783150408545073287
            1438741209321515184
            10047419486152048166
            2047317070637311557
        ]
        """
        k0 = seed
        k1 = seed_1 if seed_1 is not None else seed
        k2 = seed_2 if seed_2 is not None else seed
        k3 = seed_3 if seed_3 is not None else seed
        return wrap_s(self._df.hash_rows(k0, k1, k2, k3))

    def interpolate(self) -> DataFrame:
        """
        Interpolate intermediate values. The interpolation method is linear.

        Nulls at the beginning and end of the series remain null.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, None, 9, 10],
        ...         "bar": [6, 7, 9, None],
        ...         "baz": [1, None, None, 9],
        ...     }
        ... )
        >>> df.interpolate()
        shape: (4, 3)
        ┌──────┬──────┬──────────┐
        │ foo  ┆ bar  ┆ baz      │
        │ ---  ┆ ---  ┆ ---      │
        │ f64  ┆ f64  ┆ f64      │
        ╞══════╪══════╪══════════╡
        │ 1.0  ┆ 6.0  ┆ 1.0      │
        │ 5.0  ┆ 7.0  ┆ 3.666667 │
        │ 9.0  ┆ 9.0  ┆ 6.333333 │
        │ 10.0 ┆ null ┆ 9.0      │
        └──────┴──────┴──────────┘
        """
        return self.select(F.col("*").interpolate())

    def is_empty(self) -> bool:
        """
        Returns `True` if the DataFrame contains no rows.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        >>> df.is_empty()
        False
        >>> df.filter(pl.col("foo") > 99).is_empty()
        True
        """
        return self._df.is_empty()

    def to_struct(self, name: str = "") -> Series:
        """
        Convert a `DataFrame` to a `Series` of type `Struct`.

        Parameters
        ----------
        name
            Name for the struct Series

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5],
        ...         "b": ["one", "two", "three", "four", "five"],
        ...     }
        ... )
        >>> df.to_struct("nums")
        shape: (5,)
        Series: 'nums' [struct[2]]
        [
            {1,"one"}
            {2,"two"}
            {3,"three"}
            {4,"four"}
            {5,"five"}
        ]
        """
        return wrap_s(self._df.to_struct(name, []))

    def unnest(
        self,
        columns: ColumnNameOrSelector | Collection[ColumnNameOrSelector],
        *more_columns: ColumnNameOrSelector,
        separator: str | None = None,
    ) -> DataFrame:
        """
        Decompose struct columns into separate columns for each of their fields.

        The new columns will be inserted into the dataframe at the location of the
        struct column.

        Parameters
        ----------
        columns
            Name of the struct column(s) that should be unnested.
        *more_columns
            Additional columns to unnest, specified as positional arguments.
        separator
            Rename output column names as combination of the struct column name,
            name separator and field name.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "before": ["foo", "bar"],
        ...         "t_a": [1, 2],
        ...         "t_b": ["a", "b"],
        ...         "t_c": [True, None],
        ...         "t_d": [[1, 2], [3]],
        ...         "after": ["baz", "womp"],
        ...     }
        ... ).select("before", pl.struct(pl.col("^t_.$")).alias("t_struct"), "after")
        >>> df
        shape: (2, 3)
        ┌────────┬─────────────────────┬───────┐
        │ before ┆ t_struct            ┆ after │
        │ ---    ┆ ---                 ┆ ---   │
        │ str    ┆ struct[4]           ┆ str   │
        ╞════════╪═════════════════════╪═══════╡
        │ foo    ┆ {1,"a",true,[1, 2]} ┆ baz   │
        │ bar    ┆ {2,"b",null,[3]}    ┆ womp  │
        └────────┴─────────────────────┴───────┘
        >>> df.unnest("t_struct")
        shape: (2, 6)
        ┌────────┬─────┬─────┬──────┬───────────┬───────┐
        │ before ┆ t_a ┆ t_b ┆ t_c  ┆ t_d       ┆ after │
        │ ---    ┆ --- ┆ --- ┆ ---  ┆ ---       ┆ ---   │
        │ str    ┆ i64 ┆ str ┆ bool ┆ list[i64] ┆ str   │
        ╞════════╪═════╪═════╪══════╪═══════════╪═══════╡
        │ foo    ┆ 1   ┆ a   ┆ true ┆ [1, 2]    ┆ baz   │
        │ bar    ┆ 2   ┆ b   ┆ null ┆ [3]       ┆ womp  │
        └────────┴─────┴─────┴──────┴───────────┴───────┘
        >>> df = pl.DataFrame(
        ...     {
        ...         "before": ["foo", "bar"],
        ...         "t_a": [1, 2],
        ...         "t_b": ["a", "b"],
        ...         "t_c": [True, None],
        ...         "t_d": [[1, 2], [3]],
        ...         "after": ["baz", "womp"],
        ...     }
        ... ).select(
        ...     "before",
        ...     pl.struct(pl.col("^t_.$").name.map(lambda t: t[2:])).alias("t"),
        ...     "after",
        ... )
        >>> df.unnest("t", separator="::")
        shape: (2, 6)
        ┌────────┬──────┬──────┬──────┬───────────┬───────┐
        │ before ┆ t::a ┆ t::b ┆ t::c ┆ t::d      ┆ after │
        │ ---    ┆ ---  ┆ ---  ┆ ---  ┆ ---       ┆ ---   │
        │ str    ┆ i64  ┆ str  ┆ bool ┆ list[i64] ┆ str   │
        ╞════════╪══════╪══════╪══════╪═══════════╪═══════╡
        │ foo    ┆ 1    ┆ a    ┆ true ┆ [1, 2]    ┆ baz   │
        │ bar    ┆ 2    ┆ b    ┆ null ┆ [3]       ┆ womp  │
        └────────┴──────┴──────┴──────┴───────────┴───────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .unnest(columns, *more_columns, separator=separator)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def corr(self, **kwargs: Any) -> DataFrame:
        """
        Return pairwise Pearson product-moment correlation coefficients between columns.

        See numpy `corrcoef` for more information:
        https://numpy.org/doc/stable/reference/generated/numpy.corrcoef.html

        Notes
        -----
        This functionality requires numpy to be installed.

        Parameters
        ----------
        **kwargs
            Keyword arguments are passed to numpy `corrcoef`.

        Examples
        --------
        >>> df = pl.DataFrame({"foo": [1, 2, 3], "bar": [3, 2, 1], "ham": [7, 8, 9]})
        >>> df.corr()
        shape: (3, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ f64  ┆ f64  ┆ f64  │
        ╞══════╪══════╪══════╡
        │ 1.0  ┆ -1.0 ┆ 1.0  │
        │ -1.0 ┆ 1.0  ┆ -1.0 │
        │ 1.0  ┆ -1.0 ┆ 1.0  │
        └──────┴──────┴──────┘
        """
        correlation_matrix = np.corrcoef(self.to_numpy(), rowvar=False, **kwargs)
        if self.width == 1:
            correlation_matrix = np.array([correlation_matrix])
        return DataFrame(correlation_matrix, schema=self.columns)

    def merge_sorted(self, other: DataFrame, key: str) -> DataFrame:
        """
        Take two sorted DataFrames and merge them by the sorted key.

        The output of this operation will also be sorted.
        It is the callers responsibility that the frames
        are sorted in ascending order by that key otherwise
        the output will not make sense.

        The schemas of both DataFrames must be equal.

        Parameters
        ----------
        other
            Other DataFrame that must be merged
        key
            Key that is sorted.

        Examples
        --------
        >>> df0 = pl.DataFrame(
        ...     {"name": ["steve", "elise", "bob"], "age": [42, 44, 18]}
        ... ).sort("age")
        >>> df0
        shape: (3, 2)
        ┌───────┬─────┐
        │ name  ┆ age │
        │ ---   ┆ --- │
        │ str   ┆ i64 │
        ╞═══════╪═════╡
        │ bob   ┆ 18  │
        │ steve ┆ 42  │
        │ elise ┆ 44  │
        └───────┴─────┘
        >>> df1 = pl.DataFrame(
        ...     {"name": ["anna", "megan", "steve", "thomas"], "age": [21, 33, 42, 20]}
        ... ).sort("age")
        >>> df1
        shape: (4, 2)
        ┌────────┬─────┐
        │ name   ┆ age │
        │ ---    ┆ --- │
        │ str    ┆ i64 │
        ╞════════╪═════╡
        │ thomas ┆ 20  │
        │ anna   ┆ 21  │
        │ megan  ┆ 33  │
        │ steve  ┆ 42  │
        └────────┴─────┘
        >>> df0.merge_sorted(df1, key="age")
        shape: (7, 2)
        ┌────────┬─────┐
        │ name   ┆ age │
        │ ---    ┆ --- │
        │ str    ┆ i64 │
        ╞════════╪═════╡
        │ bob    ┆ 18  │
        │ thomas ┆ 20  │
        │ anna   ┆ 21  │
        │ megan  ┆ 33  │
        │ steve  ┆ 42  │
        │ steve  ┆ 42  │
        │ elise  ┆ 44  │
        └────────┴─────┘

        Notes
        -----
        No guarantee is given over the output row order when the key is equal
        between the both dataframes.

        The key must be sorted in ascending order.
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        require_same_type(self, other)

        return (
            self.lazy()
            .merge_sorted(other.lazy(), key)
            .collect(optimizations=QueryOptFlags._eager())
        )

    def set_sorted(
        self,
        column: str,
        *,
        descending: bool = False,
    ) -> DataFrame:
        """
        Flag a column as sorted.

        This can speed up future operations.

        Parameters
        ----------
        column
            Column that is sorted
        descending
            Whether the column is sorted in descending order.

        Warnings
        --------
        This can lead to incorrect results if the data is NOT sorted!!
        Use with care!

        """
        # NOTE: Only accepts 1 column on purpose! User think they are sorted by
        # the combined multicolumn values.
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .set_sorted(column, descending=descending)
            .collect(optimizations=QueryOptFlags._eager())
        )

    @unstable()
    def update(
        self,
        other: DataFrame,
        on: str | Sequence[str] | None = None,
        how: Literal["left", "inner", "full"] = "left",
        *,
        left_on: str | Sequence[str] | None = None,
        right_on: str | Sequence[str] | None = None,
        include_nulls: bool = False,
        maintain_order: MaintainOrderJoin | None = "left",
    ) -> DataFrame:
        """
        Update the values in this `DataFrame` with the values in `other`.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        other
            DataFrame that will be used to update the values
        on
            Column names that will be joined on. If set to `None` (default),
            the implicit row index of each frame is used as a join key.
        how : {'left', 'inner', 'full'}
            * 'left' will keep all rows from the left table; rows may be duplicated
              if multiple rows in the right frame match the left row's key.
            * 'inner' keeps only those rows where the key exists in both frames.
            * 'full' will update existing rows where the key matches while also
              adding any new rows contained in the given frame.
        left_on
           Join column(s) of the left DataFrame.
        right_on
           Join column(s) of the right DataFrame.
        include_nulls
            Overwrite values in the left frame with null values from the right frame.
            If set to `False` (default), null values in the right frame are ignored.
        maintain_order : {'none', 'left', 'right', 'left_right', 'right_left'}
            Which order of rows from the inputs to preserve. See :func:`~DataFrame.join`
            for details. Unlike `join` this function preserves the left order by
            default.

        Notes
        -----
        This is syntactic sugar for a left/inner join that preserves the order
        of the left `DataFrame` by default, with an optional coalesce when
        `include_nulls = False`.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "A": [1, 2, 3, 4],
        ...         "B": [400, 500, 600, 700],
        ...     }
        ... )
        >>> df
        shape: (4, 2)
        ┌─────┬─────┐
        │ A   ┆ B   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 400 │
        │ 2   ┆ 500 │
        │ 3   ┆ 600 │
        │ 4   ┆ 700 │
        └─────┴─────┘
        >>> new_df = pl.DataFrame(
        ...     {
        ...         "B": [-66, None, -99],
        ...         "C": [5, 3, 1],
        ...     }
        ... )

        Update `df` values with the non-null values in `new_df`, by row index:

        >>> df.update(new_df)
        shape: (4, 2)
        ┌─────┬─────┐
        │ A   ┆ B   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ -66 │
        │ 2   ┆ 500 │
        │ 3   ┆ -99 │
        │ 4   ┆ 700 │
        └─────┴─────┘

        Update `df` values with the non-null values in `new_df`, by row index,
        but only keeping those rows that are common to both frames:

        >>> df.update(new_df, how="inner")
        shape: (3, 2)
        ┌─────┬─────┐
        │ A   ┆ B   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ -66 │
        │ 2   ┆ 500 │
        │ 3   ┆ -99 │
        └─────┴─────┘

        Update `df` values with the non-null values in `new_df`, using a full
        outer join strategy that defines explicit join columns in each frame:

        >>> df.update(new_df, left_on=["A"], right_on=["C"], how="full")
        shape: (5, 2)
        ┌─────┬─────┐
        │ A   ┆ B   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ -99 │
        │ 2   ┆ 500 │
        │ 3   ┆ 600 │
        │ 4   ┆ 700 │
        │ 5   ┆ -66 │
        └─────┴─────┘

        Update `df` values including null values in `new_df`, using a full outer
        join strategy that defines explicit join columns in each frame:

        >>> df.update(new_df, left_on="A", right_on="C", how="full", include_nulls=True)
        shape: (5, 2)
        ┌─────┬──────┐
        │ A   ┆ B    │
        │ --- ┆ ---  │
        │ i64 ┆ i64  │
        ╞═════╪══════╡
        │ 1   ┆ -99  │
        │ 2   ┆ 500  │
        │ 3   ┆ null │
        │ 4   ┆ 700  │
        │ 5   ┆ -66  │
        └─────┴──────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        require_same_type(self, other)
        return (
            self.lazy()
            .update(
                other.lazy(),
                on,
                how,
                left_on=left_on,
                right_on=right_on,
                include_nulls=include_nulls,
                maintain_order=maintain_order,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    def count(self) -> DataFrame:
        """
        Return the number of non-null elements for each column.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {"a": [1, 2, 3, 4], "b": [1, 2, 1, None], "c": [None, None, None, None]}
        ... )
        >>> df.count()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ u32 ┆ u32 ┆ u32 │
        ╞═════╪═════╪═════╡
        │ 4   ┆ 3   ┆ 0   │
        └─────┴─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return self.lazy().count().collect(optimizations=QueryOptFlags._eager())

    @deprecated(
        "`DataFrame.melt` is deprecated; use `DataFrame.unpivot` instead, with "
        "`index` instead of `id_vars` and `on` instead of `value_vars`"
    )
    def melt(
        self,
        id_vars: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        value_vars: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        variable_name: str | None = None,
        value_name: str | None = None,
    ) -> DataFrame:
        """
        Unpivot a DataFrame from wide to long format.

        Optionally leaves identifiers set.

        This function is useful to massage a DataFrame into a format where one or more
        columns are identifier variables (id_vars) while all other columns, considered
        measured variables (value_vars), are "unpivoted" to the row axis leaving just
        two non-identifier columns, 'variable' and 'value'.

        .. deprecated:: 1.0.0
            Use the :meth:`.unpivot` method instead.

        Parameters
        ----------
        id_vars
            Column(s) or selector(s) to use as identifier variables.
        value_vars
            Column(s) or selector(s) to use as values variables; if `value_vars`
            is empty all columns that are not in `id_vars` will be used.
        variable_name
            Name to give to the `variable` column. Defaults to "variable"
        value_name
            Name to give to the `value` column. Defaults to "value"
        """
        return self.unpivot(
            index=id_vars,
            on=value_vars,
            variable_name=variable_name,
            value_name=value_name,
        )

    @unstable()
    def match_to_schema(
        self,
        schema: SchemaDict | Schema,
        *,
        missing_columns: Literal["insert", "raise"]
        | Mapping[str, Literal["insert", "raise"] | Expr] = "raise",
        missing_struct_fields: Literal["insert", "raise"]
        | Mapping[str, Literal["insert", "raise"]] = "raise",
        extra_columns: Literal["ignore", "raise"] = "raise",
        extra_struct_fields: Literal["ignore", "raise"]
        | Mapping[str, Literal["ignore", "raise"]] = "raise",
        integer_cast: Literal["upcast", "forbid"]
        | Mapping[str, Literal["upcast", "forbid"]] = "forbid",
        float_cast: Literal["upcast", "forbid"]
        | Mapping[str, Literal["upcast", "forbid"]] = "forbid",
    ) -> DataFrame:
        """
        Match or evolve the schema of a LazyFrame into a specific schema.

        By default, match_to_schema returns an error if the input schema does not
        exactly match the target schema. It also allows columns to be freely reordered,
        with additional coercion rules available through optional parameters.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        schema
            Target schema to match or evolve to.
        missing_columns
            Raise of insert missing columns from the input with respect to the `schema`.

            This can also be an expression per column with what to insert if it is
            missing.
        missing_struct_fields
            Raise of insert missing struct fields from the input with respect to the
            `schema`.
        extra_columns
            Raise of ignore extra columns from the input with respect to the `schema`.
        extra_struct_fields
            Raise of ignore extra struct fields from the input with respect to the
            `schema`.
        integer_cast
            Forbid of upcast for integer columns from the input to the respective column
            in `schema`.
        float_cast
            Forbid of upcast for float columns from the input to the respective column
            in `schema`.

        Examples
        --------
        Ensuring the schema matches

        >>> df = pl.DataFrame({"a": [1, 2, 3], "b": ["A", "B", "C"]})
        >>> df.match_to_schema({"a": pl.Int64, "b": pl.String})
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ A   │
        │ 2   ┆ B   │
        │ 3   ┆ C   │
        └─────┴─────┘
        >>> df.match_to_schema({"a": pl.Int64})  # doctest: +SKIP
        polars.exceptions.SchemaError: extra columns in `match_to_schema`: "b"

        Adding missing columns

        >>> (
        ...     pl.DataFrame({"a": [1, 2, 3]}).match_to_schema(
        ...         {"a": pl.Int64, "b": pl.String},
        ...         missing_columns="insert",
        ...     )
        ... )
        shape: (3, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ i64 ┆ str  │
        ╞═════╪══════╡
        │ 1   ┆ null │
        │ 2   ┆ null │
        │ 3   ┆ null │
        └─────┴──────┘
        >>> (
        ...     pl.DataFrame({"a": [1, 2, 3]}).match_to_schema(
        ...         {"a": pl.Int64, "b": pl.String},
        ...         missing_columns={"b": pl.col.a.cast(pl.String)},
        ...     )
        ... )
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ 1   │
        │ 2   ┆ 2   │
        │ 3   ┆ 3   │
        └─────┴─────┘

        Removing extra columns

        >>> (
        ...     pl.DataFrame({"a": [1, 2, 3], "b": ["A", "B", "C"]}).match_to_schema(
        ...         {"a": pl.Int64},
        ...         extra_columns="ignore",
        ...     )
        ... )
        shape: (3, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 1   │
        │ 2   │
        │ 3   │
        └─────┘

        Upcasting integers and floats

        >>> (
        ...     pl.DataFrame(
        ...         {"a": [1, 2, 3], "b": [1.0, 2.0, 3.0]},
        ...         schema={"a": pl.Int32, "b": pl.Float32},
        ...     ).match_to_schema(
        ...         {"a": pl.Int64, "b": pl.Float64},
        ...         integer_cast="upcast",
        ...         float_cast="upcast",
        ...     )
        ... )
        shape: (3, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ f64 │
        ╞═════╪═════╡
        │ 1   ┆ 1.0 │
        │ 2   ┆ 2.0 │
        │ 3   ┆ 3.0 │
        └─────┴─────┘
        """
        from polars.lazyframe.opt_flags import QueryOptFlags

        return (
            self.lazy()
            .match_to_schema(
                schema=schema,
                missing_columns=missing_columns,
                missing_struct_fields=missing_struct_fields,
                extra_columns=extra_columns,
                extra_struct_fields=extra_struct_fields,
                integer_cast=integer_cast,
                float_cast=float_cast,
            )
            .collect(optimizations=QueryOptFlags._eager())
        )

    def _to_metadata(
        self,
        columns: None | str | list[str] = None,
        stats: None | str | list[str] = None,
    ) -> DataFrame:
        """
        Get all runtime metadata for each column.

        This is unstable and is meant for debugging purposes.

        Parameters
        ----------
        columns
            Column(s) to show the information for
        stats
            Statistics to show
        """
        df = self

        if columns is not None:
            if isinstance(columns, str):
                columns = [columns]

            df = df.select(columns)

        md = self._from_pydf(df._df._to_metadata())

        if stats is not None:
            if isinstance(stats, str):
                stats = [stats]

            if "column_name" not in stats:
                stats = ["column_name"] + stats

            md = md.select(stats)

        return md

    def _row_encode(
        self,
        *,
        unordered: bool = False,
        descending: list[bool] | None = None,
        nulls_last: list[bool] | None = None,
    ) -> Series:
        """
        Row encode the given DataFrame.

        This is an internal function not meant for outside consumption and can
        be changed or removed at any point in time.

        fields have order:
        - descending
        - nulls_last
        - no_order
        """
        return self.select_seq(
            F._row_encode(
                F.all(),
                unordered=unordered,
                descending=descending,
                nulls_last=nulls_last,
            )
        ).to_series()


def _prepare_other_arg(other: Any, length: int | None = None) -> Series:
    # if not a series create singleton series such that it will broadcast
    value = other
    if not isinstance(other, pl.Series):
        if isinstance(other, str):
            pass
        elif isinstance(other, Sequence):
            msg = "operation not supported"
            raise TypeError(msg)
        other = pl.Series("", [other])

    if length is not None:
        if length > 1:
            other = other.extend_constant(value=value, n=length - 1)
        elif length == 0:
            other = other.slice(0, 0)

    return other
