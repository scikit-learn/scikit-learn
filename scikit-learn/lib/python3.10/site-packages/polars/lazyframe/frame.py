from __future__ import annotations

import contextlib
import io
import os
import warnings
from collections.abc import Collection, Iterable, Iterator, Mapping
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, time, timedelta
from functools import lru_cache, partial, reduce
from io import BytesIO, StringIO
from operator import and_
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    NoReturn,
    TypeVar,
    overload,
)

import polars._reexport as pl
import polars.selectors as cs
from polars import functions as F
from polars._dependencies import (
    _PYARROW_AVAILABLE,
    import_optional,
    subprocess,
)
from polars._dependencies import polars_cloud as pc
from polars._dependencies import pyarrow as pa
from polars._typing import (
    ParquetMetadata,
    PartitioningScheme,
)
from polars._utils.async_ import _AioDataFrameResult, _GeventDataFrameResult
from polars._utils.convert import negate_duration_string, parse_as_duration_string
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    deprecate_streaming_parameter,
    deprecated,
    issue_deprecation_warning,
)
from polars._utils.parquet import wrap_parquet_metadata_callback
from polars._utils.parse import (
    parse_into_expression,
    parse_into_list_of_expressions,
)
from polars._utils.parse.expr import parse_list_into_selector
from polars._utils.serde import serialize_polars_object
from polars._utils.slice import LazyPolarsSlice
from polars._utils.unstable import issue_unstable_warning, unstable
from polars._utils.various import (
    _is_generator,
    display_dot_graph,
    extend_bool,
    find_stacklevel,
    is_bool_sequence,
    is_sequence,
    issue_warning,
    normalize_filepath,
    parse_percentiles,
    qualified_type_name,
    require_same_type,
)
from polars._utils.wrap import wrap_df, wrap_expr
from polars.datatypes import (
    DTYPE_TEMPORAL_UNITS,
    N_INFER_DEFAULT,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Duration,
    Enum,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    Null,
    Object,
    String,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Unknown,
    is_polars_dtype,
    parse_into_datatype_expr,
    parse_into_dtype,
)
from polars.datatypes.group import DataTypeGroup
from polars.exceptions import PerformanceWarning
from polars.interchange.protocol import CompatLevel
from polars.lazyframe.engine_config import GPUEngine
from polars.lazyframe.group_by import LazyGroupBy
from polars.lazyframe.in_process import InProcessQuery
from polars.lazyframe.opt_flags import DEFAULT_QUERY_OPT_FLAGS, forward_old_opt_flags
from polars.schema import Schema
from polars.selectors import by_dtype, expand_selector

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyLazyFrame, get_engine_affinity

if TYPE_CHECKING:
    import sys
    from collections.abc import Awaitable, Iterator, Sequence
    from concurrent.futures import Future
    from io import IOBase
    from typing import IO, Literal

    from polars.lazyframe.opt_flags import QueryOptFlags

    with contextlib.suppress(ImportError):  # Module not available when building docs
        from polars._plr import PyExpr, PyPartitioning, PySelector

    with contextlib.suppress(ImportError):  # Module not available when building docs
        import polars._plr as plr

    from polars import DataFrame, DataType, Expr
    from polars._dependencies import numpy as np
    from polars._typing import (
        AsofJoinStrategy,
        ClosedInterval,
        ColumnNameOrSelector,
        CsvQuoteStyle,
        EngineType,
        ExplainFormat,
        FillNullStrategy,
        FrameInitTypes,
        IntoExpr,
        IntoExprColumn,
        IpcCompression,
        JoinStrategy,
        JoinValidation,
        Label,
        MaintainOrderJoin,
        Orientation,
        ParquetMetadata,
        PlanStage,
        PolarsDataType,
        PythonDataType,
        QuantileMethod,
        SchemaDefinition,
        SchemaDict,
        SerializationFormat,
        StartBy,
        SyncOnCloseMethod,
        UniqueKeepStrategy,
    )
    from polars.io.cloud import CredentialProviderFunction
    from polars.io.parquet import ParquetFieldOverwrites

    if sys.version_info >= (3, 10):
        from typing import Concatenate, ParamSpec
    else:
        from typing_extensions import Concatenate, ParamSpec

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004

    T = TypeVar("T")
    P = ParamSpec("P")


_COLLECT_BATCHES_POOL = ThreadPoolExecutor(thread_name_prefix="pl_col_batch_")


def _select_engine(engine: EngineType) -> EngineType:
    return get_engine_affinity() if engine == "auto" else engine


def _to_sink_target(
    path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
) -> str | Path | IO[bytes] | IO[str] | PyPartitioning:
    if isinstance(path, (str, Path)):
        return normalize_filepath(path)
    elif isinstance(path, io.IOBase):
        return path  # type: ignore[return-value]
    elif isinstance(path, PartitioningScheme):
        return path._py_partitioning
    elif callable(getattr(path, "write", None)):
        # This allows for custom writers
        return path
    else:
        msg = f"`path` argument has invalid type {qualified_type_name(path)!r}, and cannot be turned into a sink target"
        raise TypeError(msg)


def _gpu_engine_callback(
    engine: EngineType,
    *,
    streaming: bool,
    background: bool,
    new_streaming: bool,
    _eager: bool,
) -> Callable[[Any, int | None], None] | None:
    is_gpu = (is_config_obj := isinstance(engine, GPUEngine)) or engine == "gpu"
    if not (
        is_config_obj or engine in ("auto", "cpu", "in-memory", "streaming", "gpu")
    ):
        msg = f"Invalid engine argument {engine=}"
        raise ValueError(msg)
    if (streaming or background or new_streaming) and is_gpu:
        issue_warning(
            "GPU engine does not support streaming or background collection, "
            "disabling GPU engine.",
            category=UserWarning,
        )
        is_gpu = False
    if _eager:
        # Don't run on GPU in _eager mode (but don't warn)
        is_gpu = False

    if not is_gpu:
        return None
    cudf_polars = import_optional(
        "cudf_polars",
        err_prefix="GPU engine requested, but required package",
        install_message=(
            "Please install using the command "
            "`pip install cudf-polars-cu12` "
            "(CUDA 12 is required for RAPIDS cuDF v25.08 and later). "
            "If your system has a CUDA 11 driver, install with "
            "`pip install cudf-polars-cu11==25.06` "
        ),
    )
    if not is_config_obj:
        engine = GPUEngine()
    return partial(cudf_polars.execute_with_cudf, config=engine)


class LazyFrame:
    """
    Representation of a Lazy computation graph/query against a DataFrame.

    This allows for whole-query optimisation in addition to parallelism, and
    is the preferred (and highest-performance) mode of operation for polars.

    Parameters
    ----------
    data : dict, Sequence, ndarray, Series, or pandas.DataFrame
        Two-dimensional data in various forms; dict input must contain Sequences,
        Generators, or a `range`. Sequence may contain Series or other Sequences.
    schema : Sequence of str, (str,DataType) pairs, or a {str:DataType,} dict
        The LazyFrame schema may be declared in several ways:

        * As a dict of {name:type} pairs; if type is None, it will be auto-inferred.
        * As a list of column names; in this case types are automatically inferred.
        * As a list of (name,type) pairs; this is equivalent to the dictionary form.

        If you supply a list of column names that does not match the names in the
        underlying data, the names given here will overwrite them. The number
        of names given in the schema should match the underlying data dimensions.
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
    Initialising `LazyFrame(...)` directly is equivalent to `DataFrame(...).lazy()`.

    Examples
    --------
    Constructing a LazyFrame directly from a dictionary:

    >>> data = {"a": [1, 2], "b": [3, 4]}
    >>> lf = pl.LazyFrame(data)
    >>> lf.collect()
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    │ 2   ┆ 4   │
    └─────┴─────┘

    Notice that the dtypes are automatically inferred as Polars Int64:

    >>> lf.collect_schema().dtypes()
    [Int64, Int64]

    To specify a more detailed/specific frame schema you can supply the `schema`
    parameter with a dictionary of (name,dtype) pairs...

    >>> data = {"col1": [0, 2], "col2": [3, 7]}
    >>> lf2 = pl.LazyFrame(data, schema={"col1": pl.Float32, "col2": pl.Int64})
    >>> lf2.collect()
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
    >>> lf3 = pl.LazyFrame(data, schema=[("col1", pl.Float32), ("col2", pl.Int64)])
    >>> lf3.collect()
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
    >>> lf4 = pl.LazyFrame(data)
    >>> lf4.collect()
    shape: (2, 2)
    ┌──────┬──────┐
    │ col1 ┆ col2 │
    │ ---  ┆ ---  │
    │ f32  ┆ i64  │
    ╞══════╪══════╡
    │ 1.0  ┆ 3    │
    │ 2.0  ┆ 4    │
    └──────┴──────┘

    Constructing a LazyFrame from a numpy ndarray, specifying column names:

    >>> import numpy as np
    >>> data = np.array([(1, 2), (3, 4)], dtype=np.int64)
    >>> lf5 = pl.LazyFrame(data, schema=["a", "b"], orient="col")
    >>> lf5.collect()
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    │ 2   ┆ 4   │
    └─────┴─────┘

    Constructing a LazyFrame from a list of lists, row orientation specified:

    >>> data = [[1, 2, 3], [4, 5, 6]]
    >>> lf6 = pl.LazyFrame(data, schema=["a", "b", "c"], orient="row")
    >>> lf6.collect()
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

    _ldf: PyLazyFrame
    _accessors: ClassVar[set[str]] = set()

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
        from polars.dataframe import DataFrame

        self._ldf = (
            DataFrame(
                data=data,
                schema=schema,
                schema_overrides=schema_overrides,
                strict=strict,
                orient=orient,
                infer_schema_length=infer_schema_length,
                nan_to_null=nan_to_null,
            )
            .lazy()
            ._ldf
        )

    @classmethod
    def _from_pyldf(cls, ldf: PyLazyFrame) -> LazyFrame:
        self = cls.__new__(cls)
        self._ldf = ldf
        return self

    def __getstate__(self) -> bytes:
        return self.serialize()

    def __setstate__(self, state: bytes) -> None:
        self._ldf = self.deserialize(BytesIO(state))._ldf

    @classmethod
    def _scan_python_function(
        cls,
        schema: pa.schema | SchemaDict | Callable[[], SchemaDict],
        scan_fn: Any,
        *,
        pyarrow: bool = False,
        validate_schema: bool = False,
        is_pure: bool = False,
    ) -> LazyFrame:
        self = cls.__new__(cls)
        if isinstance(schema, Mapping):
            self._ldf = PyLazyFrame.scan_from_python_function_pl_schema(
                list(schema.items()),
                scan_fn,
                pyarrow=pyarrow,
                validate_schema=validate_schema,
                is_pure=is_pure,
            )
        elif _PYARROW_AVAILABLE and isinstance(schema, pa.Schema):
            self._ldf = PyLazyFrame.scan_from_python_function_arrow_schema(
                list(schema),
                scan_fn,
                pyarrow=pyarrow,
                validate_schema=validate_schema,
                is_pure=is_pure,
            )
        else:
            self._ldf = PyLazyFrame.scan_from_python_function_schema_function(
                schema, scan_fn, validate_schema=validate_schema, is_pure=is_pure
            )
        return self

    @classmethod
    def deserialize(
        cls, source: str | Path | IOBase, *, format: SerializationFormat = "binary"
    ) -> LazyFrame:
        """
        Read a logical plan from a file to construct a LazyFrame.

        Parameters
        ----------
        source
            Path to a file or a file-like object (by file-like object, we refer to
            objects that have a `read()` method, such as a file handler (e.g.
            via builtin `open` function) or `BytesIO`).
        format
            The format with which the LazyFrame was serialized. Options:

            - `"binary"`: Deserialize from binary format (bytes). This is the default.
            - `"json"`: Deserialize from JSON format (string).

        Warnings
        --------
        This function uses :mod:`pickle` if the logical plan contains Python UDFs,
        and as such inherits the security implications. Deserializing can execute
        arbitrary code, so it should only be attempted on trusted data.

        See Also
        --------
        LazyFrame.serialize

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        >>> import io
        >>> lf = pl.LazyFrame({"a": [1, 2, 3]}).sum()
        >>> bytes = lf.serialize()
        >>> pl.LazyFrame.deserialize(io.BytesIO(bytes)).collect()
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 6   │
        └─────┘
        """
        if isinstance(source, StringIO):
            source = BytesIO(source.getvalue().encode())
        elif isinstance(source, (str, Path)):
            source = normalize_filepath(source)

        if format == "binary":
            deserializer = PyLazyFrame.deserialize_binary
        elif format == "json":
            deserializer = PyLazyFrame.deserialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return cls._from_pyldf(deserializer(source))

    @property
    def columns(self) -> list[str]:
        """
        Get the column names.

        Returns
        -------
        list of str
            A list containing the name of each column in order.

        Warnings
        --------
        Determining the column names of a LazyFrame requires resolving its schema,
        which is a potentially expensive operation.
        Using :meth:`collect_schema` is the idiomatic way of resolving the schema.
        This property exists only for symmetry with the DataFrame class.

        See Also
        --------
        collect_schema
        Schema.names

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... ).select("foo", "bar")
        >>> lf.columns  # doctest: +SKIP
        ['foo', 'bar']
        """
        issue_warning(
            "Determining the column names of a LazyFrame requires resolving its schema,"
            " which is a potentially expensive operation. Use `LazyFrame.collect_schema().names()`"
            " to get the column names without this warning.",
            category=PerformanceWarning,
        )
        return self.collect_schema().names()

    @property
    def dtypes(self) -> list[DataType]:
        """
        Get the column data types.

        Returns
        -------
        list of DataType
            A list containing the data type of each column in order.

        Warnings
        --------
        Determining the data types of a LazyFrame requires resolving its schema,
        which is a potentially expensive operation.
        Using :meth:`collect_schema` is the idiomatic way to resolve the schema.
        This property exists only for symmetry with the DataFrame class.

        See Also
        --------
        collect_schema
        Schema.dtypes

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.dtypes  # doctest: +SKIP
        [Int64, Float64, String]
        """
        issue_warning(
            "Determining the data types of a LazyFrame requires resolving its schema,"
            " which is a potentially expensive operation. Use `LazyFrame.collect_schema().dtypes()`"
            " to get the data types without this warning.",
            category=PerformanceWarning,
        )
        return self.collect_schema().dtypes()

    @property
    def schema(self) -> Schema:
        """
        Get an ordered mapping of column names to their data type.

        Warnings
        --------
        Resolving the schema of a LazyFrame is a potentially expensive operation.
        Using :meth:`collect_schema` is the idiomatic way to resolve the schema.
        This property exists only for symmetry with the DataFrame class.

        See Also
        --------
        collect_schema
        Schema

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.schema  # doctest: +SKIP
        Schema({'foo': Int64, 'bar': Float64, 'ham': String})
        """
        issue_warning(
            "Resolving the schema of a LazyFrame is a potentially expensive operation."
            " Use `LazyFrame.collect_schema()` to get the schema without this warning.",
            category=PerformanceWarning,
        )
        return self.collect_schema()

    @property
    def width(self) -> int:
        """
        Get the number of columns.

        Returns
        -------
        int

        Warnings
        --------
        Determining the width of a LazyFrame requires resolving its schema,
        which is a potentially expensive operation.
        Using :meth:`collect_schema` is the idiomatic way to resolve the schema.
        This property exists only for symmetry with the DataFrame class.

        See Also
        --------
        collect_schema
        Schema.len

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [4, 5, 6],
        ...     }
        ... )
        >>> lf.width  # doctest: +SKIP
        2
        """
        issue_warning(
            "determining the width of a LazyFrame requires resolving its schema,"
            " which is a potentially expensive operation. Use `LazyFrame.collect_schema().len()`"
            " to get the width without this warning.",
            category=PerformanceWarning,
        )
        return self.collect_schema().len()

    def __bool__(self) -> NoReturn:
        msg = (
            "the truth value of a LazyFrame is ambiguous"
            "\n\nLazyFrames cannot be used in boolean context with and/or/not operators."
        )
        raise TypeError(msg)

    def _comparison_error(self, operator: str) -> NoReturn:
        msg = f'"{operator!r}" comparison not supported for LazyFrame objects'
        raise TypeError(msg)

    def __eq__(self, other: object) -> NoReturn:
        self._comparison_error("==")

    def __ne__(self, other: object) -> NoReturn:
        self._comparison_error("!=")

    def __gt__(self, other: Any) -> NoReturn:
        self._comparison_error(">")

    def __lt__(self, other: Any) -> NoReturn:
        self._comparison_error("<")

    def __ge__(self, other: Any) -> NoReturn:
        self._comparison_error(">=")

    def __le__(self, other: Any) -> NoReturn:
        self._comparison_error("<=")

    def __contains__(self, key: str) -> bool:
        return key in self.collect_schema()

    def __copy__(self) -> LazyFrame:
        return self.clone()

    def __deepcopy__(self, memo: None = None) -> LazyFrame:
        return self.clone()

    def __getitem__(self, item: slice) -> LazyFrame:
        """
        Support slice syntax, returning a new LazyFrame.

        All other forms of subscripting are currently unsupported here; use `select`,
        `filter`, or other standard methods instead.

        Notes
        -----
        LazyFrame is designed primarily for efficient computation and does not know
        its own length so, unlike DataFrame, certain slice patterns (such as those
        requiring negative stop/step) may not be supported.

        Examples
        --------
        >>> lf = pl.LazyFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        >>> lf[:2].collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 4   │
        │ 2   ┆ 5   │
        └─────┴─────┘
        >>> lf[::2].collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 4   │
        │ 3   ┆ 6   │
        └─────┴─────┘
        """
        if not isinstance(item, slice):
            msg = (
                "LazyFrame is not subscriptable (aside from slicing)"
                "\n\nUse `select()` or `filter()` instead."
            )
            raise TypeError(msg)
        return LazyPolarsSlice(self).apply(item)

    def __str__(self) -> str:
        return f"""\
naive plan: (run LazyFrame.explain(optimized=True) to see the optimized plan)

{self.explain(optimized=False)}\
"""

    def __repr__(self) -> str:
        # don't expose internal/private classpath
        return f"<{self.__class__.__name__} at 0x{id(self):X}>"

    def _repr_html_(self) -> str:
        try:
            dot = self._ldf.to_dot(optimized=False)
            svg = subprocess.check_output(
                ["dot", "-Nshape=box", "-Tsvg"], input=f"{dot}".encode()
            )
            return (
                "<h4>NAIVE QUERY PLAN</h4><p>run <b>LazyFrame.show_graph()</b> to see"
                f" the optimized version</p>{svg.decode()}"
            )
        except Exception:
            insert = self.explain(optimized=False).replace("\n", "<p></p>")

            return f"""\
<i>naive plan: (run <b>LazyFrame.explain(optimized=True)</b> to see the optimized plan)</i>
    <p></p>
    <div>{insert}</div>\
"""

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
        Serialize the logical plan of this LazyFrame to a file or string in JSON format.

        Parameters
        ----------
        file
            File path to which the result should be written. If set to `None`
            (default), the output is returned as a string instead.
        format
            The format in which to serialize. Options:

            - `"binary"`: Serialize to binary format (bytes). This is the default.
            - `"json"`: Serialize to JSON format (string) (deprecated).

        See Also
        --------
        LazyFrame.deserialize

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        Serialize the logical plan into a binary representation.

        >>> lf = pl.LazyFrame({"a": [1, 2, 3]}).sum()
        >>> bytes = lf.serialize()

        The bytes can later be deserialized back into a LazyFrame.

        >>> import io
        >>> pl.LazyFrame.deserialize(io.BytesIO(bytes)).collect()
        shape: (1, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 6   │
        └─────┘
        """
        if format == "binary":
            serializer = self._ldf.serialize_binary
        elif format == "json":
            msg = "'json' serialization format of LazyFrame is deprecated"
            warnings.warn(
                msg,
                stacklevel=find_stacklevel(),
            )
            serializer = self._ldf.serialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return serialize_polars_object(serializer, file, format)

    def pipe(
        self,
        function: Callable[Concatenate[LazyFrame, P], T],
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

        See Also
        --------
        pipe_with_schema

        Examples
        --------
        >>> def cast_str_to_int(lf: pl.LazyFrame, col_name: str) -> pl.LazyFrame:
        ...     return lf.with_columns(pl.col(col_name).cast(pl.Int64))
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": ["10", "20", "30", "40"],
        ...     }
        ... )
        >>> lf.pipe(cast_str_to_int, col_name="b").collect()
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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "b": [1, 2],
        ...         "a": [3, 4],
        ...     }
        ... )
        >>> lf.collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ b   ┆ a   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 3   │
        │ 2   ┆ 4   │
        └─────┴─────┘
        >>> lf.pipe(lambda lf: lf.select(sorted(lf.collect_schema()))).collect()
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

    @unstable()
    def pipe_with_schema(
        self,
        function: Callable[[LazyFrame, Schema], LazyFrame],
    ) -> LazyFrame:
        """
        Allows to alter the lazy frame during the plan stage with the resolved schema.

        In contrast to `pipe`, this method does not execute `function` immediately but
        only during the plan stage. This allows to use the resolved schema of the input
        to dynamically alter the lazy frame. This also means that any exceptions raised
        by `function` will only be emitted during the plan stage.

        .. warning::
            This functionality is considered **unstable**. It may be changed at any
            point without it being considered a breaking change.

        Parameters
        ----------
        function
            Callable; will receive the frame as the first parameter and the resolved
            schema as the second parameter.

        See Also
        --------
        pipe

        Examples
        --------
        >>> def cast_to_float_if_necessary(
        ...     lf: pl.LazyFrame, schema: pl.Schema
        ... ) -> pl.LazyFrame:
        ...     required_casts = [
        ...         pl.col(name).cast(pl.Float64)
        ...         for name, dtype in schema.items()
        ...         if not dtype.is_float()
        ...     ]
        ...     return lf.with_columns(required_casts)
        >>> lf = pl.LazyFrame(
        ...     {"a": [1.0, 2.0], "b": ["1.0", "2.5"], "c": [2.0, 3.0]},
        ...     schema={"a": pl.Float64, "b": pl.String, "c": pl.Float32},
        ... )
        >>> lf.pipe_with_schema(cast_to_float_if_necessary).collect()
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ f64 ┆ f64 ┆ f32 │
        ╞═════╪═════╪═════╡
        │ 1.0 ┆ 1.0 ┆ 2.0 │
        │ 2.0 ┆ 2.5 ┆ 3.0 │
        └─────┴─────┴─────┘
        """
        return self._from_pyldf(
            self._ldf.pipe_with_schema(
                lambda lf_and_schema: function(
                    self._from_pyldf(lf_and_schema[0]),
                    lf_and_schema[1],
                )._ldf
            )
        )

    def describe(
        self,
        percentiles: Sequence[float] | float | None = (0.25, 0.50, 0.75),
        *,
        interpolation: QuantileMethod = "nearest",
    ) -> DataFrame:
        """
        Creates a summary of statistics for a LazyFrame, returning a DataFrame.

        Parameters
        ----------
        percentiles
            One or more percentiles to include in the summary statistics.
            All values must be in the range `[0, 1]`.

        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method used when calculating percentiles.

        Returns
        -------
        DataFrame

        Notes
        -----
        The median is included by default as the 50% percentile.

        Warnings
        --------
        * This method does *not* maintain the laziness of the frame, and will `collect`
          the final result. This could potentially be an expensive operation.
        * We do not guarantee the output of `describe` to be stable. It will show
          statistics that we deem informative, and may be updated in the future.
          Using `describe` programmatically (versus interactive exploration) is
          not recommended for this reason.

        Examples
        --------
        >>> from datetime import date, time
        >>> lf = pl.LazyFrame(
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

        >>> lf.describe()
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
        ...     lf.describe(
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
        from polars.convert import from_dict

        schema = self.collect_schema()

        if not schema:
            msg = "cannot describe a LazyFrame that has no columns"
            raise TypeError(msg)

        # create list of metrics
        metrics = ["count", "null_count", "mean", "std", "min"]
        if quantiles := parse_percentiles(percentiles):
            metrics.extend(f"{q * 100:g}%" for q in quantiles)
        metrics.append("max")

        @lru_cache
        def skip_minmax(dt: PolarsDataType) -> bool:
            return dt.is_nested() or dt in (Categorical, Enum, Null, Object, Unknown)

        # determine which columns will produce std/mean/percentile/etc
        # statistics in a single pass over the frame schema
        has_numeric_result, sort_cols = set(), set()
        metric_exprs: list[Expr] = []
        null = F.lit(None)

        for c, dtype in schema.items():
            is_numeric = dtype.is_numeric()
            is_temporal = not is_numeric and dtype.is_temporal()

            # counts
            count_exprs = [
                F.col(c).count().name.prefix("count:"),
                F.col(c).null_count().name.prefix("null_count:"),
            ]
            # mean
            mean_expr = (
                F.col(c).mean()
                if is_temporal or is_numeric or dtype == Boolean
                else null
            )

            # standard deviation, min, max
            expr_std = F.col(c).std() if is_numeric else null
            min_expr = F.col(c).min() if not skip_minmax(dtype) else null
            max_expr = F.col(c).max() if not skip_minmax(dtype) else null

            # percentiles
            pct_exprs = []
            for p in quantiles:
                if is_numeric or is_temporal:
                    pct_expr = (
                        F.col(c).to_physical().quantile(p, interpolation).cast(dtype)
                        if is_temporal
                        else F.col(c).quantile(p, interpolation)
                    )
                    sort_cols.add(c)
                else:
                    pct_expr = null
                pct_exprs.append(pct_expr.alias(f"{p}:{c}"))

            if is_numeric or dtype.is_nested() or dtype in (Null, Boolean):
                has_numeric_result.add(c)

            # add column expressions (in end-state 'metrics' list order)
            metric_exprs.extend(
                [
                    *count_exprs,
                    mean_expr.alias(f"mean:{c}"),
                    expr_std.alias(f"std:{c}"),
                    min_expr.alias(f"min:{c}"),
                    *pct_exprs,
                    max_expr.alias(f"max:{c}"),
                ]
            )

        # calculate requested metrics in parallel, then collect the result
        df_metrics = (
            (
                # if more than one quantile, sort the relevant columns to make them O(1)
                # TODO: drop sort once we have efficient retrieval of multiple quantiles
                self.with_columns(F.col(c).sort() for c in sort_cols)
                if sort_cols
                else self
            )
            .select(*metric_exprs)
            .collect()
        )

        # reshape wide result
        n_metrics = len(metrics)
        column_metrics = [
            df_metrics.row(0)[(n * n_metrics) : (n + 1) * n_metrics]
            for n in range(schema.len())
        ]
        summary = dict(zip(schema, column_metrics))

        # cast by column type (numeric/bool -> float), (other -> string)
        for c in schema:
            summary[c] = [  # type: ignore[assignment]
                (
                    None
                    if (v is None or isinstance(v, dict))
                    else (float(v) if (c in has_numeric_result) else str(v))
                )
                for v in summary[c]
            ]

        # return results as a DataFrame
        df_summary = from_dict(summary)
        df_summary.insert_column(0, pl.Series("statistic", metrics))
        return df_summary

    @deprecate_streaming_parameter()
    @forward_old_opt_flags()
    def explain(
        self,
        *,
        format: ExplainFormat = "plain",
        optimized: bool = True,
        type_coercion: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        streaming: bool = False,
        engine: EngineType = "auto",
        tree_format: bool | None = None,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> str:
        """
        Create a string representation of the query plan.

        Different optimizations can be turned on or off.

        Parameters
        ----------
        format : {'plain', 'tree'}
            The format to use for displaying the logical plan.
        optimized
            Return an optimized query plan. Defaults to `True`.
            If this is set to `True` the subsequent
            optimization flags control which optimizations
            run.
        type_coercion
            Do type coercion optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        predicate_pushdown
            Do predicate pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        projection_pushdown
            Do projection pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        simplify_expression
            Run simplify expressions optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        slice_pushdown
            Slice pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subplan_elim
            Will try to cache branching subplans that occur on self-joins or unions.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subexpr_elim
            Common subexpressions will be cached and reused.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        cluster_with_columns
            Combine sequential independent calls to with_columns

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        collapse_joins
            Collapse a join and filters into a faster join

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        streaming
            Unused parameter, kept for backward compatibility.

           ... deprecated:: 1.30.0
                Use the `engine` parameter instead.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query
            is run using the polars in-memory engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars in-memory
            engine. If set to `"gpu"`, the GPU engine is used. Fine-grained
            control over the GPU engine, for example which device to use
            on a system with multiple devices, is possible by providing a
            :class:`~.GPUEngine` object with configuration options.

            .. note::
               GPU mode is considered **unstable**. Not all queries will run
               successfully on the GPU, however, they should fall back transparently
               to the default engine if execution is not supported.

               Running with `POLARS_VERBOSE=1` will provide information if a query
               falls back (and why).

            .. note::
               The GPU engine does not support streaming, if streaming
               is enabled then GPU execution is switched off.
        optimizations
            The optimization passes done during query optimization.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        tree_format
            Format the output as a tree.

            .. deprecated:: 0.20.30
                Use `format="tree"` instead.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [1, 2, 3, 4, 5, 6],
        ...         "c": [6, 5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> lf.group_by("a", maintain_order=True).agg(pl.all().sum()).sort(
        ...     "a"
        ... ).explain()  # doctest: +SKIP
        """
        if tree_format is not None:
            issue_deprecation_warning(
                "the `tree_format` parameter for `LazyFrame.explain` is deprecated"
                " Use the `format` parameter instead.",
                version="0.20.30",
            )
            if tree_format:
                format = "tree"

        engine = _select_engine(engine)

        if engine == "streaming":
            issue_unstable_warning("streaming mode is considered unstable.")

        if optimized:
            optimizations = optimizations.__copy__()
            optimizations._pyoptflags.streaming = engine == "streaming"
            ldf = self._ldf.with_optimizations(optimizations._pyoptflags)
            if format == "tree":
                return ldf.describe_optimized_plan_tree()
            else:
                return ldf.describe_optimized_plan()

        if format == "tree":
            return self._ldf.describe_plan_tree()
        else:
            return self._ldf.describe_plan()

    @deprecate_streaming_parameter()
    @forward_old_opt_flags()
    def show_graph(
        self,
        *,
        optimized: bool = True,
        show: bool = True,
        output_path: str | Path | None = None,
        raw_output: bool = False,
        figsize: tuple[float, float] = (16.0, 12.0),
        type_coercion: bool = True,
        _type_check: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        engine: EngineType = "auto",
        plan_stage: PlanStage = "ir",
        _check_order: bool = True,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> str | None:
        """
        Show a plot of the query plan.

        Note that Graphviz must be installed to render the visualization (if not
        already present, you can download it here: `<https://graphviz.org/download>`_).

        Parameters
        ----------
        optimized
            Optimize the query plan.
        show
            Show the figure.
        output_path
            Write the figure to disk.
        raw_output
            Return dot syntax. This cannot be combined with `show` and/or `output_path`.
        figsize
            Passed to matplotlib if `show == True`.
        type_coercion
            Do type coercion optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        predicate_pushdown
            Do predicate pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        projection_pushdown
            Do projection pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        simplify_expression
            Run simplify expressions optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        slice_pushdown
            Slice pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subplan_elim
            Will try to cache branching subplans that occur on self-joins or unions.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subexpr_elim
            Common subexpressions will be cached and reused.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        cluster_with_columns
            Combine sequential independent calls to with_columns.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        collapse_joins
            Collapse a join and filters into a faster join.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query
            is run using the polars in-memory engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars in-memory
            engine. If set to `"gpu"`, the GPU engine is used. Fine-grained
            control over the GPU engine, for example which device to use
            on a system with multiple devices, is possible by providing a
            :class:`~.GPUEngine` object with configuration options.

            .. note::
               GPU mode is considered **unstable**. Not all queries will run
               successfully on the GPU, however, they should fall back transparently
               to the default engine if execution is not supported.

               Running with `POLARS_VERBOSE=1` will provide information if a query
               falls back (and why).

            .. note::
               The GPU engine does not support streaming, if streaming
               is enabled then GPU execution is switched off.
        plan_stage : {'ir', 'physical'}
            Select the stage to display. Currently only the streaming engine has a
            separate physical stage, for the other engines both IR and physical are the
            same.
        optimizations
            The set of the optimizations considered during query optimization.


        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [1, 2, 3, 4, 5, 6],
        ...         "c": [6, 5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> lf.group_by("a", maintain_order=True).agg(pl.all().sum()).sort(
        ...     "a"
        ... ).show_graph()  # doctest: +SKIP
        """
        engine = _select_engine(engine)

        if engine == "streaming":
            issue_unstable_warning("streaming mode is considered unstable.")

        optimizations = optimizations.__copy__()
        optimizations._pyoptflags.streaming = engine == "streaming"
        _ldf = self._ldf.with_optimizations(optimizations._pyoptflags)

        if plan_stage == "ir":
            dot = _ldf.to_dot(optimized)
        elif plan_stage == "physical":
            if engine == "streaming":
                dot = _ldf.to_dot_streaming_phys(optimized)
            else:
                dot = _ldf.to_dot(optimized)
        else:
            error_msg = f"invalid plan stage '{plan_stage}'"
            raise TypeError(error_msg)

        return display_dot_graph(
            dot=dot,
            show=show,
            output_path=output_path,
            raw_output=raw_output,
            figsize=figsize,
        )

    def inspect(self, fmt: str = "{}") -> LazyFrame:
        """
        Inspect a node in the computation graph.

        Print the value that this node in the computation graph evaluates to and pass on
        the value.

        Examples
        --------
        >>> lf = pl.LazyFrame({"foo": [1, 1, -2, 3]})
        >>> (
        ...     lf.with_columns(pl.col("foo").cum_sum().alias("bar"))
        ...     .inspect()  # print the node before the filter
        ...     .filter(pl.col("bar") == pl.col("foo"))
        ... )
        <LazyFrame at ...>
        """

        def inspect(s: DataFrame) -> DataFrame:
            print(fmt.format(s))
            return s

        return self.map_batches(
            inspect, predicate_pushdown=True, projection_pushdown=True
        )

    def sort(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        *more_by: IntoExpr,
        descending: bool | Sequence[bool] = False,
        nulls_last: bool | Sequence[bool] = False,
        maintain_order: bool = False,
        multithreaded: bool = True,
    ) -> LazyFrame:
        """
        Sort the LazyFrame by the given columns.

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
        maintain_order
            Whether the order should be maintained if elements are equal.
            Note that if `true` streaming is not possible and performance might be
            worse since this requires a stable search.
        multithreaded
            Sort using multiple threads.

        Examples
        --------
        Pass a single column name to sort by that column.

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, None],
        ...         "b": [6.0, 5.0, 4.0],
        ...         "c": ["a", "c", "b"],
        ...     }
        ... )
        >>> lf.sort("a").collect()
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

        >>> lf.sort(pl.col("a") + pl.col("b") * 2, nulls_last=True).collect()
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

        >>> lf.sort(["c", "a"], descending=True).collect()
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

        >>> lf.sort("c", "a", descending=[False, True]).collect()
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
        # Fast path for sorting by a single existing column
        if (
            isinstance(by, str)
            and not more_by
            and isinstance(descending, bool)
            and isinstance(nulls_last, bool)
        ):
            return self._from_pyldf(
                self._ldf.sort(
                    by, descending, nulls_last, maintain_order, multithreaded
                )
            )

        by = parse_into_list_of_expressions(by, *more_by)
        descending = extend_bool(descending, len(by), "descending", "by")
        nulls_last = extend_bool(nulls_last, len(by), "nulls_last", "by")

        return self._from_pyldf(
            self._ldf.sort_by_exprs(
                by, descending, nulls_last, maintain_order, multithreaded
            )
        )

    def sql(self, query: str, *, table_name: str = "self") -> LazyFrame:
        """
        Execute a SQL query against the LazyFrame.

        .. versionadded:: 0.20.23

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

        See Also
        --------
        SQLContext

        Examples
        --------
        >>> lf1 = pl.LazyFrame({"a": [1, 2, 3], "b": [6, 7, 8], "c": ["z", "y", "x"]})
        >>> lf2 = pl.LazyFrame({"a": [3, 2, 1], "d": [125, -654, 888]})

        Query the LazyFrame using SQL:

        >>> lf1.sql("SELECT c, b FROM self WHERE a > 1").collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ c   ┆ b   │
        │ --- ┆ --- │
        │ str ┆ i64 │
        ╞═════╪═════╡
        │ y   ┆ 7   │
        │ x   ┆ 8   │
        └─────┴─────┘

        Apply SQL transforms (aliasing "self" to "frame") then filter
        natively (you can freely mix SQL and native operations):

        >>> lf1.sql(
        ...     query='''
        ...         SELECT
        ...             a,
        ...             (a % 2 == 0) AS a_is_even,
        ...             (b::float4 / 2) AS "b/2",
        ...             CONCAT_WS(':', c, c, c) AS c_c_c
        ...         FROM frame
        ...         ORDER BY a
        ...     ''',
        ...     table_name="frame",
        ... ).filter(~pl.col("c_c_c").str.starts_with("x")).collect()
        shape: (2, 4)
        ┌─────┬───────────┬─────┬───────┐
        │ a   ┆ a_is_even ┆ b/2 ┆ c_c_c │
        │ --- ┆ ---       ┆ --- ┆ ---   │
        │ i64 ┆ bool      ┆ f32 ┆ str   │
        ╞═════╪═══════════╪═════╪═══════╡
        │ 1   ┆ false     ┆ 3.0 ┆ z:z:z │
        │ 2   ┆ true      ┆ 3.5 ┆ y:y:y │
        └─────┴───────────┴─────┴───────┘
        """
        from polars.sql import SQLContext

        issue_unstable_warning(
            "`sql` is considered **unstable** (although it is close to being considered stable)."
        )
        with SQLContext(register_globals=False, eager=False) as ctx:
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
    ) -> LazyFrame:
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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [2, 1, 1, 3, 2, 1],
        ...     }
        ... )

        Get the rows which contain the 4 largest values in column b.

        >>> lf.top_k(4, by="b").collect()
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

        >>> lf.top_k(4, by=["b", "a"]).collect()
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
        by = parse_into_list_of_expressions(by)
        reverse = extend_bool(reverse, len(by), "reverse", "by")
        return self._from_pyldf(self._ldf.top_k(k, by=by, reverse=reverse))

    @deprecate_renamed_parameter("descending", "reverse", version="1.0.0")
    def bottom_k(
        self,
        k: int,
        *,
        by: IntoExpr | Iterable[IntoExpr],
        reverse: bool | Sequence[bool] = False,
    ) -> LazyFrame:
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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [2, 1, 1, 3, 2, 1],
        ...     }
        ... )

        Get the rows which contain the 4 smallest values in column b.

        >>> lf.bottom_k(4, by="b").collect()
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

        >>> lf.bottom_k(4, by=["a", "b"]).collect()
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
        by = parse_into_list_of_expressions(by)
        reverse = extend_bool(reverse, len(by), "reverse", "by")
        return self._from_pyldf(self._ldf.bottom_k(k, by=by, reverse=reverse))

    @forward_old_opt_flags()
    def profile(
        self,
        *,
        type_coercion: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        no_optimization: bool = False,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        show_plot: bool = False,
        truncate_nodes: int = 0,
        figsize: tuple[int, int] = (18, 8),
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
        **_kwargs: Any,
    ) -> tuple[DataFrame, DataFrame]:
        """
        Profile a LazyFrame.

        This will run the query and return a tuple
        containing the materialized DataFrame and a DataFrame that
        contains profiling information of each node that is executed.

        The units of the timings are microseconds.

        Parameters
        ----------
        type_coercion
            Do type coercion optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        predicate_pushdown
            Do predicate pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        projection_pushdown
            Do projection pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        simplify_expression
            Run simplify expressions optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        no_optimization
            Turn off (certain) optimizations.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        slice_pushdown
            Slice pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subplan_elim
            Will try to cache branching subplans that occur on self-joins or unions.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subexpr_elim
            Common subexpressions will be cached and reused.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        cluster_with_columns
            Combine sequential independent calls to with_columns

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        collapse_joins
            Collapse a join and filters into a faster join

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        show_plot
            Show a gantt chart of the profiling result
        truncate_nodes
            Truncate the label lengths in the gantt chart to this number of
            characters.
        figsize
            matplotlib figsize of the profiling plot
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query
            is run using the polars in-memory engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars in-memory
            engine. If set to `"gpu"`, the GPU engine is used. Fine-grained
            control over the GPU engine, for example which device to use
            on a system with multiple devices, is possible by providing a
            :class:`~.GPUEngine` object with configuration options.

            .. note::
               GPU mode is considered **unstable**. Not all queries will run
               successfully on the GPU, however, they should fall back transparently
               to the default engine if execution is not supported.

               Running with `POLARS_VERBOSE=1` will provide information if a query
               falls back (and why).

            .. note::
               The GPU engine does not support streaming, if streaming
               is enabled then GPU execution is switched off.
        optimizations
            The optimization passes done during query optimization.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.


        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [1, 2, 3, 4, 5, 6],
        ...         "c": [6, 5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> lf.group_by("a", maintain_order=True).agg(pl.all().sum()).sort(
        ...     "a"
        ... ).profile()  # doctest: +SKIP
        (shape: (3, 3)
         ┌─────┬─────┬─────┐
         │ a   ┆ b   ┆ c   │
         │ --- ┆ --- ┆ --- │
         │ str ┆ i64 ┆ i64 │
         ╞═════╪═════╪═════╡
         │ a   ┆ 4   ┆ 10  │
         │ b   ┆ 11  ┆ 10  │
         │ c   ┆ 6   ┆ 1   │
         └─────┴─────┴─────┘,
         shape: (3, 3)
         ┌─────────────────────────┬───────┬──────┐
         │ node                    ┆ start ┆ end  │
         │ ---                     ┆ ---   ┆ ---  │
         │ str                     ┆ u64   ┆ u64  │
         ╞═════════════════════════╪═══════╪══════╡
         │ optimization            ┆ 0     ┆ 5    │
         │ group_by_partitioned(a) ┆ 5     ┆ 470  │
         │ sort(a)                 ┆ 475   ┆ 1964 │
         └─────────────────────────┴───────┴──────┘)
        """
        for k in _kwargs:
            if k not in (  # except "private" kwargs
                "post_opt_callback",
            ):
                error_msg = f"profile() got an unexpected keyword argument '{k}'"
                raise TypeError(error_msg)
        engine = _select_engine(engine)

        optimizations = optimizations.__copy__()
        ldf = self._ldf.with_optimizations(optimizations._pyoptflags)

        callback = _gpu_engine_callback(
            engine,
            streaming=False,
            background=False,
            new_streaming=False,
            _eager=False,
        )
        if _kwargs.get("post_opt_callback") is not None:
            # Only for testing
            callback = _kwargs.get("post_opt_callback")
        df_py, timings_py = ldf.profile(callback)
        (df, timings) = wrap_df(df_py), wrap_df(timings_py)

        if show_plot:
            import_optional(
                "matplotlib",
                err_suffix="should be installed to show profiling plots",
            )
            import matplotlib.pyplot as plt

            _fig, ax = plt.subplots(1, figsize=figsize)

            max_val = timings["end"][-1]
            timings_ = timings.reverse()

            if max_val > 1e9:
                unit = "s"
                timings_ = timings_.with_columns(F.col(["start", "end"]) / 1_000_000)
            elif max_val > 1e6:
                unit = "ms"
                timings_ = timings_.with_columns(F.col(["start", "end"]) / 1000)
            else:
                unit = "us"
            if truncate_nodes > 0:
                timings_ = timings_.with_columns(
                    F.col("node").str.slice(0, truncate_nodes) + "..."
                )

            max_in_unit = timings_["end"][0]
            ax.barh(
                timings_["node"],
                width=timings_["end"] - timings_["start"],
                left=timings_["start"],
            )

            plt.title("Profiling result")
            ax.set_xlabel(f"node duration in [{unit}], total {max_in_unit}{unit}")
            ax.set_ylabel("nodes")
            plt.show()

        return df, timings

    @overload
    def collect(
        self,
        *,
        type_coercion: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        no_optimization: bool = False,
        engine: EngineType = "auto",
        background: Literal[True],
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> InProcessQuery: ...

    @overload
    def collect(
        self,
        *,
        type_coercion: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        no_optimization: bool = False,
        engine: EngineType = "auto",
        background: Literal[False] = False,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> DataFrame: ...

    @deprecate_streaming_parameter()
    @forward_old_opt_flags()
    def collect(
        self,
        *,
        type_coercion: bool = True,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        simplify_expression: bool = True,
        slice_pushdown: bool = True,
        comm_subplan_elim: bool = True,
        comm_subexpr_elim: bool = True,
        cluster_with_columns: bool = True,
        collapse_joins: bool = True,
        no_optimization: bool = False,
        engine: EngineType = "auto",
        background: bool = False,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
        **_kwargs: Any,
    ) -> DataFrame | InProcessQuery:
        """
        Materialize this LazyFrame into a DataFrame.

        By default, all query optimizations are enabled. Individual optimizations may
        be disabled by setting the corresponding parameter to `False`.

        Parameters
        ----------
        type_coercion
            Do type coercion optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        predicate_pushdown
            Do predicate pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        projection_pushdown
            Do projection pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        simplify_expression
            Run simplify expressions optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        slice_pushdown
            Slice pushdown optimization.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subplan_elim
            Will try to cache branching subplans that occur on self-joins or unions.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        comm_subexpr_elim
            Common subexpressions will be cached and reused.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        cluster_with_columns
            Combine sequential independent calls to with_columns

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        collapse_joins
            Collapse a join and filters into a faster join

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        no_optimization
            Turn off (certain) optimizations.

            .. deprecated:: 1.30.0
                Use the `optimizations` parameters.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query
            is run using the polars in-memory engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars in-memory
            engine. If set to `"gpu"`, the GPU engine is used. Fine-grained
            control over the GPU engine, for example which device to use
            on a system with multiple devices, is possible by providing a
            :class:`~.GPUEngine` object with configuration options.

            .. note::
               GPU mode is considered **unstable**. Not all queries will run
               successfully on the GPU, however, they should fall back transparently
               to the default engine if execution is not supported.

               Running with `POLARS_VERBOSE=1` will provide information if a query
               falls back (and why).

            .. note::
               The GPU engine does not support streaming, or running in the
               background. If either are enabled, then GPU execution is switched off.
        background
            Run the query in the background and get a handle to the query.
            This handle can be used to fetch the result or cancel the query.

            .. warning::
                Background mode is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        optimizations
            The optimization passes done during query optimization.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        DataFrame

        See Also
        --------
        explain : Print the query plan that is evaluated with collect.
        profile : Collect the LazyFrame and time each node in the computation graph.
        polars.collect_all : Collect multiple LazyFrames at the same time.
        polars.Config.set_streaming_chunk_size : Set the size of streaming batches.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [1, 2, 3, 4, 5, 6],
        ...         "c": [6, 5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> lf.group_by("a").agg(pl.all().sum()).collect()  # doctest: +SKIP
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 4   ┆ 10  │
        │ b   ┆ 11  ┆ 10  │
        │ c   ┆ 6   ┆ 1   │
        └─────┴─────┴─────┘

        Collect in streaming mode

        >>> lf.group_by("a").agg(pl.all().sum()).collect(
        ...     engine="streaming"
        ... )  # doctest: +SKIP
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 4   ┆ 10  │
        │ b   ┆ 11  ┆ 10  │
        │ c   ┆ 6   ┆ 1   │
        └─────┴─────┴─────┘

        Collect in GPU mode

        >>> lf.group_by("a").agg(pl.all().sum()).collect(engine="gpu")  # doctest: +SKIP
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 11  ┆ 10  │
        │ a   ┆ 4   ┆ 10  │
        │ c   ┆ 6   ┆ 1   │
        └─────┴─────┴─────┘

        With control over the device used

        >>> lf.group_by("a").agg(pl.all().sum()).collect(
        ...     engine=pl.GPUEngine(device=1)
        ... )  # doctest: +SKIP
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ b   ┆ 11  ┆ 10  │
        │ a   ┆ 4   ┆ 10  │
        │ c   ┆ 6   ┆ 1   │
        └─────┴─────┴─────┘
        """
        for k in _kwargs:
            if k not in (  # except "private" kwargs
                "new_streaming",
                "post_opt_callback",
            ):
                error_msg = f"collect() got an unexpected keyword argument '{k}'"
                raise TypeError(error_msg)

        engine = _select_engine(engine)

        new_streaming = (
            _kwargs.get("new_streaming", False) or get_engine_affinity() == "streaming"
        )

        if new_streaming:
            engine = "streaming"

        if engine == "streaming":
            issue_unstable_warning("streaming mode is considered unstable.")

        callback = _gpu_engine_callback(
            engine,
            streaming=False,
            background=background,
            new_streaming=new_streaming,
            _eager=optimizations._pyoptflags.eager,
        )

        if isinstance(engine, GPUEngine):
            engine = "gpu"

        ldf = self._ldf.with_optimizations(optimizations._pyoptflags)
        if background:
            issue_unstable_warning("background mode is considered unstable.")
            return InProcessQuery(ldf.collect_concurrently())

        # Only for testing purposes
        callback = _kwargs.get("post_opt_callback", callback)
        return wrap_df(ldf.collect(engine, callback))

    @overload
    def collect_async(
        self,
        *,
        gevent: Literal[True],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> _GeventDataFrameResult[DataFrame]: ...

    @overload
    def collect_async(
        self,
        *,
        gevent: Literal[False] = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> Awaitable[DataFrame]: ...

    @deprecate_streaming_parameter()
    def collect_async(
        self,
        *,
        gevent: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> Awaitable[DataFrame] | _GeventDataFrameResult[DataFrame]:
        """
        Collect DataFrame asynchronously in thread pool.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Collects into a DataFrame (like :func:`collect`) but, instead of returning
        a DataFrame directly, it is scheduled to be collected inside a thread pool,
        while this method returns almost instantly.

        This can be useful if you use `gevent` or `asyncio` and want to release
        control to other greenlets/tasks while LazyFrames are being collected.

        Parameters
        ----------
        gevent
            Return wrapper to `gevent.event.AsyncResult` instead of Awaitable
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query
            is run using the polars in-memory engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars in-memory
            engine.

            .. note::
               The GPU engine does not support async, or running in the
               background. If either are enabled, then GPU execution is switched off.
        optimizations
            The optimization passes done during query optimization.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        If `gevent=False` (default) then returns an awaitable.

        If `gevent=True` then returns wrapper that has a
        `.get(block=True, timeout=None)` method.

        See Also
        --------
        polars.collect_all : Collect multiple LazyFrames at the same time.
        polars.collect_all_async : Collect multiple LazyFrames at the same time lazily.

        Notes
        -----
        In case of error `set_exception` is used on
        `asyncio.Future`/`gevent.event.AsyncResult` and will be reraised by them.

        Examples
        --------
        >>> import asyncio
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "b", "c"],
        ...         "b": [1, 2, 3, 4, 5, 6],
        ...         "c": [6, 5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> async def main():
        ...     return await (
        ...         lf.group_by("a", maintain_order=True)
        ...         .agg(pl.all().sum())
        ...         .collect_async()
        ...     )
        >>> asyncio.run(main())
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ a   ┆ 4   ┆ 10  │
        │ b   ┆ 11  ┆ 10  │
        │ c   ┆ 6   ┆ 1   │
        └─────┴─────┴─────┘
        """
        engine = _select_engine(engine)

        if engine == "streaming":
            issue_unstable_warning("streaming mode is considered unstable.")

        ldf = self._ldf.with_optimizations(optimizations._pyoptflags)

        result: _GeventDataFrameResult[DataFrame] | _AioDataFrameResult[DataFrame] = (
            _GeventDataFrameResult() if gevent else _AioDataFrameResult()
        )
        ldf.collect_with_callback(engine, result._callback)
        return result

    def collect_schema(self) -> Schema:
        """
        Resolve the schema of this LazyFrame.

        Examples
        --------
        Determine the schema.

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.collect_schema()
        Schema({'foo': Int64, 'bar': Float64, 'ham': String})

        Access various properties of the schema.

        >>> schema = lf.collect_schema()
        >>> schema["bar"]
        Float64
        >>> schema.names()
        ['foo', 'bar', 'ham']
        >>> schema.dtypes()
        [Int64, Float64, String]
        >>> schema.len()
        3
        """
        return Schema(self._ldf.collect_schema(), check_dtypes=False)

    @overload
    def sink_parquet(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: str = "zstd",
        compression_level: int | None = None,
        statistics: bool | str | dict[str, bool] = True,
        row_group_size: int | None = None,
        data_page_size: int | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[False] = ...,
        field_overwrites: ParquetFieldOverwrites
        | Sequence[ParquetFieldOverwrites]
        | Mapping[str, ParquetFieldOverwrites]
        | None = None,
        engine: EngineType = "auto",
        metadata: ParquetMetadata | None = None,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> None: ...

    @overload
    def sink_parquet(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: str = "zstd",
        compression_level: int | None = None,
        statistics: bool | str | dict[str, bool] = True,
        row_group_size: int | None = None,
        data_page_size: int | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[True],
        field_overwrites: ParquetFieldOverwrites
        | Sequence[ParquetFieldOverwrites]
        | Mapping[str, ParquetFieldOverwrites]
        | None = None,
        engine: EngineType = "auto",
        metadata: ParquetMetadata | None = None,
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame: ...

    def sink_parquet(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: str = "zstd",
        compression_level: int | None = None,
        statistics: bool | str | dict[str, bool] = True,
        row_group_size: int | None = None,
        data_page_size: int | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        metadata: ParquetMetadata | None = None,
        mkdir: bool = False,
        lazy: bool = False,
        field_overwrites: ParquetFieldOverwrites
        | Sequence[ParquetFieldOverwrites]
        | Mapping[str, ParquetFieldOverwrites]
        | None = None,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame | None:
        """
        Evaluate the query in streaming mode and write to a Parquet file.

        This allows streaming results that are larger than RAM to be written to disk.

        Parameters
        ----------
        path
            File path to which the file should be written.
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
            - "full": calculate and write all available statistics.
            - `{ "statistic-key": True / False, ... }`. Available keys:

              - "min": column minimum value (default: `True`)
              - "max": column maximum value (default: `True`)
              - "distinct_count": number of unique column values (default: `False`)
              - "null_count": number of null values in column (default: `True`)
        row_group_size
            Size of the row groups in number of rows.
            If None (default), the chunks of the `DataFrame` are
            used. Writing in smaller chunks may reduce memory pressure and improve
            writing speeds.
        data_page_size
            Size limit of individual data pages.
            If not set defaults to 1024 * 1024 bytes
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
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
        sync_on_close: { None, 'data', 'all' }
            Sync to disk when before closing a file.

            * `None` does not sync.
            * `data` syncs the file contents.
            * `all` syncs the file contents and metadata.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
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
        lazy: bool
            Wait to start execution until `collect` is called.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        field_overwrites
            Property overwrites for individual Parquet fields.

            This allows more control over the writing process to the granularity of a
            Parquet field.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.
        optimizations
            The optimization passes done during query optimization.

            This has no effect if `lazy` is set to `True`.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> lf.sink_parquet("out.parquet")  # doctest: +SKIP

        Sink to a `BytesIO` object.

        >>> import io
        >>> buf = io.BytesIO()  # doctest: +SKIP
        >>> pl.LazyFrame({"x": [1, 2, 1]}).sink_parquet(buf)  # doctest: +SKIP

        Split into a hive-partitioning style partition:

        >>> pl.LazyFrame({"x": [1, 2, 1], "y": [3, 4, 5]}).sink_parquet(
        ...     pl.PartitionByKey("./out/", by="x"),
        ...     mkdir=True
        ... )  # doctest: +SKIP

        See Also
        --------
        PartitionByKey
        """
        engine = _select_engine(engine)
        if metadata is not None:
            msg = "`metadata` parameter is considered experimental"
            issue_unstable_warning(msg)

        if isinstance(statistics, bool) and statistics:
            statistics = {
                "min": True,
                "max": True,
                "distinct_count": False,
                "null_count": True,
            }
        elif isinstance(statistics, bool) and not statistics:
            statistics = {}
        elif statistics == "full":
            statistics = {
                "min": True,
                "max": True,
                "distinct_count": True,
                "null_count": True,
            }

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )

        credential_provider_builder = _init_credential_provider_builder(
            credential_provider, path, storage_options, "sink_parquet"
        )
        del credential_provider

        if storage_options:
            storage_options = list(storage_options.items())  # type: ignore[assignment]
        else:
            # Handle empty dict input
            storage_options = None

        target = _to_sink_target(path)
        sink_options = {
            "sync_on_close": sync_on_close or "none",
            "maintain_order": maintain_order,
            "mkdir": mkdir,
        }

        if isinstance(metadata, dict):
            if metadata:
                metadata = list(metadata.items())  # type: ignore[assignment]
            else:
                # Handle empty dict input
                metadata = None
        elif callable(metadata):
            metadata = wrap_parquet_metadata_callback(metadata)  # type: ignore[assignment]

        # Convert the field overwrites into something that can be ingested by Rust.
        field_overwrites_dicts: list[dict[str, Any]] = []
        if field_overwrites is not None:
            import collections

            from polars.io.parquet.field_overwrites import (
                ParquetFieldOverwrites,
                _parquet_field_overwrites_dict_to_dict_list,
                _parquet_field_overwrites_to_dict,
            )

            if isinstance(field_overwrites, ParquetFieldOverwrites):
                field_overwrites_dicts = [
                    _parquet_field_overwrites_to_dict(field_overwrites)
                ]
            elif isinstance(field_overwrites, collections.abc.Mapping):
                field_overwrites_dicts = _parquet_field_overwrites_dict_to_dict_list(
                    dict(field_overwrites)
                )
            elif isinstance(field_overwrites, collections.abc.Sequence):
                field_overwrites_dicts = [
                    _parquet_field_overwrites_to_dict(c) for c in field_overwrites
                ]
            else:
                msg = f"field_overwrites got the wrong type {type(field_overwrites)}"
                raise TypeError(msg)

        ldf_py = self._ldf.sink_parquet(
            target=target,
            compression=compression,
            compression_level=compression_level,
            statistics=statistics,
            row_group_size=row_group_size,
            data_page_size=data_page_size,
            cloud_options=storage_options,
            credential_provider=credential_provider_builder,
            retries=retries,
            sink_options=sink_options,
            metadata=metadata,
            field_overwrites=field_overwrites_dicts,
        )

        if not lazy:
            ldf_py = ldf_py.with_optimizations(optimizations._pyoptflags)
            ldf = LazyFrame._from_pyldf(ldf_py)
            ldf.collect(engine=engine)
            return None
        return LazyFrame._from_pyldf(ldf_py)

    @overload
    def sink_ipc(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: IpcCompression | None = "uncompressed",
        compat_level: CompatLevel | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[False] = ...,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> None: ...

    @overload
    def sink_ipc(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: IpcCompression | None = "uncompressed",
        compat_level: CompatLevel | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[True],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame: ...

    def sink_ipc(
        self,
        path: str | Path | IO[bytes] | PartitioningScheme,
        *,
        compression: IpcCompression | None = "uncompressed",
        compat_level: CompatLevel | None = None,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame | None:
        """
        Evaluate the query in streaming mode and write to an IPC file.

        This allows streaming results that are larger than RAM to be written to disk.

        Parameters
        ----------
        path
            File path to which the file should be written.
        compression : {'uncompressed', 'lz4', 'zstd'}
            Choose "zstd" for good compression performance.
            Choose "lz4" for fast compression/decompression.
        compat_level
            Use a specific compatibility level
            when exporting Polars' internal data structures.
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
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
        sync_on_close: { None, 'data', 'all' }
            Sync to disk when before closing a file.

            * `None` does not sync.
            * `data` syncs the file contents.
            * `all` syncs the file contents and metadata.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        mkdir: bool
            Recursively create all the directories in the path.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        lazy: bool
            Wait to start execution until `collect` is called.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.

            .. note::
               The GPU engine is currently not supported.
        optimizations
            The optimization passes done during query optimization.

            This has no effect if `lazy` is set to `True`.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> lf.sink_ipc("out.arrow")  # doctest: +SKIP

        Sink to a `BytesIO` object.

        >>> import io
        >>> buf = io.BytesIO()  # doctest: +SKIP
        >>> pl.LazyFrame({"x": [1, 2, 1]}).sink_ipc(buf)  # doctest: +SKIP

        Split into a hive-partitioning style partition:

        >>> pl.LazyFrame({"x": [1, 2, 1], "y": [3, 4, 5]}).sink_ipc(
        ...     pl.PartitionByKey("./out/", by="x"),
        ...     mkdir=True
        ... )  # doctest: +SKIP

        See Also
        --------
        PartitionByKey
        """
        engine = _select_engine(engine)

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )

        credential_provider_builder = _init_credential_provider_builder(
            credential_provider, path, storage_options, "sink_ipc"
        )
        del credential_provider

        if storage_options:
            storage_options = list(storage_options.items())  # type: ignore[assignment]
        else:
            # Handle empty dict input
            storage_options = None

        target = _to_sink_target(path)
        sink_options = {
            "sync_on_close": sync_on_close or "none",
            "maintain_order": maintain_order,
            "mkdir": mkdir,
        }

        compat_level_py: int | bool
        if compat_level is None:
            compat_level_py = True
        elif isinstance(compat_level, CompatLevel):
            compat_level_py = compat_level._version
        else:
            msg = f"`compat_level` has invalid type: {qualified_type_name(compat_level)!r}"
            raise TypeError(msg)

        if compression is None:
            compression = "uncompressed"

        ldf_py = self._ldf.sink_ipc(
            target=target,
            compression=compression,
            compat_level=compat_level_py,
            cloud_options=storage_options,
            credential_provider=credential_provider_builder,
            retries=retries,
            sink_options=sink_options,
        )

        if not lazy:
            ldf_py = ldf_py.with_optimizations(optimizations._pyoptflags)
            ldf = LazyFrame._from_pyldf(ldf_py)
            ldf.collect(engine=engine)
            return None
        return LazyFrame._from_pyldf(ldf_py)

    @overload
    def sink_csv(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
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
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[False] = ...,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> None: ...

    @overload
    def sink_csv(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
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
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[True],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame: ...

    def sink_csv(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
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
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame | None:
        """
        Evaluate the query in streaming mode and write to a CSV file.

        This allows streaming results that are larger than RAM to be written to disk.

        Parameters
        ----------
        path
            File path to which the file should be written.
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
            Use a comma as the decimal separator instead of a point. Floats will be
            encapsulated in quotes if necessary; set the field separator to override.
        null_value
            A string representing null values (defaulting to the empty string).
        quote_style : {'necessary', 'always', 'non_numeric', 'never'}
            Determines the quoting strategy used.

            - necessary (default): This puts quotes around fields only when necessary.
              They are necessary when fields contain a quote,
              delimiter or record terminator.
              Quotes are also necessary when writing an empty record
              (which is indistinguishable from a record with one empty field).
              This is the default.
            - always: This puts quotes around every field. Always.
            - never: This never puts quotes around fields, even if that results in
              invalid CSV data (e.g.: by not quoting strings containing the
              separator).
            - non_numeric: This puts quotes around all fields that are non-numeric.
              Namely, when writing a field that does not parse as a valid float
              or integer, then quotes will be used even if they aren`t strictly
              necessary.
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
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
        sync_on_close: { None, 'data', 'all' }
            Sync to disk when before closing a file.

            * `None` does not sync.
            * `data` syncs the file contents.
            * `all` syncs the file contents and metadata.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        mkdir: bool
            Recursively create all the directories in the path.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        lazy: bool
            Wait to start execution until `collect` is called.

            .. warning::
                This functionality is considered **unstable**. It may be changed at any
                point without it being considered a breaking change.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.
        optimizations
            The optimization passes done during query optimization.

            This has no effect if `lazy` is set to `True`.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> lf.sink_csv("out.csv")  # doctest: +SKIP

        Sink to a `BytesIO` object.

        >>> import io
        >>> buf = io.BytesIO()  # doctest: +SKIP
        >>> pl.LazyFrame({"x": [1, 2, 1]}).sink_csv(buf)  # doctest: +SKIP

        Split into a hive-partitioning style partition:

        >>> pl.LazyFrame({"x": [1, 2, 1], "y": [3, 4, 5]}).sink_csv(
        ...     pl.PartitionByKey("./out/", by="x"),
        ...     mkdir=True
        ... )  # doctest: +SKIP

        See Also
        --------
        PartitionByKey
        """
        from polars.io.csv._utils import _check_arg_is_1byte

        _check_arg_is_1byte("separator", separator, can_be_empty=False)
        _check_arg_is_1byte("quote_char", quote_char, can_be_empty=False)
        if not null_value:
            null_value = None
        engine = _select_engine(engine)

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )

        credential_provider_builder = _init_credential_provider_builder(
            credential_provider, path, storage_options, "sink_csv"
        )
        del credential_provider

        if storage_options:
            storage_options = list(storage_options.items())  # type: ignore[assignment]
        else:
            # Handle empty dict input
            storage_options = None

        target = _to_sink_target(path)
        sink_options = {
            "sync_on_close": sync_on_close or "none",
            "maintain_order": maintain_order,
            "mkdir": mkdir,
        }

        ldf_py = self._ldf.sink_csv(
            target=target,
            include_bom=include_bom,
            include_header=include_header,
            separator=ord(separator),
            line_terminator=line_terminator,
            quote_char=ord(quote_char),
            batch_size=batch_size,
            datetime_format=datetime_format,
            date_format=date_format,
            time_format=time_format,
            float_scientific=float_scientific,
            float_precision=float_precision,
            decimal_comma=decimal_comma,
            null_value=null_value,
            quote_style=quote_style,
            cloud_options=storage_options,
            credential_provider=credential_provider_builder,
            retries=retries,
            sink_options=sink_options,
        )

        if not lazy:
            ldf_py = ldf_py.with_optimizations(optimizations._pyoptflags)
            ldf = LazyFrame._from_pyldf(ldf_py)
            ldf.collect(engine=engine)
            return None
        return LazyFrame._from_pyldf(ldf_py)

    @overload
    def sink_ndjson(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
        *,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[False] = ...,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> None: ...

    @overload
    def sink_ndjson(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
        *,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: Literal[True],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame: ...

    def sink_ndjson(
        self,
        path: str | Path | IO[bytes] | IO[str] | PartitioningScheme,
        *,
        maintain_order: bool = True,
        storage_options: dict[str, Any] | None = None,
        credential_provider: CredentialProviderFunction
        | Literal["auto"]
        | None = "auto",
        retries: int = 2,
        sync_on_close: SyncOnCloseMethod | None = None,
        mkdir: bool = False,
        lazy: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> LazyFrame | None:
        """
        Evaluate the query in streaming mode and write to an NDJSON file.

        This allows streaming results that are larger than RAM to be written to disk.

        Parameters
        ----------
        path
            File path to which the file should be written.
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
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
        sync_on_close: { None, 'data', 'all' }
            Sync to disk when before closing a file.

            * `None` does not sync.
            * `data` syncs the file contents.
            * `all` syncs the file contents and metadata.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        mkdir: bool
            Recursively create all the directories in the path.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        lazy: bool
            Wait to start execution until `collect` is called.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.
        optimizations
            The optimization passes done during query optimization.

            This has no effect if `lazy` is set to `True`.

            .. warning::
                This functionality is considered **unstable**. It may be changed
                at any point without it being considered a breaking change.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> lf.sink_ndjson("out.ndjson")  # doctest: +SKIP

        Sink to a `BytesIO` object.

        >>> import io
        >>> buf = io.BytesIO()  # doctest: +SKIP
        >>> pl.LazyFrame({"x": [1, 2, 1]}).sink_ndjson(buf)  # doctest: +SKIP

        Split into a hive-partitioning style partition:

        >>> pl.LazyFrame({"x": [1, 2, 1], "y": [3, 4, 5]}).sink_ndjson(
        ...     pl.PartitionByKey("./out/", by="x"),
        ...     mkdir=True
        ... )  # doctest: +SKIP

        See Also
        --------
        PartitionByKey
        """
        engine = _select_engine(engine)

        from polars.io.cloud.credential_provider._builder import (
            _init_credential_provider_builder,
        )

        credential_provider_builder = _init_credential_provider_builder(
            credential_provider, path, storage_options, "sink_ndjson"
        )
        del credential_provider

        if storage_options:
            storage_options = list(storage_options.items())  # type: ignore[assignment]
        else:
            # Handle empty dict input
            storage_options = None

        target = _to_sink_target(path)
        sink_options = {
            "sync_on_close": sync_on_close or "none",
            "maintain_order": maintain_order,
            "mkdir": mkdir,
        }

        ldf_py = self._ldf.sink_json(
            target=target,
            cloud_options=storage_options,
            credential_provider=credential_provider_builder,
            retries=retries,
            sink_options=sink_options,
        )

        if not lazy:
            ldf_py = ldf_py.with_optimizations(optimizations._pyoptflags)
            ldf = LazyFrame._from_pyldf(ldf_py)
            ldf.collect(engine=engine)
            return None
        return LazyFrame._from_pyldf(ldf_py)

    @overload
    def sink_batches(
        self,
        function: Callable[[DataFrame], bool | None],
        *,
        chunk_size: int | None = None,
        maintain_order: bool = True,
        lazy: Literal[False],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> None: ...
    @overload
    def sink_batches(
        self,
        function: Callable[[DataFrame], bool | None],
        *,
        chunk_size: int | None = None,
        maintain_order: bool = True,
        lazy: Literal[True],
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> pl.LazyFrame: ...
    @unstable()
    def sink_batches(
        self,
        function: Callable[[DataFrame], bool | None],
        *,
        chunk_size: int | None = None,
        maintain_order: bool = True,
        lazy: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> pl.LazyFrame | None:
        """
        Evaluate the query and call a user-defined function for every ready batch.

        This allows streaming results that are larger than RAM in certain cases.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. warning::
            This method is much slower than native sinks. Only use it if you cannot
            implement your logic otherwise.

        Parameters
        ----------
        function
            Function to run with a batch that is ready. If the function returns
            `True`, this signals that no more results are needed, allowing for
            early stopping.
        chunk_size
            The number of rows that are buffered before the callback is called.
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.
        lazy: bool
            Wait to start execution until `collect` is called.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.
        optimizations
            The optimization passes done during query optimization.

            This has no effect if `lazy` is set to `True`.

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> lf.sink_batches(lambda df: print(df))  # doctest: +SKIP
        """

        def _wrap(pydf: plr.PyDataFrame) -> bool:
            df = wrap_df(pydf)
            return bool(function(df))

        ldf = self._ldf.sink_batches(
            function=_wrap,
            maintain_order=maintain_order,
            chunk_size=chunk_size,
        )

        if not lazy:
            ldf = ldf.with_optimizations(optimizations._pyoptflags)
            lf = LazyFrame._from_pyldf(ldf)
            lf.collect(engine=engine)
            return None
        return LazyFrame._from_pyldf(ldf)

    @unstable()
    def collect_batches(
        self,
        *,
        chunk_size: int | None = None,
        maintain_order: bool = True,
        lazy: bool = False,
        engine: EngineType = "auto",
        optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    ) -> Iterator[DataFrame]:
        """
        Evaluate the query in streaming mode and get a generator that returns chunks.

        This allows streaming results that are larger than RAM to be written to disk.

        The query will always be fully executed unless `stop` is called, so you should
        call next until all chunks have been seen.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. warning::
            This method is much slower than native sinks. Only use it if you cannot
            implement your logic otherwise.

        Parameters
        ----------
        chunk_size
            The number of rows that are buffered before a chunk is given.
        maintain_order
            Maintain the order in which data is processed.
            Setting this to `False` will be slightly faster.
        lazy
            Start the query when first requesting a batch.
        engine
            Select the engine used to process the query, optional.
            At the moment, if set to `"auto"` (default), the query is run
            using the polars streaming engine. Polars will also
            attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
            environment variable. If it cannot run the query using the
            selected engine, the query is run using the polars streaming
            engine.
        optimizations
            The optimization passes done during query optimization.

        Examples
        --------
        >>> lf = pl.scan_csv("/path/to/my_larger_than_ram_file.csv")  # doctest: +SKIP
        >>> for df in lf.collect_batches():
        ...     print(df)  # doctest: +SKIP
        """
        from queue import Queue

        class BatchCollector:
            def __init__(
                self,
                *,
                lf: pl.LazyFrame,
                chunk_size: int | None,
                maintain_order: bool,
                lazy: bool,
                engine: EngineType,
                optimizations: QueryOptFlags,
            ) -> None:
                class SharedState:
                    def __init__(self) -> None:
                        self.queue: Queue[pl.DataFrame | None] = Queue(maxsize=1)
                        self.stopped = False

                self._lf = lf
                self._chunk_size = chunk_size
                self._maintain_order = maintain_order
                self._engine = engine
                self._optimizations = optimizations
                self._shared = SharedState()
                self._fut: Future[None] | None = None

                if not lazy:
                    self._start()

            def _start(self) -> None:
                if self._fut is not None:
                    return

                # Make sure we don't capture self which would cause __del__
                # to not get called.
                shared = self._shared
                chunk_size = self._chunk_size
                maintain_order = self._maintain_order
                engine = self._engine
                optimizations = self._optimizations
                lf = self._lf

                def task() -> None:
                    def put_batch_in_queue(df: DataFrame) -> bool | None:
                        if shared.stopped:
                            return True
                        shared.queue.put(df)
                        return shared.stopped

                    try:
                        lf.sink_batches(
                            put_batch_in_queue,
                            chunk_size=chunk_size,
                            maintain_order=maintain_order,
                            engine=engine,
                            optimizations=optimizations,
                            lazy=False,
                        )
                    finally:
                        shared.queue.put(None)  # Signal the end of batches.

                self._fut = _COLLECT_BATCHES_POOL.submit(task)

            def __iter__(self) -> BatchCollector:
                return self

            def __next__(self) -> DataFrame:
                if self._shared.stopped:
                    raise StopIteration

                self._start()
                df = self._shared.queue.get()
                if df is None:
                    self._shared.stopped = True
                    raise StopIteration

                return df

            def __del__(self) -> None:
                if not self._shared.stopped:
                    # Signal to stop and unblock sink_batches task.
                    self._shared.stopped = True
                    while self._shared.queue.get() is not None:
                        pass
                if self._fut is not None:
                    self._fut.result()

        return BatchCollector(
            lf=self,
            chunk_size=chunk_size,
            maintain_order=maintain_order,
            lazy=lazy,
            engine=engine,
            optimizations=optimizations,
        )

    @deprecated(
        "`LazyFrame.fetch` is deprecated; use `LazyFrame.collect` "
        "instead, in conjunction with a call to `head`."
    )
    def fetch(
        self,
        n_rows: int = 500,
        **kwargs: Any,
    ) -> DataFrame:
        """
        Collect a small number of rows for debugging purposes.

        .. deprecated:: 1.0
            Use :meth:`collect` instead, in conjunction with a call to :meth:`head`.`

        Notes
        -----
        This is similar to a :func:`collect` operation, but it overwrites the number of
        rows read by *every* scan operation. Be aware that `fetch` does not guarantee
        the final number of rows in the DataFrame. Filters, join operations and fewer
        rows being available in the scanned data will all influence the final number
        of rows (joins are especially susceptible to this, and may return no data
        at all if `n_rows` is too small as the join keys may not be present).

        Warnings
        --------
        This is strictly a utility function that can help to debug queries using a
        smaller number of rows, and should *not* be used in production code.
        """
        return self.head(n_rows).collect(**kwargs)

    def lazy(self) -> LazyFrame:
        """
        Return lazy representation, i.e. itself.

        Useful for writing code that expects either a :class:`DataFrame` or
        :class:`LazyFrame`. On LazyFrame this is a no-op, and returns the same object.

        Returns
        -------
        LazyFrame

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [None, 2, 3, 4],
        ...         "b": [0.5, None, 2.5, 13],
        ...         "c": [True, True, False, None],
        ...     }
        ... )
        >>> lf.lazy()
        <LazyFrame at ...>
        """
        return self

    def cache(self) -> LazyFrame:
        """
        Cache the result once the execution of the physical plan hits this node.

        It is not recommended using this as the optimizer likely can do a better job.
        """
        return self._from_pyldf(self._ldf.cache())

    def cast(
        self,
        dtypes: (
            Mapping[
                ColumnNameOrSelector | PolarsDataType, PolarsDataType | PythonDataType
            ]
            | PolarsDataType
            | pl.DataTypeExpr
        ),
        *,
        strict: bool = True,
    ) -> LazyFrame:
        """
        Cast LazyFrame column(s) to the specified dtype(s).

        Parameters
        ----------
        dtypes
            Mapping of column names (or selector) to dtypes, or a single dtype
            to which all columns will be cast.
        strict
            Throw an error if a cast could not be done (for instance, due to an
            overflow).

        Examples
        --------
        >>> from datetime import date
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": [date(2020, 1, 2), date(2021, 3, 4), date(2022, 5, 6)],
        ...     }
        ... )

        Cast specific frame columns to the specified dtypes:

        >>> lf.cast({"foo": pl.Float32, "bar": pl.UInt8}).collect()
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

        >>> lf.cast({pl.Date: pl.Datetime}).collect()
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
        >>> lf.cast({cs.numeric(): pl.UInt32, cs.temporal(): pl.String}).collect()
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

        >>> lf.cast(pl.String).collect().to_dict(as_series=False)
        {'foo': ['1', '2', '3'],
         'bar': ['6.0', '7.0', '8.0'],
         'ham': ['2020-01-02', '2021-03-04', '2022-05-06']}
        """
        if not isinstance(dtypes, Mapping):
            dtypes = parse_into_datatype_expr(dtypes)
            return self._from_pyldf(self._ldf.cast_all(dtypes._pydatatype_expr, strict))

        cast_map = {}
        for c, dtype in dtypes.items():
            if (is_polars_dtype(c) or isinstance(c, DataTypeGroup)) or (
                isinstance(c, Collection) and all(is_polars_dtype(x) for x in c)
            ):
                c = by_dtype(c)  # type: ignore[arg-type]

            dtype = parse_into_dtype(dtype)
            cast_map.update(
                {c: dtype}
                if isinstance(c, str)
                else dict.fromkeys(expand_selector(self, c), dtype)  # type: ignore[arg-type]
            )

        return self._from_pyldf(self._ldf.cast(cast_map, strict))

    def clear(self, n: int = 0) -> LazyFrame:
        """
        Create an empty copy of the current LazyFrame, with zero to 'n' rows.

        Returns a copy with an identical schema but no data.

        Parameters
        ----------
        n
            Number of (empty) rows to return in the cleared frame.

        See Also
        --------
        clone : Cheap deepcopy/clone.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [None, 2, 3, 4],
        ...         "b": [0.5, None, 2.5, 13],
        ...         "c": [True, True, False, None],
        ...     }
        ... )
        >>> lf.clear().collect()
        shape: (0, 3)
        ┌─────┬─────┬──────┐
        │ a   ┆ b   ┆ c    │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ f64 ┆ bool │
        ╞═════╪═════╪══════╡
        └─────┴─────┴──────┘

        >>> lf.clear(2).collect()
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
        return pl.DataFrame(schema=self.collect_schema()).clear(n).lazy()

    def clone(self) -> LazyFrame:
        """
        Create a copy of this LazyFrame.

        This is a cheap operation that does not copy data.

        See Also
        --------
        clear : Create an empty copy of the current LazyFrame, with identical
            schema but no data.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [None, 2, 3, 4],
        ...         "b": [0.5, None, 2.5, 13],
        ...         "c": [True, True, False, None],
        ...     }
        ... )
        >>> lf.clone()
        <LazyFrame at ...>
        """
        return self._from_pyldf(self._ldf.clone())

    def _filter(
        self,
        *,
        predicates: tuple[
            IntoExprColumn
            | Iterable[IntoExprColumn]
            | bool
            | list[bool]
            | np.ndarray[Any, Any],
            ...,
        ],
        constraints: dict[str, Any],
        invert: bool = False,
    ) -> LazyFrame:
        """Common code for filter/remove ops."""
        all_predicates: list[pl.Expr] = []
        boolean_masks = []

        for p in predicates:
            # quick exit/skip conditions
            if (p is False and invert) or (p is True and not invert):
                continue  # ignore; doesn't filter/remove anything
            if (p is True and invert) or (p is False and not invert):
                return self.clear()  # discard all rows

            if _is_generator(p):
                p = tuple(p)

            # note: identify masks separately from predicates
            if is_bool_sequence(p, include_series=True):
                boolean_masks.append(pl.Series(p, dtype=Boolean))
            elif (
                (is_seq := is_sequence(p))
                and any(not isinstance(x, pl.Expr) for x in p)
            ) or (
                not is_seq
                and not isinstance(p, pl.Expr)
                and not (isinstance(p, str) and p in self.collect_schema())
            ):
                err = (
                    f"Series(…, dtype={p.dtype})"
                    if isinstance(p, pl.Series)
                    else repr(p)
                )
                msg = f"invalid predicate for `filter`: {err}"
                raise TypeError(msg)
            else:
                all_predicates.extend(
                    wrap_expr(x) for x in parse_into_list_of_expressions(p)
                )

        # unpack equality constraints from kwargs
        all_predicates.extend(
            F.col(name).eq(value) for name, value in constraints.items()
        )
        if not (all_predicates or boolean_masks):
            msg = "at least one predicate or constraint must be provided"
            raise TypeError(msg)

        # if multiple predicates, combine as 'horizontal' expression
        combined_predicate = (
            (
                F.all_horizontal(*all_predicates)
                if len(all_predicates) > 1
                else all_predicates[0]
            )
            if all_predicates
            else None
        )

        # apply reduced boolean mask first, if applicable, then predicates
        if boolean_masks:
            mask_expr = F.lit(reduce(and_, boolean_masks))
            combined_predicate = (
                mask_expr
                if combined_predicate is None
                else mask_expr & combined_predicate
            )

        if combined_predicate is None:
            return self._from_pyldf(self._ldf)

        filter_method = self._ldf.remove if invert else self._ldf.filter
        return self._from_pyldf(filter_method(combined_predicate._pyexpr))

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
    ) -> LazyFrame:
        """
        Filter rows in the LazyFrame based on a predicate expression.

        The original order of the remaining rows is preserved.

        Rows where the filter predicate does not evaluate to True are discarded
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
        False. As a result, these rows will not be retained. Ensure that null values
        are handled appropriately to avoid unexpected behaviour (see examples below).

        See Also
        --------
        remove

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3, None, 4, None, 0],
        ...         "bar": [6, 7, 8, None, None, 9, 0],
        ...         "ham": ["a", "b", "c", None, "d", "e", "f"],
        ...     }
        ... )

        Filter on one condition:

        >>> lf.filter(pl.col("foo") > 1).collect()
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

        Filter on multiple conditions:

        >>> lf.filter((pl.col("foo") < 3) & (pl.col("ham") == "a")).collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘

        Provide multiple filters using `*args` syntax:

        >>> lf.filter(
        ...     pl.col("foo") == 1,
        ...     pl.col("ham") == "a",
        ... ).collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘

        Provide multiple filters using `**kwargs` syntax:

        >>> lf.filter(foo=1, ham="a").collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        └─────┴─────┴─────┘

        Filter on an OR condition:

        >>> lf.filter(
        ...     (pl.col("foo") == 1) | (pl.col("ham") == "c"),
        ... ).collect()
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6   ┆ a   │
        │ 3   ┆ 8   ┆ c   │
        └─────┴─────┴─────┘

        Filter by comparing two columns against each other

        >>> lf.filter(
        ...     pl.col("foo") == pl.col("bar"),
        ... ).collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 0   ┆ 0   ┆ f   │
        └─────┴─────┴─────┘

        >>> lf.filter(
        ...     pl.col("foo") != pl.col("bar"),
        ... ).collect()
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

        Notice how the row with `None` values is filtered out; using `ne_missing`
        ensures that null values compare equal, and we get similar behaviour to Pandas:

        >>> lf.filter(
        ...     pl.col("foo").ne_missing(pl.col("bar")),
        ... ).collect()
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
        if not constraints:
            # early-exit conditions (exclude/include all rows)
            if not predicates or (len(predicates) == 1 and predicates[0] is True):
                return self.clone()
            if len(predicates) == 1 and predicates[0] is False:
                return self.clear()

        return self._filter(
            predicates=predicates,
            constraints=constraints,
            invert=False,
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
    ) -> LazyFrame:
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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [2, 3, None, 4, 0],
        ...         "bar": [5, 6, None, None, 0],
        ...         "ham": ["a", "b", None, "c", "d"],
        ...     }
        ... )

        Remove rows matching a condition:

        >>> lf.remove(
        ...     pl.col("bar") >= 5,
        ... ).collect()
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

        >>> lf.remove(
        ...     (pl.col("foo") >= 0) & (pl.col("bar") >= 0),
        ... ).collect()
        shape: (2, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        │ 4    ┆ null ┆ c    │
        └──────┴──────┴──────┘

        >>> lf.remove(
        ...     (pl.col("foo") >= 0) | (pl.col("bar") >= 0),
        ... ).collect()
        shape: (1, 3)
        ┌──────┬──────┬──────┐
        │ foo  ┆ bar  ┆ ham  │
        │ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╡
        │ null ┆ null ┆ null │
        └──────┴──────┴──────┘

        Provide multiple constraints using `*args` syntax:

        >>> lf.remove(
        ...     pl.col("ham").is_not_null(),
        ...     pl.col("bar") >= 0,
        ... ).collect()
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

        >>> lf.remove(foo=0, bar=0).collect()
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

        Remove rows by comparing two columns against each other; in this case, we
        remove rows where the two columns are not equal (using `ne_missing` to
        ensure that null values compare equal):

        >>> lf.remove(
        ...     pl.col("foo").ne_missing(pl.col("bar")),
        ... ).collect()
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
        if not constraints:
            # early-exit conditions (exclude/include all rows)
            if not predicates or (len(predicates) == 1 and predicates[0] is True):
                return self.clear()
            if len(predicates) == 1 and predicates[0] is False:
                return self.clone()

        return self._filter(
            predicates=predicates,
            constraints=constraints,
            invert=True,
        )

    def select(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> LazyFrame:
        """
        Select columns from this LazyFrame.

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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.select("foo").collect()
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

        >>> lf.select(["foo", "bar"]).collect()
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

        >>> lf.select(pl.col("foo"), pl.col("bar") + 1).collect()
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

        >>> lf.select(
        ...     threshold=pl.when(pl.col("foo") > 2).then(10).otherwise(0)
        ... ).collect()
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
        structify = bool(int(os.environ.get("POLARS_AUTO_STRUCTIFY", 0)))

        pyexprs = parse_into_list_of_expressions(
            *exprs, **named_exprs, __structify=structify
        )
        return self._from_pyldf(self._ldf.select(pyexprs))

    def select_seq(
        self, *exprs: IntoExpr | Iterable[IntoExpr], **named_exprs: IntoExpr
    ) -> LazyFrame:
        """
        Select columns from this LazyFrame.

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
        structify = bool(int(os.environ.get("POLARS_AUTO_STRUCTIFY", 0)))

        pyexprs = parse_into_list_of_expressions(
            *exprs, **named_exprs, __structify=structify
        )
        return self._from_pyldf(self._ldf.select_seq(pyexprs))

    def group_by(
        self,
        *by: IntoExpr | Iterable[IntoExpr],
        maintain_order: bool = False,
        **named_by: IntoExpr,
    ) -> LazyGroupBy:
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
            Setting this to `True` blocks the possibility
            to run on the streaming engine.
        **named_by
            Additional columns to group by, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        Examples
        --------
        Group by one column and call `agg` to compute the grouped sum of another
        column.

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["a", "b", "a", "b", "c"],
        ...         "b": [1, 2, 1, 3, 3],
        ...         "c": [5, 4, 3, 2, 1],
        ...     }
        ... )
        >>> lf.group_by("a").agg(pl.col("b").sum()).collect()  # doctest: +IGNORE_RESULT
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

        >>> lf.group_by("a", maintain_order=True).agg(pl.col("c")).collect()
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

        >>> lf.group_by(["a", "b"]).agg(pl.max("c")).collect()  # doctest: +SKIP
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

        >>> lf.group_by("a", pl.col("b") // 2).agg(
        ...     pl.col("c").mean()
        ... ).collect()  # doctest: +SKIP
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
        exprs = parse_into_list_of_expressions(*by, **named_by)
        lgb = self._ldf.group_by(exprs, maintain_order)
        return LazyGroupBy(lgb)

    @deprecate_renamed_parameter("by", "group_by", version="0.20.14")
    def rolling(
        self,
        index_column: IntoExpr,
        *,
        period: str | timedelta,
        offset: str | timedelta | None = None,
        closed: ClosedInterval = "right",
        group_by: IntoExpr | Iterable[IntoExpr] | None = None,
    ) -> LazyGroupBy:
        """
        Create rolling groups based on a temporal or integer column.

        Different from a `group_by_dynamic` the windows are now determined by the
        individual values and are not of constant intervals. For constant intervals
        use :func:`LazyFrame.group_by_dynamic`.

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

            In case of a rolling group by on indices, dtype needs to be one of
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
        LazyGroupBy
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
        >>> df = pl.LazyFrame({"dt": dates, "a": [3, 7, 5, 9, 2, 1]}).with_columns(
        ...     pl.col("dt").str.strptime(pl.Datetime).set_sorted()
        ... )
        >>> out = (
        ...     df.rolling(index_column="dt", period="2d")
        ...     .agg(
        ...         pl.sum("a").alias("sum_a"),
        ...         pl.min("a").alias("min_a"),
        ...         pl.max("a").alias("max_a"),
        ...     )
        ...     .collect()
        ... )
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
        """
        index_column_py = parse_into_expression(index_column)
        if offset is None:
            offset = negate_duration_string(parse_as_duration_string(period))

        pyexprs_by = (
            parse_into_list_of_expressions(group_by) if group_by is not None else []
        )
        period = parse_as_duration_string(period)
        offset = parse_as_duration_string(offset)

        lgb = self._ldf.rolling(index_column_py, period, offset, closed, pyexprs_by)
        return LazyGroupBy(lgb)

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
    ) -> LazyGroupBy:
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
        LazyGroupBy
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
        >>> lf = pl.LazyFrame(
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
        >>> lf.collect()
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

        >>> lf.group_by_dynamic("time", every="1h", closed="right").agg(
        ...     pl.col("n")
        ... ).collect()
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

        >>> lf.group_by_dynamic(
        ...     "time", every="1h", include_boundaries=True, closed="right"
        ... ).agg(pl.col("n").mean()).collect()
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

        >>> lf.group_by_dynamic("time", every="1h", closed="left").agg(
        ...     pl.col("n")
        ... ).collect()
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

        >>> lf.group_by_dynamic("time", every="1h", closed="both").agg(
        ...     pl.col("n")
        ... ).collect()
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

        >>> lf = lf.with_columns(groups=pl.Series(["a", "a", "a", "b", "b", "a", "a"]))
        >>> lf.collect()
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
        >>> lf.group_by_dynamic(
        ...     "time",
        ...     every="1h",
        ...     closed="both",
        ...     group_by="groups",
        ...     include_boundaries=True,
        ... ).agg(pl.col("n")).collect()
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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "idx": pl.int_range(0, 6, eager=True),
        ...         "A": ["A", "A", "B", "B", "B", "C"],
        ...     }
        ... )
        >>> lf.group_by_dynamic(
        ...     "idx",
        ...     every="2i",
        ...     period="3i",
        ...     include_boundaries=True,
        ...     closed="right",
        ... ).agg(pl.col("A").alias("A_agg_list")).collect()
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
        index_column_py = parse_into_expression(index_column)
        if offset is None:
            offset = "0ns"

        if period is None:
            period = every

        period = parse_as_duration_string(period)
        offset = parse_as_duration_string(offset)
        every = parse_as_duration_string(every)

        pyexprs_by = (
            parse_into_list_of_expressions(group_by) if group_by is not None else []
        )
        lgb = self._ldf.group_by_dynamic(
            index_column_py,
            every,
            period,
            offset,
            label,
            include_boundaries,
            closed,
            pyexprs_by,
            start_by,
        )
        return LazyGroupBy(lgb)

    def join_asof(
        self,
        other: LazyFrame,
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
    ) -> LazyFrame:
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

            A "nearest" search selects the last row in the right DataFrame whose value
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
            Join on these columns before doing asof join.
        by_left
            Join on these columns before doing asof join.
        by_right
            Join on these columns before doing asof join.
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

            - True: -> Always coalesce join columns.
            - False: -> Never coalesce join columns.

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
        >>> gdp = pl.LazyFrame(
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
        >>> gdp.collect()
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

        >>> population = pl.LazyFrame(
        ...     {
        ...         "date": [date(2016, 3, 1), date(2018, 8, 1), date(2019, 1, 1)],
        ...         "population": [82.19, 82.66, 83.12],
        ...     }
        ... ).sort("date")
        >>> population.collect()
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

        >>> population.join_asof(gdp, on="date", strategy="backward").collect()
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

        >>> population.join_asof(
        ...     gdp, on="date", strategy="backward", coalesce=False
        ... ).collect()
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

        >>> population.join_asof(gdp, on="date", strategy="forward").collect()
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

        >>> population.join_asof(gdp, on="date", strategy="nearest").collect()
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
        >>> gdp2 = pl.LazyFrame(
        ...     {
        ...         "country": ["Germany"] * 5 + ["Netherlands"] * 5,
        ...         "date": pl.concat([gdp_dates, gdp_dates]),
        ...         "gdp": [4164, 4411, 4566, 4696, 4827, 784, 833, 914, 910, 909],
        ...     }
        ... ).sort("country", "date")
        >>>
        >>> gdp2.collect()
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
        >>> pop2 = pl.LazyFrame(
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
        >>> pop2.collect()
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
        >>> pop2.join_asof(gdp2, by="country", on="date", strategy="nearest").collect()
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

        if isinstance(on, (str, pl.Expr)):
            left_on = on
            right_on = on

        if left_on is None or right_on is None:
            msg = "you should pass the column to join on as an argument"
            raise ValueError(msg)

        if by is not None:
            by_left_ = [by] if isinstance(by, str) else by
            by_right_ = by_left_
        elif (by_left is not None) or (by_right is not None):
            by_left_ = [by_left] if isinstance(by_left, str) else by_left  # type: ignore[assignment]
            by_right_ = [by_right] if isinstance(by_right, str) else by_right  # type: ignore[assignment]

        else:
            # no by
            by_left_ = None
            by_right_ = None

        tolerance_str: str | None = None
        tolerance_num: float | int | None = None
        if isinstance(tolerance, str):
            tolerance_str = tolerance
        elif isinstance(tolerance, timedelta):
            tolerance_str = parse_as_duration_string(tolerance)
        else:
            tolerance_num = tolerance

        if not isinstance(left_on, pl.Expr):
            left_on = F.col(left_on)
        if not isinstance(right_on, pl.Expr):
            right_on = F.col(right_on)

        return self._from_pyldf(
            self._ldf.join_asof(
                other._ldf,
                left_on._pyexpr,
                right_on._pyexpr,
                by_left_,
                by_right_,
                allow_parallel,
                force_parallel,
                suffix,
                strategy,
                tolerance_num,
                tolerance_str,
                coalesce=coalesce,
                allow_eq=allow_exact_matches,
                check_sortedness=check_sortedness,
            )
        )

    @deprecate_renamed_parameter("join_nulls", "nulls_equal", version="1.24")
    def join(
        self,
        other: LazyFrame,
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
        allow_parallel: bool = True,
        force_parallel: bool = False,
    ) -> LazyFrame:
        """
        Add a join operation to the Logical Plan.

        .. versionchanged:: 1.24
            The `join_nulls` parameter was renamed `nulls_equal`.

        Parameters
        ----------
        other
            Lazy DataFrame to join with.
        on
            Name(s) of the join columns in both DataFrames. If set, `left_on` and
            `right_on` should be None. This should not be specified if `how='cross'`.
        how : {'inner','left', 'right', 'full', 'semi', 'anti', 'cross'}
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
            Join column of the left DataFrame.
        right_on
            Join column of the right DataFrame.
        suffix
            Suffix to append to columns with a duplicate name.
        validate: {'m:m', 'm:1', '1:m', '1:1'}
            Checks if join is of specified type.

            .. list-table ::
               :header-rows: 0

               * - **m:m**
                 - *(Default)* Many-to-many. Does not result in checks.
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

        allow_parallel
            Allow the physical plan to optionally evaluate the computation of both
            DataFrames up to the join in parallel.
        force_parallel
            Force the physical plan to evaluate the computation of both DataFrames up to
            the join in parallel.

        See Also
        --------
        join_asof

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> other_lf = pl.LazyFrame(
        ...     {
        ...         "apple": ["x", "y", "z"],
        ...         "ham": ["a", "b", "d"],
        ...     }
        ... )
        >>> lf.join(other_lf, on="ham").collect()
        shape: (2, 4)
        ┌─────┬─────┬─────┬───────┐
        │ foo ┆ bar ┆ ham ┆ apple │
        │ --- ┆ --- ┆ --- ┆ ---   │
        │ i64 ┆ f64 ┆ str ┆ str   │
        ╞═════╪═════╪═════╪═══════╡
        │ 1   ┆ 6.0 ┆ a   ┆ x     │
        │ 2   ┆ 7.0 ┆ b   ┆ y     │
        └─────┴─────┴─────┴───────┘
        >>> lf.join(other_lf, on="ham", how="full").collect()
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
        >>> lf.join(other_lf, on="ham", how="left", coalesce=True).collect()
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
        >>> lf.join(other_lf, on="ham", how="semi").collect()
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 6.0 ┆ a   │
        │ 2   ┆ 7.0 ┆ b   │
        └─────┴─────┴─────┘
        >>> lf.join(other_lf, on="ham", how="anti").collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ f64 ┆ str │
        ╞═════╪═════╪═════╡
        │ 3   ┆ 8.0 ┆ c   │
        └─────┴─────┴─────┘

        >>> lf.join(other_lf, how="cross").collect()
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
        """
        require_same_type(self, other)

        if maintain_order is None:
            maintain_order = "none"

        uses_on = on is not None
        uses_left_on = left_on is not None
        uses_right_on = right_on is not None
        uses_lr_on = uses_left_on or uses_right_on
        if uses_on and uses_lr_on:
            msg = "cannot use 'on' in conjunction with 'left_on' or 'right_on'"
            raise ValueError(msg)
        elif uses_left_on != uses_right_on:
            msg = "'left_on' requires corresponding 'right_on'"
            raise ValueError(msg)

        if how == "outer":
            how = "full"
            issue_deprecation_warning(
                "use of `how='outer'` should be replaced with `how='full'`.",
                version="0.20.29",
            )
        elif how == "outer_coalesce":  # type: ignore[comparison-overlap]
            coalesce = True
            how = "full"
            issue_deprecation_warning(
                "use of `how='outer_coalesce'` should be replaced with `how='full', coalesce=True`.",
                version="0.20.29",
            )
        elif how == "cross":
            if uses_on or uses_lr_on:
                msg = "cross join should not pass join keys"
                raise ValueError(msg)
            return self._from_pyldf(
                self._ldf.join(
                    other._ldf,
                    [],
                    [],
                    allow_parallel,
                    force_parallel,
                    nulls_equal,
                    how,
                    suffix,
                    validate,
                    maintain_order,
                    coalesce=None,
                )
            )

        if uses_on:
            pyexprs = parse_into_list_of_expressions(on)
            pyexprs_left = pyexprs
            pyexprs_right = pyexprs
        elif uses_lr_on:
            pyexprs_left = parse_into_list_of_expressions(left_on)
            pyexprs_right = parse_into_list_of_expressions(right_on)
        else:
            msg = "must specify `on` OR `left_on` and `right_on`"
            raise ValueError(msg)

        return self._from_pyldf(
            self._ldf.join(
                other._ldf,
                pyexprs_left,
                pyexprs_right,
                allow_parallel,
                force_parallel,
                nulls_equal,
                how,
                suffix,
                validate,
                maintain_order,
                coalesce,
            )
        )

    @unstable()
    def join_where(
        self,
        other: LazyFrame,
        *predicates: Expr | Iterable[Expr],
        suffix: str = "_right",
    ) -> LazyFrame:
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
        Join two lazyframes together based on two predicates which get AND-ed together.

        >>> east = pl.LazyFrame(
        ...     {
        ...         "id": [100, 101, 102],
        ...         "dur": [120, 140, 160],
        ...         "rev": [12, 14, 16],
        ...         "cores": [2, 8, 4],
        ...     }
        ... )
        >>> west = pl.LazyFrame(
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
        ... ).collect()
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
        ... ).collect()
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

        pyexprs = parse_into_list_of_expressions(*predicates)

        return self._from_pyldf(
            self._ldf.join_where(
                other._ldf,
                pyexprs,
                suffix,
            )
        )

    def with_columns(
        self,
        *exprs: IntoExpr | Iterable[IntoExpr],
        **named_exprs: IntoExpr,
    ) -> LazyFrame:
        """
        Add columns to this LazyFrame.

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
        LazyFrame
            A new LazyFrame with the columns added.

        Notes
        -----
        Creating a new LazyFrame using this method does not create a new copy of
        existing data.

        Examples
        --------
        Pass an expression to add it as a new column.

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [0.5, 4, 10, 13],
        ...         "c": [True, True, False, True],
        ...     }
        ... )
        >>> lf.with_columns((pl.col("a") ** 2).alias("a^2")).collect()
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

        >>> lf.with_columns(pl.col("a").cast(pl.Float64)).collect()
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

        >>> lf.with_columns(
        ...     (pl.col("a") ** 2).alias("a^2"),
        ...     (pl.col("b") / 2).alias("b/2"),
        ...     (pl.col("c").not_()).alias("not c"),
        ... ).collect()
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

        >>> lf.with_columns(
        ...     [
        ...         (pl.col("a") ** 2).alias("a^2"),
        ...         (pl.col("b") / 2).alias("b/2"),
        ...         (pl.col("c").not_()).alias("not c"),
        ...     ]
        ... ).collect()
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

        >>> lf.with_columns(
        ...     ab=pl.col("a") * pl.col("b"),
        ...     not_c=pl.col("c").not_(),
        ... ).collect()
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
        structify = bool(int(os.environ.get("POLARS_AUTO_STRUCTIFY", 0)))

        pyexprs = parse_into_list_of_expressions(
            *exprs, **named_exprs, __structify=structify
        )
        return self._from_pyldf(self._ldf.with_columns(pyexprs))

    def with_columns_seq(
        self,
        *exprs: IntoExpr | Iterable[IntoExpr],
        **named_exprs: IntoExpr,
    ) -> LazyFrame:
        """
        Add columns to this LazyFrame.

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
        LazyFrame
            A new LazyFrame with the columns added.

        See Also
        --------
        with_columns
        """
        structify = bool(int(os.environ.get("POLARS_AUTO_STRUCTIFY", 0)))

        pyexprs = parse_into_list_of_expressions(
            *exprs, **named_exprs, __structify=structify
        )
        return self._from_pyldf(self._ldf.with_columns_seq(pyexprs))

    @deprecated(
        "`LazyFrame.with_context` is deprecated; "
        "use `pl.concat(..., how='horizontal')` instead."
    )
    def with_context(self, other: Self | list[Self]) -> LazyFrame:
        """
        Add an external context to the computation graph.

        .. deprecated:: 1.0.0
            Use :func:`concat` instead, with `how='horizontal'`

        This allows expressions to also access columns from DataFrames
        that are not part of this one.

        Parameters
        ----------
        other
            Lazy DataFrame to join with.

        Examples
        --------
        >>> lf = pl.LazyFrame({"a": [1, 2, 3], "b": ["a", "c", None]})
        >>> lf_other = pl.LazyFrame({"c": ["foo", "ham"]})
        >>> lf.with_context(lf_other).select(  # doctest: +SKIP
        ...     pl.col("b") + pl.col("c").first()
        ... ).collect()
        shape: (3, 1)
        ┌──────┐
        │ b    │
        │ ---  │
        │ str  │
        ╞══════╡
        │ afoo │
        │ cfoo │
        │ null │
        └──────┘

        Fill nulls with the median from another DataFrame:

        >>> train_lf = pl.LazyFrame(
        ...     {"feature_0": [-1.0, 0, 1], "feature_1": [-1.0, 0, 1]}
        ... )
        >>> test_lf = pl.LazyFrame(
        ...     {"feature_0": [-1.0, None, 1], "feature_1": [-1.0, 0, 1]}
        ... )
        >>> test_lf.with_context(  # doctest: +SKIP
        ...     train_lf.select(pl.all().name.suffix("_train"))
        ... ).select(
        ...     pl.col("feature_0").fill_null(pl.col("feature_0_train").median())
        ... ).collect()
        shape: (3, 1)
        ┌───────────┐
        │ feature_0 │
        │ ---       │
        │ f64       │
        ╞═══════════╡
        │ -1.0      │
        │ 0.0       │
        │ 1.0       │
        └───────────┘
        """
        if not isinstance(other, list):
            other = [other]

        return self._from_pyldf(self._ldf.with_context([lf._ldf for lf in other]))

    def drop(
        self,
        *columns: ColumnNameOrSelector | Iterable[ColumnNameOrSelector],
        strict: bool = True,
    ) -> LazyFrame:
        """
        Remove columns from the DataFrame.

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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6.0, 7.0, 8.0],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.drop("ham").collect()
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

        Drop multiple columns by passing a selector.

        >>> import polars.selectors as cs
        >>> lf.drop(cs.numeric()).collect()
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

        >>> lf.drop("foo", "ham").collect()
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
        selectors: list[ColumnNameOrSelector] = []
        for c in columns:
            if isinstance(c, Iterable) and not isinstance(c, str):
                selectors += c
            else:
                selectors += [c]

        drop_cols = parse_list_into_selector(selectors, strict=strict)
        return self._from_pyldf(self._ldf.drop(columns=drop_cols._pyselector))

    def rename(
        self, mapping: Mapping[str, str] | Callable[[str], str], *, strict: bool = True
    ) -> LazyFrame:
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

        Notes
        -----
        If existing names are swapped (e.g. 'A' points to 'B' and 'B' points to 'A'),
        polars will block projection and predicate pushdowns at this node.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, 7, 8],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.rename({"foo": "apple"}).collect()
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
        >>> lf.rename(lambda column_name: "c" + column_name[1:]).collect()
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
        if callable(mapping):
            return self.select(F.all().name.map(mapping))
        else:
            existing = list(mapping.keys())
            new = list(mapping.values())
            return self._from_pyldf(self._ldf.rename(existing, new, strict))

    def reverse(self) -> LazyFrame:
        """
        Reverse the DataFrame.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "key": ["a", "b", "c"],
        ...         "val": [1, 2, 3],
        ...     }
        ... )
        >>> lf.reverse().collect()
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
        return self._from_pyldf(self._ldf.reverse())

    def shift(
        self, n: int | IntoExprColumn = 1, *, fill_value: IntoExpr | None = None
    ) -> LazyFrame:
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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [5, 6, 7, 8],
        ...     }
        ... )
        >>> lf.shift().collect()
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

        >>> lf.shift(-2).collect()
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

        >>> lf.shift(-2, fill_value=100).collect()
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
        if fill_value is not None:
            fill_value_py = parse_into_expression(fill_value, str_as_lit=True)
        else:
            fill_value_py = None
        n_py = parse_into_expression(n)
        return self._from_pyldf(self._ldf.shift(n_py, fill_value_py))

    def slice(self, offset: int, length: int | None = None) -> LazyFrame:
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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["x", "y", "z"],
        ...         "b": [1, 3, 5],
        ...         "c": [2, 4, 6],
        ...     }
        ... )
        >>> lf.slice(1, 2).collect()
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ str ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ y   ┆ 3   ┆ 4   │
        │ z   ┆ 5   ┆ 6   │
        └─────┴─────┴─────┘
        """
        if length and length < 0:
            msg = f"negative slice lengths ({length!r}) are invalid for LazyFrame"
            raise ValueError(msg)
        return self._from_pyldf(self._ldf.slice(offset, length))

    def limit(self, n: int = 5) -> LazyFrame:
        """
        Get the first `n` rows.

        Alias for :func:`LazyFrame.head`.

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5, 6],
        ...         "b": [7, 8, 9, 10, 11, 12],
        ...     }
        ... )
        >>> lf.limit().collect()
        shape: (5, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 7   │
        │ 2   ┆ 8   │
        │ 3   ┆ 9   │
        │ 4   ┆ 10  │
        │ 5   ┆ 11  │
        └─────┴─────┘
        >>> lf.limit(2).collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 7   │
        │ 2   ┆ 8   │
        └─────┴─────┘
        """
        return self.head(n)

    def head(self, n: int = 5) -> LazyFrame:
        """
        Get the first `n` rows.

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5, 6],
        ...         "b": [7, 8, 9, 10, 11, 12],
        ...     }
        ... )
        >>> lf.head().collect()
        shape: (5, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 7   │
        │ 2   ┆ 8   │
        │ 3   ┆ 9   │
        │ 4   ┆ 10  │
        │ 5   ┆ 11  │
        └─────┴─────┘
        >>> lf.head(2).collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 7   │
        │ 2   ┆ 8   │
        └─────┴─────┘
        """
        return self.slice(0, n)

    def tail(self, n: int = 5) -> LazyFrame:
        """
        Get the last `n` rows.

        Parameters
        ----------
        n
            Number of rows to return.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4, 5, 6],
        ...         "b": [7, 8, 9, 10, 11, 12],
        ...     }
        ... )
        >>> lf.tail().collect()
        shape: (5, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 2   ┆ 8   │
        │ 3   ┆ 9   │
        │ 4   ┆ 10  │
        │ 5   ┆ 11  │
        │ 6   ┆ 12  │
        └─────┴─────┘
        >>> lf.tail(2).collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 5   ┆ 11  │
        │ 6   ┆ 12  │
        └─────┴─────┘
        """
        return self._from_pyldf(self._ldf.tail(n))

    def last(self) -> LazyFrame:
        """
        Get the last row of the DataFrame.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 5, 3],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> lf.last().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 3   ┆ 6   │
        └─────┴─────┘
        """
        return self.tail(1)

    def first(self) -> LazyFrame:
        """
        Get the first row of the DataFrame.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> lf.first().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 2   │
        └─────┴─────┘
        """
        return self.slice(0, 1)

    @deprecated(
        "`LazyFrame.approx_n_unique` is deprecated; "
        "use `select(pl.all().approx_n_unique())` instead."
    )
    def approx_n_unique(self) -> LazyFrame:
        """
        Approximate count of unique values.

        .. deprecated:: 0.20.11
            Use `select(pl.all().approx_n_unique())` instead.

        This is done using the HyperLogLog++ algorithm for cardinality estimation.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.approx_n_unique().collect()  # doctest: +SKIP
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ u32 ┆ u32 │
        ╞═════╪═════╡
        │ 4   ┆ 2   │
        └─────┴─────┘
        """
        return self.select(F.all().approx_n_unique())

    def with_row_index(self, name: str = "index", offset: int = 0) -> LazyFrame:
        """
        Add a row index as the first column in the LazyFrame.

        Parameters
        ----------
        name
            Name of the index column.
        offset
            Start the index at this offset. Cannot be negative.

        Warnings
        --------
        Using this function can have a negative effect on query performance.
        This may, for instance, block predicate pushdown optimization.

        Notes
        -----
        The resulting column does not have any special properties. It is a regular
        column of type `UInt32` (or `UInt64` in `polars[rt64]`).

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> lf.with_row_index().collect()
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
        >>> lf.with_row_index("id", offset=1000).collect()
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

        >>> lf.select(
        ...     pl.int_range(pl.len(), dtype=pl.UInt32).alias("index"),
        ...     pl.all(),
        ... ).collect()
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
            return self._from_pyldf(self._ldf.with_row_index(name, offset))
        except OverflowError:
            issue = "negative" if offset < 0 else "greater than the maximum index value"
            msg = f"`offset` input for `with_row_index` cannot be {issue}, got {offset}"
            raise ValueError(msg) from None

    @deprecated(
        "`LazyFrame.with_row_count` is deprecated; use `LazyFrame.with_row_index` instead."
        " Note that the default column name has changed from 'row_nr' to 'index'."
    )
    def with_row_count(self, name: str = "row_nr", offset: int = 0) -> LazyFrame:
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
            Start the row count at this offset.

        Warnings
        --------
        This can have a negative effect on query performance.
        This may, for instance, block predicate pushdown optimization.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 3, 5],
        ...         "b": [2, 4, 6],
        ...     }
        ... )
        >>> lf.with_row_count().collect()  # doctest: +SKIP
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

    def gather_every(self, n: int, offset: int = 0) -> LazyFrame:
        """
        Take every nth row in the LazyFrame and return as a new LazyFrame.

        Parameters
        ----------
        n
            Gather every *n*-th row.
        offset
            Starting index.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [5, 6, 7, 8],
        ...     }
        ... )
        >>> lf.gather_every(2).collect()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 5   │
        │ 3   ┆ 7   │
        └─────┴─────┘
        >>> lf.gather_every(2, offset=1).collect()
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

    def fill_null(
        self,
        value: Any | Expr | None = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
        *,
        matches_supertype: bool = True,
    ) -> LazyFrame:
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
            Fill all matching supertypes of the fill `value` literal.

        See Also
        --------
        fill_nan

        Notes
        -----
        A null value is not the same as a NaN value.
        To fill NaN values, use :func:`fill_nan`.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, None, 4],
        ...         "b": [0.5, 4, None, 13],
        ...     }
        ... )
        >>> lf.fill_null(99).collect()
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
        >>> lf.fill_null(strategy="forward").collect()
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

        >>> lf.fill_null(strategy="max").collect()
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

        >>> lf.fill_null(strategy="zero").collect()
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
        from polars import Decimal

        dtypes: Sequence[PolarsDataType] | None

        if value is not None:
            if isinstance(value, pl.Expr):
                dtypes = None
            elif isinstance(value, bool):
                dtypes = [Boolean]
            elif matches_supertype and isinstance(value, (int, float)):
                dtypes = [
                    Int8,
                    Int16,
                    Int32,
                    Int64,
                    Int128,
                    UInt8,
                    UInt16,
                    UInt32,
                    UInt64,
                    Float32,
                    Float64,
                    Decimal,
                ]
            elif isinstance(value, int):
                dtypes = [Int64]
            elif isinstance(value, float):
                dtypes = [Float64]
            elif isinstance(value, datetime):
                dtypes = [Datetime] + [Datetime(u) for u in DTYPE_TEMPORAL_UNITS]
            elif isinstance(value, timedelta):
                dtypes = [Duration] + [Duration(u) for u in DTYPE_TEMPORAL_UNITS]
            elif isinstance(value, date):
                dtypes = [Date]
            elif isinstance(value, time):
                dtypes = [Time]
            elif isinstance(value, str):
                dtypes = [String, Categorical]
            else:
                # fallback; anything not explicitly handled above
                dtypes = None

            if dtypes:
                return self.with_columns(
                    F.col(dtypes).fill_null(value, strategy, limit)
                )

        return self.select(F.all().fill_null(value, strategy, limit))

    def fill_nan(self, value: int | float | Expr | None) -> LazyFrame:
        """
        Fill floating point NaN values.

        Parameters
        ----------
        value
            Value used to fill NaN values.

        See Also
        --------
        fill_null

        Notes
        -----
        A NaN value is not the same as a null value.
        To fill null values, use :func:`fill_null`.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1.5, 2, float("nan"), 4],
        ...         "b": [0.5, 4, float("nan"), 13],
        ...     }
        ... )
        >>> lf.fill_nan(99).collect()
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
        if not isinstance(value, pl.Expr):
            value = F.lit(value)
        return self._from_pyldf(self._ldf.fill_nan(value._pyexpr))

    def std(self, ddof: int = 1) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their standard deviation value.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.std().collect()
        shape: (1, 2)
        ┌──────────┬─────┐
        │ a        ┆ b   │
        │ ---      ┆ --- │
        │ f64      ┆ f64 │
        ╞══════════╪═════╡
        │ 1.290994 ┆ 0.5 │
        └──────────┴─────┘
        >>> lf.std(ddof=0).collect()
        shape: (1, 2)
        ┌──────────┬──────────┐
        │ a        ┆ b        │
        │ ---      ┆ ---      │
        │ f64      ┆ f64      │
        ╞══════════╪══════════╡
        │ 1.118034 ┆ 0.433013 │
        └──────────┴──────────┘
        """
        return self._from_pyldf(self._ldf.std(ddof))

    def var(self, ddof: int = 1) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their variance value.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.var().collect()
        shape: (1, 2)
        ┌──────────┬──────┐
        │ a        ┆ b    │
        │ ---      ┆ ---  │
        │ f64      ┆ f64  │
        ╞══════════╪══════╡
        │ 1.666667 ┆ 0.25 │
        └──────────┴──────┘
        >>> lf.var(ddof=0).collect()
        shape: (1, 2)
        ┌──────┬────────┐
        │ a    ┆ b      │
        │ ---  ┆ ---    │
        │ f64  ┆ f64    │
        ╞══════╪════════╡
        │ 1.25 ┆ 0.1875 │
        └──────┴────────┘
        """
        return self._from_pyldf(self._ldf.var(ddof))

    def max(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their maximum value.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.max().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 4   ┆ 2   │
        └─────┴─────┘
        """
        return self._from_pyldf(self._ldf.max())

    def min(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their minimum value.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.min().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 1   │
        └─────┴─────┘
        """
        return self._from_pyldf(self._ldf.min())

    def sum(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their sum value.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.sum().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 10  ┆ 5   │
        └─────┴─────┘
        """
        return self._from_pyldf(self._ldf.sum())

    def mean(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their mean value.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.mean().collect()
        shape: (1, 2)
        ┌─────┬──────┐
        │ a   ┆ b    │
        │ --- ┆ ---  │
        │ f64 ┆ f64  │
        ╞═════╪══════╡
        │ 2.5 ┆ 1.25 │
        └─────┴──────┘
        """
        return self._from_pyldf(self._ldf.mean())

    def median(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their median value.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.median().collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ f64 ┆ f64 │
        ╞═════╪═════╡
        │ 2.5 ┆ 1.0 │
        └─────┴─────┘
        """
        return self._from_pyldf(self._ldf.median())

    def null_count(self) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame as the sum of their null value count.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, None, 3],
        ...         "bar": [6, 7, None],
        ...         "ham": ["a", "b", "c"],
        ...     }
        ... )
        >>> lf.null_count().collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ u32 ┆ u32 ┆ u32 │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 1   ┆ 0   │
        └─────┴─────┴─────┘
        """
        return self._from_pyldf(self._ldf.null_count())

    def quantile(
        self,
        quantile: float | Expr,
        interpolation: QuantileMethod = "nearest",
    ) -> LazyFrame:
        """
        Aggregate the columns in the LazyFrame to their quantile value.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [1, 2, 3, 4],
        ...         "b": [1, 2, 1, 1],
        ...     }
        ... )
        >>> lf.quantile(0.7).collect()
        shape: (1, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ f64 ┆ f64 │
        ╞═════╪═════╡
        │ 3.0 ┆ 1.0 │
        └─────┴─────┘
        """  # noqa: W505
        quantile_py = parse_into_expression(quantile)
        return self._from_pyldf(self._ldf.quantile(quantile_py, interpolation))

    def explode(
        self,
        columns: ColumnNameOrSelector | Iterable[ColumnNameOrSelector],
        *more_columns: ColumnNameOrSelector,
    ) -> LazyFrame:
        """
        Explode the DataFrame to long format by exploding the given columns.

        Parameters
        ----------
        columns
            Column names, expressions, or a selector defining them. The underlying
            columns being exploded must be of the `List` or `Array` data type.
        *more_columns
            Additional names of columns to explode, specified as positional arguments.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "letters": ["a", "a", "b", "c"],
        ...         "numbers": [[1], [2, 3], [4, 5], [6, 7, 8]],
        ...     }
        ... )
        >>> lf.explode("numbers").collect()
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
        subset = parse_list_into_selector(columns) | parse_list_into_selector(  # type: ignore[arg-type]
            more_columns
        )
        return self._from_pyldf(self._ldf.explode(subset=subset._pyselector))

    def unique(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
        *,
        keep: UniqueKeepStrategy = "any",
        maintain_order: bool = False,
    ) -> LazyFrame:
        """
        Drop duplicate rows from this DataFrame.

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
        LazyFrame
            LazyFrame with unique rows.

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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3, 1],
        ...         "bar": ["a", "a", "a", "a"],
        ...         "ham": ["b", "b", "b", "b"],
        ...     }
        ... )
        >>> lf.unique(maintain_order=True).collect()
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
        >>> lf.unique(subset=["bar", "ham"], maintain_order=True).collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ foo ┆ bar ┆ ham │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ str ┆ str │
        ╞═════╪═════╪═════╡
        │ 1   ┆ a   ┆ b   │
        └─────┴─────┴─────┘
        >>> lf.unique(keep="last", maintain_order=True).collect()
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
        selector_subset: PySelector | None = None
        if subset is not None:
            selector_subset = parse_list_into_selector(subset)._pyselector
        return self._from_pyldf(self._ldf.unique(maintain_order, selector_subset, keep))

    def drop_nans(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
    ) -> LazyFrame:
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
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [-20.5, float("nan"), 80.0],
        ...         "bar": [float("nan"), 110.0, 25.5],
        ...         "ham": ["xxx", "yyy", None],
        ...     }
        ... )

        The default behavior of this method is to drop rows where any single
        value in the row is NaN:

        >>> lf.drop_nans().collect()
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

        >>> lf.drop_nans(subset=["bar"]).collect()
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

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [float("nan"), float("nan"), float("nan"), float("nan")],
        ...         "b": [10.0, 2.5, float("nan"), 5.25],
        ...         "c": [65.75, float("nan"), float("nan"), 10.5],
        ...     }
        ... )
        >>> lf.filter(~pl.all_horizontal(pl.all().is_nan())).collect()
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
        selector_subset: PySelector | None = None
        if subset is not None:
            selector_subset = parse_list_into_selector(subset)._pyselector
        return self._from_pyldf(self._ldf.drop_nans(subset=selector_subset))

    def drop_nulls(
        self,
        subset: ColumnNameOrSelector | Collection[ColumnNameOrSelector] | None = None,
    ) -> LazyFrame:
        """
        Drop all rows that contain one or more null values.

        The original order of the remaining rows is preserved.

        See Also
        --------
        drop_nans

        Notes
        -----
        A null value is not the same as a NaN value.
        To drop NaN values, use :func:`drop_nans`.


        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, 2, 3],
        ...         "bar": [6, None, 8],
        ...         "ham": ["a", "b", None],
        ...     }
        ... )

        The default behavior of this method is to drop rows where any single
        value in the row is null:

        >>> lf.drop_nulls().collect()
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
        >>> lf.drop_nulls(subset=cs.integer()).collect()
        shape: (2, 3)
        ┌─────┬─────┬──────┐
        │ foo ┆ bar ┆ ham  │
        │ --- ┆ --- ┆ ---  │
        │ i64 ┆ i64 ┆ str  │
        ╞═════╪═════╪══════╡
        │ 1   ┆ 6   ┆ a    │
        │ 3   ┆ 8   ┆ null │
        └─────┴─────┴──────┘

        Dropping a row only if *all* values are null requires a different formulation:

        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": [None, None, None, None],
        ...         "b": [1, 2, None, 1],
        ...         "c": [1, None, None, 1],
        ...     }
        ... )
        >>> lf.filter(~pl.all_horizontal(pl.all().is_null())).collect()
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
        """
        selector_subset: PySelector | None = None
        if subset is not None:
            selector_subset = parse_list_into_selector(subset)._pyselector
        return self._from_pyldf(self._ldf.drop_nulls(subset=selector_subset))

    def unpivot(
        self,
        on: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        *,
        index: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        variable_name: str | None = None,
        value_name: str | None = None,
        streamable: bool = True,
    ) -> LazyFrame:
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
        streamable
            deprecated

        Notes
        -----
        If you're coming from pandas, this is similar to `pandas.DataFrame.melt`,
        but with `index` replacing `id_vars` and `on` replacing `value_vars`.
        In other frameworks, you might know this operation as `pivot_longer`.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "a": ["x", "y", "z"],
        ...         "b": [1, 3, 5],
        ...         "c": [2, 4, 6],
        ...     }
        ... )
        >>> import polars.selectors as cs
        >>> lf.unpivot(cs.numeric(), index="a").collect()
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
        if not streamable:
            issue_deprecation_warning(
                "the `streamable` parameter for `LazyFrame.unpivot` is deprecated"
                "This parameter has no effect",
                version="1.5.0",
            )

        selector_on: pl.Selector = (
            cs.empty() if on is None else parse_list_into_selector(on)
        )
        selector_index: pl.Selector = (
            cs.empty() if index is None else parse_list_into_selector(index)
        )

        return self._from_pyldf(
            self._ldf.unpivot(
                selector_on._pyselector,
                selector_index._pyselector,
                value_name,
                variable_name,
            )
        )

    def map_batches(
        self,
        function: Callable[[DataFrame], DataFrame],
        *,
        predicate_pushdown: bool = True,
        projection_pushdown: bool = True,
        slice_pushdown: bool = True,
        no_optimizations: bool = False,
        schema: None | SchemaDict = None,
        validate_output_schema: bool = True,
        streamable: bool = False,
    ) -> LazyFrame:
        """
        Apply a custom function.

        It is important that the function returns a Polars DataFrame.

        Parameters
        ----------
        function
            Lambda/ function to apply.
        predicate_pushdown
            Allow predicate pushdown optimization to pass this node.
        projection_pushdown
            Allow projection pushdown optimization to pass this node.
        slice_pushdown
            Allow slice pushdown optimization to pass this node.
        no_optimizations
            Turn off all optimizations past this point.
        schema
            Output schema of the function, if set to `None` we assume that the schema
            will remain unchanged by the applied function.
        validate_output_schema
            It is paramount that polars' schema is correct. This flag will ensure that
            the output schema of this function will be checked with the expected schema.
            Setting this to `False` will not do this check, but may lead to hard to
            debug bugs.
        streamable
            Whether the function that is given is eligible to be running with the
            streaming engine. That means that the function must produce the same result
            when it is executed in batches or when it is be executed on the full
            dataset.

        Warnings
        --------
        The `schema` of a `LazyFrame` must always be correct. It is up to the caller
        of this function to ensure that this invariant is upheld.

        It is important that the optimization flags are correct. If the custom function
        for instance does an aggregation of a column, `predicate_pushdown` should not
        be allowed, as this prunes rows and will influence your aggregation results.

        Notes
        -----
        A UDF passed to `map_batches` must be pure, meaning that it cannot modify or
        depend on state other than its arguments.

        Examples
        --------
        >>> lf = (  # doctest: +SKIP
        ...     pl.LazyFrame(
        ...         {
        ...             "a": pl.int_range(-100_000, 0, eager=True),
        ...             "b": pl.int_range(0, 100_000, eager=True),
        ...         }
        ...     )
        ...     .map_batches(lambda x: 2 * x, streamable=True)
        ...     .collect(engine="streaming")
        ... )
        shape: (100_000, 2)
        ┌─────────┬────────┐
        │ a       ┆ b      │
        │ ---     ┆ ---    │
        │ i64     ┆ i64    │
        ╞═════════╪════════╡
        │ -200000 ┆ 0      │
        │ -199998 ┆ 2      │
        │ -199996 ┆ 4      │
        │ -199994 ┆ 6      │
        │ …       ┆ …      │
        │ -8      ┆ 199992 │
        │ -6      ┆ 199994 │
        │ -4      ┆ 199996 │
        │ -2      ┆ 199998 │
        └─────────┴────────┘
        """
        if no_optimizations:
            predicate_pushdown = False
            projection_pushdown = False
            slice_pushdown = False

        return self._from_pyldf(
            self._ldf.map_batches(
                function,
                predicate_pushdown,
                projection_pushdown,
                slice_pushdown,
                streamable=streamable,
                schema=schema,
                validate_output=validate_output_schema,
            )
        )

    def interpolate(self) -> LazyFrame:
        """
        Interpolate intermediate values. The interpolation method is linear.

        Nulls at the beginning and end of the series remain null.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "foo": [1, None, 9, 10],
        ...         "bar": [6, 7, 9, None],
        ...         "baz": [1, None, None, 9],
        ...     }
        ... )
        >>> lf.interpolate().collect()
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

    def unnest(
        self,
        columns: ColumnNameOrSelector | Collection[ColumnNameOrSelector],
        *more_columns: ColumnNameOrSelector,
        separator: str | None = None,
    ) -> LazyFrame:
        """
        Decompose struct columns into separate columns for each of their fields.

        The new columns will be inserted into the DataFrame at the location of the
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
        >>> df = pl.LazyFrame(
        ...     {
        ...         "before": ["foo", "bar"],
        ...         "t_a": [1, 2],
        ...         "t_b": ["a", "b"],
        ...         "t_c": [True, None],
        ...         "t_d": [[1, 2], [3]],
        ...         "after": ["baz", "womp"],
        ...     }
        ... ).select("before", pl.struct(pl.col("^t_.$")).alias("t_struct"), "after")
        >>> df.collect()
        shape: (2, 3)
        ┌────────┬─────────────────────┬───────┐
        │ before ┆ t_struct            ┆ after │
        │ ---    ┆ ---                 ┆ ---   │
        │ str    ┆ struct[4]           ┆ str   │
        ╞════════╪═════════════════════╪═══════╡
        │ foo    ┆ {1,"a",true,[1, 2]} ┆ baz   │
        │ bar    ┆ {2,"b",null,[3]}    ┆ womp  │
        └────────┴─────────────────────┴───────┘
        >>> df.unnest("t_struct").collect()
        shape: (2, 6)
        ┌────────┬─────┬─────┬──────┬───────────┬───────┐
        │ before ┆ t_a ┆ t_b ┆ t_c  ┆ t_d       ┆ after │
        │ ---    ┆ --- ┆ --- ┆ ---  ┆ ---       ┆ ---   │
        │ str    ┆ i64 ┆ str ┆ bool ┆ list[i64] ┆ str   │
        ╞════════╪═════╪═════╪══════╪═══════════╪═══════╡
        │ foo    ┆ 1   ┆ a   ┆ true ┆ [1, 2]    ┆ baz   │
        │ bar    ┆ 2   ┆ b   ┆ null ┆ [3]       ┆ womp  │
        └────────┴─────┴─────┴──────┴───────────┴───────┘
        >>> df = pl.LazyFrame(
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
        >>> df.unnest("t", separator="::").collect()
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
        subset = parse_list_into_selector(columns) | parse_list_into_selector(
            more_columns
        )
        return self._from_pyldf(self._ldf.unnest(subset._pyselector, separator))

    def merge_sorted(self, other: LazyFrame, key: str) -> LazyFrame:
        """
        Take two sorted DataFrames and merge them by the sorted key.

        The output of this operation will also be sorted.
        It is the callers responsibility that the frames
        are sorted in ascending order by that key otherwise
        the output will not make sense.

        The schemas of both LazyFrames must be equal.

        Parameters
        ----------
        other
            Other DataFrame that must be merged
        key
            Key that is sorted.

        Examples
        --------
        >>> df0 = pl.LazyFrame(
        ...     {"name": ["steve", "elise", "bob"], "age": [42, 44, 18]}
        ... ).sort("age")
        >>> df0.collect()
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
        >>> df1 = pl.LazyFrame(
        ...     {"name": ["anna", "megan", "steve", "thomas"], "age": [21, 33, 42, 20]}
        ... ).sort("age")
        >>> df1.collect()
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
        >>> df0.merge_sorted(df1, key="age").collect()
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
        require_same_type(self, other)
        return self._from_pyldf(self._ldf.merge_sorted(other._ldf, key))

    def set_sorted(
        self,
        column: str,
        *more_columns: str,
        descending: bool | list[bool] = False,
        nulls_last: bool | list[bool] = False,
    ) -> LazyFrame:
        """
        Flag a column as sorted.

        This can speed up future operations.

        Parameters
        ----------
        column
            Column that is sorted
        more_columns
            Columns that are sorted over after `column`.
        descending
            Whether the column is sorted in descending order.
        nulls_last
            Whether the nulls are at the end.

        Warnings
        --------
        This can lead to incorrect results if the data is NOT sorted!!
        Use with care!

        """
        # NOTE: Only accepts 1 column on purpose! User think they are sorted by
        # the combined multicolumn values.
        if not isinstance(column, str):
            msg = "expected a 'str' for argument 'column' in 'set_sorted'"
            raise TypeError(msg)

        ds: list[bool]
        nl: list[bool]
        if isinstance(descending, bool):
            ds = [descending]
        else:
            ds = descending
        if isinstance(nulls_last, bool):
            nl = [nulls_last]
        else:
            nl = nulls_last

        return self._from_pyldf(
            self._ldf.hint_sorted(
                [column] + list(more_columns), descending=ds, nulls_last=nl
            )
        )

    @unstable()
    def update(
        self,
        other: LazyFrame,
        on: str | Sequence[str] | None = None,
        how: Literal["left", "inner", "full"] = "left",
        *,
        left_on: str | Sequence[str] | None = None,
        right_on: str | Sequence[str] | None = None,
        include_nulls: bool = False,
        maintain_order: MaintainOrderJoin | None = "left",
    ) -> LazyFrame:
        """
        Update the values in this `LazyFrame` with the values in `other`.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        other
            LazyFrame that will be used to update the values
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
            Which order of rows from the inputs to preserve. See :func:`~LazyFrame.join`
            for details. Unlike `join` this function preserves the left order by
            default.

        Notes
        -----
        This is syntactic sugar for a left/inner join that preserves the order
        of the left `DataFrame` by default, with an optional coalesce when
        `include_nulls = False`.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {
        ...         "A": [1, 2, 3, 4],
        ...         "B": [400, 500, 600, 700],
        ...     }
        ... )
        >>> lf.collect()
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
        >>> new_lf = pl.LazyFrame(
        ...     {
        ...         "B": [-66, None, -99],
        ...         "C": [5, 3, 1],
        ...     }
        ... )

        Update `df` values with the non-null values in `new_df`, by row index:

        >>> lf.update(new_lf).collect()
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

        >>> lf.update(new_lf, how="inner").collect()
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

        >>> lf.update(new_lf, left_on=["A"], right_on=["C"], how="full").collect()
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

        Update `df` values including null values in `new_df`, using a full
        outer join strategy that defines explicit join columns in each frame:

        >>> lf.update(
        ...     new_lf, left_on="A", right_on="C", how="full", include_nulls=True
        ... ).collect()
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
        require_same_type(self, other)
        if how in ("outer", "outer_coalesce"):
            how = "full"
            issue_deprecation_warning(
                "use of `how='outer'` should be replaced with `how='full'`.",
                version="0.20.29",
            )

        if how not in ("left", "inner", "full"):
            msg = f"`how` must be one of {{'left', 'inner', 'full'}}; found {how!r}"
            raise ValueError(msg)

        row_index_used = False
        if on is None:
            if left_on is None and right_on is None:
                # no keys provided--use row index
                row_index_used = True
                row_index_name = "__POLARS_ROW_INDEX"
                self = self.with_row_index(row_index_name)
                other = other.with_row_index(row_index_name)
                left_on = right_on = [row_index_name]
            else:
                # one of left or right is missing, raise error
                if left_on is None:
                    msg = "missing join columns for left frame"
                    raise ValueError(msg)
                if right_on is None:
                    msg = "missing join columns for right frame"
                    raise ValueError(msg)
        else:
            # move on into left/right_on to simplify logic
            left_on = right_on = on

        if isinstance(left_on, str):
            left_on = [left_on]
        if isinstance(right_on, str):
            right_on = [right_on]

        left_schema = self.collect_schema()
        for name in left_on:
            if name not in left_schema:
                msg = f"left join column {name!r} not found"
                raise ValueError(msg)
        right_schema = other.collect_schema()
        for name in right_on:
            if name not in right_schema:
                msg = f"right join column {name!r} not found"
                raise ValueError(msg)

        # no need to join if *only* join columns are in other (inner/left update only)
        if how != "full" and len(right_schema) == len(right_on):
            if row_index_used:
                return self.drop(row_index_name)
            return self

        # only use non-idx right columns present in left frame
        right_other = set(right_schema).intersection(left_schema) - set(right_on)

        # When include_nulls is True, we need to distinguish records after the join that
        # were originally null in the right frame, as opposed to records that were null
        # because the key was missing from the right frame.
        # Add a validity column to track whether row was matched or not.
        if include_nulls:
            validity = ("__POLARS_VALIDITY",)
            other = other.with_columns(F.lit(True).alias(validity[0]))
        else:
            validity = ()  # type: ignore[assignment]

        tmp_name = "__POLARS_RIGHT"
        drop_columns = [*(f"{name}{tmp_name}" for name in right_other), *validity]
        result = (
            self.join(
                other.select(*right_on, *right_other, *validity),
                left_on=left_on,
                right_on=right_on,
                how=how,
                suffix=tmp_name,
                coalesce=True,
                maintain_order=maintain_order,
            )
            .with_columns(
                (
                    # use left value only when right value failed to join
                    F.when(F.col(validity).is_null())
                    .then(F.col(name))
                    .otherwise(F.col(f"{name}{tmp_name}"))
                    if include_nulls
                    else F.coalesce([f"{name}{tmp_name}", F.col(name)])
                ).alias(name)
                for name in right_other
            )
            .drop(drop_columns)
        )
        if row_index_used:
            result = result.drop(row_index_name)

        return self._from_pyldf(result._ldf)

    def count(self) -> LazyFrame:
        """
        Return the number of non-null elements for each column.

        Examples
        --------
        >>> lf = pl.LazyFrame(
        ...     {"a": [1, 2, 3, 4], "b": [1, 2, 1, None], "c": [None, None, None, None]}
        ... )
        >>> lf.count().collect()
        shape: (1, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ u32 ┆ u32 ┆ u32 │
        ╞═════╪═════╪═════╡
        │ 4   ┆ 3   ┆ 0   │
        └─────┴─────┴─────┘
        """
        return self._from_pyldf(self._ldf.count())

    @deprecated(
        "`LazyFrame.melt` is deprecated; use `LazyFrame.unpivot` instead, with "
        "`index` instead of `id_vars` and `on` instead of `value_vars`"
    )
    def melt(
        self,
        id_vars: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        value_vars: ColumnNameOrSelector | Sequence[ColumnNameOrSelector] | None = None,
        variable_name: str | None = None,
        value_name: str | None = None,
        *,
        streamable: bool = True,
    ) -> LazyFrame:
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
        streamable
            Allow this node to run in the streaming engine.
            If this runs in streaming, the output of the unpivot operation
            will not have a stable ordering.
        """
        return self.unpivot(
            index=id_vars,
            on=value_vars,
            variable_name=variable_name,
            value_name=value_name,
            streamable=streamable,
        )

    @unstable()
    def remote(
        self,
        context: pc.ComputeContext | None = None,
        plan_type: pc._typing.PlanTypePreference = "dot",
    ) -> pc.LazyFrameRemote:
        """
        Run a query remotely on Polars Cloud.

        This allows you to run Polars remotely on
        one or more workers via several strategies
        for distributed compute.

        Read more in the `Announcement post <https://pola.rs/posts/polars-cloud-what-we-are-building/>`_

        Parameters
        ----------
        context
            Compute context in which queries are executed.
            If none given, it will take the default context.
        plan_type: {'plain', 'dot'}
            Whether to give a dot diagram of a plain text
            version of logical plan.

        Examples
        --------
        Run a query on a cloud instance.

        >>> lf = pl.LazyFrame([1, 2, 3]).sum()
        >>> in_progress = lf.remote().collect()  # doctest: +SKIP
        >>> # do some other work
        >>> in_progress.await_result()  # doctest: +SKIP
        shape: (1, 1)
        ┌──────────┐
        │ column_0 │
        │ ---      │
        │ i64      │
        ╞══════════╡
        │ 6        │
        └──────────┘

        Run a query distributed.

        >>> lf = (
        ...     pl.scan_parquet("s3://my_bucket/").group_by("key").agg(pl.sum("values"))
        ... )
        >>> in_progress = lf.remote().distributed().collect()  # doctest: +SKIP
        >>> in_progress.await_result()  # doctest: +SKIP
        shape: (1, 1)
        ┌──────────┐
        │ column_0 │
        │ ---      │
        │ i64      │
        ╞══════════╡
        │ 6        │
        └──────────┘

        """
        return pc.LazyFrameRemote(lf=self, context=context, plan_type=plan_type)

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
    ) -> LazyFrame:
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

        >>> lf = pl.LazyFrame({"a": [1, 2, 3], "b": ["A", "B", "C"]})
        >>> lf.match_to_schema({"a": pl.Int64, "b": pl.String}).collect()
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
        >>> (lf.match_to_schema({"a": pl.Int64}).collect())  # doctest: +SKIP
        polars.exceptions.SchemaError: extra columns in `match_to_schema`: "b"

        Adding missing columns

        >>> (
        ...     pl.LazyFrame({"a": [1, 2, 3]})
        ...     .match_to_schema(
        ...         {"a": pl.Int64, "b": pl.String},
        ...         missing_columns="insert",
        ...     )
        ...     .collect()
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
        ...     pl.LazyFrame({"a": [1, 2, 3]})
        ...     .match_to_schema(
        ...         {"a": pl.Int64, "b": pl.String},
        ...         missing_columns={"b": pl.col.a.cast(pl.String)},
        ...     )
        ...     .collect()
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
        ...     pl.LazyFrame({"a": [1, 2, 3], "b": ["A", "B", "C"]})
        ...     .match_to_schema(
        ...         {"a": pl.Int64},
        ...         extra_columns="ignore",
        ...     )
        ...     .collect()
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
        ...     pl.LazyFrame(
        ...         {"a": [1, 2, 3], "b": [1.0, 2.0, 3.0]},
        ...         schema={"a": pl.Int32, "b": pl.Float32},
        ...     )
        ...     .match_to_schema(
        ...         {"a": pl.Int64, "b": pl.Float64},
        ...         integer_cast="upcast",
        ...         float_cast="upcast",
        ...     )
        ...     .collect()
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
        from polars import Expr

        def prepare_missing_columns(
            value: Literal["insert", "raise"] | Expr,
        ) -> Literal["insert", "raise"] | PyExpr:
            if isinstance(value, Expr):
                return value._pyexpr
            return value

        schema_prep: Schema
        if isinstance(schema, Mapping):
            schema_prep = Schema(schema)
        else:
            schema_prep = schema

        missing_columns_pyexpr: (
            Literal["insert", "raise"] | dict[str, Literal["insert", "raise"] | PyExpr]
        )
        if isinstance(missing_columns, Mapping):
            missing_columns_pyexpr = {
                key: prepare_missing_columns(value)
                for key, value in missing_columns.items()
            }
        elif isinstance(missing_columns, Expr):
            missing_columns_pyexpr = prepare_missing_columns(missing_columns)
        else:
            missing_columns_pyexpr = missing_columns

        return LazyFrame._from_pyldf(
            self._ldf.match_to_schema(
                schema=schema_prep,
                missing_columns=missing_columns_pyexpr,
                missing_struct_fields=missing_struct_fields,
                extra_columns=extra_columns,
                extra_struct_fields=extra_struct_fields,
                integer_cast=integer_cast,
                float_cast=float_cast,
            )
        )

    def _to_metadata(
        self,
        columns: None | str | list[str] = None,
        stats: None | str | list[str] = None,
    ) -> DataFrame:
        """
        Get all runtime metadata for each column.

        This is unstable and is meant for debugging purposes.
        """
        lf = self

        if columns is not None:
            if isinstance(columns, str):
                columns = [columns]

            lf = lf.select(columns)

        return lf.collect()._to_metadata(stats=stats)
