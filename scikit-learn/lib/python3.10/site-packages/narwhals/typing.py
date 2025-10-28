from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Literal, Protocol, TypeVar, Union

from narwhals._compliant import CompliantDataFrame, CompliantLazyFrame, CompliantSeries
from narwhals._native import (
    IntoDataFrame,
    IntoDataFrameT,
    IntoFrame,
    IntoFrameT,
    IntoLazyFrame,
    IntoLazyFrameT,
    IntoSeries,
    IntoSeriesT,
)
from narwhals._typing import Backend, EagerAllowed, IntoBackend, LazyAllowed

if TYPE_CHECKING:
    import datetime as dt
    import os
    from collections.abc import Sequence
    from decimal import Decimal
    from types import ModuleType

    import numpy as np
    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import TypeAlias

    from narwhals import dtypes
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.expr import Expr
    from narwhals.schema import Schema
    from narwhals.series import Series

    class SupportsNativeNamespace(Protocol):
        def __native_namespace__(self) -> ModuleType: ...

    # ruff: noqa: N802
    class DTypes(Protocol):
        @property
        def Decimal(self) -> type[dtypes.Decimal]: ...
        @property
        def Int128(self) -> type[dtypes.Int128]: ...
        @property
        def Int64(self) -> type[dtypes.Int64]: ...
        @property
        def Int32(self) -> type[dtypes.Int32]: ...
        @property
        def Int16(self) -> type[dtypes.Int16]: ...
        @property
        def Int8(self) -> type[dtypes.Int8]: ...
        @property
        def UInt128(self) -> type[dtypes.UInt128]: ...
        @property
        def UInt64(self) -> type[dtypes.UInt64]: ...
        @property
        def UInt32(self) -> type[dtypes.UInt32]: ...
        @property
        def UInt16(self) -> type[dtypes.UInt16]: ...
        @property
        def UInt8(self) -> type[dtypes.UInt8]: ...
        @property
        def Float64(self) -> type[dtypes.Float64]: ...
        @property
        def Float32(self) -> type[dtypes.Float32]: ...
        @property
        def String(self) -> type[dtypes.String]: ...
        @property
        def Boolean(self) -> type[dtypes.Boolean]: ...
        @property
        def Object(self) -> type[dtypes.Object]: ...
        @property
        def Categorical(self) -> type[dtypes.Categorical]: ...
        @property
        def Enum(self) -> type[dtypes.Enum]: ...
        @property
        def Datetime(self) -> type[dtypes.Datetime]: ...
        @property
        def Duration(self) -> type[dtypes.Duration]: ...
        @property
        def Date(self) -> type[dtypes.Date]: ...
        @property
        def Field(self) -> type[dtypes.Field]: ...
        @property
        def Struct(self) -> type[dtypes.Struct]: ...
        @property
        def List(self) -> type[dtypes.List]: ...
        @property
        def Array(self) -> type[dtypes.Array]: ...
        @property
        def Unknown(self) -> type[dtypes.Unknown]: ...
        @property
        def Time(self) -> type[dtypes.Time]: ...
        @property
        def Binary(self) -> type[dtypes.Binary]: ...


IntoExpr: TypeAlias = Union["Expr", str, "Series[Any]"]
"""Anything which can be converted to an expression.

Use this to mean "either a Narwhals expression, or something which can be converted
into one". For example, `exprs` in `DataFrame.select` is typed to accept `IntoExpr`,
as it can either accept a `nw.Expr` (e.g. `df.select(nw.col('a'))`) or a string
which will be interpreted as a `nw.Expr`, e.g. `df.select('a')`.
"""

Frame: TypeAlias = Union["DataFrame[Any]", "LazyFrame[Any]"]
"""Narwhals DataFrame or Narwhals LazyFrame.

Use this if your function can work with either and your function doesn't care
about its backend.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import Frame
    >>> @nw.narwhalify
    ... def agnostic_columns(df: Frame) -> list[str]:
    ...     return df.columns
"""

FrameT = TypeVar("FrameT", "DataFrame[Any]", "LazyFrame[Any]")
"""TypeVar bound to Narwhals DataFrame or Narwhals LazyFrame.

Use this if your function accepts either `nw.DataFrame` or `nw.LazyFrame` and returns
an object of the same kind.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import FrameT
    >>> @nw.narwhalify
    ... def agnostic_func(df: FrameT) -> FrameT:
    ...     return df.with_columns(c=nw.col("a") + 1)
"""

DataFrameT = TypeVar("DataFrameT", bound="DataFrame[Any]")
"""TypeVar bound to Narwhals DataFrame.

Use this if your function can accept a Narwhals DataFrame and returns a Narwhals
DataFrame backed by the same backend.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import DataFrameT
    >>> @nw.narwhalify
    >>> def func(df: DataFrameT) -> DataFrameT:
    ...     return df.with_columns(c=df["a"] + 1)
"""

LazyFrameT = TypeVar("LazyFrameT", bound="LazyFrame[Any]")
SeriesT = TypeVar("SeriesT", bound="Series[Any]")
DTypeBackend: TypeAlias = 'Literal["pyarrow", "numpy_nullable"] | None'
SizeUnit: TypeAlias = Literal[
    "b",
    "kb",
    "mb",
    "gb",
    "tb",
    "bytes",
    "kilobytes",
    "megabytes",
    "gigabytes",
    "terabytes",
]

TimeUnit: TypeAlias = Literal["ns", "us", "ms", "s"]

AsofJoinStrategy: TypeAlias = Literal["backward", "forward", "nearest"]
"""Join strategy.

- *"backward"*: Selects the last row in the right DataFrame whose `on` key
    is less than or equal to the left's key.
- *"forward"*: Selects the first row in the right DataFrame whose `on` key
    is greater than or equal to the left's key.
- *"nearest"*: Search selects the last row in the right DataFrame whose value
    is nearest to the left's key.
"""

ClosedInterval: TypeAlias = Literal["left", "right", "none", "both"]
"""Define which sides of the interval are closed (inclusive)."""

ConcatMethod: TypeAlias = Literal["horizontal", "vertical", "diagonal"]
"""Concatenating strategy.

- *"vertical"*: Concatenate vertically. Column names must match.
- *"horizontal"*: Concatenate horizontally. If lengths don't match, then
    missing rows are filled with null values.
- *"diagonal"*: Finds a union between the column schemas and fills missing
    column values with null.
"""

FillNullStrategy: TypeAlias = Literal["forward", "backward"]
"""Strategy used to fill null values."""

JoinStrategy: TypeAlias = Literal["inner", "left", "full", "cross", "semi", "anti"]
"""Join strategy.

- *"inner"*: Returns rows that have matching values in both tables.
- *"left"*: Returns all rows from the left table, and the matched rows from
    the right table.
- *"full"*: Returns all rows in both dataframes, with the `suffix` appended to
    the right join keys.
- *"cross"*: Returns the Cartesian product of rows from both tables.
- *"semi"*: Filter rows that have a match in the right table.
- *"anti"*: Filter rows that do not have a match in the right table.
"""

PivotAgg: TypeAlias = Literal[
    "min", "max", "first", "last", "sum", "mean", "median", "len"
]
"""A predefined aggregate function string."""

RankMethod: TypeAlias = Literal["average", "min", "max", "dense", "ordinal"]
"""The method used to assign ranks to tied elements.

- *"average"*: The average of the ranks that would have been assigned to
    all the tied values is assigned to each value.
- *"min"*: The minimum of the ranks that would have been assigned to all
    the tied values is assigned to each value. (This is also referred to
    as "competition" ranking.)
- *"max"*: The maximum of the ranks that would have been assigned to all
    the tied values is assigned to each value.
- *"dense"*: Like "min", but the rank of the next highest element is
    assigned the rank immediately after those assigned to the tied elements.
- *"ordinal"*: All values are given a distinct rank, corresponding to the
    order that the values occur in the Series.
"""

RollingInterpolationMethod: TypeAlias = Literal[
    "nearest", "higher", "lower", "midpoint", "linear"
]
"""Interpolation method."""

UniqueKeepStrategy: TypeAlias = Literal["any", "first", "last", "none"]
"""Which of the duplicate rows to keep.

- *"any"*: Does not give any guarantee of which row is kept.
    This allows more optimizations.
- *"none"*: Don't keep duplicate rows.
- *"first"*: Keep first unique row.
- *"last"*: Keep last unique row.
"""

ModeKeepStrategy: TypeAlias = Literal["any", "all"]
"""Which of the mode's to keep.

- *"any"*: Does not give any guarantee of which mode is kept.
- *"all"*: Keeps all the mode's.
"""

_ShapeT = TypeVar("_ShapeT", bound="tuple[int, ...]")
_NDArray: TypeAlias = "np.ndarray[_ShapeT, Any]"
_1DArray: TypeAlias = "_NDArray[tuple[int]]"
_1DArrayInt: TypeAlias = "np.ndarray[tuple[int], np.dtype[np.integer[Any]]]"
_2DArray: TypeAlias = "_NDArray[tuple[int, int]]"  # noqa: PYI047
_AnyDArray: TypeAlias = "_NDArray[tuple[int, ...]]"  # noqa: PYI047
_NumpyScalar: TypeAlias = "np.generic[Any]"
Into1DArray: TypeAlias = "_1DArray | _NumpyScalar"
"""A 1-dimensional `numpy.ndarray` or scalar that can be converted into one."""

PandasLikeDType: TypeAlias = "pd.api.extensions.ExtensionDtype | np.dtype[Any]"


NumericLiteral: TypeAlias = "int | float | Decimal"
TemporalLiteral: TypeAlias = "dt.date | dt.datetime | dt.time | dt.timedelta"
NonNestedLiteral: TypeAlias = (
    "NumericLiteral | TemporalLiteral | str | bool | bytes | None"
)
PythonLiteral: TypeAlias = "NonNestedLiteral | list[Any] | tuple[Any, ...]"

NonNestedDType: TypeAlias = "dtypes.NumericType | dtypes.TemporalType | dtypes.String | dtypes.Boolean | dtypes.Binary | dtypes.Categorical | dtypes.Unknown | dtypes.Object"
"""Any Narwhals DType that does not have required arguments."""

IntoDType: TypeAlias = "dtypes.DType | type[NonNestedDType]"
"""Anything that can be converted into a Narwhals DType.

Examples:
    >>> import polars as pl
    >>> import narwhals as nw
    >>> df_native = pl.DataFrame({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})
    >>> df = nw.from_native(df_native)
    >>> df.select(
    ...     nw.col("a").cast(nw.Int32),
    ...     nw.col("b").cast(nw.String()).str.split(".").cast(nw.List(nw.Int8)),
    ... )
    ┌──────────────────┐
    |Narwhals DataFrame|
    |------------------|
    |shape: (3, 2)     |
    |┌─────┬──────────┐|
    |│ a   ┆ b        │|
    |│ --- ┆ ---      │|
    |│ i32 ┆ list[i8] │|
    |╞═════╪══════════╡|
    |│ 1   ┆ [4, 0]   │|
    |│ 2   ┆ [5, 0]   │|
    |│ 3   ┆ [6, 0]   │|
    |└─────┴──────────┘|
    └──────────────────┘
"""


# TODO @dangotbanned: fix this?
# Constructor allows tuples, but we don't support that *everywhere* yet
IntoSchema: TypeAlias = "Mapping[str, dtypes.DType] | Schema"
"""Anything that can be converted into a Narwhals Schema.

Defined by column names and their associated *instantiated* Narwhals DType.

Examples:
    >>> import narwhals as nw
    >>> import pyarrow as pa
    >>> data = {"a": [1, 2, 3], "b": [None, "hi", "howdy"], "c": [2.1, 2.0, None]}
    >>> nw.DataFrame.from_dict(
    ...     data,
    ...     schema={"a": nw.UInt8(), "b": nw.String(), "c": nw.Float32()},
    ...     backend="pyarrow",
    ... )
    ┌────────────────────────┐
    |   Narwhals DataFrame   |
    |------------------------|
    |pyarrow.Table           |
    |a: uint8                |
    |b: string               |
    |c: float                |
    |----                    |
    |a: [[1,2,3]]            |
    |b: [[null,"hi","howdy"]]|
    |c: [[2.1,2,null]]       |
    └────────────────────────┘
"""

IntoArrowSchema: TypeAlias = "pa.Schema | Mapping[str, pa.DataType]"
IntoPolarsSchema: TypeAlias = "pl.Schema | Mapping[str, pl.DataType]"
IntoPandasSchema: TypeAlias = Mapping[str, PandasLikeDType]

FileSource: TypeAlias = "str | os.PathLike[str]"
"""Path to a file.

Either a string or an object that implements [`__fspath__`], such as [`pathlib.Path`].

[`__fspath__`]: https://docs.python.org/3/library/os.html#os.PathLike
[`pathlib.Path`]: https://docs.python.org/3/library/pathlib.html#pathlib.Path
"""


# Annotations for `__getitem__` methods
_T = TypeVar("_T")
_Slice: TypeAlias = "slice[_T, Any, Any] | slice[Any, _T, Any] | slice[None, None, _T]"
_SliceNone: TypeAlias = "slice[None, None, None]"
# Index/column positions
SingleIndexSelector: TypeAlias = int
_SliceIndex: TypeAlias = "_Slice[int] | _SliceNone"
"""E.g. `[1:]` or `[:3]` or `[::2]`."""
SizedMultiIndexSelector: TypeAlias = "Sequence[int] | _T | _1DArrayInt"
MultiIndexSelector: TypeAlias = "_SliceIndex | SizedMultiIndexSelector[_T]"
# Labels/column names
SingleNameSelector: TypeAlias = str
_SliceName: TypeAlias = "_Slice[str] | _SliceNone"
SizedMultiNameSelector: TypeAlias = "Sequence[str] | _T | _1DArray"
MultiNameSelector: TypeAlias = "_SliceName | SizedMultiNameSelector[_T]"
# Mixed selectors
SingleColSelector: TypeAlias = "SingleIndexSelector | SingleNameSelector"
MultiColSelector: TypeAlias = "MultiIndexSelector[_T] | MultiNameSelector[_T]"


__all__ = [
    "Backend",
    "CompliantDataFrame",
    "CompliantLazyFrame",
    "CompliantSeries",
    "DataFrameT",
    "EagerAllowed",
    "Frame",
    "FrameT",
    "IntoBackend",
    "IntoDataFrame",
    "IntoDataFrameT",
    "IntoExpr",
    "IntoFrame",
    "IntoFrameT",
    "IntoLazyFrame",
    "IntoLazyFrameT",
    "IntoSeries",
    "IntoSeriesT",
    "LazyAllowed",
]
