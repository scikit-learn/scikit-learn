"""The home for *mostly* [structural] counterparts to [nominal] native types.

If you find yourself being yelled at by a typechecker and ended up here - **do not fear!**

We have 5 funky flavors, which tackle two different problem spaces.

How do we describe [Native types] when ...
- ... **wrapping in** a [Narwhals type]?
- ... **matching to** an [`Implementation`]?

## Wrapping in a Narwhals type
[//]: # (TODO @dangotbanned: Replace `Thing` with a better name)

The following examples use the placeholder type `Thing` which represents one of:
- `DataFrame`: (Eager) 2D data structure representing data as a table with rows and columns.
- `LazyFrame`: (Lazy) Computation graph/query against a DataFrame/database.
- `Series`: 1D data structure representing a single column.

Our goal is to **wrap** a *partially-unknown* native object **in** a [generic class]:

    def wrapping_in_df(native: IntoDataFrameT) -> DataFrame[IntoDataFrameT]: ...
    def wrapping_in_lf(native: IntoLazyFrameT) -> LazyFrame[IntoLazyFrameT]: ...
    def wrapping_in_ser(native: IntoSeriesT) -> Series[IntoSeriesT]: ...

### (1) `Native<Thing>`
Minimal [`Protocol`]s that are [assignable to] *almost any* supported native type of that group:

    class NativeThing(Protocol):
        def something_common(self, *args: Any, **kwargs: Any) -> Any: ...

Note:
    This group is primarily a building block for more useful types.

### (2) `Into<Thing>`
*Publicly* exported [`TypeAlias`]s of **(1)**:

    IntoThing: TypeAlias = NativeThing

**But**, occasionally, there'll be an edge-case which we can spell like:

    IntoThing: TypeAlias = Union[<type that does not fit the protocol>, NativeThing]

Tip:
    Reach for these when there **isn't a need to preserve** the original native type.

### (3) `Into<Thing>T`
*Publicly* exported [`TypeVar`]s, bound to **(2)**:

    IntoThingT = TypeVar("IntoThingT", bound=IntoThing)

Important:
    In most situations, you'll want to use these as they **do preserve** the original native type.

Putting it all together, we can now add a *narwhals-level* wrapper:

    class Thing(Generic[IntoThingT]):
        def to_native(self) -> IntoThingT: ...

## Matching to an `Implementation`
This problem differs as we need to *create* a relationship between *otherwise-unrelated* types.

Comparing the problems side-by-side, we can more clearly see this difference:

    def wrapping_in_df(native: IntoDataFrameT) -> DataFrame[IntoDataFrameT]: ...
    def matching_to_polars(native: pl.DataFrame) -> Literal[Implementation.POLARS]: ...

### (4) `Native<Backend>`
If we want to describe a set of specific types and **match** them in [`@overload`s], then these the tools we need.

For common and easily-installed backends, [`TypeAlias`]s are composed of the native type(s):

    NativePolars: TypeAlias = pl.DataFrame | pl.LazyFrame | pl.Series

Otherwise, we need to define a [`Protocol`] which the native type(s) can **match** against *when* installed:

    class NativeDask(NativeLazyFrame, Protocol):
        _partition_type: type[pd.DataFrame]

Tip:
    The goal is to be as minimal as possible, while still being *specific-enough* to **not match** something else.

Important:
    See [ibis#9276 comment] for a more *in-depth* example that doesn't fit here 😄

### (5) `is_native_<backend>`
[Type guards] for **(4)**, *similar* to those found in `nw.dependencies`.

They differ by checking **all** native types/protocols in a single-call and using ``Native<Backend>`` aliases.

[structural]: https://typing.python.org/en/latest/spec/glossary.html#term-structural
[nominal]: https://typing.python.org/en/latest/spec/glossary.html#term-nominal
[Native types]: https://narwhals-dev.github.io/narwhals/how_it_works/#polars-and-other-implementations
[Narwhals type]: https://narwhals-dev.github.io/narwhals/api-reference/dataframe/
[`Implementation`]: https://narwhals-dev.github.io/narwhals/api-reference/implementation/
[`Protocol`]: https://typing.python.org/en/latest/spec/protocol.html
[assignable to]: https://typing.python.org/en/latest/spec/glossary.html#term-assignable
[`TypeAlias`]: https://mypy.readthedocs.io/en/stable/kinds_of_types.html#type-aliases
[`TypeVar`]: https://mypy.readthedocs.io/en/stable/generics.html#type-variables-with-upper-bounds
[generic class]: https://docs.python.org/3/library/typing.html#user-defined-generic-types
[`@overload`s]: https://typing.python.org/en/latest/spec/overload.html
[ibis#9276 comment]: https://github.com/ibis-project/ibis/issues/9276#issuecomment-3292016818
[Type guards]: https://typing.python.org/en/latest/spec/narrowing.html
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Iterable, Sized
from typing import TYPE_CHECKING, Any, Protocol, TypeVar, Union, cast

from narwhals.dependencies import (
    get_cudf,
    get_modin,
    get_pandas,
    get_polars,
    get_pyarrow,
    is_dask_dataframe,
    is_duckdb_relation,
    is_ibis_table,
    is_pyspark_connect_dataframe,
    is_pyspark_dataframe,
    is_sqlframe_dataframe,
)

if TYPE_CHECKING:
    import duckdb
    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from sqlframe.base.dataframe import BaseDataFrame as _BaseDataFrame
    from typing_extensions import Self, TypeAlias, TypeIs

    SQLFrameDataFrame = _BaseDataFrame[Any, Any, Any, Any, Any]
    T = TypeVar("T")
    _Guard: TypeAlias = "Callable[[Any], TypeIs[T]]"
    Incomplete: TypeAlias = Any

__all__ = [
    "IntoDataFrame",
    "IntoDataFrameT",
    "IntoFrame",
    "IntoFrameT",
    "IntoLazyFrame",
    "IntoLazyFrameT",
    "IntoSeries",
    "IntoSeriesT",
    "NativeAny",
    "NativeArrow",
    "NativeCuDF",
    "NativeDask",
    "NativeDataFrame",
    "NativeDuckDB",
    "NativeFrame",
    "NativeIbis",
    "NativeKnown",
    "NativeLazyFrame",
    "NativeModin",
    "NativePandas",
    "NativePandasLike",
    "NativePandasLikeDataFrame",
    "NativePandasLikeSeries",
    "NativePolars",
    "NativePySpark",
    "NativePySparkConnect",
    "NativeSQLFrame",
    "NativeSeries",
    "NativeSparkLike",
    "NativeUnknown",
    "is_native_arrow",
    "is_native_cudf",
    "is_native_dask",
    "is_native_duckdb",
    "is_native_ibis",
    "is_native_modin",
    "is_native_pandas",
    "is_native_pandas_like",
    "is_native_polars",
    "is_native_pyspark",
    "is_native_pyspark_connect",
    "is_native_spark_like",
    "is_native_sqlframe",
]


# All dataframes supported by Narwhals have a
# `columns` property. Their similarities don't extend
# _that_ much further unfortunately...
class NativeFrame(Protocol):
    @property
    def columns(self) -> Any: ...
    def join(self, *args: Any, **kwargs: Any) -> Any: ...


class NativeDataFrame(Sized, NativeFrame, Protocol): ...


class NativeLazyFrame(NativeFrame, Protocol):
    def explain(self, *args: Any, **kwargs: Any) -> Any: ...


class NativeSeries(Sized, Iterable[Any], Protocol):
    def filter(self, *args: Any, **kwargs: Any) -> Any: ...


class _BasePandasLike(Sized, Protocol):
    index: Any
    """`mypy` doesn't like the asymmetric `property` setter in `pandas`."""

    def __getitem__(self, key: Any, /) -> Any: ...
    def __mul__(self, other: float | Collection[float] | Self, /) -> Self: ...
    def __floordiv__(self, other: float | Collection[float] | Self, /) -> Self: ...
    @property
    def loc(self) -> Any: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    def set_axis(self, labels: Any, *, axis: Any = ..., copy: bool = ...) -> Self: ...
    def copy(self, deep: bool = ...) -> Self: ...  # noqa: FBT001
    def rename(self, *args: Any, **kwds: Any) -> Self | Incomplete:
        """`mypy` & `pyright` disagree on overloads.

        `Incomplete` used to fix [more important issue](https://github.com/narwhals-dev/narwhals/pull/3016#discussion_r2296139744).
        """


class _BasePandasLikeFrame(NativeDataFrame, _BasePandasLike, Protocol): ...


class _BasePandasLikeSeries(NativeSeries, _BasePandasLike, Protocol):
    def where(self, cond: Any, other: Any = ..., /) -> Self | Incomplete: ...


class NativeDask(NativeLazyFrame, Protocol):
    _partition_type: type[pd.DataFrame]


class _CuDFDataFrame(_BasePandasLikeFrame, Protocol):
    def to_pylibcudf(self, *args: Any, **kwds: Any) -> Any: ...


class _CuDFSeries(_BasePandasLikeSeries, Protocol):
    def to_pylibcudf(self, *args: Any, **kwds: Any) -> Any: ...


class NativeIbis(NativeFrame, Protocol):
    def sql(self, *args: Any, **kwds: Any) -> Any: ...
    def __pyarrow_result__(self, *args: Any, **kwds: Any) -> Any: ...
    def __pandas_result__(self, *args: Any, **kwds: Any) -> Any: ...
    def __polars_result__(self, *args: Any, **kwds: Any) -> Any: ...


class _ModinDataFrame(_BasePandasLikeFrame, Protocol):
    _pandas_class: type[pd.DataFrame]


class _ModinSeries(_BasePandasLikeSeries, Protocol):
    _pandas_class: type[pd.Series[Any]]


class _PySparkDataFrame(NativeLazyFrame, Protocol):
    def dropDuplicatesWithinWatermark(self, *arg: Any, **kwargs: Any) -> Any: ...  # noqa: N802


NativePolars: TypeAlias = "pl.DataFrame | pl.LazyFrame | pl.Series"
NativeArrow: TypeAlias = "pa.Table | pa.ChunkedArray[Any]"
NativeDuckDB: TypeAlias = "duckdb.DuckDBPyRelation"
NativePandas: TypeAlias = "pd.DataFrame | pd.Series[Any]"
NativeModin: TypeAlias = "_ModinDataFrame | _ModinSeries"
NativeCuDF: TypeAlias = "_CuDFDataFrame | _CuDFSeries"
NativePandasLikeSeries: TypeAlias = "pd.Series[Any] | _CuDFSeries | _ModinSeries"
NativePandasLikeDataFrame: TypeAlias = "pd.DataFrame | _CuDFDataFrame | _ModinDataFrame"
NativePandasLike: TypeAlias = "NativePandasLikeDataFrame | NativePandasLikeSeries"
NativeSQLFrame: TypeAlias = "_BaseDataFrame[Any, Any, Any, Any, Any]"
NativePySpark: TypeAlias = _PySparkDataFrame
NativePySparkConnect: TypeAlias = _PySparkDataFrame
NativeSparkLike: TypeAlias = "NativeSQLFrame | NativePySpark | NativePySparkConnect"
NativeKnown: TypeAlias = "NativePolars | NativeArrow | NativePandasLike | NativeSparkLike | NativeDuckDB | NativeDask | NativeIbis"
NativeUnknown: TypeAlias = "NativeDataFrame | NativeSeries | NativeLazyFrame"
NativeAny: TypeAlias = "NativeKnown | NativeUnknown"

IntoDataFrame: TypeAlias = NativeDataFrame
"""Anything which can be converted to a Narwhals DataFrame.

Use this if your function accepts a narwhalifiable object but doesn't care about its backend.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoDataFrame
    >>> def agnostic_shape(df_native: IntoDataFrame) -> tuple[int, int]:
    ...     df = nw.from_native(df_native, eager_only=True)
    ...     return df.shape
"""

IntoLazyFrame: TypeAlias = Union[NativeLazyFrame, NativeIbis]
IntoFrame: TypeAlias = Union[IntoDataFrame, IntoLazyFrame]
"""Anything which can be converted to a Narwhals DataFrame or LazyFrame.

Use this if your function can accept an object which can be converted to either
`nw.DataFrame` or `nw.LazyFrame` and it doesn't care about its backend.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoFrame
    >>> def agnostic_columns(df_native: IntoFrame) -> list[str]:
    ...     df = nw.from_native(df_native)
    ...     return df.collect_schema().names()
"""

IntoSeries: TypeAlias = NativeSeries
"""Anything which can be converted to a Narwhals Series.

Use this if your function can accept an object which can be converted to `nw.Series`
and it doesn't care about its backend.

Examples:
    >>> from typing import Any
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoSeries
    >>> def agnostic_to_list(s_native: IntoSeries) -> list[Any]:
    ...     s = nw.from_native(s_native)
    ...     return s.to_list()
"""

IntoFrameT = TypeVar("IntoFrameT", bound=IntoFrame)
"""TypeVar bound to object convertible to Narwhals DataFrame or Narwhals LazyFrame.

Use this if your function accepts an object which is convertible to `nw.DataFrame`
or `nw.LazyFrame` and returns an object of the same type.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoFrameT
    >>> def agnostic_func(df_native: IntoFrameT) -> IntoFrameT:
    ...     df = nw.from_native(df_native)
    ...     return df.with_columns(c=nw.col("a") + 1).to_native()
"""

IntoDataFrameT = TypeVar("IntoDataFrameT", bound=IntoDataFrame)
"""TypeVar bound to object convertible to Narwhals DataFrame.

Use this if your function accepts an object which can be converted to `nw.DataFrame`
and returns an object of the same class.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoDataFrameT
    >>> def agnostic_func(df_native: IntoDataFrameT) -> IntoDataFrameT:
    ...     df = nw.from_native(df_native, eager_only=True)
    ...     return df.with_columns(c=df["a"] + 1).to_native()
"""

IntoLazyFrameT = TypeVar("IntoLazyFrameT", bound=IntoLazyFrame)
IntoSeriesT = TypeVar("IntoSeriesT", bound=IntoSeries)
"""TypeVar bound to object convertible to Narwhals Series.

Use this if your function accepts an object which can be converted to `nw.Series`
and returns an object of the same class.

Examples:
    >>> import narwhals as nw
    >>> from narwhals.typing import IntoSeriesT
    >>> def agnostic_abs(s_native: IntoSeriesT) -> IntoSeriesT:
    ...     s = nw.from_native(s_native, series_only=True)
    ...     return s.abs().to_native()
"""


def is_native_polars(obj: Any) -> TypeIs[NativePolars]:
    return (pl := get_polars()) is not None and isinstance(
        obj, (pl.DataFrame, pl.Series, pl.LazyFrame)
    )


def is_native_arrow(obj: Any) -> TypeIs[NativeArrow]:
    return (pa := get_pyarrow()) is not None and isinstance(
        obj, (pa.Table, pa.ChunkedArray)
    )


is_native_dask = cast("_Guard[NativeDask]", is_dask_dataframe)
is_native_duckdb: _Guard[NativeDuckDB] = is_duckdb_relation
is_native_sqlframe: _Guard[NativeSQLFrame] = is_sqlframe_dataframe
is_native_pyspark = cast("_Guard[NativePySpark]", is_pyspark_dataframe)
is_native_pyspark_connect = cast(
    "_Guard[NativePySparkConnect]", is_pyspark_connect_dataframe
)
is_native_ibis = cast("_Guard[NativeIbis]", is_ibis_table)


def is_native_pandas(obj: Any) -> TypeIs[NativePandas]:
    return (pd := get_pandas()) is not None and isinstance(obj, (pd.DataFrame, pd.Series))


def is_native_modin(obj: Any) -> TypeIs[NativeModin]:
    return (mpd := get_modin()) is not None and isinstance(
        obj, (mpd.DataFrame, mpd.Series)
    )


def is_native_cudf(obj: Any) -> TypeIs[NativeCuDF]:
    return (cudf := get_cudf()) is not None and isinstance(
        obj, (cudf.DataFrame, cudf.Series)
    )  # pragma: no cover


def is_native_pandas_like(obj: Any) -> TypeIs[NativePandasLike]:
    return is_native_pandas(obj) or is_native_cudf(obj) or is_native_modin(obj)


def is_native_spark_like(obj: Any) -> TypeIs[NativeSparkLike]:
    return (
        is_native_sqlframe(obj)
        or is_native_pyspark(obj)
        or is_native_pyspark_connect(obj)
    )
