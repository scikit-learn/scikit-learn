from __future__ import annotations

from types import ModuleType
from typing import TYPE_CHECKING, Literal, Union

from narwhals._typing_compat import TypeVar
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from typing_extensions import TypeAlias

# `str` aliases
_Polars: TypeAlias = Literal["polars"]
_Arrow: TypeAlias = Literal["pyarrow"]
_Dask: TypeAlias = Literal["dask"]
_DuckDB: TypeAlias = Literal["duckdb"]
_Pandas: TypeAlias = Literal["pandas"]
_Modin: TypeAlias = Literal["modin"]
_CuDF: TypeAlias = Literal["cudf"]
_PySpark: TypeAlias = Literal["pyspark"]
_SQLFrame: TypeAlias = Literal["sqlframe"]
_PySparkConnect: TypeAlias = Literal["pyspark[connect]"]
_Ibis: TypeAlias = Literal["ibis"]
_PandasLike: TypeAlias = Literal[_Pandas, _CuDF, _Modin]
_SparkLike: TypeAlias = Literal[_PySpark, _SQLFrame, _PySparkConnect]
_EagerOnly: TypeAlias = Literal[_PandasLike, _Arrow]
_EagerAllowed: TypeAlias = Literal[_Polars, _EagerOnly]
_LazyOnly: TypeAlias = Literal[_SparkLike, _Dask, _DuckDB, _Ibis]
_LazyAllowed: TypeAlias = Literal[_Polars, _LazyOnly]

# `Implementation` aliases
_PandasImpl: TypeAlias = Literal[Implementation.PANDAS]
_ModinImpl: TypeAlias = Literal[Implementation.MODIN]
_CuDFImpl: TypeAlias = Literal[Implementation.CUDF]
_PySparkImpl: TypeAlias = Literal[Implementation.PYSPARK]
_SQLFrameImpl: TypeAlias = Literal[Implementation.SQLFRAME]
_PySparkConnectImpl: TypeAlias = Literal[Implementation.PYSPARK_CONNECT]
_PolarsImpl: TypeAlias = Literal[Implementation.POLARS]
_ArrowImpl: TypeAlias = Literal[Implementation.PYARROW]
_DaskImpl: TypeAlias = Literal[Implementation.DASK]
_DuckDBImpl: TypeAlias = Literal[Implementation.DUCKDB]
_IbisImpl: TypeAlias = Literal[Implementation.IBIS]
_PandasLikeImpl: TypeAlias = Literal[_PandasImpl, _CuDFImpl, _ModinImpl]
_SparkLikeImpl: TypeAlias = Literal[_PySparkImpl, _SQLFrameImpl, _PySparkConnectImpl]
_EagerOnlyImpl: TypeAlias = Literal[_PandasLikeImpl, _ArrowImpl]
_EagerAllowedImpl: TypeAlias = Literal[_EagerOnlyImpl, _PolarsImpl]  # noqa: PYI047
_LazyOnlyImpl: TypeAlias = Literal[_SparkLikeImpl, _DaskImpl, _DuckDBImpl, _IbisImpl]
_LazyAllowedImpl: TypeAlias = Literal[_LazyOnlyImpl, _PolarsImpl]  # noqa: PYI047

# NOTE: Temporary aliases for gaps in `LazyFrame.collect`, `DataFrame.lazy`, see:
# - https://github.com/narwhals-dev/narwhals/pull/2971#discussion_r2277137003
# - https://github.com/narwhals-dev/narwhals/pull/3002#issuecomment-3194267667
_LazyFrameCollectImpl: TypeAlias = Literal[_PandasImpl, _PolarsImpl, _ArrowImpl]  # noqa: PYI047

# `str | Implementation` aliases
Pandas: TypeAlias = Literal[_Pandas, _PandasImpl]
CuDF: TypeAlias = Literal[_CuDF, _CuDFImpl]
Modin: TypeAlias = Literal[_Modin, _ModinImpl]
PySpark: TypeAlias = Literal[_PySpark, _PySparkImpl]
SQLFrame: TypeAlias = Literal[_SQLFrame, _SQLFrameImpl]
PySparkConnect: TypeAlias = Literal[_PySparkConnect, _PySparkConnectImpl]
Polars: TypeAlias = Literal[_Polars, _PolarsImpl]
Arrow: TypeAlias = Literal[_Arrow, _ArrowImpl]
Dask: TypeAlias = Literal[_Dask, _DaskImpl]
DuckDB: TypeAlias = Literal[_DuckDB, _DuckDBImpl]
Ibis: TypeAlias = Literal[_Ibis, _IbisImpl]
PandasLike: TypeAlias = Literal[_PandasLike, _PandasLikeImpl]
SparkLike: TypeAlias = Literal[_SparkLike, _SparkLikeImpl]
EagerOnly: TypeAlias = Literal[PandasLike, Arrow]
LazyOnly: TypeAlias = Literal[SparkLike, Dask, DuckDB, Ibis]
EagerAllowed: TypeAlias = Literal[EagerOnly, Polars]
"""A string name or [`narwhals.Implementation`][] of an eager backend.

- A string name, one of: `"cudf"`, `"modin"`, `"pandas"`, `"pyarrow"`, `"polars"`.
- An Implementation, such as: `Implementation.CUDF`, `Implementation.MODIN`, ...
"""

LazyAllowed: TypeAlias = Literal[LazyOnly, Polars]
"""A string name or [`narwhals.Implementation`][] of a lazy backend.

- A string name, such as: `"duckdb"`, `"ibis"`, `"dask"`, `"sqlframe"`, ...
- An Implementation, such as: `Implementation.POLARS`, `Implementation.PYSPARK`, ...
"""

BackendName: TypeAlias = Literal[_EagerAllowed, _LazyAllowed]
Backend: TypeAlias = Literal[EagerAllowed, LazyAllowed]
"""A string name or [`narwhals.Implementation`][] of a supported
backend (either eager or lazy).

- A string name, such as: `"duckdb"`, `"ibis"`, `"pandas"`, `"pyarrow"`, `"polars"`, ...
- An Implementation, such as: `Implementation.DASK`, `Implementation.PYSPARK`, ...
"""

BackendT = TypeVar("BackendT", bound=Backend)
IntoBackend: TypeAlias = Union[BackendT, ModuleType]
"""Anything that can be converted into a [`narwhals.Implementation`][].

`backend` can be specified in three ways.

Examples:
    A string backend name, such as: `"pandas"`, `"pyarrow"`, `"modin"`, `"cudf"`

    >>> import pandas as pd
    >>> import narwhals as nw
    >>>
    >>> data = {"c": [5, 2], "d": [1, 4]}
    >>> nw.DataFrame.from_dict(data, backend="pandas")
    ┌──────────────────┐
    |Narwhals DataFrame|
    |------------------|
    |        c  d      |
    |     0  5  1      |
    |     1  2  4      |
    └──────────────────┘

    An Implementation, such as: `Implementation.POLARS`, `Implementation.DUCKDB`, `Implementation.PYSPARK`

    >>> import narwhals as nw
    >>> nw.read_parquet("file.parquet", backend=nw.Implementation.PYARROW)
    ┌──────────────────┐
    |Narwhals DataFrame|
    |------------------|
    |  pyarrow.Table   |
    |  a: int64        |
    |  b: int64        |
    |  ----            |
    |  a: [[1,2]]      |
    |  b: [[4,5]]      |
    └──────────────────┘

    A python module, such as `dask`, `ibis`, `sqlframe`

    >>> import numpy as np
    >>> import polars as pl
    >>> import narwhals as nw
    >>>
    >>> arr = np.arange(5, 10)
    >>> nw.Series.from_numpy("arr", arr, dtype=nw.Int8, backend=pl)
    ┌──────────────────┐
    | Narwhals Series  |
    |------------------|
    |shape: (5,)       |
    |Series: 'arr' [i8]|
    |[                 |
    |        5         |
    |        6         |
    |        7         |
    |        8         |
    |        9         |
    |]                 |
    └──────────────────┘
"""

IntoBackendAny: TypeAlias = IntoBackend[Backend]
IntoBackendEager: TypeAlias = IntoBackend[EagerAllowed]
IntoBackendLazy: TypeAlias = IntoBackend[LazyAllowed]
