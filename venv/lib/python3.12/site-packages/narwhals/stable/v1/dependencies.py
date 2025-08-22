from __future__ import annotations

import sys
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import cudf
    import dask.dataframe as dd
    import ibis
    import modin.pandas as mpd
    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import TypeIs


from narwhals.dependencies import (
    IMPORT_HOOKS,
    get_cudf,
    get_dask_dataframe,
    get_ibis,
    get_modin,
    get_numpy,
    get_pandas,
    get_polars,
    get_pyarrow,
    is_into_dataframe,
    is_into_series,
    is_narwhals_dataframe,
    is_narwhals_lazyframe,
    is_narwhals_series,
    is_numpy_array,
    is_pandas_index,
)


def is_pandas_dataframe(df: Any) -> TypeIs[pd.DataFrame]:
    """Check whether `df` is a pandas DataFrame without importing pandas."""
    return ((pd := get_pandas()) is not None and isinstance(df, pd.DataFrame)) or any(
        (mod := sys.modules.get(module_name, None)) is not None
        and isinstance(df, mod.pandas.DataFrame)
        for module_name in IMPORT_HOOKS
    )


def is_pandas_series(ser: Any) -> TypeIs[pd.Series[Any]]:
    """Check whether `ser` is a pandas Series without importing pandas."""
    return ((pd := get_pandas()) is not None and isinstance(ser, pd.Series)) or any(
        (mod := sys.modules.get(module_name, None)) is not None
        and isinstance(ser, mod.pandas.Series)
        for module_name in IMPORT_HOOKS
    )


def is_modin_dataframe(df: Any) -> TypeIs[mpd.DataFrame]:
    """Check whether `df` is a modin DataFrame without importing modin."""
    return (mpd := get_modin()) is not None and isinstance(df, mpd.DataFrame)


def is_modin_series(ser: Any) -> TypeIs[mpd.Series]:
    """Check whether `ser` is a modin Series without importing modin."""
    return (mpd := get_modin()) is not None and isinstance(ser, mpd.Series)


def is_cudf_dataframe(df: Any) -> TypeIs[cudf.DataFrame]:
    """Check whether `df` is a cudf DataFrame without importing cudf."""
    return (cudf := get_cudf()) is not None and isinstance(df, cudf.DataFrame)


def is_cudf_series(ser: Any) -> TypeIs[cudf.Series[Any]]:
    """Check whether `ser` is a cudf Series without importing cudf."""
    return (cudf := get_cudf()) is not None and isinstance(ser, cudf.Series)


def is_dask_dataframe(df: Any) -> TypeIs[dd.DataFrame]:
    """Check whether `df` is a Dask DataFrame without importing Dask."""
    return (dd := get_dask_dataframe()) is not None and isinstance(df, dd.DataFrame)


def is_ibis_table(df: Any) -> TypeIs[ibis.Table]:
    """Check whether `df` is a Ibis Table without importing Ibis."""
    return (ibis := get_ibis()) is not None and isinstance(df, ibis.expr.types.Table)


def is_polars_dataframe(df: Any) -> TypeIs[pl.DataFrame]:
    """Check whether `df` is a Polars DataFrame without importing Polars."""
    return (pl := get_polars()) is not None and isinstance(df, pl.DataFrame)


def is_polars_lazyframe(df: Any) -> TypeIs[pl.LazyFrame]:
    """Check whether `df` is a Polars LazyFrame without importing Polars."""
    return (pl := get_polars()) is not None and isinstance(df, pl.LazyFrame)


def is_polars_series(ser: Any) -> TypeIs[pl.Series]:
    """Check whether `ser` is a Polars Series without importing Polars."""
    return (pl := get_polars()) is not None and isinstance(ser, pl.Series)


def is_pyarrow_chunked_array(ser: Any) -> TypeIs[pa.ChunkedArray[Any]]:
    """Check whether `ser` is a PyArrow ChunkedArray without importing PyArrow."""
    return (pa := get_pyarrow()) is not None and isinstance(ser, pa.ChunkedArray)


def is_pyarrow_table(df: Any) -> TypeIs[pa.Table]:
    """Check whether `df` is a PyArrow Table without importing PyArrow."""
    return (pa := get_pyarrow()) is not None and isinstance(df, pa.Table)


def is_pandas_like_dataframe(df: Any) -> bool:
    """Check whether `df` is a pandas-like DataFrame without doing any imports.

    By "pandas-like", we mean: pandas, Modin, cuDF.
    """
    return is_pandas_dataframe(df) or is_modin_dataframe(df) or is_cudf_dataframe(df)


def is_pandas_like_series(ser: Any) -> bool:
    """Check whether `ser` is a pandas-like Series without doing any imports.

    By "pandas-like", we mean: pandas, Modin, cuDF.
    """
    return is_pandas_series(ser) or is_modin_series(ser) or is_cudf_series(ser)


__all__ = [
    "get_cudf",
    "get_ibis",
    "get_modin",
    "get_numpy",
    "get_pandas",
    "get_polars",
    "get_pyarrow",
    "is_cudf_dataframe",
    "is_cudf_series",
    "is_dask_dataframe",
    "is_ibis_table",
    "is_into_dataframe",
    "is_into_series",
    "is_modin_dataframe",
    "is_modin_series",
    "is_narwhals_dataframe",
    "is_narwhals_lazyframe",
    "is_narwhals_series",
    "is_numpy_array",
    "is_pandas_dataframe",
    "is_pandas_index",
    "is_pandas_like_dataframe",
    "is_pandas_like_series",
    "is_pandas_series",
    "is_polars_dataframe",
    "is_polars_lazyframe",
    "is_polars_series",
    "is_pyarrow_chunked_array",
    "is_pyarrow_table",
]
