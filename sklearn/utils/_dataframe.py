"""Functions to determine if an object is a dataframe or series."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import sys


def is_df_or_series(X):
    """Return True if the X is a dataframe or series."""
    return is_pandas_df_or_series(X) or is_polars_df_or_series(X) or is_pyarrow_data(X)


def is_pandas_df_or_series(X):
    """Return True if the X is a pandas dataframe or series."""
    try:
        pd = sys.modules["pandas"]
    except KeyError:
        return False
    return isinstance(X, (pd.DataFrame, pd.Series))


def is_pandas_df(X):
    """Return True if the X is a pandas dataframe."""
    try:
        pd = sys.modules["pandas"]
    except KeyError:
        return False
    return isinstance(X, pd.DataFrame)


def is_pyarrow_data(X):
    """Return True if the X is a pyarrow Table, RecordBatch, Array or ChunkedArray."""
    try:
        pa = sys.modules["pyarrow"]
    except KeyError:
        return False
    return isinstance(X, (pa.Table, pa.RecordBatch, pa.Array, pa.ChunkedArray))


def is_polars_df_or_series(X):
    """Return True if the X is a polars dataframe or series."""
    try:
        pl = sys.modules["polars"]
    except KeyError:
        return False
    return isinstance(X, (pl.DataFrame, pl.Series))


def is_polars_df(X):
    """Return True if the X is a polars dataframe."""
    try:
        pl = sys.modules["polars"]
    except KeyError:
        return False
    return isinstance(X, pl.DataFrame)
