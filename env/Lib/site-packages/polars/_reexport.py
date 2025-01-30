"""Re-export Polars functionality to avoid cyclical imports."""

from polars.dataframe import DataFrame
from polars.expr import Expr, When
from polars.lazyframe import LazyFrame
from polars.schema import Schema
from polars.series import Series

__all__ = [
    "DataFrame",
    "Expr",
    "LazyFrame",
    "Schema",
    "Series",
    "When",
]
