from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar, Union

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

if TYPE_CHECKING:
    from typing_extensions import TypeAlias

    from narwhals.stable.v2 import DataFrame, Expr, LazyFrame, Series


IntoExpr: TypeAlias = Union["Expr", str, "Series[Any]"]
"""Anything which can be converted to an expression.

Use this to mean "either a Narwhals expression, or something
which can be converted into one". For example, `exprs` in `DataFrame.select` is
typed to accept `IntoExpr`, as it can either accept a `nw.Expr`
(e.g. `df.select(nw.col('a'))`) or a string which will be interpreted as a
`nw.Expr`, e.g. `df.select('a')`.
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


__all__ = [
    "DataFrameT",
    "Frame",
    "FrameT",
    "IntoDataFrame",
    "IntoDataFrameT",
    "IntoExpr",
    "IntoFrame",
    "IntoFrameT",
    "IntoLazyFrame",
    "IntoLazyFrameT",
    "IntoSeries",
    "IntoSeriesT",
]
