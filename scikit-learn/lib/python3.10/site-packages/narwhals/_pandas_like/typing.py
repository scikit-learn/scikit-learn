from __future__ import annotations  # pragma: no cover

from typing import TYPE_CHECKING  # pragma: no cover

from narwhals._typing_compat import TypeVar

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd
    from typing_extensions import TypeAlias

    from narwhals._native import (
        NativePandasLikeDataFrame,
        _CuDFDataFrame,
        _CuDFSeries,
        _ModinDataFrame,
        _ModinSeries,
    )
    from narwhals._pandas_like.expr import PandasLikeExpr
    from narwhals._pandas_like.series import PandasLikeSeries

    IntoPandasLikeExpr: TypeAlias = "PandasLikeExpr | PandasLikeSeries"

NativeSeriesT = TypeVar(
    "NativeSeriesT",
    "pd.Series[Any]",
    "_CuDFSeries",
    "_ModinSeries",
    default="pd.Series[Any]",
)
NativeDataFrameT = TypeVar(
    "NativeDataFrameT", bound="NativePandasLikeDataFrame", default="pd.DataFrame"
)
NativeNDFrameT = TypeVar(
    "NativeNDFrameT",
    "pd.DataFrame",
    "pd.Series[Any]",
    "_CuDFDataFrame",
    "_CuDFSeries",
    "_ModinDataFrame",
    "_ModinSeries",
)
