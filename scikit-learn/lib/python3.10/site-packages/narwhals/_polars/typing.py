from __future__ import annotations  # pragma: no cover

from typing import TYPE_CHECKING  # pragma: no cover

if TYPE_CHECKING:
    from typing import TypeVar

    from narwhals._polars.dataframe import PolarsDataFrame, PolarsLazyFrame

    FrameT = TypeVar("FrameT", PolarsDataFrame, PolarsLazyFrame)
