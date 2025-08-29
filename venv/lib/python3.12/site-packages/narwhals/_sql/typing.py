from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

if TYPE_CHECKING:
    from narwhals._sql.dataframe import SQLLazyFrame
    from narwhals._sql.expr import SQLExpr

    SQLExprAny = SQLExpr[Any, Any]
    SQLLazyFrameAny = SQLLazyFrame[Any, Any, Any]

SQLExprT = TypeVar("SQLExprT", bound="SQLExprAny")
SQLExprT_contra = TypeVar("SQLExprT_contra", bound="SQLExprAny", contravariant=True)
SQLLazyFrameT = TypeVar("SQLLazyFrameT", bound="SQLLazyFrameAny")
