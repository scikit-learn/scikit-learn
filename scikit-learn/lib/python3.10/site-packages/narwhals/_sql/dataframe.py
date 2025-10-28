from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol

from narwhals._compliant.dataframe import CompliantLazyFrame
from narwhals._compliant.typing import (
    CompliantExprT_contra,
    NativeExprT,
    NativeLazyFrameT,
)
from narwhals._translate import ToNarwhalsT_co
from narwhals._utils import check_columns_exist

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self, TypeAlias

    from narwhals._compliant.window import WindowInputs
    from narwhals._sql.expr import SQLExpr
    from narwhals.exceptions import ColumnNotFoundError

    Incomplete: TypeAlias = Any


class SQLLazyFrame(
    CompliantLazyFrame[CompliantExprT_contra, NativeLazyFrameT, ToNarwhalsT_co],
    Protocol[CompliantExprT_contra, NativeLazyFrameT, ToNarwhalsT_co],
):
    def _evaluate_window_expr(
        self,
        expr: SQLExpr[Self, NativeExprT],
        /,
        window_inputs: WindowInputs[NativeExprT],
    ) -> NativeExprT:
        result = expr.window_function(self, window_inputs)
        assert len(result) == 1  # debug assertion  # noqa: S101
        return result[0]

    def _evaluate_expr(self, expr: CompliantExprT_contra, /) -> Any:
        result = expr(self)
        assert len(result) == 1  # debug assertion  # noqa: S101
        return result[0]

    def _check_columns_exist(self, subset: Sequence[str]) -> ColumnNotFoundError | None:
        return check_columns_exist(subset, available=self.columns)
