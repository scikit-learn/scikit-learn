from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from narwhals._compliant.typing import NativeExprT
from narwhals._compliant.when_then import CompliantThen, CompliantWhen
from narwhals._sql.typing import SQLExprT, SQLLazyFrameT

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from narwhals._compliant.typing import WindowFunction
    from narwhals._compliant.when_then import IntoExpr
    from narwhals._compliant.window import WindowInputs
    from narwhals._utils import _LimitedContext


class SQLWhen(
    CompliantWhen[SQLLazyFrameT, NativeExprT, SQLExprT],
    Protocol[SQLLazyFrameT, NativeExprT, SQLExprT],
):
    @property
    def _then(self) -> type[SQLThen[SQLLazyFrameT, NativeExprT, SQLExprT]]: ...

    def __call__(self, df: SQLLazyFrameT) -> Sequence[NativeExprT]:
        is_expr = self._condition._is_expr
        when = df.__narwhals_namespace__()._when
        lit = df.__narwhals_namespace__()._lit
        condition = df._evaluate_expr(self._condition)
        then_ = self._then_value
        then = df._evaluate_expr(then_) if is_expr(then_) else lit(then_)
        other_ = self._otherwise_value
        if other_ is None:
            result = when(condition, then)
        else:
            otherwise = df._evaluate_expr(other_) if is_expr(other_) else lit(other_)
            result = when(condition, then).otherwise(otherwise)
        return [result]

    @classmethod
    def from_expr(cls, condition: SQLExprT, /, *, context: _LimitedContext) -> Self:
        obj = cls.__new__(cls)
        obj._condition = condition
        obj._then_value = None
        obj._otherwise_value = None
        obj._implementation = context._implementation
        obj._version = context._version
        return obj

    def _window_function(
        self, df: SQLLazyFrameT, window_inputs: WindowInputs[NativeExprT]
    ) -> Sequence[NativeExprT]:
        when = df.__narwhals_namespace__()._when
        lit = df.__narwhals_namespace__()._lit
        is_expr = self._condition._is_expr
        condition = self._condition.window_function(df, window_inputs)[0]
        then_ = self._then_value
        then = (
            then_.window_function(df, window_inputs)[0] if is_expr(then_) else lit(then_)
        )

        other_ = self._otherwise_value
        if other_ is None:
            result = when(condition, then)
        else:
            other = (
                other_.window_function(df, window_inputs)[0]
                if is_expr(other_)
                else lit(other_)
            )
            result = when(condition, then).otherwise(other)
        return [result]


class SQLThen(
    CompliantThen[
        SQLLazyFrameT,
        NativeExprT,
        SQLExprT,
        SQLWhen[SQLLazyFrameT, NativeExprT, SQLExprT],
    ],
    Protocol[SQLLazyFrameT, NativeExprT, SQLExprT],
):
    _window_function: WindowFunction[SQLLazyFrameT, NativeExprT] | None

    @classmethod
    def from_when(
        cls,
        when: SQLWhen[SQLLazyFrameT, NativeExprT, SQLExprT],
        then: IntoExpr[NativeExprT, SQLExprT],
        /,
    ) -> Self:
        when._then_value = then
        obj = cls.__new__(cls)
        obj._call = when
        obj._window_function = when._window_function
        obj._when_value = when
        obj._evaluate_output_names = getattr(
            then, "_evaluate_output_names", lambda _df: ["literal"]
        )
        obj._alias_output_names = getattr(then, "_alias_output_names", None)
        obj._implementation = when._implementation
        obj._version = when._version
        return obj
