from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, TypeVar, cast

from narwhals._compliant.expr import CompliantExpr
from narwhals._compliant.typing import (
    CompliantExprAny,
    CompliantFrameAny,
    CompliantSeriesOrNativeExprAny,
    EagerDataFrameT,
    EagerExprT,
    EagerSeriesT,
    LazyExprAny,
    NativeSeriesT,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self, TypeAlias

    from narwhals._compliant.typing import EvalSeries
    from narwhals._utils import Implementation, Version, _LimitedContext
    from narwhals.typing import NonNestedLiteral


__all__ = ["CompliantThen", "CompliantWhen", "EagerWhen"]

ExprT = TypeVar("ExprT", bound=CompliantExprAny)
LazyExprT = TypeVar("LazyExprT", bound=LazyExprAny)
SeriesT = TypeVar("SeriesT", bound=CompliantSeriesOrNativeExprAny)
FrameT = TypeVar("FrameT", bound=CompliantFrameAny)

Scalar: TypeAlias = Any
"""A native literal value."""

IntoExpr: TypeAlias = "SeriesT | ExprT | NonNestedLiteral | Scalar"
"""Anything that is convertible into a `CompliantExpr`."""


class CompliantWhen(Protocol[FrameT, SeriesT, ExprT]):
    _condition: ExprT
    _then_value: IntoExpr[SeriesT, ExprT]
    _otherwise_value: IntoExpr[SeriesT, ExprT] | None
    _implementation: Implementation
    _version: Version

    @property
    def _then(self) -> type[CompliantThen[FrameT, SeriesT, ExprT, Self]]: ...
    def __call__(self, compliant_frame: FrameT, /) -> Sequence[SeriesT]: ...
    def then(
        self, value: IntoExpr[SeriesT, ExprT], /
    ) -> CompliantThen[FrameT, SeriesT, ExprT, Self]:
        return self._then.from_when(self, value)

    @classmethod
    def from_expr(cls, condition: ExprT, /, *, context: _LimitedContext) -> Self:
        obj = cls.__new__(cls)
        obj._condition = condition
        obj._then_value = None
        obj._otherwise_value = None
        obj._implementation = context._implementation
        obj._version = context._version
        return obj


WhenT_contra = TypeVar(
    "WhenT_contra", bound=CompliantWhen[Any, Any, Any], contravariant=True
)


class CompliantThen(
    CompliantExpr[FrameT, SeriesT], Protocol[FrameT, SeriesT, ExprT, WhenT_contra]
):
    _call: EvalSeries[FrameT, SeriesT]
    _when_value: CompliantWhen[FrameT, SeriesT, ExprT]
    _implementation: Implementation
    _version: Version

    @classmethod
    def from_when(cls, when: WhenT_contra, then: IntoExpr[SeriesT, ExprT], /) -> Self:
        when._then_value = then
        obj = cls.__new__(cls)
        obj._call = when
        obj._when_value = when
        obj._evaluate_output_names = getattr(
            then, "_evaluate_output_names", lambda _df: ["literal"]
        )
        obj._alias_output_names = getattr(then, "_alias_output_names", None)
        obj._implementation = when._implementation
        obj._version = when._version
        return obj

    def otherwise(self, otherwise: IntoExpr[SeriesT, ExprT], /) -> ExprT:
        self._when_value._otherwise_value = otherwise
        return cast("ExprT", self)


class EagerWhen(
    CompliantWhen[EagerDataFrameT, EagerSeriesT, EagerExprT],
    Protocol[EagerDataFrameT, EagerSeriesT, EagerExprT, NativeSeriesT],
):
    def _if_then_else(
        self,
        when: NativeSeriesT,
        then: NativeSeriesT,
        otherwise: NativeSeriesT | NonNestedLiteral | Scalar,
        /,
    ) -> NativeSeriesT: ...

    def __call__(self, df: EagerDataFrameT, /) -> Sequence[EagerSeriesT]:
        is_expr = self._condition._is_expr
        when: EagerSeriesT = self._condition(df)[0]
        then: EagerSeriesT
        align = when._align_full_broadcast

        if is_expr(self._then_value):
            then = self._then_value(df)[0]
        else:
            then = when.alias("literal")._from_scalar(self._then_value)
            then._broadcast = True

        if is_expr(self._otherwise_value):
            otherwise = self._otherwise_value(df)[0]
            when, then, otherwise = align(when, then, otherwise)
            result = self._if_then_else(when.native, then.native, otherwise.native)
        else:
            when, then = align(when, then)
            result = self._if_then_else(when.native, then.native, self._otherwise_value)
        return [then._with_native(result)]
