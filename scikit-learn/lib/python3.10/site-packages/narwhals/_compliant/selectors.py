"""Almost entirely complete, generic `selectors` implementation."""

from __future__ import annotations

import re
from functools import partial
from typing import TYPE_CHECKING, Protocol, TypeVar, overload

from narwhals._compliant.expr import CompliantExpr
from narwhals._utils import (
    _parse_time_unit_and_time_zone,
    dtype_matches_time_unit_and_time_zone,
    get_column_names,
    is_compliant_dataframe,
    zip_strict,
)

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable, Iterator, Sequence
    from datetime import timezone

    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.expr import NativeExpr
    from narwhals._compliant.typing import (
        CompliantDataFrameAny,
        CompliantExprAny,
        CompliantFrameAny,
        CompliantLazyFrameAny,
        CompliantSeriesAny,
        CompliantSeriesOrNativeExprAny,
        EvalNames,
        EvalSeries,
    )
    from narwhals._utils import Implementation, Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.typing import TimeUnit

__all__ = [
    "CompliantSelector",
    "CompliantSelectorNamespace",
    "EagerSelectorNamespace",
    "LazySelectorNamespace",
]


SeriesOrExprT = TypeVar("SeriesOrExprT", bound="CompliantSeriesOrNativeExprAny")
SeriesT = TypeVar("SeriesT", bound="CompliantSeriesAny")
ExprT = TypeVar("ExprT", bound="NativeExpr")
FrameT = TypeVar("FrameT", bound="CompliantFrameAny")
DataFrameT = TypeVar("DataFrameT", bound="CompliantDataFrameAny")
LazyFrameT = TypeVar("LazyFrameT", bound="CompliantLazyFrameAny")
SelectorOrExpr: TypeAlias = (
    "CompliantSelector[FrameT, SeriesOrExprT] | CompliantExpr[FrameT, SeriesOrExprT]"
)


class CompliantSelectorNamespace(Protocol[FrameT, SeriesOrExprT]):
    # NOTE: `narwhals`
    _implementation: Implementation
    _version: Version

    @property
    def _selector(self) -> type[CompliantSelector[FrameT, SeriesOrExprT]]: ...
    @classmethod
    def from_namespace(cls, context: _LimitedContext, /) -> Self:
        obj = cls.__new__(cls)
        obj._implementation = context._implementation
        obj._version = context._version
        return obj

    def _iter_columns(self, df: FrameT, /) -> Iterator[SeriesOrExprT]: ...
    def _iter_schema(self, df: FrameT, /) -> Iterator[tuple[str, DType]]: ...
    def _iter_columns_dtypes(
        self, df: FrameT, /
    ) -> Iterator[tuple[SeriesOrExprT, DType]]: ...
    def _iter_columns_names(self, df: FrameT, /) -> Iterator[tuple[SeriesOrExprT, str]]:
        yield from zip_strict(self._iter_columns(df), df.columns)

    def _is_dtype(
        self: CompliantSelectorNamespace[FrameT, SeriesOrExprT], dtype: type[DType], /
    ) -> CompliantSelector[FrameT, SeriesOrExprT]:
        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            return [
                ser for ser, tp in self._iter_columns_dtypes(df) if isinstance(tp, dtype)
            ]

        def names(df: FrameT) -> Sequence[str]:
            return [name for name, tp in self._iter_schema(df) if isinstance(tp, dtype)]

        return self._selector.from_callables(series, names, context=self)

    # NOTE: `polars`
    def by_dtype(
        self, dtypes: Collection[DType | type[DType]]
    ) -> CompliantSelector[FrameT, SeriesOrExprT]:
        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            return [ser for ser, tp in self._iter_columns_dtypes(df) if tp in dtypes]

        def names(df: FrameT) -> Sequence[str]:
            return [name for name, tp in self._iter_schema(df) if tp in dtypes]

        return self._selector.from_callables(series, names, context=self)

    def matches(self, pattern: str) -> CompliantSelector[FrameT, SeriesOrExprT]:
        p = re.compile(pattern)

        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            if (
                is_compliant_dataframe(df)
                and not self._implementation.is_duckdb()
                and not self._implementation.is_ibis()
            ):
                return [df.get_column(col) for col in df.columns if p.search(col)]

            return [ser for ser, name in self._iter_columns_names(df) if p.search(name)]

        def names(df: FrameT) -> Sequence[str]:
            return [col for col in df.columns if p.search(col)]

        return self._selector.from_callables(series, names, context=self)

    def numeric(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            return [ser for ser, tp in self._iter_columns_dtypes(df) if tp.is_numeric()]

        def names(df: FrameT) -> Sequence[str]:
            return [name for name, tp in self._iter_schema(df) if tp.is_numeric()]

        return self._selector.from_callables(series, names, context=self)

    def categorical(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        return self._is_dtype(self._version.dtypes.Categorical)

    def string(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        return self._is_dtype(self._version.dtypes.String)

    def boolean(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        return self._is_dtype(self._version.dtypes.Boolean)

    def all(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            return list(self._iter_columns(df))

        return self._selector.from_callables(series, get_column_names, context=self)

    def datetime(
        self,
        time_unit: TimeUnit | Iterable[TimeUnit] | None,
        time_zone: str | timezone | Iterable[str | timezone | None] | None,
    ) -> CompliantSelector[FrameT, SeriesOrExprT]:
        time_units, time_zones = _parse_time_unit_and_time_zone(time_unit, time_zone)
        matches = partial(
            dtype_matches_time_unit_and_time_zone,
            dtypes=self._version.dtypes,
            time_units=time_units,
            time_zones=time_zones,
        )

        def series(df: FrameT) -> Sequence[SeriesOrExprT]:
            return [ser for ser, tp in self._iter_columns_dtypes(df) if matches(tp)]

        def names(df: FrameT) -> Sequence[str]:
            return [name for name, tp in self._iter_schema(df) if matches(tp)]

        return self._selector.from_callables(series, names, context=self)


class EagerSelectorNamespace(
    CompliantSelectorNamespace[DataFrameT, SeriesT], Protocol[DataFrameT, SeriesT]
):
    def _iter_schema(self, df: DataFrameT, /) -> Iterator[tuple[str, DType]]:
        for ser in self._iter_columns(df):
            yield ser.name, ser.dtype

    def _iter_columns(self, df: DataFrameT, /) -> Iterator[SeriesT]:
        yield from df.iter_columns()

    def _iter_columns_dtypes(self, df: DataFrameT, /) -> Iterator[tuple[SeriesT, DType]]:
        for ser in self._iter_columns(df):
            yield ser, ser.dtype


class LazySelectorNamespace(
    CompliantSelectorNamespace[LazyFrameT, ExprT], Protocol[LazyFrameT, ExprT]
):
    def _iter_schema(self, df: LazyFrameT) -> Iterator[tuple[str, DType]]:
        yield from df.schema.items()

    def _iter_columns(self, df: LazyFrameT) -> Iterator[ExprT]:
        yield from df._iter_columns()

    def _iter_columns_dtypes(self, df: LazyFrameT, /) -> Iterator[tuple[ExprT, DType]]:
        yield from zip_strict(self._iter_columns(df), df.schema.values())


class CompliantSelector(
    CompliantExpr[FrameT, SeriesOrExprT], Protocol[FrameT, SeriesOrExprT]
):
    _call: EvalSeries[FrameT, SeriesOrExprT]
    _function_name: str
    _implementation: Implementation
    _version: Version

    @classmethod
    def from_callables(
        cls,
        call: EvalSeries[FrameT, SeriesOrExprT],
        evaluate_output_names: EvalNames[FrameT],
        *,
        context: _LimitedContext,
    ) -> Self:
        obj = cls.__new__(cls)
        obj._call = call
        obj._evaluate_output_names = evaluate_output_names
        obj._alias_output_names = None
        obj._implementation = context._implementation
        obj._version = context._version
        return obj

    @property
    def selectors(self) -> CompliantSelectorNamespace[FrameT, SeriesOrExprT]:
        return self.__narwhals_namespace__().selectors

    def _to_expr(self) -> CompliantExpr[FrameT, SeriesOrExprT]: ...

    def _is_selector(
        self, other: Self | CompliantExpr[FrameT, SeriesOrExprT]
    ) -> TypeIs[CompliantSelector[FrameT, SeriesOrExprT]]:
        return isinstance(other, type(self))

    @overload
    def __sub__(self, other: Self) -> Self: ...
    @overload
    def __sub__(
        self, other: CompliantExpr[FrameT, SeriesOrExprT]
    ) -> CompliantExpr[FrameT, SeriesOrExprT]: ...
    def __sub__(
        self, other: SelectorOrExpr[FrameT, SeriesOrExprT]
    ) -> SelectorOrExpr[FrameT, SeriesOrExprT]:
        if self._is_selector(other):

            def series(df: FrameT) -> Sequence[SeriesOrExprT]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [
                    x
                    for x, name in zip_strict(self(df), lhs_names)
                    if name not in rhs_names
                ]

            def names(df: FrameT) -> Sequence[str]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [x for x in lhs_names if x not in rhs_names]

            return self.selectors._selector.from_callables(series, names, context=self)
        return self._to_expr() - other

    @overload
    def __or__(self, other: Self) -> Self: ...
    @overload
    def __or__(
        self, other: CompliantExpr[FrameT, SeriesOrExprT]
    ) -> CompliantExpr[FrameT, SeriesOrExprT]: ...
    def __or__(
        self, other: SelectorOrExpr[FrameT, SeriesOrExprT]
    ) -> SelectorOrExpr[FrameT, SeriesOrExprT]:
        if self._is_selector(other):

            def series(df: FrameT) -> Sequence[SeriesOrExprT]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [
                    *(
                        x
                        for x, name in zip_strict(self(df), lhs_names)
                        if name not in rhs_names
                    ),
                    *other(df),
                ]

            def names(df: FrameT) -> Sequence[str]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [*(x for x in lhs_names if x not in rhs_names), *rhs_names]

            return self.selectors._selector.from_callables(series, names, context=self)
        return self._to_expr() | other

    @overload
    def __and__(self, other: Self) -> Self: ...
    @overload
    def __and__(
        self, other: CompliantExpr[FrameT, SeriesOrExprT]
    ) -> CompliantExpr[FrameT, SeriesOrExprT]: ...
    def __and__(
        self, other: SelectorOrExpr[FrameT, SeriesOrExprT]
    ) -> SelectorOrExpr[FrameT, SeriesOrExprT]:
        if self._is_selector(other):

            def series(df: FrameT) -> Sequence[SeriesOrExprT]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [
                    x for x, name in zip_strict(self(df), lhs_names) if name in rhs_names
                ]

            def names(df: FrameT) -> Sequence[str]:
                lhs_names, rhs_names = _eval_lhs_rhs(df, self, other)
                return [x for x in lhs_names if x in rhs_names]

            return self.selectors._selector.from_callables(series, names, context=self)
        return self._to_expr() & other

    def __invert__(self) -> CompliantSelector[FrameT, SeriesOrExprT]:
        return self.selectors.all() - self


def _eval_lhs_rhs(
    df: CompliantFrameAny, lhs: CompliantExprAny, rhs: CompliantExprAny
) -> tuple[Sequence[str], Sequence[str]]:
    return lhs._evaluate_output_names(df), rhs._evaluate_output_names(df)
