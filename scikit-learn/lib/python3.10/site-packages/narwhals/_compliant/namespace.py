from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Protocol, overload

from narwhals._compliant.typing import (
    CompliantExprT,
    CompliantFrameT,
    CompliantLazyFrameT,
    DepthTrackingExprT,
    EagerDataFrameT,
    EagerExprT,
    EagerSeriesT,
    LazyExprT,
    NativeFrameT,
    NativeSeriesT,
)
from narwhals._expression_parsing import is_expr, is_series
from narwhals._utils import (
    exclude_column_names,
    get_column_names,
    passthrough_column_names,
)
from narwhals.dependencies import is_numpy_array, is_numpy_array_2d

if TYPE_CHECKING:
    from collections.abc import Container, Iterable, Sequence

    from typing_extensions import TypeAlias, TypeIs

    from narwhals._compliant.selectors import CompliantSelectorNamespace
    from narwhals._compliant.when_then import CompliantWhen, EagerWhen
    from narwhals._utils import Implementation, Version
    from narwhals.expr import Expr
    from narwhals.series import Series
    from narwhals.typing import (
        ConcatMethod,
        Into1DArray,
        IntoDType,
        IntoSchema,
        NonNestedLiteral,
        _1DArray,
        _2DArray,
    )

    Incomplete: TypeAlias = Any

__all__ = [
    "CompliantNamespace",
    "DepthTrackingNamespace",
    "EagerNamespace",
    "LazyNamespace",
]


class CompliantNamespace(Protocol[CompliantFrameT, CompliantExprT]):
    # NOTE: `narwhals`
    _implementation: Implementation
    _version: Version

    @property
    def _expr(self) -> type[CompliantExprT]: ...
    def parse_into_expr(
        self, data: Expr | NonNestedLiteral | Any, /, *, str_as_lit: bool
    ) -> CompliantExprT | NonNestedLiteral:
        if is_expr(data):
            expr = data._to_compliant_expr(self)
            assert isinstance(expr, self._expr)  # noqa: S101
            return expr
        if isinstance(data, str) and not str_as_lit:
            return self.col(data)
        return data

    # NOTE: `polars`
    def all(self) -> CompliantExprT:
        return self._expr.from_column_names(get_column_names, context=self)

    def col(self, *column_names: str) -> CompliantExprT:
        return self._expr.from_column_names(
            passthrough_column_names(column_names), context=self
        )

    def exclude(self, excluded_names: Container[str]) -> CompliantExprT:
        return self._expr.from_column_names(
            partial(exclude_column_names, names=excluded_names), context=self
        )

    def nth(self, *column_indices: int) -> CompliantExprT:
        return self._expr.from_column_indices(*column_indices, context=self)

    def len(self) -> CompliantExprT: ...
    def lit(self, value: NonNestedLiteral, dtype: IntoDType | None) -> CompliantExprT: ...
    def all_horizontal(
        self, *exprs: CompliantExprT, ignore_nulls: bool
    ) -> CompliantExprT: ...
    def any_horizontal(
        self, *exprs: CompliantExprT, ignore_nulls: bool
    ) -> CompliantExprT: ...
    def sum_horizontal(self, *exprs: CompliantExprT) -> CompliantExprT: ...
    def mean_horizontal(self, *exprs: CompliantExprT) -> CompliantExprT: ...
    def min_horizontal(self, *exprs: CompliantExprT) -> CompliantExprT: ...
    def max_horizontal(self, *exprs: CompliantExprT) -> CompliantExprT: ...
    def concat(
        self, items: Iterable[CompliantFrameT], *, how: ConcatMethod
    ) -> CompliantFrameT: ...
    def when(
        self, predicate: CompliantExprT
    ) -> CompliantWhen[CompliantFrameT, Incomplete, CompliantExprT]: ...
    def concat_str(
        self, *exprs: CompliantExprT, separator: str, ignore_nulls: bool
    ) -> CompliantExprT: ...
    @property
    def selectors(self) -> CompliantSelectorNamespace[Any, Any]: ...
    def coalesce(self, *exprs: CompliantExprT) -> CompliantExprT: ...
    # NOTE: typing this accurately requires 2x more `TypeVar`s
    def from_native(self, data: Any, /) -> Any: ...
    def is_native(self, obj: Any, /) -> TypeIs[Any]:
        """Return `True` if `obj` can be passed to `from_native`."""
        ...


class DepthTrackingNamespace(
    CompliantNamespace[CompliantFrameT, DepthTrackingExprT],
    Protocol[CompliantFrameT, DepthTrackingExprT],
):
    def all(self) -> DepthTrackingExprT:
        return self._expr.from_column_names(
            get_column_names, function_name="all", context=self
        )

    def col(self, *column_names: str) -> DepthTrackingExprT:
        return self._expr.from_column_names(
            passthrough_column_names(column_names), function_name="col", context=self
        )

    def exclude(self, excluded_names: Container[str]) -> DepthTrackingExprT:
        return self._expr.from_column_names(
            partial(exclude_column_names, names=excluded_names),
            function_name="exclude",
            context=self,
        )


class LazyNamespace(
    CompliantNamespace[CompliantLazyFrameT, LazyExprT],
    Protocol[CompliantLazyFrameT, LazyExprT, NativeFrameT],
):
    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @property
    def _lazyframe(self) -> type[CompliantLazyFrameT]: ...
    def is_native(self, obj: Any, /) -> TypeIs[NativeFrameT]:
        return self._lazyframe._is_native(obj)

    def from_native(self, data: NativeFrameT | Any, /) -> CompliantLazyFrameT:
        if self._lazyframe._is_native(data):
            return self._lazyframe.from_native(data, context=self)
        msg = f"Unsupported type: {type(data).__name__!r}"  # pragma: no cover
        raise TypeError(msg)


class EagerNamespace(
    DepthTrackingNamespace[EagerDataFrameT, EagerExprT],
    Protocol[EagerDataFrameT, EagerSeriesT, EagerExprT, NativeFrameT, NativeSeriesT],
):
    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @property
    def _dataframe(self) -> type[EagerDataFrameT]: ...
    @property
    def _series(self) -> type[EagerSeriesT]: ...
    def when(
        self, predicate: EagerExprT
    ) -> EagerWhen[EagerDataFrameT, EagerSeriesT, EagerExprT, NativeSeriesT]: ...

    def is_native(self, obj: Any, /) -> TypeIs[NativeFrameT | NativeSeriesT]:
        return self._dataframe._is_native(obj) or self._series._is_native(obj)

    @overload
    def from_native(self, data: NativeFrameT, /) -> EagerDataFrameT: ...
    @overload
    def from_native(self, data: NativeSeriesT, /) -> EagerSeriesT: ...
    def from_native(
        self, data: NativeFrameT | NativeSeriesT | Any, /
    ) -> EagerDataFrameT | EagerSeriesT:
        if self._dataframe._is_native(data):
            return self._dataframe.from_native(data, context=self)
        if self._series._is_native(data):
            return self._series.from_native(data, context=self)
        msg = f"Unsupported type: {type(data).__name__!r}"
        raise TypeError(msg)

    def parse_into_expr(
        self,
        data: Expr | Series[NativeSeriesT] | _1DArray | NonNestedLiteral,
        /,
        *,
        str_as_lit: bool,
    ) -> EagerExprT | NonNestedLiteral:
        if not (is_series(data) or is_numpy_array(data)):
            return super().parse_into_expr(data, str_as_lit=str_as_lit)
        return self._expr._from_series(
            data._compliant_series
            if is_series(data)
            else self._series.from_numpy(data, context=self)
        )

    @overload
    def from_numpy(self, data: Into1DArray, /, schema: None = ...) -> EagerSeriesT: ...

    @overload
    def from_numpy(
        self, data: _2DArray, /, schema: IntoSchema | Sequence[str] | None
    ) -> EagerDataFrameT: ...

    def from_numpy(
        self,
        data: Into1DArray | _2DArray,
        /,
        schema: IntoSchema | Sequence[str] | None = None,
    ) -> EagerDataFrameT | EagerSeriesT:
        if is_numpy_array_2d(data):
            return self._dataframe.from_numpy(data, schema=schema, context=self)
        return self._series.from_numpy(data, context=self)

    def _concat_diagonal(self, dfs: Sequence[NativeFrameT], /) -> NativeFrameT: ...
    def _concat_horizontal(
        self, dfs: Sequence[NativeFrameT | Any], /
    ) -> NativeFrameT: ...
    def _concat_vertical(self, dfs: Sequence[NativeFrameT], /) -> NativeFrameT: ...
    def concat(
        self, items: Iterable[EagerDataFrameT], *, how: ConcatMethod
    ) -> EagerDataFrameT:
        dfs = [item.native for item in items]
        if how == "horizontal":
            native = self._concat_horizontal(dfs)
        elif how == "vertical":
            native = self._concat_vertical(dfs)
        elif how == "diagonal":
            native = self._concat_diagonal(dfs)
        else:  # pragma: no cover
            raise NotImplementedError
        return self._dataframe.from_native(native, context=self)
