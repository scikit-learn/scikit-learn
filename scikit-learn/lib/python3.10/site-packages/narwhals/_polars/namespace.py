from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import polars as pl

from narwhals._expression_parsing import is_expr, is_series
from narwhals._polars.expr import PolarsExpr
from narwhals._polars.series import PolarsSeries
from narwhals._polars.utils import extract_args_kwargs, narwhals_to_native_dtype
from narwhals._utils import Implementation, requires, zip_strict
from narwhals.dependencies import is_numpy_array_2d
from narwhals.dtypes import DType

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from datetime import timezone

    from typing_extensions import TypeIs

    from narwhals._compliant import CompliantSelectorNamespace, CompliantWhen
    from narwhals._polars.dataframe import Method, PolarsDataFrame, PolarsLazyFrame
    from narwhals._polars.typing import FrameT
    from narwhals._utils import Version, _LimitedContext
    from narwhals.expr import Expr
    from narwhals.series import Series
    from narwhals.typing import (
        Into1DArray,
        IntoDType,
        IntoSchema,
        NonNestedLiteral,
        TimeUnit,
        _1DArray,
        _2DArray,
    )


class PolarsNamespace:
    all: Method[PolarsExpr]
    coalesce: Method[PolarsExpr]
    col: Method[PolarsExpr]
    exclude: Method[PolarsExpr]
    sum_horizontal: Method[PolarsExpr]
    min_horizontal: Method[PolarsExpr]
    max_horizontal: Method[PolarsExpr]

    when: Method[CompliantWhen[PolarsDataFrame, PolarsSeries, PolarsExpr]]

    _implementation: Implementation = Implementation.POLARS
    _version: Version

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def __getattr__(self, attr: str) -> Any:
        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            return self._expr(getattr(pl, attr)(*pos, **kwds), version=self._version)

        return func

    @property
    def _dataframe(self) -> type[PolarsDataFrame]:
        from narwhals._polars.dataframe import PolarsDataFrame

        return PolarsDataFrame

    @property
    def _lazyframe(self) -> type[PolarsLazyFrame]:
        from narwhals._polars.dataframe import PolarsLazyFrame

        return PolarsLazyFrame

    @property
    def _expr(self) -> type[PolarsExpr]:
        return PolarsExpr

    @property
    def _series(self) -> type[PolarsSeries]:
        return PolarsSeries

    def parse_into_expr(
        self,
        data: Expr | NonNestedLiteral | Series[pl.Series] | _1DArray,
        /,
        *,
        str_as_lit: bool,
    ) -> PolarsExpr | None:
        if data is None:
            # NOTE: To avoid `pl.lit(None)` failing this `None` check
            # https://github.com/pola-rs/polars/blob/58dd8e5770f16a9bef9009a1c05f00e15a5263c7/py-polars/polars/expr/expr.py#L2870-L2872
            return data
        if is_expr(data):
            expr = data._to_compliant_expr(self)
            assert isinstance(expr, self._expr)  # noqa: S101
            return expr
        if isinstance(data, str) and not str_as_lit:
            return self.col(data)
        return self.lit(data.to_native() if is_series(data) else data, None)

    def is_native(self, obj: Any) -> TypeIs[pl.DataFrame | pl.LazyFrame | pl.Series]:
        return isinstance(obj, (pl.DataFrame, pl.LazyFrame, pl.Series))

    @overload
    def from_native(self, data: pl.DataFrame, /) -> PolarsDataFrame: ...
    @overload
    def from_native(self, data: pl.LazyFrame, /) -> PolarsLazyFrame: ...
    @overload
    def from_native(self, data: pl.Series, /) -> PolarsSeries: ...
    def from_native(
        self, data: pl.DataFrame | pl.LazyFrame | pl.Series | Any, /
    ) -> PolarsDataFrame | PolarsLazyFrame | PolarsSeries:
        if self._dataframe._is_native(data):
            return self._dataframe.from_native(data, context=self)
        if self._series._is_native(data):
            return self._series.from_native(data, context=self)
        if self._lazyframe._is_native(data):
            return self._lazyframe.from_native(data, context=self)
        msg = f"Unsupported type: {type(data).__name__!r}"  # pragma: no cover
        raise TypeError(msg)  # pragma: no cover

    @overload
    def from_numpy(self, data: Into1DArray, /, schema: None = ...) -> PolarsSeries: ...

    @overload
    def from_numpy(
        self, data: _2DArray, /, schema: IntoSchema | Sequence[str] | None
    ) -> PolarsDataFrame: ...

    def from_numpy(
        self,
        data: Into1DArray | _2DArray,
        /,
        schema: IntoSchema | Sequence[str] | None = None,
    ) -> PolarsDataFrame | PolarsSeries:
        if is_numpy_array_2d(data):
            return self._dataframe.from_numpy(data, schema=schema, context=self)
        return self._series.from_numpy(data, context=self)  # pragma: no cover

    @requires.backend_version(
        (1, 0, 0), "Please use `col` for columns selection instead."
    )
    def nth(self, *indices: int) -> PolarsExpr:
        return self._expr(pl.nth(*indices), version=self._version)

    def len(self) -> PolarsExpr:
        if self._backend_version < (0, 20, 5):
            return self._expr(pl.count().alias("len"), self._version)
        return self._expr(pl.len(), self._version)

    def all_horizontal(self, *exprs: PolarsExpr, ignore_nulls: bool) -> PolarsExpr:
        it = (expr.fill_null(True) for expr in exprs) if ignore_nulls else iter(exprs)
        return self._expr(pl.all_horizontal(*(expr.native for expr in it)), self._version)

    def any_horizontal(self, *exprs: PolarsExpr, ignore_nulls: bool) -> PolarsExpr:
        it = (expr.fill_null(False) for expr in exprs) if ignore_nulls else iter(exprs)
        return self._expr(pl.any_horizontal(*(expr.native for expr in it)), self._version)

    def concat(
        self,
        items: Iterable[FrameT],
        *,
        how: Literal["vertical", "horizontal", "diagonal"],
    ) -> PolarsDataFrame | PolarsLazyFrame:
        result = pl.concat((item.native for item in items), how=how)
        if isinstance(result, pl.DataFrame):
            return self._dataframe(result, version=self._version)
        return self._lazyframe.from_native(result, context=self)

    def lit(self, value: Any, dtype: IntoDType | None) -> PolarsExpr:
        if dtype is not None:
            return self._expr(
                pl.lit(value, dtype=narwhals_to_native_dtype(dtype, self._version)),
                version=self._version,
            )
        return self._expr(pl.lit(value), version=self._version)

    def mean_horizontal(self, *exprs: PolarsExpr) -> PolarsExpr:
        if self._backend_version < (0, 20, 8):
            return self._expr(
                pl.sum_horizontal(e._native_expr for e in exprs)
                / pl.sum_horizontal(1 - e.is_null()._native_expr for e in exprs),
                version=self._version,
            )

        return self._expr(
            pl.mean_horizontal(e._native_expr for e in exprs), version=self._version
        )

    def concat_str(
        self, *exprs: PolarsExpr, separator: str, ignore_nulls: bool
    ) -> PolarsExpr:
        pl_exprs: list[pl.Expr] = [expr._native_expr for expr in exprs]

        if self._backend_version < (0, 20, 6):
            null_mask = [expr.is_null() for expr in pl_exprs]
            sep = pl.lit(separator)

            if not ignore_nulls:
                null_mask_result = pl.any_horizontal(*null_mask)
                output_expr = pl.reduce(
                    lambda x, y: x.cast(pl.String()) + sep + y.cast(pl.String()),  # type: ignore[arg-type,return-value]
                    pl_exprs,
                )
                result = pl.when(~null_mask_result).then(output_expr)
            else:
                init_value, *values = [
                    pl.when(nm).then(pl.lit("")).otherwise(expr.cast(pl.String()))
                    for expr, nm in zip_strict(pl_exprs, null_mask)
                ]
                separators = [
                    pl.when(~nm).then(sep).otherwise(pl.lit("")) for nm in null_mask[:-1]
                ]

                result = pl.fold(  # type: ignore[assignment]
                    acc=init_value,
                    function=operator.add,
                    exprs=[s + v for s, v in zip_strict(separators, values)],
                )

            return self._expr(result, version=self._version)

        return self._expr(
            pl.concat_str(pl_exprs, separator=separator, ignore_nulls=ignore_nulls),
            version=self._version,
        )

    # NOTE: Implementation is too different to annotate correctly (vs other `*SelectorNamespace`)
    # 1. Others have lots of private stuff for code reuse
    #    i. None of that is useful here
    # 2. We don't have a `PolarsSelector` abstraction, and just use `PolarsExpr`
    @property
    def selectors(self) -> CompliantSelectorNamespace[PolarsDataFrame, PolarsSeries]:
        return cast(
            "CompliantSelectorNamespace[PolarsDataFrame, PolarsSeries]",
            PolarsSelectorNamespace(self),
        )


class PolarsSelectorNamespace:
    _implementation = Implementation.POLARS

    def __init__(self, context: _LimitedContext, /) -> None:
        self._version = context._version

    def by_dtype(self, dtypes: Iterable[DType]) -> PolarsExpr:
        native_dtypes = [
            narwhals_to_native_dtype(dtype, self._version).__class__
            if isinstance(dtype, type) and issubclass(dtype, DType)
            else narwhals_to_native_dtype(dtype, self._version)
            for dtype in dtypes
        ]
        return PolarsExpr(pl.selectors.by_dtype(native_dtypes), version=self._version)

    def matches(self, pattern: str) -> PolarsExpr:
        return PolarsExpr(pl.selectors.matches(pattern=pattern), version=self._version)

    def numeric(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.numeric(), version=self._version)

    def boolean(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.boolean(), version=self._version)

    def string(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.string(), version=self._version)

    def categorical(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.categorical(), version=self._version)

    def all(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.all(), version=self._version)

    def datetime(
        self,
        time_unit: TimeUnit | Iterable[TimeUnit] | None,
        time_zone: str | timezone | Iterable[str | timezone | None] | None,
    ) -> PolarsExpr:
        return PolarsExpr(
            pl.selectors.datetime(time_unit=time_unit, time_zone=time_zone),  # type: ignore[arg-type]
            version=self._version,
        )
