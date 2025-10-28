from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Literal, cast

from narwhals._expression_parsing import ExprKind, ExprMetadata
from narwhals._spark_like.expr_dt import SparkLikeExprDateTimeNamespace
from narwhals._spark_like.expr_list import SparkLikeExprListNamespace
from narwhals._spark_like.expr_str import SparkLikeExprStringNamespace
from narwhals._spark_like.expr_struct import SparkLikeExprStructNamespace
from narwhals._spark_like.utils import (
    import_functions,
    import_native_dtypes,
    import_window,
    narwhals_to_native_dtype,
    true_divide,
)
from narwhals._sql.expr import SQLExpr
from narwhals._utils import (
    Implementation,
    Version,
    extend_bool,
    not_implemented,
    zip_strict,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, Sequence

    from sqlframe.base.column import Column
    from sqlframe.base.window import Window, WindowSpec
    from typing_extensions import Self, TypeAlias

    from narwhals._compliant import WindowInputs
    from narwhals._compliant.typing import (
        AliasNames,
        EvalNames,
        EvalSeries,
        WindowFunction,
    )
    from narwhals._spark_like.dataframe import SparkLikeLazyFrame
    from narwhals._spark_like.namespace import SparkLikeNamespace
    from narwhals._utils import _LimitedContext
    from narwhals.typing import FillNullStrategy, IntoDType, NonNestedLiteral, RankMethod

    NativeRankMethod: TypeAlias = Literal["rank", "dense_rank", "row_number"]
    SparkWindowFunction = WindowFunction[SparkLikeLazyFrame, Column]
    SparkWindowInputs = WindowInputs[Column]


class SparkLikeExpr(SQLExpr["SparkLikeLazyFrame", "Column"]):
    def __init__(
        self,
        call: EvalSeries[SparkLikeLazyFrame, Column],
        window_function: SparkWindowFunction | None = None,
        *,
        evaluate_output_names: EvalNames[SparkLikeLazyFrame],
        alias_output_names: AliasNames | None,
        version: Version,
        implementation: Implementation,
    ) -> None:
        self._call = call
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._version = version
        self._implementation = implementation
        self._metadata: ExprMetadata | None = None
        self._window_function: SparkWindowFunction | None = window_function

    _REMAP_RANK_METHOD: ClassVar[Mapping[RankMethod, NativeRankMethod]] = {
        "min": "rank",
        "max": "rank",
        "average": "rank",
        "dense": "dense_rank",
        "ordinal": "row_number",
    }

    def _count_star(self) -> Column:
        return self._F.count("*")

    def _window_expression(
        self,
        expr: Column,
        partition_by: Sequence[str | Column] = (),
        order_by: Sequence[str | Column] = (),
        rows_start: int | None = None,
        rows_end: int | None = None,
        *,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> Column:
        window = self.partition_by(*partition_by)
        if order_by:
            window = window.orderBy(
                *self._sort(*order_by, descending=descending, nulls_last=nulls_last)
            )
        if rows_start is not None and rows_end is not None:
            window = window.rowsBetween(rows_start, rows_end)
        elif rows_end is not None:
            window = window.rowsBetween(self._Window.unboundedPreceding, rows_end)
        elif rows_start is not None:  # pragma: no cover
            window = window.rowsBetween(rows_start, self._Window.unboundedFollowing)
        return expr.over(window)

    def _first(self, expr: Column, *order_by: str) -> Column:
        # Docs say it's non-deterministic, with no way to specify order.
        msg = "`first` is not supported for PySpark."
        raise NotImplementedError(msg)

    def _last(self, expr: Column, *order_by: str) -> Column:  # pragma: no cover
        # Docs say it's non-deterministic, with no way to specify order.
        msg = "`last` is not supported for PySpark."
        raise NotImplementedError(msg)

    def broadcast(self, kind: Literal[ExprKind.AGGREGATION, ExprKind.LITERAL]) -> Self:
        if kind is ExprKind.LITERAL:
            return self
        return self.over([self._F.lit(1)], [])

    @property
    def _F(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import functions

            return functions
        return import_functions(self._implementation)

    @property
    def _native_dtypes(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import types

            return types
        return import_native_dtypes(self._implementation)

    @property
    def _Window(self) -> type[Window]:
        if TYPE_CHECKING:
            from sqlframe.base.window import Window

            return Window
        return import_window(self._implementation)

    def _sort(
        self,
        *cols: Column | str,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> Iterator[Column]:
        F = self._F
        n = len(cols)
        descending = extend_bool(descending or False, n)
        nulls_last = extend_bool(nulls_last or False, n)
        mapping = {
            (False, False): F.asc_nulls_first,
            (False, True): F.asc_nulls_last,
            (True, False): F.desc_nulls_first,
            (True, True): F.desc_nulls_last,
        }
        yield from (
            mapping[(_desc, _nulls_last)](col)
            for col, _desc, _nulls_last in zip_strict(cols, descending, nulls_last)
        )

    def partition_by(self, *cols: Column | str) -> WindowSpec:
        """Wraps `Window().partitionBy`, with default and `WindowInputs` handling."""
        return self._Window.partitionBy(*cols or [self._F.lit(1)])

    def __narwhals_namespace__(self) -> SparkLikeNamespace:  # pragma: no cover
        from narwhals._spark_like.namespace import SparkLikeNamespace

        return SparkLikeNamespace(
            version=self._version, implementation=self._implementation
        )

    @classmethod
    def _alias_native(cls, expr: Column, name: str) -> Column:
        return expr.alias(name)

    @classmethod
    def from_column_names(
        cls: type[Self],
        evaluate_column_names: EvalNames[SparkLikeLazyFrame],
        /,
        *,
        context: _LimitedContext,
    ) -> Self:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            return [df._F.col(col_name) for col_name in evaluate_column_names(df)]

        return cls(
            func,
            evaluate_output_names=evaluate_column_names,
            alias_output_names=None,
            version=context._version,
            implementation=context._implementation,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            columns = df.columns
            return [df._F.col(columns[i]) for i in column_indices]

        return cls(
            func,
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            version=context._version,
            implementation=context._implementation,
        )

    def __truediv__(self, other: SparkLikeExpr) -> Self:
        def _truediv(expr: Column, other: Column) -> Column:
            return true_divide(self._F, expr, other)

        return self._with_binary(_truediv, other)

    def __rtruediv__(self, other: SparkLikeExpr) -> Self:
        def _rtruediv(expr: Column, other: Column) -> Column:
            return true_divide(self._F, other, expr)

        return self._with_binary(_rtruediv, other).alias("literal")

    def __floordiv__(self, other: SparkLikeExpr) -> Self:
        def _floordiv(expr: Column, other: Column) -> Column:
            F = self._F
            return F.when(
                other != F.lit(0), F.floor(true_divide(F, expr, other))
            ).otherwise(F.lit(None))

        return self._with_binary(_floordiv, other)

    def __rfloordiv__(self, other: SparkLikeExpr) -> Self:
        def _rfloordiv(expr: Column, other: Column) -> Column:
            F = self._F
            return F.when(
                expr != F.lit(0), F.floor(true_divide(F, other, expr))
            ).otherwise(F.lit(None))

        return self._with_binary(_rfloordiv, other).alias("literal")

    def __invert__(self) -> Self:
        invert = cast("Callable[..., Column]", operator.invert)
        return self._with_elementwise(invert)

    def cast(self, dtype: IntoDType) -> Self:
        def func(df: SparkLikeLazyFrame) -> Sequence[Column]:
            spark_dtype = narwhals_to_native_dtype(
                dtype, self._version, self._native_dtypes, df.native.sparkSession
            )
            return [expr.cast(spark_dtype) for expr in self(df)]

        def window_f(
            df: SparkLikeLazyFrame, inputs: SparkWindowInputs
        ) -> Sequence[Column]:
            spark_dtype = narwhals_to_native_dtype(
                dtype, self._version, self._native_dtypes, df.native.sparkSession
            )
            return [expr.cast(spark_dtype) for expr in self.window_function(df, inputs)]

        return self.__class__(
            func,
            window_f,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
            implementation=self._implementation,
        )

    def median(self) -> Self:
        def _median(expr: Column) -> Column:
            if self._implementation in {
                Implementation.PYSPARK,
                Implementation.PYSPARK_CONNECT,
            } and Implementation.PYSPARK._backend_version() < (3, 4):  # pragma: no cover
                # Use percentile_approx with default accuracy parameter (10000)
                return self._F.percentile_approx(expr.cast("double"), 0.5)

            return self._F.median(expr)

        return self._with_callable(_median)

    def null_count(self) -> Self:
        def _null_count(expr: Column) -> Column:
            return self._F.count_if(self._F.isnull(expr))

        return self._with_callable(_null_count)

    def std(self, *, ddof: int) -> Self:
        F = self._F
        if ddof == 0:
            return self._with_callable(F.stddev_pop)
        if ddof == 1:
            return self._with_callable(F.stddev_samp)

        def func(expr: Column) -> Column:
            n_rows = F.count(expr)
            return F.stddev_samp(expr) * F.sqrt((n_rows - 1) / (n_rows - ddof))

        return self._with_callable(func)

    def var(self, *, ddof: int) -> Self:
        F = self._F
        if ddof == 0:
            return self._with_callable(F.var_pop)
        if ddof == 1:
            return self._with_callable(F.var_samp)

        def func(expr: Column) -> Column:
            n_rows = F.count(expr)
            return F.var_samp(expr) * (n_rows - 1) / (n_rows - ddof)

        return self._with_callable(func)

    def is_finite(self) -> Self:
        def _is_finite(expr: Column) -> Column:
            # A value is finite if it's not NaN, and not infinite, while NULLs should be
            # preserved
            is_finite_condition = (
                ~self._F.isnan(expr)
                & (expr != self._F.lit(float("inf")))
                & (expr != self._F.lit(float("-inf")))
            )
            return self._F.when(~self._F.isnull(expr), is_finite_condition).otherwise(
                None
            )

        return self._with_elementwise(_is_finite)

    def is_in(self, other: Sequence[Any]) -> Self:
        def _is_in(expr: Column) -> Column:
            return expr.isin(other) if other else self._F.lit(False)

        return self._with_elementwise(_is_in)

    def len(self) -> Self:
        def _len(_expr: Column) -> Column:
            # Use count(*) to count all rows including nulls
            return self._F.count("*")

        return self._with_callable(_len)

    def skew(self) -> Self:
        return self._with_callable(self._F.skewness)

    def kurtosis(self) -> Self:
        return self._with_callable(self._F.kurtosis)

    def is_nan(self) -> Self:
        def _is_nan(expr: Column) -> Column:
            return self._F.when(self._F.isnull(expr), None).otherwise(self._F.isnan(expr))

        return self._with_elementwise(_is_nan)

    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self:
        if strategy is not None:

            def _fill_with_strategy(
                df: SparkLikeLazyFrame, inputs: SparkWindowInputs
            ) -> Sequence[Column]:
                fn = self._F.last_value if strategy == "forward" else self._F.first_value
                if strategy == "forward":
                    start = self._Window.unboundedPreceding if limit is None else -limit
                    end = self._Window.currentRow
                else:
                    start = self._Window.currentRow
                    end = self._Window.unboundedFollowing if limit is None else limit
                return [
                    fn(expr, ignoreNulls=True).over(
                        self.partition_by(*inputs.partition_by)
                        .orderBy(*self._sort(*inputs.order_by))
                        .rowsBetween(start, end)
                    )
                    for expr in self(df)
                ]

            return self._with_window_function(_fill_with_strategy)

        def _fill_constant(expr: Column, value: Column) -> Column:
            return self._F.ifnull(expr, value)

        return self._with_elementwise(_fill_constant, value=value)

    @property
    def str(self) -> SparkLikeExprStringNamespace:
        return SparkLikeExprStringNamespace(self)

    @property
    def dt(self) -> SparkLikeExprDateTimeNamespace:
        return SparkLikeExprDateTimeNamespace(self)

    @property
    def list(self) -> SparkLikeExprListNamespace:
        return SparkLikeExprListNamespace(self)

    @property
    def struct(self) -> SparkLikeExprStructNamespace:
        return SparkLikeExprStructNamespace(self)

    quantile = not_implemented()
