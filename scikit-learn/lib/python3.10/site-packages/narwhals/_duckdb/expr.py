from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Callable, Literal, cast

from duckdb import CoalesceOperator, StarExpression

from narwhals._duckdb.expr_dt import DuckDBExprDateTimeNamespace
from narwhals._duckdb.expr_list import DuckDBExprListNamespace
from narwhals._duckdb.expr_str import DuckDBExprStringNamespace
from narwhals._duckdb.expr_struct import DuckDBExprStructNamespace
from narwhals._duckdb.utils import (
    DeferredTimeZone,
    F,
    col,
    generate_order_by_sql,
    lit,
    narwhals_to_native_dtype,
    sql_expression,
    when,
    window_expression,
)
from narwhals._expression_parsing import ExprKind, ExprMetadata
from narwhals._sql.expr import SQLExpr
from narwhals._utils import Implementation, Version, extend_bool

if TYPE_CHECKING:
    from collections.abc import Sequence

    from duckdb import Expression
    from typing_extensions import Self

    from narwhals._compliant import WindowInputs
    from narwhals._compliant.typing import (
        AliasNames,
        EvalNames,
        EvalSeries,
        WindowFunction,
    )
    from narwhals._duckdb.dataframe import DuckDBLazyFrame
    from narwhals._duckdb.namespace import DuckDBNamespace
    from narwhals._utils import _LimitedContext
    from narwhals.typing import (
        FillNullStrategy,
        IntoDType,
        NonNestedLiteral,
        RollingInterpolationMethod,
    )

    DuckDBWindowFunction = WindowFunction[DuckDBLazyFrame, Expression]
    DuckDBWindowInputs = WindowInputs[Expression]


class DuckDBExpr(SQLExpr["DuckDBLazyFrame", "Expression"]):
    _implementation = Implementation.DUCKDB

    def __init__(
        self,
        call: EvalSeries[DuckDBLazyFrame, Expression],
        window_function: DuckDBWindowFunction | None = None,
        *,
        evaluate_output_names: EvalNames[DuckDBLazyFrame],
        alias_output_names: AliasNames | None,
        version: Version,
        implementation: Implementation = Implementation.DUCKDB,
    ) -> None:
        self._call = call
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._version = version
        self._metadata: ExprMetadata | None = None
        self._window_function: DuckDBWindowFunction | None = window_function

    def _count_star(self) -> Expression:
        return F("count", StarExpression())

    def _window_expression(
        self,
        expr: Expression,
        partition_by: Sequence[str | Expression] = (),
        order_by: Sequence[str | Expression] = (),
        rows_start: int | None = None,
        rows_end: int | None = None,
        *,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> Expression:
        return window_expression(
            expr,
            partition_by,
            order_by,
            rows_start,
            rows_end,
            descending=descending,
            nulls_last=nulls_last,
        )

    def _first_last(
        self, function: str, expr: Expression, order_by: Sequence[str], /
    ) -> Expression:
        # https://github.com/duckdb/duckdb/discussions/19252
        flags = extend_bool(False, len(order_by))
        order_by_sql = generate_order_by_sql(
            *order_by, descending=flags, nulls_last=flags
        )
        return sql_expression(f"{function}({expr} {order_by_sql})")

    def _first(self, expr: Expression, *order_by: str) -> Expression:
        return self._first_last("first", expr, order_by)

    def _last(self, expr: Expression, *order_by: str) -> Expression:
        return self._first_last("last", expr, order_by)

    def __narwhals_namespace__(self) -> DuckDBNamespace:  # pragma: no cover
        from narwhals._duckdb.namespace import DuckDBNamespace

        return DuckDBNamespace(version=self._version)

    def broadcast(self, kind: Literal[ExprKind.AGGREGATION, ExprKind.LITERAL]) -> Self:
        if kind is ExprKind.LITERAL:
            return self
        if self._backend_version < (1, 3):
            msg = "At least version 1.3 of DuckDB is required for binary operations between aggregates and columns."
            raise NotImplementedError(msg)
        return self.over([lit(1)], [])

    @classmethod
    def from_column_names(
        cls,
        evaluate_column_names: EvalNames[DuckDBLazyFrame],
        /,
        *,
        context: _LimitedContext,
    ) -> Self:
        def func(df: DuckDBLazyFrame) -> list[Expression]:
            return [col(name) for name in evaluate_column_names(df)]

        return cls(
            func,
            evaluate_output_names=evaluate_column_names,
            alias_output_names=None,
            version=context._version,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: DuckDBLazyFrame) -> list[Expression]:
            columns = df.columns
            return [col(columns[i]) for i in column_indices]

        return cls(
            func,
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            version=context._version,
        )

    @classmethod
    def _alias_native(cls, expr: Expression, name: str) -> Expression:
        return expr.alias(name)

    def __invert__(self) -> Self:
        invert = cast("Callable[..., Expression]", operator.invert)
        return self._with_elementwise(invert)

    def skew(self) -> Self:
        W = self._window_expression  # noqa: N806

        def func(expr: Expression) -> Expression:
            count = F("count", expr)
            # Adjust population skewness by correction factor to get sample skewness
            sample_skewness = (
                F("skewness", expr)
                * (count - lit(2))
                / F("sqrt", count * (count - lit(1)))
            )
            return when(count == lit(0), lit(None)).otherwise(
                when(count == lit(1), lit(float("nan"))).otherwise(
                    when(count == lit(2), lit(0.0)).otherwise(sample_skewness)
                )
            )

        def window_f(df: DuckDBLazyFrame, inputs: DuckDBWindowInputs) -> list[Expression]:
            ret = []
            for expr in self(df):
                count = W(F("count", expr), inputs.partition_by)
                # Adjust population skewness by correction factor to get sample skewness
                sample_skewness = (
                    W(F("skewness", expr), inputs.partition_by)
                    * (count - lit(2))
                    / F("sqrt", count * (count - lit(1)))
                )
                ret.append(
                    when(count == lit(0), lit(None)).otherwise(
                        when(count == lit(1), lit(float("nan"))).otherwise(
                            when(count == lit(2), lit(0.0)).otherwise(sample_skewness)
                        )
                    )
                )
            return ret

        return self._with_callable(func, window_f)

    def kurtosis(self) -> Self:
        return self._with_callable(lambda expr: F("kurtosis_pop", expr))

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> Self:
        def func(expr: Expression) -> Expression:
            if interpolation == "linear":
                return F("quantile_cont", expr, lit(quantile))
            msg = "Only linear interpolation methods are supported for DuckDB quantile."
            raise NotImplementedError(msg)

        return self._with_callable(func)

    def len(self) -> Self:
        return self._with_callable(lambda _expr: F("count"))

    def std(self, *, ddof: int) -> Self:
        if ddof == 0:
            return self._with_callable(lambda expr: F("stddev_pop", expr))
        if ddof == 1:
            return self._with_callable(lambda expr: F("stddev_samp", expr))

        def _std(expr: Expression) -> Expression:
            n_samples = F("count", expr)
            return (
                F("stddev_pop", expr)
                * F("sqrt", n_samples)
                / (F("sqrt", (n_samples - lit(ddof))))
            )

        return self._with_callable(_std)

    def var(self, *, ddof: int) -> Self:
        if ddof == 0:
            return self._with_callable(lambda expr: F("var_pop", expr))
        if ddof == 1:
            return self._with_callable(lambda expr: F("var_samp", expr))

        def _var(expr: Expression) -> Expression:
            n_samples = F("count", expr)
            return F("var_pop", expr) * n_samples / (n_samples - lit(ddof))

        return self._with_callable(_var)

    def null_count(self) -> Self:
        return self._with_callable(lambda expr: F("sum", expr.isnull().cast("int")))

    def is_nan(self) -> Self:
        return self._with_elementwise(lambda expr: F("isnan", expr))

    def is_finite(self) -> Self:
        return self._with_elementwise(lambda expr: F("isfinite", expr))

    def is_in(self, other: Sequence[Any]) -> Self:
        return self._with_elementwise(lambda expr: F("contains", lit(other), expr))

    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self:
        if strategy is not None:
            if self._backend_version < (1, 3):  # pragma: no cover
                msg = f"`fill_null` with `strategy={strategy}` is only available in 'duckdb>=1.3.0'."
                raise NotImplementedError(msg)

            def _fill_with_strategy(
                df: DuckDBLazyFrame, inputs: DuckDBWindowInputs
            ) -> Sequence[Expression]:
                fill_func = "last_value" if strategy == "forward" else "first_value"
                rows_start, rows_end = (
                    (-limit if limit is not None else None, 0)
                    if strategy == "forward"
                    else (0, limit)
                )
                return [
                    window_expression(
                        F(fill_func, expr),
                        inputs.partition_by,
                        inputs.order_by,
                        rows_start=rows_start,
                        rows_end=rows_end,
                        ignore_nulls=True,
                    )
                    for expr in self(df)
                ]

            return self._with_window_function(_fill_with_strategy)

        def _fill_constant(expr: Expression, value: Any) -> Expression:
            return CoalesceOperator(expr, value)

        return self._with_elementwise(_fill_constant, value=value)

    def cast(self, dtype: IntoDType) -> Self:
        def func(df: DuckDBLazyFrame) -> list[Expression]:
            tz = DeferredTimeZone(df.native)
            native_dtype = narwhals_to_native_dtype(dtype, self._version, tz)
            return [expr.cast(native_dtype) for expr in self(df)]

        def window_f(df: DuckDBLazyFrame, inputs: DuckDBWindowInputs) -> list[Expression]:
            tz = DeferredTimeZone(df.native)
            native_dtype = narwhals_to_native_dtype(dtype, self._version, tz)
            return [expr.cast(native_dtype) for expr in self.window_function(df, inputs)]

        return self.__class__(
            func,
            window_f,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    @property
    def str(self) -> DuckDBExprStringNamespace:
        return DuckDBExprStringNamespace(self)

    @property
    def dt(self) -> DuckDBExprDateTimeNamespace:
        return DuckDBExprDateTimeNamespace(self)

    @property
    def list(self) -> DuckDBExprListNamespace:
        return DuckDBExprListNamespace(self)

    @property
    def struct(self) -> DuckDBExprStructNamespace:
        return DuckDBExprStructNamespace(self)
