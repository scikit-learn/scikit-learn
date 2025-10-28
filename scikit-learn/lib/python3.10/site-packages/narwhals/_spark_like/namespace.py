from __future__ import annotations

import operator
from functools import reduce
from typing import TYPE_CHECKING, Any

from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._spark_like.dataframe import SparkLikeLazyFrame
from narwhals._spark_like.expr import SparkLikeExpr
from narwhals._spark_like.selectors import SparkLikeSelectorNamespace
from narwhals._spark_like.utils import (
    import_functions,
    import_native_dtypes,
    narwhals_to_native_dtype,
    true_divide,
)
from narwhals._sql.namespace import SQLNamespace
from narwhals._sql.when_then import SQLThen, SQLWhen
from narwhals._utils import zip_strict

if TYPE_CHECKING:
    from collections.abc import Iterable

    from sqlframe.base.column import Column

    from narwhals._spark_like.dataframe import SQLFrameDataFrame  # noqa: F401
    from narwhals._utils import Implementation, Version
    from narwhals.typing import ConcatMethod, IntoDType, NonNestedLiteral, PythonLiteral

# Adjust slight SQL vs PySpark differences
FUNCTION_REMAPPINGS = {
    "starts_with": "startswith",
    "ends_with": "endswith",
    "trim": "btrim",
    "str_split": "split",
    "regexp_matches": "regexp",
}


class SparkLikeNamespace(
    SQLNamespace[SparkLikeLazyFrame, SparkLikeExpr, "SQLFrameDataFrame", "Column"]
):
    def __init__(self, *, version: Version, implementation: Implementation) -> None:
        self._version = version
        self._implementation = implementation

    @property
    def selectors(self) -> SparkLikeSelectorNamespace:
        return SparkLikeSelectorNamespace.from_namespace(self)

    @property
    def _expr(self) -> type[SparkLikeExpr]:
        return SparkLikeExpr

    @property
    def _lazyframe(self) -> type[SparkLikeLazyFrame]:
        return SparkLikeLazyFrame

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

    def _function(self, name: str, *args: Column | PythonLiteral) -> Column:
        return getattr(self._F, FUNCTION_REMAPPINGS.get(name, name))(*args)

    def _lit(self, value: Any) -> Column:
        return self._F.lit(value)

    def _when(
        self, condition: Column, value: Column, otherwise: Column | None = None
    ) -> Column:
        if otherwise is None:
            return self._F.when(condition, value)
        return self._F.when(condition, value).otherwise(otherwise)

    def _coalesce(self, *exprs: Column) -> Column:
        return self._F.coalesce(*exprs)

    def lit(self, value: NonNestedLiteral, dtype: IntoDType | None) -> SparkLikeExpr:
        def _lit(df: SparkLikeLazyFrame) -> list[Column]:
            column = df._F.lit(value)
            if dtype:
                native_dtype = narwhals_to_native_dtype(
                    dtype, self._version, df._native_dtypes, df.native.sparkSession
                )
                column = column.cast(native_dtype)

            return [column]

        return self._expr(
            call=_lit,
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
            implementation=self._implementation,
        )

    def len(self) -> SparkLikeExpr:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            return [df._F.count("*")]

        return self._expr(
            func,
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
            implementation=self._implementation,
        )

    def mean_horizontal(self, *exprs: SparkLikeExpr) -> SparkLikeExpr:
        def func(cols: Iterable[Column]) -> Column:
            cols = tuple(cols)
            F = exprs[0]._F
            numerator = reduce(
                operator.add, (self._F.coalesce(col, self._F.lit(0)) for col in cols)
            )
            denominator = reduce(
                operator.add,
                (col.isNotNull().cast(self._native_dtypes.IntegerType()) for col in cols),
            )
            return true_divide(F, numerator, denominator)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def concat(
        self, items: Iterable[SparkLikeLazyFrame], *, how: ConcatMethod
    ) -> SparkLikeLazyFrame:
        dfs = [item._native_frame for item in items]
        if how == "vertical":
            cols_0 = dfs[0].columns
            for i, df in enumerate(dfs[1:], start=1):
                cols_current = df.columns
                if not ((len(cols_current) == len(cols_0)) and (cols_current == cols_0)):
                    msg = (
                        "unable to vstack, column names don't match:\n"
                        f"   - dataframe 0: {cols_0}\n"
                        f"   - dataframe {i}: {cols_current}\n"
                    )
                    raise TypeError(msg)

            return SparkLikeLazyFrame(
                native_dataframe=reduce(lambda x, y: x.union(y), dfs),
                version=self._version,
                implementation=self._implementation,
            )

        if how == "diagonal":
            return SparkLikeLazyFrame(
                native_dataframe=reduce(
                    lambda x, y: x.unionByName(y, allowMissingColumns=True), dfs
                ),
                version=self._version,
                implementation=self._implementation,
            )
        raise NotImplementedError

    def concat_str(
        self, *exprs: SparkLikeExpr, separator: str, ignore_nulls: bool
    ) -> SparkLikeExpr:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            cols = [s for _expr in exprs for s in _expr(df)]
            cols_casted = [s.cast(df._native_dtypes.StringType()) for s in cols]
            null_mask = [df._F.isnull(s) for s in cols]

            if not ignore_nulls:
                null_mask_result = reduce(operator.or_, null_mask)
                result = df._F.when(
                    ~null_mask_result,
                    reduce(
                        lambda x, y: df._F.format_string(f"%s{separator}%s", x, y),
                        cols_casted,
                    ),
                ).otherwise(df._F.lit(None))
            else:
                init_value, *values = [
                    df._F.when(~nm, col).otherwise(df._F.lit(""))
                    for col, nm in zip_strict(cols_casted, null_mask)
                ]

                separators = (
                    df._F.when(nm, df._F.lit("")).otherwise(df._F.lit(separator))
                    for nm in null_mask[:-1]
                )
                result = reduce(
                    lambda x, y: df._F.format_string("%s%s", x, y),
                    (
                        df._F.format_string("%s%s", s, v)
                        for s, v in zip_strict(separators, values)
                    ),
                    init_value,
                )

            return [result]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
            implementation=self._implementation,
        )

    def when(self, predicate: SparkLikeExpr) -> SparkLikeWhen:
        return SparkLikeWhen.from_expr(predicate, context=self)


class SparkLikeWhen(SQLWhen[SparkLikeLazyFrame, "Column", SparkLikeExpr]):
    @property
    def _then(self) -> type[SparkLikeThen]:
        return SparkLikeThen


class SparkLikeThen(
    SQLThen[SparkLikeLazyFrame, "Column", SparkLikeExpr], SparkLikeExpr
): ...
