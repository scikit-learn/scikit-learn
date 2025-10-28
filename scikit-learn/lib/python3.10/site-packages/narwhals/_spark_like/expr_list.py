from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace

if TYPE_CHECKING:
    from sqlframe.base.column import Column

    from narwhals._spark_like.expr import SparkLikeExpr
    from narwhals.typing import NonNestedLiteral


class SparkLikeExprListNamespace(
    LazyExprNamespace["SparkLikeExpr"], ListNamespace["SparkLikeExpr"]
):
    def len(self) -> SparkLikeExpr:
        return self.compliant._with_elementwise(self.compliant._F.array_size)

    def unique(self) -> SparkLikeExpr:
        return self.compliant._with_elementwise(self.compliant._F.array_distinct)

    def contains(self, item: NonNestedLiteral) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            return F.array_contains(expr, F.lit(item))

        return self.compliant._with_elementwise(func)

    def get(self, index: int) -> SparkLikeExpr:
        def _get(expr: Column) -> Column:
            return expr.getItem(index)

        return self.compliant._with_elementwise(_get)
