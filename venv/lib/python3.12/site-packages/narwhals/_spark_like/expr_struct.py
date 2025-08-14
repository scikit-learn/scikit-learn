from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StructNamespace

if TYPE_CHECKING:
    from sqlframe.base.column import Column

    from narwhals._spark_like.expr import SparkLikeExpr


class SparkLikeExprStructNamespace(
    LazyExprNamespace["SparkLikeExpr"], StructNamespace["SparkLikeExpr"]
):
    def field(self, name: str) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            return expr.getField(name)

        return self.compliant._with_elementwise(func).alias(name)
