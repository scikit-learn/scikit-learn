from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace
from narwhals._duckdb.utils import F, lit, when

if TYPE_CHECKING:
    from duckdb import Expression

    from narwhals._duckdb.expr import DuckDBExpr
    from narwhals.typing import NonNestedLiteral


class DuckDBExprListNamespace(
    LazyExprNamespace["DuckDBExpr"], ListNamespace["DuckDBExpr"]
):
    def len(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("len", expr))

    def unique(self) -> DuckDBExpr:
        def func(expr: Expression) -> Expression:
            expr_distinct = F("list_distinct", expr)
            return when(
                F("array_position", expr, lit(None)).isnotnull(),
                F("list_append", expr_distinct, lit(None)),
            ).otherwise(expr_distinct)

        return self.compliant._with_callable(func)

    def contains(self, item: NonNestedLiteral) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("list_contains", expr, lit(item))
        )
