from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace
from narwhals._duckdb.utils import F, lit, when
from narwhals._utils import requires

if TYPE_CHECKING:
    from duckdb import Expression

    from narwhals._duckdb.expr import DuckDBExpr
    from narwhals.typing import NonNestedLiteral


class DuckDBExprListNamespace(
    LazyExprNamespace["DuckDBExpr"], ListNamespace["DuckDBExpr"]
):
    def len(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("len", expr))

    @requires.backend_version((1, 3))  # bugged before 1.3
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

    def get(self, index: int) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("list_extract", expr, lit(index + 1))
        )
