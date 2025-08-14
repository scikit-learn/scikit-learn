from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StructNamespace
from narwhals._duckdb.utils import F, lit

if TYPE_CHECKING:
    from narwhals._duckdb.expr import DuckDBExpr


class DuckDBExprStructNamespace(
    LazyExprNamespace["DuckDBExpr"], StructNamespace["DuckDBExpr"]
):
    def field(self, name: str) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("struct_extract", expr, lit(name))
        ).alias(name)
