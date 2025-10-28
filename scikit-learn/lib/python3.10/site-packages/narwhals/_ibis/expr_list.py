from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace

if TYPE_CHECKING:
    import ibis.expr.types as ir

    from narwhals._ibis.expr import IbisExpr
    from narwhals.typing import NonNestedLiteral


class IbisExprListNamespace(LazyExprNamespace["IbisExpr"], ListNamespace["IbisExpr"]):
    def len(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.length())

    def unique(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.unique())

    def contains(self, item: NonNestedLiteral) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.contains(item))

    def get(self, index: int) -> IbisExpr:
        def _get(expr: ir.ArrayColumn) -> ir.Column:
            return expr[index]

        return self.compliant._with_callable(_get)
