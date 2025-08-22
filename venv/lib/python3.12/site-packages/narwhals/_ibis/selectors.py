from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import CompliantSelector, LazySelectorNamespace
from narwhals._ibis.expr import IbisExpr

if TYPE_CHECKING:
    import ibis.expr.types as ir  # noqa: F401

    from narwhals._ibis.dataframe import IbisLazyFrame  # noqa: F401
    from narwhals._ibis.expr import IbisWindowFunction


class IbisSelectorNamespace(LazySelectorNamespace["IbisLazyFrame", "ir.Value"]):
    @property
    def _selector(self) -> type[IbisSelector]:
        return IbisSelector


class IbisSelector(  # type: ignore[misc]
    CompliantSelector["IbisLazyFrame", "ir.Value"], IbisExpr
):
    _window_function: IbisWindowFunction | None = None

    def _to_expr(self) -> IbisExpr:
        return IbisExpr(
            self._call,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )
