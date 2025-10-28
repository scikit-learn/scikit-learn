from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._arrow.expr import ArrowExpr
from narwhals._compliant import CompliantSelector, EagerSelectorNamespace

if TYPE_CHECKING:
    from narwhals._arrow.dataframe import ArrowDataFrame  # noqa: F401
    from narwhals._arrow.series import ArrowSeries  # noqa: F401
    from narwhals._compliant.typing import ScalarKwargs


class ArrowSelectorNamespace(EagerSelectorNamespace["ArrowDataFrame", "ArrowSeries"]):
    @property
    def _selector(self) -> type[ArrowSelector]:
        return ArrowSelector


class ArrowSelector(CompliantSelector["ArrowDataFrame", "ArrowSeries"], ArrowExpr):  # type: ignore[misc]
    _depth: int = 0
    _scalar_kwargs: ScalarKwargs = {}  # noqa: RUF012
    _function_name: str = "selector"

    def _to_expr(self) -> ArrowExpr:
        return ArrowExpr(
            self._call,
            depth=self._depth,
            function_name=self._function_name,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )
