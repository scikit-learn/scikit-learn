from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from narwhals._compliant.group_by import CompliantGroupBy, ParseKeysGroupBy
from narwhals._compliant.typing import CompliantLazyFrameT, NativeExprT_co
from narwhals._sql.typing import SQLExprT_contra
from narwhals._utils import zip_strict

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator


class SQLGroupBy(
    ParseKeysGroupBy[CompliantLazyFrameT, SQLExprT_contra],
    CompliantGroupBy[CompliantLazyFrameT, SQLExprT_contra],
    Protocol[CompliantLazyFrameT, SQLExprT_contra, NativeExprT_co],
):
    _keys: list[str]
    _output_key_names: list[str]

    def _evaluate_expr(self, expr: SQLExprT_contra, /) -> Iterator[NativeExprT_co]:
        output_names = expr._evaluate_output_names(self.compliant)
        aliases = (
            expr._alias_output_names(output_names)
            if expr._alias_output_names
            else output_names
        )
        native_exprs = expr(self.compliant)
        if expr._is_multi_output_unnamed():
            exclude = {*self._keys, *self._output_key_names}
            for native_expr, name, alias in zip_strict(
                native_exprs, output_names, aliases
            ):
                if name not in exclude:
                    yield expr._alias_native(native_expr, alias)
        else:
            for native_expr, alias in zip_strict(native_exprs, aliases):
                yield expr._alias_native(native_expr, alias)

    def _evaluate_exprs(
        self, exprs: Iterable[SQLExprT_contra], /
    ) -> Iterator[NativeExprT_co]:
        for expr in exprs:
            yield from self._evaluate_expr(expr)
