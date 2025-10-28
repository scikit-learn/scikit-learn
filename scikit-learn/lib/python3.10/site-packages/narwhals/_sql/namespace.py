from __future__ import annotations

import operator
from functools import reduce
from typing import TYPE_CHECKING, Any, Protocol

from narwhals._compliant import LazyNamespace
from narwhals._compliant.typing import NativeExprT, NativeFrameT
from narwhals._sql.typing import SQLExprT, SQLLazyFrameT

if TYPE_CHECKING:
    from collections.abc import Iterable

    from narwhals.typing import PythonLiteral


class SQLNamespace(
    LazyNamespace[SQLLazyFrameT, SQLExprT, NativeFrameT],
    Protocol[SQLLazyFrameT, SQLExprT, NativeFrameT, NativeExprT],
):
    def _function(self, name: str, *args: NativeExprT | PythonLiteral) -> NativeExprT: ...
    def _lit(self, value: Any) -> NativeExprT: ...
    def _when(
        self,
        condition: NativeExprT,
        value: NativeExprT,
        otherwise: NativeExprT | None = None,
    ) -> NativeExprT: ...
    def _coalesce(self, *exprs: NativeExprT) -> NativeExprT: ...

    # Horizontal functions
    def any_horizontal(self, *exprs: SQLExprT, ignore_nulls: bool) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            if ignore_nulls:
                cols = (self._coalesce(col, self._lit(False)) for col in cols)
            return reduce(operator.or_, cols)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def all_horizontal(self, *exprs: SQLExprT, ignore_nulls: bool) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            if ignore_nulls:
                cols = (self._coalesce(col, self._lit(True)) for col in cols)
            return reduce(operator.and_, cols)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def max_horizontal(self, *exprs: SQLExprT) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            return self._function("greatest", *cols)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def min_horizontal(self, *exprs: SQLExprT) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            return self._function("least", *cols)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def sum_horizontal(self, *exprs: SQLExprT) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            return reduce(
                operator.add, (self._coalesce(col, self._lit(0)) for col in cols)
            )

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    # Other
    def coalesce(self, *exprs: SQLExprT) -> SQLExprT:
        def func(cols: Iterable[NativeExprT]) -> NativeExprT:
            return self._coalesce(*cols)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)
