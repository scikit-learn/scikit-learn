from __future__ import annotations

import polars._reexport as pl


class DataTypeExprListNameSpace:
    """Namespace for list datatype expressions."""

    _accessor = "list"

    def __init__(self, expr: pl.DataTypeExpr) -> None:
        self._pydatatype_expr = expr._pydatatype_expr

    def inner_dtype(self) -> pl.DataTypeExpr:
        """Get the inner DataType of list."""
        return pl.DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.list_inner_dtype()
        )
