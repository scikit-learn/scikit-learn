from __future__ import annotations

import polars._reexport as pl


class DataTypeExprArrNameSpace:
    """Namespace for arr datatype expressions."""

    _accessor = "arr"

    def __init__(self, expr: pl.DataTypeExpr) -> None:
        self._pydatatype_expr = expr._pydatatype_expr

    def inner_dtype(self) -> pl.DataTypeExpr:
        """Get the inner DataType of array."""
        return pl.DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.arr_inner_dtype()
        )

    def width(self) -> pl.Expr:
        """
        Get the array width.

        Examples
        --------
        >>> pl.select(pl.Array(pl.Int8, (1, 2, 3)).to_dtype_expr().arr.width())
        shape: (1, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ u32     │
        ╞═════════╡
        │ 1       │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(self._pydatatype_expr.arr_width())

    def shape(self) -> pl.Expr:
        """
        Get the array shape.

        Examples
        --------
        >>> pl.select(pl.Array(pl.Int8, (1, 2, 3)).to_dtype_expr().arr.shape())
        shape: (3, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ u32     │
        ╞═════════╡
        │ 1       │
        │ 2       │
        │ 3       │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(self._pydatatype_expr.arr_shape())
