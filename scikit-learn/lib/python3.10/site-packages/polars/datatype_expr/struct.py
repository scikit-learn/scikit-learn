from __future__ import annotations

import polars._reexport as pl
from polars._utils.various import qualified_type_name


class DataTypeExprStructNameSpace:
    """Namespace for struct datatype expressions."""

    _accessor = "struct"

    def __init__(self, expr: pl.DataTypeExpr) -> None:
        self._pydatatype_expr = expr._pydatatype_expr

    def __getitem__(self, item: str | int) -> pl.DataTypeExpr:
        if isinstance(item, str):
            return self.field_dtype(item)
        elif isinstance(item, int):
            return pl.DataTypeExpr._from_pydatatype_expr(
                self._pydatatype_expr.struct_field_dtype_by_index(item)
            )
        else:
            msg = f"expected type 'int | str', got {qualified_type_name(item)!r} ({item!r})"
            raise TypeError(msg)

    def field_dtype(self, field_name: str) -> pl.DataTypeExpr:
        """

        Get the DataType of field with a specific field name.

        Notes
        -----
        The `struct` namespace has implemented `__getitem__` so you can also access
        fields by index:

        >>> (
        ...     pl.Struct({"x": pl.Int64, "y": pl.String})
        ...     .to_dtype_expr()
        ...     .struct[1]
        ...     .collect_dtype({})
        ... )
        String
        """
        return pl.DataTypeExpr._from_pydatatype_expr(
            self._pydatatype_expr.struct_field_dtype_by_name(field_name)
        )

    def field_names(self) -> pl.Expr:
        """
        Get the field names in a struct as a list.

        Examples
        --------
        >>> pl.select(
        ...     pl.Struct({"x": pl.Int64, "y": pl.String})
        ...     .to_dtype_expr()
        ...     .struct.field_names()
        ... )
        shape: (2, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ str     │
        ╞═════════╡
        │ x       │
        │ y       │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(self._pydatatype_expr.struct_field_names())
