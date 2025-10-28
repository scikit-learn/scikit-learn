from __future__ import annotations

from typing import TYPE_CHECKING

import polars._reexport as pl
from polars import functions as F
from polars._utils.unstable import unstable
from polars._utils.various import qualified_type_name

if TYPE_CHECKING:
    from collections.abc import Mapping

    from polars import Expr
    from polars._typing import PolarsDataType


@unstable()
def dtype_of(col_or_expr: str | Expr) -> pl.DataTypeExpr:
    """
    Get a lazily evaluated :class:`DataType` of a column or expression.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    Examples
    --------
    >>> def inspect(expr: pl.Expr) -> pl.Expr:
    ...     def print_and_return(s: pl.Series) -> pl.Series:
    ...         print(s)
    ...         return s
    ...
    ...     return expr.map_batches(
    ...         print_and_return,
    ...         # Clarify that the expression returns the same datatype as the input
    ...         # datatype.
    ...         return_dtype=pl.dtype_of(expr),
    ...     )
    >>> df = pl.DataFrame(
    ...     {
    ...         "UserID": [1, 2, 3, 4, 5],
    ...         "Name": ["Alice", "Bob", "Charlie", "Diana", "Ethan"],
    ...     }
    ... )
    >>> df.select(inspect(pl.col("Name")))
    shape: (5,)
    Series: 'Name' [str]
    [
        "Alice"
        "Bob"
        "Charlie"
        "Diana"
        "Ethan"
    ]
    shape: (5, 1)
    ┌─────────┐
    │ Name    │
    │ ---     │
    │ str     │
    ╞═════════╡
    │ Alice   │
    │ Bob     │
    │ Charlie │
    │ Diana   │
    │ Ethan   │
    └─────────┘
    """
    from polars._plr import PyDataTypeExpr

    e: Expr
    if isinstance(col_or_expr, str):
        e = F.col(col_or_expr)
    else:
        e = col_or_expr

    return pl.DataTypeExpr._from_pydatatype_expr(PyDataTypeExpr.of_expr(e._pyexpr))


@unstable()
def self_dtype() -> pl.DataTypeExpr:
    """
    Get the dtype of `self` in `map_elements` and `map_batches`.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """
    from polars._plr import PyDataTypeExpr

    return pl.DataTypeExpr._from_pydatatype_expr(PyDataTypeExpr.self_dtype())


@unstable()
def struct_with_fields(
    mapping: Mapping[str, PolarsDataType | pl.DataTypeExpr],
) -> pl.DataTypeExpr:
    """
    Create a new datatype expression that represents a Struct datatype.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """
    from polars._plr import PyDataTypeExpr

    def preprocess(dtype_expr: PolarsDataType | pl.DataTypeExpr) -> PyDataTypeExpr:
        if isinstance(dtype_expr, pl.DataType):
            return dtype_expr.to_dtype_expr()._pydatatype_expr
        if isinstance(dtype_expr, pl.DataTypeClass):
            return dtype_expr.to_dtype_expr()._pydatatype_expr
        elif isinstance(dtype_expr, pl.DataTypeExpr):
            return dtype_expr._pydatatype_expr
        else:
            msg = f"mapping item must be a datatype or datatype expression; found {qualified_type_name(dtype_expr)!r}"
            raise TypeError(msg)

    fields = [(name, preprocess(dtype_expr)) for (name, dtype_expr) in mapping.items()]

    return pl.DataTypeExpr._from_pydatatype_expr(
        PyDataTypeExpr.struct_with_fields(fields)
    )
