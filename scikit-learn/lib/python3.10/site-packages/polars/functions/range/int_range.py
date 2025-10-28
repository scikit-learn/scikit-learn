from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, overload

from polars import functions as F
from polars._utils.parse import parse_into_expression
from polars._utils.wrap import wrap_expr, wrap_s
from polars.datatypes import Int64
from polars.datatypes._parse import parse_into_datatype_expr

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

if TYPE_CHECKING:
    from typing import Literal

    from polars import DataTypeExpr, Expr, Series
    from polars._typing import IntoExprColumn, PolarsIntegerType


@overload
def arange(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def arange(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def arange(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: bool,
) -> Expr | Series: ...


def arange(
    start: int | IntoExprColumn = 0,
    end: int | IntoExprColumn | None = None,
    step: int = 1,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = Int64,
    eager: bool = False,
) -> Expr | Series:
    """
    Generate a range of integers.

    Alias for :func:`int_range`.

    Parameters
    ----------
    start
        Lower bound of the range (inclusive).
    end
        Upper bound of the range (exclusive).
    step
        Step size of the range.
    dtype
        Data type of the range. Defaults to `Int64`.
    eager
        Evaluate immediately and return a `Series`.
        If set to `False` (default), return an expression instead.

    Returns
    -------
    Expr or Series
        Column of integer data type `dtype`.

    See Also
    --------
    int_range : Generate a range of integers.
    int_ranges : Generate a range of integers for each row of the input columns.

    Examples
    --------
    >>> pl.arange(0, 3, eager=True)
    shape: (3,)
    Series: 'literal' [i64]
    [
            0
            1
            2
    ]
    """
    return int_range(start, end, step, dtype=dtype, eager=eager)


@overload
def int_range(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def int_range(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def int_range(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: bool,
) -> Expr | Series: ...


def int_range(
    start: int | IntoExprColumn = 0,
    end: int | IntoExprColumn | None = None,
    step: int = 1,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = Int64,
    eager: bool = False,
) -> Expr | Series:
    """
    Generate a range of integers.

    Parameters
    ----------
    start
        Start of the range (inclusive). Defaults to 0.
    end
        End of the range (exclusive). If set to `None` (default),
        the value of `start` is used and `start` is set to `0`.
    step
        Step size of the range.
    dtype
        Data type of the range.
    eager
        Evaluate immediately and return a `Series`.
        If set to `False` (default), return an expression instead.

    Returns
    -------
    Expr or Series
        Column of integer data type `dtype`.

    See Also
    --------
    int_ranges : Generate a range of integers for each row of the input columns.

    Examples
    --------
    >>> pl.int_range(0, 3, eager=True)
    shape: (3,)
    Series: 'literal' [i64]
    [
            0
            1
            2
    ]

    `end` can be omitted for a shorter syntax.

    >>> pl.int_range(3, eager=True)
    shape: (3,)
    Series: 'literal' [i64]
    [
            0
            1
            2
    ]

    Generate an index column by using `int_range` in conjunction with :func:`len`.

    >>> df = pl.DataFrame({"a": [1, 3, 5], "b": [2, 4, 6]})
    >>> df.select(
    ...     pl.int_range(pl.len(), dtype=pl.UInt32).alias("index"),
    ...     pl.all(),
    ... )
    shape: (3, 3)
    ┌───────┬─────┬─────┐
    │ index ┆ a   ┆ b   │
    │ ---   ┆ --- ┆ --- │
    │ u32   ┆ i64 ┆ i64 │
    ╞═══════╪═════╪═════╡
    │ 0     ┆ 1   ┆ 2   │
    │ 1     ┆ 3   ┆ 4   │
    │ 2     ┆ 5   ┆ 6   │
    └───────┴─────┴─────┘
    """
    if end is None:
        end = start
        start = 0

    dtype_expr = parse_into_datatype_expr(dtype)
    if isinstance(start, int) and isinstance(end, int) and eager:
        return wrap_s(
            plr.eager_int_range(start, end, step, dtype_expr._pydatatype_expr)
        )

    start_pyexpr = parse_into_expression(start)
    end_pyexpr = parse_into_expression(end)
    result = wrap_expr(
        plr.int_range(start_pyexpr, end_pyexpr, step, dtype_expr._pydatatype_expr)
    )

    if eager:
        return F.select(result).to_series()

    return result


@overload
def int_ranges(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int | IntoExprColumn = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def int_ranges(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int | IntoExprColumn = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def int_ranges(
    start: int | IntoExprColumn = ...,
    end: int | IntoExprColumn | None = ...,
    step: int | IntoExprColumn = ...,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = ...,
    eager: bool,
) -> Expr | Series: ...


def int_ranges(
    start: int | IntoExprColumn = 0,
    end: int | IntoExprColumn | None = None,
    step: int | IntoExprColumn = 1,
    *,
    dtype: PolarsIntegerType | DataTypeExpr = Int64,
    eager: bool = False,
) -> Expr | Series:
    """
    Generate a range of integers for each row of the input columns.

    Parameters
    ----------
    start
        Start of the range (inclusive). Defaults to 0.
    end
        End of the range (exclusive). If set to `None` (default),
        the value of `start` is used and `start` is set to `0`.
    step
        Step size of the range.
    dtype
        Integer data type of the ranges. Defaults to `Int64`.
    eager
        Evaluate immediately and return a `Series`.
        If set to `False` (default), return an expression instead.

    Returns
    -------
    Expr or Series
        Column of data type `List(dtype)`.

    See Also
    --------
    int_range : Generate a single range of integers.

    Examples
    --------
    >>> df = pl.DataFrame({"start": [1, -1], "end": [3, 2]})
    >>> df.with_columns(int_range=pl.int_ranges("start", "end"))
    shape: (2, 3)
    ┌───────┬─────┬────────────┐
    │ start ┆ end ┆ int_range  │
    │ ---   ┆ --- ┆ ---        │
    │ i64   ┆ i64 ┆ list[i64]  │
    ╞═══════╪═════╪════════════╡
    │ 1     ┆ 3   ┆ [1, 2]     │
    │ -1    ┆ 2   ┆ [-1, 0, 1] │
    └───────┴─────┴────────────┘

    `end` can be omitted for a shorter syntax.

    >>> df.select("end", int_range=pl.int_ranges("end"))
    shape: (2, 2)
    ┌─────┬───────────┐
    │ end ┆ int_range │
    │ --- ┆ ---       │
    │ i64 ┆ list[i64] │
    ╞═════╪═══════════╡
    │ 3   ┆ [0, 1, 2] │
    │ 2   ┆ [0, 1]    │
    └─────┴───────────┘
    """
    if end is None:
        end = start
        start = 0

    dtype_expr = parse_into_datatype_expr(dtype)
    start_pyexpr = parse_into_expression(start)
    end_pyexpr = parse_into_expression(end)
    step_pyexpr = parse_into_expression(step)
    result = wrap_expr(
        plr.int_ranges(
            start_pyexpr, end_pyexpr, step_pyexpr, dtype_expr._pydatatype_expr
        )
    )

    if eager:
        return F.select(result).to_series()

    return result
