from __future__ import annotations

import contextlib
from decimal import Decimal as D
from functools import lru_cache
from typing import TYPE_CHECKING, Any, overload

from polars import functions as F
from polars._utils.parse import parse_into_expression
from polars._utils.wrap import wrap_expr
from polars.datatypes import (
    Array,
    Boolean,
    Decimal,
    Float64,
    List,
    Utf8,
)
from polars.datatypes.group import FLOAT_DTYPES, INTEGER_DTYPES

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr


if TYPE_CHECKING:
    from typing import Literal

    from polars import Expr, Series
    from polars._typing import IntoExpr, PolarsDataType


# create a lookup of dtypes that have a reasonable one/zero mapping; for
# anything more elaborate should use `repeat`
@lru_cache(16)
def _one_or_zero_by_dtype(value: int, dtype: PolarsDataType) -> Any:
    if dtype in INTEGER_DTYPES:
        return value
    elif dtype in FLOAT_DTYPES:
        return float(value)
    elif dtype == Boolean:
        return bool(value)
    elif dtype == Utf8:
        return str(value)
    elif isinstance(dtype, Decimal):
        return D(value)
    elif isinstance(dtype, (List, Array)):
        arr_width = getattr(dtype, "size", 1)
        return [_one_or_zero_by_dtype(value, dtype.inner)] * arr_width
    return None


@overload
def repeat(
    value: IntoExpr | None,
    n: int | Expr,
    *,
    dtype: PolarsDataType | None = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def repeat(
    value: IntoExpr | None,
    n: int | Expr,
    *,
    dtype: PolarsDataType | None = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def repeat(
    value: IntoExpr | None,
    n: int | Expr,
    *,
    dtype: PolarsDataType | None = ...,
    eager: bool,
) -> Expr | Series: ...


def repeat(
    value: IntoExpr | None,
    n: int | Expr,
    *,
    dtype: PolarsDataType | None = None,
    eager: bool = False,
) -> Expr | Series:
    """
    Construct a column of length `n` filled with the given value.

    Parameters
    ----------
    value
        Value to repeat.
    n
        Length of the resulting column.
    dtype
        Data type of the resulting column. If set to `None` (default), data type is
        inferred from the given value. Defaults to Int32 for integer values, unless
        Int64 is required to fit the given value. Defaults to Float64 for float values.
    eager
        Evaluate immediately and return a `Series`. If set to `False` (default),
        return an expression instead.

    Notes
    -----
    If you want to construct a column in lazy mode and do not need a pre-determined
    length, use :func:`lit` instead.

    See Also
    --------
    lit

    Examples
    --------
    Construct a column with a repeated value in a lazy context.

    >>> pl.select(pl.repeat("z", n=3)).to_series()
    shape: (3,)
    Series: 'repeat' [str]
    [
            "z"
            "z"
            "z"
    ]

    Generate a Series directly by setting `eager=True`.

    >>> pl.repeat(3, n=3, dtype=pl.Int8, eager=True)
    shape: (3,)
    Series: 'repeat' [i8]
    [
            3
            3
            3
    ]
    """
    if isinstance(n, int):
        n = F.lit(n)
    value = parse_into_expression(value, str_as_lit=True, dtype=dtype)
    expr = wrap_expr(plr.repeat(value, n._pyexpr, dtype))
    if eager:
        return F.select(expr).to_series()
    return expr


@overload
def ones(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def ones(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: Literal[True],
) -> Series: ...


@overload
def ones(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: bool,
) -> Expr | Series: ...


def ones(
    n: int | Expr,
    dtype: PolarsDataType = Float64,
    *,
    eager: bool = False,
) -> Expr | Series:
    """
    Construct a column of length `n` filled with ones.

    This is syntactic sugar for the `repeat` function.

    Parameters
    ----------
    n
        Length of the resulting column.
    dtype
        Data type of the resulting column. Defaults to Float64.
    eager
        Evaluate immediately and return a `Series`. If set to `False`,
        return an expression instead.

    Notes
    -----
    If you want to construct a column in lazy mode and do not need a pre-determined
    length, use :func:`lit` instead.

    See Also
    --------
    repeat
    lit

    Examples
    --------
    >>> pl.ones(3, pl.Int8, eager=True)
    shape: (3,)
    Series: 'ones' [i8]
    [
        1
        1
        1
    ]
    """
    if (one := _one_or_zero_by_dtype(1, dtype)) is None:
        msg = f"invalid dtype for `ones`; found {dtype}"
        raise TypeError(msg)

    return repeat(one, n=n, dtype=dtype, eager=eager).alias("ones")


@overload
def zeros(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def zeros(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: Literal[True],
) -> Series: ...


@overload
def zeros(
    n: int | Expr,
    dtype: PolarsDataType = ...,
    *,
    eager: bool,
) -> Expr | Series: ...


def zeros(
    n: int | Expr,
    dtype: PolarsDataType = Float64,
    *,
    eager: bool = False,
) -> Expr | Series:
    """
    Construct a column of length `n` filled with zeros.

    This is syntactic sugar for the `repeat` function.

    Parameters
    ----------
    n
        Length of the resulting column.
    dtype
        Data type of the resulting column. Defaults to Float64.
    eager
        Evaluate immediately and return a `Series`. If set to `False`,
        return an expression instead.

    Notes
    -----
    If you want to construct a column in lazy mode and do not need a pre-determined
    length, use :func:`lit` instead.

    See Also
    --------
    repeat
    lit

    Examples
    --------
    >>> pl.zeros(3, pl.Int8, eager=True)
    shape: (3,)
    Series: 'zeros' [i8]
    [
        0
        0
        0
    ]
    """
    if (zero := _one_or_zero_by_dtype(0, dtype)) is None:
        msg = f"invalid dtype for `zeros`; found {dtype}"
        raise TypeError(msg)

    return repeat(zero, n=n, dtype=dtype, eager=eager).alias("zeros")
