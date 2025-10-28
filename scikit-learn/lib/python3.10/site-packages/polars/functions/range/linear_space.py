from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, overload

from polars import functions as F
from polars._utils.parse import parse_into_expression
from polars._utils.unstable import unstable
from polars._utils.wrap import wrap_expr

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

from typing import Literal

if TYPE_CHECKING:
    from polars import Expr, Series
    from polars._typing import (
        ClosedInterval,
        IntoExpr,
        IntoExprColumn,
        NumericLiteral,
        TemporalLiteral,
    )


@overload
def linear_space(
    start: NumericLiteral | TemporalLiteral | IntoExpr,
    end: NumericLiteral | TemporalLiteral | IntoExpr,
    num_samples: int | IntoExpr,
    *,
    closed: ClosedInterval = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def linear_space(
    start: NumericLiteral | TemporalLiteral | IntoExpr,
    end: NumericLiteral | TemporalLiteral | IntoExpr,
    num_samples: int | IntoExpr,
    *,
    closed: ClosedInterval = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def linear_space(
    start: NumericLiteral | TemporalLiteral | IntoExpr,
    end: NumericLiteral | TemporalLiteral | IntoExpr,
    num_samples: int | IntoExpr,
    *,
    closed: ClosedInterval = ...,
    eager: bool,
) -> Expr | Series: ...


@unstable()
def linear_space(
    start: NumericLiteral | TemporalLiteral | IntoExpr,
    end: NumericLiteral | TemporalLiteral | IntoExpr,
    num_samples: int | IntoExpr,
    *,
    closed: ClosedInterval = "both",
    eager: bool = False,
) -> Expr | Series:
    """
    Create sequence of evenly-spaced points.

    Parameters
    ----------
    start
        Lower bound of the range.
    end
        Upper bound of the range.
    num_samples
        Number of samples in the output sequence.
    closed : {'both', 'left', 'right', 'none'}
        Define which sides of the interval are closed (inclusive).
    eager
        Evaluate immediately and return a `Series`.
        If set to `False` (default), return an expression instead.

    .. warning::
        This functionality is experimental. It may be changed at any point without it
        being considered a breaking change.

    Notes
    -----
    `linear_space` works with numeric and temporal dtypes. When the `start` and `end`
    parameters are `Date` dtypes, the output sequence consists of equally-spaced
    `Datetime` elements with millisecond precision.

    Returns
    -------
    Expr or Series
        Column of data type `:class:Time`.

    Examples
    --------
    >>> pl.linear_space(start=0, end=1, num_samples=3, eager=True)
    shape: (3,)
    Series: 'literal' [f64]
    [
            0.0
            0.5
            1.0
    ]
    >>> pl.linear_space(start=0, end=1, num_samples=3, closed="left", eager=True)
    shape: (3,)
    Series: 'literal' [f64]
    [
            0.0
            0.333333
            0.666667
    ]
    >>> pl.linear_space(start=0, end=1, num_samples=3, closed="right", eager=True)
    shape: (3,)
    Series: 'literal' [f64]
    [
            0.333333
            0.666667
            1.0
    ]
    >>> pl.linear_space(start=0, end=1, num_samples=3, closed="none", eager=True)
    shape: (3,)
    Series: 'literal' [f64]
    [
            0.25
            0.5
            0.75
    ]
    >>> from datetime import time
    >>> pl.linear_space(
    ...     start=time(hour=1), end=time(hour=12), num_samples=3, eager=True
    ... )
    shape: (3,)
    Series: 'literal' [time]
    [
            01:00:00
            06:30:00
            12:00:00
    ]

    `Date` endpoints generate a sequence of `Datetime` values:

    >>> from datetime import date
    >>> pl.linear_space(
    ...     start=date(2025, 1, 1),
    ...     end=date(2025, 2, 1),
    ...     num_samples=3,
    ...     closed="right",
    ...     eager=True,
    ... )
    shape: (3,)
    Series: 'literal' [datetime[μs]]
    [
            2025-01-11 08:00:00
            2025-01-21 16:00:00
            2025-02-01 00:00:00
    ]

    When `eager=False` (default), an expression is produced. You can generate a sequence
    using the length of the dataframe:

    >>> df = pl.DataFrame({"a": [1, 2, 3, 4, 5]})
    >>> df.with_columns(pl.linear_space(0, 1, pl.len()).alias("ls"))
    shape: (5, 2)
    ┌─────┬──────┐
    │ a   ┆ ls   │
    │ --- ┆ ---  │
    │ i64 ┆ f64  │
    ╞═════╪══════╡
    │ 1   ┆ 0.0  │
    │ 2   ┆ 0.25 │
    │ 3   ┆ 0.5  │
    │ 4   ┆ 0.75 │
    │ 5   ┆ 1.0  │
    └─────┴──────┘
    """
    start_pyexpr = parse_into_expression(start)
    end_pyexpr = parse_into_expression(end)
    num_samples_pyexpr = parse_into_expression(num_samples)
    result = wrap_expr(
        plr.linear_space(start_pyexpr, end_pyexpr, num_samples_pyexpr, closed)
    )

    if eager:
        return F.select(result).to_series()

    return result


@overload
def linear_spaces(
    start: NumericLiteral | TemporalLiteral | IntoExprColumn,
    end: NumericLiteral | TemporalLiteral | IntoExprColumn,
    num_samples: int | IntoExprColumn,
    *,
    closed: ClosedInterval = ...,
    as_array: bool = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def linear_spaces(
    start: NumericLiteral | TemporalLiteral | IntoExprColumn,
    end: NumericLiteral | TemporalLiteral | IntoExprColumn,
    num_samples: int | IntoExprColumn,
    *,
    closed: ClosedInterval = ...,
    as_array: bool = ...,
    eager: Literal[True],
) -> Series: ...


@overload
def linear_spaces(
    start: NumericLiteral | TemporalLiteral | IntoExprColumn,
    end: NumericLiteral | TemporalLiteral | IntoExprColumn,
    num_samples: int | IntoExprColumn,
    *,
    closed: ClosedInterval = ...,
    as_array: bool = ...,
    eager: bool,
) -> Expr | Series: ...


def linear_spaces(
    start: NumericLiteral | TemporalLiteral | IntoExprColumn,
    end: NumericLiteral | TemporalLiteral | IntoExprColumn,
    num_samples: int | IntoExprColumn,
    *,
    closed: ClosedInterval = "both",
    as_array: bool = False,
    eager: bool = False,
) -> Expr | Series:
    """
    Generate a sequence of evenly-spaced values for each row between `start` and `end`.

    The number of values in each sequence is determined by `num_samples`.

    Parameters
    ----------
    start
        Lower bound of the range.
    end
        Upper bound of the range.
    num_samples
        Number of samples in the output sequence.
    closed : {'both', 'left', 'right', 'none'}
        Define which sides of the interval are closed (inclusive).
    as_array
        Return result as a fixed-length `Array`. `num_samples` must be a constant.
    eager
        Evaluate immediately and return a `Series`.
        If set to `False` (default), return an expression instead.

    .. warning::
        This functionality is experimental. It may be changed at any point without it
        being considered a breaking change.

    Returns
    -------
    Expr or Series
        Column of data type `List(dtype)`.

    See Also
    --------
    linear_space : Generate a single sequence of linearly-spaced values.

    Examples
    --------
    >>> df = pl.DataFrame({"start": [1, -1], "end": [3, 2], "num_samples": [4, 5]})
    >>> df.with_columns(ls=pl.linear_spaces("start", "end", "num_samples"))
    shape: (2, 4)
    ┌───────┬─────┬─────────────┬────────────────────────┐
    │ start ┆ end ┆ num_samples ┆ ls                     │
    │ ---   ┆ --- ┆ ---         ┆ ---                    │
    │ i64   ┆ i64 ┆ i64         ┆ list[f64]              │
    ╞═══════╪═════╪═════════════╪════════════════════════╡
    │ 1     ┆ 3   ┆ 4           ┆ [1.0, 1.666667, … 3.0] │
    │ -1    ┆ 2   ┆ 5           ┆ [-1.0, -0.25, … 2.0]   │
    └───────┴─────┴─────────────┴────────────────────────┘
    >>> df.with_columns(ls=pl.linear_spaces("start", "end", 3, as_array=True))
    shape: (2, 4)
    ┌───────┬─────┬─────────────┬──────────────────┐
    │ start ┆ end ┆ num_samples ┆ ls               │
    │ ---   ┆ --- ┆ ---         ┆ ---              │
    │ i64   ┆ i64 ┆ i64         ┆ array[f64, 3]    │
    ╞═══════╪═════╪═════════════╪══════════════════╡
    │ 1     ┆ 3   ┆ 4           ┆ [1.0, 2.0, 3.0]  │
    │ -1    ┆ 2   ┆ 5           ┆ [-1.0, 0.5, 2.0] │
    └───────┴─────┴─────────────┴──────────────────┘
    """
    start_pyexpr = parse_into_expression(start)
    end_pyexpr = parse_into_expression(end)
    num_samples_pyexpr = parse_into_expression(num_samples)
    result = wrap_expr(
        plr.linear_spaces(
            start_pyexpr, end_pyexpr, num_samples_pyexpr, closed, as_array
        )
    )

    if eager:
        return F.select(result).to_series()

    return result
