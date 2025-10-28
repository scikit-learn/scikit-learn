from __future__ import annotations

import contextlib
import warnings
from typing import TYPE_CHECKING, Any, Callable, overload

import polars._reexport as pl
import polars.functions as F
import polars.selectors as cs
from polars._dependencies import _check_for_numpy
from polars._dependencies import numpy as np
from polars._utils.async_ import _AioDataFrameResult, _GeventDataFrameResult
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    deprecate_streaming_parameter,
    deprecated,
    issue_deprecation_warning,
)
from polars._utils.parse import (
    parse_into_expression,
    parse_into_list_of_expressions,
)
from polars._utils.unstable import issue_unstable_warning, unstable
from polars._utils.various import extend_bool, qualified_type_name
from polars._utils.wrap import wrap_df, wrap_expr, wrap_s
from polars.datatypes import DTYPE_TEMPORAL_UNITS, Date, Datetime, Int64
from polars.datatypes._parse import parse_into_datatype_expr
from polars.lazyframe.opt_flags import (
    DEFAULT_QUERY_OPT_FLAGS,
    forward_old_opt_flags,
)
from polars.meta.index_type import get_index_type

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

if TYPE_CHECKING:
    import sys
    from collections.abc import Awaitable, Collection, Iterable, Sequence
    from typing import Literal

    from polars import DataFrame, Expr, LazyFrame, Series
    from polars._typing import (
        CorrelationMethod,
        EngineType,
        EpochTimeUnit,
        IntoExpr,
        PolarsDataType,
        QuantileMethod,
    )
    from polars.lazyframe.opt_flags import (
        QueryOptFlags,
    )

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


def field(name: str | list[str]) -> Expr:
    """
    Select a field in the current `struct.with_fields` scope.

    Parameters
    ----------
    name
        Name of the field(s) to select.

    Examples
    --------
    >>> df = pl.DataFrame({"a": [{"x": 5, "y": 2}, {"x": 3, "y": 4}]})
    >>> df.select(pl.col("a").struct.with_fields(pl.field("x") ** 2))
    shape: (2, 1)
    ┌───────────┐
    │ a         │
    │ ---       │
    │ struct[2] │
    ╞═══════════╡
    │ {25,2}    │
    │ {9,4}     │
    └───────────┘
    """
    if isinstance(name, str):
        name = [name]
    return wrap_expr(plr.field(name))


def element() -> Expr:
    """
    Alias for an element being evaluated in an `eval` or `filter` expression.

    Examples
    --------
    A horizontal rank computation by taking the elements of a list

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...     }
    ... )
    >>> df.with_columns(
    ...     pl.concat_list(["a", "b"]).list.eval(pl.element().rank()).alias("rank")
    ... )
    shape: (3, 3)
    ┌─────┬─────┬────────────┐
    │ a   ┆ b   ┆ rank       │
    │ --- ┆ --- ┆ ---        │
    │ i64 ┆ i64 ┆ list[f64]  │
    ╞═════╪═════╪════════════╡
    │ 1   ┆ 4   ┆ [1.0, 2.0] │
    │ 8   ┆ 5   ┆ [2.0, 1.0] │
    │ 3   ┆ 2   ┆ [2.0, 1.0] │
    └─────┴─────┴────────────┘

    A mathematical operation on array elements

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...     }
    ... )
    >>> df.with_columns(
    ...     pl.concat_list(["a", "b"]).list.eval(pl.element() * 2).alias("a_b_doubled")
    ... )
    shape: (3, 3)
    ┌─────┬─────┬─────────────┐
    │ a   ┆ b   ┆ a_b_doubled │
    │ --- ┆ --- ┆ ---         │
    │ i64 ┆ i64 ┆ list[i64]   │
    ╞═════╪═════╪═════════════╡
    │ 1   ┆ 4   ┆ [2, 8]      │
    │ 8   ┆ 5   ┆ [16, 10]    │
    │ 3   ┆ 2   ┆ [6, 4]      │
    └─────┴─────┴─────────────┘

    A filter operation on list elements

    >>> import polars as pl
    >>> df = pl.DataFrame({"a": [1, 8, 3], "b": [4, 5, 2]})
    >>> df.with_columns(
    ...     evens=pl.concat_list("a", "b").list.filter(pl.element() % 2 == 0)
    ... )
    shape: (3, 3)
    ┌─────┬─────┬───────────┐
    │ a   ┆ b   ┆ evens     │
    │ --- ┆ --- ┆ ---       │
    │ i64 ┆ i64 ┆ list[i64] │
    ╞═════╪═════╪═══════════╡
    │ 1   ┆ 4   ┆ [4]       │
    │ 8   ┆ 5   ┆ [8]       │
    │ 3   ┆ 2   ┆ [2]       │
    └─────┴─────┴───────────┘
    """
    return F.col("")


def count(*columns: str) -> Expr:
    """
    Return the number of non-null values in the column.

    This function is syntactic sugar for `col(columns).count()`.

    Calling this function without any arguments returns the number of rows in the
    context. **This way of using the function is deprecated.** Please use :func:`len`
    instead.

    Parameters
    ----------
    *columns
        One or more column names.

    Returns
    -------
    Expr
        Expression of data type :class:`UInt32`.

    See Also
    --------
    Expr.count

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, None],
    ...         "b": [3, None, None],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.count("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 2   │
    └─────┘

    Return the number of non-null values in multiple columns.

    >>> df.select(pl.count("b", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ b   ┆ c   │
    │ --- ┆ --- │
    │ u32 ┆ u32 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    └─────┴─────┘

    Return the number of rows in a context. **This way of using the function is
    deprecated.** Please use :func:`len` instead.

    >>> df.select(pl.count())  # doctest: +SKIP
    shape: (1, 1)
    ┌───────┐
    │ count │
    │ ---   │
    │ u32   │
    ╞═══════╡
    │ 3     │
    └───────┘
    """
    if not columns:
        issue_deprecation_warning(
            "`pl.count()` is deprecated. Please use `pl.len()` instead.",
            version="0.20.5",
        )
        return F.len().alias("count")
    return F.col(*columns).count()


def cum_count(*columns: str, reverse: bool = False) -> Expr:
    """
    Return the cumulative count of the non-null values in the column.

    This function is syntactic sugar for `col(columns).cum_count()`.

    Parameters
    ----------
    *columns
        Name(s) of the columns to use.
    reverse
        Reverse the operation.

    Examples
    --------
    >>> df = pl.DataFrame({"a": [1, 2, None], "b": [3, None, None]})
    >>> df.with_columns(
    ...     ca=pl.cum_count("a"),
    ...     cb=pl.cum_count("b"),
    ... )
    shape: (3, 4)
    ┌──────┬──────┬─────┬─────┐
    │ a    ┆ b    ┆ ca  ┆ cb  │
    │ ---  ┆ ---  ┆ --- ┆ --- │
    │ i64  ┆ i64  ┆ u32 ┆ u32 │
    ╞══════╪══════╪═════╪═════╡
    │ 1    ┆ 3    ┆ 1   ┆ 1   │
    │ 2    ┆ null ┆ 2   ┆ 1   │
    │ null ┆ null ┆ 2   ┆ 1   │
    └──────┴──────┴─────┴─────┘
    """
    return F.col(*columns).cum_count(reverse=reverse)


def implode(*columns: str) -> Expr:
    """
    Aggregate all column values into a list.

    This function is syntactic sugar for `pl.col(name).implode()`.

    Parameters
    ----------
    *columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [9, 8, 7],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.implode("a"))
    shape: (1, 1)
    ┌───────────┐
    │ a         │
    │ ---       │
    │ list[i64] │
    ╞═══════════╡
    │ [1, 2, 3] │
    └───────────┘
    >>> df.select(pl.implode("b", "c"))
    shape: (1, 2)
    ┌───────────┬───────────────────────┐
    │ b         ┆ c                     │
    │ ---       ┆ ---                   │
    │ list[i64] ┆ list[str]             │
    ╞═══════════╪═══════════════════════╡
    │ [9, 8, 7] ┆ ["foo", "bar", "foo"] │
    └───────────┴───────────────────────┘

    """
    return F.col(*columns).implode()


def std(column: str, ddof: int = 1) -> Expr:
    """
    Get the standard deviation.

    This function is syntactic sugar for `pl.col(column).std(ddof)`.

    Parameters
    ----------
    column
        Column name.
    ddof
        “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
        where N represents the number of elements.
        By default ddof is 1.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.std("a"))
    shape: (1, 1)
    ┌──────────┐
    │ a        │
    │ ---      │
    │ f64      │
    ╞══════════╡
    │ 3.605551 │
    └──────────┘
    >>> df["a"].std()
    3.605551275463989
    """
    return F.col(column).std(ddof)


def var(column: str, ddof: int = 1) -> Expr:
    """
    Get the variance.

    This function is syntactic sugar for `pl.col(column).var(ddof)`.

    Parameters
    ----------
    column
        Column name.
    ddof
        “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
        where N represents the number of elements.
        By default ddof is 1.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     },
    ... )
    >>> df.select(pl.var("a"))
    shape: (1, 1)
    ┌──────┐
    │ a    │
    │ ---  │
    │ f64  │
    ╞══════╡
    │ 13.0 │
    └──────┘
    >>> df["a"].var()
    13.0
    """
    return F.col(column).var(ddof)


def mean(*columns: str) -> Expr:
    """
    Get the mean value.

    This function is syntactic sugar for `pl.col(columns).mean()`.

    Parameters
    ----------
    *columns
        One or more column names.

    See Also
    --------
    mean_horizontal

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.mean("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 4.0 │
    └─────┘
    >>> df.select(pl.mean("a", "b"))
    shape: (1, 2)
    ┌─────┬──────────┐
    │ a   ┆ b        │
    │ --- ┆ ---      │
    │ f64 ┆ f64      │
    ╞═════╪══════════╡
    │ 4.0 ┆ 3.666667 │
    └─────┴──────────┘

    """
    return F.col(*columns).mean()


def median(*columns: str) -> Expr:
    """
    Get the median value.

    This function is syntactic sugar for `pl.col(columns).median()`.

    Parameters
    ----------
    columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.median("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 3.0 │
    └─────┘
    >>> df.select(pl.median("a", "b"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ f64 ┆ f64 │
    ╞═════╪═════╡
    │ 3.0 ┆ 4.0 │
    └─────┴─────┘

    """
    return F.col(*columns).median()


def n_unique(*columns: str) -> Expr:
    """
    Count unique values.

    This function is syntactic sugar for `pl.col(columns).n_unique()`.

    Parameters
    ----------
    columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 1],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.n_unique("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 2   │
    └─────┘
    >>> df.select(pl.n_unique("b", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ b   ┆ c   │
    │ --- ┆ --- │
    │ u32 ┆ u32 │
    ╞═════╪═════╡
    │ 3   ┆ 2   │
    └─────┴─────┘

    """
    return F.col(*columns).n_unique()


def approx_n_unique(*columns: str) -> Expr:
    """
    Approximate count of unique values.

    This function is syntactic sugar for `pl.col(columns).approx_n_unique()`, and
    uses the HyperLogLog++ algorithm for cardinality estimation.

    Parameters
    ----------
    columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 1],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.approx_n_unique("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 2   │
    └─────┘
    >>> df.select(pl.approx_n_unique("b", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ b   ┆ c   │
    │ --- ┆ --- │
    │ u32 ┆ u32 │
    ╞═════╪═════╡
    │ 3   ┆ 2   │
    └─────┴─────┘

    """
    return F.col(*columns).approx_n_unique()


def first(*columns: str) -> Expr:
    """
    Get the first column or value.

    This function has different behavior depending on the presence of `columns`
    values. If none given (the default), returns an expression that takes the first
    column of the context; otherwise, takes the first value of the given column(s).

    Parameters
    ----------
    *columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "baz"],
    ...     }
    ... )

    Return the first column:

    >>> df.select(pl.first())
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 8   │
    │ 3   │
    └─────┘

    Return the first value for the given column(s):

    >>> df.select(pl.first("b"))
    shape: (1, 1)
    ┌─────┐
    │ b   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 4   │
    └─────┘
    >>> df.select(pl.first("a", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ c   │
    │ --- ┆ --- │
    │ i64 ┆ str │
    ╞═════╪═════╡
    │ 1   ┆ foo │
    └─────┴─────┘

    """
    if not columns:
        return cs.first().as_expr()

    return F.col(*columns).first()


def last(*columns: str) -> Expr:
    """
    Get the last column or value.

    This function has different behavior depending on the presence of `columns`
    values. If none given (the default), returns an expression that takes the last
    column of the context; otherwise, takes the last value of the given column(s).

    Parameters
    ----------
    *columns
        One or more column names.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "baz"],
    ...     }
    ... )

    Return the last column:

    >>> df.select(pl.last())
    shape: (3, 1)
    ┌─────┐
    │ c   │
    │ --- │
    │ str │
    ╞═════╡
    │ foo │
    │ bar │
    │ baz │
    └─────┘

    Return the last value for the given column(s):

    >>> df.select(pl.last("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 3   │
    └─────┘
    >>> df.select(pl.last("b", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ b   ┆ c   │
    │ --- ┆ --- │
    │ i64 ┆ str │
    ╞═════╪═════╡
    │ 2   ┆ baz │
    └─────┴─────┘

    """
    if not columns:
        return cs.last().as_expr()

    return F.col(*columns).last()


def nth(*indices: int | Sequence[int], strict: bool = True) -> Expr:
    """
    Get the nth column(s) of the context.

    Parameters
    ----------
    indices
        One or more indices representing the columns to retrieve.
    strict
        By default, all specified indices must be valid; if any index is out of bounds,
        an error is raised. If set to `False`, out-of-bounds indices are ignored.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "baz"],
    ...     }
    ... )
    >>> df.select(pl.nth(1))
    shape: (3, 1)
    ┌─────┐
    │ b   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 4   │
    │ 5   │
    │ 2   │
    └─────┘
    >>> df.select(pl.nth(2, 0))
    shape: (3, 2)
    ┌─────┬─────┐
    │ c   ┆ a   │
    │ --- ┆ --- │
    │ str ┆ i64 │
    ╞═════╪═════╡
    │ foo ┆ 1   │
    │ bar ┆ 8   │
    │ baz ┆ 3   │
    └─────┴─────┘
    """
    return cs.by_index(*indices, require_all=strict).as_expr()


def head(column: str, n: int = 10) -> Expr:
    """
    Get the first `n` rows.

    This function is syntactic sugar for `pl.col(column).head(n)`.

    Parameters
    ----------
    column
        Column name.
    n
        Number of rows to return.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.head("a"))
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 8   │
    │ 3   │
    └─────┘
    >>> df.select(pl.head("a", 2))
    shape: (2, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 8   │
    └─────┘
    """
    return F.col(column).head(n)


def tail(column: str, n: int = 10) -> Expr:
    """
    Get the last `n` rows.

    This function is syntactic sugar for `pl.col(column).tail(n)`.

    Parameters
    ----------
    column
        Column name.
    n
        Number of rows to return.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.tail("a"))
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 8   │
    │ 3   │
    └─────┘
    >>> df.select(pl.tail("a", 2))
    shape: (2, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 8   │
    │ 3   │
    └─────┘
    """
    return F.col(column).tail(n)


@overload
def corr(
    a: IntoExpr,
    b: IntoExpr,
    *,
    method: CorrelationMethod = ...,
    ddof: int | None = ...,
    propagate_nans: bool = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def corr(
    a: IntoExpr,
    b: IntoExpr,
    *,
    method: CorrelationMethod = ...,
    ddof: int | None = ...,
    propagate_nans: bool = ...,
    eager: Literal[True],
) -> Series: ...


def corr(
    a: IntoExpr,
    b: IntoExpr,
    *,
    method: CorrelationMethod = "pearson",
    ddof: int | None = None,
    propagate_nans: bool = False,
    eager: bool = False,
) -> Expr | Series:
    """
    Compute the Pearson's or Spearman rank correlation between two columns.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    ddof
        Has no effect, do not use.

        .. deprecated:: 1.17.0

    method : {'pearson', 'spearman'}
        Correlation method.
    propagate_nans
        If `True` any `NaN` encountered will lead to `NaN` in the output.
        Defaults to `False` where `NaN` are regarded as larger than any finite number
        and thus lead to the highest rank.
    eager
        Evaluate immediately and return a `Series`; this requires that at least one
        of the given arguments is a `Series`. If set to `False` (default), return
        an expression instead.

    Examples
    --------
    Pearson's correlation:

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.corr("a", "b"))
    shape: (1, 1)
    ┌──────────┐
    │ a        │
    │ ---      │
    │ f64      │
    ╞══════════╡
    │ 0.544705 │
    └──────────┘

    Spearman rank correlation:

    >>> df.select(pl.corr("a", "b", method="spearman"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 0.5 │
    └─────┘

    Eager evaluation:

    >>> s1 = pl.Series("a", [1, 8, 3])
    >>> s2 = pl.Series("b", [4, 5, 2])
    >>> pl.corr(s1, s2, eager=True)
    shape: (1,)
    Series: 'a' [f64]
    [
        0.544705
    ]
    >>> pl.corr(s1, s2, method="spearman", eager=True)
    shape: (1,)
    Series: 'a' [f64]
    [
        0.5
    ]
    """
    if ddof is not None:
        issue_deprecation_warning(
            "the `ddof` parameter has no effect. Do not use it.",
            version="1.17.0",
        )

    if eager:
        if not (isinstance(a, pl.Series) or isinstance(b, pl.Series)):
            msg = "expected at least one Series in 'corr' inputs if 'eager=True'"
            raise ValueError(msg)

        frame = pl.DataFrame([e for e in (a, b) if isinstance(e, pl.Series)])
        exprs = ((e.name if isinstance(e, pl.Series) else e) for e in (a, b))
        return frame.select(
            corr(*exprs, eager=False, method=method, propagate_nans=propagate_nans)
        ).to_series()
    else:
        a_pyexpr = parse_into_expression(a)
        b_pyexpr = parse_into_expression(b)

        if method == "pearson":
            return wrap_expr(plr.pearson_corr(a_pyexpr, b_pyexpr))
        elif method == "spearman":
            return wrap_expr(plr.spearman_rank_corr(a_pyexpr, b_pyexpr, propagate_nans))
        else:
            msg = f"method must be one of {{'pearson', 'spearman'}}, got {method!r}"
            raise ValueError(msg)


@overload
def cov(
    a: IntoExpr,
    b: IntoExpr,
    *,
    ddof: int = ...,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def cov(
    a: IntoExpr,
    b: IntoExpr,
    *,
    ddof: int = ...,
    eager: Literal[True],
) -> Series: ...


def cov(
    a: IntoExpr,
    b: IntoExpr,
    *,
    ddof: int = 1,
    eager: bool = False,
) -> Expr | Series:
    """
    Compute the covariance between two columns/ expressions.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    ddof
        "Delta Degrees of Freedom": the divisor used in the calculation is N - ddof,
        where N represents the number of elements.
        By default ddof is 1.
    eager
        Evaluate immediately and return a `Series`; this requires that at least one
        of the given arguments is a `Series`. If set to `False` (default), return
        an expression instead.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     },
    ... )

    >>> df.select(
    ...     x=pl.cov("a", "b"),
    ...     y=pl.cov("a", "b", ddof=2),
    ... )
    shape: (1, 2)
    ┌─────┬─────┐
    │ x   ┆ y   │
    │ --- ┆ --- │
    │ f64 ┆ f64 │
    ╞═════╪═════╡
    │ 3.0 ┆ 6.0 │
    └─────┴─────┘

    Eager evaluation:

    >>> s1 = pl.Series("a", [1, 8, 3])
    >>> s2 = pl.Series("b", [4, 5, 2])
    >>> pl.cov(s1, s2, eager=True)
    shape: (1,)
    Series: 'a' [f64]
    [
        3.0
    ]
    """
    if eager:
        if not (isinstance(a, pl.Series) or isinstance(b, pl.Series)):
            msg = "expected at least one Series in 'cov' inputs if 'eager=True'"
            raise ValueError(msg)

        frame = pl.DataFrame([e for e in (a, b) if isinstance(e, pl.Series)])
        exprs = ((e.name if isinstance(e, pl.Series) else e) for e in (a, b))
        return frame.select(cov(*exprs, eager=False, ddof=ddof)).to_series()
    else:
        a_pyexpr = parse_into_expression(a)
        b_pyexpr = parse_into_expression(b)
        return wrap_expr(plr.cov(a_pyexpr, b_pyexpr, ddof))


class _map_batches_wrapper:
    def __init__(
        self,
        function: Callable[[Sequence[Series]], Series | Any],
        *,
        returns_scalar: bool,
    ) -> None:
        self.function = function
        self.returns_scalar = returns_scalar

    def __call__(
        self, sl: list[plr.PySeries], *args: Any, **kwargs: Any
    ) -> plr.PySeries:
        return_dtype = kwargs["return_dtype"]
        slp = [wrap_s(s) for s in sl]

        # ufunc and numba don't expect return_dtype
        try:
            rv = self.function(slp, *args, **kwargs)
        except TypeError as e:
            if "unexpected keyword argument 'return_dtype'" in e.args[0]:
                kwargs.pop("return_dtype")
                rv = self.function(slp, *args, **kwargs)
            else:
                raise

        if _check_for_numpy(rv) and isinstance(rv, np.ndarray):
            rv = pl.Series(rv, dtype=return_dtype)

        if isinstance(rv, pl.Series):
            return rv._s
        elif self.returns_scalar:
            return pl.Series([rv], dtype=return_dtype)._s
        else:
            msg = f"`map` with `returns_scalar=False` must return a Series; found {qualified_type_name(rv)!r}.\n\nIf `returns_scalar` is set to `True`, a returned value can be a scalar value."
            raise TypeError(msg)


def map_batches(
    exprs: Sequence[str | Expr],
    function: Callable[[Sequence[Series]], Series | Any],
    return_dtype: PolarsDataType | pl.DataTypeExpr | None = None,
    *,
    is_elementwise: bool = False,
    returns_scalar: bool = False,
) -> Expr:
    """
    Map a custom function over multiple columns/expressions.

    Produces a single Series result.

    .. warning::
        This method is much slower than the native expressions API.
        Only use it if you cannot implement your logic otherwise.

    Parameters
    ----------
    exprs
        Expression(s) representing the input Series to the function.
    function
        Function to apply over the input.
    return_dtype
        Datatype of the output Series.

        It is recommended to set this whenever possible. If this is `None`, it tries
        to infer the datatype by calling the function with dummy data and looking at
        the output.
    is_elementwise
        Set to true if the operations is elementwise for better performance
        and optimization.

        An elementwise operations has unit or equal length for all inputs
        and can be ran sequentially on slices without results being affected.
    returns_scalar
        If the function returns a scalar, by default it will be wrapped in
        a list in the output, since the assumption is that the function
        always returns something Series-like. If you want to keep the
        result as a scalar, set this argument to True.

    Notes
    -----
    A UDF passed to `map_batches` must be pure, meaning that it cannot modify
    or depend on state other than its arguments. We may call the function
    with arbitrary input data.

    Returns
    -------
    Expr
        Expression with the data type given by `return_dtype`.

    Examples
    --------
    >>> def test_func(a, b, c):
    ...     return a + b + c
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3, 4],
    ...         "b": [4, 5, 6, 7],
    ...     }
    ... )
    >>>
    >>> df.with_columns(
    ...     (
    ...         pl.struct(["a", "b"]).map_batches(
    ...             lambda x: test_func(x.struct.field("a"), x.struct.field("b"), 1)
    ...         )
    ...     ).alias("a+b+c")
    ... )
    shape: (4, 3)
    ┌─────┬─────┬───────┐
    │ a   ┆ b   ┆ a+b+c │
    │ --- ┆ --- ┆ ---   │
    │ i64 ┆ i64 ┆ i64   │
    ╞═════╪═════╪═══════╡
    │ 1   ┆ 4   ┆ 6     │
    │ 2   ┆ 5   ┆ 8     │
    │ 3   ┆ 6   ┆ 10    │
    │ 4   ┆ 7   ┆ 12    │
    └─────┴─────┴───────┘
    """
    pyexprs = parse_into_list_of_expressions(exprs)

    return_dtype_expr = (
        parse_into_datatype_expr(return_dtype)._pydatatype_expr
        if return_dtype is not None
        else None
    )

    return wrap_expr(
        plr.map_expr(
            pyexprs,
            _map_batches_wrapper(function, returns_scalar=returns_scalar),
            return_dtype_expr,
            is_elementwise=is_elementwise,
            returns_scalar=returns_scalar,
        )
    )


def map_groups(
    exprs: Sequence[str | Expr],
    function: Callable[[Sequence[Series]], Series | Any],
    return_dtype: PolarsDataType | pl.DataTypeExpr | None = None,
    *,
    is_elementwise: bool = False,
    returns_scalar: bool = False,
) -> Expr:
    """
    Apply a custom/user-defined function (UDF) in a GroupBy context.

    .. warning::
        This method is much slower than the native expressions API.
        Only use it if you cannot implement your logic otherwise.

    Parameters
    ----------
    exprs
        Expression(s) representing the input Series to the function.
    function
        Function to apply over the input; should be of type Callable[[Series], Series].
    return_dtype
        Datatype of the output Series.

        It is recommended to set this whenever possible. If this is `None`, it tries
        to infer the datatype by calling the function with dummy data and looking at
        the output.
    is_elementwise
        Set to true if the operations is elementwise for better performance
        and optimization.

        An elementwise operations has unit or equal length for all inputs
        and can be ran sequentially on slices without results being affected.
    returns_scalar
        If the function returns a single scalar as output.

    Notes
    -----
    A UDF passed to `map_batches` must be pure, meaning that it cannot modify
    or depend on state other than its arguments. Polars may call the function
    with arbitrary input data.

    Returns
    -------
    Expr
        Expression with the data type given by `return_dtype`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "group": [1, 1, 2],
    ...         "a": [1, 3, 3],
    ...         "b": [5, 6, 7],
    ...     }
    ... )
    >>> df
    shape: (3, 3)
    ┌───────┬─────┬─────┐
    │ group ┆ a   ┆ b   │
    │ ---   ┆ --- ┆ --- │
    │ i64   ┆ i64 ┆ i64 │
    ╞═══════╪═════╪═════╡
    │ 1     ┆ 1   ┆ 5   │
    │ 1     ┆ 3   ┆ 6   │
    │ 2     ┆ 3   ┆ 7   │
    └───────┴─────┴─────┘
    >>> (
    ...     df.group_by("group").agg(
    ...         pl.map_groups(
    ...             exprs=["a", "b"],
    ...             function=lambda list_of_series: list_of_series[0]
    ...             / list_of_series[0].sum()
    ...             + list_of_series[1],
    ...             return_dtype=pl.Float64,
    ...         ).alias("my_custom_aggregation")
    ...     )
    ... ).sort("group")
    shape: (2, 2)
    ┌───────┬───────────────────────┐
    │ group ┆ my_custom_aggregation │
    │ ---   ┆ ---                   │
    │ i64   ┆ list[f64]             │
    ╞═══════╪═══════════════════════╡
    │ 1     ┆ [5.25, 6.75]          │
    │ 2     ┆ [8.0]                 │
    └───────┴───────────────────────┘

    The output for group `1` can be understood as follows:

    - group `1` contains Series `'a': [1, 3]` and `'b': [5, 6]`
    - applying the function to those lists of Series, one gets the output
      `[1 / 4 + 5, 3 / 4 + 6]`, i.e. `[5.25, 6.75]`
    """
    return map_batches(
        exprs,
        function,
        return_dtype,
        is_elementwise=is_elementwise,
        returns_scalar=returns_scalar,
    )


def _row_encode(
    exprs: pl.Selector | pl.Expr | Sequence[str | pl.Expr],
    *,
    unordered: bool = False,
    descending: list[bool] | None = None,
    nulls_last: list[bool] | None = None,
) -> Expr:
    if isinstance(exprs, pl.Selector):
        exprs = [exprs.as_expr()]
    elif isinstance(exprs, pl.Expr):
        exprs = [exprs]

    pyexprs = parse_into_list_of_expressions(exprs)

    if unordered:
        assert descending is None
        assert nulls_last is None

        result = plr.PyExpr.row_encode_unordered(pyexprs)
    else:
        result = plr.PyExpr.row_encode_ordered(pyexprs, descending, nulls_last)

    return wrap_expr(result)


def _wrap_acc_lamba(
    function: Callable[[Series, Series], Series],
) -> Callable[[tuple[plr.PySeries, plr.PySeries]], plr.PySeries]:
    def wrapper(t: tuple[plr.PySeries, plr.PySeries]) -> plr.PySeries:
        a, b = t
        return function(wrap_s(a), wrap_s(b))._s

    return wrapper


def fold(
    acc: IntoExpr,
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
    *,
    returns_scalar: bool = False,
    return_dtype: pl.DataTypeExpr | PolarsDataType | None = None,
) -> Expr:
    """
    Accumulate over multiple columns horizontally/ row wise with a left fold.

    Parameters
    ----------
    acc
        Accumulator Expression. This is the value that will be initialized when the fold
        starts. For a sum this could for instance be lit(0).
    function
        Function to apply over the accumulator and the value.
        Fn(acc, value) -> new_value
    exprs
        Expressions to aggregate over. May also be a wildcard expression.
    returns_scalar
        Whether or not `function` applied returns a scalar. This must be set correctly
        by the user.
    return_dtype
        Output datatype.
        If not set, the dtype will be inferred based on the dtype
        of the accumulator.

    Notes
    -----
    If you simply want the first encountered expression as accumulator,
    consider using `reduce`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [3, 4, 5],
    ...         "c": [5, 6, 7],
    ...     }
    ... )
    >>> df
    shape: (3, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 5   │
    │ 2   ┆ 4   ┆ 6   │
    │ 3   ┆ 5   ┆ 7   │
    └─────┴─────┴─────┘

    Horizontally sum over all columns and add 1.

    >>> df.select(
    ...     pl.fold(
    ...         acc=pl.lit(1), function=lambda acc, x: acc + x, exprs=pl.col("*")
    ...     ).alias("sum"),
    ... )
    shape: (3, 1)
    ┌─────┐
    │ sum │
    │ --- │
    │ i32 │
    ╞═════╡
    │ 10  │
    │ 13  │
    │ 16  │
    └─────┘

    You can also apply a condition/predicate on all columns:

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [0, 1, 2],
    ...     }
    ... )
    >>> df
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 0   │
    │ 2   ┆ 1   │
    │ 3   ┆ 2   │
    └─────┴─────┘

    >>> df.filter(
    ...     pl.fold(
    ...         acc=pl.lit(True),
    ...         function=lambda acc, x: acc & x,
    ...         exprs=pl.col("*") > 1,
    ...     )
    ... )
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 3   ┆ 2   │
    └─────┴─────┘
    """
    # in case of col("*")
    pyacc = parse_into_expression(acc, str_as_lit=True)
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    rt: plr.PyDataTypeExpr | None = None
    if return_dtype is not None:
        rt = parse_into_datatype_expr(return_dtype)._pydatatype_expr

    pyexprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.fold(
            pyacc,
            _wrap_acc_lamba(function),
            pyexprs,
            returns_scalar=returns_scalar,
            return_dtype=rt,
        )
    )


def reduce(
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
    *,
    returns_scalar: bool = False,
    return_dtype: pl.DataTypeExpr | PolarsDataType | None = None,
) -> Expr:
    """
    Accumulate over multiple columns horizontally/ row wise with a left fold.

    Parameters
    ----------
    function
        Function to apply over the accumulator and the value.
        Fn(acc, value) -> new_value
    exprs
        Expressions to aggregate over. May also be a wildcard expression.
    returns_scalar
        Whether or not `function` applied returns a scalar. This must be set correctly
        by the user.
    return_dtype
        Output datatype.
        If not set, the dtype will be inferred based on the dtype of the input
        expressions.

    Notes
    -----
    See `fold` for the version with an explicit accumulator.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [0, 1, 2],
    ...     }
    ... )
    >>> df
    shape: (3, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 0   │
    │ 2   ┆ 1   │
    │ 3   ┆ 2   │
    └─────┴─────┘

    Horizontally sum over all columns.

    >>> df.select(
    ...     pl.reduce(function=lambda acc, x: acc + x, exprs=pl.col("*")).alias("sum")
    ... )
    shape: (3, 1)
    ┌─────┐
    │ sum │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 3   │
    │ 5   │
    └─────┘
    """
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    rt: plr.PyDataTypeExpr | None = None
    if return_dtype is not None:
        rt = parse_into_datatype_expr(return_dtype)._pydatatype_expr

    pyexprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.reduce(
            _wrap_acc_lamba(function),
            pyexprs,
            returns_scalar=returns_scalar,
            return_dtype=rt,
        )
    )


def cum_fold(
    acc: IntoExpr,
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
    *,
    returns_scalar: bool = False,
    return_dtype: pl.DataTypeExpr | PolarsDataType | None = None,
    include_init: bool = False,
) -> Expr:
    """
    Cumulatively fold horizontally across columns with a left fold.

    Every cumulative result is added as a separate field in a Struct column.

    Parameters
    ----------
    acc
        Accumulator expression. This is the value that will be initialized when the fold
        starts. For a sum this could for instance be lit(0).
    function
        Function to apply over the accumulator and the value.
        Fn(acc, value) -> new_value
    exprs
        Expressions to aggregate over. May also be a wildcard expression.
    returns_scalar
        Whether or not `function` applied returns a scalar. This must be set correctly
        by the user.
    return_dtype
        Output datatype.
        If not set, the dtype will be inferred based on the dtype of the accumulator.
    include_init
        Include the initial accumulator state as struct field.

    Notes
    -----
    If you simply want the first encountered expression as accumulator,
    consider using :func:`cum_reduce`.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [3, 4, 5],
    ...         "c": [5, 6, 7],
    ...     }
    ... )
    >>> df.with_columns(
    ...     pl.cum_fold(acc=pl.lit(1), function=lambda acc, x: acc + x, exprs=pl.all())
    ... )
    shape: (3, 4)
    ┌─────┬─────┬─────┬───────────┐
    │ a   ┆ b   ┆ c   ┆ cum_fold  │
    │ --- ┆ --- ┆ --- ┆ ---       │
    │ i64 ┆ i64 ┆ i64 ┆ struct[3] │
    ╞═════╪═════╪═════╪═══════════╡
    │ 1   ┆ 3   ┆ 5   ┆ {2,5,10}  │
    │ 2   ┆ 4   ┆ 6   ┆ {3,7,13}  │
    │ 3   ┆ 5   ┆ 7   ┆ {4,9,16}  │
    └─────┴─────┴─────┴───────────┘
    """
    # in case of col("*")
    pyacc = parse_into_expression(acc, str_as_lit=True)
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    rt: plr.PyDataTypeExpr | None = None
    if return_dtype is not None:
        rt = parse_into_datatype_expr(return_dtype)._pydatatype_expr

    pyexprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.cum_fold(
            pyacc,
            _wrap_acc_lamba(function),
            pyexprs,
            returns_scalar=returns_scalar,
            return_dtype=rt,
            include_init=include_init,
        ).alias("cum_fold")
    )


def cum_reduce(
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
    *,
    returns_scalar: bool = False,
    return_dtype: pl.DataTypeExpr | PolarsDataType | None = None,
) -> Expr:
    """
    Cumulatively reduce horizontally across columns with a left fold.

    Every cumulative result is added as a separate field in a Struct column.

    Parameters
    ----------
    function
        Function to apply over the accumulator and the value.
        Fn(acc, value) -> new_value
    exprs
        Expressions to aggregate over. May also be a wildcard expression.
    returns_scalar
        Whether or not `function` applied returns a scalar. This must be set correctly
        by the user.
    return_dtype
        Output datatype.
        If not set, the dtype will be inferred based on the dtype of the input
        expressions.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [3, 4, 5],
    ...         "c": [5, 6, 7],
    ...     }
    ... )
    >>> df.with_columns(pl.cum_reduce(function=lambda acc, x: acc + x, exprs=pl.all()))
    shape: (3, 4)
    ┌─────┬─────┬─────┬────────────┐
    │ a   ┆ b   ┆ c   ┆ cum_reduce │
    │ --- ┆ --- ┆ --- ┆ ---        │
    │ i64 ┆ i64 ┆ i64 ┆ struct[3]  │
    ╞═════╪═════╪═════╪════════════╡
    │ 1   ┆ 3   ┆ 5   ┆ {1,4,9}    │
    │ 2   ┆ 4   ┆ 6   ┆ {2,6,12}   │
    │ 3   ┆ 5   ┆ 7   ┆ {3,8,15}   │
    └─────┴─────┴─────┴────────────┘
    """
    # in case of col("*")
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    rt: plr.PyDataTypeExpr | None = None
    if return_dtype is not None:
        rt = parse_into_datatype_expr(return_dtype)._pydatatype_expr

    pyexprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.cum_reduce(
            _wrap_acc_lamba(function),
            pyexprs,
            returns_scalar=returns_scalar,
            return_dtype=rt,
        ).alias("cum_reduce")
    )


def arctan2(y: str | Expr, x: str | Expr) -> Expr:
    """
    Compute two argument arctan in radians.

    Returns the angle (in radians) in the plane between the
    positive x-axis and the ray from the origin to (x,y).

    Parameters
    ----------
    y
        Column name or Expression.
    x
        Column name or Expression.

    Examples
    --------
    >>> c = (2**0.5) / 2
    >>> df = pl.DataFrame(
    ...     {
    ...         "y": [c, -c, c, -c],
    ...         "x": [c, c, -c, -c],
    ...     }
    ... )
    >>> df.with_columns(pl.arctan2("y", "x").alias("atan2"))
    shape: (4, 3)
    ┌───────────┬───────────┬───────────┐
    │ y         ┆ x         ┆ atan2     │
    │ ---       ┆ ---       ┆ ---       │
    │ f64       ┆ f64       ┆ f64       │
    ╞═══════════╪═══════════╪═══════════╡
    │ 0.707107  ┆ 0.707107  ┆ 0.785398  │
    │ -0.707107 ┆ 0.707107  ┆ -0.785398 │
    │ 0.707107  ┆ -0.707107 ┆ 2.356194  │
    │ -0.707107 ┆ -0.707107 ┆ -2.356194 │
    └───────────┴───────────┴───────────┘
    """
    if isinstance(y, str):
        y = F.col(y)
    if isinstance(x, str):
        x = F.col(x)
    if not hasattr(x, "_pyexpr"):
        msg = f"`arctan2` expected a `str` or `Expr` got a `{qualified_type_name(x)}`"
        raise TypeError(msg)
    if not hasattr(y, "_pyexpr"):
        msg = f"`arctan2` expected a `str` or `Expr` got a `{qualified_type_name(y)}`"
        raise TypeError(msg)

    return wrap_expr(plr.arctan2(y._pyexpr, x._pyexpr))


@deprecated("`arctan2d` is deprecated; use `arctan2` followed by `.degrees()` instead.")
def arctan2d(y: str | Expr, x: str | Expr) -> Expr:
    """
    Compute two argument arctan in degrees.

    .. deprecated:: 1.0.0
        Use `arctan2` followed by :meth:`Expr.degrees` instead.

    Returns the angle (in degrees) in the plane between the positive x-axis
    and the ray from the origin to (x,y).

    Parameters
    ----------
    y
        Column name or Expression.
    x
        Column name or Expression.

    Examples
    --------
    >>> c = (2**0.5) / 2
    >>> df = pl.DataFrame(
    ...     {
    ...         "y": [c, -c, c, -c],
    ...         "x": [c, c, -c, -c],
    ...     }
    ... )
    >>> df.select(  # doctest: +SKIP
    ...     pl.arctan2d("y", "x").alias("atan2d"),
    ...     pl.arctan2("y", "x").alias("atan2"),
    ... )
    shape: (4, 2)
    ┌────────┬───────────┐
    │ atan2d ┆ atan2     │
    │ ---    ┆ ---       │
    │ f64    ┆ f64       │
    ╞════════╪═══════════╡
    │ 45.0   ┆ 0.785398  │
    │ -45.0  ┆ -0.785398 │
    │ 135.0  ┆ 2.356194  │
    │ -135.0 ┆ -2.356194 │
    └────────┴───────────┘
    """
    return arctan2(y, x).degrees()


def exclude(
    columns: str | PolarsDataType | Collection[str] | Collection[PolarsDataType],
    *more_columns: str | PolarsDataType,
) -> Expr:
    """
    Represent all columns except for the given columns.

    Syntactic sugar for `pl.all().exclude(columns)`.

    Parameters
    ----------
    columns
        The name or datatype of the column(s) to exclude. Accepts regular expression
        input. Regular expressions should start with `^` and end with `$`.
    *more_columns
        Additional names or datatypes of columns to exclude, specified as positional
        arguments.

    Examples
    --------
    Exclude by column name(s):

    >>> df = pl.DataFrame(
    ...     {
    ...         "aa": [1, 2, 3],
    ...         "ba": ["a", "b", None],
    ...         "cc": [None, 2.5, 1.5],
    ...     }
    ... )
    >>> df.select(pl.exclude("ba"))
    shape: (3, 2)
    ┌─────┬──────┐
    │ aa  ┆ cc   │
    │ --- ┆ ---  │
    │ i64 ┆ f64  │
    ╞═════╪══════╡
    │ 1   ┆ null │
    │ 2   ┆ 2.5  │
    │ 3   ┆ 1.5  │
    └─────┴──────┘

    Exclude by regex, e.g. removing all columns whose names end with the letter "a":

    >>> df.select(pl.exclude("^.*a$"))
    shape: (3, 1)
    ┌──────┐
    │ cc   │
    │ ---  │
    │ f64  │
    ╞══════╡
    │ null │
    │ 2.5  │
    │ 1.5  │
    └──────┘

    Exclude by dtype(s), e.g. removing all columns of type Int64 or Float64:

    >>> df.select(pl.exclude([pl.Int64, pl.Float64]))
    shape: (3, 1)
    ┌──────┐
    │ ba   │
    │ ---  │
    │ str  │
    ╞══════╡
    │ a    │
    │ b    │
    │ null │
    └──────┘

    """
    return F.col("*").exclude(columns, *more_columns)


def groups(column: str) -> Expr:
    """
    Syntactic sugar for `pl.col("foo").agg_groups()`.

    .. deprecated:: 1.35
        Use `df.with_row_index().group_by(...).agg(pl.col('index'))` instead.
        This method will be removed in Polars 2.0.
    """
    warnings.warn(
        "pl.groups() is deprecated and will be removed in Polars 2.0. "
        "Use df.with_row_index().group_by(...).agg(pl.col('index')) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return F.col(column).agg_groups()


def quantile(
    column: str,
    quantile: float | Expr,
    interpolation: QuantileMethod = "nearest",
) -> Expr:
    """
    Syntactic sugar for `pl.col("foo").quantile(..)`.

    Parameters
    ----------
    column
        Column name.
    quantile
        Quantile between 0.0 and 1.0.
    interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
        Interpolation method.
    """
    return F.col(column).quantile(quantile, interpolation)


def arg_sort_by(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    descending: bool | Sequence[bool] = False,
    nulls_last: bool | Sequence[bool] = False,
    multithreaded: bool = True,
    maintain_order: bool = False,
) -> Expr:
    """
    Return the row indices that would sort the column(s).

    Parameters
    ----------
    exprs
        Column(s) to arg sort by. Accepts expression input. Strings are parsed as column
        names.
    *more_exprs
        Additional columns to arg sort by, specified as positional arguments.
    descending
        Sort in descending order. When sorting by multiple columns, can be specified
        per column by passing a sequence of booleans.
    nulls_last
        Place null values last.
    multithreaded
        Sort using multiple threads.
    maintain_order
        Whether the order should be maintained if elements are equal.

    See Also
    --------
    Expr.gather: Take values by index.
    Expr.rank : Get the rank of each row.

    Examples
    --------
    Pass a single column name to compute the arg sort by that column.

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [0, 1, 1, 0],
    ...         "b": [3, 2, 3, 2],
    ...         "c": [1, 2, 3, 4],
    ...     }
    ... )
    >>> df.select(pl.arg_sort_by("a"))
    shape: (4, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 0   │
    │ 3   │
    │ 1   │
    │ 2   │
    └─────┘

    Compute the arg sort by multiple columns by either passing a list of columns, or by
    specifying each column as a positional argument.

    >>> df.select(pl.arg_sort_by(["a", "b"], descending=True))
    shape: (4, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 2   │
    │ 1   │
    │ 0   │
    │ 3   │
    └─────┘

    Use gather to apply the arg sort to other columns.

    >>> df.select(pl.col("c").gather(pl.arg_sort_by("a")))
    shape: (4, 1)
    ┌─────┐
    │ c   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 4   │
    │ 2   │
    │ 3   │
    └─────┘
    """
    exprs = parse_into_list_of_expressions(exprs, *more_exprs)
    descending = extend_bool(descending, len(exprs), "descending", "exprs")
    nulls_last = extend_bool(nulls_last, len(exprs), "nulls_last", "exprs")
    return wrap_expr(
        plr.arg_sort_by(exprs, descending, nulls_last, multithreaded, maintain_order)
    )


@deprecate_streaming_parameter()
@forward_old_opt_flags()
def collect_all(
    lazy_frames: Iterable[LazyFrame],
    *,
    type_coercion: bool = True,
    predicate_pushdown: bool = True,
    projection_pushdown: bool = True,
    simplify_expression: bool = True,
    no_optimization: bool = False,
    slice_pushdown: bool = True,
    comm_subplan_elim: bool = True,
    comm_subexpr_elim: bool = True,
    cluster_with_columns: bool = True,
    collapse_joins: bool = True,
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
    engine: EngineType = "auto",
) -> list[DataFrame]:
    """
    Collect multiple LazyFrames at the same time.

    This can run all the computation graphs in parallel or combined.

    Common Subplan Elimination is applied on the combined plan, meaning
    that diverging queries will run only once.

    Parameters
    ----------
    lazy_frames
        A list of LazyFrames to collect.
    type_coercion
        Do type coercion optimization.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    predicate_pushdown
        Do predicate pushdown optimization.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    projection_pushdown
        Do projection pushdown optimization.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    simplify_expression
        Run simplify expressions optimization.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    no_optimization
        Turn off optimizations.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    slice_pushdown
        Slice pushdown optimization.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    comm_subplan_elim
        Will try to cache branching subplans that occur on self-joins or unions.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    comm_subexpr_elim
        Common subexpressions will be cached and reused.

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    cluster_with_columns
        Combine sequential independent calls to with_columns

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    collapse_joins
        Collapse a join and filters into a faster join

        .. deprecated:: 1.30.0
            Use the `optimizations` parameters.
    optimizations
        The optimization passes done during query optimization.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    engine
        Select the engine used to process the query, optional.
        At the moment, if set to `"auto"` (default), the query
        is run using the polars in-memory engine. Polars will also
        attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
        environment variable. If it cannot run the query using the
        selected engine, the query is run using the polars in-memory
        engine.

        .. note::
           The GPU engine does not support async, or running in the
           background. If either are enabled, then GPU execution is switched off.

    Returns
    -------
    list of DataFrames
        The collected DataFrames, returned in the same order as the input LazyFrames.

    """
    if engine == "streaming":
        issue_unstable_warning("streaming mode is considered unstable.")

    lfs = [lf._ldf for lf in lazy_frames]
    out = plr.collect_all(lfs, engine, optimizations._pyoptflags)

    # wrap the pydataframes into dataframe
    result = [wrap_df(pydf) for pydf in out]

    return result


@overload
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: Literal[True],
    engine: EngineType = "auto",
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
) -> _GeventDataFrameResult[list[DataFrame]]: ...


@overload
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: Literal[False] = False,
    engine: EngineType = "auto",
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
) -> Awaitable[list[DataFrame]]: ...


@unstable()
@deprecate_streaming_parameter()
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: bool = False,
    engine: EngineType = "auto",
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
) -> Awaitable[list[DataFrame]] | _GeventDataFrameResult[list[DataFrame]]:
    """
    Collect multiple LazyFrames at the same time asynchronously in thread pool.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    Collects into a list of DataFrame (like :func:`polars.collect_all`),
    but instead of returning them directly, they are scheduled to be collected
    inside thread pool, while this method returns almost instantly.

    May be useful if you use gevent or asyncio and want to release control to other
    greenlets/tasks while LazyFrames are being collected.

    Parameters
    ----------
    lazy_frames
        A list of LazyFrames to collect.
    gevent
        Return wrapper to `gevent.event.AsyncResult` instead of Awaitable
    optimizations
        The optimization passes done during query optimization.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    engine
        Select the engine used to process the query, optional.
        At the moment, if set to `"auto"` (default), the query
        is run using the polars in-memory engine. Polars will also
        attempt to use the engine set by the `POLARS_ENGINE_AFFINITY`
        environment variable. If it cannot run the query using the
        selected engine, the query is run using the polars in-memory
        engine.

        .. note::
           The GPU engine does not support async, or running in the
           background. If either are enabled, then GPU execution is switched off.

    See Also
    --------
    polars.collect_all : Collect multiple LazyFrames at the same time.
    LazyFrame.collect_async : To collect single frame.

    Notes
    -----
    In case of error `set_exception` is used on
    `asyncio.Future`/`gevent.event.AsyncResult` and will be reraised by them.

    Returns
    -------
    If `gevent=False` (default) then returns awaitable.

    If `gevent=True` then returns wrapper that has
    `.get(block=True, timeout=None)` method.
    """
    if engine == "streaming":
        issue_unstable_warning("streaming mode is considered unstable.")

    result: (
        _GeventDataFrameResult[list[DataFrame]] | _AioDataFrameResult[list[DataFrame]]
    ) = _GeventDataFrameResult() if gevent else _AioDataFrameResult()
    lfs = [lf._ldf for lf in lazy_frames]
    plr.collect_all_with_callback(
        lfs, engine, optimizations._pyoptflags, result._callback_all
    )
    return result


@unstable()
def explain_all(
    lazy_frames: Iterable[LazyFrame],
    *,
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
) -> str:
    """
    Explain multiple LazyFrames as if passed to `collect_all`.

    Common Subplan Elimination is applied on the combined plan, meaning
    that diverging queries will run only once.

    Parameters
    ----------
    lazy_frames
        A list of LazyFrames to collect.
    optimizations
        The optimization passes done during query optimization.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

    Returns
    -------
    Explained plan.
    """
    lfs = [lf._ldf for lf in lazy_frames]
    return plr.explain_all(lfs, optimizations._pyoptflags)


@overload
def select(
    *exprs: IntoExpr | Iterable[IntoExpr],
    eager: Literal[True] = ...,
    **named_exprs: IntoExpr,
) -> DataFrame: ...


@overload
def select(
    *exprs: IntoExpr | Iterable[IntoExpr],
    eager: Literal[False],
    **named_exprs: IntoExpr,
) -> LazyFrame: ...


def select(
    *exprs: IntoExpr | Iterable[IntoExpr], eager: bool = True, **named_exprs: IntoExpr
) -> DataFrame | LazyFrame:
    """
    Run polars expressions without a context.

    This is syntactic sugar for running `df.select` on an empty DataFrame
    (or LazyFrame if eager=False).

    Parameters
    ----------
    *exprs
        Column(s) to select, specified as positional arguments.
        Accepts expression input. Strings are parsed as column names,
        other non-expression inputs are parsed as literals.
    eager
        Evaluate immediately and return a `DataFrame` (default); if set to `False`,
        return a `LazyFrame` instead.
    **named_exprs
        Additional columns to select, specified as keyword arguments.
        The columns will be renamed to the keyword used.

    Returns
    -------
    DataFrame or LazyFrame

    Examples
    --------
    >>> foo = pl.Series("foo", [1, 2, 3])
    >>> bar = pl.Series("bar", [3, 2, 1])
    >>> pl.select(min=pl.min_horizontal(foo, bar))
    shape: (3, 1)
    ┌─────┐
    │ min │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 2   │
    │ 1   │
    └─────┘

    >>> pl.select(pl.int_range(0, 100_000, 2).alias("n"), eager=False).filter(
    ...     pl.col("n") % 22_500 == 0
    ... ).collect()
    shape: (5, 1)
    ┌───────┐
    │ n     │
    │ ---   │
    │ i64   │
    ╞═══════╡
    │ 0     │
    │ 22500 │
    │ 45000 │
    │ 67500 │
    │ 90000 │
    └───────┘
    """
    empty_frame = pl.DataFrame() if eager else pl.LazyFrame()
    return empty_frame.select(*exprs, **named_exprs)


@overload
def arg_where(condition: Expr | Series, *, eager: Literal[False] = ...) -> Expr: ...


@overload
def arg_where(condition: Expr | Series, *, eager: Literal[True]) -> Series: ...


def arg_where(condition: Expr | Series, *, eager: bool = False) -> Expr | Series:
    """
    Return indices where `condition` evaluates `True`.

    Parameters
    ----------
    condition
        Boolean expression to evaluate
    eager
        Evaluate immediately and return a `Series`; this requires that the given
        condition is itself a `Series`. If set to `False` (default), return
        an expression instead.

    See Also
    --------
    Series.arg_true : Return indices where Series is True

    Examples
    --------
    >>> df = pl.DataFrame({"a": [1, 2, 3, 4, 5]})
    >>> df.select(
    ...     [
    ...         pl.arg_where(pl.col("a") % 2 == 0),
    ...     ]
    ... ).to_series()
    shape: (2,)
    Series: 'a' [u32]
    [
        1
        3
    ]
    """
    if eager:
        if not isinstance(condition, pl.Series):
            msg = (
                "expected Series in 'arg_where' if 'eager=True', got"
                f" {type(condition).__name__!r}"
            )
            raise ValueError(msg)
        return condition.to_frame().select(arg_where(F.col(condition.name))).to_series()
    else:
        condition_pyexpr = parse_into_expression(condition)
        return wrap_expr(plr.arg_where(condition_pyexpr))


@overload
def coalesce(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    eager: Literal[False] = ...,
) -> Expr: ...


@overload
def coalesce(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    eager: Literal[True],
) -> Series: ...


@overload
def coalesce(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    eager: bool,
) -> Expr | Series: ...


def coalesce(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    eager: bool = False,
) -> Expr | Series:
    """
    Folds the columns from left to right, keeping the first non-null value.

    Parameters
    ----------
    exprs
        Columns to coalesce. Accepts expression input. Strings are parsed as column
        names, other non-expression inputs are parsed as literals.
    *more_exprs
        Additional columns to coalesce, specified as positional arguments.
    eager
        Evaluate immediately and return a `Series`; this requires that at least one
        of the given arguments is a `Series`. If set to `False` (default), return
        an expression instead.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, None, None, None],
    ...         "b": [1, 2, None, None],
    ...         "c": [5, None, 3, None],
    ...     }
    ... )

    >>> df.with_columns(pl.coalesce("a", "b", "c", 10).alias("d"))
    shape: (4, 4)
    ┌──────┬──────┬──────┬─────┐
    │ a    ┆ b    ┆ c    ┆ d   │
    │ ---  ┆ ---  ┆ ---  ┆ --- │
    │ i64  ┆ i64  ┆ i64  ┆ i64 │
    ╞══════╪══════╪══════╪═════╡
    │ 1    ┆ 1    ┆ 5    ┆ 1   │
    │ null ┆ 2    ┆ null ┆ 2   │
    │ null ┆ null ┆ 3    ┆ 3   │
    │ null ┆ null ┆ null ┆ 10  │
    └──────┴──────┴──────┴─────┘

    >>> df.with_columns(pl.coalesce(pl.col(["a", "b", "c"]), 10.0).alias("d"))
    shape: (4, 4)
    ┌──────┬──────┬──────┬──────┐
    │ a    ┆ b    ┆ c    ┆ d    │
    │ ---  ┆ ---  ┆ ---  ┆ ---  │
    │ i64  ┆ i64  ┆ i64  ┆ f64  │
    ╞══════╪══════╪══════╪══════╡
    │ 1    ┆ 1    ┆ 5    ┆ 1.0  │
    │ null ┆ 2    ┆ null ┆ 2.0  │
    │ null ┆ null ┆ 3    ┆ 3.0  │
    │ null ┆ null ┆ null ┆ 10.0 │
    └──────┴──────┴──────┴──────┘

    >>> s1 = pl.Series("a", [None, 2, None])
    >>> s2 = pl.Series("b", [1, None, 3])
    >>> pl.coalesce(s1, s2, eager=True)
    shape: (3,)
    Series: 'a' [i64]
    [
        1
        2
        3
    ]
    """
    if eager:
        exprs = [exprs, *more_exprs]
        if not (series := [e for e in exprs if isinstance(e, pl.Series)]):
            msg = "expected at least one Series in 'coalesce' if 'eager=True'"
            raise ValueError(msg)

        exprs = [(e.name if isinstance(e, pl.Series) else e) for e in exprs]
        return pl.DataFrame(series).select(coalesce(exprs, eager=False)).to_series()
    else:
        exprs = parse_into_list_of_expressions(exprs, *more_exprs)
        return wrap_expr(plr.coalesce(exprs))


@overload
def from_epoch(column: str | Expr, time_unit: EpochTimeUnit = ...) -> Expr: ...


@overload
def from_epoch(
    column: Series | Sequence[int], time_unit: EpochTimeUnit = ...
) -> Series: ...


def from_epoch(
    column: str | Expr | Series | Sequence[int], time_unit: EpochTimeUnit = "s"
) -> Expr | Series:
    """
    Utility function that parses an epoch timestamp (or Unix time) to Polars Date(time).

    Depending on the `time_unit` provided, this function will return a different dtype:

    - time_unit="d" returns pl.Date
    - time_unit="s" returns pl.Datetime["us"] (pl.Datetime's default)
    - time_unit="ms" returns pl.Datetime["ms"]
    - time_unit="us" returns pl.Datetime["us"]
    - time_unit="ns" returns pl.Datetime["ns"]

    Parameters
    ----------
    column
        Series or expression to parse integers to pl.Datetime.
    time_unit
        The unit of time of the timesteps since epoch time.

    Examples
    --------
    >>> df = pl.DataFrame({"timestamp": [1666683077, 1666683099]}).lazy()
    >>> df.select(pl.from_epoch(pl.col("timestamp"), time_unit="s")).collect()
    shape: (2, 1)
    ┌─────────────────────┐
    │ timestamp           │
    │ ---                 │
    │ datetime[μs]        │
    ╞═════════════════════╡
    │ 2022-10-25 07:31:17 │
    │ 2022-10-25 07:31:39 │
    └─────────────────────┘

    The function can also be used in an eager context by passing a Series.

    >>> s = pl.Series([12345, 12346])
    >>> pl.from_epoch(s, time_unit="d")
    shape: (2,)
    Series: '' [date]
    [
            2003-10-20
            2003-10-21
    ]
    """
    if isinstance(column, str):
        column = F.col(column)
    elif not isinstance(column, (pl.Series, pl.Expr)):
        column = pl.Series(column)  # Sequence input handled by Series constructor

    if time_unit == "d":
        return column.cast(Date)
    elif time_unit == "s":
        return (column.cast(Int64) * 1_000_000).cast(Datetime("us"))
    elif time_unit in DTYPE_TEMPORAL_UNITS:
        return column.cast(Datetime(time_unit))
    else:
        msg = f"`time_unit` must be one of {{'ns', 'us', 'ms', 's', 'd'}}, got {time_unit!r}"
        raise ValueError(msg)


@deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
def rolling_cov(
    a: str | Expr,
    b: str | Expr,
    *,
    window_size: int,
    min_samples: int | None = None,
    ddof: int = 1,
) -> Expr:
    """
    Compute the rolling covariance between two columns/ expressions.

    The window at a given row includes the row itself and the
    `window_size - 1` elements before it.

    .. versionchanged:: 1.21.0
        The `min_periods` parameter was renamed `min_samples`.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    window_size
        The length of the window.
    min_samples
        The number of values in the window that should be non-null before computing
        a result. If None, it will be set equal to window size.
    ddof
        Delta degrees of freedom. The divisor used in calculations
        is `N - ddof`, where `N` represents the number of elements.
    """
    if min_samples is None:
        min_samples = window_size
    if isinstance(a, str):
        a = F.col(a)
    if isinstance(b, str):
        b = F.col(b)
    return wrap_expr(
        plr.rolling_cov(a._pyexpr, b._pyexpr, window_size, min_samples, ddof)
    )


@deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
def rolling_corr(
    a: str | Expr,
    b: str | Expr,
    *,
    window_size: int,
    min_samples: int | None = None,
    ddof: int = 1,
) -> Expr:
    """
    Compute the rolling correlation between two columns/ expressions.

    The window at a given row includes the row itself and the
    `window_size - 1` elements before it.

    .. versionchanged:: 1.21.0
        The `min_periods` parameter was renamed `min_samples`.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    window_size
        The length of the window.
    min_samples
        The number of values in the window that should be non-null before computing
        a result. If None, it will be set equal to window size.
    ddof
        Delta degrees of freedom. The divisor used in calculations
        is `N - ddof`, where `N` represents the number of elements.
    """
    if min_samples is None:
        min_samples = window_size
    if isinstance(a, str):
        a = F.col(a)
    if isinstance(b, str):
        b = F.col(b)
    return wrap_expr(
        plr.rolling_corr(a._pyexpr, b._pyexpr, window_size, min_samples, ddof)
    )


@overload
def sql_expr(sql: str) -> Expr:  # type: ignore[overload-overlap]
    ...


@overload
def sql_expr(sql: Sequence[str]) -> list[Expr]: ...


def sql_expr(sql: str | Sequence[str]) -> Expr | list[Expr]:
    """
    Parse one or more SQL expressions to Polars expression(s).

    Parameters
    ----------
    sql
        One or more SQL expressions.

    Examples
    --------
    Parse a single SQL expression:

    >>> df = pl.DataFrame({"a": [2, 1]})
    >>> expr = pl.sql_expr("MAX(a)")
    >>> df.select(expr)
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 2   │
    └─────┘

    Parse multiple SQL expressions:

    >>> df.with_columns(
    ...     *pl.sql_expr(["POWER(a,a) AS a_a", "CAST(a AS TEXT) AS a_txt"]),
    ... )
    shape: (2, 3)
    ┌─────┬─────┬───────┐
    │ a   ┆ a_a ┆ a_txt │
    │ --- ┆ --- ┆ ---   │
    │ i64 ┆ i64 ┆ str   │
    ╞═════╪═════╪═══════╡
    │ 2   ┆ 4   ┆ 2     │
    │ 1   ┆ 1   ┆ 1     │
    └─────┴─────┴───────┘
    """
    if isinstance(sql, str):
        return wrap_expr(plr.sql_expr(sql))
    else:
        return [wrap_expr(plr.sql_expr(q)) for q in sql]


@unstable()
def row_index(name: str = "index") -> pl.Expr:
    """
    Generates a sequence of integers.

    The length of the returned sequence will match the context length, and the
    datatype will match the one returned by `get_index_dtype()`.

    .. versionadded:: 1.32.0

    If you would like to generate sequences with custom offsets / length /
    step size / datatypes, it is recommended to use `int_range` instead.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    Parameters
    ----------
    name
        Name of the returned column.

    Returns
    -------
    Expr
        Column of integers.

    See Also
    --------
    int_range : Generate a range of integers.

    Examples
    --------
    >>> df = pl.DataFrame({"x": ["A", "A", "B", "B", "B"]})
    >>> df.with_columns(pl.row_index(), pl.row_index("another_index"))
    shape: (5, 3)
    ┌─────┬───────┬───────────────┐
    │ x   ┆ index ┆ another_index │
    │ --- ┆ ---   ┆ ---           │
    │ str ┆ u32   ┆ u32           │
    ╞═════╪═══════╪═══════════════╡
    │ A   ┆ 0     ┆ 0             │
    │ A   ┆ 1     ┆ 1             │
    │ B   ┆ 2     ┆ 2             │
    │ B   ┆ 3     ┆ 3             │
    │ B   ┆ 4     ┆ 4             │
    └─────┴───────┴───────────────┘
    >>> df.group_by("x").agg(pl.row_index()).sort("x")
    shape: (2, 2)
    ┌─────┬───────────┐
    │ x   ┆ index     │
    │ --- ┆ ---       │
    │ str ┆ list[u32] │
    ╞═════╪═══════════╡
    │ A   ┆ [0, 1]    │
    │ B   ┆ [0, 1, 2] │
    └─────┴───────────┘
    >>> df.select(pl.row_index())
    shape: (5, 1)
    ┌───────┐
    │ index │
    │ ---   │
    │ u32   │
    ╞═══════╡
    │ 0     │
    │ 1     │
    │ 2     │
    │ 3     │
    │ 4     │
    └───────┘
    """
    # Notes
    # * Dispatching to `int_range` means that we cannot accept an offset
    #   parameter, as unlike `DataFrame.with_row_index()`, `int_range` will simply
    #   truncate instead of raising an error.
    return F.int_range(
        F.len(),
        dtype=get_index_type(),
    ).alias(name)
