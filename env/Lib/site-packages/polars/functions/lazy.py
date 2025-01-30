from __future__ import annotations

import contextlib
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable, overload

import polars._reexport as pl
import polars.functions as F
from polars._utils.async_ import _AioDataFrameResult, _GeventDataFrameResult
from polars._utils.deprecation import deprecate_function, issue_deprecation_warning
from polars._utils.parse import (
    parse_into_expression,
    parse_into_list_of_expressions,
)
from polars._utils.unstable import issue_unstable_warning, unstable
from polars._utils.various import extend_bool
from polars._utils.wrap import wrap_df, wrap_expr
from polars.datatypes import DTYPE_TEMPORAL_UNITS, Date, Datetime, Int64

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr

if TYPE_CHECKING:
    from collections.abc import Awaitable, Collection, Iterable
    from typing import Literal

    from polars import DataFrame, Expr, LazyFrame, Series
    from polars._typing import (
        CorrelationMethod,
        EpochTimeUnit,
        IntoExpr,
        PolarsDataType,
        RollingInterpolationMethod,
    )


def field(name: str | list[str]) -> Expr:
    """
    Select a field in the current `struct.with_fields` scope.

    name
        Name of the field(s) to select.
    """
    if isinstance(name, str):
        name = [name]
    return wrap_expr(plr.field(name))


def element() -> Expr:
    """
    Alias for an element being evaluated in an `eval` expression.

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
    >>> df.select(pl.cum_count("a"))
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ u32 │
    ╞═════╡
    │ 1   │
    │ 2   │
    │ 2   │
    └─────┘
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
        return wrap_expr(plr.first())

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
        return wrap_expr(plr.last())

    return F.col(*columns).last()


def nth(*indices: int | Sequence[int]) -> Expr:
    """
    Get the nth column(s) of the context.

    Parameters
    ----------
    indices
        One or more indices representing the columns to retrieve.

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
    if len(indices) == 1 and isinstance(indices[0], Sequence):
        indices = indices[0]  # type: ignore[assignment]

    return wrap_expr(plr.index_cols(indices))


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


def corr(
    a: IntoExpr,
    b: IntoExpr,
    *,
    method: CorrelationMethod = "pearson",
    ddof: int | None = None,
    propagate_nans: bool = False,
) -> Expr:
    """
    Compute the Pearson's or Spearman rank correlation correlation between two columns.

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

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.corr("a", "b", method="spearman"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 0.5 │
    └─────┘
    """
    if ddof is not None:
        issue_deprecation_warning(
            "The `ddof` parameter has no effect. Do not use it.",
            version="1.17.0",
        )

    a = parse_into_expression(a)
    b = parse_into_expression(b)

    if method == "pearson":
        return wrap_expr(plr.pearson_corr(a, b))
    elif method == "spearman":
        return wrap_expr(plr.spearman_rank_corr(a, b, propagate_nans))
    else:
        msg = f"method must be one of {{'pearson', 'spearman'}}, got {method!r}"
        raise ValueError(msg)


def cov(a: IntoExpr, b: IntoExpr, ddof: int = 1) -> Expr:
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

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     },
    ... )
    >>> df.select(pl.cov("a", "b"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 3.0 │
    └─────┘
    """
    a = parse_into_expression(a)
    b = parse_into_expression(b)
    return wrap_expr(plr.cov(a, b, ddof))


def map_batches(
    exprs: Sequence[str] | Sequence[Expr],
    function: Callable[[Sequence[Series]], Series],
    return_dtype: PolarsDataType | None = None,
) -> Expr:
    """
    Map a custom function over multiple columns/expressions.

    Produces a single Series result.

    Parameters
    ----------
    exprs
        Expression(s) representing the input Series to the function.
    function
        Function to apply over the input.
    return_dtype
        dtype of the output Series.

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
    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.map_mul(
            exprs, function, return_dtype, map_groups=False, returns_scalar=False
        )
    )


def map_groups(
    exprs: Sequence[str | Expr],
    function: Callable[[Sequence[Series]], Series | Any],
    return_dtype: PolarsDataType | None = None,
    *,
    returns_scalar: bool = True,
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
        dtype of the output Series.
    returns_scalar
        If the function returns a single scalar as output.

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
    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(
        plr.map_mul(
            exprs,
            function,
            return_dtype,
            map_groups=True,
            returns_scalar=returns_scalar,
        )
    )


def fold(
    acc: IntoExpr,
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
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
    │ i64 │
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
    acc = parse_into_expression(acc, str_as_lit=True)
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(plr.fold(acc, function, exprs))


def reduce(
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
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
    # in case of col("*")
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(plr.reduce(function, exprs))


def cum_fold(
    acc: IntoExpr,
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
    *,
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
    acc = parse_into_expression(acc, str_as_lit=True)
    if isinstance(exprs, pl.Expr):
        exprs = [exprs]

    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(plr.cum_fold(acc, function, exprs, include_init).alias("cum_fold"))


def cum_reduce(
    function: Callable[[Series, Series], Series],
    exprs: Sequence[Expr | str] | Expr,
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

    exprs = parse_into_list_of_expressions(exprs)
    return wrap_expr(plr.cum_reduce(function, exprs).alias("cum_reduce"))


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
    return wrap_expr(plr.arctan2(y._pyexpr, x._pyexpr))


@deprecate_function("Use `arctan2` followed by `.degrees()` instead.", version="1.0.0")
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
    """Syntactic sugar for `pl.col("foo").agg_groups()`."""
    return F.col(column).agg_groups()


def quantile(
    column: str,
    quantile: float | Expr,
    interpolation: RollingInterpolationMethod = "nearest",
) -> Expr:
    """
    Syntactic sugar for `pl.col("foo").quantile(..)`.

    Parameters
    ----------
    column
        Column name.
    quantile
        Quantile between 0.0 and 1.0.
    interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear'}
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


def collect_all(
    lazy_frames: Iterable[LazyFrame],
    *,
    type_coercion: bool = True,
    _type_check: bool = True,
    predicate_pushdown: bool = True,
    projection_pushdown: bool = True,
    simplify_expression: bool = True,
    no_optimization: bool = False,
    slice_pushdown: bool = True,
    comm_subplan_elim: bool = True,
    comm_subexpr_elim: bool = True,
    cluster_with_columns: bool = True,
    collapse_joins: bool = True,
    streaming: bool = False,
    _check_order: bool = True,
) -> list[DataFrame]:
    """
    Collect multiple LazyFrames at the same time.

    This runs all the computation graphs in parallel on the Polars threadpool.

    Parameters
    ----------
    lazy_frames
        A list of LazyFrames to collect.
    type_coercion
        Do type coercion optimization.
    predicate_pushdown
        Do predicate pushdown optimization.
    projection_pushdown
        Do projection pushdown optimization.
    simplify_expression
        Run simplify expressions optimization.
    no_optimization
        Turn off optimizations.
    slice_pushdown
        Slice pushdown optimization.
    comm_subplan_elim
        Will try to cache branching subplans that occur on self-joins or unions.
    comm_subexpr_elim
        Common subexpressions will be cached and reused.
    cluster_with_columns
        Combine sequential independent calls to with_columns
    collapse_joins
        Collapse a join and filters into a faster join
    streaming
        Process the query in batches to handle larger-than-memory data.
        If set to `False` (default), the entire query is processed in a single
        batch.

        .. warning::
            Streaming mode is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. note::
            Use :func:`explain` to see if Polars can process the query in streaming
            mode.

    Returns
    -------
    list of DataFrames
        The collected DataFrames, returned in the same order as the input LazyFrames.
    """
    if no_optimization:
        predicate_pushdown = False
        projection_pushdown = False
        slice_pushdown = False
        comm_subplan_elim = False
        comm_subexpr_elim = False
        cluster_with_columns = False
        collapse_joins = False

    if streaming:
        issue_unstable_warning("Streaming mode is considered unstable.")
        comm_subplan_elim = False

    prepared = []

    for lf in lazy_frames:
        type_check = _type_check
        ldf = lf._ldf.optimization_toggle(
            type_coercion,
            type_check,
            predicate_pushdown,
            projection_pushdown,
            simplify_expression,
            slice_pushdown,
            comm_subplan_elim,
            comm_subexpr_elim,
            cluster_with_columns,
            collapse_joins,
            streaming,
            _eager=False,
            _check_order=_check_order,
            new_streaming=False,
        )
        prepared.append(ldf)

    out = plr.collect_all(prepared)

    # wrap the pydataframes into dataframe
    result = [wrap_df(pydf) for pydf in out]

    return result


@overload
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: Literal[True],
    type_coercion: bool = True,
    _type_check: bool = True,
    predicate_pushdown: bool = True,
    projection_pushdown: bool = True,
    simplify_expression: bool = True,
    no_optimization: bool = True,
    slice_pushdown: bool = True,
    comm_subplan_elim: bool = True,
    comm_subexpr_elim: bool = True,
    cluster_with_columns: bool = True,
    collapse_joins: bool = True,
    streaming: bool = True,
) -> _GeventDataFrameResult[list[DataFrame]]: ...


@overload
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: Literal[False] = False,
    type_coercion: bool = True,
    _type_check: bool = True,
    predicate_pushdown: bool = True,
    projection_pushdown: bool = True,
    simplify_expression: bool = True,
    no_optimization: bool = False,
    slice_pushdown: bool = True,
    comm_subplan_elim: bool = True,
    comm_subexpr_elim: bool = True,
    cluster_with_columns: bool = True,
    collapse_joins: bool = True,
    streaming: bool = False,
) -> Awaitable[list[DataFrame]]: ...


@unstable()
def collect_all_async(
    lazy_frames: Iterable[LazyFrame],
    *,
    gevent: bool = False,
    type_coercion: bool = True,
    _type_check: bool = True,
    predicate_pushdown: bool = True,
    projection_pushdown: bool = True,
    simplify_expression: bool = True,
    no_optimization: bool = False,
    slice_pushdown: bool = True,
    comm_subplan_elim: bool = True,
    comm_subexpr_elim: bool = True,
    cluster_with_columns: bool = True,
    collapse_joins: bool = True,
    streaming: bool = False,
    _check_order: bool = True,
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
    type_coercion
        Do type coercion optimization.
    predicate_pushdown
        Do predicate pushdown optimization.
    projection_pushdown
        Do projection pushdown optimization.
    simplify_expression
        Run simplify expressions optimization.
    no_optimization
        Turn off (certain) optimizations.
    slice_pushdown
        Slice pushdown optimization.
    comm_subplan_elim
        Will try to cache branching subplans that occur on self-joins or unions.
    comm_subexpr_elim
        Common subexpressions will be cached and reused.
    cluster_with_columns
        Combine sequential independent calls to with_columns
    collapse_joins
        Collapse a join and filters into a faster join
    streaming
        Process the query in batches to handle larger-than-memory data.
        If set to `False` (default), the entire query is processed in a single
        batch.

        .. warning::
            Streaming mode is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. note::
            Use :func:`explain` to see if Polars can process the query in streaming
            mode.

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
    if no_optimization:
        predicate_pushdown = False
        projection_pushdown = False
        slice_pushdown = False
        comm_subplan_elim = False
        comm_subexpr_elim = False
        cluster_with_columns = False
        collapse_joins = False

    if streaming:
        issue_unstable_warning("Streaming mode is considered unstable.")
        comm_subplan_elim = False

    prepared = []

    for lf in lazy_frames:
        type_check = _type_check
        ldf = lf._ldf.optimization_toggle(
            type_coercion,
            type_check,
            predicate_pushdown,
            projection_pushdown,
            simplify_expression,
            slice_pushdown,
            comm_subplan_elim,
            comm_subexpr_elim,
            cluster_with_columns,
            collapse_joins,
            streaming,
            _eager=False,
            _check_order=_check_order,
            new_streaming=False,
        )
        prepared.append(ldf)

    result: (
        _GeventDataFrameResult[list[DataFrame]] | _AioDataFrameResult[list[DataFrame]]
    ) = _GeventDataFrameResult() if gevent else _AioDataFrameResult()
    plr.collect_all_with_callback(prepared, result._callback_all)
    return result


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


@overload
def arg_where(condition: Expr | Series, *, eager: bool) -> Expr | Series: ...


def arg_where(condition: Expr | Series, *, eager: bool = False) -> Expr | Series:
    """
    Return indices where `condition` evaluates `True`.

    Parameters
    ----------
    condition
        Boolean expression to evaluate
    eager
        Evaluate immediately and return a `Series`. If set to `False` (default),
        return an expression instead.

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
                "expected 'Series' in 'arg_where' if 'eager=True', got"
                f" {type(condition).__name__!r}"
            )
            raise ValueError(msg)
        return condition.to_frame().select(arg_where(F.col(condition.name))).to_series()
    else:
        condition = parse_into_expression(condition)
        return wrap_expr(plr.arg_where(condition))


def coalesce(exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr) -> Expr:
    """
    Folds the columns from left to right, keeping the first non-null value.

    Parameters
    ----------
    exprs
        Columns to coalesce. Accepts expression input. Strings are parsed as column
        names, other non-expression inputs are parsed as literals.
    *more_exprs
        Additional columns to coalesce, specified as positional arguments.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, None, None, None],
    ...         "b": [1, 2, None, None],
    ...         "c": [5, None, 3, None],
    ...     }
    ... )
    >>> df.with_columns(pl.coalesce(["a", "b", "c", 10]).alias("d"))
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
    """
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


@unstable()
def rolling_cov(
    a: str | Expr,
    b: str | Expr,
    *,
    window_size: int,
    min_periods: int | None = None,
    ddof: int = 1,
) -> Expr:
    """
    Compute the rolling covariance between two columns/ expressions.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    The window at a given row includes the row itself and the
    `window_size - 1` elements before it.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    window_size
        The length of the window.
    min_periods
        The number of values in the window that should be non-null before computing
        a result. If None, it will be set equal to window size.
    ddof
        Delta degrees of freedom. The divisor used in calculations
        is `N - ddof`, where `N` represents the number of elements.
    """
    if min_periods is None:
        min_periods = window_size
    if isinstance(a, str):
        a = F.col(a)
    if isinstance(b, str):
        b = F.col(b)
    return wrap_expr(
        plr.rolling_cov(a._pyexpr, b._pyexpr, window_size, min_periods, ddof)
    )


@unstable()
def rolling_corr(
    a: str | Expr,
    b: str | Expr,
    *,
    window_size: int,
    min_periods: int | None = None,
    ddof: int = 1,
) -> Expr:
    """
    Compute the rolling correlation between two columns/ expressions.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    The window at a given row includes the row itself and the
    `window_size - 1` elements before it.

    Parameters
    ----------
    a
        Column name or Expression.
    b
        Column name or Expression.
    window_size
        The length of the window.
    min_periods
        The number of values in the window that should be non-null before computing
        a result. If None, it will be set equal to window size.
    ddof
        Delta degrees of freedom. The divisor used in calculations
        is `N - ddof`, where `N` represents the number of elements.
    """
    if min_periods is None:
        min_periods = window_size
    if isinstance(a, str):
        a = F.col(a)
    if isinstance(b, str):
        b = F.col(b)
    return wrap_expr(
        plr.rolling_corr(a._pyexpr, b._pyexpr, window_size, min_periods, ddof)
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
