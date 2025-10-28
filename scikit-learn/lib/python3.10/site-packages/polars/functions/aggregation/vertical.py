from __future__ import annotations

from typing import TYPE_CHECKING

import polars.functions as F

if TYPE_CHECKING:
    from polars import Expr


def all(*names: str, ignore_nulls: bool = True) -> Expr:
    """
    Either return an expression representing all columns, or evaluate a bitwise AND operation.

    If no arguments are passed, this function is syntactic sugar for `col("*")`.
    Otherwise, this function is syntactic sugar for `col(names).all()`.

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.
    ignore_nulls
        * If set to `True` (default), null values are ignored. If there
          are no non-null values, the output is `True`.
        * If set to `False`, `Kleene logic`_ is used to deal with nulls:
          if the column contains any null values and no `False` values,
          the output is null.

        .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

    See Also
    --------
    all_horizontal

    Examples
    --------
    Selecting all columns.

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [True, False, True],
    ...         "b": [False, False, False],
    ...     }
    ... )
    >>> df.select(pl.all().sum())
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ u32 ┆ u32 │
    ╞═════╪═════╡
    │ 2   ┆ 0   │
    └─────┴─────┘

    Evaluate bitwise AND for a column.

    >>> df.select(pl.all("a"))
    shape: (1, 1)
    ┌───────┐
    │ a     │
    │ ---   │
    │ bool  │
    ╞═══════╡
    │ false │
    └───────┘
    """  # noqa: W505
    if not names:
        return F.col("*")

    return F.col(*names).all(ignore_nulls=ignore_nulls)


def any(*names: str, ignore_nulls: bool = True) -> Expr | bool | None:
    """
    Evaluate a bitwise OR operation.

    Syntactic sugar for `col(names).any()`.

    See Also
    --------
    any_horizontal

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.
    ignore_nulls
        * If set to `True` (default), null values are ignored. If there
          are no non-null values, the output is `False`.
        * If set to `False`, `Kleene logic`_ is used to deal with nulls:
          if the column contains any null values and no `True` values,
          the output is null.

        .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [True, False, True],
    ...         "b": [False, False, False],
    ...     }
    ... )
    >>> df.select(pl.any("a"))
    shape: (1, 1)
    ┌──────┐
    │ a    │
    │ ---  │
    │ bool │
    ╞══════╡
    │ true │
    └──────┘
    """
    return F.col(*names).any(ignore_nulls=ignore_nulls)


def max(*names: str) -> Expr:
    """
    Get the maximum value.

    Syntactic sugar for `col(names).max()`.

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.

    See Also
    --------
    max_horizontal

    Examples
    --------
    Get the maximum value of a column.

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.max("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 8   │
    └─────┘

    Get the maximum value of multiple columns.

    >>> df.select(pl.max("^a|b$"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 8   ┆ 5   │
    └─────┴─────┘
    >>> df.select(pl.max("a", "b"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 8   ┆ 5   │
    └─────┴─────┘
    """
    return F.col(*names).max()


def min(*names: str) -> Expr:
    """
    Get the minimum value.

    Syntactic sugar for `col(names).min()`.

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.

    See Also
    --------
    min_horizontal

    Examples
    --------
    Get the minimum value of a column.

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 8, 3],
    ...         "b": [4, 5, 2],
    ...         "c": ["foo", "bar", "foo"],
    ...     }
    ... )
    >>> df.select(pl.min("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    └─────┘

    Get the minimum value of multiple columns.

    >>> df.select(pl.min("^a|b$"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 2   │
    └─────┴─────┘
    >>> df.select(pl.min("a", "b"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 2   │
    └─────┴─────┘
    """
    return F.col(*names).min()


def sum(*names: str) -> Expr:
    """
    Sum all values.

    Syntactic sugar for `col(name).sum()`.

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.

    Notes
    -----
    If there are no non-null values, then the output is `0`.
    If you would prefer empty sums to return `None`, you can
    use `pl.when(pl.col(name).count()>0).then(pl.sum(name))` instead
    of `pl.sum(name)`.

    See Also
    --------
    sum_horizontal

    Examples
    --------
    Sum a column.

    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2],
    ...         "b": [3, 4],
    ...         "c": [5, 6],
    ...     }
    ... )
    >>> df.select(pl.sum("a"))
    shape: (1, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 3   │
    └─────┘

    Sum multiple columns.

    >>> df.select(pl.sum("a", "c"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ a   ┆ c   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 3   ┆ 11  │
    └─────┴─────┘
    >>> df.select(pl.sum("^.*[bc]$"))
    shape: (1, 2)
    ┌─────┬─────┐
    │ b   ┆ c   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 7   ┆ 11  │
    └─────┴─────┘
    """
    return F.col(*names).sum()


def cum_sum(*names: str) -> Expr:
    """
    Cumulatively sum all values.

    Syntactic sugar for `col(names).cum_sum()`.

    Parameters
    ----------
    *names
        Name(s) of the columns to use in the aggregation.

    See Also
    --------
    cumsum_horizontal

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "a": [1, 2, 3],
    ...         "b": [4, 5, 6],
    ...     }
    ... )
    >>> df.select(pl.cum_sum("a"))
    shape: (3, 1)
    ┌─────┐
    │ a   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 3   │
    │ 6   │
    └─────┘
    """
    return F.col(*names).cum_sum()
