from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, Any

import polars._reexport as pl
from polars._utils.parse import parse_predicates_constraints_into_expression

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr

if TYPE_CHECKING:
    from collections.abc import Iterable

    from polars._typing import IntoExprColumn


def when(
    *predicates: IntoExprColumn | Iterable[IntoExprColumn] | bool,
    **constraints: Any,
) -> pl.When:
    """
    Start a `when-then-otherwise` expression.

    Expression similar to an `if-else` statement in Python. Always initiated by a
    `pl.when(<condition>).then(<value if condition>)`., and optionally followed by
    chaining one or more `.when(<condition>).then(<value>)` statements.

    Chained when-then operations should be read as Python `if, elif, ... elif` blocks,
    not as `if, if, ... if`, i.e. the first condition that evaluates to True will be
    picked.

    If none of the conditions are `True`, an optional `.otherwise(<value if all
    statements are false>)` can be appended at the end. If not appended, and none
    of the conditions are `True`, `None` will be returned.

    Parameters
    ----------
    predicates
        Condition(s) that must be met in order to apply the subsequent statement.
        Accepts one or more boolean expressions, which are implicitly combined with
        `&`. String input is parsed as a column name.
    constraints
        Apply conditions as `col_name = value` keyword arguments that are treated as
        equality matches, such as `x = 123`. As with the predicates parameter, multiple
        conditions are implicitly combined using `&`.

    Warnings
    --------
    Polars computes all expressions passed to `when-then-otherwise` in parallel and
    filters afterwards. This means each expression must be valid on its own, regardless
    of the conditions in the `when-then-otherwise` chain.

    Examples
    --------
    Below we add a column with the value 1, where column "foo" > 2 and the value -1
    where it isn't.

    >>> df = pl.DataFrame({"foo": [1, 3, 4], "bar": [3, 4, 0]})
    >>> df.with_columns(pl.when(pl.col("foo") > 2).then(1).otherwise(-1).alias("val"))
    shape: (3, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ val │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i32 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ -1  │
    │ 3   ┆ 4   ┆ 1   │
    │ 4   ┆ 0   ┆ 1   │
    └─────┴─────┴─────┘

    Or with multiple when-then operations chained:

    >>> df.with_columns(
    ...     pl.when(pl.col("foo") > 2)
    ...     .then(1)
    ...     .when(pl.col("bar") > 2)
    ...     .then(4)
    ...     .otherwise(-1)
    ...     .alias("val")
    ... )
    shape: (3, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ val │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i32 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 4   │
    │ 3   ┆ 4   ┆ 1   │
    │ 4   ┆ 0   ┆ 1   │
    └─────┴─────┴─────┘

    Note how in the example above for the second row in the DataFrame,
    where `foo=3` and `bar=4`, the first `when` evaluates to `True`, and therefore
    the second `when`, which is also `True`, is not evaluated.

    The `otherwise` at the end is optional. If left out, any rows where none
    of the `when` expressions evaluate to True, are set to `null`:

    >>> df.with_columns(pl.when(pl.col("foo") > 2).then(1).alias("val"))
    shape: (3, 3)
    ┌─────┬─────┬──────┐
    │ foo ┆ bar ┆ val  │
    │ --- ┆ --- ┆ ---  │
    │ i64 ┆ i64 ┆ i32  │
    ╞═════╪═════╪══════╡
    │ 1   ┆ 3   ┆ null │
    │ 3   ┆ 4   ┆ 1    │
    │ 4   ┆ 0   ┆ 1    │
    └─────┴─────┴──────┘

    Pass multiple predicates, each of which must be met:

    >>> df.with_columns(
    ...     val=pl.when(
    ...         pl.col("bar") > 0,
    ...         pl.col("foo") % 2 != 0,
    ...     )
    ...     .then(99)
    ...     .otherwise(-1)
    ... )
    shape: (3, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ val │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i32 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 99  │
    │ 3   ┆ 4   ┆ 99  │
    │ 4   ┆ 0   ┆ -1  │
    └─────┴─────┴─────┘

    Pass conditions as keyword arguments:

    >>> df.with_columns(val=pl.when(foo=4, bar=0).then(99).otherwise(-1))
    shape: (3, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ val │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i32 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ -1  │
    │ 3   ┆ 4   ┆ -1  │
    │ 4   ┆ 0   ┆ 99  │
    └─────┴─────┴─────┘
    """
    condition = parse_predicates_constraints_into_expression(*predicates, **constraints)
    return pl.When(plr.when(condition))
