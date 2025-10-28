from __future__ import annotations

import contextlib
import sys
from collections.abc import Iterable
from datetime import datetime, timedelta
from typing import TYPE_CHECKING

import polars._reexport as pl
from polars._utils.wrap import wrap_expr
from polars.datatypes import (
    Datetime,
    Duration,
    is_polars_dtype,
    parse_into_dtype,
)
from polars.datatypes.group import (
    DATETIME_DTYPES,
    DURATION_DTYPES,
    FLOAT_DTYPES,
    INTEGER_DTYPES,
)

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

if TYPE_CHECKING:
    from polars._typing import PolarsDataType, PythonDataType
    from polars.expr.expr import Expr

    if not sys.version_info >= (3, 11):
        from typing import Any

__all__ = ["col"]


def _create_col(
    name: (
        str
        | PolarsDataType
        | PythonDataType
        | Iterable[str]
        | Iterable[PolarsDataType | PythonDataType]
    ),
    *more_names: str | PolarsDataType | PythonDataType,
) -> Expr:
    """Create one or more column expressions representing column(s) in a DataFrame."""
    dtypes: list[PolarsDataType]
    if more_names:
        if isinstance(name, str):
            names_str = [name]
            names_str.extend(more_names)  # type: ignore[arg-type]
            return pl.Selector._by_name(names_str, strict=True).as_expr()
        elif is_polars_dtype(name):
            dtypes = [name]
            dtypes.extend(more_names)  # type: ignore[arg-type]
            return pl.Selector._by_dtype(dtypes).as_expr()  # type: ignore[arg-type]
        else:
            msg = (
                "invalid input for `col`"
                f"\n\nExpected `str` or `DataType`, got {type(name).__name__!r}."
            )
            raise TypeError(msg)

    if isinstance(name, str):
        return wrap_expr(plr.col(name))
    elif is_polars_dtype(name):
        dtypes = _polars_dtype_match(name)
        return pl.Selector._by_dtype(dtypes).as_expr()  # type: ignore[arg-type]
    elif isinstance(name, type):
        dtypes = _python_dtype_match(name)
        return pl.Selector._by_dtype(dtypes).as_expr()  # type: ignore[arg-type]
    elif isinstance(name, Iterable):
        names = list(name)
        if not names:
            return pl.Selector._by_name(names, strict=True).as_expr()  # type: ignore[arg-type]

        item = names[0]
        if isinstance(item, str):
            return pl.Selector._by_name(names, strict=True).as_expr()  # type: ignore[arg-type]
        elif is_polars_dtype(item):
            dtypes = []
            for nm in names:
                dtypes.extend(_polars_dtype_match(nm))  # type: ignore[arg-type]
            return pl.Selector._by_dtype(dtypes).as_expr()  # type: ignore[arg-type]
        elif isinstance(item, type):
            dtypes = []
            for nm in names:
                dtypes.extend(_python_dtype_match(nm))  # type: ignore[arg-type]
            return pl.Selector._by_dtype(dtypes).as_expr()  # type: ignore[arg-type]
        else:
            msg = (
                "invalid input for `col`"
                "\n\nExpected iterable of type `str` or `DataType`,"
                f" got iterable of type {type(item).__name__!r}."
            )
            raise TypeError(msg)
    else:
        msg = (
            "invalid input for `col`"
            f"\n\nExpected `str` or `DataType`, got {type(name).__name__!r}."
        )
        raise TypeError(msg)


def _python_dtype_match(tp: PythonDataType) -> list[PolarsDataType]:
    if tp is int:
        return list(INTEGER_DTYPES)
    elif tp is float:
        return list(FLOAT_DTYPES)
    elif tp is datetime:
        return list(DATETIME_DTYPES)
    elif tp is timedelta:
        return list(DURATION_DTYPES)
    return [parse_into_dtype(tp)]


def _polars_dtype_match(tp: PolarsDataType) -> list[PolarsDataType]:
    if Datetime.is_(tp):
        return list(DATETIME_DTYPES)
    elif Duration.is_(tp):
        return list(DURATION_DTYPES)
    return [tp]


class Col:
    """
    Create Polars column expressions.

    Notes
    -----
    An instance of this class is exported under the name `col`. It can be used as
    though it were a function by calling, for example, `pl.col("foo")`.
    See the :func:`__call__` method for further documentation.

    This helper class enables an alternative syntax for creating a column expression
    through attribute lookup. For example `col.foo` creates an expression equal to
    `col("foo")`. See the :func:`__getattr__` method for further documentation.

    The function call syntax is considered the idiomatic way of constructing a column
    expression. The alternative attribute syntax can be useful for quick prototyping as
    it can save some keystrokes, but has drawbacks in both expressiveness and
    readability.

    Examples
    --------
    >>> from polars import col
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": [1, 2],
    ...         "bar": [3, 4],
    ...     }
    ... )

    Create a new column expression using the standard syntax:

    >>> df.with_columns(baz=(col("foo") * col("bar")) / 2)
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ baz │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ f64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 1.5 │
    │ 2   ┆ 4   ┆ 4.0 │
    └─────┴─────┴─────┘

    Use attribute lookup to create a new column expression:

    >>> df.with_columns(baz=(col.foo + col.bar))
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ baz │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 4   │
    │ 2   ┆ 4   ┆ 6   │
    └─────┴─────┴─────┘
    """

    def __call__(
        self,
        name: (
            str
            | PolarsDataType
            | PythonDataType
            | Iterable[str]
            | Iterable[PolarsDataType | PythonDataType]
        ),
        *more_names: str | PolarsDataType | PythonDataType,
    ) -> Expr:
        """
        Create one or more expressions representing columns in a DataFrame.

        Parameters
        ----------
        name
            The name or datatype of the column(s) to represent.
            Accepts regular expression input; regular expressions
            should start with `^` and end with `$`.
        *more_names
            Additional names or datatypes of columns to represent,
            specified as positional arguments.

        See Also
        --------
        first
        last
        nth

        Examples
        --------
        Pass a single column name to represent that column.

        >>> df = pl.DataFrame(
        ...     {
        ...         "ham": [1, 2],
        ...         "hamburger": [11, 22],
        ...         "foo": [2, 1],
        ...         "bar": ["a", "b"],
        ...     }
        ... )
        >>> df.select(pl.col("foo"))
        shape: (2, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 2   │
        │ 1   │
        └─────┘

        Use dot syntax to save keystrokes for quick prototyping.

        >>> from polars import col as c
        >>> df.select(c.foo + c.ham)
        shape: (2, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 3   │
        │ 3   │
        └─────┘

        Use the wildcard `*` to represent all columns.

        >>> df.select(pl.col("*"))
        shape: (2, 4)
        ┌─────┬───────────┬─────┬─────┐
        │ ham ┆ hamburger ┆ foo ┆ bar │
        │ --- ┆ ---       ┆ --- ┆ --- │
        │ i64 ┆ i64       ┆ i64 ┆ str │
        ╞═════╪═══════════╪═════╪═════╡
        │ 1   ┆ 11        ┆ 2   ┆ a   │
        │ 2   ┆ 22        ┆ 1   ┆ b   │
        └─────┴───────────┴─────┴─────┘
        >>> df.select(pl.col("*").exclude("ham"))
        shape: (2, 3)
        ┌───────────┬─────┬─────┐
        │ hamburger ┆ foo ┆ bar │
        │ ---       ┆ --- ┆ --- │
        │ i64       ┆ i64 ┆ str │
        ╞═══════════╪═════╪═════╡
        │ 11        ┆ 2   ┆ a   │
        │ 22        ┆ 1   ┆ b   │
        └───────────┴─────┴─────┘

        Regular expression input is supported.

        >>> df.select(pl.col("^ham.*$"))
        shape: (2, 2)
        ┌─────┬───────────┐
        │ ham ┆ hamburger │
        │ --- ┆ ---       │
        │ i64 ┆ i64       │
        ╞═════╪═══════════╡
        │ 1   ┆ 11        │
        │ 2   ┆ 22        │
        └─────┴───────────┘

        Multiple columns can be represented by passing a list of names.

        >>> df.select(pl.col(["hamburger", "foo"]))
        shape: (2, 2)
        ┌───────────┬─────┐
        │ hamburger ┆ foo │
        │ ---       ┆ --- │
        │ i64       ┆ i64 │
        ╞═══════════╪═════╡
        │ 11        ┆ 2   │
        │ 22        ┆ 1   │
        └───────────┴─────┘

        Or use positional arguments to represent multiple columns in the same way.

        >>> df.select(pl.col("hamburger", "foo"))
        shape: (2, 2)
        ┌───────────┬─────┐
        │ hamburger ┆ foo │
        │ ---       ┆ --- │
        │ i64       ┆ i64 │
        ╞═══════════╪═════╡
        │ 11        ┆ 2   │
        │ 22        ┆ 1   │
        └───────────┴─────┘

        Easily select all columns that match a certain data type by passing that
        datatype.

        >>> df.select(pl.col(pl.String))
        shape: (2, 1)
        ┌─────┐
        │ bar │
        │ --- │
        │ str │
        ╞═════╡
        │ a   │
        │ b   │
        └─────┘
        >>> df.select(pl.col(pl.Int64, pl.Float64))
        shape: (2, 3)
        ┌─────┬───────────┬─────┐
        │ ham ┆ hamburger ┆ foo │
        │ --- ┆ ---       ┆ --- │
        │ i64 ┆ i64       ┆ i64 │
        ╞═════╪═══════════╪═════╡
        │ 1   ┆ 11        ┆ 2   │
        │ 2   ┆ 22        ┆ 1   │
        └─────┴───────────┴─────┘
        """
        return _create_col(name, *more_names)

    def __getattr__(self, name: str) -> Expr:
        """
        Create a column expression using attribute syntax.

        Note that this syntax does not support passing data
        types or multiple column names.

        Parameters
        ----------
        name
            The name of the column to represent.

        Examples
        --------
        >>> from polars import col as c
        >>> df = pl.DataFrame(
        ...     {
        ...         "foo": [1, 2],
        ...         "bar": [3, 4],
        ...     }
        ... )
        >>> df.select(c.foo + c.bar)
        shape: (2, 1)
        ┌─────┐
        │ foo │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 4   │
        │ 6   │
        └─────┘
        """
        # For autocomplete to work with IPython
        if name.startswith("__wrapped__"):
            return getattr(type(self), name)

        return _create_col(name)

    if not sys.version_info >= (3, 11):

        def __getstate__(self) -> Any:
            return self.__dict__

        def __setstate__(self, state: Any) -> None:
            self.__dict__ = state


col: Col = Col()
