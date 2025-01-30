from __future__ import annotations

from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from polars import Expr


class ExprNameNameSpace:
    """Namespace for expressions that operate on expression names."""

    _accessor = "name"

    def __init__(self, expr: Expr) -> None:
        self._from_pyexpr = expr._from_pyexpr
        self._pyexpr = expr._pyexpr

    def keep(self) -> Expr:
        """
        Keep the original root name of the expression.

        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.
        Consider using `.name.map` for advanced renaming.

        See Also
        --------
        Expr.alias
        map

        Examples
        --------
        Prevent errors due to potential duplicate column names.

        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2],
        ...         "b": [3, 4],
        ...     }
        ... )
        >>> df.select((pl.lit(10) / pl.all()).name.keep())
        shape: (2, 2)
        ┌──────┬──────────┐
        │ a    ┆ b        │
        │ ---  ┆ ---      │
        │ f64  ┆ f64      │
        ╞══════╪══════════╡
        │ 10.0 ┆ 3.333333 │
        │ 5.0  ┆ 2.5      │
        └──────┴──────────┘

        Undo an alias operation.

        >>> df.with_columns((pl.col("a") * 9).alias("c").name.keep())
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ b   │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 9   ┆ 3   │
        │ 18  ┆ 4   │
        └─────┴─────┘
        """
        return self._from_pyexpr(self._pyexpr.name_keep())

    def map(self, function: Callable[[str], str]) -> Expr:
        """
        Rename the output of an expression by mapping a function over the root name.

        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.


        Parameters
        ----------
        function
            Function that maps a root name to a new name.

        See Also
        --------
        keep
        prefix
        suffix

        Examples
        --------
        Remove a common suffix and convert to lower case.

        >>> df = pl.DataFrame(
        ...     {
        ...         "A_reverse": [3, 2, 1],
        ...         "B_reverse": ["z", "y", "x"],
        ...     }
        ... )
        >>> df.with_columns(
        ...     pl.all().reverse().name.map(lambda c: c.rstrip("_reverse").lower())
        ... )
        shape: (3, 4)
        ┌───────────┬───────────┬─────┬─────┐
        │ A_reverse ┆ B_reverse ┆ a   ┆ b   │
        │ ---       ┆ ---       ┆ --- ┆ --- │
        │ i64       ┆ str       ┆ i64 ┆ str │
        ╞═══════════╪═══════════╪═════╪═════╡
        │ 3         ┆ z         ┆ 1   ┆ x   │
        │ 2         ┆ y         ┆ 2   ┆ y   │
        │ 1         ┆ x         ┆ 3   ┆ z   │
        └───────────┴───────────┴─────┴─────┘
        """
        return self._from_pyexpr(self._pyexpr.name_map(function))

    def prefix(self, prefix: str) -> Expr:
        """
        Add a prefix to the root column name of the expression.

        Parameters
        ----------
        prefix
            Prefix to add to the root column name.


        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.
        Consider using `.name.map` for advanced renaming.

        See Also
        --------
        suffix

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.with_columns(pl.all().reverse().name.prefix("reverse_"))
        shape: (3, 4)
        ┌─────┬─────┬───────────┬───────────┐
        │ a   ┆ b   ┆ reverse_a ┆ reverse_b │
        │ --- ┆ --- ┆ ---       ┆ ---       │
        │ i64 ┆ str ┆ i64       ┆ str       │
        ╞═════╪═════╪═══════════╪═══════════╡
        │ 1   ┆ x   ┆ 3         ┆ z         │
        │ 2   ┆ y   ┆ 2         ┆ y         │
        │ 3   ┆ z   ┆ 1         ┆ x         │
        └─────┴─────┴───────────┴───────────┘
        """
        return self._from_pyexpr(self._pyexpr.name_prefix(prefix))

    def suffix(self, suffix: str) -> Expr:
        """
        Add a suffix to the root column name of the expression.

        Parameters
        ----------
        suffix
            Suffix to add to the root column name.

        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.
        Consider using `.name.map` for advanced renaming.

        See Also
        --------
        prefix

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "a": [1, 2, 3],
        ...         "b": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.with_columns(pl.all().reverse().name.suffix("_reverse"))
        shape: (3, 4)
        ┌─────┬─────┬───────────┬───────────┐
        │ a   ┆ b   ┆ a_reverse ┆ b_reverse │
        │ --- ┆ --- ┆ ---       ┆ ---       │
        │ i64 ┆ str ┆ i64       ┆ str       │
        ╞═════╪═════╪═══════════╪═══════════╡
        │ 1   ┆ x   ┆ 3         ┆ z         │
        │ 2   ┆ y   ┆ 2         ┆ y         │
        │ 3   ┆ z   ┆ 1         ┆ x         │
        └─────┴─────┴───────────┴───────────┘
        """
        return self._from_pyexpr(self._pyexpr.name_suffix(suffix))

    def to_lowercase(self) -> Expr:
        """
        Make the root column name lowercase.

        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.
        Consider using `.name.map` for advanced renaming.

        See Also
        --------
        prefix
        suffix
        to_uppercase

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "ColX": [1, 2, 3],
        ...         "ColY": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.with_columns(pl.all().name.to_lowercase())
        shape: (3, 4)
        ┌──────┬──────┬──────┬──────┐
        │ ColX ┆ ColY ┆ colx ┆ coly │
        │ ---  ┆ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ str  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╪══════╡
        │ 1    ┆ x    ┆ 1    ┆ x    │
        │ 2    ┆ y    ┆ 2    ┆ y    │
        │ 3    ┆ z    ┆ 3    ┆ z    │
        └──────┴──────┴──────┴──────┘
        """
        return self._from_pyexpr(self._pyexpr.name_to_lowercase())

    def to_uppercase(self) -> Expr:
        """
        Make the root column name uppercase.

        Notes
        -----
        This will undo any previous renaming operations on the expression.

        Due to implementation constraints, this method can only be called as the last
        expression in a chain. Only one name operation per expression will work.
        Consider using `.name.map` for advanced renaming.

        See Also
        --------
        prefix
        suffix
        to_lowercase

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "ColX": [1, 2, 3],
        ...         "ColY": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.with_columns(pl.all().name.to_uppercase())
        shape: (3, 4)
        ┌──────┬──────┬──────┬──────┐
        │ ColX ┆ ColY ┆ COLX ┆ COLY │
        │ ---  ┆ ---  ┆ ---  ┆ ---  │
        │ i64  ┆ str  ┆ i64  ┆ str  │
        ╞══════╪══════╪══════╪══════╡
        │ 1    ┆ x    ┆ 1    ┆ x    │
        │ 2    ┆ y    ┆ 2    ┆ y    │
        │ 3    ┆ z    ┆ 3    ┆ z    │
        └──────┴──────┴──────┴──────┘
        """
        return self._from_pyexpr(self._pyexpr.name_to_uppercase())

    def map_fields(self, function: Callable[[str], str]) -> Expr:
        """
        Rename fields of a struct by mapping a function over the field name(s).

        Notes
        -----
        This only takes effect for struct columns.

        Parameters
        ----------
        function
            Function that maps a field name to a new name.

        See Also
        --------
        prefix_fields
        suffix_fields

        Examples
        --------
        >>> df = pl.DataFrame({"x": {"a": 1, "b": 2}})
        >>> df.select(pl.col("x").name.map_fields(lambda x: x.upper())).schema
        Schema({'x': Struct({'A': Int64, 'B': Int64})})
        """
        return self._from_pyexpr(self._pyexpr.name_map_fields(function))

    def prefix_fields(self, prefix: str) -> Expr:
        """
        Add a prefix to all field names of a struct.

        Notes
        -----
        This only takes effect for struct columns.

        Parameters
        ----------
        prefix
            Prefix to add to the field name.

        See Also
        --------
        map_fields
        suffix_fields

        Examples
        --------
        >>> df = pl.DataFrame({"x": {"a": 1, "b": 2}})
        >>> df.select(pl.col("x").name.prefix_fields("prefix_")).schema
        Schema({'x': Struct({'prefix_a': Int64, 'prefix_b': Int64})})
        """
        return self._from_pyexpr(self._pyexpr.name_prefix_fields(prefix))

    def suffix_fields(self, suffix: str) -> Expr:
        """
        Add a suffix to all field names of a struct.

        Notes
        -----
        This only takes effect for struct columns.

        Parameters
        ----------
        suffix
            Suffix to add to the field name.

        See Also
        --------
        map_fields
        prefix_fields

        Examples
        --------
        >>> df = pl.DataFrame({"x": {"a": 1, "b": 2}})
        >>> df.select(pl.col("x").name.suffix_fields("_suffix")).schema
        Schema({'x': Struct({'a_suffix': Int64, 'b_suffix': Int64})})
        """
        return self._from_pyexpr(self._pyexpr.name_suffix_fields(suffix))
