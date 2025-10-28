from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from polars._utils.wrap import wrap_expr

if TYPE_CHECKING:
    from polars import Expr


class ExprNameNameSpace:
    """Namespace for expressions that operate on expression names."""

    _accessor = "name"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def keep(self) -> Expr:
        """
        Keep the original root name of the expression.

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
        return wrap_expr(self._pyexpr.name_keep())

    def map(self, function: Callable[[str], str]) -> Expr:
        """
        Rename the output of an expression by mapping a function over the root name.

        Parameters
        ----------
        function
            Function that maps a root name to a new name.

        See Also
        --------
        keep
        prefix
        suffix
        replace

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
        ...     pl.all()
        ...     .reverse()
        ...     .name.map(lambda c: c.removesuffix("_reverse").lower())
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
        return wrap_expr(self._pyexpr.name_map(function))

    def prefix(self, prefix: str) -> Expr:
        """
        Add a prefix to the root column name of the expression.

        Parameters
        ----------
        prefix
            Prefix to add to the root column name.

        See Also
        --------
        suffix
        map

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
        return wrap_expr(self._pyexpr.name_prefix(prefix))

    def suffix(self, suffix: str) -> Expr:
        """
        Add a suffix to the root column name of the expression.

        Parameters
        ----------
        suffix
            Suffix to add to the root column name.

        See Also
        --------
        prefix
        map

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
        return wrap_expr(self._pyexpr.name_suffix(suffix))

    def to_lowercase(self) -> Expr:
        """
        Make the root column name lowercase.

        See Also
        --------
        prefix
        suffix
        to_uppercase
        map

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
        return wrap_expr(self._pyexpr.name_to_lowercase())

    def to_uppercase(self) -> Expr:
        """
        Make the root column name uppercase.

        See Also
        --------
        prefix
        suffix
        to_lowercase
        map

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
        return wrap_expr(self._pyexpr.name_to_uppercase())

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
        return wrap_expr(self._pyexpr.name_map_fields(function))

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
        return wrap_expr(self._pyexpr.name_prefix_fields(prefix))

    def replace(self, pattern: str, value: str, *, literal: bool = False) -> Expr:
        r"""
        Replace matching regex/literal substring in the name with a new value.

        Parameters
        ----------
        pattern
            A valid regular expression pattern, compatible with the `regex crate
            <https://docs.rs/regex/latest/regex/>`_.
        value
            String that will replace the matched substring.
        literal
            Treat `pattern` as a literal string, not a regex.

        Notes
        -----
        * To modify regular expression behaviour (such as case-sensitivity) with flags,
          use the inline `(?iLmsuxU)` syntax. See the regex crate's section on
          `grouping and flags <https://docs.rs/regex/latest/regex/#grouping-and-flags>`_
          for additional information about the use of inline expression modifiers.

        * The dollar sign (`$`) is a special character related to capture groups; if you
          want to replace some target pattern with characters that include a literal `$`
          you should escape it by doubling it up as `$$`, or set `literal=True` if you
          do not need a full regular expression pattern match. Otherwise, you will be
          referencing a (potentially non-existent) capture group.

        See Also
        --------
        Expr.str.replace

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "n_foo": [1, 2, 3],
        ...         "n_bar": ["x", "y", "z"],
        ...     }
        ... )
        >>> df.select(pl.all().name.replace(r"^n_", "col_"))
        shape: (3, 2)
        ┌─────────┬─────────┐
        │ col_foo ┆ col_bar │
        │ ---     ┆ ---     │
        │ i64     ┆ str     │
        ╞═════════╪═════════╡
        │ 1       ┆ x       │
        │ 2       ┆ y       │
        │ 3       ┆ z       │
        └─────────┴─────────┘
        >>> df.select(pl.all().name.replace(r"(a|e|i|o|u)", "@")).schema
        Schema({'n_f@@': Int64, 'n_b@r': String})

        Apply case-insensitive string replacement using the `(?i)` flag.

        >>> pl.DataFrame({"Foo": [1], "faz": [2]}).select(
        ...     pl.all().name.replace(r"(?i)^f", "b")
        ... )
        shape: (1, 2)
        ┌─────┬─────┐
        │ boo ┆ baz │
        │ --- ┆ --- │
        │ i64 ┆ i64 │
        ╞═════╪═════╡
        │ 1   ┆ 2   │
        └─────┴─────┘

        Capture groups are supported. Use `$1` or `${1}` in the `value` string to refer
        to the first capture group in the pattern, `$2` or `${2}` to refer to the
        second capture group, and so on. You can also use named capture groups.

        >>> df = pl.DataFrame({"x_1": [1], "x_2": [2], "group_id": ["xyz"]})
        >>> df.select(pl.all().name.replace(r"_(\d+)$", ":$1"))
        shape: (1, 3)
        ┌─────┬─────┬──────────┐
        │ x:1 ┆ x:2 ┆ group_id │
        │ --- ┆ --- ┆ ---      │
        │ i64 ┆ i64 ┆ str      │
        ╞═════╪═════╪══════════╡
        │ 1   ┆ 2   ┆ xyz      │
        └─────┴─────┴──────────┘

        The `${1}` form is used to disambiguate the group reference from surrounding
        text.

        >>> df = pl.DataFrame({"hat": [1], "hut": [2]}).with_row_index()
        >>> df.with_columns(pl.all().name.replace(r"^h(.)t", "s$1m"))  # doctest: +SKIP
        # ComputeError: the name 's' passed to `LazyFrame.with_columns` is duplicate

        >>> df.with_columns(pl.all().name.replace(r"^h(.)t", "s${1}m"))
        shape: (1, 5)
        ┌───────┬─────┬─────┬─────┬─────┐
        │ index ┆ hat ┆ hut ┆ sam ┆ sum │
        │ ---   ┆ --- ┆ --- ┆ --- ┆ --- │
        │ u32   ┆ i64 ┆ i64 ┆ i64 ┆ i64 │
        ╞═══════╪═════╪═════╪═════╪═════╡
        │ 0     ┆ 1   ┆ 2   ┆ 1   ┆ 2   │
        └───────┴─────┴─────┴─────┴─────┘
        """
        return wrap_expr(self._pyexpr.name_replace(pattern, value, literal))

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
        return wrap_expr(self._pyexpr.name_suffix_fields(suffix))
