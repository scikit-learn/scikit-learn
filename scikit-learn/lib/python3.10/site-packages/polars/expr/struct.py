from __future__ import annotations

import os
from typing import TYPE_CHECKING

from polars._utils.parse import parse_into_list_of_expressions
from polars._utils.various import qualified_type_name
from polars._utils.wrap import wrap_expr

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from polars import Expr
    from polars._typing import IntoExpr


class ExprStructNameSpace:
    """Namespace for struct related expressions."""

    _accessor = "struct"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def __getitem__(self, item: str | int) -> Expr:
        if isinstance(item, str):
            return self.field(item)
        elif isinstance(item, int):
            return wrap_expr(self._pyexpr.struct_field_by_index(item))
        else:
            msg = f"expected type 'int | str', got {qualified_type_name(item)!r} ({item!r})"
            raise TypeError(msg)

    def field(self, name: str | list[str], *more_names: str) -> Expr:
        """
        Retrieve one or multiple `Struct` field(s) as a new Series.

        Parameters
        ----------
        name
            Name of the struct field to retrieve.
        *more_names
            Additional struct field names.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "aaa": [1, 2],
        ...         "bbb": ["ab", "cd"],
        ...         "ccc": [True, None],
        ...         "ddd": [[1, 2], [3]],
        ...     }
        ... ).select(pl.struct("aaa", "bbb", "ccc", "ddd").alias("struct_col"))
        >>> df
        shape: (2, 1)
        ┌──────────────────────┐
        │ struct_col           │
        │ ---                  │
        │ struct[4]            │
        ╞══════════════════════╡
        │ {1,"ab",true,[1, 2]} │
        │ {2,"cd",null,[3]}    │
        └──────────────────────┘

        Retrieve struct field(s) as Series:

        >>> df.select(pl.col("struct_col").struct.field("bbb"))
        shape: (2, 1)
        ┌─────┐
        │ bbb │
        │ --- │
        │ str │
        ╞═════╡
        │ ab  │
        │ cd  │
        └─────┘

        >>> df.select(
        ...     pl.col("struct_col").struct.field("bbb"),
        ...     pl.col("struct_col").struct.field("ddd"),
        ... )
        shape: (2, 2)
        ┌─────┬───────────┐
        │ bbb ┆ ddd       │
        │ --- ┆ ---       │
        │ str ┆ list[i64] │
        ╞═════╪═══════════╡
        │ ab  ┆ [1, 2]    │
        │ cd  ┆ [3]       │
        └─────┴───────────┘

        Use wildcard expansion:

        >>> df.select(pl.col("struct_col").struct.field("*"))
        shape: (2, 4)
        ┌─────┬─────┬──────┬───────────┐
        │ aaa ┆ bbb ┆ ccc  ┆ ddd       │
        │ --- ┆ --- ┆ ---  ┆ ---       │
        │ i64 ┆ str ┆ bool ┆ list[i64] │
        ╞═════╪═════╪══════╪═══════════╡
        │ 1   ┆ ab  ┆ true ┆ [1, 2]    │
        │ 2   ┆ cd  ┆ null ┆ [3]       │
        └─────┴─────┴──────┴───────────┘

        Retrieve multiple fields by name:

        >>> df.select(pl.col("struct_col").struct.field("aaa", "bbb"))
        shape: (2, 2)
        ┌─────┬─────┐
        │ aaa ┆ bbb │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ ab  │
        │ 2   ┆ cd  │
        └─────┴─────┘

        Retrieve multiple fields by regex expansion:

        >>> df.select(pl.col("struct_col").struct.field("^a.*|b.*$"))
        shape: (2, 2)
        ┌─────┬─────┐
        │ aaa ┆ bbb │
        │ --- ┆ --- │
        │ i64 ┆ str │
        ╞═════╪═════╡
        │ 1   ┆ ab  │
        │ 2   ┆ cd  │
        └─────┴─────┘

        Notes
        -----
        The `struct` namespace has implemented `__getitem__`
        so you can also access fields by index:

        >>> df.select(pl.col("struct_col").struct[1])
        shape: (2, 1)
        ┌─────┐
        │ bbb │
        │ --- │
        │ str │
        ╞═════╡
        │ ab  │
        │ cd  │
        └─────┘
        """
        if more_names:
            name = [*([name] if isinstance(name, str) else name), *more_names]
        if isinstance(name, list):
            return wrap_expr(self._pyexpr.struct_multiple_fields(name))

        return wrap_expr(self._pyexpr.struct_field_by_name(name))

    def unnest(self) -> Expr:
        """
        Expand the struct into its individual fields.

        Alias for `Expr.struct.field("*")`.

        >>> df = pl.DataFrame(
        ...     {
        ...         "aaa": [1, 2],
        ...         "bbb": ["ab", "cd"],
        ...         "ccc": [True, None],
        ...         "ddd": [[1, 2], [3]],
        ...     }
        ... ).select(pl.struct("aaa", "bbb", "ccc", "ddd").alias("struct_col"))
        >>> df
        shape: (2, 1)
        ┌──────────────────────┐
        │ struct_col           │
        │ ---                  │
        │ struct[4]            │
        ╞══════════════════════╡
        │ {1,"ab",true,[1, 2]} │
        │ {2,"cd",null,[3]}    │
        └──────────────────────┘
        >>> df.select(pl.col("struct_col").struct.unnest())
        shape: (2, 4)
        ┌─────┬─────┬──────┬───────────┐
        │ aaa ┆ bbb ┆ ccc  ┆ ddd       │
        │ --- ┆ --- ┆ ---  ┆ ---       │
        │ i64 ┆ str ┆ bool ┆ list[i64] │
        ╞═════╪═════╪══════╪═══════════╡
        │ 1   ┆ ab  ┆ true ┆ [1, 2]    │
        │ 2   ┆ cd  ┆ null ┆ [3]       │
        └─────┴─────┴──────┴───────────┘
        """
        return self.field("*")

    def rename_fields(self, names: Sequence[str]) -> Expr:
        """
        Rename the fields of the struct.

        Parameters
        ----------
        names
            New names, given in the same order as the struct's fields.

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "aaa": [1, 2],
        ...         "bbb": ["ab", "cd"],
        ...         "ccc": [True, None],
        ...         "ddd": [[1, 2], [3]],
        ...     }
        ... ).select(pl.struct("aaa", "bbb", "ccc", "ddd").alias("struct_col"))
        >>> df
        shape: (2, 1)
        ┌──────────────────────┐
        │ struct_col           │
        │ ---                  │
        │ struct[4]            │
        ╞══════════════════════╡
        │ {1,"ab",true,[1, 2]} │
        │ {2,"cd",null,[3]}    │
        └──────────────────────┘

        >>> df.unnest("struct_col")
        shape: (2, 4)
        ┌─────┬─────┬──────┬───────────┐
        │ aaa ┆ bbb ┆ ccc  ┆ ddd       │
        │ --- ┆ --- ┆ ---  ┆ ---       │
        │ i64 ┆ str ┆ bool ┆ list[i64] │
        ╞═════╪═════╪══════╪═══════════╡
        │ 1   ┆ ab  ┆ true ┆ [1, 2]    │
        │ 2   ┆ cd  ┆ null ┆ [3]       │
        └─────┴─────┴──────┴───────────┘

        Rename fields:

        >>> df = df.select(
        ...     pl.col("struct_col").struct.rename_fields(["www", "xxx", "yyy", "zzz"])
        ... )
        >>> df.unnest("struct_col")
        shape: (2, 4)
        ┌─────┬─────┬──────┬───────────┐
        │ www ┆ xxx ┆ yyy  ┆ zzz       │
        │ --- ┆ --- ┆ ---  ┆ ---       │
        │ i64 ┆ str ┆ bool ┆ list[i64] │
        ╞═════╪═════╪══════╪═══════════╡
        │ 1   ┆ ab  ┆ true ┆ [1, 2]    │
        │ 2   ┆ cd  ┆ null ┆ [3]       │
        └─────┴─────┴──────┴───────────┘

        Following a rename, the previous field names (obviously) cannot be referenced:

        >>> df.select(pl.col("struct_col").struct.field("aaa"))  # doctest: +SKIP
        StructFieldNotFoundError: aaa
        """
        return wrap_expr(self._pyexpr.struct_rename_fields(names))

    def json_encode(self) -> Expr:
        """
        Convert this struct to a string column with json values.

        Examples
        --------
        >>> pl.DataFrame(
        ...     {"a": [{"a": [1, 2], "b": [45]}, {"a": [9, 1, 3], "b": None}]}
        ... ).with_columns(pl.col("a").struct.json_encode().alias("encoded"))
        shape: (2, 2)
        ┌──────────────────┬────────────────────────┐
        │ a                ┆ encoded                │
        │ ---              ┆ ---                    │
        │ struct[2]        ┆ str                    │
        ╞══════════════════╪════════════════════════╡
        │ {[1, 2],[45]}    ┆ {"a":[1,2],"b":[45]}   │
        │ {[9, 1, 3],null} ┆ {"a":[9,1,3],"b":null} │
        └──────────────────┴────────────────────────┘
        """
        return wrap_expr(self._pyexpr.struct_json_encode())

    def with_fields(
        self,
        *exprs: IntoExpr | Iterable[IntoExpr],
        **named_exprs: IntoExpr,
    ) -> Expr:
        """
        Add or overwrite fields of this struct.

        This is similar to `with_columns` on `DataFrame`.

        .. versionadded:: 0.20.27

        Examples
        --------
        >>> df = pl.DataFrame(
        ...     {
        ...         "coords": [{"x": 1, "y": 4}, {"x": 4, "y": 9}, {"x": 9, "y": 16}],
        ...         "multiply": [10, 2, 3],
        ...     }
        ... )
        >>> df
        shape: (3, 2)
        ┌───────────┬──────────┐
        │ coords    ┆ multiply │
        │ ---       ┆ ---      │
        │ struct[2] ┆ i64      │
        ╞═══════════╪══════════╡
        │ {1,4}     ┆ 10       │
        │ {4,9}     ┆ 2        │
        │ {9,16}    ┆ 3        │
        └───────────┴──────────┘
        >>> df = df.with_columns(
        ...     pl.col("coords").struct.with_fields(
        ...         pl.field("x").sqrt(),
        ...         y_mul=pl.field("y") * pl.col("multiply"),
        ...     )
        ... )
        >>> df
        shape: (3, 2)
        ┌─────────────┬──────────┐
        │ coords      ┆ multiply │
        │ ---         ┆ ---      │
        │ struct[3]   ┆ i64      │
        ╞═════════════╪══════════╡
        │ {1.0,4,40}  ┆ 10       │
        │ {2.0,9,18}  ┆ 2        │
        │ {3.0,16,48} ┆ 3        │
        └─────────────┴──────────┘
        >>> df.unnest("coords")
        shape: (3, 4)
        ┌─────┬─────┬───────┬──────────┐
        │ x   ┆ y   ┆ y_mul ┆ multiply │
        │ --- ┆ --- ┆ ---   ┆ ---      │
        │ f64 ┆ i64 ┆ i64   ┆ i64      │
        ╞═════╪═════╪═══════╪══════════╡
        │ 1.0 ┆ 4   ┆ 40    ┆ 10       │
        │ 2.0 ┆ 9   ┆ 18    ┆ 2        │
        │ 3.0 ┆ 16  ┆ 48    ┆ 3        │
        └─────┴─────┴───────┴──────────┘

        Parameters
        ----------
        *exprs
            Field(s) to add, specified as positional arguments.
            Accepts expression input. Strings are parsed as column names, other
            non-expression inputs are parsed as literals.
        **named_exprs
            Additional fields to add, specified as keyword arguments.
            The columns will be renamed to the keyword used.

        See Also
        --------
        field
        """
        structify = bool(int(os.environ.get("POLARS_AUTO_STRUCTIFY", 0)))

        pyexprs = parse_into_list_of_expressions(
            *exprs, **named_exprs, __structify=structify
        )

        return wrap_expr(self._pyexpr.struct_with_fields(pyexprs))
