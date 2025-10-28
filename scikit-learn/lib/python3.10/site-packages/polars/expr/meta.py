from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import polars._reexport as pl
from polars._utils.deprecation import deprecated
from polars._utils.serde import serialize_polars_object
from polars._utils.various import display_dot_graph
from polars._utils.wrap import wrap_expr
from polars.exceptions import ComputeError

if TYPE_CHECKING:
    import sys
    from io import IOBase
    from pathlib import Path

    from polars import Expr
    from polars._typing import SchemaDict, SerializationFormat

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


class ExprMetaNameSpace:
    """Namespace for expressions on a meta level."""

    _accessor = "meta"

    def __init__(self, expr: Expr) -> None:
        self._pyexpr = expr._pyexpr

    def __eq__(self, other: ExprMetaNameSpace | Expr) -> bool:  # type: ignore[override]
        return self._pyexpr.meta_eq(other._pyexpr)

    def __ne__(self, other: ExprMetaNameSpace | Expr) -> bool:  # type: ignore[override]
        return not self == other

    def eq(self, other: ExprMetaNameSpace | Expr) -> bool:
        """
        Indicate if this expression is the same as another expression.

        Examples
        --------
        >>> foo_bar = pl.col("foo").alias("bar")
        >>> foo = pl.col("foo")
        >>> foo_bar.meta.eq(foo)
        False
        >>> foo_bar2 = pl.col("foo").alias("bar")
        >>> foo_bar.meta.eq(foo_bar2)
        True
        """
        return self._pyexpr.meta_eq(other._pyexpr)

    def ne(self, other: ExprMetaNameSpace | Expr) -> bool:
        """
        Indicate if this expression is NOT the same as another expression.

        Examples
        --------
        >>> foo_bar = pl.col("foo").alias("bar")
        >>> foo = pl.col("foo")
        >>> foo_bar.meta.ne(foo)
        True
        >>> foo_bar2 = pl.col("foo").alias("bar")
        >>> foo_bar.meta.ne(foo_bar2)
        False
        """
        return not self.eq(other)

    def has_multiple_outputs(self) -> bool:
        """
        Indicate if this expression expands into multiple expressions.

        Examples
        --------
        >>> e = pl.col(["a", "b"]).name.suffix("_foo")
        >>> e.meta.has_multiple_outputs()
        True
        """
        return self._pyexpr.meta_has_multiple_outputs()

    def is_column(self) -> bool:
        r"""
        Indicate if this expression is a basic (non-regex) unaliased column.

        Examples
        --------
        >>> e = pl.col("foo")
        >>> e.meta.is_column()
        True
        >>> e = pl.col("foo") * pl.col("bar")
        >>> e.meta.is_column()
        False
        >>> e = pl.col(r"^col.*\d+$")
        >>> e.meta.is_column()
        False
        """
        return self._pyexpr.meta_is_column()

    def is_regex_projection(self) -> bool:
        """
        Indicate if this expression expands to columns that match a regex pattern.

        Examples
        --------
        >>> e = pl.col("^.*$").name.prefix("foo_")
        >>> e.meta.is_regex_projection()
        True
        """
        return self._pyexpr.meta_is_regex_projection()

    def is_column_selection(self, *, allow_aliasing: bool = False) -> bool:
        """
        Indicate if this expression only selects columns (optionally with aliasing).

        This can include bare columns, columns matched by regex or dtype, selectors
        and exclude ops, and (optionally) column/expression aliasing.

        .. versionadded:: 0.20.30

        Parameters
        ----------
        allow_aliasing
            If False (default), any aliasing is not considered to be column selection.
            Set True to allow for column selection that also includes aliasing.

        Examples
        --------
        >>> import polars.selectors as cs
        >>> e = pl.col("foo")
        >>> e.meta.is_column_selection()
        True
        >>> e = pl.col("foo").alias("bar")
        >>> e.meta.is_column_selection()
        False
        >>> e.meta.is_column_selection(allow_aliasing=True)
        True
        >>> e = pl.col("foo") * pl.col("bar")
        >>> e.meta.is_column_selection()
        False
        >>> e = cs.starts_with("foo")
        >>> e.meta.is_column_selection()
        True
        >>> e = cs.starts_with("foo").exclude("foo!")
        >>> e.meta.is_column_selection()
        True
        """
        return self._pyexpr.meta_is_column_selection(allow_aliasing)

    def is_literal(self, *, allow_aliasing: bool = False) -> bool:
        """
        Indicate if this expression is a literal value (optionally aliased).

        .. versionadded:: 1.14

        Parameters
        ----------
        allow_aliasing
            If False (default), only a bare literal will match.
            Set True to also allow for aliased literals.

        Examples
        --------
        >>> from datetime import datetime
        >>> e = pl.lit(123)
        >>> e.meta.is_literal()
        True
        >>> e = pl.lit(987.654321).alias("foo")
        >>> e.meta.is_literal()
        False
        >>> e = pl.lit(datetime.now()).alias("bar")
        >>> e.meta.is_literal(allow_aliasing=True)
        True
        """
        return self._pyexpr.meta_is_literal(allow_aliasing)

    @overload
    def output_name(self, *, raise_if_undetermined: Literal[True] = True) -> str: ...

    @overload
    def output_name(self, *, raise_if_undetermined: Literal[False]) -> str | None: ...

    def output_name(self, *, raise_if_undetermined: bool = True) -> str | None:
        """
        Get the column name that this expression would produce.

        It may not always be possible to determine the output name as that can depend
        on the schema of the context; in that case this will raise `ComputeError` if
        `raise_if_undetermined` is True (the default), or `None` otherwise.

        Examples
        --------
        >>> e = pl.col("foo") * pl.col("bar")
        >>> e.meta.output_name()
        'foo'
        >>> e_filter = pl.col("foo").filter(pl.col("bar") == 13)
        >>> e_filter.meta.output_name()
        'foo'
        >>> e_sum_over = pl.sum("foo").over("groups")
        >>> e_sum_over.meta.output_name()
        'foo'
        >>> e_sum_slice = pl.sum("foo").slice(pl.len() - 10, pl.col("bar"))
        >>> e_sum_slice.meta.output_name()
        'foo'
        >>> pl.len().meta.output_name()
        'len'
        """
        try:
            return self._pyexpr.meta_output_name()
        except ComputeError:
            if not raise_if_undetermined:
                return None
            raise

    def pop(self, *, schema: SchemaDict | None = None) -> list[Expr]:
        """
        Pop the latest expression and return the input(s) of the popped expression.

        Returns
        -------
        list of Expr
            A list of expressions which in most cases will have a unit length.
            This is not the case when an expression has multiple inputs.
            For instance in a `fold` expression.

        Examples
        --------
        >>> e = pl.col("foo") + pl.col("bar")
        >>> first = e.meta.pop()[0]
        >>> first.meta == pl.col("bar")
        True
        >>> first.meta == pl.col("foo")
        False
        """
        return [wrap_expr(e) for e in self._pyexpr.meta_pop(schema)]

    def root_names(self) -> list[str]:
        """
        Get a list with the root column name.

        Examples
        --------
        >>> e = pl.col("foo") * pl.col("bar")
        >>> e.meta.root_names()
        ['foo', 'bar']
        >>> e_filter = pl.col("foo").filter(pl.col("bar") == 13)
        >>> e_filter.meta.root_names()
        ['foo', 'bar']
        >>> e_sum_over = pl.sum("foo").over("groups")
        >>> e_sum_over.meta.root_names()
        ['foo', 'groups']
        >>> e_sum_slice = pl.sum("foo").slice(pl.len() - 10, pl.col("bar"))
        >>> e_sum_slice.meta.root_names()
        ['foo', 'bar']
        """
        return self._pyexpr.meta_root_names()

    def undo_aliases(self) -> Expr:
        """
        Undo any renaming operation like `alias` or `name.keep`.

        Examples
        --------
        >>> e = pl.col("foo").alias("bar")
        >>> e.meta.undo_aliases().meta == pl.col("foo")
        True
        >>> e = pl.col("foo").sum().over("bar")
        >>> e.name.keep().meta.undo_aliases().meta == e
        True
        """
        return wrap_expr(self._pyexpr.meta_undo_aliases())

    def as_selector(self) -> pl.Selector:
        """
        Try to turn this expression in a selector.

        Raises if the underlying expressions is not a column or selector.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
        """
        return pl.Selector._from_pyselector(self._pyexpr.into_selector())

    @overload
    def serialize(
        self, file: None = ..., *, format: Literal["binary"] = ...
    ) -> bytes: ...

    @overload
    def serialize(self, file: None = ..., *, format: Literal["json"]) -> str: ...

    @overload
    def serialize(
        self, file: IOBase | str | Path, *, format: SerializationFormat = ...
    ) -> None: ...

    def serialize(
        self,
        file: IOBase | str | Path | None = None,
        *,
        format: SerializationFormat = "binary",
    ) -> bytes | str | None:
        r"""
        Serialize this expression to a file or string in JSON format.

        Parameters
        ----------
        file
            File path to which the result should be written. If set to `None`
            (default), the output is returned as a string instead.
        format
            The format in which to serialize. Options:

            - `"binary"`: Serialize to binary format (bytes). This is the default.
            - `"json"`: Serialize to JSON format (string).

        See Also
        --------
        Expr.deserialize

        Notes
        -----
        Serialization is not stable across Polars versions: a LazyFrame serialized
        in one Polars version may not be deserializable in another Polars version.

        Examples
        --------
        Serialize the expression into a binary representation.

        >>> expr = pl.col("foo").sum().over("bar")
        >>> bytes = expr.meta.serialize()
        >>> type(bytes)
        <class 'bytes'>

        The bytes can later be deserialized back into an `Expr` object.

        >>> import io
        >>> pl.Expr.deserialize(io.BytesIO(bytes))
        <Expr ['col("foo").sum().over([col("ba…'] at ...>
        """
        if format == "binary":
            serializer = self._pyexpr.serialize_binary
        elif format == "json":
            serializer = self._pyexpr.serialize_json
        else:
            msg = f"`format` must be one of {{'binary', 'json'}}, got {format!r}"
            raise ValueError(msg)

        return serialize_polars_object(serializer, file, format)

    @overload
    def write_json(self, file: None = ...) -> str: ...

    @overload
    def write_json(self, file: IOBase | str | Path) -> None: ...

    @deprecated("`meta.write_json` was renamed; use `meta.serialize` instead")
    def write_json(self, file: IOBase | str | Path | None = None) -> str | None:
        """
        Write expression to json.

        .. deprecated:: 0.20.11
            This method has been renamed to :meth:`serialize`.
        """
        return self.serialize(file, format="json")

    @overload
    def tree_format(
        self,
        *,
        return_as_string: Literal[False] = ...,
        schema: None | SchemaDict = None,
    ) -> None: ...

    @overload
    def tree_format(
        self, *, return_as_string: Literal[True], schema: None | SchemaDict = None
    ) -> str: ...

    def tree_format(  # noqa: D417 (TODO: document schema parameter)
        self, *, return_as_string: bool = False, schema: None | SchemaDict = None
    ) -> str | None:
        """
        Format the expression as a tree.

        Parameters
        ----------
        return_as_string:
            If True, return as string rather than printing to stdout.

        Examples
        --------
        >>> e = (pl.col("foo") * pl.col("bar")).sum().over(pl.col("ham")) / 2
        >>> e.meta.tree_format(return_as_string=True)  # doctest: +SKIP
        """
        s = self._pyexpr.meta_tree_format(schema)
        if return_as_string:
            return s
        else:
            print(s)
            return None

    def show_graph(  # noqa: D417 (TODO: document schema parameter)
        self,
        *,
        show: bool = True,
        output_path: str | Path | None = None,
        raw_output: bool = False,
        figsize: tuple[float, float] = (16.0, 12.0),
        schema: None | SchemaDict = None,
    ) -> str | None:
        """
        Format the expression as a Graphviz graph.

        Note that Graphviz must be installed to render the visualization (if not
        already present, you can download it here: `<https://graphviz.org/download>`_).

        Parameters
        ----------
        show
            Show the figure.
        output_path
            Write the figure to disk.
        raw_output
            Return dot syntax. This cannot be combined with `show` and/or `output_path`.
        figsize
            Passed to matplotlib if `show == True`.

        Examples
        --------
        >>> e = (pl.col("foo") * pl.col("bar")).sum().over(pl.col("ham")) / 2
        >>> e.meta.show_graph()  # doctest: +SKIP
        """
        dot = self._pyexpr.meta_show_graph(schema)
        return display_dot_graph(
            dot=dot,
            show=show,
            output_path=output_path,
            raw_output=raw_output,
            figsize=figsize,
        )
