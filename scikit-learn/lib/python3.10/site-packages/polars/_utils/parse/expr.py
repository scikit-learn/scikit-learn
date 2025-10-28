from __future__ import annotations

import contextlib
from collections.abc import Collection, Iterable, Mapping
from typing import TYPE_CHECKING, Any

import polars._reexport as pl
from polars import functions as F
from polars._utils.various import qualified_type_name
from polars.exceptions import ComputeError

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

if TYPE_CHECKING:
    from polars import Expr
    from polars._plr import PyExpr
    from polars._typing import ColumnNameOrSelector, IntoExpr, PolarsDataType


def parse_into_expression(
    input: IntoExpr,
    *,
    str_as_lit: bool = False,
    list_as_series: bool = False,
    structify: bool = False,
    dtype: PolarsDataType | None = None,
) -> PyExpr:
    """
    Parse a single input into an expression.

    Parameters
    ----------
    input
        The input to be parsed as an expression.
    str_as_lit
        Interpret string input as a string literal. If set to `False` (default),
        strings are parsed as column names.
    list_as_series
        Interpret list input as a Series literal. If set to `False` (default),
        lists are parsed as list literals.
    structify
        Convert multi-column expressions to a single struct expression.
    dtype
        If the input is expected to resolve to a literal with a known dtype, pass
        this to the `lit` constructor.

    Returns
    -------
    PyExpr
    """
    if isinstance(input, pl.Expr):
        expr = input
        if structify:
            expr = _structify_expression(expr)

    elif isinstance(input, str) and not str_as_lit:
        expr = F.col(input)
    elif isinstance(input, list) and list_as_series:
        expr = F.lit(pl.Series(input), dtype=dtype)
    else:
        expr = F.lit(input, dtype=dtype)

    return expr._pyexpr


def _structify_expression(expr: Expr) -> Expr:
    unaliased_expr = expr.meta.undo_aliases()
    if unaliased_expr.meta.has_multiple_outputs():
        try:
            expr_name = expr.meta.output_name()
        except ComputeError:
            expr = F.struct(expr)
        else:
            expr = F.struct(unaliased_expr).alias(expr_name)
    return expr


def parse_into_list_of_expressions(
    *inputs: IntoExpr | Iterable[IntoExpr],
    __structify: bool = False,
    **named_inputs: IntoExpr,
) -> list[PyExpr]:
    """
    Parse multiple inputs into a list of expressions.

    Parameters
    ----------
    *inputs
        Inputs to be parsed as expressions, specified as positional arguments.
    **named_inputs
        Additional inputs to be parsed as expressions, specified as keyword arguments.
        The expressions will be renamed to the keyword used.
    __structify
        Convert multi-column expressions to a single struct expression.

    Returns
    -------
    list of PyExpr
    """
    exprs = _parse_positional_inputs(inputs, structify=__structify)  # type: ignore[arg-type]
    if named_inputs:
        named_exprs = _parse_named_inputs(named_inputs, structify=__structify)
        exprs.extend(named_exprs)

    return exprs


def parse_into_selector(
    i: ColumnNameOrSelector,
    *,
    strict: bool = True,
) -> pl.Selector:
    if isinstance(i, str):
        import polars.selectors as cs

        return cs.by_name([i], require_all=strict)
    elif isinstance(i, pl.Selector):
        return i
    elif isinstance(i, pl.Expr):
        return i.meta.as_selector()
    else:
        msg = f"cannot turn {qualified_type_name(i)!r} into selector"
        raise TypeError(msg)


def parse_list_into_selector(
    inputs: ColumnNameOrSelector | Collection[ColumnNameOrSelector],
    *,
    strict: bool = True,
) -> pl.Selector:
    if isinstance(inputs, Collection) and not isinstance(inputs, str):
        import polars.selectors as cs

        columns = list(filter(lambda i: isinstance(i, str), inputs))
        selector = cs.by_name(columns, require_all=strict)  # type: ignore[arg-type]

        if len(columns) == len(inputs):
            return selector

        # A bit cleaner
        if len(columns) == 0:
            selector = cs.empty()

        for i in inputs:
            selector |= parse_into_selector(i, strict=strict)
        return selector
    else:
        return parse_into_selector(inputs, strict=strict)


def _parse_positional_inputs(
    inputs: tuple[IntoExpr, ...] | tuple[Iterable[IntoExpr]],
    *,
    structify: bool = False,
) -> list[PyExpr]:
    inputs_iter = _parse_inputs_as_iterable(inputs)
    return [parse_into_expression(e, structify=structify) for e in inputs_iter]


def _parse_inputs_as_iterable(
    inputs: tuple[Any, ...] | tuple[Iterable[Any]],
) -> Iterable[Any]:
    if not inputs:
        return []

    # Ensures that the outermost element cannot be a Dictionary (as an iterable)
    if len(inputs) == 1 and isinstance(inputs[0], Mapping):
        msg = (
            "Cannot pass a dictionary as a single positional argument.\n"
            "If you merely want the *keys*, use:\n"
            "  • df.method(*your_dict.keys())\n"
            "If you need the key value pairs, use one of:\n"
            "  • unpack as keywords:    df.method(**your_dict)\n"
            "  • build expressions:     df.method(expr.alias(k) for k, expr in your_dict.items())"
        )
        raise TypeError(msg)

    # Treat elements of a single iterable as separate inputs
    if len(inputs) == 1 and _is_iterable(inputs[0]):
        return inputs[0]

    return inputs


def _is_iterable(input: Any | Iterable[Any]) -> bool:
    return isinstance(input, Iterable) and not isinstance(
        input, (str, bytes, pl.Series)
    )


def _parse_named_inputs(
    named_inputs: dict[str, IntoExpr], *, structify: bool = False
) -> Iterable[PyExpr]:
    for name, input in named_inputs.items():
        yield parse_into_expression(input, structify=structify).alias(name)


def parse_predicates_constraints_into_expression(
    *predicates: IntoExpr | Iterable[IntoExpr],
    **constraints: Any,
) -> PyExpr:
    """
    Parse predicates and constraints into a single expression.

    The result is an AND-reduction of all inputs.

    Parameters
    ----------
    *predicates
        Predicates to be parsed, specified as positional arguments.
    **constraints
        Constraints to be parsed, specified as keyword arguments.
        These will be converted to predicates of the form "keyword equals input value".

    Returns
    -------
    PyExpr
    """
    all_predicates = _parse_positional_inputs(predicates)  # type: ignore[arg-type]

    if constraints:
        constraint_predicates = _parse_constraints(constraints)
        all_predicates.extend(constraint_predicates)

    return _combine_predicates(all_predicates)


def _parse_constraints(constraints: dict[str, IntoExpr]) -> Iterable[PyExpr]:
    for name, value in constraints.items():
        yield F.col(name).eq(value)._pyexpr


def _combine_predicates(predicates: list[PyExpr]) -> PyExpr:
    if not predicates:
        msg = "at least one predicate or constraint must be provided"
        raise TypeError(msg)

    if len(predicates) == 1:
        return predicates[0]

    return plr.all_horizontal(predicates)
