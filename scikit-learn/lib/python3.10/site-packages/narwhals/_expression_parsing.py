# Utilities for expression parsing
# Useful for backends which don't have any concept of expressions, such
# and pandas or PyArrow.
# ! Any change to this module will trigger the pyspark and pyspark-connect tests in CI
from __future__ import annotations

from enum import Enum, auto
from itertools import chain
from typing import TYPE_CHECKING, Any, Callable, Literal, TypeVar

from narwhals._utils import is_compliant_expr, zip_strict
from narwhals.dependencies import is_narwhals_series, is_numpy_array, is_numpy_array_1d
from narwhals.exceptions import InvalidOperationError, MultiOutputExpressionError

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Never, TypeIs

    from narwhals._compliant import CompliantExpr, CompliantFrameT
    from narwhals._compliant.typing import (
        AliasNames,
        CompliantExprAny,
        CompliantFrameAny,
        CompliantNamespaceAny,
        EvalNames,
    )
    from narwhals.expr import Expr
    from narwhals.series import Series
    from narwhals.typing import IntoExpr, NonNestedLiteral, _1DArray

    T = TypeVar("T")


def is_expr(obj: Any) -> TypeIs[Expr]:
    """Check whether `obj` is a Narwhals Expr."""
    from narwhals.expr import Expr

    return isinstance(obj, Expr)


def is_series(obj: Any) -> TypeIs[Series[Any]]:
    """Check whether `obj` is a Narwhals Expr."""
    from narwhals.series import Series

    return isinstance(obj, Series)


def is_into_expr_eager(obj: Any) -> TypeIs[Expr | Series[Any] | str | _1DArray]:
    from narwhals.expr import Expr
    from narwhals.series import Series

    return isinstance(obj, (Series, Expr, str)) or is_numpy_array_1d(obj)


def combine_evaluate_output_names(
    *exprs: CompliantExpr[CompliantFrameT, Any],
) -> EvalNames[CompliantFrameT]:
    # Follow left-hand-rule for naming. E.g. `nw.sum_horizontal(expr1, expr2)` takes the
    # first name of `expr1`.
    if not is_compliant_expr(exprs[0]):  # pragma: no cover
        msg = f"Safety assertion failed, expected expression, got: {type(exprs[0])}. Please report a bug."
        raise AssertionError(msg)

    def evaluate_output_names(df: CompliantFrameT) -> Sequence[str]:
        return exprs[0]._evaluate_output_names(df)[:1]

    return evaluate_output_names


def combine_alias_output_names(*exprs: CompliantExprAny) -> AliasNames | None:
    # Follow left-hand-rule for naming. E.g. `nw.sum_horizontal(expr1.alias(alias), expr2)` takes the
    # aliasing function of `expr1` and apply it to the first output name of `expr1`.
    if exprs[0]._alias_output_names is None:
        return None

    def alias_output_names(names: Sequence[str]) -> Sequence[str]:
        return exprs[0]._alias_output_names(names)[:1]  # type: ignore[misc]

    return alias_output_names


def evaluate_output_names_and_aliases(
    expr: CompliantExprAny, df: CompliantFrameAny, exclude: Sequence[str]
) -> tuple[Sequence[str], Sequence[str]]:
    output_names = expr._evaluate_output_names(df)
    aliases = (
        output_names
        if expr._alias_output_names is None
        else expr._alias_output_names(output_names)
    )
    if exclude:
        assert expr._metadata is not None  # noqa: S101
        if expr._metadata.expansion_kind.is_multi_unnamed():
            output_names, aliases = zip_strict(
                *[
                    (x, alias)
                    for x, alias in zip_strict(output_names, aliases)
                    if x not in exclude
                ]
            )
    return output_names, aliases


class ExprKind(Enum):
    """Describe which kind of expression we are dealing with."""

    LITERAL = auto()
    """e.g. `nw.lit(1)`"""

    AGGREGATION = auto()
    """Reduces to a single value, not affected by row order, e.g. `nw.col('a').mean()`"""

    ORDERABLE_AGGREGATION = auto()
    """Reduces to a single value, affected by row order, e.g. `nw.col('a').arg_max()`"""

    ELEMENTWISE = auto()
    """Preserves length, can operate without context for surrounding rows, e.g. `nw.col('a').abs()`."""

    ORDERABLE_WINDOW = auto()
    """Depends on the rows around it and on their order, e.g. `diff`."""

    WINDOW = auto()
    """Depends on the rows around it and possibly their order, e.g. `rank`."""

    FILTRATION = auto()
    """Changes length, not affected by row order, e.g. `drop_nulls`."""

    ORDERABLE_FILTRATION = auto()
    """Changes length, affected by row order, e.g. `tail`."""

    OVER = auto()
    """Results from calling `.over` on expression."""

    UNKNOWN = auto()
    """Based on the information we have, we can't determine the ExprKind."""

    @property
    def is_scalar_like(self) -> bool:
        return self in {ExprKind.LITERAL, ExprKind.AGGREGATION}

    @property
    def is_orderable_window(self) -> bool:
        return self in {ExprKind.ORDERABLE_WINDOW, ExprKind.ORDERABLE_AGGREGATION}

    @classmethod
    def from_expr(cls, obj: Expr) -> ExprKind:
        meta = obj._metadata
        if meta.is_literal:
            return ExprKind.LITERAL
        if meta.is_scalar_like:
            return ExprKind.AGGREGATION
        if meta.is_elementwise:
            return ExprKind.ELEMENTWISE
        return ExprKind.UNKNOWN

    @classmethod
    def from_into_expr(
        cls, obj: IntoExpr | NonNestedLiteral | _1DArray, *, str_as_lit: bool
    ) -> ExprKind:
        if is_expr(obj):
            return cls.from_expr(obj)
        if (
            is_narwhals_series(obj)
            or is_numpy_array(obj)
            or (isinstance(obj, str) and not str_as_lit)
        ):
            return ExprKind.ELEMENTWISE
        return ExprKind.LITERAL


def is_scalar_like(
    obj: ExprKind,
) -> TypeIs[Literal[ExprKind.LITERAL, ExprKind.AGGREGATION]]:
    return obj.is_scalar_like


class ExpansionKind(Enum):
    """Describe what kind of expansion the expression performs."""

    SINGLE = auto()
    """e.g. `nw.col('a'), nw.sum_horizontal(nw.all())`"""

    MULTI_NAMED = auto()
    """e.g. `nw.col('a', 'b')`"""

    MULTI_UNNAMED = auto()
    """e.g. `nw.all()`, nw.nth(0, 1)"""

    def is_multi_unnamed(self) -> bool:
        return self is ExpansionKind.MULTI_UNNAMED

    def is_multi_output(self) -> bool:
        return self in {ExpansionKind.MULTI_NAMED, ExpansionKind.MULTI_UNNAMED}

    def __and__(self, other: ExpansionKind) -> Literal[ExpansionKind.MULTI_UNNAMED]:
        if self is ExpansionKind.MULTI_UNNAMED and other is ExpansionKind.MULTI_UNNAMED:
            # e.g. nw.selectors.all() - nw.selectors.numeric().
            return ExpansionKind.MULTI_UNNAMED
        # Don't attempt anything more complex, keep it simple and raise in the face of ambiguity.
        msg = f"Unsupported ExpansionKind combination, got {self} and {other}, please report a bug."  # pragma: no cover
        raise AssertionError(msg)  # pragma: no cover


class ExprMetadata:
    """Expression metadata.

    Parameters:
        expansion_kind: What kind of expansion the expression performs.
        has_windows: Whether it already contains window functions.
        is_elementwise: Whether it can operate row-by-row without context
            of the other rows around it.
        is_literal: Whether it is just a literal wrapped in an expression.
        is_scalar_like: Whether it is a literal or an aggregation.
        last_node: The ExprKind of the last node.
        n_orderable_ops: The number of order-dependent operations. In the
            lazy case, this number must be `0` by the time the expression
            is evaluated.
        preserves_length: Whether the expression preserves the input length.
    """

    __slots__ = (
        "expansion_kind",
        "has_windows",
        "is_elementwise",
        "is_literal",
        "is_scalar_like",
        "last_node",
        "n_orderable_ops",
        "preserves_length",
    )

    def __init__(
        self,
        expansion_kind: ExpansionKind,
        last_node: ExprKind,
        *,
        has_windows: bool = False,
        n_orderable_ops: int = 0,
        preserves_length: bool = True,
        is_elementwise: bool = True,
        is_scalar_like: bool = False,
        is_literal: bool = False,
    ) -> None:
        if is_literal:
            assert is_scalar_like  # noqa: S101  # debug assertion
        if is_elementwise:
            assert preserves_length  # noqa: S101  # debug assertion
        self.expansion_kind: ExpansionKind = expansion_kind
        self.last_node: ExprKind = last_node
        self.has_windows: bool = has_windows
        self.n_orderable_ops: int = n_orderable_ops
        self.is_elementwise: bool = is_elementwise
        self.preserves_length: bool = preserves_length
        self.is_scalar_like: bool = is_scalar_like
        self.is_literal: bool = is_literal

    def __init_subclass__(cls, /, *args: Any, **kwds: Any) -> Never:  # pragma: no cover
        msg = f"Cannot subclass {cls.__name__!r}"
        raise TypeError(msg)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"ExprMetadata(\n"
            f"  expansion_kind: {self.expansion_kind},\n"
            f"  last_node: {self.last_node},\n"
            f"  has_windows: {self.has_windows},\n"
            f"  n_orderable_ops: {self.n_orderable_ops},\n"
            f"  is_elementwise: {self.is_elementwise},\n"
            f"  preserves_length: {self.preserves_length},\n"
            f"  is_scalar_like: {self.is_scalar_like},\n"
            f"  is_literal: {self.is_literal},\n"
            ")"
        )

    @property
    def is_filtration(self) -> bool:
        return not self.preserves_length and not self.is_scalar_like

    def with_aggregation(self) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply aggregations to scalar-like expressions."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.AGGREGATION,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=True,
            is_literal=False,
        )

    def with_orderable_aggregation(self) -> ExprMetadata:
        # Deprecated, used only in stable.v1.
        if self.is_scalar_like:  # pragma: no cover
            msg = "Can't apply aggregations to scalar-like expressions."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.ORDERABLE_AGGREGATION,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=True,
            is_literal=False,
        )

    def with_elementwise_op(self) -> ExprMetadata:
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.ELEMENTWISE,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=self.preserves_length,
            is_elementwise=self.is_elementwise,
            is_scalar_like=self.is_scalar_like,
            is_literal=self.is_literal,
        )

    def with_window(self) -> ExprMetadata:
        # Window function which may (but doesn't have to) be used with `over(order_by=...)`.
        if self.is_scalar_like:
            msg = "Can't apply window (e.g. `rank`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.WINDOW,
            has_windows=self.has_windows,
            # The function isn't order-dependent (but, users can still use `order_by` if they wish!),
            # so we don't increment `n_orderable_ops`.
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=self.preserves_length,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    def with_orderable_window(self) -> ExprMetadata:
        # Window function which must be used with `over(order_by=...)`.
        if self.is_scalar_like:
            msg = "Can't apply orderable window (e.g. `diff`, `shift`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.ORDERABLE_WINDOW,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=self.preserves_length,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    def with_ordered_over(self) -> ExprMetadata:
        if self.has_windows:
            msg = "Cannot nest `over` statements."
            raise InvalidOperationError(msg)
        if self.is_elementwise or self.is_filtration:
            msg = (
                "Cannot use `over` on expressions which are elementwise\n"
                "(e.g. `abs`) or which change length (e.g. `drop_nulls`)."
            )
            raise InvalidOperationError(msg)
        n_orderable_ops = self.n_orderable_ops
        if not n_orderable_ops and self.last_node is not ExprKind.WINDOW:
            msg = (
                "Cannot use `order_by` in `over` on expression which isn't orderable.\n"
                "If your expression is orderable, then make sure that `over(order_by=...)`\n"
                "comes immediately after the order-dependent expression.\n\n"
                "Hint: instead of\n"
                "  - `(nw.col('price').diff() + 1).over(order_by='date')`\n"
                "write:\n"
                "  + `nw.col('price').diff().over(order_by='date') + 1`\n"
            )
            raise InvalidOperationError(msg)
        if self.last_node.is_orderable_window:
            n_orderable_ops -= 1
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.OVER,
            has_windows=True,
            n_orderable_ops=n_orderable_ops,
            preserves_length=True,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    def with_partitioned_over(self) -> ExprMetadata:
        if self.has_windows:
            msg = "Cannot nest `over` statements."
            raise InvalidOperationError(msg)
        if self.is_elementwise or self.is_filtration:
            msg = (
                "Cannot use `over` on expressions which are elementwise\n"
                "(e.g. `abs`) or which change length (e.g. `drop_nulls`)."
            )
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.OVER,
            has_windows=True,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=True,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    def with_filtration(self) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply filtration (e.g. `drop_nulls`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.FILTRATION,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    def with_orderable_filtration(self) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply filtration (e.g. `drop_nulls`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            ExprKind.ORDERABLE_FILTRATION,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
        )

    @staticmethod
    def aggregation() -> ExprMetadata:
        return ExprMetadata(
            ExpansionKind.SINGLE,
            ExprKind.AGGREGATION,
            is_elementwise=False,
            preserves_length=False,
            is_scalar_like=True,
        )

    @staticmethod
    def literal() -> ExprMetadata:
        return ExprMetadata(
            ExpansionKind.SINGLE,
            ExprKind.LITERAL,
            is_elementwise=False,
            preserves_length=False,
            is_literal=True,
            is_scalar_like=True,
        )

    @staticmethod
    def selector_single() -> ExprMetadata:
        # e.g. `nw.col('a')`, `nw.nth(0)`
        return ExprMetadata(ExpansionKind.SINGLE, ExprKind.ELEMENTWISE)

    @staticmethod
    def selector_multi_named() -> ExprMetadata:
        # e.g. `nw.col('a', 'b')`
        return ExprMetadata(ExpansionKind.MULTI_NAMED, ExprKind.ELEMENTWISE)

    @staticmethod
    def selector_multi_unnamed() -> ExprMetadata:
        # e.g. `nw.all()`
        return ExprMetadata(ExpansionKind.MULTI_UNNAMED, ExprKind.ELEMENTWISE)

    @classmethod
    def from_binary_op(cls, lhs: Expr, rhs: IntoExpr, /) -> ExprMetadata:
        # We may be able to allow multi-output rhs in the future:
        # https://github.com/narwhals-dev/narwhals/issues/2244.
        return combine_metadata(
            lhs, rhs, str_as_lit=True, allow_multi_output=False, to_single_output=False
        )

    @classmethod
    def from_horizontal_op(cls, *exprs: IntoExpr) -> ExprMetadata:
        return combine_metadata(
            *exprs, str_as_lit=False, allow_multi_output=True, to_single_output=True
        )


def combine_metadata(
    *args: IntoExpr | object | None,
    str_as_lit: bool,
    allow_multi_output: bool,
    to_single_output: bool,
) -> ExprMetadata:
    """Combine metadata from `args`.

    Arguments:
        args: Arguments, maybe expressions, literals, or Series.
        str_as_lit: Whether to interpret strings as literals or as column names.
        allow_multi_output: Whether to allow multi-output inputs.
        to_single_output: Whether the result is always single-output, regardless
            of the inputs (e.g. `nw.sum_horizontal`).
    """
    n_filtrations = 0
    result_expansion_kind = ExpansionKind.SINGLE
    result_has_windows = False
    result_n_orderable_ops = 0
    # result preserves length if at least one input does
    result_preserves_length = False
    # result is elementwise if all inputs are elementwise
    result_is_elementwise = True
    # result is scalar-like if all inputs are scalar-like
    result_is_scalar_like = True
    # result is literal if all inputs are literal
    result_is_literal = True

    for i, arg in enumerate(args):
        if (isinstance(arg, str) and not str_as_lit) or is_series(arg):
            result_preserves_length = True
            result_is_scalar_like = False
            result_is_literal = False
        elif is_expr(arg):
            metadata = arg._metadata
            if metadata.expansion_kind.is_multi_output():
                expansion_kind = metadata.expansion_kind
                if i > 0 and not allow_multi_output:
                    # Left-most argument is always allowed to be multi-output.
                    msg = (
                        "Multi-output expressions (e.g. nw.col('a', 'b'), nw.all()) "
                        "are not supported in this context."
                    )
                    raise MultiOutputExpressionError(msg)
                if not to_single_output:
                    result_expansion_kind = (
                        result_expansion_kind & expansion_kind
                        if i > 0
                        else expansion_kind
                    )

            result_has_windows |= metadata.has_windows
            result_n_orderable_ops += metadata.n_orderable_ops
            result_preserves_length |= metadata.preserves_length
            result_is_elementwise &= metadata.is_elementwise
            result_is_scalar_like &= metadata.is_scalar_like
            result_is_literal &= metadata.is_literal
            n_filtrations += int(metadata.is_filtration)

    if n_filtrations > 1:
        msg = "Length-changing expressions can only be used in isolation, or followed by an aggregation"
        raise InvalidOperationError(msg)
    if result_preserves_length and n_filtrations:
        msg = "Cannot combine length-changing expressions with length-preserving ones or aggregations"
        raise InvalidOperationError(msg)

    return ExprMetadata(
        result_expansion_kind,
        # n-ary operations align positionally, and so the last node is elementwise.
        ExprKind.ELEMENTWISE,
        has_windows=result_has_windows,
        n_orderable_ops=result_n_orderable_ops,
        preserves_length=result_preserves_length,
        is_elementwise=result_is_elementwise,
        is_scalar_like=result_is_scalar_like,
        is_literal=result_is_literal,
    )


def check_expressions_preserve_length(*args: IntoExpr, function_name: str) -> None:
    # Raise if any argument in `args` isn't length-preserving.
    # For Series input, we don't raise (yet), we let such checks happen later,
    # as this function works lazily and so can't evaluate lengths.
    from narwhals.series import Series

    if not all(
        (is_expr(x) and x._metadata.preserves_length) or isinstance(x, (str, Series))
        for x in args
    ):
        msg = f"Expressions which aggregate or change length cannot be passed to '{function_name}'."
        raise InvalidOperationError(msg)


def all_exprs_are_scalar_like(*args: IntoExpr, **kwargs: IntoExpr) -> bool:
    # Raise if any argument in `args` isn't an aggregation or literal.
    # For Series input, we don't raise (yet), we let such checks happen later,
    # as this function works lazily and so can't evaluate lengths.
    exprs = chain(args, kwargs.values())
    return all(is_expr(x) and x._metadata.is_scalar_like for x in exprs)


def apply_n_ary_operation(
    plx: CompliantNamespaceAny,
    n_ary_function: Callable[..., CompliantExprAny],
    *comparands: IntoExpr | NonNestedLiteral | _1DArray,
    str_as_lit: bool,
) -> CompliantExprAny:
    parse = plx.parse_into_expr
    compliant_exprs = (parse(into, str_as_lit=str_as_lit) for into in comparands)
    kinds = [
        ExprKind.from_into_expr(comparand, str_as_lit=str_as_lit)
        for comparand in comparands
    ]

    broadcast = any(not kind.is_scalar_like for kind in kinds)
    compliant_exprs = (
        compliant_expr.broadcast(kind)
        if broadcast and is_compliant_expr(compliant_expr) and is_scalar_like(kind)
        else compliant_expr
        for compliant_expr, kind in zip_strict(compliant_exprs, kinds)
    )
    return n_ary_function(*compliant_exprs)
