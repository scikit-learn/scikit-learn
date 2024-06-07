"""Constant folding of expressions.

For example, 3 + 5 can be constant folded into 8.
"""

from __future__ import annotations

from typing import Final, Union

from mypy.nodes import (
    ComplexExpr,
    Expression,
    FloatExpr,
    IntExpr,
    NameExpr,
    OpExpr,
    StrExpr,
    UnaryExpr,
    Var,
)

# All possible result types of constant folding
ConstantValue = Union[int, bool, float, complex, str]
CONST_TYPES: Final = (int, bool, float, complex, str)


def constant_fold_expr(expr: Expression, cur_mod_id: str) -> ConstantValue | None:
    """Return the constant value of an expression for supported operations.

    Among other things, support int arithmetic and string
    concatenation. For example, the expression 3 + 5 has the constant
    value 8.

    Also bind simple references to final constants defined in the
    current module (cur_mod_id). Binding to references is best effort
    -- we don't bind references to other modules. Mypyc trusts these
    to be correct in compiled modules, so that it can replace a
    constant expression (or a reference to one) with the statically
    computed value. We don't want to infer constant values based on
    stubs, in particular, as these might not match the implementation
    (due to version skew, for example).

    Return None if unsuccessful.
    """
    if isinstance(expr, IntExpr):
        return expr.value
    if isinstance(expr, StrExpr):
        return expr.value
    if isinstance(expr, FloatExpr):
        return expr.value
    if isinstance(expr, ComplexExpr):
        return expr.value
    elif isinstance(expr, NameExpr):
        if expr.name == "True":
            return True
        elif expr.name == "False":
            return False
        node = expr.node
        if (
            isinstance(node, Var)
            and node.is_final
            and node.fullname.rsplit(".", 1)[0] == cur_mod_id
        ):
            value = node.final_value
            if isinstance(value, (CONST_TYPES)):
                return value
    elif isinstance(expr, OpExpr):
        left = constant_fold_expr(expr.left, cur_mod_id)
        right = constant_fold_expr(expr.right, cur_mod_id)
        if left is not None and right is not None:
            return constant_fold_binary_op(expr.op, left, right)
    elif isinstance(expr, UnaryExpr):
        value = constant_fold_expr(expr.expr, cur_mod_id)
        if value is not None:
            return constant_fold_unary_op(expr.op, value)
    return None


def constant_fold_binary_op(
    op: str, left: ConstantValue, right: ConstantValue
) -> ConstantValue | None:
    if isinstance(left, int) and isinstance(right, int):
        return constant_fold_binary_int_op(op, left, right)

    # Float and mixed int/float arithmetic.
    if isinstance(left, float) and isinstance(right, float):
        return constant_fold_binary_float_op(op, left, right)
    elif isinstance(left, float) and isinstance(right, int):
        return constant_fold_binary_float_op(op, left, right)
    elif isinstance(left, int) and isinstance(right, float):
        return constant_fold_binary_float_op(op, left, right)

    # String concatenation and multiplication.
    if op == "+" and isinstance(left, str) and isinstance(right, str):
        return left + right
    elif op == "*" and isinstance(left, str) and isinstance(right, int):
        return left * right
    elif op == "*" and isinstance(left, int) and isinstance(right, str):
        return left * right

    # Complex construction.
    if op == "+" and isinstance(left, (int, float)) and isinstance(right, complex):
        return left + right
    elif op == "+" and isinstance(left, complex) and isinstance(right, (int, float)):
        return left + right
    elif op == "-" and isinstance(left, (int, float)) and isinstance(right, complex):
        return left - right
    elif op == "-" and isinstance(left, complex) and isinstance(right, (int, float)):
        return left - right

    return None


def constant_fold_binary_int_op(op: str, left: int, right: int) -> int | float | None:
    if op == "+":
        return left + right
    if op == "-":
        return left - right
    elif op == "*":
        return left * right
    elif op == "/":
        if right != 0:
            return left / right
    elif op == "//":
        if right != 0:
            return left // right
    elif op == "%":
        if right != 0:
            return left % right
    elif op == "&":
        return left & right
    elif op == "|":
        return left | right
    elif op == "^":
        return left ^ right
    elif op == "<<":
        if right >= 0:
            return left << right
    elif op == ">>":
        if right >= 0:
            return left >> right
    elif op == "**":
        if right >= 0:
            ret = left**right
            assert isinstance(ret, int)
            return ret
    return None


def constant_fold_binary_float_op(
    op: str, left: int | float, right: int | float
) -> float | None:
    assert not (isinstance(left, int) and isinstance(right, int)), (op, left, right)
    if op == "+":
        return left + right
    elif op == "-":
        return left - right
    elif op == "*":
        return left * right
    elif op == "/":
        if right != 0:
            return left / right
    elif op == "//":
        if right != 0:
            return left // right
    elif op == "%":
        if right != 0:
            return left % right
    elif op == "**":
        if (left < 0 and isinstance(right, int)) or left > 0:
            try:
                ret = left**right
            except OverflowError:
                return None
            else:
                assert isinstance(ret, float), ret
                return ret

    return None


def constant_fold_unary_op(op: str, value: ConstantValue) -> int | float | None:
    if op == "-" and isinstance(value, (int, float)):
        return -value
    elif op == "~" and isinstance(value, int):
        return ~value
    elif op == "+" and isinstance(value, (int, float)):
        return value
    return None
