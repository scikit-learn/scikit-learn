"""Convert tagged int primitive ops to lower-level ops."""

from __future__ import annotations

from typing import NamedTuple

from mypyc.ir.ops import Assign, BasicBlock, Branch, ComparisonOp, Register, Value
from mypyc.ir.rtypes import bool_rprimitive, is_short_int_rprimitive
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.lower.registry import lower_primitive_op
from mypyc.primitives.int_ops import int_equal_, int_less_than_
from mypyc.primitives.registry import CFunctionDescription


# Description for building int comparison ops
#
# Fields:
#   binary_op_variant: identify which IntOp to use when operands are short integers
#   c_func_description: the C function to call when operands are tagged integers
#   c_func_negated: whether to negate the C function call's result
#   c_func_swap_operands: whether to swap lhs and rhs when call the function
class IntComparisonOpDescription(NamedTuple):
    binary_op_variant: int
    c_func_description: CFunctionDescription
    c_func_negated: bool
    c_func_swap_operands: bool


# Provide mapping from textual op to short int's op variant and boxed int's description.
# Note that these are not complete implementations and require extra IR.
int_comparison_op_mapping: dict[str, IntComparisonOpDescription] = {
    "==": IntComparisonOpDescription(ComparisonOp.EQ, int_equal_, False, False),
    "!=": IntComparisonOpDescription(ComparisonOp.NEQ, int_equal_, True, False),
    "<": IntComparisonOpDescription(ComparisonOp.SLT, int_less_than_, False, False),
    "<=": IntComparisonOpDescription(ComparisonOp.SLE, int_less_than_, True, True),
    ">": IntComparisonOpDescription(ComparisonOp.SGT, int_less_than_, False, True),
    ">=": IntComparisonOpDescription(ComparisonOp.SGE, int_less_than_, True, False),
}


def compare_tagged(
    self: LowLevelIRBuilder, lhs: Value, rhs: Value, op: str, line: int
) -> Value:
    """Compare two tagged integers using given operator (value context)."""
    # generate fast binary logic ops on short ints
    if (
        is_short_int_rprimitive(lhs.type) or is_short_int_rprimitive(rhs.type)
    ) and op in (
        "==",
        "!=",
    ):
        quick = True
    else:
        quick = is_short_int_rprimitive(lhs.type) and is_short_int_rprimitive(rhs.type)
    if quick:
        return self.comparison_op(lhs, rhs, int_comparison_op_mapping[op][0], line)
    op_type, c_func_desc, negate_result, swap_op = int_comparison_op_mapping[op]
    result = Register(bool_rprimitive)
    short_int_block, int_block, out = BasicBlock(), BasicBlock(), BasicBlock()
    check_lhs = self.check_tagged_short_int(lhs, line, negated=True)
    if op in ("==", "!="):
        self.add(Branch(check_lhs, int_block, short_int_block, Branch.BOOL))
    else:
        # for non-equality logical ops (less/greater than, etc.), need to check both sides
        short_lhs = BasicBlock()
        self.add(Branch(check_lhs, int_block, short_lhs, Branch.BOOL))
        self.activate_block(short_lhs)
        check_rhs = self.check_tagged_short_int(rhs, line, negated=True)
        self.add(Branch(check_rhs, int_block, short_int_block, Branch.BOOL))
    self.activate_block(int_block)
    if swap_op:
        args = [rhs, lhs]
    else:
        args = [lhs, rhs]
    call = self.call_c(c_func_desc, args, line)
    if negate_result:
        # TODO: introduce UnaryIntOp?
        call_result = self.unary_op(call, "not", line)
    else:
        call_result = call
    self.add(Assign(result, call_result, line))
    self.goto(out)
    self.activate_block(short_int_block)
    eq = self.comparison_op(lhs, rhs, op_type, line)
    self.add(Assign(result, eq, line))
    self.goto_and_activate(out)
    return result


@lower_primitive_op("int_eq")
def lower_int_eq(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], "==", line)


@lower_primitive_op("int_ne")
def lower_int_ne(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], "!=", line)


@lower_primitive_op("int_lt")
def lower_int_lt(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], "<", line)


@lower_primitive_op("int_le")
def lower_int_le(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], "<=", line)


@lower_primitive_op("int_gt")
def lower_int_gt(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], ">", line)


@lower_primitive_op("int_ge")
def lower_int_ge(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    return compare_tagged(builder, args[0], args[1], ">=", line)
