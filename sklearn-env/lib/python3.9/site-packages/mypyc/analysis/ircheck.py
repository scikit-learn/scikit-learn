"""Utilities for checking that internal ir is valid and consistent."""
from typing import List, Union
from mypyc.ir.pprint import format_func
from mypyc.ir.ops import (
    OpVisitor, BasicBlock, Op, ControlOp, Goto, Branch, Return, Unreachable,
    Assign, AssignMulti, LoadErrorValue, LoadLiteral, GetAttr, SetAttr, LoadStatic,
    InitStatic, TupleGet, TupleSet, IncRef, DecRef, Call, MethodCall, Cast,
    Box, Unbox, RaiseStandardError, CallC, Truncate, LoadGlobal, IntOp, ComparisonOp,
    LoadMem, SetMem, GetElementPtr, LoadAddress, KeepAlive
)
from mypyc.ir.func_ir import FuncIR


class FnError(object):
    def __init__(self, source: Union[Op, BasicBlock], desc: str) -> None:
        self.source = source
        self.desc = desc

    def __eq__(self, other: object) -> bool:
        return isinstance(other, FnError) and self.source == other.source and \
            self.desc == other.desc

    def __repr__(self) -> str:
        return f"FnError(source={self.source}, desc={self.desc})"


def check_func_ir(fn: FuncIR) -> List[FnError]:
    """Applies validations to a given function ir and returns a list of errors found."""
    errors = []

    for block in fn.blocks:
        if not block.terminated:
            errors.append(FnError(
                source=block.ops[-1] if block.ops else block,
                desc="Block not terminated",
            ))

    op_checker = OpChecker(fn)
    for block in fn.blocks:
        for op in block.ops:
            op.accept(op_checker)

    return errors + op_checker.errors


class IrCheckException(Exception):
    pass


def assert_func_ir_valid(fn: FuncIR) -> None:
    errors = check_func_ir(fn)
    if errors:
        raise IrCheckException("Internal error: Generated invalid IR: \n" + "\n".join(
            format_func(fn, [(e.source, e.desc) for e in errors])),
        )


class OpChecker(OpVisitor[None]):
    def __init__(self, parent_fn: FuncIR) -> None:
        self.parent_fn = parent_fn
        self.errors: List[FnError] = []

    def fail(self, source: Op, desc: str) -> None:
        self.errors.append(FnError(source=source, desc=desc))

    def check_control_op_targets(self, op: ControlOp) -> None:
        for target in op.targets():
            if target not in self.parent_fn.blocks:
                self.fail(source=op, desc=f"Invalid control operation target: {target.label}")

    def visit_goto(self, op: Goto) -> None:
        self.check_control_op_targets(op)

    def visit_branch(self, op: Branch) -> None:
        self.check_control_op_targets(op)

    def visit_return(self, op: Return) -> None:
        pass

    def visit_unreachable(self, op: Unreachable) -> None:
        pass

    def visit_assign(self, op: Assign) -> None:
        pass

    def visit_assign_multi(self, op: AssignMulti) -> None:
        pass

    def visit_load_error_value(self, op: LoadErrorValue) -> None:
        pass

    def visit_load_literal(self, op: LoadLiteral) -> None:
        pass

    def visit_get_attr(self, op: GetAttr) -> None:
        pass

    def visit_set_attr(self, op: SetAttr) -> None:
        pass

    def visit_load_static(self, op: LoadStatic) -> None:
        pass

    def visit_init_static(self, op: InitStatic) -> None:
        pass

    def visit_tuple_get(self, op: TupleGet) -> None:
        pass

    def visit_tuple_set(self, op: TupleSet) -> None:
        pass

    def visit_inc_ref(self, op: IncRef) -> None:
        pass

    def visit_dec_ref(self, op: DecRef) -> None:
        pass

    def visit_call(self, op: Call) -> None:
        pass

    def visit_method_call(self, op: MethodCall) -> None:
        pass

    def visit_cast(self, op: Cast) -> None:
        pass

    def visit_box(self, op: Box) -> None:
        pass

    def visit_unbox(self, op: Unbox) -> None:
        pass

    def visit_raise_standard_error(self, op: RaiseStandardError) -> None:
        pass

    def visit_call_c(self, op: CallC) -> None:
        pass

    def visit_truncate(self, op: Truncate) -> None:
        pass

    def visit_load_global(self, op: LoadGlobal) -> None:
        pass

    def visit_int_op(self, op: IntOp) -> None:
        pass

    def visit_comparison_op(self, op: ComparisonOp) -> None:
        pass

    def visit_load_mem(self, op: LoadMem) -> None:
        pass

    def visit_set_mem(self, op: SetMem) -> None:
        pass

    def visit_get_element_ptr(self, op: GetElementPtr) -> None:
        pass

    def visit_load_address(self, op: LoadAddress) -> None:
        pass

    def visit_keep_alive(self, op: KeepAlive) -> None:
        pass
