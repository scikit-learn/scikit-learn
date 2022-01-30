import unittest
from typing import List

from mypyc.analysis.ircheck import check_func_ir, FnError
from mypyc.ir.rtypes import none_rprimitive
from mypyc.ir.ops import BasicBlock, Op, Return, Integer, Goto
from mypyc.ir.func_ir import FuncIR, FuncDecl, FuncSignature
from mypyc.ir.pprint import format_func


def assert_has_error(fn: FuncIR, error: FnError) -> None:
    errors = check_func_ir(fn)
    assert errors == [error]


def assert_no_errors(fn: FuncIR) -> None:
    assert not check_func_ir(fn)


NONE_VALUE = Integer(0, rtype=none_rprimitive)


class TestIrcheck(unittest.TestCase):
    def setUp(self) -> None:
        self.label = 0

    def basic_block(self, ops: List[Op]) -> BasicBlock:
        self.label += 1
        block = BasicBlock(self.label)
        block.ops = ops
        return block

    def func_decl(self, name: str) -> FuncDecl:
        return FuncDecl(name=name, class_name=None, module_name="module", sig=FuncSignature(
            args=[], ret_type=none_rprimitive,
        ))

    def test_valid_fn(self) -> None:
        assert_no_errors(FuncIR(
            decl=self.func_decl(name="func_1"),
            arg_regs=[],
            blocks=[self.basic_block(ops=[
                Return(value=NONE_VALUE),
            ])],
        ))

    def test_block_not_terminated_empty_block(self) -> None:
        block = self.basic_block([])
        fn = FuncIR(
            decl=self.func_decl(name="func_1"),
            arg_regs=[],
            blocks=[block],
        )
        assert_has_error(fn, FnError(source=block, desc="Block not terminated"))

    def test_valid_goto(self) -> None:
        block_1 = self.basic_block([Return(value=NONE_VALUE)])
        block_2 = self.basic_block([Goto(label=block_1)])
        fn = FuncIR(
            decl=self.func_decl(name="func_1"),
            arg_regs=[],
            blocks=[block_1, block_2],
        )
        assert_no_errors(fn)

    def test_invalid_goto(self) -> None:
        block_1 = self.basic_block([Return(value=NONE_VALUE)])
        goto = Goto(label=block_1)
        block_2 = self.basic_block([goto])
        fn = FuncIR(
            decl=self.func_decl(name="func_1"),
            arg_regs=[],
            # block_1 omitted
            blocks=[block_2],
        )
        assert_has_error(fn, FnError(source=goto, desc="Invalid control operation target: 1"))

    def test_pprint(self) -> None:
        block_1 = self.basic_block([Return(value=NONE_VALUE)])
        goto = Goto(label=block_1)
        block_2 = self.basic_block([goto])
        fn = FuncIR(
            decl=self.func_decl(name="func_1"),
            arg_regs=[],
            # block_1 omitted
            blocks=[block_2],
        )
        errors = [(goto, "Invalid control operation target: 1")]
        formatted = format_func(fn, errors)
        assert formatted == [
            "def func_1():",
            "L0:",
            "    goto L1",
            "  ERR: Invalid control operation target: 1",
        ]
