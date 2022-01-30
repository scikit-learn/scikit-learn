"""Insert checks for uninitialized values."""

from typing import List

from mypyc.analysis.dataflow import (
    get_cfg,
    cleanup_cfg,
    analyze_must_defined_regs,
    AnalysisDict
)
from mypyc.ir.ops import (
    BasicBlock, Op, Branch, Value, RaiseStandardError, Unreachable, Register,
    LoadAddress, Assign, LoadErrorValue
)
from mypyc.ir.func_ir import FuncIR, all_values


def insert_uninit_checks(ir: FuncIR) -> None:
    # Remove dead blocks from the CFG, which helps avoid spurious
    # checks due to unused error handling blocks.
    cleanup_cfg(ir.blocks)

    cfg = get_cfg(ir.blocks)
    must_defined = analyze_must_defined_regs(
        ir.blocks,
        cfg,
        set(ir.arg_regs),
        all_values(ir.arg_regs, ir.blocks))

    ir.blocks = split_blocks_at_uninits(ir.blocks, must_defined.before)


def split_blocks_at_uninits(blocks: List[BasicBlock],
                            pre_must_defined: 'AnalysisDict[Value]') -> List[BasicBlock]:
    new_blocks: List[BasicBlock] = []

    init_registers = []
    init_registers_set = set()

    # First split blocks on ops that may raise.
    for block in blocks:
        ops = block.ops
        block.ops = []
        cur_block = block
        new_blocks.append(cur_block)

        for i, op in enumerate(ops):
            defined = pre_must_defined[block, i]
            for src in op.unique_sources():
                # If a register operand is not guaranteed to be
                # initialized is an operand to something other than a
                # check that it is defined, insert a check.

                # Note that for register operand in a LoadAddress op,
                # we should be able to use it without initialization
                # as we may need to use its address to update itself
                if (isinstance(src, Register) and src not in defined
                        and not (isinstance(op, Branch) and op.op == Branch.IS_ERROR)
                        and not isinstance(op, LoadAddress)):
                    new_block, error_block = BasicBlock(), BasicBlock()
                    new_block.error_handler = error_block.error_handler = cur_block.error_handler
                    new_blocks += [error_block, new_block]

                    if src not in init_registers_set:
                        init_registers.append(src)
                        init_registers_set.add(src)

                    cur_block.ops.append(Branch(src,
                                                true_label=error_block,
                                                false_label=new_block,
                                                op=Branch.IS_ERROR,
                                                line=op.line))
                    raise_std = RaiseStandardError(
                        RaiseStandardError.UNBOUND_LOCAL_ERROR,
                        'local variable "{}" referenced before assignment'.format(src.name),
                        op.line)
                    error_block.ops.append(raise_std)
                    error_block.ops.append(Unreachable())
                    cur_block = new_block
            cur_block.ops.append(op)

    if init_registers:
        new_ops: List[Op] = []
        for reg in init_registers:
            err = LoadErrorValue(reg.type, undefines=True)
            new_ops.append(err)
            new_ops.append(Assign(reg, err))
        new_blocks[0].ops[0:0] = new_ops

    return new_blocks
