"""Simple copy propagation optimization.

Example input:

    x = f()
    y = x

The register x is redundant and we can directly assign its value to y:

    y = f()

This can optimize away registers that are assigned to once.
"""

from __future__ import annotations

from mypyc.ir.func_ir import FuncIR
from mypyc.ir.ops import (
    Assign,
    AssignMulti,
    LoadAddress,
    LoadErrorValue,
    Register,
    Value,
)
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.options import CompilerOptions
from mypyc.sametype import is_same_type
from mypyc.transform.ir_transform import IRTransform


def do_copy_propagation(fn: FuncIR, options: CompilerOptions) -> None:
    """Perform copy propagation optimization for fn."""

    # Anything with an assignment count >1 will not be optimized
    # here, as it would be require data flow analysis and we want to
    # keep this simple and fast, at least until we've made data flow
    # analysis much faster.
    counts: dict[Value, int] = {}
    replacements: dict[Value, Value] = {}
    for arg in fn.arg_regs:
        # Arguments are always assigned to initially
        counts[arg] = 1

    for block in fn.blocks:
        for op in block.ops:
            if isinstance(op, Assign):
                c = counts.get(op.dest, 0)
                counts[op.dest] = c + 1
                # Does this look like a supported assignment?
                # TODO: Something needs LoadErrorValue assignments to be preserved?
                if (
                    c == 0
                    and is_same_type(op.dest.type, op.src.type)
                    and not isinstance(op.src, LoadErrorValue)
                ):
                    replacements[op.dest] = op.src
                elif c == 1:
                    # Too many assignments -- don't replace this one
                    replacements.pop(op.dest, 0)
            elif isinstance(op, AssignMulti):
                # Copy propagation not supported for AssignMulti destinations
                counts[op.dest] = 2
                replacements.pop(op.dest, 0)
            elif isinstance(op, LoadAddress):
                # We don't support taking the address of an arbitrary Value,
                # so we'll need to preserve the operands of LoadAddress.
                if isinstance(op.src, Register):
                    counts[op.src] = 2
                    replacements.pop(op.src, 0)

    # Follow chains of propagation with more than one assignment.
    for src, dst in list(replacements.items()):
        if counts.get(dst, 0) > 1:
            # Not supported
            del replacements[src]
        else:
            while dst in replacements:
                dst = replacements[dst]
                if counts.get(dst, 0) > 1:
                    # Not supported
                    del replacements[src]
        if src in replacements:
            replacements[src] = dst

    builder = LowLevelIRBuilder(None, options)
    transform = CopyPropagationTransform(builder, replacements)
    transform.transform_blocks(fn.blocks)
    fn.blocks = builder.blocks


class CopyPropagationTransform(IRTransform):
    def __init__(self, builder: LowLevelIRBuilder, map: dict[Value, Value]) -> None:
        super().__init__(builder)
        self.op_map.update(map)
        self.removed = set(map)

    def visit_assign(self, op: Assign) -> Value | None:
        if op.dest in self.removed:
            return None
        return self.add(op)
