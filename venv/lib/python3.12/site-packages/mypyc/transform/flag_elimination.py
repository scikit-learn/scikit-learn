"""Bool register elimination optimization.

Example input:

  L1:
    r0 = f()
    b = r0
    goto L3
  L2:
    r1 = g()
    b = r1
    goto L3
  L3:
    if b goto L4 else goto L5

The register b is redundant and we replace the assignments with two copies of
the branch in L3:

  L1:
    r0 = f()
    if r0 goto L4 else goto L5
  L2:
    r1 = g()
    if r1 goto L4 else goto L5

This helps generate simpler IR for tagged integers comparisons, for example.
"""

from __future__ import annotations

from mypyc.ir.func_ir import FuncIR
from mypyc.ir.ops import Assign, BasicBlock, Branch, Goto, Register, Unreachable
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.options import CompilerOptions
from mypyc.transform.ir_transform import IRTransform


def do_flag_elimination(fn: FuncIR, options: CompilerOptions) -> None:
    # Find registers that are used exactly once as source, and in a branch.
    counts: dict[Register, int] = {}
    branches: dict[Register, Branch] = {}
    labels: dict[Register, BasicBlock] = {}
    for block in fn.blocks:
        for i, op in enumerate(block.ops):
            for src in op.sources():
                if isinstance(src, Register):
                    counts[src] = counts.get(src, 0) + 1
            if i == 0 and isinstance(op, Branch) and isinstance(op.value, Register):
                branches[op.value] = op
                labels[op.value] = block

    # Based on these we can find the candidate registers.
    candidates: set[Register] = {
        r for r in branches if counts.get(r, 0) == 1 and r not in fn.arg_regs
    }

    # Remove candidates with invalid assignments.
    for block in fn.blocks:
        for i, op in enumerate(block.ops):
            if isinstance(op, Assign) and op.dest in candidates:
                next_op = block.ops[i + 1]
                if not (isinstance(next_op, Goto) and next_op.label is labels[op.dest]):
                    # Not right
                    candidates.remove(op.dest)

    builder = LowLevelIRBuilder(None, options)
    transform = FlagEliminationTransform(
        builder, {x: y for x, y in branches.items() if x in candidates}
    )
    transform.transform_blocks(fn.blocks)
    fn.blocks = builder.blocks


class FlagEliminationTransform(IRTransform):
    def __init__(
        self, builder: LowLevelIRBuilder, branch_map: dict[Register, Branch]
    ) -> None:
        super().__init__(builder)
        self.branch_map = branch_map
        self.branches = set(branch_map.values())

    def visit_assign(self, op: Assign) -> None:
        old_branch = self.branch_map.get(op.dest)
        if old_branch:
            # Replace assignment with a copy of the old branch, which is in a
            # separate basic block. The old branch will be deletecd in visit_branch.
            new_branch = Branch(
                op.src,
                old_branch.true,
                old_branch.false,
                old_branch.op,
                old_branch.line,
                rare=old_branch.rare,
            )
            new_branch.negated = old_branch.negated
            new_branch.traceback_entry = old_branch.traceback_entry
            self.add(new_branch)
        else:
            self.add(op)

    def visit_goto(self, op: Goto) -> None:
        # This is a no-op if basic block already terminated
        self.builder.goto(op.label)

    def visit_branch(self, op: Branch) -> None:
        if op in self.branches:
            # This branch is optimized away
            self.add(Unreachable())
        else:
            self.add(op)
