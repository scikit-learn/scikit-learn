"""Transform IR to lower-level ops.

Higher-level ops are used in earlier compiler passes, as they make
various analyses, optimizations and transforms easier to implement.
Later passes use lower-level ops, as they are easier to generate code
from, and they help with lower-level optimizations.

Lowering of various primitive ops is implemented in the mypyc.lower
package.
"""

from mypyc.ir.func_ir import FuncIR
from mypyc.ir.ops import PrimitiveOp, Value
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.lower.registry import lowering_registry
from mypyc.options import CompilerOptions
from mypyc.transform.ir_transform import IRTransform


def lower_ir(ir: FuncIR, options: CompilerOptions) -> None:
    builder = LowLevelIRBuilder(None, options)
    visitor = LoweringVisitor(builder)
    visitor.transform_blocks(ir.blocks)
    ir.blocks = builder.blocks


class LoweringVisitor(IRTransform):
    def visit_primitive_op(self, op: PrimitiveOp) -> Value:
        # The lowering implementation functions of various primitive ops are stored
        # in a registry, which is populated using function decorators. The name
        # of op (such as "int_eq") is used as the key.
        lower_fn = lowering_registry[op.desc.name]
        return lower_fn(self.builder, op.args, op.line)
