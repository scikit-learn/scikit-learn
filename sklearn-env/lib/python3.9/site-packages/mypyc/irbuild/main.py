"""Transform a mypy AST to the IR form (Intermediate Representation).

For example, consider a function like this:

   def f(x: int) -> int:
       return x * 2 + 1

It would be translated to something that conceptually looks like this:

   r0 = 2
   r1 = 1
   r2 = x * r0 :: int
   r3 = r2 + r1 :: int
   return r3

This module deals with the module-level IR transformation logic and
putting it all together. The actual IR is implemented in mypyc.ir.

For the core of the IR transform implementation, look at build_ir()
below, mypyc.irbuild.builder, and mypyc.irbuild.visitor.
"""

from mypy.backports import OrderedDict
from typing import List, Dict, Callable, Any, TypeVar, cast

from mypy.nodes import MypyFile, Expression, ClassDef
from mypy.types import Type
from mypy.state import strict_optional_set
from mypy.build import Graph

from mypyc.common import TOP_LEVEL_NAME
from mypyc.errors import Errors
from mypyc.options import CompilerOptions
from mypyc.ir.rtypes import none_rprimitive
from mypyc.ir.module_ir import ModuleIR, ModuleIRs
from mypyc.ir.func_ir import FuncIR, FuncDecl, FuncSignature
from mypyc.irbuild.prebuildvisitor import PreBuildVisitor
from mypyc.irbuild.vtable import compute_vtable
from mypyc.irbuild.prepare import build_type_map, find_singledispatch_register_impls
from mypyc.irbuild.builder import IRBuilder
from mypyc.irbuild.visitor import IRBuilderVisitor
from mypyc.irbuild.mapper import Mapper


# The stubs for callable contextmanagers are busted so cast it to the
# right type...
F = TypeVar('F', bound=Callable[..., Any])
strict_optional_dec = cast(Callable[[F], F], strict_optional_set(True))


@strict_optional_dec  # Turn on strict optional for any type manipulations we do
def build_ir(modules: List[MypyFile],
             graph: Graph,
             types: Dict[Expression, Type],
             mapper: 'Mapper',
             options: CompilerOptions,
             errors: Errors) -> ModuleIRs:
    """Build IR for a set of modules that have been type-checked by mypy."""

    build_type_map(mapper, modules, graph, types, options, errors)
    singledispatch_info = find_singledispatch_register_impls(modules, errors)

    result: ModuleIRs = OrderedDict()

    # Generate IR for all modules.
    class_irs = []

    for module in modules:
        # First pass to determine free symbols.
        pbv = PreBuildVisitor(errors, module, singledispatch_info.decorators_to_remove)
        module.accept(pbv)

        # Construct and configure builder objects (cyclic runtime dependency).
        visitor = IRBuilderVisitor()
        builder = IRBuilder(
            module.fullname, types, graph, errors, mapper, pbv, visitor, options,
            singledispatch_info.singledispatch_impls,
        )
        visitor.builder = builder

        # Second pass does the bulk of the work.
        transform_mypy_file(builder, module)
        module_ir = ModuleIR(
            module.fullname,
            list(builder.imports),
            builder.functions,
            builder.classes,
            builder.final_names
        )
        result[module.fullname] = module_ir
        class_irs.extend(builder.classes)

    # Compute vtables.
    for cir in class_irs:
        if cir.is_ext_class:
            compute_vtable(cir)

    return result


def transform_mypy_file(builder: IRBuilder, mypyfile: MypyFile) -> None:
    """Generate IR for a single module."""

    if mypyfile.fullname in ('typing', 'abc'):
        # These module are special; their contents are currently all
        # built-in primitives.
        return

    builder.set_module(mypyfile.fullname, mypyfile.path)

    classes = [node for node in mypyfile.defs if isinstance(node, ClassDef)]

    # Collect all classes.
    for cls in classes:
        ir = builder.mapper.type_to_ir[cls.info]
        builder.classes.append(ir)

    builder.enter('<top level>')

    # Make sure we have a builtins import
    builder.gen_import('builtins', -1)

    # Generate ops.
    for node in mypyfile.defs:
        builder.accept(node)

    builder.maybe_add_implicit_return()

    # Generate special function representing module top level.
    args, _, blocks, ret_type, _ = builder.leave()
    sig = FuncSignature([], none_rprimitive)
    func_ir = FuncIR(FuncDecl(TOP_LEVEL_NAME, None, builder.module_name, sig), args, blocks,
                     traceback_name="<module>")
    builder.functions.append(func_ir)
