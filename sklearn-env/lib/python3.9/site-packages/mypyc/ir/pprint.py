"""Utilities for pretty-printing IR in a human-readable form."""

from collections import defaultdict
from typing import Any, Dict, List, Union, Sequence, Tuple

from typing_extensions import Final

from mypyc.common import short_name
from mypyc.ir.ops import (
    Goto, Branch, Return, Unreachable, Assign, Integer, LoadErrorValue, GetAttr, SetAttr,
    LoadStatic, InitStatic, TupleGet, TupleSet, IncRef, DecRef, Call, MethodCall, Cast, Box, Unbox,
    RaiseStandardError, CallC, Truncate, LoadGlobal, IntOp, ComparisonOp, LoadMem, SetMem,
    GetElementPtr, LoadAddress, Register, Value, OpVisitor, BasicBlock, ControlOp, LoadLiteral,
    AssignMulti, KeepAlive, Op
)
from mypyc.ir.func_ir import FuncIR, all_values_full
from mypyc.ir.module_ir import ModuleIRs
from mypyc.ir.rtypes import is_bool_rprimitive, is_int_rprimitive, RType

ErrorSource = Union[BasicBlock, Op]


class IRPrettyPrintVisitor(OpVisitor[str]):
    """Internal visitor that pretty-prints ops."""

    def __init__(self, names: Dict[Value, str]) -> None:
        # This should contain a name for all values that are shown as
        # registers in the output. This is not just for Register
        # instances -- all Ops that produce values need (generated) names.
        self.names = names

    def visit_goto(self, op: Goto) -> str:
        return self.format('goto %l', op.label)

    branch_op_names: Final = {
        Branch.BOOL: ('%r', 'bool'),
        Branch.IS_ERROR: ('is_error(%r)', ''),
    }

    def visit_branch(self, op: Branch) -> str:
        fmt, typ = self.branch_op_names[op.op]
        if op.negated:
            fmt = 'not {}'.format(fmt)

        cond = self.format(fmt, op.value)
        tb = ''
        if op.traceback_entry:
            tb = ' (error at %s:%d)' % op.traceback_entry
        fmt = 'if {} goto %l{} else goto %l'.format(cond, tb)
        if typ:
            fmt += ' :: {}'.format(typ)
        return self.format(fmt, op.true, op.false)

    def visit_return(self, op: Return) -> str:
        return self.format('return %r', op.value)

    def visit_unreachable(self, op: Unreachable) -> str:
        return "unreachable"

    def visit_assign(self, op: Assign) -> str:
        return self.format('%r = %r', op.dest, op.src)

    def visit_assign_multi(self, op: AssignMulti) -> str:
        return self.format('%r = [%s]',
                           op.dest,
                           ', '.join(self.format('%r', v) for v in op.src))

    def visit_load_error_value(self, op: LoadErrorValue) -> str:
        return self.format('%r = <error> :: %s', op, op.type)

    def visit_load_literal(self, op: LoadLiteral) -> str:
        prefix = ''
        # For values that have a potential unboxed representation, make
        # it explicit that this is a Python object.
        if isinstance(op.value, int):
            prefix = 'object '
        return self.format('%r = %s%s', op, prefix, repr(op.value))

    def visit_get_attr(self, op: GetAttr) -> str:
        return self.format('%r = %r.%s', op, op.obj, op.attr)

    def visit_set_attr(self, op: SetAttr) -> str:
        return self.format('%r.%s = %r; %r = is_error', op.obj, op.attr, op.src, op)

    def visit_load_static(self, op: LoadStatic) -> str:
        ann = '  ({})'.format(repr(op.ann)) if op.ann else ''
        name = op.identifier
        if op.module_name is not None:
            name = '{}.{}'.format(op.module_name, name)
        return self.format('%r = %s :: %s%s', op, name, op.namespace, ann)

    def visit_init_static(self, op: InitStatic) -> str:
        name = op.identifier
        if op.module_name is not None:
            name = '{}.{}'.format(op.module_name, name)
        return self.format('%s = %r :: %s', name, op.value, op.namespace)

    def visit_tuple_get(self, op: TupleGet) -> str:
        return self.format('%r = %r[%d]', op, op.src, op.index)

    def visit_tuple_set(self, op: TupleSet) -> str:
        item_str = ', '.join(self.format('%r', item) for item in op.items)
        return self.format('%r = (%s)', op, item_str)

    def visit_inc_ref(self, op: IncRef) -> str:
        s = self.format('inc_ref %r', op.src)
        # TODO: Remove bool check (it's unboxed)
        if is_bool_rprimitive(op.src.type) or is_int_rprimitive(op.src.type):
            s += ' :: {}'.format(short_name(op.src.type.name))
        return s

    def visit_dec_ref(self, op: DecRef) -> str:
        s = self.format('%sdec_ref %r', 'x' if op.is_xdec else '', op.src)
        # TODO: Remove bool check (it's unboxed)
        if is_bool_rprimitive(op.src.type) or is_int_rprimitive(op.src.type):
            s += ' :: {}'.format(short_name(op.src.type.name))
        return s

    def visit_call(self, op: Call) -> str:
        args = ', '.join(self.format('%r', arg) for arg in op.args)
        # TODO: Display long name?
        short_name = op.fn.shortname
        s = '%s(%s)' % (short_name, args)
        if not op.is_void:
            s = self.format('%r = ', op) + s
        return s

    def visit_method_call(self, op: MethodCall) -> str:
        args = ', '.join(self.format('%r', arg) for arg in op.args)
        s = self.format('%r.%s(%s)', op.obj, op.method, args)
        if not op.is_void:
            s = self.format('%r = ', op) + s
        return s

    def visit_cast(self, op: Cast) -> str:
        return self.format('%r = cast(%s, %r)', op, op.type, op.src)

    def visit_box(self, op: Box) -> str:
        return self.format('%r = box(%s, %r)', op, op.src.type, op.src)

    def visit_unbox(self, op: Unbox) -> str:
        return self.format('%r = unbox(%s, %r)', op, op.type, op.src)

    def visit_raise_standard_error(self, op: RaiseStandardError) -> str:
        if op.value is not None:
            if isinstance(op.value, str):
                return self.format('%r = raise %s(%s)', op, op.class_name, repr(op.value))
            elif isinstance(op.value, Value):
                return self.format('%r = raise %s(%r)', op, op.class_name, op.value)
            else:
                assert False, 'value type must be either str or Value'
        else:
            return self.format('%r = raise %s', op, op.class_name)

    def visit_call_c(self, op: CallC) -> str:
        args_str = ', '.join(self.format('%r', arg) for arg in op.args)
        if op.is_void:
            return self.format('%s(%s)', op.function_name, args_str)
        else:
            return self.format('%r = %s(%s)', op, op.function_name, args_str)

    def visit_truncate(self, op: Truncate) -> str:
        return self.format("%r = truncate %r: %t to %t", op, op.src, op.src_type, op.type)

    def visit_load_global(self, op: LoadGlobal) -> str:
        ann = '  ({})'.format(repr(op.ann)) if op.ann else ''
        return self.format('%r = load_global %s :: static%s', op, op.identifier, ann)

    def visit_int_op(self, op: IntOp) -> str:
        return self.format('%r = %r %s %r', op, op.lhs, IntOp.op_str[op.op], op.rhs)

    def visit_comparison_op(self, op: ComparisonOp) -> str:
        if op.op in (ComparisonOp.SLT, ComparisonOp.SGT, ComparisonOp.SLE, ComparisonOp.SGE):
            sign_format = " :: signed"
        elif op.op in (ComparisonOp.ULT, ComparisonOp.UGT, ComparisonOp.ULE, ComparisonOp.UGE):
            sign_format = " :: unsigned"
        else:
            sign_format = ""
        return self.format('%r = %r %s %r%s', op, op.lhs, ComparisonOp.op_str[op.op],
                           op.rhs, sign_format)

    def visit_load_mem(self, op: LoadMem) -> str:
        return self.format("%r = load_mem %r :: %t*", op, op.src, op.type)

    def visit_set_mem(self, op: SetMem) -> str:
        return self.format("set_mem %r, %r :: %t*", op.dest, op.src, op.dest_type)

    def visit_get_element_ptr(self, op: GetElementPtr) -> str:
        return self.format("%r = get_element_ptr %r %s :: %t", op, op.src, op.field, op.src_type)

    def visit_load_address(self, op: LoadAddress) -> str:
        if isinstance(op.src, Register):
            return self.format("%r = load_address %r", op, op.src)
        else:
            return self.format("%r = load_address %s", op, op.src)

    def visit_keep_alive(self, op: KeepAlive) -> str:
        return self.format('keep_alive %s' % ', '.join(self.format('%r', v)
                                                       for v in op.src))

    # Helpers

    def format(self, fmt: str, *args: Any) -> str:
        """Helper for formatting strings.

        These format sequences are supported in fmt:

          %s: arbitrary object converted to string using str()
          %r: name of IR value/register
          %d: int
          %f: float
          %l: BasicBlock (formatted as label 'Ln')
          %t: RType
        """
        result = []
        i = 0
        arglist = list(args)
        while i < len(fmt):
            n = fmt.find('%', i)
            if n < 0:
                n = len(fmt)
            result.append(fmt[i:n])
            if n < len(fmt):
                typespec = fmt[n + 1]
                arg = arglist.pop(0)
                if typespec == 'r':
                    # Register/value
                    assert isinstance(arg, Value)
                    if isinstance(arg, Integer):
                        result.append(str(arg.value))
                    else:
                        result.append(self.names[arg])
                elif typespec == 'd':
                    # Integer
                    result.append('%d' % arg)
                elif typespec == 'f':
                    # Float
                    result.append('%f' % arg)
                elif typespec == 'l':
                    # Basic block (label)
                    assert isinstance(arg, BasicBlock)
                    result.append('L%s' % arg.label)
                elif typespec == 't':
                    # RType
                    assert isinstance(arg, RType)
                    result.append(arg.name)
                elif typespec == 's':
                    # String
                    result.append(str(arg))
                else:
                    raise ValueError('Invalid format sequence %{}'.format(typespec))
                i = n + 2
            else:
                i = n
        return ''.join(result)


def format_registers(func_ir: FuncIR,
                     names: Dict[Value, str]) -> List[str]:
    result = []
    i = 0
    regs = all_values_full(func_ir.arg_regs, func_ir.blocks)
    while i < len(regs):
        i0 = i
        group = [names[regs[i0]]]
        while i + 1 < len(regs) and regs[i + 1].type == regs[i0].type:
            i += 1
            group.append(names[regs[i]])
        i += 1
        result.append('%s :: %s' % (', '.join(group), regs[i0].type))
    return result


def format_blocks(blocks: List[BasicBlock],
                  names: Dict[Value, str],
                  source_to_error: Dict[ErrorSource, List[str]]) -> List[str]:
    """Format a list of IR basic blocks into a human-readable form."""
    # First label all of the blocks
    for i, block in enumerate(blocks):
        block.label = i

    handler_map: Dict[BasicBlock, List[BasicBlock]] = {}
    for b in blocks:
        if b.error_handler:
            handler_map.setdefault(b.error_handler, []).append(b)

    visitor = IRPrettyPrintVisitor(names)

    lines = []
    for i, block in enumerate(blocks):
        handler_msg = ''
        if block in handler_map:
            labels = sorted('L%d' % b.label for b in handler_map[block])
            handler_msg = ' (handler for {})'.format(', '.join(labels))

        lines.append('L%d:%s' % (block.label, handler_msg))
        if block in source_to_error:
            for error in source_to_error[block]:
                lines.append(f"  ERR: {error}")
        ops = block.ops
        if (isinstance(ops[-1], Goto) and i + 1 < len(blocks)
                and ops[-1].label == blocks[i + 1]
                and not source_to_error.get(ops[-1], [])):
            # Hide the last goto if it just goes to the next basic block,
            # and there are no assocatiated errors with the op.
            ops = ops[:-1]
        for op in ops:
            line = '    ' + op.accept(visitor)
            lines.append(line)
            if op in source_to_error:
                for error in source_to_error[op]:
                    lines.append(f"  ERR: {error}")

        if not isinstance(block.ops[-1], (Goto, Branch, Return, Unreachable)):
            # Each basic block needs to exit somewhere.
            lines.append('    [MISSING BLOCK EXIT OPCODE]')
    return lines


def format_func(fn: FuncIR, errors: Sequence[Tuple[ErrorSource, str]] = ()) -> List[str]:
    lines = []
    cls_prefix = fn.class_name + '.' if fn.class_name else ''
    lines.append('def {}{}({}):'.format(cls_prefix, fn.name,
                                        ', '.join(arg.name for arg in fn.args)))
    names = generate_names_for_ir(fn.arg_regs, fn.blocks)
    for line in format_registers(fn, names):
        lines.append('    ' + line)

    source_to_error = defaultdict(list)
    for source, error in errors:
        source_to_error[source].append(error)

    code = format_blocks(fn.blocks, names, source_to_error)
    lines.extend(code)
    return lines


def format_modules(modules: ModuleIRs) -> List[str]:
    ops = []
    for module in modules.values():
        for fn in module.functions:
            ops.extend(format_func(fn))
            ops.append('')
    return ops


def generate_names_for_ir(args: List[Register], blocks: List[BasicBlock]) -> Dict[Value, str]:
    """Generate unique names for IR values.

    Give names such as 'r5' to temp values in IR which are useful when
    pretty-printing or generating C. Ensure generated names are unique.
    """
    names: Dict[Value, str] = {}
    used_names = set()

    temp_index = 0

    for arg in args:
        names[arg] = arg.name
        used_names.add(arg.name)

    for block in blocks:
        for op in block.ops:
            values = []

            for source in op.sources():
                if source not in names:
                    values.append(source)

            if isinstance(op, (Assign, AssignMulti)):
                values.append(op.dest)
            elif isinstance(op, ControlOp) or op.is_void:
                continue
            elif op not in names:
                values.append(op)

            for value in values:
                if value in names:
                    continue
                if isinstance(value, Register) and value.name:
                    name = value.name
                elif isinstance(value, Integer):
                    continue
                else:
                    name = 'r%d' % temp_index
                    temp_index += 1

                # Append _2, _3, ... if needed to make the name unique.
                if name in used_names:
                    n = 2
                    while True:
                        candidate = '%s_%d' % (name, n)
                        if candidate not in used_names:
                            name = candidate
                            break
                        n += 1

                names[value] = name
                used_names.add(name)

    return names
