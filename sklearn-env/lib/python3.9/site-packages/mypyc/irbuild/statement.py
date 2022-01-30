"""Transform mypy statement ASTs to mypyc IR (Intermediate Representation).

The top-level AST transformation logic is implemented in mypyc.irbuild.visitor
and mypyc.irbuild.builder.

A few statements are transformed in mypyc.irbuild.function (yield, for example).
"""

from typing import Optional, List, Tuple, Sequence, Callable
import importlib.util

from mypy.nodes import (
    Block, ExpressionStmt, ReturnStmt, AssignmentStmt, OperatorAssignmentStmt, IfStmt, WhileStmt,
    ForStmt, BreakStmt, ContinueStmt, RaiseStmt, TryStmt, WithStmt, AssertStmt, DelStmt,
    Expression, StrExpr, TempNode, Lvalue, Import, ImportFrom, ImportAll, TupleExpr, ListExpr,
    StarExpr
)

from mypyc.ir.ops import (
    Assign, Unreachable, RaiseStandardError, LoadErrorValue, BasicBlock, TupleGet, Value, Register,
    Branch, NO_TRACEBACK_LINE_NO
)
from mypyc.ir.rtypes import RInstance, exc_rtuple
from mypyc.primitives.generic_ops import py_delattr_op
from mypyc.primitives.misc_ops import type_op, import_from_op
from mypyc.primitives.exc_ops import (
    raise_exception_op, reraise_exception_op, error_catch_op, exc_matches_op, restore_exc_info_op,
    get_exc_value_op, keep_propagating_op, get_exc_info_op
)
from mypyc.irbuild.targets import (
    AssignmentTarget, AssignmentTargetRegister, AssignmentTargetIndex, AssignmentTargetAttr,
    AssignmentTargetTuple
)
from mypyc.irbuild.nonlocalcontrol import (
    ExceptNonlocalControl, FinallyNonlocalControl, TryFinallyNonlocalControl
)
from mypyc.irbuild.for_helpers import for_loop_helper
from mypyc.irbuild.builder import IRBuilder

GenFunc = Callable[[], None]


def transform_block(builder: IRBuilder, block: Block) -> None:
    if not block.is_unreachable:
        for stmt in block.body:
            builder.accept(stmt)
    # Raise a RuntimeError if we hit a non-empty unreachable block.
    # Don't complain about empty unreachable blocks, since mypy inserts
    # those after `if MYPY`.
    elif block.body:
        builder.add(RaiseStandardError(RaiseStandardError.RUNTIME_ERROR,
                                       'Reached allegedly unreachable code!',
                                       block.line))
        builder.add(Unreachable())


def transform_expression_stmt(builder: IRBuilder, stmt: ExpressionStmt) -> None:
    if isinstance(stmt.expr, StrExpr):
        # Docstring. Ignore
        return
    # ExpressionStmts do not need to be coerced like other Expressions.
    stmt.expr.accept(builder.visitor)


def transform_return_stmt(builder: IRBuilder, stmt: ReturnStmt) -> None:
    if stmt.expr:
        retval = builder.accept(stmt.expr)
    else:
        retval = builder.builder.none()
    retval = builder.coerce(retval, builder.ret_types[-1], stmt.line)
    builder.nonlocal_control[-1].gen_return(builder, retval, stmt.line)


def transform_assignment_stmt(builder: IRBuilder, stmt: AssignmentStmt) -> None:
    lvalues = stmt.lvalues
    assert len(lvalues) >= 1
    builder.disallow_class_assignments(lvalues, stmt.line)
    first_lvalue = lvalues[0]
    if stmt.type and isinstance(stmt.rvalue, TempNode):
        # This is actually a variable annotation without initializer. Don't generate
        # an assignment but we need to call get_assignment_target since it adds a
        # name binding as a side effect.
        builder.get_assignment_target(first_lvalue, stmt.line)
        return

    # Special case multiple assignments like 'x, y = e1, e2'.
    if (isinstance(first_lvalue, (TupleExpr, ListExpr))
            and isinstance(stmt.rvalue, (TupleExpr, ListExpr))
            and len(first_lvalue.items) == len(stmt.rvalue.items)
            and all(is_simple_lvalue(item) for item in first_lvalue.items)
            and len(lvalues) == 1):
        temps = []
        for right in stmt.rvalue.items:
            rvalue_reg = builder.accept(right)
            temp = Register(rvalue_reg.type)
            builder.assign(temp, rvalue_reg, stmt.line)
            temps.append(temp)
        for (left, temp) in zip(first_lvalue.items, temps):
            assignment_target = builder.get_assignment_target(left)
            builder.assign(assignment_target, temp, stmt.line)
        return

    line = stmt.rvalue.line
    rvalue_reg = builder.accept(stmt.rvalue)
    if builder.non_function_scope() and stmt.is_final_def:
        builder.init_final_static(first_lvalue, rvalue_reg)
    for lvalue in lvalues:
        target = builder.get_assignment_target(lvalue)
        builder.assign(target, rvalue_reg, line)


def is_simple_lvalue(expr: Expression) -> bool:
    return not isinstance(expr, (StarExpr, ListExpr, TupleExpr))


def transform_operator_assignment_stmt(builder: IRBuilder, stmt: OperatorAssignmentStmt) -> None:
    """Operator assignment statement such as x += 1"""
    builder.disallow_class_assignments([stmt.lvalue], stmt.line)
    target = builder.get_assignment_target(stmt.lvalue)
    target_value = builder.read(target, stmt.line)
    rreg = builder.accept(stmt.rvalue)
    # the Python parser strips the '=' from operator assignment statements, so re-add it
    op = stmt.op + '='
    res = builder.binary_op(target_value, rreg, op, stmt.line)
    # usually operator assignments are done in-place
    # but when target doesn't support that we need to manually assign
    builder.assign(target, res, res.line)


def transform_import(builder: IRBuilder, node: Import) -> None:
    if node.is_mypy_only:
        return
    globals = builder.load_globals_dict()
    for node_id, as_name in node.ids:
        builder.gen_import(node_id, node.line)

        # Update the globals dict with the appropriate module:
        # * For 'import foo.bar as baz' we add 'foo.bar' with the name 'baz'
        # * For 'import foo.bar' we add 'foo' with the name 'foo'
        # Typically we then ignore these entries and access things directly
        # via the module static, but we will use the globals version for modules
        # that mypy couldn't find, since it doesn't analyze module references
        # from those properly.

        # TODO: Don't add local imports to the global namespace

        # Miscompiling imports inside of functions, like below in import from.
        if as_name:
            name = as_name
            base = node_id
        else:
            base = name = node_id.split('.')[0]

        obj = builder.get_module(base, node.line)

        builder.gen_method_call(
            globals, '__setitem__', [builder.load_str(name), obj],
            result_type=None, line=node.line)


def transform_import_from(builder: IRBuilder, node: ImportFrom) -> None:
    if node.is_mypy_only:
        return

    module_state = builder.graph[builder.module_name]
    if module_state.ancestors is not None and module_state.ancestors:
        module_package = module_state.ancestors[0]
    elif builder.module_path.endswith("__init__.py"):
        module_package = builder.module_name
    else:
        module_package = ''

    id = importlib.util.resolve_name('.' * node.relative + node.id, module_package)

    globals = builder.load_globals_dict()
    imported_names = [name for name, _ in node.names]
    module = builder.gen_import_from(id, globals, imported_names, node.line)

    # Copy everything into our module's dict.
    # Note that we miscompile import from inside of functions here,
    # since that case *shouldn't* load it into the globals dict.
    # This probably doesn't matter much and the code runs basically right.
    for name, maybe_as_name in node.names:
        as_name = maybe_as_name or name
        obj = builder.call_c(import_from_op,
                             [module, builder.load_str(id),
                              builder.load_str(name), builder.load_str(as_name)],
                             node.line)
        builder.gen_method_call(
            globals, '__setitem__', [builder.load_str(as_name), obj],
            result_type=None, line=node.line)


def transform_import_all(builder: IRBuilder, node: ImportAll) -> None:
    if node.is_mypy_only:
        return
    builder.gen_import(node.id, node.line)


def transform_if_stmt(builder: IRBuilder, stmt: IfStmt) -> None:
    if_body, next = BasicBlock(), BasicBlock()
    else_body = BasicBlock() if stmt.else_body else next

    # If statements are normalized
    assert len(stmt.expr) == 1

    builder.process_conditional(stmt.expr[0], if_body, else_body)
    builder.activate_block(if_body)
    builder.accept(stmt.body[0])
    builder.goto(next)
    if stmt.else_body:
        builder.activate_block(else_body)
        builder.accept(stmt.else_body)
        builder.goto(next)
    builder.activate_block(next)


def transform_while_stmt(builder: IRBuilder, s: WhileStmt) -> None:
    body, next, top, else_block = BasicBlock(), BasicBlock(), BasicBlock(), BasicBlock()
    normal_loop_exit = else_block if s.else_body is not None else next

    builder.push_loop_stack(top, next)

    # Split block so that we get a handle to the top of the loop.
    builder.goto_and_activate(top)
    builder.process_conditional(s.expr, body, normal_loop_exit)

    builder.activate_block(body)
    builder.accept(s.body)
    # Add branch to the top at the end of the body.
    builder.goto(top)

    builder.pop_loop_stack()

    if s.else_body is not None:
        builder.activate_block(else_block)
        builder.accept(s.else_body)
        builder.goto(next)

    builder.activate_block(next)


def transform_for_stmt(builder: IRBuilder, s: ForStmt) -> None:
    if s.is_async:
        builder.error('async for is unimplemented', s.line)

    def body() -> None:
        builder.accept(s.body)

    def else_block() -> None:
        assert s.else_body is not None
        builder.accept(s.else_body)

    for_loop_helper(builder, s.index, s.expr, body,
                    else_block if s.else_body else None, s.line)


def transform_break_stmt(builder: IRBuilder, node: BreakStmt) -> None:
    builder.nonlocal_control[-1].gen_break(builder, node.line)


def transform_continue_stmt(builder: IRBuilder, node: ContinueStmt) -> None:
    builder.nonlocal_control[-1].gen_continue(builder, node.line)


def transform_raise_stmt(builder: IRBuilder, s: RaiseStmt) -> None:
    if s.expr is None:
        builder.call_c(reraise_exception_op, [], NO_TRACEBACK_LINE_NO)
        builder.add(Unreachable())
        return

    exc = builder.accept(s.expr)
    builder.call_c(raise_exception_op, [exc], s.line)
    builder.add(Unreachable())


def transform_try_except(builder: IRBuilder,
                         body: GenFunc,
                         handlers: Sequence[
                             Tuple[Optional[Expression], Optional[Expression], GenFunc]],
                         else_body: Optional[GenFunc],
                         line: int) -> None:
    """Generalized try/except/else handling that takes functions to gen the bodies.

    The point of this is to also be able to support with."""
    assert handlers, "try needs except"

    except_entry, exit_block, cleanup_block = BasicBlock(), BasicBlock(), BasicBlock()
    double_except_block = BasicBlock()
    # If there is an else block, jump there after the try, otherwise just leave
    else_block = BasicBlock() if else_body else exit_block

    # Compile the try block with an error handler
    builder.builder.push_error_handler(except_entry)
    builder.goto_and_activate(BasicBlock())
    body()
    builder.goto(else_block)
    builder.builder.pop_error_handler()

    # The error handler catches the error and then checks it
    # against the except clauses. We compile the error handler
    # itself with an error handler so that it can properly restore
    # the *old* exc_info if an exception occurs.
    # The exception chaining will be done automatically when the
    # exception is raised, based on the exception in exc_info.
    builder.builder.push_error_handler(double_except_block)
    builder.activate_block(except_entry)
    old_exc = builder.maybe_spill(builder.call_c(error_catch_op, [], line))
    # Compile the except blocks with the nonlocal control flow overridden to clear exc_info
    builder.nonlocal_control.append(
        ExceptNonlocalControl(builder.nonlocal_control[-1], old_exc))

    # Process the bodies
    for type, var, handler_body in handlers:
        next_block = None
        if type:
            next_block, body_block = BasicBlock(), BasicBlock()
            matches = builder.call_c(
                exc_matches_op, [builder.accept(type)], type.line
            )
            builder.add(Branch(matches, body_block, next_block, Branch.BOOL))
            builder.activate_block(body_block)
        if var:
            target = builder.get_assignment_target(var)
            builder.assign(
                target,
                builder.call_c(get_exc_value_op, [], var.line),
                var.line
            )
        handler_body()
        builder.goto(cleanup_block)
        if next_block:
            builder.activate_block(next_block)

    # Reraise the exception if needed
    if next_block:
        builder.call_c(reraise_exception_op, [], NO_TRACEBACK_LINE_NO)
        builder.add(Unreachable())

    builder.nonlocal_control.pop()
    builder.builder.pop_error_handler()

    # Cleanup for if we leave except through normal control flow:
    # restore the saved exc_info information and continue propagating
    # the exception if it exists.
    builder.activate_block(cleanup_block)
    builder.call_c(restore_exc_info_op, [builder.read(old_exc)], line)
    builder.goto(exit_block)

    # Cleanup for if we leave except through a raised exception:
    # restore the saved exc_info information and continue propagating
    # the exception.
    builder.activate_block(double_except_block)
    builder.call_c(restore_exc_info_op, [builder.read(old_exc)], line)
    builder.call_c(keep_propagating_op, [], NO_TRACEBACK_LINE_NO)
    builder.add(Unreachable())

    # If present, compile the else body in the obvious way
    if else_body:
        builder.activate_block(else_block)
        else_body()
        builder.goto(exit_block)

    builder.activate_block(exit_block)


def transform_try_except_stmt(builder: IRBuilder, t: TryStmt) -> None:
    def body() -> None:
        builder.accept(t.body)

    # Work around scoping woes
    def make_handler(body: Block) -> GenFunc:
        return lambda: builder.accept(body)

    handlers = [(type, var, make_handler(body))
                for type, var, body in zip(t.types, t.vars, t.handlers)]
    else_body = (lambda: builder.accept(t.else_body)) if t.else_body else None
    transform_try_except(builder, body, handlers, else_body, t.line)


def try_finally_try(builder: IRBuilder,
                    err_handler: BasicBlock,
                    return_entry: BasicBlock,
                    main_entry: BasicBlock,
                    try_body: GenFunc) -> Optional[Register]:
    # Compile the try block with an error handler
    control = TryFinallyNonlocalControl(return_entry)
    builder.builder.push_error_handler(err_handler)

    builder.nonlocal_control.append(control)
    builder.goto_and_activate(BasicBlock())
    try_body()
    builder.goto(main_entry)
    builder.nonlocal_control.pop()
    builder.builder.pop_error_handler()

    return control.ret_reg


def try_finally_entry_blocks(builder: IRBuilder,
                             err_handler: BasicBlock,
                             return_entry: BasicBlock,
                             main_entry: BasicBlock,
                             finally_block: BasicBlock,
                             ret_reg: Optional[Register]) -> Value:
    old_exc = Register(exc_rtuple)

    # Entry block for non-exceptional flow
    builder.activate_block(main_entry)
    if ret_reg:
        builder.add(
            Assign(
                ret_reg,
                builder.add(LoadErrorValue(builder.ret_types[-1]))
            )
        )
    builder.goto(return_entry)

    builder.activate_block(return_entry)
    builder.add(Assign(old_exc, builder.add(LoadErrorValue(exc_rtuple))))
    builder.goto(finally_block)

    # Entry block for errors
    builder.activate_block(err_handler)
    if ret_reg:
        builder.add(
            Assign(
                ret_reg,
                builder.add(LoadErrorValue(builder.ret_types[-1]))
            )
        )
    builder.add(Assign(old_exc, builder.call_c(error_catch_op, [], -1)))
    builder.goto(finally_block)

    return old_exc


def try_finally_body(
        builder: IRBuilder,
        finally_block: BasicBlock,
        finally_body: GenFunc,
        ret_reg: Optional[Value],
        old_exc: Value) -> Tuple[BasicBlock, FinallyNonlocalControl]:
    cleanup_block = BasicBlock()
    # Compile the finally block with the nonlocal control flow overridden to restore exc_info
    builder.builder.push_error_handler(cleanup_block)
    finally_control = FinallyNonlocalControl(
        builder.nonlocal_control[-1], ret_reg, old_exc)
    builder.nonlocal_control.append(finally_control)
    builder.activate_block(finally_block)
    finally_body()
    builder.nonlocal_control.pop()

    return cleanup_block, finally_control


def try_finally_resolve_control(builder: IRBuilder,
                                cleanup_block: BasicBlock,
                                finally_control: FinallyNonlocalControl,
                                old_exc: Value,
                                ret_reg: Optional[Value]) -> BasicBlock:
    """Resolve the control flow out of a finally block.

    This means returning if there was a return, propagating
    exceptions, break/continue (soon), or just continuing on.
    """
    reraise, rest = BasicBlock(), BasicBlock()
    builder.add(Branch(old_exc, rest, reraise, Branch.IS_ERROR))

    # Reraise the exception if there was one
    builder.activate_block(reraise)
    builder.call_c(reraise_exception_op, [], NO_TRACEBACK_LINE_NO)
    builder.add(Unreachable())
    builder.builder.pop_error_handler()

    # If there was a return, keep returning
    if ret_reg:
        builder.activate_block(rest)
        return_block, rest = BasicBlock(), BasicBlock()
        builder.add(Branch(ret_reg, rest, return_block, Branch.IS_ERROR))

        builder.activate_block(return_block)
        builder.nonlocal_control[-1].gen_return(builder, ret_reg, -1)

    # TODO: handle break/continue
    builder.activate_block(rest)
    out_block = BasicBlock()
    builder.goto(out_block)

    # If there was an exception, restore again
    builder.activate_block(cleanup_block)
    finally_control.gen_cleanup(builder, -1)
    builder.call_c(keep_propagating_op, [], NO_TRACEBACK_LINE_NO)
    builder.add(Unreachable())

    return out_block


def transform_try_finally_stmt(builder: IRBuilder,
                               try_body: GenFunc,
                               finally_body: GenFunc) -> None:
    """Generalized try/finally handling that takes functions to gen the bodies.

    The point of this is to also be able to support with."""
    # Finally is a big pain, because there are so many ways that
    # exits can occur. We emit 10+ basic blocks for every finally!

    err_handler, main_entry, return_entry, finally_block = (
        BasicBlock(), BasicBlock(), BasicBlock(), BasicBlock())

    # Compile the body of the try
    ret_reg = try_finally_try(
        builder, err_handler, return_entry, main_entry, try_body)

    # Set up the entry blocks for the finally statement
    old_exc = try_finally_entry_blocks(
        builder, err_handler, return_entry, main_entry, finally_block, ret_reg)

    # Compile the body of the finally
    cleanup_block, finally_control = try_finally_body(
        builder, finally_block, finally_body, ret_reg, old_exc)

    # Resolve the control flow out of the finally block
    out_block = try_finally_resolve_control(
        builder, cleanup_block, finally_control, old_exc, ret_reg)

    builder.activate_block(out_block)


def transform_try_stmt(builder: IRBuilder, t: TryStmt) -> None:
    # Our compilation strategy for try/except/else/finally is to
    # treat try/except/else and try/finally as separate language
    # constructs that we compile separately. When we have a
    # try/except/else/finally, we treat the try/except/else as the
    # body of a try/finally block.
    if t.finally_body:
        def transform_try_body() -> None:
            if t.handlers:
                transform_try_except_stmt(builder, t)
            else:
                builder.accept(t.body)
        body = t.finally_body

        transform_try_finally_stmt(builder, transform_try_body, lambda: builder.accept(body))
    else:
        transform_try_except_stmt(builder, t)


def get_sys_exc_info(builder: IRBuilder) -> List[Value]:
    exc_info = builder.call_c(get_exc_info_op, [], -1)
    return [builder.add(TupleGet(exc_info, i, -1)) for i in range(3)]


def transform_with(builder: IRBuilder,
                   expr: Expression,
                   target: Optional[Lvalue],
                   body: GenFunc,
                   line: int) -> None:
    # This is basically a straight transcription of the Python code in PEP 343.
    # I don't actually understand why a bunch of it is the way it is.
    # We could probably optimize the case where the manager is compiled by us,
    # but that is not our common case at all, so.
    mgr_v = builder.accept(expr)
    typ = builder.call_c(type_op, [mgr_v], line)
    exit_ = builder.maybe_spill(builder.py_get_attr(typ, '__exit__', line))
    value = builder.py_call(
        builder.py_get_attr(typ, '__enter__', line), [mgr_v], line
    )
    mgr = builder.maybe_spill(mgr_v)
    exc = builder.maybe_spill_assignable(builder.true())

    def try_body() -> None:
        if target:
            builder.assign(builder.get_assignment_target(target), value, line)
        body()

    def except_body() -> None:
        builder.assign(exc, builder.false(), line)
        out_block, reraise_block = BasicBlock(), BasicBlock()
        builder.add_bool_branch(
            builder.py_call(builder.read(exit_),
                            [builder.read(mgr)] + get_sys_exc_info(builder), line),
            out_block,
            reraise_block
        )
        builder.activate_block(reraise_block)
        builder.call_c(reraise_exception_op, [], NO_TRACEBACK_LINE_NO)
        builder.add(Unreachable())
        builder.activate_block(out_block)

    def finally_body() -> None:
        out_block, exit_block = BasicBlock(), BasicBlock()
        builder.add(
            Branch(builder.read(exc), exit_block, out_block, Branch.BOOL)
        )
        builder.activate_block(exit_block)
        none = builder.none_object()
        builder.py_call(
            builder.read(exit_), [builder.read(mgr), none, none, none], line
        )
        builder.goto_and_activate(out_block)

    transform_try_finally_stmt(
        builder,
        lambda: transform_try_except(builder,
                                     try_body,
                                     [(None, None, except_body)],
                                     None,
                                     line),
        finally_body
    )


def transform_with_stmt(builder: IRBuilder, o: WithStmt) -> None:
    if o.is_async:
        builder.error('async with is unimplemented', o.line)

    # Generate separate logic for each expr in it, left to right
    def generate(i: int) -> None:
        if i >= len(o.expr):
            builder.accept(o.body)
        else:
            transform_with(builder, o.expr[i], o.target[i], lambda: generate(i + 1), o.line)

    generate(0)


def transform_assert_stmt(builder: IRBuilder, a: AssertStmt) -> None:
    if builder.options.strip_asserts:
        return
    cond = builder.accept(a.expr)
    ok_block, error_block = BasicBlock(), BasicBlock()
    builder.add_bool_branch(cond, ok_block, error_block)
    builder.activate_block(error_block)
    if a.msg is None:
        # Special case (for simpler generated code)
        builder.add(RaiseStandardError(RaiseStandardError.ASSERTION_ERROR, None, a.line))
    elif isinstance(a.msg, StrExpr):
        # Another special case
        builder.add(RaiseStandardError(RaiseStandardError.ASSERTION_ERROR, a.msg.value,
                                       a.line))
    else:
        # The general case -- explicitly construct an exception instance
        message = builder.accept(a.msg)
        exc_type = builder.load_module_attr_by_fullname('builtins.AssertionError', a.line)
        exc = builder.py_call(exc_type, [message], a.line)
        builder.call_c(raise_exception_op, [exc], a.line)
    builder.add(Unreachable())
    builder.activate_block(ok_block)


def transform_del_stmt(builder: IRBuilder, o: DelStmt) -> None:
    transform_del_item(builder, builder.get_assignment_target(o.expr), o.line)


def transform_del_item(builder: IRBuilder, target: AssignmentTarget, line: int) -> None:
    if isinstance(target, AssignmentTargetIndex):
        builder.gen_method_call(
            target.base,
            '__delitem__',
            [target.index],
            result_type=None,
            line=line
        )
    elif isinstance(target, AssignmentTargetAttr):
        if isinstance(target.obj_type, RInstance):
            cl = target.obj_type.class_ir
            if not cl.is_deletable(target.attr):
                builder.error('"{}" cannot be deleted'.format(target.attr), line)
                builder.note(
                    'Using "__deletable__ = ' +
                    '[\'<attr>\']" in the class body enables "del obj.<attr>"', line)
        key = builder.load_str(target.attr)
        builder.call_c(py_delattr_op, [target.obj, key], line)
    elif isinstance(target, AssignmentTargetRegister):
        # Delete a local by assigning an error value to it, which will
        # prompt the insertion of uninit checks.
        builder.add(Assign(target.register,
                           builder.add(LoadErrorValue(target.type, undefines=True))))
    elif isinstance(target, AssignmentTargetTuple):
        for subtarget in target.items:
            transform_del_item(builder, subtarget, line)
