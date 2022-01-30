"""Generate IR for generator functions.

A generator function is represented by a class that implements the
generator protocol and keeps track of the generator state, including
local variables.

The top-level logic for dealing with generator functions is in
mypyc.irbuild.function.
"""

from typing import List

from mypy.nodes import Var, ARG_OPT

from mypyc.common import SELF_NAME, NEXT_LABEL_ATTR_NAME, ENV_ATTR_NAME
from mypyc.ir.ops import (
    BasicBlock, Call, Return, Goto, Integer, SetAttr, Unreachable, RaiseStandardError,
    Value, Register
)
from mypyc.ir.rtypes import RInstance, int_rprimitive, object_rprimitive
from mypyc.ir.func_ir import FuncIR, FuncDecl, FuncSignature, RuntimeArg
from mypyc.ir.class_ir import ClassIR
from mypyc.primitives.exc_ops import raise_exception_with_tb_op
from mypyc.irbuild.env_class import (
    add_args_to_env, load_outer_env, load_env_registers, finalize_env_class
)
from mypyc.irbuild.builder import IRBuilder, gen_arg_defaults
from mypyc.irbuild.context import FuncInfo, GeneratorClass


def gen_generator_func(builder: IRBuilder) -> None:
    setup_generator_class(builder)
    load_env_registers(builder)
    gen_arg_defaults(builder)
    finalize_env_class(builder)
    builder.add(Return(instantiate_generator_class(builder)))


def instantiate_generator_class(builder: IRBuilder) -> Value:
    fitem = builder.fn_info.fitem
    generator_reg = builder.add(Call(builder.fn_info.generator_class.ir.ctor, [], fitem.line))

    # Get the current environment register. If the current function is nested, then the
    # generator class gets instantiated from the callable class' '__call__' method, and hence
    # we use the callable class' environment register. Otherwise, we use the original
    # function's environment register.
    if builder.fn_info.is_nested:
        curr_env_reg = builder.fn_info.callable_class.curr_env_reg
    else:
        curr_env_reg = builder.fn_info.curr_env_reg

    # Set the generator class' environment attribute to point at the environment class
    # defined in the current scope.
    builder.add(SetAttr(generator_reg, ENV_ATTR_NAME, curr_env_reg, fitem.line))

    # Set the generator class' environment class' NEXT_LABEL_ATTR_NAME attribute to 0.
    zero = Integer(0)
    builder.add(SetAttr(curr_env_reg, NEXT_LABEL_ATTR_NAME, zero, fitem.line))
    return generator_reg


def setup_generator_class(builder: IRBuilder) -> ClassIR:
    name = '{}_gen'.format(builder.fn_info.namespaced_name())

    generator_class_ir = ClassIR(name, builder.module_name, is_generated=True)
    generator_class_ir.attributes[ENV_ATTR_NAME] = RInstance(builder.fn_info.env_class)
    generator_class_ir.mro = [generator_class_ir]

    builder.classes.append(generator_class_ir)
    builder.fn_info.generator_class = GeneratorClass(generator_class_ir)
    return generator_class_ir


def create_switch_for_generator_class(builder: IRBuilder) -> None:
    builder.add(Goto(builder.fn_info.generator_class.switch_block))
    block = BasicBlock()
    builder.fn_info.generator_class.continuation_blocks.append(block)
    builder.activate_block(block)


def populate_switch_for_generator_class(builder: IRBuilder) -> None:
    cls = builder.fn_info.generator_class
    line = builder.fn_info.fitem.line

    builder.activate_block(cls.switch_block)
    for label, true_block in enumerate(cls.continuation_blocks):
        false_block = BasicBlock()
        comparison = builder.binary_op(
            cls.next_label_reg, Integer(label), '==', line
        )
        builder.add_bool_branch(comparison, true_block, false_block)
        builder.activate_block(false_block)

    builder.add(RaiseStandardError(RaiseStandardError.STOP_ITERATION, None, line))
    builder.add(Unreachable())


def add_raise_exception_blocks_to_generator_class(builder: IRBuilder, line: int) -> None:
    """Add error handling blocks to a generator class.

    Generates blocks to check if error flags are set while calling the
    helper method for generator functions, and raises an exception if
    those flags are set.
    """
    cls = builder.fn_info.generator_class
    assert cls.exc_regs is not None
    exc_type, exc_val, exc_tb = cls.exc_regs

    # Check to see if an exception was raised.
    error_block = BasicBlock()
    ok_block = BasicBlock()
    comparison = builder.translate_is_op(exc_type, builder.none_object(), 'is not', line)
    builder.add_bool_branch(comparison, error_block, ok_block)

    builder.activate_block(error_block)
    builder.call_c(raise_exception_with_tb_op, [exc_type, exc_val, exc_tb], line)
    builder.add(Unreachable())
    builder.goto_and_activate(ok_block)


def add_methods_to_generator_class(builder: IRBuilder,
                                   fn_info: FuncInfo,
                                   sig: FuncSignature,
                                   arg_regs: List[Register],
                                   blocks: List[BasicBlock],
                                   is_coroutine: bool) -> None:
    helper_fn_decl = add_helper_to_generator_class(builder, arg_regs, blocks, sig, fn_info)
    add_next_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_send_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_iter_to_generator_class(builder, fn_info)
    add_throw_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_close_to_generator_class(builder, fn_info)
    if is_coroutine:
        add_await_to_generator_class(builder, fn_info)


def add_helper_to_generator_class(builder: IRBuilder,
                                  arg_regs: List[Register],
                                  blocks: List[BasicBlock],
                                  sig: FuncSignature,
                                  fn_info: FuncInfo) -> FuncDecl:
    """Generates a helper method for a generator class, called by '__next__' and 'throw'."""
    sig = FuncSignature((RuntimeArg(SELF_NAME, object_rprimitive),
                         RuntimeArg('type', object_rprimitive),
                         RuntimeArg('value', object_rprimitive),
                         RuntimeArg('traceback', object_rprimitive),
                         RuntimeArg('arg', object_rprimitive)
                         ), sig.ret_type)
    helper_fn_decl = FuncDecl('__mypyc_generator_helper__', fn_info.generator_class.ir.name,
                              builder.module_name, sig)
    helper_fn_ir = FuncIR(helper_fn_decl, arg_regs, blocks,
                          fn_info.fitem.line, traceback_name=fn_info.fitem.name)
    fn_info.generator_class.ir.methods['__mypyc_generator_helper__'] = helper_fn_ir
    builder.functions.append(helper_fn_ir)
    return helper_fn_decl


def add_iter_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__iter__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, '__iter__', object_rprimitive, fn_info):
        builder.add(Return(builder.self()))


def add_next_to_generator_class(builder: IRBuilder,
                                fn_info: FuncInfo,
                                fn_decl: FuncDecl,
                                sig: FuncSignature) -> None:
    """Generates the '__next__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, '__next__',
                              object_rprimitive, fn_info):
        none_reg = builder.none_object()
        # Call the helper function with error flags set to Py_None, and return that result.
        result = builder.add(Call(fn_decl,
                                  [builder.self(), none_reg, none_reg, none_reg, none_reg],
                                  fn_info.fitem.line))
        builder.add(Return(result))


def add_send_to_generator_class(builder: IRBuilder,
                                fn_info: FuncInfo,
                                fn_decl: FuncDecl,
                                sig: FuncSignature) -> None:
    """Generates the 'send' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, 'send', object_rprimitive, fn_info):
        arg = builder.add_argument('arg', object_rprimitive)
        none_reg = builder.none_object()
        # Call the helper function with error flags set to Py_None, and return that result.
        result = builder.add(Call(
            fn_decl,
            [builder.self(), none_reg, none_reg, none_reg, builder.read(arg)],
            fn_info.fitem.line))
        builder.add(Return(result))


def add_throw_to_generator_class(builder: IRBuilder,
                                 fn_info: FuncInfo,
                                 fn_decl: FuncDecl,
                                 sig: FuncSignature) -> None:
    """Generates the 'throw' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, 'throw',
                              object_rprimitive, fn_info):
        typ = builder.add_argument('type', object_rprimitive)
        val = builder.add_argument('value', object_rprimitive, ARG_OPT)
        tb = builder.add_argument('traceback', object_rprimitive, ARG_OPT)

        # Because the value and traceback arguments are optional and hence
        # can be NULL if not passed in, we have to assign them Py_None if
        # they are not passed in.
        none_reg = builder.none_object()
        builder.assign_if_null(val, lambda: none_reg, builder.fn_info.fitem.line)
        builder.assign_if_null(tb, lambda: none_reg, builder.fn_info.fitem.line)

        # Call the helper function using the arguments passed in, and return that result.
        result = builder.add(Call(
            fn_decl,
            [builder.self(), builder.read(typ), builder.read(val), builder.read(tb), none_reg],
            fn_info.fitem.line))
        builder.add(Return(result))


def add_close_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__close__' method for a generator class."""
    # TODO: Currently this method just triggers a runtime error.
    #       We should fill this out (https://github.com/mypyc/mypyc/issues/790).
    with builder.enter_method(fn_info.generator_class.ir, 'close', object_rprimitive, fn_info):
        builder.add(RaiseStandardError(RaiseStandardError.RUNTIME_ERROR,
                                       'close method on generator classes unimplemented',
                                       fn_info.fitem.line))
        builder.add(Unreachable())


def add_await_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__await__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, '__await__', object_rprimitive, fn_info):
        builder.add(Return(builder.self()))


def setup_env_for_generator_class(builder: IRBuilder) -> None:
    """Populates the environment for a generator class."""
    fitem = builder.fn_info.fitem
    cls = builder.fn_info.generator_class
    self_target = builder.add_self_to_env(cls.ir)

    # Add the type, value, and traceback variables to the environment.
    exc_type = builder.add_local(Var('type'), object_rprimitive, is_arg=True)
    exc_val = builder.add_local(Var('value'), object_rprimitive, is_arg=True)
    exc_tb = builder.add_local(Var('traceback'), object_rprimitive, is_arg=True)
    # TODO: Use the right type here instead of object?
    exc_arg = builder.add_local(Var('arg'), object_rprimitive, is_arg=True)

    cls.exc_regs = (exc_type, exc_val, exc_tb)
    cls.send_arg_reg = exc_arg

    cls.self_reg = builder.read(self_target, fitem.line)
    cls.curr_env_reg = load_outer_env(builder, cls.self_reg, builder.symtables[-1])

    # Define a variable representing the label to go to the next time
    # the '__next__' function of the generator is called, and add it
    # as an attribute to the environment class.
    cls.next_label_target = builder.add_var_to_env_class(
        Var(NEXT_LABEL_ATTR_NAME),
        int_rprimitive,
        cls,
        reassign=False
    )

    # Add arguments from the original generator function to the
    # environment of the generator class.
    add_args_to_env(builder, local=False, base=cls, reassign=False)

    # Set the next label register for the generator class.
    cls.next_label_reg = builder.read(cls.next_label_target, fitem.line)
