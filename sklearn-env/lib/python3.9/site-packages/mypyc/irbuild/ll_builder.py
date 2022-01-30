"""A "low-level" IR builder class.

LowLevelIRBuilder provides core abstractions we use for constructing
IR as well as a number of higher-level ones (accessing attributes,
calling functions and methods, and coercing between types, for
example). The core principle of the low-level IR builder is that all
of its facilities operate solely on the IR level and not the AST
level---it has *no knowledge* of mypy types or expressions.
"""

from typing import (
    Callable, List, Tuple, Optional, Sequence
)

from typing_extensions import Final

from mypy.nodes import ArgKind, ARG_POS, ARG_STAR, ARG_STAR2
from mypy.operators import op_methods
from mypy.types import AnyType, TypeOfAny
from mypy.checkexpr import map_actuals_to_formals

from mypyc.ir.ops import (
    BasicBlock, Op, Integer, Value, Register, Assign, Branch, Goto, Call, Box, Unbox, Cast,
    GetAttr, LoadStatic, MethodCall, CallC, Truncate, LoadLiteral, AssignMulti,
    RaiseStandardError, Unreachable, LoadErrorValue,
    NAMESPACE_TYPE, NAMESPACE_MODULE, NAMESPACE_STATIC, IntOp, GetElementPtr,
    LoadMem, ComparisonOp, LoadAddress, TupleGet, KeepAlive, ERR_NEVER, ERR_FALSE, SetMem
)
from mypyc.ir.rtypes import (
    RType, RUnion, RInstance, RArray, optional_value_type, int_rprimitive, float_rprimitive,
    bool_rprimitive, list_rprimitive, str_rprimitive, is_none_rprimitive, object_rprimitive,
    c_pyssize_t_rprimitive, is_short_int_rprimitive, is_tagged, PyVarObject, short_int_rprimitive,
    is_list_rprimitive, is_tuple_rprimitive, is_dict_rprimitive, is_set_rprimitive, PySetObject,
    none_rprimitive, RTuple, is_bool_rprimitive, is_str_rprimitive, c_int_rprimitive,
    pointer_rprimitive, PyObject, PyListObject, bit_rprimitive, is_bit_rprimitive,
    object_pointer_rprimitive, c_size_t_rprimitive, dict_rprimitive, bytes_rprimitive,
    is_bytes_rprimitive
)
from mypyc.ir.func_ir import FuncDecl, FuncSignature
from mypyc.ir.class_ir import ClassIR, all_concrete_classes
from mypyc.common import (
    FAST_ISINSTANCE_MAX_SUBCLASSES, MAX_LITERAL_SHORT_INT, MIN_LITERAL_SHORT_INT, PLATFORM_SIZE,
    use_vectorcall, use_method_vectorcall
)
from mypyc.primitives.registry import (
    method_call_ops, CFunctionDescription,
    binary_ops, unary_ops, ERR_NEG_INT
)
from mypyc.primitives.bytes_ops import bytes_compare
from mypyc.primitives.list_ops import (
    list_extend_op, new_list_op, list_build_op
)
from mypyc.primitives.tuple_ops import (
    list_tuple_op, new_tuple_op, new_tuple_with_length_op
)
from mypyc.primitives.dict_ops import (
    dict_update_in_display_op, dict_new_op, dict_build_op, dict_ssize_t_size_op
)
from mypyc.primitives.generic_ops import (
    py_getattr_op, py_call_op, py_call_with_kwargs_op, py_method_call_op,
    py_vectorcall_op, py_vectorcall_method_op,
    generic_len_op, generic_ssize_t_len_op
)
from mypyc.primitives.misc_ops import (
    none_object_op, fast_isinstance_op, bool_op
)
from mypyc.primitives.int_ops import int_comparison_op_mapping
from mypyc.primitives.exc_ops import err_occurred_op, keep_propagating_op
from mypyc.primitives.str_ops import (
    unicode_compare, str_check_if_true, str_ssize_t_size_op
)
from mypyc.primitives.set_ops import new_set_op
from mypyc.rt_subtype import is_runtime_subtype
from mypyc.subtype import is_subtype
from mypyc.sametype import is_same_type
from mypyc.irbuild.mapper import Mapper
from mypyc.options import CompilerOptions
from mypyc.irbuild.util import concrete_arg_kind


DictEntry = Tuple[Optional[Value], Value]

# If the number of items is less than the threshold when initializing
# a list, we would inline the generate IR using SetMem and expanded
# for-loop. Otherwise, we would call `list_build_op` for larger lists.
# TODO: The threshold is a randomly chosen number which needs further
#       study on real-world projects for a better balance.
LIST_BUILDING_EXPANSION_THRESHOLD = 10

# From CPython
PY_VECTORCALL_ARGUMENTS_OFFSET: Final = 1 << (PLATFORM_SIZE * 8 - 1)


class LowLevelIRBuilder:
    def __init__(
        self,
        current_module: str,
        mapper: Mapper,
        options: CompilerOptions,
    ) -> None:
        self.current_module = current_module
        self.mapper = mapper
        self.options = options
        self.args: List[Register] = []
        self.blocks: List[BasicBlock] = []
        # Stack of except handler entry blocks
        self.error_handlers: List[Optional[BasicBlock]] = [None]

    # Basic operations

    def add(self, op: Op) -> Value:
        """Add an op."""
        assert not self.blocks[-1].terminated, "Can't add to finished block"
        self.blocks[-1].ops.append(op)
        return op

    def goto(self, target: BasicBlock) -> None:
        """Add goto to a basic block."""
        if not self.blocks[-1].terminated:
            self.add(Goto(target))

    def activate_block(self, block: BasicBlock) -> None:
        """Add a basic block and make it the active one (target of adds)."""
        if self.blocks:
            assert self.blocks[-1].terminated

        block.error_handler = self.error_handlers[-1]
        self.blocks.append(block)

    def goto_and_activate(self, block: BasicBlock) -> None:
        """Add goto a block and make it the active block."""
        self.goto(block)
        self.activate_block(block)

    def push_error_handler(self, handler: Optional[BasicBlock]) -> None:
        self.error_handlers.append(handler)

    def pop_error_handler(self) -> Optional[BasicBlock]:
        return self.error_handlers.pop()

    def self(self) -> Register:
        """Return reference to the 'self' argument.

        This only works in a method.
        """
        return self.args[0]

    # Type conversions

    def box(self, src: Value) -> Value:
        if src.type.is_unboxed:
            return self.add(Box(src))
        else:
            return src

    def unbox_or_cast(self, src: Value, target_type: RType, line: int) -> Value:
        if target_type.is_unboxed:
            return self.add(Unbox(src, target_type, line))
        else:
            return self.add(Cast(src, target_type, line))

    def coerce(self, src: Value, target_type: RType, line: int, force: bool = False) -> Value:
        """Generate a coercion/cast from one type to other (only if needed).

        For example, int -> object boxes the source int; int -> int emits nothing;
        object -> int unboxes the object. All conversions preserve object value.

        If force is true, always generate an op (even if it is just an assignment) so
        that the result will have exactly target_type as the type.

        Returns the register with the converted value (may be same as src).
        """
        if src.type.is_unboxed and not target_type.is_unboxed:
            return self.box(src)
        if ((src.type.is_unboxed and target_type.is_unboxed)
                and not is_runtime_subtype(src.type, target_type)):
            # To go from one unboxed type to another, we go through a boxed
            # in-between value, for simplicity.
            tmp = self.box(src)
            return self.unbox_or_cast(tmp, target_type, line)
        if ((not src.type.is_unboxed and target_type.is_unboxed)
                or not is_subtype(src.type, target_type)):
            return self.unbox_or_cast(src, target_type, line)
        elif force:
            tmp = Register(target_type)
            self.add(Assign(tmp, src))
            return tmp
        return src

    def coerce_nullable(self, src: Value, target_type: RType, line: int) -> Value:
        """Generate a coercion from a potentially null value."""
        if (
            src.type.is_unboxed == target_type.is_unboxed
            and (
                (target_type.is_unboxed and is_runtime_subtype(src.type, target_type))
                or (not target_type.is_unboxed and is_subtype(src.type, target_type))
            )
        ):
            return src

        target = Register(target_type)

        valid, invalid, out = BasicBlock(), BasicBlock(), BasicBlock()
        self.add(Branch(src, invalid, valid, Branch.IS_ERROR))

        self.activate_block(valid)
        coerced = self.coerce(src, target_type, line)
        self.add(Assign(target, coerced, line))
        self.goto(out)

        self.activate_block(invalid)
        error = self.add(LoadErrorValue(target_type))
        self.add(Assign(target, error, line))

        self.goto_and_activate(out)
        return target

    # Attribute access

    def get_attr(self, obj: Value, attr: str, result_type: RType, line: int) -> Value:
        """Get a native or Python attribute of an object."""
        if (isinstance(obj.type, RInstance) and obj.type.class_ir.is_ext_class
                and obj.type.class_ir.has_attr(attr)):
            return self.add(GetAttr(obj, attr, line))
        elif isinstance(obj.type, RUnion):
            return self.union_get_attr(obj, obj.type, attr, result_type, line)
        else:
            return self.py_get_attr(obj, attr, line)

    def union_get_attr(self,
                       obj: Value,
                       rtype: RUnion,
                       attr: str,
                       result_type: RType,
                       line: int) -> Value:
        """Get an attribute of an object with a union type."""

        def get_item_attr(value: Value) -> Value:
            return self.get_attr(value, attr, result_type, line)

        return self.decompose_union_helper(obj, rtype, result_type, get_item_attr, line)

    def py_get_attr(self, obj: Value, attr: str, line: int) -> Value:
        """Get a Python attribute (slow).

        Prefer get_attr() which generates optimized code for native classes.
        """
        key = self.load_str(attr)
        return self.call_c(py_getattr_op, [obj, key], line)

    # isinstance() checks

    def isinstance_helper(self, obj: Value, class_irs: List[ClassIR], line: int) -> Value:
        """Fast path for isinstance() that checks against a list of native classes."""
        if not class_irs:
            return self.false()
        ret = self.isinstance_native(obj, class_irs[0], line)
        for class_ir in class_irs[1:]:
            def other() -> Value:
                return self.isinstance_native(obj, class_ir, line)
            ret = self.shortcircuit_helper('or', bool_rprimitive, lambda: ret, other, line)
        return ret

    def get_type_of_obj(self, obj: Value, line: int) -> Value:
        ob_type_address = self.add(GetElementPtr(obj, PyObject, 'ob_type', line))
        ob_type = self.add(LoadMem(object_rprimitive, ob_type_address))
        self.add(KeepAlive([obj]))
        return ob_type

    def type_is_op(self, obj: Value, type_obj: Value, line: int) -> Value:
        typ = self.get_type_of_obj(obj, line)
        return self.add(ComparisonOp(typ, type_obj, ComparisonOp.EQ, line))

    def isinstance_native(self, obj: Value, class_ir: ClassIR, line: int) -> Value:
        """Fast isinstance() check for a native class.

        If there are three or fewer concrete (non-trait) classes among the class
        and all its children, use even faster type comparison checks `type(obj)
        is typ`.
        """
        concrete = all_concrete_classes(class_ir)
        if concrete is None or len(concrete) > FAST_ISINSTANCE_MAX_SUBCLASSES + 1:
            return self.call_c(fast_isinstance_op,
                               [obj, self.get_native_type(class_ir)],
                               line)
        if not concrete:
            # There can't be any concrete instance that matches this.
            return self.false()
        type_obj = self.get_native_type(concrete[0])
        ret = self.type_is_op(obj, type_obj, line)
        for c in concrete[1:]:
            def other() -> Value:
                return self.type_is_op(obj, self.get_native_type(c), line)
            ret = self.shortcircuit_helper('or', bool_rprimitive, lambda: ret, other, line)
        return ret

    # Calls

    def _construct_varargs(self,
                           args: Sequence[Tuple[Value, ArgKind, Optional[str]]],
                           line: int,
                           *,
                           has_star: bool,
                           has_star2: bool) -> Tuple[Optional[Value], Optional[Value]]:
        """Construct *args and **kwargs from a collection of arguments

        This is pretty complicated, and almost all of the complication here stems from
        one of two things (but mostly the second):
          * The handling of ARG_STAR/ARG_STAR2. We want to create as much of the args/kwargs
            values in one go as we can, so we collect values until our hand is forced, and
            then we emit creation of the list/tuple, and expand it from there if needed.

          * Support potentially nullable argument values. This has very narrow applicability,
            as this will never be done by our compiled Python code, but is critically used
            by gen_glue_method when generating glue methods to mediate between the function
            signature of a parent class and its subclasses.

            For named-only arguments, this is quite simple: if it is
            null, don't put it in the dict.

            For positional-or-named arguments, things are much more complicated.
              * First, anything that was passed as a positional arg
                must be forwarded along as a positional arg. It *must
                not* be converted to a named arg. This is because mypy
                does not enforce that positional-or-named arguments
                have the same name in subclasses, and it is not
                uncommon for code to have different names in
                subclasses (a bunch of mypy's visitors do this, for
                example!). This is arguably a bug in both mypy and code doing
                this, and they ought to be using positional-only arguments, but
                positional-only arguments are new and ugly.

              * On the flip side, we're willing to accept the
                infelicity of sometimes turning an argument that was
                passed by keyword into a positional argument. It's wrong,
                but it's very marginal, and avoiding it would require passing
                a bitmask of which arguments were named with every function call,
                or something similar.
                (See some discussion of this in testComplicatedArgs)

            Thus, our strategy for positional-or-named arguments is to
            always pass them as positional, except in the one
            situation where we can not, and where we can be absolutely
            sure they were passed by name: when an *earlier*
            positional argument was missing its value.

            This means that if we have a method `f(self, x: int=..., y: object=...)`:
              * x and y present:      args=(x, y), kwargs={}
              * x present, y missing: args=(x,),   kwargs={}
              * x missing, y present: args=(),     kwargs={'y': y}

            To implement this, when we have multiple optional
            positional arguments, we maintain a flag in a register
            that tracks whether an argument has been missing, and for
            each such optional argument (except the first), we check
            the flag to determine whether to append the argument to
            the *args list or add it to the **kwargs dict. What a
            mess!

            This is what really makes everything here such a tangle;
            otherwise the *args and **kwargs code could be separated.

        The arguments has_star and has_star2 indicate whether the target function
        takes an ARG_STAR and ARG_STAR2 argument, respectively.
        (These will always be true when making a pycall, and be based
        on the actual target signature for a native call.)
        """

        star_result: Optional[Value] = None
        star2_result: Optional[Value] = None
        # We aggregate values that need to go into *args and **kwargs
        # in these lists. Once all arguments are processed (in the
        # happiest case), or we encounter an ARG_STAR/ARG_STAR2 or a
        # nullable arg, then we create the list and/or dict.
        star_values: List[Value] = []
        star2_keys: List[Value] = []
        star2_values: List[Value] = []

        seen_empty_reg: Optional[Register] = None

        for value, kind, name in args:
            if kind == ARG_STAR:
                if star_result is None:
                    star_result = self.new_list_op(star_values, line)
                self.call_c(list_extend_op, [star_result, value], line)
            elif kind == ARG_STAR2:
                if star2_result is None:
                    star2_result = self._create_dict(star2_keys, star2_values, line)

                self.call_c(
                    dict_update_in_display_op,
                    [star2_result, value],
                    line=line
                )
            else:
                nullable = kind.is_optional()
                maybe_pos = kind.is_positional() and has_star
                maybe_named = kind.is_named() or (kind.is_optional() and name and has_star2)

                # If the argument is nullable, we need to create the
                # relevant args/kwargs objects so that we can
                # conditionally modify them.
                if nullable:
                    if maybe_pos and star_result is None:
                        star_result = self.new_list_op(star_values, line)
                    if maybe_named and star2_result is None:
                        star2_result = self._create_dict(star2_keys, star2_values, line)

                # Easy cases: just collect the argument.
                if maybe_pos and star_result is None:
                    star_values.append(value)
                    continue

                if maybe_named and star2_result is None:
                    assert name is not None
                    key = self.load_str(name)
                    star2_keys.append(key)
                    star2_values.append(value)
                    continue

                # OK, anything that is nullable or *after* a nullable arg needs to be here
                # TODO: We could try harder to avoid creating basic blocks in the common case
                new_seen_empty_reg = seen_empty_reg

                out = BasicBlock()
                if nullable:
                    # If this is the first nullable positional arg we've seen, create
                    # a register to track whether anything has been null.
                    # (We won't *check* the register until the next argument, though.)
                    if maybe_pos and not seen_empty_reg:
                        new_seen_empty_reg = Register(bool_rprimitive)
                        self.add(Assign(new_seen_empty_reg, self.false(), line))

                    skip = BasicBlock() if maybe_pos else out
                    keep = BasicBlock()
                    self.add(Branch(value, skip, keep, Branch.IS_ERROR))
                    self.activate_block(keep)

                # If this could be positional or named and we /might/ have seen a missing
                # positional arg, then we need to compile *both* a positional and named
                # version! What a pain!
                if maybe_pos and maybe_named and seen_empty_reg:
                    pos_block, named_block = BasicBlock(), BasicBlock()
                    self.add(Branch(seen_empty_reg, named_block, pos_block, Branch.BOOL))
                else:
                    pos_block = named_block = BasicBlock()
                    self.goto(pos_block)

                if maybe_pos:
                    self.activate_block(pos_block)
                    assert star_result
                    self.translate_special_method_call(
                        star_result, 'append', [value], result_type=None, line=line)
                    self.goto(out)

                if maybe_named and (not maybe_pos or seen_empty_reg):
                    self.activate_block(named_block)
                    assert name is not None
                    key = self.load_str(name)
                    assert star2_result
                    self.translate_special_method_call(
                        star2_result, '__setitem__', [key, value], result_type=None, line=line)
                    self.goto(out)

                if nullable and maybe_pos and new_seen_empty_reg:
                    assert skip is not out
                    self.activate_block(skip)
                    self.add(Assign(new_seen_empty_reg, self.true(), line))
                    self.goto(out)

                self.activate_block(out)

                seen_empty_reg = new_seen_empty_reg

        assert not (star_result or star_values) or has_star
        assert not (star2_result or star2_values) or has_star2
        if has_star:
            # If we managed to make it this far without creating a
            # *args list, then we can directly create a
            # tuple. Otherwise create the tuple from the list.
            if star_result is None:
                star_result = self.new_tuple(star_values, line)
            else:
                star_result = self.call_c(list_tuple_op, [star_result], line)
        if has_star2 and star2_result is None:
            star2_result = self._create_dict(star2_keys, star2_values, line)

        return star_result, star2_result

    def py_call(self,
                function: Value,
                arg_values: List[Value],
                line: int,
                arg_kinds: Optional[List[ArgKind]] = None,
                arg_names: Optional[Sequence[Optional[str]]] = None) -> Value:
        """Call a Python function (non-native and slow).

        Use py_call_op or py_call_with_kwargs_op for Python function call.
        """
        if use_vectorcall(self.options.capi_version):
            # More recent Python versions support faster vectorcalls.
            result = self._py_vector_call(function, arg_values, line, arg_kinds, arg_names)
            if result is not None:
                return result

        # If all arguments are positional, we can use py_call_op.
        if arg_kinds is None or all(kind == ARG_POS for kind in arg_kinds):
            return self.call_c(py_call_op, [function] + arg_values, line)

        # Otherwise fallback to py_call_with_kwargs_op.
        assert arg_names is not None

        pos_args_tuple, kw_args_dict = self._construct_varargs(
            list(zip(arg_values, arg_kinds, arg_names)), line, has_star=True, has_star2=True
        )
        assert pos_args_tuple and kw_args_dict

        return self.call_c(
            py_call_with_kwargs_op, [function, pos_args_tuple, kw_args_dict], line)

    def _py_vector_call(self,
                        function: Value,
                        arg_values: List[Value],
                        line: int,
                        arg_kinds: Optional[List[ArgKind]] = None,
                        arg_names: Optional[Sequence[Optional[str]]] = None) -> Optional[Value]:
        """Call function using the vectorcall API if possible.

        Return the return value if successful. Return None if a non-vectorcall
        API should be used instead.
        """
        # We can do this if all args are positional or named (no *args or **kwargs, not optional).
        if arg_kinds is None or all(not kind.is_star() and not kind.is_optional()
                                    for kind in arg_kinds):
            if arg_values:
                # Create a C array containing all arguments as boxed values.
                array = Register(RArray(object_rprimitive, len(arg_values)))
                coerced_args = [self.coerce(arg, object_rprimitive, line) for arg in arg_values]
                self.add(AssignMulti(array, coerced_args))
                arg_ptr = self.add(LoadAddress(object_pointer_rprimitive, array))
            else:
                arg_ptr = Integer(0, object_pointer_rprimitive)
            num_pos = num_positional_args(arg_values, arg_kinds)
            keywords = self._vectorcall_keywords(arg_names)
            value = self.call_c(py_vectorcall_op, [function,
                                                   arg_ptr,
                                                   Integer(num_pos, c_size_t_rprimitive),
                                                   keywords],
                                line)
            if arg_values:
                # Make sure arguments won't be freed until after the call.
                # We need this because RArray doesn't support automatic
                # memory management.
                self.add(KeepAlive(coerced_args))
            return value
        return None

    def _vectorcall_keywords(self, arg_names: Optional[Sequence[Optional[str]]]) -> Value:
        """Return a reference to a tuple literal with keyword argument names.

        Return null pointer if there are no keyword arguments.
        """
        if arg_names:
            kw_list = [name for name in arg_names if name is not None]
            if kw_list:
                return self.add(LoadLiteral(tuple(kw_list), object_rprimitive))
        return Integer(0, object_rprimitive)

    def py_method_call(self,
                       obj: Value,
                       method_name: str,
                       arg_values: List[Value],
                       line: int,
                       arg_kinds: Optional[List[ArgKind]],
                       arg_names: Optional[Sequence[Optional[str]]]) -> Value:
        """Call a Python method (non-native and slow)."""
        if use_method_vectorcall(self.options.capi_version):
            # More recent Python versions support faster vectorcalls.
            result = self._py_vector_method_call(
                obj, method_name, arg_values, line, arg_kinds, arg_names)
            if result is not None:
                return result

        if arg_kinds is None or all(kind == ARG_POS for kind in arg_kinds):
            # Use legacy method call API
            method_name_reg = self.load_str(method_name)
            return self.call_c(py_method_call_op, [obj, method_name_reg] + arg_values, line)
        else:
            # Use py_call since it supports keyword arguments (and vectorcalls).
            method = self.py_get_attr(obj, method_name, line)
            return self.py_call(method, arg_values, line, arg_kinds=arg_kinds, arg_names=arg_names)

    def _py_vector_method_call(self,
                               obj: Value,
                               method_name: str,
                               arg_values: List[Value],
                               line: int,
                               arg_kinds: Optional[List[ArgKind]],
                               arg_names: Optional[Sequence[Optional[str]]]) -> Optional[Value]:
        """Call method using the vectorcall API if possible.

        Return the return value if successful. Return None if a non-vectorcall
        API should be used instead.
        """
        if arg_kinds is None or all(not kind.is_star() and not kind.is_optional()
                                    for kind in arg_kinds):
            method_name_reg = self.load_str(method_name)
            array = Register(RArray(object_rprimitive, len(arg_values) + 1))
            self_arg = self.coerce(obj, object_rprimitive, line)
            coerced_args = [self_arg] + [self.coerce(arg, object_rprimitive, line)
                                         for arg in arg_values]
            self.add(AssignMulti(array, coerced_args))
            arg_ptr = self.add(LoadAddress(object_pointer_rprimitive, array))
            num_pos = num_positional_args(arg_values, arg_kinds)
            keywords = self._vectorcall_keywords(arg_names)
            value = self.call_c(py_vectorcall_method_op,
                                [method_name_reg,
                                 arg_ptr,
                                 Integer((num_pos + 1) | PY_VECTORCALL_ARGUMENTS_OFFSET,
                                         c_size_t_rprimitive),
                                 keywords],
                                line)
            # Make sure arguments won't be freed until after the call.
            # We need this because RArray doesn't support automatic
            # memory management.
            self.add(KeepAlive(coerced_args))
            return value
        return None

    def call(self,
             decl: FuncDecl,
             args: Sequence[Value],
             arg_kinds: List[ArgKind],
             arg_names: Sequence[Optional[str]],
             line: int) -> Value:
        """Call a native function."""
        # Normalize args to positionals.
        args = self.native_args_to_positional(
            args, arg_kinds, arg_names, decl.sig, line)
        return self.add(Call(decl, args, line))

    def native_args_to_positional(self,
                                  args: Sequence[Value],
                                  arg_kinds: List[ArgKind],
                                  arg_names: Sequence[Optional[str]],
                                  sig: FuncSignature,
                                  line: int) -> List[Value]:
        """Prepare arguments for a native call.

        Given args/kinds/names and a target signature for a native call, map
        keyword arguments to their appropriate place in the argument list,
        fill in error values for unspecified default arguments,
        package arguments that will go into *args/**kwargs into a tuple/dict,
        and coerce arguments to the appropriate type.
        """

        sig_arg_kinds = [arg.kind for arg in sig.args]
        sig_arg_names = [arg.name for arg in sig.args]
        concrete_kinds = [concrete_arg_kind(arg_kind) for arg_kind in arg_kinds]
        formal_to_actual = map_actuals_to_formals(concrete_kinds,
                                                  arg_names,
                                                  sig_arg_kinds,
                                                  sig_arg_names,
                                                  lambda n: AnyType(TypeOfAny.special_form))

        # First scan for */** and construct those
        has_star = has_star2 = False
        star_arg_entries = []
        for lst, arg in zip(formal_to_actual, sig.args):
            if arg.kind.is_star():
                star_arg_entries.extend([(args[i], arg_kinds[i], arg_names[i]) for i in lst])
            has_star = has_star or arg.kind == ARG_STAR
            has_star2 = has_star2 or arg.kind == ARG_STAR2

        star_arg, star2_arg = self._construct_varargs(
            star_arg_entries, line, has_star=has_star, has_star2=has_star2
        )

        # Flatten out the arguments, loading error values for default
        # arguments, constructing tuples/dicts for star args, and
        # coercing everything to the expected type.
        output_args = []
        for lst, arg in zip(formal_to_actual, sig.args):
            if arg.kind == ARG_STAR:
                assert star_arg
                output_arg = star_arg
            elif arg.kind == ARG_STAR2:
                assert star2_arg
                output_arg = star2_arg
            elif not lst:
                output_arg = self.add(LoadErrorValue(arg.type, is_borrowed=True))
            else:
                base_arg = args[lst[0]]

                if arg_kinds[lst[0]].is_optional():
                    output_arg = self.coerce_nullable(base_arg, arg.type, line)
                else:
                    output_arg = self.coerce(base_arg, arg.type, line)

            output_args.append(output_arg)

        return output_args

    def gen_method_call(self,
                        base: Value,
                        name: str,
                        arg_values: List[Value],
                        result_type: Optional[RType],
                        line: int,
                        arg_kinds: Optional[List[ArgKind]] = None,
                        arg_names: Optional[List[Optional[str]]] = None) -> Value:
        """Generate either a native or Python method call."""
        # If we have *args, then fallback to Python method call.
        if arg_kinds is not None and any(kind.is_star() for kind in arg_kinds):
            return self.py_method_call(base, name, arg_values, base.line, arg_kinds, arg_names)

        # If the base type is one of ours, do a MethodCall
        if (isinstance(base.type, RInstance) and base.type.class_ir.is_ext_class
                and not base.type.class_ir.builtin_base):
            if base.type.class_ir.has_method(name):
                decl = base.type.class_ir.method_decl(name)
                if arg_kinds is None:
                    assert arg_names is None, "arg_kinds not present but arg_names is"
                    arg_kinds = [ARG_POS for _ in arg_values]
                    arg_names = [None for _ in arg_values]
                else:
                    assert arg_names is not None, "arg_kinds present but arg_names is not"

                # Normalize args to positionals.
                assert decl.bound_sig
                arg_values = self.native_args_to_positional(
                    arg_values, arg_kinds, arg_names, decl.bound_sig, line)
                return self.add(MethodCall(base, name, arg_values, line))
            elif base.type.class_ir.has_attr(name):
                function = self.add(GetAttr(base, name, line))
                return self.py_call(function, arg_values, line,
                                    arg_kinds=arg_kinds, arg_names=arg_names)

        elif isinstance(base.type, RUnion):
            return self.union_method_call(base, base.type, name, arg_values, result_type, line,
                                          arg_kinds, arg_names)

        # Try to do a special-cased method call
        if not arg_kinds or arg_kinds == [ARG_POS] * len(arg_values):
            target = self.translate_special_method_call(base, name, arg_values, result_type, line)
            if target:
                return target

        # Fall back to Python method call
        return self.py_method_call(base, name, arg_values, line, arg_kinds, arg_names)

    def union_method_call(self,
                          base: Value,
                          obj_type: RUnion,
                          name: str,
                          arg_values: List[Value],
                          return_rtype: Optional[RType],
                          line: int,
                          arg_kinds: Optional[List[ArgKind]],
                          arg_names: Optional[List[Optional[str]]]) -> Value:
        """Generate a method call with a union type for the object."""
        # Union method call needs a return_rtype for the type of the output register.
        # If we don't have one, use object_rprimitive.
        return_rtype = return_rtype or object_rprimitive

        def call_union_item(value: Value) -> Value:
            return self.gen_method_call(value, name, arg_values, return_rtype, line,
                                        arg_kinds, arg_names)

        return self.decompose_union_helper(base, obj_type, return_rtype, call_union_item, line)

    # Loading various values

    def none(self) -> Value:
        """Load unboxed None value (type: none_rprimitive)."""
        return Integer(1, none_rprimitive)

    def true(self) -> Value:
        """Load unboxed True value (type: bool_rprimitive)."""
        return Integer(1, bool_rprimitive)

    def false(self) -> Value:
        """Load unboxed False value (type: bool_rprimitive)."""
        return Integer(0, bool_rprimitive)

    def none_object(self) -> Value:
        """Load Python None value (type: object_rprimitive)."""
        return self.add(LoadAddress(none_object_op.type, none_object_op.src, line=-1))

    def load_int(self, value: int) -> Value:
        """Load a tagged (Python) integer literal value."""
        if value > MAX_LITERAL_SHORT_INT or value < MIN_LITERAL_SHORT_INT:
            return self.add(LoadLiteral(value, int_rprimitive))
        else:
            return Integer(value)

    def load_float(self, value: float) -> Value:
        """Load a float literal value."""
        return self.add(LoadLiteral(value, float_rprimitive))

    def load_str(self, value: str) -> Value:
        """Load a str literal value.

        This is useful for more than just str literals; for example, method calls
        also require a PyObject * form for the name of the method.
        """
        return self.add(LoadLiteral(value, str_rprimitive))

    def load_bytes(self, value: bytes) -> Value:
        """Load a bytes literal value."""
        return self.add(LoadLiteral(value, bytes_rprimitive))

    def load_complex(self, value: complex) -> Value:
        """Load a complex literal value."""
        return self.add(LoadLiteral(value, object_rprimitive))

    def load_static_checked(self, typ: RType, identifier: str, module_name: Optional[str] = None,
                            namespace: str = NAMESPACE_STATIC,
                            line: int = -1,
                            error_msg: Optional[str] = None) -> Value:
        if error_msg is None:
            error_msg = 'name "{}" is not defined'.format(identifier)
        ok_block, error_block = BasicBlock(), BasicBlock()
        value = self.add(LoadStatic(typ, identifier, module_name, namespace, line=line))
        self.add(Branch(value, error_block, ok_block, Branch.IS_ERROR, rare=True))
        self.activate_block(error_block)
        self.add(RaiseStandardError(RaiseStandardError.NAME_ERROR,
                                    error_msg,
                                    line))
        self.add(Unreachable())
        self.activate_block(ok_block)
        return value

    def load_module(self, name: str) -> Value:
        return self.add(LoadStatic(object_rprimitive, name, namespace=NAMESPACE_MODULE))

    def get_native_type(self, cls: ClassIR) -> Value:
        """Load native type object."""
        fullname = '%s.%s' % (cls.module_name, cls.name)
        return self.load_native_type_object(fullname)

    def load_native_type_object(self, fullname: str) -> Value:
        module, name = fullname.rsplit('.', 1)
        return self.add(LoadStatic(object_rprimitive, name, module, NAMESPACE_TYPE))

    # Other primitive operations
    def binary_op(self, lreg: Value, rreg: Value, op: str, line: int) -> Value:
        ltype = lreg.type
        rtype = rreg.type

        # Special case tuple comparison here so that nested tuples can be supported
        if isinstance(ltype, RTuple) and isinstance(rtype, RTuple) and op in ('==', '!='):
            return self.compare_tuples(lreg, rreg, op, line)

        # Special case == and != when we can resolve the method call statically
        if op in ('==', '!='):
            value = self.translate_eq_cmp(lreg, rreg, op, line)
            if value is not None:
                return value

        # Special case various ops
        if op in ('is', 'is not'):
            return self.translate_is_op(lreg, rreg, op, line)
        # TODO: modify 'str' to use same interface as 'compare_bytes' as it avoids
        # call to PyErr_Occurred()
        if is_str_rprimitive(ltype) and is_str_rprimitive(rtype) and op in ('==', '!='):
            return self.compare_strings(lreg, rreg, op, line)
        if is_bytes_rprimitive(ltype) and is_bytes_rprimitive(rtype) and op in ('==', '!='):
            return self.compare_bytes(lreg, rreg, op, line)
        if is_tagged(ltype) and is_tagged(rtype) and op in int_comparison_op_mapping:
            return self.compare_tagged(lreg, rreg, op, line)
        if is_bool_rprimitive(ltype) and is_bool_rprimitive(rtype) and op in (
                '&', '&=', '|', '|=', '^', '^='):
            return self.bool_bitwise_op(lreg, rreg, op[0], line)
        if isinstance(rtype, RInstance) and op in ('in', 'not in'):
            return self.translate_instance_contains(rreg, lreg, op, line)

        call_c_ops_candidates = binary_ops.get(op, [])
        target = self.matching_call_c(call_c_ops_candidates, [lreg, rreg], line)
        assert target, 'Unsupported binary operation: %s' % op
        return target

    def check_tagged_short_int(self, val: Value, line: int, negated: bool = False) -> Value:
        """Check if a tagged integer is a short integer.

        Return the result of the check (value of type 'bit').
        """
        int_tag = Integer(1, c_pyssize_t_rprimitive, line)
        bitwise_and = self.int_op(c_pyssize_t_rprimitive, val, int_tag, IntOp.AND, line)
        zero = Integer(0, c_pyssize_t_rprimitive, line)
        op = ComparisonOp.NEQ if negated else ComparisonOp.EQ
        check = self.comparison_op(bitwise_and, zero, op, line)
        return check

    def compare_tagged(self, lhs: Value, rhs: Value, op: str, line: int) -> Value:
        """Compare two tagged integers using given operator (value context)."""
        # generate fast binary logic ops on short ints
        if is_short_int_rprimitive(lhs.type) and is_short_int_rprimitive(rhs.type):
            return self.comparison_op(lhs, rhs, int_comparison_op_mapping[op][0], line)
        op_type, c_func_desc, negate_result, swap_op = int_comparison_op_mapping[op]
        result = Register(bool_rprimitive)
        short_int_block, int_block, out = BasicBlock(), BasicBlock(), BasicBlock()
        check_lhs = self.check_tagged_short_int(lhs, line)
        if op in ("==", "!="):
            check = check_lhs
        else:
            # for non-equality logical ops (less/greater than, etc.), need to check both sides
            check_rhs = self.check_tagged_short_int(rhs, line)
            check = self.int_op(bit_rprimitive, check_lhs, check_rhs, IntOp.AND, line)
        self.add(Branch(check, short_int_block, int_block, Branch.BOOL))
        self.activate_block(short_int_block)
        eq = self.comparison_op(lhs, rhs, op_type, line)
        self.add(Assign(result, eq, line))
        self.goto(out)
        self.activate_block(int_block)
        if swap_op:
            args = [rhs, lhs]
        else:
            args = [lhs, rhs]
        call = self.call_c(c_func_desc, args, line)
        if negate_result:
            # TODO: introduce UnaryIntOp?
            call_result = self.unary_op(call, "not", line)
        else:
            call_result = call
        self.add(Assign(result, call_result, line))
        self.goto_and_activate(out)
        return result

    def compare_tagged_condition(self,
                                 lhs: Value,
                                 rhs: Value,
                                 op: str,
                                 true: BasicBlock,
                                 false: BasicBlock,
                                 line: int) -> None:
        """Compare two tagged integers using given operator (conditional context).

        Assume lhs and and rhs are tagged integers.

        Args:
            lhs: Left operand
            rhs: Right operand
            op: Operation, one of '==', '!=', '<', '<=', '>', '<='
            true: Branch target if comparison is true
            false: Branch target if comparison is false
        """
        is_eq = op in ("==", "!=")
        if ((is_short_int_rprimitive(lhs.type) and is_short_int_rprimitive(rhs.type))
            or (is_eq and (is_short_int_rprimitive(lhs.type) or
                           is_short_int_rprimitive(rhs.type)))):
            # We can skip the tag check
            check = self.comparison_op(lhs, rhs, int_comparison_op_mapping[op][0], line)
            self.add(Branch(check, true, false, Branch.BOOL))
            return
        op_type, c_func_desc, negate_result, swap_op = int_comparison_op_mapping[op]
        int_block, short_int_block = BasicBlock(), BasicBlock()
        check_lhs = self.check_tagged_short_int(lhs, line, negated=True)
        if is_eq or is_short_int_rprimitive(rhs.type):
            self.add(Branch(check_lhs, int_block, short_int_block, Branch.BOOL))
        else:
            # For non-equality logical ops (less/greater than, etc.), need to check both sides
            rhs_block = BasicBlock()
            self.add(Branch(check_lhs, int_block, rhs_block, Branch.BOOL))
            self.activate_block(rhs_block)
            check_rhs = self.check_tagged_short_int(rhs, line, negated=True)
            self.add(Branch(check_rhs, int_block, short_int_block, Branch.BOOL))
        # Arbitrary integers (slow path)
        self.activate_block(int_block)
        if swap_op:
            args = [rhs, lhs]
        else:
            args = [lhs, rhs]
        call = self.call_c(c_func_desc, args, line)
        if negate_result:
            self.add(Branch(call, false, true, Branch.BOOL))
        else:
            self.add(Branch(call, true, false, Branch.BOOL))
        # Short integers (fast path)
        self.activate_block(short_int_block)
        eq = self.comparison_op(lhs, rhs, op_type, line)
        self.add(Branch(eq, true, false, Branch.BOOL))

    def compare_strings(self, lhs: Value, rhs: Value, op: str, line: int) -> Value:
        """Compare two strings"""
        compare_result = self.call_c(unicode_compare, [lhs, rhs], line)
        error_constant = Integer(-1, c_int_rprimitive, line)
        compare_error_check = self.add(ComparisonOp(compare_result,
                                                    error_constant, ComparisonOp.EQ, line))
        exception_check, propagate, final_compare = BasicBlock(), BasicBlock(), BasicBlock()
        branch = Branch(compare_error_check, exception_check, final_compare, Branch.BOOL)
        branch.negated = False
        self.add(branch)
        self.activate_block(exception_check)
        check_error_result = self.call_c(err_occurred_op, [], line)
        null = Integer(0, pointer_rprimitive, line)
        compare_error_check = self.add(ComparisonOp(check_error_result,
                                                    null, ComparisonOp.NEQ, line))
        branch = Branch(compare_error_check, propagate, final_compare, Branch.BOOL)
        branch.negated = False
        self.add(branch)
        self.activate_block(propagate)
        self.call_c(keep_propagating_op, [], line)
        self.goto(final_compare)
        self.activate_block(final_compare)
        op_type = ComparisonOp.EQ if op == '==' else ComparisonOp.NEQ
        return self.add(ComparisonOp(compare_result,
                                     Integer(0, c_int_rprimitive), op_type, line))

    def compare_bytes(self, lhs: Value, rhs: Value, op: str, line: int) -> Value:
        compare_result = self.call_c(bytes_compare, [lhs, rhs], line)
        op_type = ComparisonOp.EQ if op == '==' else ComparisonOp.NEQ
        return self.add(ComparisonOp(compare_result,
                                     Integer(1, c_int_rprimitive), op_type, line))

    def compare_tuples(self,
                       lhs: Value,
                       rhs: Value,
                       op: str,
                       line: int = -1) -> Value:
        """Compare two tuples item by item"""
        # type cast to pass mypy check
        assert isinstance(lhs.type, RTuple) and isinstance(rhs.type, RTuple)
        equal = True if op == '==' else False
        result = Register(bool_rprimitive)
        # empty tuples
        if len(lhs.type.types) == 0 and len(rhs.type.types) == 0:
            self.add(Assign(result, self.true() if equal else self.false(), line))
            return result
        length = len(lhs.type.types)
        false_assign, true_assign, out = BasicBlock(), BasicBlock(), BasicBlock()
        check_blocks = [BasicBlock() for _ in range(length)]
        lhs_items = [self.add(TupleGet(lhs, i, line)) for i in range(length)]
        rhs_items = [self.add(TupleGet(rhs, i, line)) for i in range(length)]

        if equal:
            early_stop, final = false_assign, true_assign
        else:
            early_stop, final = true_assign, false_assign

        for i in range(len(lhs.type.types)):
            if i != 0:
                self.activate_block(check_blocks[i])
            lhs_item = lhs_items[i]
            rhs_item = rhs_items[i]
            compare = self.binary_op(lhs_item, rhs_item, op, line)
            # Cast to bool if necessary since most types uses comparison returning a object type
            # See generic_ops.py for more information
            if not is_bool_rprimitive(compare.type):
                compare = self.call_c(bool_op, [compare], line)
            if i < len(lhs.type.types) - 1:
                branch = Branch(compare, early_stop, check_blocks[i + 1], Branch.BOOL)
            else:
                branch = Branch(compare, early_stop, final, Branch.BOOL)
            # if op is ==, we branch on false, else branch on true
            branch.negated = equal
            self.add(branch)
        self.activate_block(false_assign)
        self.add(Assign(result, self.false(), line))
        self.goto(out)
        self.activate_block(true_assign)
        self.add(Assign(result, self.true(), line))
        self.goto_and_activate(out)
        return result

    def translate_instance_contains(self, inst: Value, item: Value, op: str, line: int) -> Value:
        res = self.gen_method_call(inst, '__contains__', [item], None, line)
        if not is_bool_rprimitive(res.type):
            res = self.call_c(bool_op, [res], line)
        if op == 'not in':
            res = self.bool_bitwise_op(res, Integer(1, rtype=bool_rprimitive), '^', line)
        return res

    def bool_bitwise_op(self, lreg: Value, rreg: Value, op: str, line: int) -> Value:
        if op == '&':
            code = IntOp.AND
        elif op == '|':
            code = IntOp.OR
        elif op == '^':
            code = IntOp.XOR
        else:
            assert False, op
        return self.add(IntOp(bool_rprimitive, lreg, rreg, code, line))

    def unary_not(self,
                  value: Value,
                  line: int) -> Value:
        mask = Integer(1, value.type, line)
        return self.int_op(value.type, value, mask, IntOp.XOR, line)

    def unary_op(self,
                 value: Value,
                 expr_op: str,
                 line: int) -> Value:
        typ = value.type
        if (is_bool_rprimitive(typ) or is_bit_rprimitive(typ)) and expr_op == 'not':
            return self.unary_not(value, line)
        if isinstance(typ, RInstance):
            if expr_op == '-':
                method = '__neg__'
            elif expr_op == '~':
                method = '__invert__'
            else:
                method = ''
            if method and typ.class_ir.has_method(method):
                return self.gen_method_call(value, method, [], None, line)
        call_c_ops_candidates = unary_ops.get(expr_op, [])
        target = self.matching_call_c(call_c_ops_candidates, [value], line)
        assert target, 'Unsupported unary operation: %s' % expr_op
        return target

    def make_dict(self, key_value_pairs: Sequence[DictEntry], line: int) -> Value:
        result: Optional[Value] = None
        keys: List[Value] = []
        values: List[Value] = []
        for key, value in key_value_pairs:
            if key is not None:
                # key:value
                if result is None:
                    keys.append(key)
                    values.append(value)
                    continue

                self.translate_special_method_call(
                    result,
                    '__setitem__',
                    [key, value],
                    result_type=None,
                    line=line)
            else:
                # **value
                if result is None:
                    result = self._create_dict(keys, values, line)

                self.call_c(
                    dict_update_in_display_op,
                    [result, value],
                    line=line
                )

        if result is None:
            result = self._create_dict(keys, values, line)

        return result

    def new_list_op_with_length(self, length: Value, line: int) -> Value:
        """This function returns an uninitialized list.

        If the length is non-zero, the caller must initialize the list, before
        it can be made visible to user code -- otherwise the list object is broken.
        You might need further initialization with `new_list_set_item_op` op.

        Args:
            length: desired length of the new list. The rtype should be
                    c_pyssize_t_rprimitive
            line: line number
        """
        return self.call_c(new_list_op, [length], line)

    def new_list_op(self, values: List[Value], line: int) -> Value:
        length: List[Value] = [Integer(len(values), c_pyssize_t_rprimitive, line)]
        if len(values) >= LIST_BUILDING_EXPANSION_THRESHOLD:
            return self.call_c(list_build_op, length + values, line)

        # If the length of the list is less than the threshold,
        # LIST_BUILDING_EXPANSION_THRESHOLD, we directly expand the
        # for-loop and inline the SetMem operation, which is faster
        # than list_build_op, however generates more code.
        result_list = self.call_c(new_list_op, length, line)
        if len(values) == 0:
            return result_list
        args = [self.coerce(item, object_rprimitive, line) for item in values]
        ob_item_ptr = self.add(GetElementPtr(result_list, PyListObject, 'ob_item', line))
        ob_item_base = self.add(LoadMem(pointer_rprimitive, ob_item_ptr, line))
        for i in range(len(values)):
            if i == 0:
                item_address = ob_item_base
            else:
                offset = Integer(PLATFORM_SIZE * i, c_pyssize_t_rprimitive, line)
                item_address = self.add(IntOp(pointer_rprimitive, ob_item_base, offset,
                                              IntOp.ADD, line))
            self.add(SetMem(object_rprimitive, item_address, args[i], line))
        self.add(KeepAlive([result_list]))
        return result_list

    def new_set_op(self, values: List[Value], line: int) -> Value:
        return self.call_c(new_set_op, values, line)

    def shortcircuit_helper(self, op: str,
                            expr_type: RType,
                            left: Callable[[], Value],
                            right: Callable[[], Value], line: int) -> Value:
        # Having actual Phi nodes would be really nice here!
        target = Register(expr_type)
        # left_body takes the value of the left side, right_body the right
        left_body, right_body, next_block = BasicBlock(), BasicBlock(), BasicBlock()
        # true_body is taken if the left is true, false_body if it is false.
        # For 'and' the value is the right side if the left is true, and for 'or'
        # it is the right side if the left is false.
        true_body, false_body = (
            (right_body, left_body) if op == 'and' else (left_body, right_body))

        left_value = left()
        self.add_bool_branch(left_value, true_body, false_body)

        self.activate_block(left_body)
        left_coerced = self.coerce(left_value, expr_type, line)
        self.add(Assign(target, left_coerced))
        self.goto(next_block)

        self.activate_block(right_body)
        right_value = right()
        right_coerced = self.coerce(right_value, expr_type, line)
        self.add(Assign(target, right_coerced))
        self.goto(next_block)

        self.activate_block(next_block)
        return target

    def add_bool_branch(self, value: Value, true: BasicBlock, false: BasicBlock) -> None:
        if is_runtime_subtype(value.type, int_rprimitive):
            zero = Integer(0, short_int_rprimitive)
            self.compare_tagged_condition(value, zero, '!=', true, false, value.line)
            return
        elif is_same_type(value.type, str_rprimitive):
            value = self.call_c(str_check_if_true, [value], value.line)
        elif (is_same_type(value.type, list_rprimitive)
                or is_same_type(value.type, dict_rprimitive)):
            length = self.builtin_len(value, value.line)
            zero = Integer(0)
            value = self.binary_op(length, zero, '!=', value.line)
        elif (isinstance(value.type, RInstance) and value.type.class_ir.is_ext_class
                and value.type.class_ir.has_method('__bool__')):
            # Directly call the __bool__ method on classes that have it.
            value = self.gen_method_call(value, '__bool__', [], bool_rprimitive, value.line)
        else:
            value_type = optional_value_type(value.type)
            if value_type is not None:
                is_none = self.translate_is_op(value, self.none_object(), 'is not', value.line)
                branch = Branch(is_none, true, false, Branch.BOOL)
                self.add(branch)
                always_truthy = False
                if isinstance(value_type, RInstance):
                    # check whether X.__bool__ is always just the default (object.__bool__)
                    if (not value_type.class_ir.has_method('__bool__')
                            and value_type.class_ir.is_method_final('__bool__')):
                        always_truthy = True

                if not always_truthy:
                    # Optional[X] where X may be falsey and requires a check
                    branch.true = BasicBlock()
                    self.activate_block(branch.true)
                    # unbox_or_cast instead of coerce because we want the
                    # type to change even if it is a subtype.
                    remaining = self.unbox_or_cast(value, value_type, value.line)
                    self.add_bool_branch(remaining, true, false)
                return
            elif not is_bool_rprimitive(value.type) and not is_bit_rprimitive(value.type):
                value = self.call_c(bool_op, [value], value.line)
        self.add(Branch(value, true, false, Branch.BOOL))

    def call_c(self,
               desc: CFunctionDescription,
               args: List[Value],
               line: int,
               result_type: Optional[RType] = None) -> Value:
        """Call function using C/native calling convention (not a Python callable)."""
        # Handle void function via singleton RVoid instance
        coerced = []
        # Coerce fixed number arguments
        for i in range(min(len(args), len(desc.arg_types))):
            formal_type = desc.arg_types[i]
            arg = args[i]
            arg = self.coerce(arg, formal_type, line)
            coerced.append(arg)
        # Reorder args if necessary
        if desc.ordering is not None:
            assert desc.var_arg_type is None
            coerced = [coerced[i] for i in desc.ordering]
        # Coerce any var_arg
        var_arg_idx = -1
        if desc.var_arg_type is not None:
            var_arg_idx = len(desc.arg_types)
            for i in range(len(desc.arg_types), len(args)):
                arg = args[i]
                arg = self.coerce(arg, desc.var_arg_type, line)
                coerced.append(arg)
        # Add extra integer constant if any
        for item in desc.extra_int_constants:
            val, typ = item
            extra_int_constant = Integer(val, typ, line)
            coerced.append(extra_int_constant)
        error_kind = desc.error_kind
        if error_kind == ERR_NEG_INT:
            # Handled with an explicit comparison
            error_kind = ERR_NEVER
        target = self.add(CallC(desc.c_function_name, coerced, desc.return_type, desc.steals,
                                desc.is_borrowed, error_kind, line, var_arg_idx))
        if desc.error_kind == ERR_NEG_INT:
            comp = ComparisonOp(target,
                                Integer(0, desc.return_type, line),
                                ComparisonOp.SGE,
                                line)
            comp.error_kind = ERR_FALSE
            self.add(comp)

        if desc.truncated_type is None:
            result = target
        else:
            truncate = self.add(Truncate(target, desc.return_type, desc.truncated_type))
            result = truncate
        if result_type and not is_runtime_subtype(result.type, result_type):
            if is_none_rprimitive(result_type):
                # Special case None return. The actual result may actually be a bool
                # and so we can't just coerce it.
                result = self.none()
            else:
                result = self.coerce(target, result_type, line)
        return result

    def matching_call_c(self,
                        candidates: List[CFunctionDescription],
                        args: List[Value],
                        line: int,
                        result_type: Optional[RType] = None) -> Optional[Value]:
        matching: Optional[CFunctionDescription] = None
        for desc in candidates:
            if len(desc.arg_types) != len(args):
                continue
            if all(is_subtype(actual.type, formal)
                   for actual, formal in zip(args, desc.arg_types)):
                if matching:
                    assert matching.priority != desc.priority, 'Ambiguous:\n1) %s\n2) %s' % (
                        matching, desc)
                    if desc.priority > matching.priority:
                        matching = desc
                else:
                    matching = desc
        if matching:
            target = self.call_c(matching, args, line, result_type)
            return target
        return None

    def int_op(self, type: RType, lhs: Value, rhs: Value, op: int, line: int) -> Value:
        return self.add(IntOp(type, lhs, rhs, op, line))

    def comparison_op(self, lhs: Value, rhs: Value, op: int, line: int) -> Value:
        return self.add(ComparisonOp(lhs, rhs, op, line))

    def builtin_len(self, val: Value, line: int, use_pyssize_t: bool = False) -> Value:
        """Generate len(val).

        Return short_int_rprimitive by default.
        Return c_pyssize_t if use_pyssize_t is true (unshifted).
        """
        typ = val.type
        size_value = None
        if (is_list_rprimitive(typ) or is_tuple_rprimitive(typ)
                or is_bytes_rprimitive(typ)):
            elem_address = self.add(GetElementPtr(val, PyVarObject, 'ob_size'))
            size_value = self.add(LoadMem(c_pyssize_t_rprimitive, elem_address))
            self.add(KeepAlive([val]))
        elif is_set_rprimitive(typ):
            elem_address = self.add(GetElementPtr(val, PySetObject, 'used'))
            size_value = self.add(LoadMem(c_pyssize_t_rprimitive, elem_address))
            self.add(KeepAlive([val]))
        elif is_dict_rprimitive(typ):
            size_value = self.call_c(dict_ssize_t_size_op, [val], line)
        elif is_str_rprimitive(typ):
            size_value = self.call_c(str_ssize_t_size_op, [val], line)

        if size_value is not None:
            if use_pyssize_t:
                return size_value
            offset = Integer(1, c_pyssize_t_rprimitive, line)
            return self.int_op(short_int_rprimitive, size_value, offset,
                               IntOp.LEFT_SHIFT, line)

        if isinstance(typ, RInstance):
            # TODO: Support use_pyssize_t
            assert not use_pyssize_t
            length = self.gen_method_call(val, '__len__', [], int_rprimitive, line)
            length = self.coerce(length, int_rprimitive, line)
            ok, fail = BasicBlock(), BasicBlock()
            self.compare_tagged_condition(length, Integer(0), '>=', ok, fail, line)
            self.activate_block(fail)
            self.add(RaiseStandardError(RaiseStandardError.VALUE_ERROR,
                                        "__len__() should return >= 0",
                                        line))
            self.add(Unreachable())
            self.activate_block(ok)
            return length

        # generic case
        if use_pyssize_t:
            return self.call_c(generic_ssize_t_len_op, [val], line)
        else:
            return self.call_c(generic_len_op, [val], line)

    def new_tuple(self, items: List[Value], line: int) -> Value:
        size: Value = Integer(len(items), c_pyssize_t_rprimitive)
        return self.call_c(new_tuple_op, [size] + items, line)

    def new_tuple_with_length(self, length: Value, line: int) -> Value:
        """This function returns an uninitialized tuple.

        If the length is non-zero, the caller must initialize the tuple, before
        it can be made visible to user code -- otherwise the tuple object is broken.
        You might need further initialization with `new_tuple_set_item_op` op.

        Args:
            length: desired length of the new tuple. The rtype should be
                    c_pyssize_t_rprimitive
            line: line number
        """
        return self.call_c(new_tuple_with_length_op, [length], line)

    # Internal helpers

    def decompose_union_helper(self,
                               obj: Value,
                               rtype: RUnion,
                               result_type: RType,
                               process_item: Callable[[Value], Value],
                               line: int) -> Value:
        """Generate isinstance() + specialized operations for union items.

        Say, for Union[A, B] generate ops resembling this (pseudocode):

            if isinstance(obj, A):
                result = <result of process_item(cast(A, obj)>
            else:
                result = <result of process_item(cast(B, obj)>

        Args:
            obj: value with a union type
            rtype: the union type
            result_type: result of the operation
            process_item: callback to generate op for a single union item (arg is coerced
                to union item type)
            line: line number
        """
        # TODO: Optimize cases where a single operation can handle multiple union items
        #     (say a method is implemented in a common base class)
        fast_items = []
        rest_items = []
        for item in rtype.items:
            if isinstance(item, RInstance):
                fast_items.append(item)
            else:
                # For everything but RInstance we fall back to C API
                rest_items.append(item)
        exit_block = BasicBlock()
        result = Register(result_type)
        for i, item in enumerate(fast_items):
            more_types = i < len(fast_items) - 1 or rest_items
            if more_types:
                # We are not at the final item so we need one more branch
                op = self.isinstance_native(obj, item.class_ir, line)
                true_block, false_block = BasicBlock(), BasicBlock()
                self.add_bool_branch(op, true_block, false_block)
                self.activate_block(true_block)
            coerced = self.coerce(obj, item, line)
            temp = process_item(coerced)
            temp2 = self.coerce(temp, result_type, line)
            self.add(Assign(result, temp2))
            self.goto(exit_block)
            if more_types:
                self.activate_block(false_block)
        if rest_items:
            # For everything else we use generic operation. Use force=True to drop the
            # union type.
            coerced = self.coerce(obj, object_rprimitive, line, force=True)
            temp = process_item(coerced)
            temp2 = self.coerce(temp, result_type, line)
            self.add(Assign(result, temp2))
            self.goto(exit_block)
        self.activate_block(exit_block)
        return result

    def translate_special_method_call(self,
                                      base_reg: Value,
                                      name: str,
                                      args: List[Value],
                                      result_type: Optional[RType],
                                      line: int) -> Optional[Value]:
        """Translate a method call which is handled nongenerically.

        These are special in the sense that we have code generated specifically for them.
        They tend to be method calls which have equivalents in C that are more direct
        than calling with the PyObject api.

        Return None if no translation found; otherwise return the target register.
        """
        call_c_ops_candidates = method_call_ops.get(name, [])
        call_c_op = self.matching_call_c(call_c_ops_candidates, [base_reg] + args,
                                         line, result_type)
        return call_c_op

    def translate_eq_cmp(self,
                         lreg: Value,
                         rreg: Value,
                         expr_op: str,
                         line: int) -> Optional[Value]:
        """Add a equality comparison operation.

        Args:
            expr_op: either '==' or '!='
        """
        ltype = lreg.type
        rtype = rreg.type
        if not (isinstance(ltype, RInstance) and ltype == rtype):
            return None

        class_ir = ltype.class_ir
        # Check whether any subclasses of the operand redefines __eq__
        # or it might be redefined in a Python parent class or by
        # dataclasses
        cmp_varies_at_runtime = (
            not class_ir.is_method_final('__eq__')
            or not class_ir.is_method_final('__ne__')
            or class_ir.inherits_python
            or class_ir.is_augmented
        )

        if cmp_varies_at_runtime:
            # We might need to call left.__eq__(right) or right.__eq__(left)
            # depending on which is the more specific type.
            return None

        if not class_ir.has_method('__eq__'):
            # There's no __eq__ defined, so just use object identity.
            identity_ref_op = 'is' if expr_op == '==' else 'is not'
            return self.translate_is_op(lreg, rreg, identity_ref_op, line)

        return self.gen_method_call(
            lreg,
            op_methods[expr_op],
            [rreg],
            ltype,
            line
        )

    def translate_is_op(self,
                        lreg: Value,
                        rreg: Value,
                        expr_op: str,
                        line: int) -> Value:
        """Create equality comparison operation between object identities

        Args:
            expr_op: either 'is' or 'is not'
        """
        op = ComparisonOp.EQ if expr_op == 'is' else ComparisonOp.NEQ
        lhs = self.coerce(lreg, object_rprimitive, line)
        rhs = self.coerce(rreg, object_rprimitive, line)
        return self.add(ComparisonOp(lhs, rhs, op, line))

    def _create_dict(self,
                     keys: List[Value],
                     values: List[Value],
                     line: int) -> Value:
        """Create a dictionary(possibly empty) using keys and values"""
        # keys and values should have the same number of items
        size = len(keys)
        if size > 0:
            size_value: Value = Integer(size, c_pyssize_t_rprimitive)
            # merge keys and values
            items = [i for t in list(zip(keys, values)) for i in t]
            return self.call_c(dict_build_op, [size_value] + items, line)
        else:
            return self.call_c(dict_new_op, [], line)


def num_positional_args(arg_values: List[Value], arg_kinds: Optional[List[ArgKind]]) -> int:
    if arg_kinds is None:
        return len(arg_values)
    num_pos = 0
    for kind in arg_kinds:
        if kind == ARG_POS:
            num_pos += 1
    return num_pos
