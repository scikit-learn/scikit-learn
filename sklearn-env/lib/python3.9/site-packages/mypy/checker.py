"""Mypy type checker."""

import itertools
import fnmatch
from contextlib import contextmanager

from typing import (
    Any, Dict, Set, List, cast, Tuple, TypeVar, Union, Optional, NamedTuple, Iterator,
    Iterable, Sequence, Mapping, Generic, AbstractSet, Callable
)
from typing_extensions import Final, TypeAlias as _TypeAlias

from mypy.backports import nullcontext
from mypy.errors import Errors, report_internal_error
from mypy.nodes import (
    SymbolTable, Statement, MypyFile, Var, Expression, Lvalue, Node,
    OverloadedFuncDef, FuncDef, FuncItem, FuncBase, TypeInfo,
    ClassDef, Block, AssignmentStmt, NameExpr, MemberExpr, IndexExpr,
    TupleExpr, ListExpr, ExpressionStmt, ReturnStmt, IfStmt,
    WhileStmt, OperatorAssignmentStmt, WithStmt, AssertStmt,
    RaiseStmt, TryStmt, ForStmt, DelStmt, CallExpr, IntExpr, StrExpr,
    UnicodeExpr, OpExpr, UnaryExpr, LambdaExpr, TempNode, SymbolTableNode,
    Context, Decorator, PrintStmt, BreakStmt, PassStmt, ContinueStmt,
    ComparisonExpr, StarExpr, EllipsisExpr, RefExpr, PromoteExpr,
    Import, ImportFrom, ImportAll, ImportBase, TypeAlias,
    ARG_POS, ARG_STAR, LITERAL_TYPE, LDEF, MDEF, GDEF,
    CONTRAVARIANT, COVARIANT, INVARIANT, TypeVarExpr, AssignmentExpr,
    is_final_node,
    ARG_NAMED)
from mypy import nodes
from mypy import operators
from mypy.literals import literal, literal_hash, Key
from mypy.typeanal import has_any_from_unimported_type, check_for_explicit_any, make_optional_type
from mypy.types import (
    Type, AnyType, CallableType, FunctionLike, Overloaded, TupleType, TypedDictType,
    Instance, NoneType, strip_type, TypeType, TypeOfAny,
    UnionType, TypeVarId, TypeVarType, PartialType, DeletedType, UninhabitedType,
    is_named_instance, union_items, TypeQuery, LiteralType,
    is_optional, remove_optional, TypeTranslator, StarType, get_proper_type, ProperType,
    get_proper_types, is_literal_type, TypeAliasType, TypeGuardedType, ParamSpecType
)
from mypy.sametypes import is_same_type
from mypy.messages import (
    MessageBuilder, make_inferred_type_note, append_invariance_notes, pretty_seq,
    format_type, format_type_bare, format_type_distinctly, SUGGESTED_TEST_FIXTURES
)
import mypy.checkexpr
from mypy.checkmember import (
    MemberContext, analyze_member_access, analyze_descriptor_access,
    type_object_type,
    analyze_decorator_or_funcbase_access,
)
from mypy.semanal_enum import ENUM_BASES, ENUM_SPECIAL_PROPS
from mypy.typeops import (
    map_type_from_supertype, bind_self, erase_to_bound, make_simplified_union,
    erase_def_to_union_or_bound, erase_to_union_or_bound, coerce_to_literal,
    try_getting_str_literals_from_type, try_getting_int_literals_from_type,
    tuple_fallback, is_singleton_type, try_expanding_sum_type_to_union,
    true_only, false_only, function_type, get_type_vars, custom_special_method,
    is_literal_type_like,
)
from mypy import message_registry
from mypy.message_registry import ErrorMessage
from mypy.subtypes import (
    is_subtype, is_equivalent, is_proper_subtype, is_more_precise,
    restrict_subtype_away, is_subtype_ignoring_tvars, is_callable_compatible,
    unify_generic_callable, find_member
)
from mypy.constraints import SUPERTYPE_OF
from mypy.maptype import map_instance_to_supertype
from mypy.typevars import fill_typevars, has_no_typevars, fill_typevars_with_any
from mypy.semanal import set_callable_name, refers_to_fullname
from mypy.mro import calculate_mro, MroError
from mypy.erasetype import erase_typevars, remove_instance_last_known_values, erase_type
from mypy.expandtype import expand_type, expand_type_by_instance
from mypy.visitor import NodeVisitor
from mypy.join import join_types
from mypy.treetransform import TransformVisitor
from mypy.binder import ConditionalTypeBinder, get_declaration
from mypy.meet import is_overlapping_erased_types, is_overlapping_types
from mypy.options import Options
from mypy.plugin import Plugin, CheckerPluginInterface
from mypy.sharedparse import BINARY_MAGIC_METHODS
from mypy.scope import Scope
from mypy import state, errorcodes as codes
from mypy.traverser import has_return_statement, all_return_statements
from mypy.errorcodes import ErrorCode
from mypy.util import is_typeshed_file

T = TypeVar('T')

DEFAULT_LAST_PASS: Final = 1  # Pass numbers start at 0

DeferredNodeType: _TypeAlias = Union[FuncDef, LambdaExpr, OverloadedFuncDef, Decorator]
FineGrainedDeferredNodeType: _TypeAlias = Union[FuncDef, MypyFile, OverloadedFuncDef]

# A node which is postponed to be processed during the next pass.
# In normal mode one can defer functions and methods (also decorated and/or overloaded)
# and lambda expressions. Nested functions can't be deferred -- only top-level functions
# and methods of classes not defined within a function can be deferred.
DeferredNode = NamedTuple(
    'DeferredNode',
    [
        ('node', DeferredNodeType),
        ('active_typeinfo', Optional[TypeInfo]),  # And its TypeInfo (for semantic analysis
                                                  # self type handling)
    ])

# Same as above, but for fine-grained mode targets. Only top-level functions/methods
# and module top levels are allowed as such.
FineGrainedDeferredNode = NamedTuple(
    'FineGrainedDeferredNode',
    [
        ('node', FineGrainedDeferredNodeType),
        ('active_typeinfo', Optional[TypeInfo]),
    ])

# Data structure returned by find_isinstance_check representing
# information learned from the truth or falsehood of a condition.  The
# dict maps nodes representing expressions like 'a[0].x' to their
# refined types under the assumption that the condition has a
# particular truth value. A value of None means that the condition can
# never have that truth value.

# NB: The keys of this dict are nodes in the original source program,
# which are compared by reference equality--effectively, being *the
# same* expression of the program, not just two identical expressions
# (such as two references to the same variable). TODO: it would
# probably be better to have the dict keyed by the nodes' literal_hash
# field instead.

TypeMap: _TypeAlias = Optional[Dict[Expression, Type]]

# An object that represents either a precise type or a type with an upper bound;
# it is important for correct type inference with isinstance.
TypeRange = NamedTuple(
    'TypeRange',
    [
        ('item', Type),
        ('is_upper_bound', bool),  # False => precise type
    ])

# Keeps track of partial types in a single scope. In fine-grained incremental
# mode partial types initially defined at the top level cannot be completed in
# a function, and we use the 'is_function' attribute to enforce this.
PartialTypeScope = NamedTuple('PartialTypeScope', [('map', Dict[Var, Context]),
                                                   ('is_function', bool),
                                                   ('is_local', bool),
                                                   ])


class TypeChecker(NodeVisitor[None], CheckerPluginInterface):
    """Mypy type checker.

    Type check mypy source files that have been semantically analyzed.

    You must create a separate instance for each source file.
    """

    # Are we type checking a stub?
    is_stub = False
    # Error message reporter
    errors: Errors
    # Utility for generating messages
    msg: MessageBuilder
    # Types of type checked nodes
    type_map: Dict[Expression, Type]

    # Helper for managing conditional types
    binder: ConditionalTypeBinder
    # Helper for type checking expressions
    expr_checker: mypy.checkexpr.ExpressionChecker

    tscope: Scope
    scope: "CheckerScope"
    # Stack of function return types
    return_types: List[Type]
    # Flags; true for dynamically typed functions
    dynamic_funcs: List[bool]
    # Stack of collections of variables with partial types
    partial_types: List[PartialTypeScope]
    # Vars for which partial type errors are already reported
    # (to avoid logically duplicate errors with different error context).
    partial_reported: Set[Var]
    globals: SymbolTable
    modules: Dict[str, MypyFile]
    # Nodes that couldn't be checked because some types weren't available. We'll run
    # another pass and try these again.
    deferred_nodes: List[DeferredNode]
    # Type checking pass number (0 = first pass)
    pass_num = 0
    # Last pass number to take
    last_pass = DEFAULT_LAST_PASS
    # Have we deferred the current function? If yes, don't infer additional
    # types during this pass within the function.
    current_node_deferred = False
    # Is this file a typeshed stub?
    is_typeshed_stub = False
    # Should strict Optional-related errors be suppressed in this file?
    suppress_none_errors = False  # TODO: Get it from options instead
    options: Options
    # Used for collecting inferred attribute types so that they can be checked
    # for consistency.
    inferred_attribute_types: Optional[Dict[Var, Type]] = None
    # Don't infer partial None types if we are processing assignment from Union
    no_partial_types: bool = False

    # The set of all dependencies (suppressed or not) that this module accesses, either
    # directly or indirectly.
    module_refs: Set[str]

    # A map from variable nodes to a snapshot of the frame ids of the
    # frames that were active when the variable was declared. This can
    # be used to determine nearest common ancestor frame of a variable's
    # declaration and the current frame, which lets us determine if it
    # was declared in a different branch of the same `if` statement
    # (if that frame is a conditional_frame).
    var_decl_frames: Dict[Var, Set[int]]

    # Plugin that provides special type checking rules for specific library
    # functions such as open(), etc.
    plugin: Plugin

    def __init__(self, errors: Errors, modules: Dict[str, MypyFile], options: Options,
                 tree: MypyFile, path: str, plugin: Plugin) -> None:
        """Construct a type checker.

        Use errors to report type check errors.
        """
        self.errors = errors
        self.modules = modules
        self.options = options
        self.tree = tree
        self.path = path
        self.msg = MessageBuilder(errors, modules)
        self.plugin = plugin
        self.expr_checker = mypy.checkexpr.ExpressionChecker(self, self.msg, self.plugin)
        self.tscope = Scope()
        self.scope = CheckerScope(tree)
        self.binder = ConditionalTypeBinder()
        self.globals = tree.names
        self.return_types = []
        self.dynamic_funcs = []
        self.partial_types = []
        self.partial_reported = set()
        self.var_decl_frames = {}
        self.deferred_nodes = []
        self.type_map = {}
        self.module_refs = set()
        self.pass_num = 0
        self.current_node_deferred = False
        self.is_stub = tree.is_stub
        self.is_typeshed_stub = is_typeshed_file(path)
        self.inferred_attribute_types = None
        if options.strict_optional_whitelist is None:
            self.suppress_none_errors = not options.show_none_errors
        else:
            self.suppress_none_errors = not any(fnmatch.fnmatch(path, pattern)
                                                for pattern
                                                in options.strict_optional_whitelist)
        # If True, process function definitions. If False, don't. This is used
        # for processing module top levels in fine-grained incremental mode.
        self.recurse_into_functions = True
        # This internal flag is used to track whether we a currently type-checking
        # a final declaration (assignment), so that some errors should be suppressed.
        # Should not be set manually, use get_final_context/enter_final_context instead.
        # NOTE: we use the context manager to avoid "threading" an additional `is_final_def`
        # argument through various `checker` and `checkmember` functions.
        self._is_final_def = False

    @property
    def type_context(self) -> List[Optional[Type]]:
        return self.expr_checker.type_context

    def reset(self) -> None:
        """Cleanup stale state that might be left over from a typechecking run.

        This allows us to reuse TypeChecker objects in fine-grained
        incremental mode.
        """
        # TODO: verify this is still actually worth it over creating new checkers
        self.partial_reported.clear()
        self.module_refs.clear()
        self.binder = ConditionalTypeBinder()
        self.type_map.clear()

        assert self.inferred_attribute_types is None
        assert self.partial_types == []
        assert self.deferred_nodes == []
        assert len(self.scope.stack) == 1
        assert self.partial_types == []

    def check_first_pass(self) -> None:
        """Type check the entire file, but defer functions with unresolved references.

        Unresolved references are forward references to variables
        whose types haven't been inferred yet.  They may occur later
        in the same file or in a different file that's being processed
        later (usually due to an import cycle).

        Deferred functions will be processed by check_second_pass().
        """
        self.recurse_into_functions = True
        with state.strict_optional_set(self.options.strict_optional):
            self.errors.set_file(self.path, self.tree.fullname, scope=self.tscope)
            with self.tscope.module_scope(self.tree.fullname):
                with self.enter_partial_types(), self.binder.top_frame_context():
                    for d in self.tree.defs:
                        if (self.binder.is_unreachable()
                                and self.should_report_unreachable_issues()
                                and not self.is_raising_or_empty(d)):
                            self.msg.unreachable_statement(d)
                            break
                        self.accept(d)

                assert not self.current_node_deferred

                all_ = self.globals.get('__all__')
                if all_ is not None and all_.type is not None:
                    all_node = all_.node
                    assert all_node is not None
                    seq_str = self.named_generic_type('typing.Sequence',
                                                      [self.named_type('builtins.str')])
                    if self.options.python_version[0] < 3:
                        seq_str = self.named_generic_type('typing.Sequence',
                                                          [self.named_type('builtins.unicode')])
                    if not is_subtype(all_.type, seq_str):
                        str_seq_s, all_s = format_type_distinctly(seq_str, all_.type)
                        self.fail(message_registry.ALL_MUST_BE_SEQ_STR.format(str_seq_s, all_s),
                                  all_node)

    def check_second_pass(self,
                          todo: Optional[Sequence[Union[DeferredNode,
                                                        FineGrainedDeferredNode]]] = None
                          ) -> bool:
        """Run second or following pass of type checking.

        This goes through deferred nodes, returning True if there were any.
        """
        self.recurse_into_functions = True
        with state.strict_optional_set(self.options.strict_optional):
            if not todo and not self.deferred_nodes:
                return False
            self.errors.set_file(self.path, self.tree.fullname, scope=self.tscope)
            with self.tscope.module_scope(self.tree.fullname):
                self.pass_num += 1
                if not todo:
                    todo = self.deferred_nodes
                else:
                    assert not self.deferred_nodes
                self.deferred_nodes = []
                done: Set[Union[DeferredNodeType, FineGrainedDeferredNodeType]] = set()
                for node, active_typeinfo in todo:
                    if node in done:
                        continue
                    # This is useful for debugging:
                    # print("XXX in pass %d, class %s, function %s" %
                    #       (self.pass_num, type_name, node.fullname or node.name))
                    done.add(node)
                    with self.tscope.class_scope(active_typeinfo) if active_typeinfo \
                            else nullcontext():
                        with self.scope.push_class(active_typeinfo) if active_typeinfo \
                                else nullcontext():
                            self.check_partial(node)
            return True

    def check_partial(self, node: Union[DeferredNodeType, FineGrainedDeferredNodeType]) -> None:
        if isinstance(node, MypyFile):
            self.check_top_level(node)
        else:
            self.recurse_into_functions = True
            if isinstance(node, LambdaExpr):
                self.expr_checker.accept(node)
            else:
                self.accept(node)

    def check_top_level(self, node: MypyFile) -> None:
        """Check only the top-level of a module, skipping function definitions."""
        self.recurse_into_functions = False
        with self.enter_partial_types():
            with self.binder.top_frame_context():
                for d in node.defs:
                    d.accept(self)

        assert not self.current_node_deferred
        # TODO: Handle __all__

    def defer_node(self, node: DeferredNodeType, enclosing_class: Optional[TypeInfo]) -> None:
        """Defer a node for processing during next type-checking pass.

        Args:
            node: function/method being deferred
            enclosing_class: for methods, the class where the method is defined
        NOTE: this can't handle nested functions/methods.
        """
        # We don't freeze the entire scope since only top-level functions and methods
        # can be deferred. Only module/class level scope information is needed.
        # Module-level scope information is preserved in the TypeChecker instance.
        self.deferred_nodes.append(DeferredNode(node, enclosing_class))

    def handle_cannot_determine_type(self, name: str, context: Context) -> None:
        node = self.scope.top_non_lambda_function()
        if self.pass_num < self.last_pass and isinstance(node, FuncDef):
            # Don't report an error yet. Just defer. Note that we don't defer
            # lambdas because they are coupled to the surrounding function
            # through the binder and the inferred type of the lambda, so it
            # would get messy.
            enclosing_class = self.scope.enclosing_class()
            self.defer_node(node, enclosing_class)
            # Set a marker so that we won't infer additional types in this
            # function. Any inferred types could be bogus, because there's at
            # least one type that we don't know.
            self.current_node_deferred = True
        else:
            self.msg.cannot_determine_type(name, context)

    def accept(self, stmt: Statement) -> None:
        """Type check a node in the given type context."""
        try:
            stmt.accept(self)
        except Exception as err:
            report_internal_error(err, self.errors.file, stmt.line, self.errors, self.options)

    def accept_loop(self, body: Statement, else_body: Optional[Statement] = None, *,
                    exit_condition: Optional[Expression] = None) -> None:
        """Repeatedly type check a loop body until the frame doesn't change.
        If exit_condition is set, assume it must be False on exit from the loop.

        Then check the else_body.
        """
        # The outer frame accumulates the results of all iterations
        with self.binder.frame_context(can_skip=False, conditional_frame=True):
            while True:
                with self.binder.frame_context(can_skip=True,
                                               break_frame=2, continue_frame=1):
                    self.accept(body)
                if not self.binder.last_pop_changed:
                    break
            if exit_condition:
                _, else_map = self.find_isinstance_check(exit_condition)
                self.push_type_map(else_map)
            if else_body:
                self.accept(else_body)

    #
    # Definitions
    #

    def visit_overloaded_func_def(self, defn: OverloadedFuncDef) -> None:
        if not self.recurse_into_functions:
            return
        with self.tscope.function_scope(defn):
            self._visit_overloaded_func_def(defn)

    def _visit_overloaded_func_def(self, defn: OverloadedFuncDef) -> None:
        num_abstract = 0
        if not defn.items:
            # In this case we have already complained about none of these being
            # valid overloads.
            return None
        if len(defn.items) == 1:
            self.fail(message_registry.MULTIPLE_OVERLOADS_REQUIRED, defn)

        if defn.is_property:
            # HACK: Infer the type of the property.
            self.visit_decorator(cast(Decorator, defn.items[0]))
        for fdef in defn.items:
            assert isinstance(fdef, Decorator)
            self.check_func_item(fdef.func, name=fdef.func.name)
            if fdef.func.is_abstract:
                num_abstract += 1
        if num_abstract not in (0, len(defn.items)):
            self.fail(message_registry.INCONSISTENT_ABSTRACT_OVERLOAD, defn)
        if defn.impl:
            defn.impl.accept(self)
        if defn.info:
            self.check_method_override(defn)
            self.check_inplace_operator_method(defn)
        if not defn.is_property:
            self.check_overlapping_overloads(defn)
        return None

    def check_overlapping_overloads(self, defn: OverloadedFuncDef) -> None:
        # At this point we should have set the impl already, and all remaining
        # items are decorators

        if self.msg.errors.file in self.msg.errors.ignored_files:
            # This is a little hacky, however, the quadratic check here is really expensive, this
            # method has no side effects, so we should skip it if we aren't going to report
            # anything. In some other places we swallow errors in stubs, but this error is very
            # useful for stubs!
            return

        # Compute some info about the implementation (if it exists) for use below
        impl_type: Optional[CallableType] = None
        if defn.impl:
            if isinstance(defn.impl, FuncDef):
                inner_type: Optional[Type] = defn.impl.type
            elif isinstance(defn.impl, Decorator):
                inner_type = defn.impl.var.type
            else:
                assert False, "Impl isn't the right type"

            # This can happen if we've got an overload with a different
            # decorator or if the implementation is untyped -- we gave up on the types.
            inner_type = get_proper_type(inner_type)
            if inner_type is not None and not isinstance(inner_type, AnyType):
                if isinstance(inner_type, CallableType):
                    impl_type = inner_type
                elif isinstance(inner_type, Instance):
                    inner_call = get_proper_type(
                        analyze_member_access(
                            name='__call__',
                            typ=inner_type,
                            context=defn.impl,
                            is_lvalue=False,
                            is_super=False,
                            is_operator=True,
                            msg=self.msg,
                            original_type=inner_type,
                            chk=self,
                        ),
                    )
                    if isinstance(inner_call, CallableType):
                        impl_type = inner_call
                if impl_type is None:
                    self.msg.not_callable(inner_type, defn.impl)

        is_descriptor_get = defn.info and defn.name == "__get__"
        for i, item in enumerate(defn.items):
            # TODO overloads involving decorators
            assert isinstance(item, Decorator)
            sig1 = self.function_type(item.func)
            assert isinstance(sig1, CallableType)

            for j, item2 in enumerate(defn.items[i + 1:]):
                assert isinstance(item2, Decorator)
                sig2 = self.function_type(item2.func)
                assert isinstance(sig2, CallableType)

                if not are_argument_counts_overlapping(sig1, sig2):
                    continue

                if overload_can_never_match(sig1, sig2):
                    self.msg.overloaded_signature_will_never_match(
                        i + 1, i + j + 2, item2.func)
                elif not is_descriptor_get:
                    # Note: we force mypy to check overload signatures in strict-optional mode
                    # so we don't incorrectly report errors when a user tries typing an overload
                    # that happens to have a 'if the argument is None' fallback.
                    #
                    # For example, the following is fine in strict-optional mode but would throw
                    # the unsafe overlap error when strict-optional is disabled:
                    #
                    #     @overload
                    #     def foo(x: None) -> int: ...
                    #     @overload
                    #     def foo(x: str) -> str: ...
                    #
                    # See Python 2's map function for a concrete example of this kind of overload.
                    with state.strict_optional_set(True):
                        if is_unsafe_overlapping_overload_signatures(sig1, sig2):
                            self.msg.overloaded_signatures_overlap(
                                i + 1, i + j + 2, item.func)

            if impl_type is not None:
                assert defn.impl is not None

                # We perform a unification step that's very similar to what
                # 'is_callable_compatible' would have done if we had set
                # 'unify_generics' to True -- the only difference is that
                # we check and see if the impl_type's return value is a
                # *supertype* of the overload alternative, not a *subtype*.
                #
                # This is to match the direction the implementation's return
                # needs to be compatible in.
                if impl_type.variables:
                    impl = unify_generic_callable(impl_type, sig1,
                                                  ignore_return=False,
                                                  return_constraint_direction=SUPERTYPE_OF)
                    if impl is None:
                        self.msg.overloaded_signatures_typevar_specific(i + 1, defn.impl)
                        continue
                else:
                    impl = impl_type

                # Prevent extra noise from inconsistent use of @classmethod by copying
                # the first arg from the method being checked against.
                if sig1.arg_types and defn.info:
                    impl = impl.copy_modified(arg_types=[sig1.arg_types[0]] + impl.arg_types[1:])

                # Is the overload alternative's arguments subtypes of the implementation's?
                if not is_callable_compatible(impl, sig1,
                                              is_compat=is_subtype_no_promote,
                                              ignore_return=True):
                    self.msg.overloaded_signatures_arg_specific(i + 1, defn.impl)

                # Is the overload alternative's return type a subtype of the implementation's?
                if not is_subtype_no_promote(sig1.ret_type, impl.ret_type):
                    self.msg.overloaded_signatures_ret_specific(i + 1, defn.impl)

    # Here's the scoop about generators and coroutines.
    #
    # There are two kinds of generators: classic generators (functions
    # with `yield` or `yield from` in the body) and coroutines
    # (functions declared with `async def`).  The latter are specified
    # in PEP 492 and only available in Python >= 3.5.
    #
    # Classic generators can be parameterized with three types:
    # - ty is the Yield type (the type of y in `yield y`)
    # - tc is the type reCeived by yield (the type of c in `c = yield`).
    # - tr is the Return type (the type of r in `return r`)
    #
    # A classic generator must define a return type that's either
    # `Generator[ty, tc, tr]`, Iterator[ty], or Iterable[ty] (or
    # object or Any).  If tc/tr are not given, both are None.
    #
    # A coroutine must define a return type corresponding to tr; the
    # other two are unconstrained.  The "external" return type (seen
    # by the caller) is Awaitable[tr].
    #
    # In addition, there's the synthetic type AwaitableGenerator: it
    # inherits from both Awaitable and Generator and can be used both
    # in `yield from` and in `await`.  This type is set automatically
    # for functions decorated with `@types.coroutine` or
    # `@asyncio.coroutine`.  Its single parameter corresponds to tr.
    #
    # PEP 525 adds a new type, the asynchronous generator, which was
    # first released in Python 3.6. Async generators are `async def`
    # functions that can also `yield` values. They can be parameterized
    # with two types, ty and tc, because they cannot return a value.
    #
    # There are several useful methods, each taking a type t and a
    # flag c indicating whether it's for a generator or coroutine:
    #
    # - is_generator_return_type(t, c) returns whether t is a Generator,
    #   Iterator, Iterable (if not c), or Awaitable (if c), or
    #   AwaitableGenerator (regardless of c).
    # - is_async_generator_return_type(t) returns whether t is an
    #   AsyncGenerator.
    # - get_generator_yield_type(t, c) returns ty.
    # - get_generator_receive_type(t, c) returns tc.
    # - get_generator_return_type(t, c) returns tr.

    def is_generator_return_type(self, typ: Type, is_coroutine: bool) -> bool:
        """Is `typ` a valid type for a generator/coroutine?

        True if `typ` is a *supertype* of Generator or Awaitable.
        Also true it it's *exactly* AwaitableGenerator (modulo type parameters).
        """
        typ = get_proper_type(typ)
        if is_coroutine:
            # This means we're in Python 3.5 or later.
            at = self.named_generic_type('typing.Awaitable', [AnyType(TypeOfAny.special_form)])
            if is_subtype(at, typ):
                return True
        else:
            any_type = AnyType(TypeOfAny.special_form)
            gt = self.named_generic_type('typing.Generator', [any_type, any_type, any_type])
            if is_subtype(gt, typ):
                return True
        return isinstance(typ, Instance) and typ.type.fullname == 'typing.AwaitableGenerator'

    def is_async_generator_return_type(self, typ: Type) -> bool:
        """Is `typ` a valid type for an async generator?

        True if `typ` is a supertype of AsyncGenerator.
        """
        try:
            any_type = AnyType(TypeOfAny.special_form)
            agt = self.named_generic_type('typing.AsyncGenerator', [any_type, any_type])
        except KeyError:
            # we're running on a version of typing that doesn't have AsyncGenerator yet
            return False
        return is_subtype(agt, typ)

    def get_generator_yield_type(self, return_type: Type, is_coroutine: bool) -> Type:
        """Given the declared return type of a generator (t), return the type it yields (ty)."""
        return_type = get_proper_type(return_type)

        if isinstance(return_type, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=return_type)
        elif (not self.is_generator_return_type(return_type, is_coroutine)
                and not self.is_async_generator_return_type(return_type)):
            # If the function doesn't have a proper Generator (or
            # Awaitable) return type, anything is permissible.
            return AnyType(TypeOfAny.from_error)
        elif not isinstance(return_type, Instance):
            # Same as above, but written as a separate branch so the typechecker can understand.
            return AnyType(TypeOfAny.from_error)
        elif return_type.type.fullname == 'typing.Awaitable':
            # Awaitable: ty is Any.
            return AnyType(TypeOfAny.special_form)
        elif return_type.args:
            # AwaitableGenerator, Generator, AsyncGenerator, Iterator, or Iterable; ty is args[0].
            ret_type = return_type.args[0]
            # TODO not best fix, better have dedicated yield token
            return ret_type
        else:
            # If the function's declared supertype of Generator has no type
            # parameters (i.e. is `object`), then the yielded values can't
            # be accessed so any type is acceptable.  IOW, ty is Any.
            # (However, see https://github.com/python/mypy/issues/1933)
            return AnyType(TypeOfAny.special_form)

    def get_generator_receive_type(self, return_type: Type, is_coroutine: bool) -> Type:
        """Given a declared generator return type (t), return the type its yield receives (tc)."""
        return_type = get_proper_type(return_type)

        if isinstance(return_type, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=return_type)
        elif (not self.is_generator_return_type(return_type, is_coroutine)
                and not self.is_async_generator_return_type(return_type)):
            # If the function doesn't have a proper Generator (or
            # Awaitable) return type, anything is permissible.
            return AnyType(TypeOfAny.from_error)
        elif not isinstance(return_type, Instance):
            # Same as above, but written as a separate branch so the typechecker can understand.
            return AnyType(TypeOfAny.from_error)
        elif return_type.type.fullname == 'typing.Awaitable':
            # Awaitable, AwaitableGenerator: tc is Any.
            return AnyType(TypeOfAny.special_form)
        elif (return_type.type.fullname in ('typing.Generator', 'typing.AwaitableGenerator')
              and len(return_type.args) >= 3):
            # Generator: tc is args[1].
            return return_type.args[1]
        elif return_type.type.fullname == 'typing.AsyncGenerator' and len(return_type.args) >= 2:
            return return_type.args[1]
        else:
            # `return_type` is a supertype of Generator, so callers won't be able to send it
            # values.  IOW, tc is None.
            return NoneType()

    def get_coroutine_return_type(self, return_type: Type) -> Type:
        return_type = get_proper_type(return_type)
        if isinstance(return_type, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=return_type)
        assert isinstance(return_type, Instance), "Should only be called on coroutine functions."
        # Note: return type is the 3rd type parameter of Coroutine.
        return return_type.args[2]

    def get_generator_return_type(self, return_type: Type, is_coroutine: bool) -> Type:
        """Given the declared return type of a generator (t), return the type it returns (tr)."""
        return_type = get_proper_type(return_type)

        if isinstance(return_type, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=return_type)
        elif not self.is_generator_return_type(return_type, is_coroutine):
            # If the function doesn't have a proper Generator (or
            # Awaitable) return type, anything is permissible.
            return AnyType(TypeOfAny.from_error)
        elif not isinstance(return_type, Instance):
            # Same as above, but written as a separate branch so the typechecker can understand.
            return AnyType(TypeOfAny.from_error)
        elif return_type.type.fullname == 'typing.Awaitable' and len(return_type.args) == 1:
            # Awaitable: tr is args[0].
            return return_type.args[0]
        elif (return_type.type.fullname in ('typing.Generator', 'typing.AwaitableGenerator')
              and len(return_type.args) >= 3):
            # AwaitableGenerator, Generator: tr is args[2].
            return return_type.args[2]
        else:
            # Supertype of Generator (Iterator, Iterable, object): tr is any.
            return AnyType(TypeOfAny.special_form)

    def visit_func_def(self, defn: FuncDef) -> None:
        if not self.recurse_into_functions:
            return
        with self.tscope.function_scope(defn):
            self._visit_func_def(defn)

    def _visit_func_def(self, defn: FuncDef) -> None:
        """Type check a function definition."""
        self.check_func_item(defn, name=defn.name)
        if defn.info:
            if not defn.is_dynamic() and not defn.is_overload and not defn.is_decorated:
                # If the definition is the implementation for an
                # overload, the legality of the override has already
                # been typechecked, and decorated methods will be
                # checked when the decorator is.
                self.check_method_override(defn)
            self.check_inplace_operator_method(defn)
        if defn.original_def:
            # Override previous definition.
            new_type = self.function_type(defn)
            if isinstance(defn.original_def, FuncDef):
                # Function definition overrides function definition.
                if not is_same_type(new_type, self.function_type(defn.original_def)):
                    self.msg.incompatible_conditional_function_def(defn)
            else:
                # Function definition overrides a variable initialized via assignment or a
                # decorated function.
                orig_type = defn.original_def.type
                if orig_type is None:
                    # XXX This can be None, as happens in
                    # test_testcheck_TypeCheckSuite.testRedefinedFunctionInTryWithElse
                    self.msg.note("Internal mypy error checking function redefinition", defn)
                    return
                if isinstance(orig_type, PartialType):
                    if orig_type.type is None:
                        # Ah this is a partial type. Give it the type of the function.
                        orig_def = defn.original_def
                        if isinstance(orig_def, Decorator):
                            var = orig_def.var
                        else:
                            var = orig_def
                        partial_types = self.find_partial_types(var)
                        if partial_types is not None:
                            var.type = new_type
                            del partial_types[var]
                    else:
                        # Trying to redefine something like partial empty list as function.
                        self.fail(message_registry.INCOMPATIBLE_REDEFINITION, defn)
                else:
                    # TODO: Update conditional type binder.
                    self.check_subtype(new_type, orig_type, defn,
                                       message_registry.INCOMPATIBLE_REDEFINITION,
                                       'redefinition with type',
                                       'original type')

    def check_func_item(self, defn: FuncItem,
                        type_override: Optional[CallableType] = None,
                        name: Optional[str] = None) -> None:
        """Type check a function.

        If type_override is provided, use it as the function type.
        """
        self.dynamic_funcs.append(defn.is_dynamic() and not type_override)

        with self.enter_partial_types(is_function=True):
            typ = self.function_type(defn)
            if type_override:
                typ = type_override.copy_modified(line=typ.line, column=typ.column)
            if isinstance(typ, CallableType):
                with self.enter_attribute_inference_context():
                    self.check_func_def(defn, typ, name)
            else:
                raise RuntimeError('Not supported')

        self.dynamic_funcs.pop()
        self.current_node_deferred = False

        if name == '__exit__':
            self.check__exit__return_type(defn)

    @contextmanager
    def enter_attribute_inference_context(self) -> Iterator[None]:
        old_types = self.inferred_attribute_types
        self.inferred_attribute_types = {}
        yield None
        self.inferred_attribute_types = old_types

    def check_func_def(self, defn: FuncItem, typ: CallableType, name: Optional[str]) -> None:
        """Type check a function definition."""
        # Expand type variables with value restrictions to ordinary types.
        expanded = self.expand_typevars(defn, typ)
        for item, typ in expanded:
            old_binder = self.binder
            self.binder = ConditionalTypeBinder()
            with self.binder.top_frame_context():
                defn.expanded.append(item)

                # We may be checking a function definition or an anonymous
                # function. In the first case, set up another reference with the
                # precise type.
                if isinstance(item, FuncDef):
                    fdef = item
                    # Check if __init__ has an invalid, non-None return type.
                    if (fdef.info and fdef.name in ('__init__', '__init_subclass__') and
                            not isinstance(get_proper_type(typ.ret_type), NoneType) and
                            not self.dynamic_funcs[-1]):
                        self.fail(message_registry.MUST_HAVE_NONE_RETURN_TYPE.format(fdef.name),
                                  item)

                    # Check validity of __new__ signature
                    if fdef.info and fdef.name == '__new__':
                        self.check___new___signature(fdef, typ)

                    self.check_for_missing_annotations(fdef)
                    if self.options.disallow_any_unimported:
                        if fdef.type and isinstance(fdef.type, CallableType):
                            ret_type = fdef.type.ret_type
                            if has_any_from_unimported_type(ret_type):
                                self.msg.unimported_type_becomes_any("Return type", ret_type, fdef)
                            for idx, arg_type in enumerate(fdef.type.arg_types):
                                if has_any_from_unimported_type(arg_type):
                                    prefix = "Argument {} to \"{}\"".format(idx + 1, fdef.name)
                                    self.msg.unimported_type_becomes_any(prefix, arg_type, fdef)
                    check_for_explicit_any(fdef.type, self.options, self.is_typeshed_stub,
                                           self.msg, context=fdef)

                if name:  # Special method names
                    if defn.info and self.is_reverse_op_method(name):
                        self.check_reverse_op_method(item, typ, name, defn)
                    elif name in ('__getattr__', '__getattribute__'):
                        self.check_getattr_method(typ, defn, name)
                    elif name == '__setattr__':
                        self.check_setattr_method(typ, defn)

                # Refuse contravariant return type variable
                if isinstance(typ.ret_type, TypeVarType):
                    if typ.ret_type.variance == CONTRAVARIANT:
                        self.fail(message_registry.RETURN_TYPE_CANNOT_BE_CONTRAVARIANT,
                                  typ.ret_type)

                # Check that Generator functions have the appropriate return type.
                if defn.is_generator:
                    if defn.is_async_generator:
                        if not self.is_async_generator_return_type(typ.ret_type):
                            self.fail(message_registry.INVALID_RETURN_TYPE_FOR_ASYNC_GENERATOR,
                                      typ)
                    else:
                        if not self.is_generator_return_type(typ.ret_type, defn.is_coroutine):
                            self.fail(message_registry.INVALID_RETURN_TYPE_FOR_GENERATOR, typ)

                    # Python 2 generators aren't allowed to return values.
                    orig_ret_type = get_proper_type(typ.ret_type)
                    if (self.options.python_version[0] == 2 and
                            isinstance(orig_ret_type, Instance) and
                            orig_ret_type.type.fullname == 'typing.Generator'):
                        if not isinstance(get_proper_type(orig_ret_type.args[2]),
                                          (NoneType, AnyType)):
                            self.fail(message_registry.INVALID_GENERATOR_RETURN_ITEM_TYPE, typ)

                # Fix the type if decorated with `@types.coroutine` or `@asyncio.coroutine`.
                if defn.is_awaitable_coroutine:
                    # Update the return type to AwaitableGenerator.
                    # (This doesn't exist in typing.py, only in typing.pyi.)
                    t = typ.ret_type
                    c = defn.is_coroutine
                    ty = self.get_generator_yield_type(t, c)
                    tc = self.get_generator_receive_type(t, c)
                    if c:
                        tr = self.get_coroutine_return_type(t)
                    else:
                        tr = self.get_generator_return_type(t, c)
                    ret_type = self.named_generic_type('typing.AwaitableGenerator',
                                                       [ty, tc, tr, t])
                    typ = typ.copy_modified(ret_type=ret_type)
                    defn.type = typ

                # Push return type.
                self.return_types.append(typ.ret_type)

                # Store argument types.
                for i in range(len(typ.arg_types)):
                    arg_type = typ.arg_types[i]
                    with self.scope.push_function(defn):
                        # We temporary push the definition to get the self type as
                        # visible from *inside* of this function/method.
                        ref_type: Optional[Type] = self.scope.active_self_type()
                    if (isinstance(defn, FuncDef) and ref_type is not None and i == 0
                            and not defn.is_static
                            and typ.arg_kinds[0] not in [nodes.ARG_STAR, nodes.ARG_STAR2]):
                        isclass = defn.is_class or defn.name in ('__new__', '__init_subclass__')
                        if isclass:
                            ref_type = mypy.types.TypeType.make_normalized(ref_type)
                        erased = get_proper_type(erase_to_bound(arg_type))
                        if not is_subtype_ignoring_tvars(ref_type, erased):
                            note = None
                            if (isinstance(erased, Instance) and erased.type.is_protocol or
                                    isinstance(erased, TypeType) and
                                    isinstance(erased.item, Instance) and
                                    erased.item.type.is_protocol):
                                # We allow the explicit self-type to be not a supertype of
                                # the current class if it is a protocol. For such cases
                                # the consistency check will be performed at call sites.
                                msg = None
                            elif typ.arg_names[i] in {'self', 'cls'}:
                                if (self.options.python_version[0] < 3
                                        and is_same_type(erased, arg_type) and not isclass):
                                    msg = message_registry.INVALID_SELF_TYPE_OR_EXTRA_ARG
                                    note = '(Hint: typically annotations omit the type for self)'
                                else:
                                    msg = message_registry.ERASED_SELF_TYPE_NOT_SUPERTYPE.format(
                                        erased, ref_type)
                            else:
                                msg = message_registry.MISSING_OR_INVALID_SELF_TYPE
                            if msg:
                                self.fail(msg, defn)
                                if note:
                                    self.note(note, defn)
                    elif isinstance(arg_type, TypeVarType):
                        # Refuse covariant parameter type variables
                        # TODO: check recursively for inner type variables
                        if (
                            arg_type.variance == COVARIANT and
                            defn.name not in ('__init__', '__new__')
                        ):
                            ctx: Context = arg_type
                            if ctx.line < 0:
                                ctx = typ
                            self.fail(message_registry.FUNCTION_PARAMETER_CANNOT_BE_COVARIANT, ctx)
                    if typ.arg_kinds[i] == nodes.ARG_STAR:
                        if not isinstance(arg_type, ParamSpecType):
                            # builtins.tuple[T] is typing.Tuple[T, ...]
                            arg_type = self.named_generic_type('builtins.tuple',
                                                               [arg_type])
                    elif typ.arg_kinds[i] == nodes.ARG_STAR2:
                        if not isinstance(arg_type, ParamSpecType):
                            arg_type = self.named_generic_type('builtins.dict',
                                                               [self.str_type(),
                                                                arg_type])
                    item.arguments[i].variable.type = arg_type

                # Type check initialization expressions.
                body_is_trivial = self.is_trivial_body(defn.body)
                self.check_default_args(item, body_is_trivial)

            # Type check body in a new scope.
            with self.binder.top_frame_context():
                with self.scope.push_function(defn):
                    # We suppress reachability warnings when we use TypeVars with value
                    # restrictions: we only want to report a warning if a certain statement is
                    # marked as being suppressed in *all* of the expansions, but we currently
                    # have no good way of doing this.
                    #
                    # TODO: Find a way of working around this limitation
                    if len(expanded) >= 2:
                        self.binder.suppress_unreachable_warnings()
                    self.accept(item.body)
                unreachable = self.binder.is_unreachable()

            if self.options.warn_no_return and not unreachable:
                if (defn.is_generator or
                        is_named_instance(self.return_types[-1], 'typing.AwaitableGenerator')):
                    return_type = self.get_generator_return_type(self.return_types[-1],
                                                                 defn.is_coroutine)
                elif defn.is_coroutine:
                    return_type = self.get_coroutine_return_type(self.return_types[-1])
                else:
                    return_type = self.return_types[-1]

                return_type = get_proper_type(return_type)
                if not isinstance(return_type, (NoneType, AnyType)) and not body_is_trivial:
                    # Control flow fell off the end of a function that was
                    # declared to return a non-None type and is not
                    # entirely pass/Ellipsis/raise NotImplementedError.
                    if isinstance(return_type, UninhabitedType):
                        # This is a NoReturn function
                        self.fail(message_registry.INVALID_IMPLICIT_RETURN, defn)
                    else:
                        self.fail(message_registry.MISSING_RETURN_STATEMENT, defn)

            self.return_types.pop()

            self.binder = old_binder

    def check_default_args(self, item: FuncItem, body_is_trivial: bool) -> None:
        for arg in item.arguments:
            if arg.initializer is None:
                continue
            if body_is_trivial and isinstance(arg.initializer, EllipsisExpr):
                continue
            name = arg.variable.name
            msg = 'Incompatible default for '
            if name.startswith('__tuple_arg_'):
                msg += "tuple argument {}".format(name[12:])
            else:
                msg += 'argument "{}"'.format(name)
            self.check_simple_assignment(
                arg.variable.type,
                arg.initializer,
                context=arg.initializer,
                msg=msg,
                lvalue_name='argument',
                rvalue_name='default',
                code=codes.ASSIGNMENT)

    def is_forward_op_method(self, method_name: str) -> bool:
        if self.options.python_version[0] == 2 and method_name == '__div__':
            return True
        else:
            return method_name in operators.reverse_op_methods

    def is_reverse_op_method(self, method_name: str) -> bool:
        if self.options.python_version[0] == 2 and method_name == '__rdiv__':
            return True
        else:
            return method_name in operators.reverse_op_method_set

    def check_for_missing_annotations(self, fdef: FuncItem) -> None:
        # Check for functions with unspecified/not fully specified types.
        def is_unannotated_any(t: Type) -> bool:
            if not isinstance(t, ProperType):
                return False
            return isinstance(t, AnyType) and t.type_of_any == TypeOfAny.unannotated

        has_explicit_annotation = (isinstance(fdef.type, CallableType)
                                   and any(not is_unannotated_any(t)
                                           for t in fdef.type.arg_types + [fdef.type.ret_type]))

        show_untyped = not self.is_typeshed_stub or self.options.warn_incomplete_stub
        check_incomplete_defs = self.options.disallow_incomplete_defs and has_explicit_annotation
        if show_untyped and (self.options.disallow_untyped_defs or check_incomplete_defs):
            if fdef.type is None and self.options.disallow_untyped_defs:
                if (not fdef.arguments or (len(fdef.arguments) == 1 and
                        (fdef.arg_names[0] == 'self' or fdef.arg_names[0] == 'cls'))):
                    self.fail(message_registry.RETURN_TYPE_EXPECTED, fdef)
                    if not has_return_statement(fdef) and not fdef.is_generator:
                        self.note('Use "-> None" if function does not return a value', fdef,
                                  code=codes.NO_UNTYPED_DEF)
                else:
                    self.fail(message_registry.FUNCTION_TYPE_EXPECTED, fdef)
            elif isinstance(fdef.type, CallableType):
                ret_type = get_proper_type(fdef.type.ret_type)
                if is_unannotated_any(ret_type):
                    self.fail(message_registry.RETURN_TYPE_EXPECTED, fdef)
                elif fdef.is_generator:
                    if is_unannotated_any(self.get_generator_return_type(ret_type,
                                                                         fdef.is_coroutine)):
                        self.fail(message_registry.RETURN_TYPE_EXPECTED, fdef)
                elif fdef.is_coroutine and isinstance(ret_type, Instance):
                    if is_unannotated_any(self.get_coroutine_return_type(ret_type)):
                        self.fail(message_registry.RETURN_TYPE_EXPECTED, fdef)
                if any(is_unannotated_any(t) for t in fdef.type.arg_types):
                    self.fail(message_registry.ARGUMENT_TYPE_EXPECTED, fdef)

    def check___new___signature(self, fdef: FuncDef, typ: CallableType) -> None:
        self_type = fill_typevars_with_any(fdef.info)
        bound_type = bind_self(typ, self_type, is_classmethod=True)
        # Check that __new__ (after binding cls) returns an instance
        # type (or any).
        if isinstance(fdef.info, TypeInfo) and fdef.info.is_metaclass():
            # This is a metaclass, so it must return a new unrelated type.
            self.check_subtype(
                bound_type.ret_type,
                self.type_type(),
                fdef,
                message_registry.INVALID_NEW_TYPE,
                'returns',
                'but must return a subtype of'
            )
        elif not isinstance(get_proper_type(bound_type.ret_type),
                          (AnyType, Instance, TupleType)):
            self.fail(
                message_registry.NON_INSTANCE_NEW_TYPE.format(
                    format_type(bound_type.ret_type)),
                fdef)
        else:
            # And that it returns a subtype of the class
            self.check_subtype(
                bound_type.ret_type,
                self_type,
                fdef,
                message_registry.INVALID_NEW_TYPE,
                'returns',
                'but must return a subtype of'
            )

    def is_trivial_body(self, block: Block) -> bool:
        """Returns 'true' if the given body is "trivial" -- if it contains just a "pass",
        "..." (ellipsis), or "raise NotImplementedError()". A trivial body may also
        start with a statement containing just a string (e.g. a docstring).

        Note: functions that raise other kinds of exceptions do not count as
        "trivial". We use this function to help us determine when it's ok to
        relax certain checks on body, but functions that raise arbitrary exceptions
        are more likely to do non-trivial work. For example:

           def halt(self, reason: str = ...) -> NoReturn:
               raise MyCustomError("Fatal error: " + reason, self.line, self.context)

        A function that raises just NotImplementedError is much less likely to be
        this complex.
        """
        body = block.body

        # Skip a docstring
        if (body and isinstance(body[0], ExpressionStmt) and
                isinstance(body[0].expr, (StrExpr, UnicodeExpr))):
            body = block.body[1:]

        if len(body) == 0:
            # There's only a docstring (or no body at all).
            return True
        elif len(body) > 1:
            return False

        stmt = body[0]

        if isinstance(stmt, RaiseStmt):
            expr = stmt.expr
            if expr is None:
                return False
            if isinstance(expr, CallExpr):
                expr = expr.callee

            return (isinstance(expr, NameExpr)
                    and expr.fullname == 'builtins.NotImplementedError')

        return (isinstance(stmt, PassStmt) or
                (isinstance(stmt, ExpressionStmt) and
                 isinstance(stmt.expr, EllipsisExpr)))

    def check_reverse_op_method(self, defn: FuncItem,
                                reverse_type: CallableType, reverse_name: str,
                                context: Context) -> None:
        """Check a reverse operator method such as __radd__."""
        # Decides whether it's worth calling check_overlapping_op_methods().

        # This used to check for some very obscure scenario.  It now
        # just decides whether it's worth calling
        # check_overlapping_op_methods().

        assert defn.info

        # First check for a valid signature
        method_type = CallableType([AnyType(TypeOfAny.special_form),
                                    AnyType(TypeOfAny.special_form)],
                                   [nodes.ARG_POS, nodes.ARG_POS],
                                   [None, None],
                                   AnyType(TypeOfAny.special_form),
                                   self.named_type('builtins.function'))
        if not is_subtype(reverse_type, method_type):
            self.msg.invalid_signature(reverse_type, context)
            return

        if reverse_name in ('__eq__', '__ne__'):
            # These are defined for all objects => can't cause trouble.
            return

        # With 'Any' or 'object' return type we are happy, since any possible
        # return value is valid.
        ret_type = get_proper_type(reverse_type.ret_type)
        if isinstance(ret_type, AnyType):
            return
        if isinstance(ret_type, Instance):
            if ret_type.type.fullname == 'builtins.object':
                return
        if reverse_type.arg_kinds[0] == ARG_STAR:
            reverse_type = reverse_type.copy_modified(arg_types=[reverse_type.arg_types[0]] * 2,
                                                      arg_kinds=[ARG_POS] * 2,
                                                      arg_names=[reverse_type.arg_names[0], "_"])
        assert len(reverse_type.arg_types) >= 2

        if self.options.python_version[0] == 2 and reverse_name == '__rdiv__':
            forward_name = '__div__'
        else:
            forward_name = operators.normal_from_reverse_op[reverse_name]
        forward_inst = get_proper_type(reverse_type.arg_types[1])
        if isinstance(forward_inst, TypeVarType):
            forward_inst = get_proper_type(forward_inst.upper_bound)
        elif isinstance(forward_inst, TupleType):
            forward_inst = tuple_fallback(forward_inst)
        elif isinstance(forward_inst, (FunctionLike, TypedDictType, LiteralType)):
            forward_inst = forward_inst.fallback
        if isinstance(forward_inst, TypeType):
            item = forward_inst.item
            if isinstance(item, Instance):
                opt_meta = item.type.metaclass_type
                if opt_meta is not None:
                    forward_inst = opt_meta
        if not (isinstance(forward_inst, (Instance, UnionType))
                and forward_inst.has_readable_member(forward_name)):
            return
        forward_base = reverse_type.arg_types[1]
        forward_type = self.expr_checker.analyze_external_member_access(forward_name, forward_base,
                                                                        context=defn)
        self.check_overlapping_op_methods(reverse_type, reverse_name, defn.info,
                                          forward_type, forward_name, forward_base,
                                          context=defn)

    def check_overlapping_op_methods(self,
                                     reverse_type: CallableType,
                                     reverse_name: str,
                                     reverse_class: TypeInfo,
                                     forward_type: Type,
                                     forward_name: str,
                                     forward_base: Type,
                                     context: Context) -> None:
        """Check for overlapping method and reverse method signatures.

        This function assumes that:

        -   The reverse method has valid argument count and kinds.
        -   If the reverse operator method accepts some argument of type
            X, the forward operator method also belong to class X.

            For example, if we have the reverse operator `A.__radd__(B)`, then the
            corresponding forward operator must have the type `B.__add__(...)`.
        """

        # Note: Suppose we have two operator methods "A.__rOP__(B) -> R1" and
        # "B.__OP__(C) -> R2". We check if these two methods are unsafely overlapping
        # by using the following algorithm:
        #
        # 1. Rewrite "B.__OP__(C) -> R1"  to "temp1(B, C) -> R1"
        #
        # 2. Rewrite "A.__rOP__(B) -> R2" to "temp2(B, A) -> R2"
        #
        # 3. Treat temp1 and temp2 as if they were both variants in the same
        #    overloaded function. (This mirrors how the Python runtime calls
        #    operator methods: we first try __OP__, then __rOP__.)
        #
        #    If the first signature is unsafely overlapping with the second,
        #    report an error.
        #
        # 4. However, if temp1 shadows temp2 (e.g. the __rOP__ method can never
        #    be called), do NOT report an error.
        #
        #    This behavior deviates from how we handle overloads -- many of the
        #    modules in typeshed seem to define __OP__ methods that shadow the
        #    corresponding __rOP__ method.
        #
        # Note: we do not attempt to handle unsafe overlaps related to multiple
        # inheritance. (This is consistent with how we handle overloads: we also
        # do not try checking unsafe overlaps due to multiple inheritance there.)

        for forward_item in union_items(forward_type):
            if isinstance(forward_item, CallableType):
                if self.is_unsafe_overlapping_op(forward_item, forward_base, reverse_type):
                    self.msg.operator_method_signatures_overlap(
                        reverse_class, reverse_name,
                        forward_base, forward_name, context)
            elif isinstance(forward_item, Overloaded):
                for item in forward_item.items:
                    if self.is_unsafe_overlapping_op(item, forward_base, reverse_type):
                        self.msg.operator_method_signatures_overlap(
                            reverse_class, reverse_name,
                            forward_base, forward_name,
                            context)
            elif not isinstance(forward_item, AnyType):
                self.msg.forward_operator_not_callable(forward_name, context)

    def is_unsafe_overlapping_op(self,
                                 forward_item: CallableType,
                                 forward_base: Type,
                                 reverse_type: CallableType) -> bool:
        # TODO: check argument kinds?
        if len(forward_item.arg_types) < 1:
            # Not a valid operator method -- can't succeed anyway.
            return False

        # Erase the type if necessary to make sure we don't have a single
        # TypeVar in forward_tweaked. (Having a function signature containing
        # just a single TypeVar can lead to unpredictable behavior.)
        forward_base_erased = forward_base
        if isinstance(forward_base, TypeVarType):
            forward_base_erased = erase_to_bound(forward_base)

        # Construct normalized function signatures corresponding to the
        # operator methods. The first argument is the left operand and the
        # second operand is the right argument -- we switch the order of
        # the arguments of the reverse method.

        forward_tweaked = forward_item.copy_modified(
            arg_types=[forward_base_erased, forward_item.arg_types[0]],
            arg_kinds=[nodes.ARG_POS] * 2,
            arg_names=[None] * 2,
        )
        reverse_tweaked = reverse_type.copy_modified(
            arg_types=[reverse_type.arg_types[1], reverse_type.arg_types[0]],
            arg_kinds=[nodes.ARG_POS] * 2,
            arg_names=[None] * 2,
        )

        reverse_base_erased = reverse_type.arg_types[0]
        if isinstance(reverse_base_erased, TypeVarType):
            reverse_base_erased = erase_to_bound(reverse_base_erased)

        if is_same_type(reverse_base_erased, forward_base_erased):
            return False
        elif is_subtype(reverse_base_erased, forward_base_erased):
            first = reverse_tweaked
            second = forward_tweaked
        else:
            first = forward_tweaked
            second = reverse_tweaked

        return is_unsafe_overlapping_overload_signatures(first, second)

    def check_inplace_operator_method(self, defn: FuncBase) -> None:
        """Check an inplace operator method such as __iadd__.

        They cannot arbitrarily overlap with __add__.
        """
        method = defn.name
        if method not in operators.inplace_operator_methods:
            return
        typ = bind_self(self.function_type(defn))
        cls = defn.info
        other_method = '__' + method[3:]
        if cls.has_readable_member(other_method):
            instance = fill_typevars(cls)
            typ2 = get_proper_type(self.expr_checker.analyze_external_member_access(
                other_method, instance, defn))
            fail = False
            if isinstance(typ2, FunctionLike):
                if not is_more_general_arg_prefix(typ, typ2):
                    fail = True
            else:
                # TODO overloads
                fail = True
            if fail:
                self.msg.signatures_incompatible(method, other_method, defn)

    def check_getattr_method(self, typ: Type, context: Context, name: str) -> None:
        if len(self.scope.stack) == 1:
            # module scope
            if name == '__getattribute__':
                self.fail(message_registry.MODULE_LEVEL_GETATTRIBUTE, context)
                return
            # __getattr__ is fine at the module level as of Python 3.7 (PEP 562). We could
            # show an error for Python < 3.7, but that would be annoying in code that supports
            # both 3.7 and older versions.
            method_type = CallableType([self.named_type('builtins.str')],
                                       [nodes.ARG_POS],
                                       [None],
                                       AnyType(TypeOfAny.special_form),
                                       self.named_type('builtins.function'))
        elif self.scope.active_class():
            method_type = CallableType([AnyType(TypeOfAny.special_form),
                                        self.named_type('builtins.str')],
                                       [nodes.ARG_POS, nodes.ARG_POS],
                                       [None, None],
                                       AnyType(TypeOfAny.special_form),
                                       self.named_type('builtins.function'))
        else:
            return
        if not is_subtype(typ, method_type):
            self.msg.invalid_signature_for_special_method(typ, context, name)

    def check_setattr_method(self, typ: Type, context: Context) -> None:
        if not self.scope.active_class():
            return
        method_type = CallableType([AnyType(TypeOfAny.special_form),
                                    self.named_type('builtins.str'),
                                    AnyType(TypeOfAny.special_form)],
                                   [nodes.ARG_POS, nodes.ARG_POS, nodes.ARG_POS],
                                   [None, None, None],
                                   NoneType(),
                                   self.named_type('builtins.function'))
        if not is_subtype(typ, method_type):
            self.msg.invalid_signature_for_special_method(typ, context, '__setattr__')

    def expand_typevars(self, defn: FuncItem,
                        typ: CallableType) -> List[Tuple[FuncItem, CallableType]]:
        # TODO use generator
        subst: List[List[Tuple[TypeVarId, Type]]] = []
        tvars = list(typ.variables) or []
        if defn.info:
            # Class type variables
            tvars += defn.info.defn.type_vars or []
        # TODO(PEP612): audit for paramspec
        for tvar in tvars:
            if isinstance(tvar, TypeVarType) and tvar.values:
                subst.append([(tvar.id, value) for value in tvar.values])
        # Make a copy of the function to check for each combination of
        # value restricted type variables. (Except when running mypyc,
        # where we need one canonical version of the function.)
        if subst and not self.options.mypyc:
            result: List[Tuple[FuncItem, CallableType]] = []
            for substitutions in itertools.product(*subst):
                mapping = dict(substitutions)
                expanded = cast(CallableType, expand_type(typ, mapping))
                result.append((expand_func(defn, mapping), expanded))
            return result
        else:
            return [(defn, typ)]

    def check_method_override(self, defn: Union[FuncDef, OverloadedFuncDef, Decorator]) -> None:
        """Check if function definition is compatible with base classes.

        This may defer the method if a signature is not available in at least one base class.
        """
        # Check against definitions in base classes.
        for base in defn.info.mro[1:]:
            if self.check_method_or_accessor_override_for_base(defn, base):
                # Node was deferred, we will have another attempt later.
                return

    def check_method_or_accessor_override_for_base(self, defn: Union[FuncDef,
                                                                     OverloadedFuncDef,
                                                                     Decorator],
                                                   base: TypeInfo) -> bool:
        """Check if method definition is compatible with a base class.

        Return True if the node was deferred because one of the corresponding
        superclass nodes is not ready.
        """
        if base:
            name = defn.name
            base_attr = base.names.get(name)
            if base_attr:
                # First, check if we override a final (always an error, even with Any types).
                if is_final_node(base_attr.node):
                    self.msg.cant_override_final(name, base.name, defn)
                # Second, final can't override anything writeable independently of types.
                if defn.is_final:
                    self.check_if_final_var_override_writable(name, base_attr.node, defn)

            # Check the type of override.
            if name not in ('__init__', '__new__', '__init_subclass__'):
                # Check method override
                # (__init__, __new__, __init_subclass__ are special).
                if self.check_method_override_for_base_with_name(defn, name, base):
                    return True
                if name in operators.inplace_operator_methods:
                    # Figure out the name of the corresponding operator method.
                    method = '__' + name[3:]
                    # An inplace operator method such as __iadd__ might not be
                    # always introduced safely if a base class defined __add__.
                    # TODO can't come up with an example where this is
                    #      necessary; now it's "just in case"
                    return self.check_method_override_for_base_with_name(defn, method,
                                                                         base)
        return False

    def check_method_override_for_base_with_name(
            self, defn: Union[FuncDef, OverloadedFuncDef, Decorator],
            name: str, base: TypeInfo) -> bool:
        """Check if overriding an attribute `name` of `base` with `defn` is valid.

        Return True if the supertype node was not analysed yet, and `defn` was deferred.
        """
        base_attr = base.names.get(name)
        if base_attr:
            # The name of the method is defined in the base class.

            # Point errors at the 'def' line (important for backward compatibility
            # of type ignores).
            if not isinstance(defn, Decorator):
                context = defn
            else:
                context = defn.func

            # Construct the type of the overriding method.
            if isinstance(defn, (FuncDef, OverloadedFuncDef)):
                typ: Type = self.function_type(defn)
                override_class_or_static = defn.is_class or defn.is_static
                override_class = defn.is_class
            else:
                assert defn.var.is_ready
                assert defn.var.type is not None
                typ = defn.var.type
                override_class_or_static = defn.func.is_class or defn.func.is_static
                override_class = defn.func.is_class
            typ = get_proper_type(typ)
            if isinstance(typ, FunctionLike) and not is_static(context):
                typ = bind_self(typ, self.scope.active_self_type(),
                                is_classmethod=override_class)
            # Map the overridden method type to subtype context so that
            # it can be checked for compatibility.
            original_type = get_proper_type(base_attr.type)
            original_node = base_attr.node
            if original_type is None:
                if self.pass_num < self.last_pass:
                    # If there are passes left, defer this node until next pass,
                    # otherwise try reconstructing the method type from available information.
                    self.defer_node(defn, defn.info)
                    return True
                elif isinstance(original_node, (FuncDef, OverloadedFuncDef)):
                    original_type = self.function_type(original_node)
                elif isinstance(original_node, Decorator):
                    original_type = self.function_type(original_node.func)
                elif isinstance(original_node, Var):
                    # Super type can define method as an attribute.
                    # See https://github.com/python/mypy/issues/10134

                    # We also check that sometimes `original_node.type` is None.
                    # This is the case when we use something like `__hash__ = None`.
                    if original_node.type is not None:
                        original_type = get_proper_type(original_node.type)
                    else:
                        original_type = NoneType()
                else:
                    assert False, str(base_attr.node)
            if isinstance(original_node, (FuncDef, OverloadedFuncDef)):
                original_class_or_static = original_node.is_class or original_node.is_static
            elif isinstance(original_node, Decorator):
                fdef = original_node.func
                original_class_or_static = fdef.is_class or fdef.is_static
            else:
                original_class_or_static = False  # a variable can't be class or static
            if isinstance(original_type, AnyType) or isinstance(typ, AnyType):
                pass
            elif isinstance(original_type, FunctionLike) and isinstance(typ, FunctionLike):
                original = self.bind_and_map_method(base_attr, original_type,
                                                    defn.info, base)
                # Check that the types are compatible.
                # TODO overloaded signatures
                self.check_override(typ,
                                    original,
                                    defn.name,
                                    name,
                                    base.name,
                                    original_class_or_static,
                                    override_class_or_static,
                                    context)
            elif is_equivalent(original_type, typ):
                # Assume invariance for a non-callable attribute here. Note
                # that this doesn't affect read-only properties which can have
                # covariant overrides.
                #
                pass
            elif (base_attr.node and not self.is_writable_attribute(base_attr.node)
                  and is_subtype(typ, original_type)):
                # If the attribute is read-only, allow covariance
                pass
            else:
                self.msg.signature_incompatible_with_supertype(
                    defn.name, name, base.name, context)
        return False

    def bind_and_map_method(self, sym: SymbolTableNode, typ: FunctionLike,
                            sub_info: TypeInfo, super_info: TypeInfo) -> FunctionLike:
        """Bind self-type and map type variables for a method.

        Arguments:
            sym: a symbol that points to method definition
            typ: method type on the definition
            sub_info: class where the method is used
            super_info: class where the method was defined
        """
        if (isinstance(sym.node, (FuncDef, OverloadedFuncDef, Decorator))
                and not is_static(sym.node)):
            if isinstance(sym.node, Decorator):
                is_class_method = sym.node.func.is_class
            else:
                is_class_method = sym.node.is_class
            bound = bind_self(typ, self.scope.active_self_type(), is_class_method)
        else:
            bound = typ
        return cast(FunctionLike, map_type_from_supertype(bound, sub_info, super_info))

    def get_op_other_domain(self, tp: FunctionLike) -> Optional[Type]:
        if isinstance(tp, CallableType):
            if tp.arg_kinds and tp.arg_kinds[0] == ARG_POS:
                return tp.arg_types[0]
            return None
        elif isinstance(tp, Overloaded):
            raw_items = [self.get_op_other_domain(it) for it in tp.items]
            items = [it for it in raw_items if it]
            if items:
                return make_simplified_union(items)
            return None
        else:
            assert False, "Need to check all FunctionLike subtypes here"

    def check_override(self, override: FunctionLike, original: FunctionLike,
                       name: str, name_in_super: str, supertype: str,
                       original_class_or_static: bool,
                       override_class_or_static: bool,
                       node: Context) -> None:
        """Check a method override with given signatures.

        Arguments:
          override:  The signature of the overriding method.
          original:  The signature of the original supertype method.
          name:      The name of the subtype. This and the next argument are
                     only used for generating error messages.
          supertype: The name of the supertype.
        """
        # Use boolean variable to clarify code.
        fail = False
        op_method_wider_note = False
        if not is_subtype(override, original, ignore_pos_arg_names=True):
            fail = True
        elif isinstance(override, Overloaded) and self.is_forward_op_method(name):
            # Operator method overrides cannot extend the domain, as
            # this could be unsafe with reverse operator methods.
            original_domain = self.get_op_other_domain(original)
            override_domain = self.get_op_other_domain(override)
            if (original_domain and override_domain and
                    not is_subtype(override_domain, original_domain)):
                fail = True
                op_method_wider_note = True
        if isinstance(original, FunctionLike) and isinstance(override, FunctionLike):
            if original_class_or_static and not override_class_or_static:
                fail = True
            elif isinstance(original, CallableType) and isinstance(override, CallableType):
                if original.type_guard is not None and override.type_guard is None:
                    fail = True

        if is_private(name):
            fail = False

        if fail:
            emitted_msg = False
            if (isinstance(override, CallableType) and
                    isinstance(original, CallableType) and
                    len(override.arg_types) == len(original.arg_types) and
                    override.min_args == original.min_args):
                # Give more detailed messages for the common case of both
                # signatures having the same number of arguments and no
                # overloads.

                # override might have its own generic function type
                # variables. If an argument or return type of override
                # does not have the correct subtyping relationship
                # with the original type even after these variables
                # are erased, then it is definitely an incompatibility.

                override_ids = override.type_var_ids()
                type_name = None
                if isinstance(override.definition, FuncDef):
                    type_name = override.definition.info.name

                def erase_override(t: Type) -> Type:
                    return erase_typevars(t, ids_to_erase=override_ids)

                for i in range(len(override.arg_types)):
                    if not is_subtype(original.arg_types[i],
                                      erase_override(override.arg_types[i])):
                        arg_type_in_super = original.arg_types[i]
                        self.msg.argument_incompatible_with_supertype(
                            i + 1,
                            name,
                            type_name,
                            name_in_super,
                            arg_type_in_super,
                            supertype,
                            node
                        )
                        emitted_msg = True

                if not is_subtype(erase_override(override.ret_type),
                                  original.ret_type):
                    self.msg.return_type_incompatible_with_supertype(
                        name, name_in_super, supertype, original.ret_type, override.ret_type, node)
                    emitted_msg = True
            elif isinstance(override, Overloaded) and isinstance(original, Overloaded):
                # Give a more detailed message in the case where the user is trying to
                # override an overload, and the subclass's overload is plausible, except
                # that the order of the variants are wrong.
                #
                # For example, if the parent defines the overload f(int) -> int and f(str) -> str
                # (in that order), and if the child swaps the two and does f(str) -> str and
                # f(int) -> int
                order = []
                for child_variant in override.items:
                    for i, parent_variant in enumerate(original.items):
                        if is_subtype(child_variant, parent_variant):
                            order.append(i)
                            break

                if len(order) == len(original.items) and order != sorted(order):
                    self.msg.overload_signature_incompatible_with_supertype(
                        name, name_in_super, supertype, override, node)
                    emitted_msg = True

            if not emitted_msg:
                # Fall back to generic incompatibility message.
                self.msg.signature_incompatible_with_supertype(
                    name, name_in_super, supertype, node, original=original, override=override)
            if op_method_wider_note:
                self.note("Overloaded operator methods can't have wider argument types"
                          " in overrides", node, code=codes.OVERRIDE)

    def check__exit__return_type(self, defn: FuncItem) -> None:
        """Generate error if the return type of __exit__ is problematic.

        If __exit__ always returns False but the return type is declared
        as bool, mypy thinks that a with statement may "swallow"
        exceptions even though this is not the case, resulting in
        invalid reachability inference.
        """
        if not defn.type or not isinstance(defn.type, CallableType):
            return

        ret_type = get_proper_type(defn.type.ret_type)
        if not has_bool_item(ret_type):
            return

        returns = all_return_statements(defn)
        if not returns:
            return

        if all(isinstance(ret.expr, NameExpr) and ret.expr.fullname == 'builtins.False'
               for ret in returns):
            self.msg.incorrect__exit__return(defn)

    def visit_class_def(self, defn: ClassDef) -> None:
        """Type check a class definition."""
        typ = defn.info
        for base in typ.mro[1:]:
            if base.is_final:
                self.fail(message_registry.CANNOT_INHERIT_FROM_FINAL.format(base.name), defn)
        with self.tscope.class_scope(defn.info), self.enter_partial_types(is_class=True):
            old_binder = self.binder
            self.binder = ConditionalTypeBinder()
            with self.binder.top_frame_context():
                with self.scope.push_class(defn.info):
                    self.accept(defn.defs)
            self.binder = old_binder
            if not (defn.info.typeddict_type or defn.info.tuple_type or defn.info.is_enum):
                # If it is not a normal class (not a special form) check class keywords.
                self.check_init_subclass(defn)
            if not defn.has_incompatible_baseclass:
                # Otherwise we've already found errors; more errors are not useful
                self.check_multiple_inheritance(typ)
            self.check_final_deletable(typ)

            if defn.decorators:
                sig: Type = type_object_type(defn.info, self.named_type)
                # Decorators are applied in reverse order.
                for decorator in reversed(defn.decorators):
                    if (isinstance(decorator, CallExpr)
                            and isinstance(decorator.analyzed, PromoteExpr)):
                        # _promote is a special type checking related construct.
                        continue

                    dec = self.expr_checker.accept(decorator)
                    temp = self.temp_node(sig, context=decorator)
                    fullname = None
                    if isinstance(decorator, RefExpr):
                        fullname = decorator.fullname

                    # TODO: Figure out how to have clearer error messages.
                    # (e.g. "class decorator must be a function that accepts a type."
                    sig, _ = self.expr_checker.check_call(dec, [temp],
                                                          [nodes.ARG_POS], defn,
                                                          callable_name=fullname)
                # TODO: Apply the sig to the actual TypeInfo so we can handle decorators
                # that completely swap out the type.  (e.g. Callable[[Type[A]], Type[B]])
        if typ.is_protocol and typ.defn.type_vars:
            self.check_protocol_variance(defn)

    def check_final_deletable(self, typ: TypeInfo) -> None:
        # These checks are only for mypyc. Only perform some checks that are easier
        # to implement here than in mypyc.
        for attr in typ.deletable_attributes:
            node = typ.names.get(attr)
            if node and isinstance(node.node, Var) and node.node.is_final:
                self.fail(message_registry.CANNOT_MAKE_DELETABLE_FINAL, node.node)

    def check_init_subclass(self, defn: ClassDef) -> None:
        """Check that keywords in a class definition are valid arguments for __init_subclass__().

        In this example:
            1   class Base:
            2       def __init_subclass__(cls, thing: int):
            3           pass
            4   class Child(Base, thing=5):
            5       def __init_subclass__(cls):
            6           pass
            7   Child()

        Base.__init_subclass__(thing=5) is called at line 4. This is what we simulate here.
        Child.__init_subclass__ is never called.
        """
        if (defn.info.metaclass_type and
                defn.info.metaclass_type.type.fullname not in ('builtins.type', 'abc.ABCMeta')):
            # We can't safely check situations when both __init_subclass__ and a custom
            # metaclass are present.
            return
        # At runtime, only Base.__init_subclass__ will be called, so
        # we skip the current class itself.
        for base in defn.info.mro[1:]:
            if '__init_subclass__' not in base.names:
                continue
            name_expr = NameExpr(defn.name)
            name_expr.node = base
            callee = MemberExpr(name_expr, '__init_subclass__')
            args = list(defn.keywords.values())
            arg_names: List[Optional[str]] = list(defn.keywords.keys())
            # 'metaclass' keyword is consumed by the rest of the type machinery,
            # and is never passed to __init_subclass__ implementations
            if 'metaclass' in arg_names:
                idx = arg_names.index('metaclass')
                arg_names.pop(idx)
                args.pop(idx)
            arg_kinds = [ARG_NAMED] * len(args)
            call_expr = CallExpr(callee, args, arg_kinds, arg_names)
            call_expr.line = defn.line
            call_expr.column = defn.column
            call_expr.end_line = defn.end_line
            self.expr_checker.accept(call_expr,
                                     allow_none_return=True,
                                     always_allow_any=True)
            # We are only interested in the first Base having __init_subclass__,
            # all other bases have already been checked.
            break

    def check_protocol_variance(self, defn: ClassDef) -> None:
        """Check that protocol definition is compatible with declared
        variances of type variables.

        Note that we also prohibit declaring protocol classes as invariant
        if they are actually covariant/contravariant, since this may break
        transitivity of subtyping, see PEP 544.
        """
        info = defn.info
        object_type = Instance(info.mro[-1], [])
        tvars = info.defn.type_vars
        for i, tvar in enumerate(tvars):
            up_args: List[Type] = [
                object_type if i == j else AnyType(TypeOfAny.special_form)
                for j, _ in enumerate(tvars)
            ]
            down_args: List[Type] = [
                UninhabitedType() if i == j else AnyType(TypeOfAny.special_form)
                for j, _ in enumerate(tvars)
            ]
            up, down = Instance(info, up_args), Instance(info, down_args)
            # TODO: add advanced variance checks for recursive protocols
            if is_subtype(down, up, ignore_declared_variance=True):
                expected = COVARIANT
            elif is_subtype(up, down, ignore_declared_variance=True):
                expected = CONTRAVARIANT
            else:
                expected = INVARIANT
            if isinstance(tvar, TypeVarType) and expected != tvar.variance:
                self.msg.bad_proto_variance(tvar.variance, tvar.name, expected, defn)

    def check_multiple_inheritance(self, typ: TypeInfo) -> None:
        """Check for multiple inheritance related errors."""
        if len(typ.bases) <= 1:
            # No multiple inheritance.
            return
        # Verify that inherited attributes are compatible.
        mro = typ.mro[1:]
        for i, base in enumerate(mro):
            # Attributes defined in both the type and base are skipped.
            # Normal checks for attribute compatibility should catch any problems elsewhere.
            non_overridden_attrs = base.names.keys() - typ.names.keys()
            for name in non_overridden_attrs:
                if is_private(name):
                    continue
                for base2 in mro[i + 1:]:
                    # We only need to check compatibility of attributes from classes not
                    # in a subclass relationship. For subclasses, normal (single inheritance)
                    # checks suffice (these are implemented elsewhere).
                    if name in base2.names and base2 not in base.mro:
                        self.check_compatibility(name, base, base2, typ)

    def determine_type_of_class_member(self, sym: SymbolTableNode) -> Optional[Type]:
        if sym.type is not None:
            return sym.type
        if isinstance(sym.node, FuncBase):
            return self.function_type(sym.node)
        if isinstance(sym.node, TypeInfo):
            # nested class
            return type_object_type(sym.node, self.named_type)
        if isinstance(sym.node, TypeVarExpr):
            # Use of TypeVars is rejected in an expression/runtime context, so
            # we don't need to check supertype compatibility for them.
            return AnyType(TypeOfAny.special_form)
        return None

    def check_compatibility(self, name: str, base1: TypeInfo,
                            base2: TypeInfo, ctx: TypeInfo) -> None:
        """Check if attribute name in base1 is compatible with base2 in multiple inheritance.

        Assume base1 comes before base2 in the MRO, and that base1 and base2 don't have
        a direct subclass relationship (i.e., the compatibility requirement only derives from
        multiple inheritance).

        This check verifies that a definition taken from base1 (and mapped to the current
        class ctx), is type compatible with the definition taken from base2 (also mapped), so
        that unsafe subclassing like this can be detected:
            class A(Generic[T]):
                def foo(self, x: T) -> None: ...

            class B:
                def foo(self, x: str) -> None: ...

            class C(B, A[int]): ...  # this is unsafe because...

            x: A[int] = C()
            x.foo  # ...runtime type is (str) -> None, while static type is (int) -> None
        """
        if name in ('__init__', '__new__', '__init_subclass__'):
            # __init__ and friends can be incompatible -- it's a special case.
            return
        first = base1.names[name]
        second = base2.names[name]
        first_type = get_proper_type(self.determine_type_of_class_member(first))
        second_type = get_proper_type(self.determine_type_of_class_member(second))

        if (isinstance(first_type, FunctionLike) and
                isinstance(second_type, FunctionLike)):
            if first_type.is_type_obj() and second_type.is_type_obj():
                # For class objects only check the subtype relationship of the classes,
                # since we allow incompatible overrides of '__init__'/'__new__'
                ok = is_subtype(left=fill_typevars_with_any(first_type.type_object()),
                                right=fill_typevars_with_any(second_type.type_object()))
            else:
                # First bind/map method types when necessary.
                first_sig = self.bind_and_map_method(first, first_type, ctx, base1)
                second_sig = self.bind_and_map_method(second, second_type, ctx, base2)
                ok = is_subtype(first_sig, second_sig, ignore_pos_arg_names=True)
        elif first_type and second_type:
            ok = is_equivalent(first_type, second_type)
            if not ok:
                second_node = base2[name].node
                if isinstance(second_node, Decorator) and second_node.func.is_property:
                    ok = is_subtype(first_type, cast(CallableType, second_type).ret_type)
        else:
            if first_type is None:
                self.msg.cannot_determine_type_in_base(name, base1.name, ctx)
            if second_type is None:
                self.msg.cannot_determine_type_in_base(name, base2.name, ctx)
            ok = True
        # Final attributes can never be overridden, but can override
        # non-final read-only attributes.
        if is_final_node(second.node):
            self.msg.cant_override_final(name, base2.name, ctx)
        if is_final_node(first.node):
            self.check_if_final_var_override_writable(name, second.node, ctx)
        # __slots__ and __deletable__ are special and the type can vary across class hierarchy.
        if name in ('__slots__', '__deletable__'):
            ok = True
        if not ok:
            self.msg.base_class_definitions_incompatible(name, base1, base2,
                                                         ctx)

    def visit_import_from(self, node: ImportFrom) -> None:
        self.check_import(node)

    def visit_import_all(self, node: ImportAll) -> None:
        self.check_import(node)

    def visit_import(self, s: Import) -> None:
        pass

    def check_import(self, node: ImportBase) -> None:
        for assign in node.assignments:
            lvalue = assign.lvalues[0]
            lvalue_type, _, __ = self.check_lvalue(lvalue)
            if lvalue_type is None:
                # TODO: This is broken.
                lvalue_type = AnyType(TypeOfAny.special_form)
            message = '{} "{}"'.format(message_registry.INCOMPATIBLE_IMPORT_OF,
                                       cast(NameExpr, assign.rvalue).name)
            self.check_simple_assignment(lvalue_type, assign.rvalue, node,
                                         msg=message, lvalue_name='local name',
                                         rvalue_name='imported name')

    #
    # Statements
    #

    def visit_block(self, b: Block) -> None:
        if b.is_unreachable:
            # This block was marked as being unreachable during semantic analysis.
            # It turns out any blocks marked in this way are *intentionally* marked
            # as unreachable -- so we don't display an error.
            self.binder.unreachable()
            return
        for s in b.body:
            if self.binder.is_unreachable():
                if self.should_report_unreachable_issues() and not self.is_raising_or_empty(s):
                    self.msg.unreachable_statement(s)
                break
            self.accept(s)

    def should_report_unreachable_issues(self) -> bool:
        return (self.in_checked_function()
                and self.options.warn_unreachable
                and not self.binder.is_unreachable_warning_suppressed())

    def is_raising_or_empty(self, s: Statement) -> bool:
        """Returns 'true' if the given statement either throws an error of some kind
        or is a no-op.

        We use this function mostly while handling the '--warn-unreachable' flag. When
        that flag is present, we normally report an error on any unreachable statement.
        But if that statement is just something like a 'pass' or a just-in-case 'assert False',
        reporting an error would be annoying.
        """
        if isinstance(s, AssertStmt) and is_false_literal(s.expr):
            return True
        elif isinstance(s, (RaiseStmt, PassStmt)):
            return True
        elif isinstance(s, ExpressionStmt):
            if isinstance(s.expr, EllipsisExpr):
                return True
            elif isinstance(s.expr, CallExpr):
                with self.expr_checker.msg.disable_errors():
                    typ = get_proper_type(self.expr_checker.accept(
                        s.expr, allow_none_return=True, always_allow_any=True))

                if isinstance(typ, UninhabitedType):
                    return True
        return False

    def visit_assignment_stmt(self, s: AssignmentStmt) -> None:
        """Type check an assignment statement.

        Handle all kinds of assignment statements (simple, indexed, multiple).
        """
        # Avoid type checking type aliases in stubs to avoid false
        # positives about modern type syntax available in stubs such
        # as X | Y.
        if not (s.is_alias_def and self.is_stub):
            with self.enter_final_context(s.is_final_def):
                self.check_assignment(s.lvalues[-1], s.rvalue, s.type is None, s.new_syntax)

        if s.is_alias_def:
            self.check_type_alias_rvalue(s)

        if (s.type is not None and
                self.options.disallow_any_unimported and
                has_any_from_unimported_type(s.type)):
            if isinstance(s.lvalues[-1], TupleExpr):
                # This is a multiple assignment. Instead of figuring out which type is problematic,
                # give a generic error message.
                self.msg.unimported_type_becomes_any("A type on this line",
                                                     AnyType(TypeOfAny.special_form), s)
            else:
                self.msg.unimported_type_becomes_any("Type of variable", s.type, s)
        check_for_explicit_any(s.type, self.options, self.is_typeshed_stub, self.msg, context=s)

        if len(s.lvalues) > 1:
            # Chained assignment (e.g. x = y = ...).
            # Make sure that rvalue type will not be reinferred.
            if s.rvalue not in self.type_map:
                self.expr_checker.accept(s.rvalue)
            rvalue = self.temp_node(self.type_map[s.rvalue], s)
            for lv in s.lvalues[:-1]:
                with self.enter_final_context(s.is_final_def):
                    self.check_assignment(lv, rvalue, s.type is None)

        self.check_final(s)
        if (s.is_final_def and s.type and not has_no_typevars(s.type)
                and self.scope.active_class() is not None):
            self.fail(message_registry.DEPENDENT_FINAL_IN_CLASS_BODY, s)

    def check_type_alias_rvalue(self, s: AssignmentStmt) -> None:
        if not (self.is_stub and isinstance(s.rvalue, OpExpr) and s.rvalue.op == '|'):
            # We do this mostly for compatibility with old semantic analyzer.
            # TODO: should we get rid of this?
            alias_type = self.expr_checker.accept(s.rvalue)
        else:
            # Avoid type checking 'X | Y' in stubs, since there can be errors
            # on older Python targets.
            alias_type = AnyType(TypeOfAny.special_form)

            def accept_items(e: Expression) -> None:
                if isinstance(e, OpExpr) and e.op == '|':
                    accept_items(e.left)
                    accept_items(e.right)
                else:
                    # Nested union types have been converted to type context
                    # in semantic analysis (such as in 'list[int | str]'),
                    # so we don't need to deal with them here.
                    self.expr_checker.accept(e)

            accept_items(s.rvalue)
        self.store_type(s.lvalues[-1], alias_type)

    def check_assignment(self, lvalue: Lvalue, rvalue: Expression, infer_lvalue_type: bool = True,
                         new_syntax: bool = False) -> None:
        """Type check a single assignment: lvalue = rvalue."""
        if isinstance(lvalue, TupleExpr) or isinstance(lvalue, ListExpr):
            self.check_assignment_to_multiple_lvalues(lvalue.items, rvalue, rvalue,
                                                      infer_lvalue_type)
        else:
            self.try_infer_partial_generic_type_from_assignment(lvalue, rvalue, '=')
            lvalue_type, index_lvalue, inferred = self.check_lvalue(lvalue)
            # If we're assigning to __getattr__ or similar methods, check that the signature is
            # valid.
            if isinstance(lvalue, NameExpr) and lvalue.node:
                name = lvalue.node.name
                if name in ('__setattr__', '__getattribute__', '__getattr__'):
                    # If an explicit type is given, use that.
                    if lvalue_type:
                        signature = lvalue_type
                    else:
                        signature = self.expr_checker.accept(rvalue)
                    if signature:
                        if name == '__setattr__':
                            self.check_setattr_method(signature, lvalue)
                        else:
                            self.check_getattr_method(signature, lvalue, name)

            # Defer PartialType's super type checking.
            if (isinstance(lvalue, RefExpr) and
                    not (isinstance(lvalue_type, PartialType) and lvalue_type.type is None)):
                if self.check_compatibility_all_supers(lvalue, lvalue_type, rvalue):
                    # We hit an error on this line; don't check for any others
                    return

            if lvalue_type:
                if isinstance(lvalue_type, PartialType) and lvalue_type.type is None:
                    # Try to infer a proper type for a variable with a partial None type.
                    rvalue_type = self.expr_checker.accept(rvalue)
                    if isinstance(get_proper_type(rvalue_type), NoneType):
                        # This doesn't actually provide any additional information -- multiple
                        # None initializers preserve the partial None type.
                        return

                    if is_valid_inferred_type(rvalue_type):
                        var = lvalue_type.var
                        partial_types = self.find_partial_types(var)
                        if partial_types is not None:
                            if not self.current_node_deferred:
                                # Partial type can't be final, so strip any literal values.
                                rvalue_type = remove_instance_last_known_values(rvalue_type)
                                inferred_type = make_simplified_union(
                                    [rvalue_type, NoneType()])
                                self.set_inferred_type(var, lvalue, inferred_type)
                            else:
                                var.type = None
                            del partial_types[var]
                            lvalue_type = var.type
                    else:
                        # Try to infer a partial type. No need to check the return value, as
                        # an error will be reported elsewhere.
                        self.infer_partial_type(lvalue_type.var, lvalue, rvalue_type)
                    # Handle None PartialType's super type checking here, after it's resolved.
                    if (isinstance(lvalue, RefExpr) and
                            self.check_compatibility_all_supers(lvalue, lvalue_type, rvalue)):
                        # We hit an error on this line; don't check for any others
                        return
                elif (is_literal_none(rvalue) and
                        isinstance(lvalue, NameExpr) and
                        isinstance(lvalue.node, Var) and
                        lvalue.node.is_initialized_in_class and
                        not new_syntax):
                    # Allow None's to be assigned to class variables with non-Optional types.
                    rvalue_type = lvalue_type
                elif (isinstance(lvalue, MemberExpr) and
                        lvalue.kind is None):  # Ignore member access to modules
                    instance_type = self.expr_checker.accept(lvalue.expr)
                    rvalue_type, lvalue_type, infer_lvalue_type = self.check_member_assignment(
                        instance_type, lvalue_type, rvalue, context=rvalue)
                else:
                    # Hacky special case for assigning a literal None
                    # to a variable defined in a previous if
                    # branch. When we detect this, we'll go back and
                    # make the type optional. This is somewhat
                    # unpleasant, and a generalization of this would
                    # be an improvement!
                    if (is_literal_none(rvalue) and
                            isinstance(lvalue, NameExpr) and
                            lvalue.kind == LDEF and
                            isinstance(lvalue.node, Var) and
                            lvalue.node.type and
                            lvalue.node in self.var_decl_frames and
                            not isinstance(get_proper_type(lvalue_type), AnyType)):
                        decl_frame_map = self.var_decl_frames[lvalue.node]
                        # Check if the nearest common ancestor frame for the definition site
                        # and the current site is the enclosing frame of an if/elif/else block.
                        has_if_ancestor = False
                        for frame in reversed(self.binder.frames):
                            if frame.id in decl_frame_map:
                                has_if_ancestor = frame.conditional_frame
                                break
                        if has_if_ancestor:
                            lvalue_type = make_optional_type(lvalue_type)
                            self.set_inferred_type(lvalue.node, lvalue, lvalue_type)

                    rvalue_type = self.check_simple_assignment(lvalue_type, rvalue, context=rvalue,
                                                               code=codes.ASSIGNMENT)

                # Special case: only non-abstract non-protocol classes can be assigned to
                # variables with explicit type Type[A], where A is protocol or abstract.
                rvalue_type = get_proper_type(rvalue_type)
                lvalue_type = get_proper_type(lvalue_type)
                if (isinstance(rvalue_type, CallableType) and rvalue_type.is_type_obj() and
                        (rvalue_type.type_object().is_abstract or
                         rvalue_type.type_object().is_protocol) and
                        isinstance(lvalue_type, TypeType) and
                        isinstance(lvalue_type.item, Instance) and
                        (lvalue_type.item.type.is_abstract or
                         lvalue_type.item.type.is_protocol)):
                    self.msg.concrete_only_assign(lvalue_type, rvalue)
                    return
                if rvalue_type and infer_lvalue_type and not isinstance(lvalue_type, PartialType):
                    # Don't use type binder for definitions of special forms, like named tuples.
                    if not (isinstance(lvalue, NameExpr) and lvalue.is_special_form):
                        self.binder.assign_type(lvalue, rvalue_type, lvalue_type, False)

            elif index_lvalue:
                self.check_indexed_assignment(index_lvalue, rvalue, lvalue)

            if inferred:
                rvalue_type = self.expr_checker.accept(rvalue)
                if not inferred.is_final:
                    rvalue_type = remove_instance_last_known_values(rvalue_type)
                self.infer_variable_type(inferred, lvalue, rvalue_type, rvalue)
            self.check_assignment_to_slots(lvalue)

    # (type, operator) tuples for augmented assignments supported with partial types
    partial_type_augmented_ops: Final = {
        ('builtins.list', '+'),
        ('builtins.set', '|'),
    }

    def try_infer_partial_generic_type_from_assignment(self,
                                                       lvalue: Lvalue,
                                                       rvalue: Expression,
                                                       op: str) -> None:
        """Try to infer a precise type for partial generic type from assignment.

        'op' is '=' for normal assignment and a binary operator ('+', ...) for
        augmented assignment.

        Example where this happens:

            x = []
            if foo():
                x = [1]  # Infer List[int] as type of 'x'
        """
        var = None
        if (isinstance(lvalue, NameExpr)
                and isinstance(lvalue.node, Var)
                and isinstance(lvalue.node.type, PartialType)):
            var = lvalue.node
        elif isinstance(lvalue, MemberExpr):
            var = self.expr_checker.get_partial_self_var(lvalue)
        if var is not None:
            typ = var.type
            assert isinstance(typ, PartialType)
            if typ.type is None:
                return
            # Return if this is an unsupported augmented assignment.
            if op != '=' and (typ.type.fullname, op) not in self.partial_type_augmented_ops:
                return
            # TODO: some logic here duplicates the None partial type counterpart
            #       inlined in check_assignment(), see #8043.
            partial_types = self.find_partial_types(var)
            if partial_types is None:
                return
            rvalue_type = self.expr_checker.accept(rvalue)
            rvalue_type = get_proper_type(rvalue_type)
            if isinstance(rvalue_type, Instance):
                if rvalue_type.type == typ.type and is_valid_inferred_type(rvalue_type):
                    var.type = rvalue_type
                    del partial_types[var]
            elif isinstance(rvalue_type, AnyType):
                var.type = fill_typevars_with_any(typ.type)
                del partial_types[var]

    def check_compatibility_all_supers(self, lvalue: RefExpr, lvalue_type: Optional[Type],
                                       rvalue: Expression) -> bool:
        lvalue_node = lvalue.node
        # Check if we are a class variable with at least one base class
        if (isinstance(lvalue_node, Var) and
                lvalue.kind in (MDEF, None) and  # None for Vars defined via self
                len(lvalue_node.info.bases) > 0):

            for base in lvalue_node.info.mro[1:]:
                tnode = base.names.get(lvalue_node.name)
                if tnode is not None:
                    if not self.check_compatibility_classvar_super(lvalue_node,
                                                                   base,
                                                                   tnode.node):
                        # Show only one error per variable
                        break

                    if not self.check_compatibility_final_super(lvalue_node,
                                                                base,
                                                                tnode.node):
                        # Show only one error per variable
                        break

            direct_bases = lvalue_node.info.direct_base_classes()
            last_immediate_base = direct_bases[-1] if direct_bases else None

            for base in lvalue_node.info.mro[1:]:
                # Only check __slots__ against the 'object'
                # If a base class defines a Tuple of 3 elements, a child of
                # this class should not be allowed to define it as a Tuple of
                # anything other than 3 elements. The exception to this rule
                # is __slots__, where it is allowed for any child class to
                # redefine it.
                if lvalue_node.name == "__slots__" and base.fullname != "builtins.object":
                    continue
                # We don't care about the type of "__deletable__".
                if lvalue_node.name == "__deletable__":
                    continue

                if is_private(lvalue_node.name):
                    continue

                base_type, base_node = self.lvalue_type_from_base(lvalue_node, base)

                if base_type:
                    assert base_node is not None
                    if not self.check_compatibility_super(lvalue,
                                                          lvalue_type,
                                                          rvalue,
                                                          base,
                                                          base_type,
                                                          base_node):
                        # Only show one error per variable; even if other
                        # base classes are also incompatible
                        return True
                    if base is last_immediate_base:
                        # At this point, the attribute was found to be compatible with all
                        # immediate parents.
                        break
        return False

    def check_compatibility_super(self, lvalue: RefExpr, lvalue_type: Optional[Type],
                                  rvalue: Expression, base: TypeInfo, base_type: Type,
                                  base_node: Node) -> bool:
        lvalue_node = lvalue.node
        assert isinstance(lvalue_node, Var)

        # Do not check whether the rvalue is compatible if the
        # lvalue had a type defined; this is handled by other
        # parts, and all we have to worry about in that case is
        # that lvalue is compatible with the base class.
        compare_node = None
        if lvalue_type:
            compare_type = lvalue_type
            compare_node = lvalue.node
        else:
            compare_type = self.expr_checker.accept(rvalue, base_type)
            if isinstance(rvalue, NameExpr):
                compare_node = rvalue.node
                if isinstance(compare_node, Decorator):
                    compare_node = compare_node.func

        base_type = get_proper_type(base_type)
        compare_type = get_proper_type(compare_type)
        if compare_type:
            if (isinstance(base_type, CallableType) and
                    isinstance(compare_type, CallableType)):
                base_static = is_node_static(base_node)
                compare_static = is_node_static(compare_node)

                # In case compare_static is unknown, also check
                # if 'definition' is set. The most common case for
                # this is with TempNode(), where we lose all
                # information about the real rvalue node (but only get
                # the rvalue type)
                if compare_static is None and compare_type.definition:
                    compare_static = is_node_static(compare_type.definition)

                # Compare against False, as is_node_static can return None
                if base_static is False and compare_static is False:
                    # Class-level function objects and classmethods become bound
                    # methods: the former to the instance, the latter to the
                    # class
                    base_type = bind_self(base_type, self.scope.active_self_type())
                    compare_type = bind_self(compare_type, self.scope.active_self_type())

                # If we are a static method, ensure to also tell the
                # lvalue it now contains a static method
                if base_static and compare_static:
                    lvalue_node.is_staticmethod = True

            return self.check_subtype(compare_type, base_type, rvalue,
                                      message_registry.INCOMPATIBLE_TYPES_IN_ASSIGNMENT,
                                      'expression has type',
                                      'base class "%s" defined the type as' % base.name,
                                      code=codes.ASSIGNMENT)
        return True

    def lvalue_type_from_base(self, expr_node: Var,
                              base: TypeInfo) -> Tuple[Optional[Type], Optional[Node]]:
        """For a NameExpr that is part of a class, walk all base classes and try
        to find the first class that defines a Type for the same name."""
        expr_name = expr_node.name
        base_var = base.names.get(expr_name)

        if base_var:
            base_node = base_var.node
            base_type = base_var.type
            if isinstance(base_node, Decorator):
                base_node = base_node.func
                base_type = base_node.type

            if base_type:
                if not has_no_typevars(base_type):
                    self_type = self.scope.active_self_type()
                    assert self_type is not None, "Internal error: base lookup outside class"
                    if isinstance(self_type, TupleType):
                        instance = tuple_fallback(self_type)
                    else:
                        instance = self_type
                    itype = map_instance_to_supertype(instance, base)
                    base_type = expand_type_by_instance(base_type, itype)

                base_type = get_proper_type(base_type)
                if isinstance(base_type, CallableType) and isinstance(base_node, FuncDef):
                    # If we are a property, return the Type of the return
                    # value, not the Callable
                    if base_node.is_property:
                        base_type = get_proper_type(base_type.ret_type)
                if isinstance(base_type, FunctionLike) and isinstance(base_node,
                                                                      OverloadedFuncDef):
                    # Same for properties with setter
                    if base_node.is_property:
                        base_type = base_type.items[0].ret_type

                return base_type, base_node

        return None, None

    def check_compatibility_classvar_super(self, node: Var,
                                           base: TypeInfo, base_node: Optional[Node]) -> bool:
        if not isinstance(base_node, Var):
            return True
        if node.is_classvar and not base_node.is_classvar:
            self.fail(message_registry.CANNOT_OVERRIDE_INSTANCE_VAR.format(base.name), node)
            return False
        elif not node.is_classvar and base_node.is_classvar:
            self.fail(message_registry.CANNOT_OVERRIDE_CLASS_VAR.format(base.name), node)
            return False
        return True

    def check_compatibility_final_super(self, node: Var,
                                        base: TypeInfo, base_node: Optional[Node]) -> bool:
        """Check if an assignment overrides a final attribute in a base class.

        This only checks situations where either a node in base class is not a variable
        but a final method, or where override is explicitly declared as final.
        In these cases we give a more detailed error message. In addition, we check that
        a final variable doesn't override writeable attribute, which is not safe.

        Other situations are checked in `check_final()`.
        """
        if not isinstance(base_node, (Var, FuncBase, Decorator)):
            return True
        if base_node.is_final and (node.is_final or not isinstance(base_node, Var)):
            # Give this error only for explicit override attempt with `Final`, or
            # if we are overriding a final method with variable.
            # Other override attempts will be flagged as assignment to constant
            # in `check_final()`.
            self.msg.cant_override_final(node.name, base.name, node)
            return False
        if node.is_final:
            if base.fullname in ENUM_BASES and node.name in ENUM_SPECIAL_PROPS:
                return True
            self.check_if_final_var_override_writable(node.name, base_node, node)
        return True

    def check_if_final_var_override_writable(self,
                                             name: str,
                                             base_node: Optional[Node],
                                             ctx: Context) -> None:
        """Check that a final variable doesn't override writeable attribute.

        This is done to prevent situations like this:
            class C:
                attr = 1
            class D(C):
                attr: Final = 2

            x: C = D()
            x.attr = 3  # Oops!
        """
        writable = True
        if base_node:
            writable = self.is_writable_attribute(base_node)
        if writable:
            self.msg.final_cant_override_writable(name, ctx)

    def get_final_context(self) -> bool:
        """Check whether we a currently checking a final declaration."""
        return self._is_final_def

    @contextmanager
    def enter_final_context(self, is_final_def: bool) -> Iterator[None]:
        """Store whether the current checked assignment is a final declaration."""
        old_ctx = self._is_final_def
        self._is_final_def = is_final_def
        try:
            yield
        finally:
            self._is_final_def = old_ctx

    def check_final(self,
                    s: Union[AssignmentStmt, OperatorAssignmentStmt, AssignmentExpr]) -> None:
        """Check if this assignment does not assign to a final attribute.

        This function performs the check only for name assignments at module
        and class scope. The assignments to `obj.attr` and `Cls.attr` are checked
        in checkmember.py.
        """
        if isinstance(s, AssignmentStmt):
            lvs = self.flatten_lvalues(s.lvalues)
        elif isinstance(s, AssignmentExpr):
            lvs = [s.target]
        else:
            lvs = [s.lvalue]
        is_final_decl = s.is_final_def if isinstance(s, AssignmentStmt) else False
        if is_final_decl and self.scope.active_class():
            lv = lvs[0]
            assert isinstance(lv, RefExpr)
            assert isinstance(lv.node, Var)
            if (lv.node.final_unset_in_class and not lv.node.final_set_in_init and
                    not self.is_stub and  # It is OK to skip initializer in stub files.
                    # Avoid extra error messages, if there is no type in Final[...],
                    # then we already reported the error about missing r.h.s.
                    isinstance(s, AssignmentStmt) and s.type is not None):
                self.msg.final_without_value(s)
        for lv in lvs:
            if isinstance(lv, RefExpr) and isinstance(lv.node, Var):
                name = lv.node.name
                cls = self.scope.active_class()
                if cls is not None:
                    # These additional checks exist to give more error messages
                    # even if the final attribute was overridden with a new symbol
                    # (which is itself an error)...
                    for base in cls.mro[1:]:
                        sym = base.names.get(name)
                        # We only give this error if base node is variable,
                        # overriding final method will be caught in
                        # `check_compatibility_final_super()`.
                        if sym and isinstance(sym.node, Var):
                            if sym.node.is_final and not is_final_decl:
                                self.msg.cant_assign_to_final(name, sym.node.info is None, s)
                                # ...but only once
                                break
                if lv.node.is_final and not is_final_decl:
                    self.msg.cant_assign_to_final(name, lv.node.info is None, s)

    def check_assignment_to_slots(self, lvalue: Lvalue) -> None:
        if not isinstance(lvalue, MemberExpr):
            return

        inst = get_proper_type(self.expr_checker.accept(lvalue.expr))
        if not isinstance(inst, Instance):
            return
        if inst.type.slots is None:
            return  # Slots do not exist, we can allow any assignment
        if lvalue.name in inst.type.slots:
            return  # We are assigning to an existing slot
        for base_info in inst.type.mro[:-1]:
            if base_info.names.get('__setattr__') is not None:
                # When type has `__setattr__` defined,
                # we can assign any dynamic value.
                # We exclude object, because it always has `__setattr__`.
                return

        definition = inst.type.get(lvalue.name)
        if definition is None:
            # We don't want to duplicate
            # `"SomeType" has no attribute "some_attr"`
            # error twice.
            return
        if self.is_assignable_slot(lvalue, definition.type):
            return

        self.fail(
            message_registry.NAME_NOT_IN_SLOTS.format(
                lvalue.name, inst.type.fullname,
            ),
            lvalue,
        )

    def is_assignable_slot(self, lvalue: Lvalue, typ: Optional[Type]) -> bool:
        if getattr(lvalue, 'node', None):
            return False  # This is a definition

        typ = get_proper_type(typ)
        if typ is None or isinstance(typ, AnyType):
            return True  # Any can be literally anything, like `@propery`
        if isinstance(typ, Instance):
            # When working with instances, we need to know if they contain
            # `__set__` special method. Like `@property` does.
            # This makes assigning to properties possible,
            # even without extra slot spec.
            return typ.type.get('__set__') is not None
        if isinstance(typ, FunctionLike):
            return True  # Can be a property, or some other magic
        if isinstance(typ, UnionType):
            return all(self.is_assignable_slot(lvalue, u) for u in typ.items)
        return False

    def check_assignment_to_multiple_lvalues(self, lvalues: List[Lvalue], rvalue: Expression,
                                             context: Context,
                                             infer_lvalue_type: bool = True) -> None:
        if isinstance(rvalue, TupleExpr) or isinstance(rvalue, ListExpr):
            # Recursively go into Tuple or List expression rhs instead of
            # using the type of rhs, because this allowed more fine grained
            # control in cases like: a, b = [int, str] where rhs would get
            # type List[object]
            rvalues: List[Expression] = []
            iterable_type: Optional[Type] = None
            last_idx: Optional[int] = None
            for idx_rval, rval in enumerate(rvalue.items):
                if isinstance(rval, StarExpr):
                    typs = get_proper_type(self.expr_checker.visit_star_expr(rval).type)
                    if isinstance(typs, TupleType):
                        rvalues.extend([TempNode(typ) for typ in typs.items])
                    elif self.type_is_iterable(typs) and isinstance(typs, Instance):
                        if (iterable_type is not None
                                and iterable_type != self.iterable_item_type(typs)):
                            self.fail(message_registry.CONTIGUOUS_ITERABLE_EXPECTED, context)
                        else:
                            if last_idx is None or last_idx + 1 == idx_rval:
                                rvalues.append(rval)
                                last_idx = idx_rval
                                iterable_type = self.iterable_item_type(typs)
                            else:
                                self.fail(message_registry.CONTIGUOUS_ITERABLE_EXPECTED, context)
                    else:
                        self.fail(message_registry.ITERABLE_TYPE_EXPECTED.format(typs),
                             context)
                else:
                    rvalues.append(rval)
            iterable_start: Optional[int] = None
            iterable_end: Optional[int] = None
            for i, rval in enumerate(rvalues):
                if isinstance(rval, StarExpr):
                    typs = get_proper_type(self.expr_checker.visit_star_expr(rval).type)
                    if self.type_is_iterable(typs) and isinstance(typs, Instance):
                        if iterable_start is None:
                            iterable_start = i
                        iterable_end = i
            if (iterable_start is not None
                    and iterable_end is not None
                    and iterable_type is not None):
                iterable_num = iterable_end - iterable_start + 1
                rvalue_needed = len(lvalues) - (len(rvalues) - iterable_num)
                if rvalue_needed > 0:
                    rvalues = rvalues[0: iterable_start] + [TempNode(iterable_type)
                        for i in range(rvalue_needed)] + rvalues[iterable_end + 1:]

            if self.check_rvalue_count_in_assignment(lvalues, len(rvalues), context):
                star_index = next((i for i, lv in enumerate(lvalues) if
                                   isinstance(lv, StarExpr)), len(lvalues))

                left_lvs = lvalues[:star_index]
                star_lv = cast(StarExpr,
                               lvalues[star_index]) if star_index != len(lvalues) else None
                right_lvs = lvalues[star_index + 1:]

                left_rvs, star_rvs, right_rvs = self.split_around_star(
                    rvalues, star_index, len(lvalues))

                lr_pairs = list(zip(left_lvs, left_rvs))
                if star_lv:
                    rv_list = ListExpr(star_rvs)
                    rv_list.set_line(rvalue.get_line())
                    lr_pairs.append((star_lv.expr, rv_list))
                lr_pairs.extend(zip(right_lvs, right_rvs))

                for lv, rv in lr_pairs:
                    self.check_assignment(lv, rv, infer_lvalue_type)
        else:
            self.check_multi_assignment(lvalues, rvalue, context, infer_lvalue_type)

    def check_rvalue_count_in_assignment(self, lvalues: List[Lvalue], rvalue_count: int,
                                         context: Context) -> bool:
        if any(isinstance(lvalue, StarExpr) for lvalue in lvalues):
            if len(lvalues) - 1 > rvalue_count:
                self.msg.wrong_number_values_to_unpack(rvalue_count,
                                                       len(lvalues) - 1, context)
                return False
        elif rvalue_count != len(lvalues):
            self.msg.wrong_number_values_to_unpack(rvalue_count, len(lvalues), context)
            return False
        return True

    def check_multi_assignment(self, lvalues: List[Lvalue],
                               rvalue: Expression,
                               context: Context,
                               infer_lvalue_type: bool = True,
                               rv_type: Optional[Type] = None,
                               undefined_rvalue: bool = False) -> None:
        """Check the assignment of one rvalue to a number of lvalues."""

        # Infer the type of an ordinary rvalue expression.
        # TODO: maybe elsewhere; redundant.
        rvalue_type = get_proper_type(rv_type or self.expr_checker.accept(rvalue))

        if isinstance(rvalue_type, UnionType):
            # If this is an Optional type in non-strict Optional code, unwrap it.
            relevant_items = rvalue_type.relevant_items()
            if len(relevant_items) == 1:
                rvalue_type = get_proper_type(relevant_items[0])

        if isinstance(rvalue_type, AnyType):
            for lv in lvalues:
                if isinstance(lv, StarExpr):
                    lv = lv.expr
                temp_node = self.temp_node(AnyType(TypeOfAny.from_another_any,
                                                   source_any=rvalue_type), context)
                self.check_assignment(lv, temp_node, infer_lvalue_type)
        elif isinstance(rvalue_type, TupleType):
            self.check_multi_assignment_from_tuple(lvalues, rvalue, rvalue_type,
                                                   context, undefined_rvalue, infer_lvalue_type)
        elif isinstance(rvalue_type, UnionType):
            self.check_multi_assignment_from_union(lvalues, rvalue, rvalue_type, context,
                                                   infer_lvalue_type)
        elif isinstance(rvalue_type, Instance) and rvalue_type.type.fullname == 'builtins.str':
            self.msg.unpacking_strings_disallowed(context)
        else:
            self.check_multi_assignment_from_iterable(lvalues, rvalue_type,
                                                      context, infer_lvalue_type)

    def check_multi_assignment_from_union(self, lvalues: List[Expression], rvalue: Expression,
                                          rvalue_type: UnionType, context: Context,
                                          infer_lvalue_type: bool) -> None:
        """Check assignment to multiple lvalue targets when rvalue type is a Union[...].
        For example:

            t: Union[Tuple[int, int], Tuple[str, str]]
            x, y = t
            reveal_type(x)  # Union[int, str]

        The idea in this case is to process the assignment for every item of the union.
        Important note: the types are collected in two places, 'union_types' contains
        inferred types for first assignments, 'assignments' contains the narrowed types
        for binder.
        """
        self.no_partial_types = True
        transposed: Tuple[List[Type], ...] = tuple([] for _ in self.flatten_lvalues(lvalues))
        # Notify binder that we want to defer bindings and instead collect types.
        with self.binder.accumulate_type_assignments() as assignments:
            for item in rvalue_type.items:
                # Type check the assignment separately for each union item and collect
                # the inferred lvalue types for each union item.
                self.check_multi_assignment(lvalues, rvalue, context,
                                            infer_lvalue_type=infer_lvalue_type,
                                            rv_type=item, undefined_rvalue=True)
                for t, lv in zip(transposed, self.flatten_lvalues(lvalues)):
                    t.append(self.type_map.pop(lv, AnyType(TypeOfAny.special_form)))
        union_types = tuple(make_simplified_union(col) for col in transposed)
        for expr, items in assignments.items():
            # Bind a union of types collected in 'assignments' to every expression.
            if isinstance(expr, StarExpr):
                expr = expr.expr

            # TODO: See todo in binder.py, ConditionalTypeBinder.assign_type
            # It's unclear why the 'declared_type' param is sometimes 'None'
            clean_items: List[Tuple[Type, Type]] = []
            for type, declared_type in items:
                assert declared_type is not None
                clean_items.append((type, declared_type))

            # TODO: fix signature of zip() in typeshed.
            types, declared_types = cast(Any, zip)(*clean_items)
            self.binder.assign_type(expr,
                                    make_simplified_union(list(types)),
                                    make_simplified_union(list(declared_types)),
                                    False)
        for union, lv in zip(union_types, self.flatten_lvalues(lvalues)):
            # Properly store the inferred types.
            _1, _2, inferred = self.check_lvalue(lv)
            if inferred:
                self.set_inferred_type(inferred, lv, union)
            else:
                self.store_type(lv, union)
        self.no_partial_types = False

    def flatten_lvalues(self, lvalues: List[Expression]) -> List[Expression]:
        res: List[Expression] = []
        for lv in lvalues:
            if isinstance(lv, (TupleExpr, ListExpr)):
                res.extend(self.flatten_lvalues(lv.items))
            if isinstance(lv, StarExpr):
                # Unwrap StarExpr, since it is unwrapped by other helpers.
                lv = lv.expr
            res.append(lv)
        return res

    def check_multi_assignment_from_tuple(self, lvalues: List[Lvalue], rvalue: Expression,
                                          rvalue_type: TupleType, context: Context,
                                          undefined_rvalue: bool,
                                          infer_lvalue_type: bool = True) -> None:
        if self.check_rvalue_count_in_assignment(lvalues, len(rvalue_type.items), context):
            star_index = next((i for i, lv in enumerate(lvalues)
                               if isinstance(lv, StarExpr)), len(lvalues))

            left_lvs = lvalues[:star_index]
            star_lv = cast(StarExpr, lvalues[star_index]) if star_index != len(lvalues) else None
            right_lvs = lvalues[star_index + 1:]

            if not undefined_rvalue:
                # Infer rvalue again, now in the correct type context.
                lvalue_type = self.lvalue_type_for_inference(lvalues, rvalue_type)
                reinferred_rvalue_type = get_proper_type(self.expr_checker.accept(rvalue,
                                                                                  lvalue_type))

                if isinstance(reinferred_rvalue_type, UnionType):
                    # If this is an Optional type in non-strict Optional code, unwrap it.
                    relevant_items = reinferred_rvalue_type.relevant_items()
                    if len(relevant_items) == 1:
                        reinferred_rvalue_type = get_proper_type(relevant_items[0])
                if isinstance(reinferred_rvalue_type, UnionType):
                    self.check_multi_assignment_from_union(lvalues, rvalue,
                                                           reinferred_rvalue_type, context,
                                                           infer_lvalue_type)
                    return
                if isinstance(reinferred_rvalue_type, AnyType):
                    # We can get Any if the current node is
                    # deferred. Doing more inference in deferred nodes
                    # is hard, so give up for now.  We can also get
                    # here if reinferring types above changes the
                    # inferred return type for an overloaded function
                    # to be ambiguous.
                    return
                assert isinstance(reinferred_rvalue_type, TupleType)
                rvalue_type = reinferred_rvalue_type

            left_rv_types, star_rv_types, right_rv_types = self.split_around_star(
                rvalue_type.items, star_index, len(lvalues))

            for lv, rv_type in zip(left_lvs, left_rv_types):
                self.check_assignment(lv, self.temp_node(rv_type, context), infer_lvalue_type)
            if star_lv:
                list_expr = ListExpr([self.temp_node(rv_type, context)
                                      for rv_type in star_rv_types])
                list_expr.set_line(context.get_line())
                self.check_assignment(star_lv.expr, list_expr, infer_lvalue_type)
            for lv, rv_type in zip(right_lvs, right_rv_types):
                self.check_assignment(lv, self.temp_node(rv_type, context), infer_lvalue_type)

    def lvalue_type_for_inference(self, lvalues: List[Lvalue], rvalue_type: TupleType) -> Type:
        star_index = next((i for i, lv in enumerate(lvalues)
                           if isinstance(lv, StarExpr)), len(lvalues))
        left_lvs = lvalues[:star_index]
        star_lv = cast(StarExpr, lvalues[star_index]) if star_index != len(lvalues) else None
        right_lvs = lvalues[star_index + 1:]
        left_rv_types, star_rv_types, right_rv_types = self.split_around_star(
            rvalue_type.items, star_index, len(lvalues))

        type_parameters: List[Type] = []

        def append_types_for_inference(lvs: List[Expression], rv_types: List[Type]) -> None:
            for lv, rv_type in zip(lvs, rv_types):
                sub_lvalue_type, index_expr, inferred = self.check_lvalue(lv)
                if sub_lvalue_type and not isinstance(sub_lvalue_type, PartialType):
                    type_parameters.append(sub_lvalue_type)
                else:  # index lvalue
                    # TODO Figure out more precise type context, probably
                    #      based on the type signature of the _set method.
                    type_parameters.append(rv_type)

        append_types_for_inference(left_lvs, left_rv_types)

        if star_lv:
            sub_lvalue_type, index_expr, inferred = self.check_lvalue(star_lv.expr)
            if sub_lvalue_type and not isinstance(sub_lvalue_type, PartialType):
                type_parameters.extend([sub_lvalue_type] * len(star_rv_types))
            else:  # index lvalue
                # TODO Figure out more precise type context, probably
                #      based on the type signature of the _set method.
                type_parameters.extend(star_rv_types)

        append_types_for_inference(right_lvs, right_rv_types)

        return TupleType(type_parameters, self.named_type('builtins.tuple'))

    def split_around_star(self, items: List[T], star_index: int,
                          length: int) -> Tuple[List[T], List[T], List[T]]:
        """Splits a list of items in three to match another list of length 'length'
        that contains a starred expression at 'star_index' in the following way:

        star_index = 2, length = 5 (i.e., [a,b,*,c,d]), items = [1,2,3,4,5,6,7]
        returns in: ([1,2], [3,4,5], [6,7])
        """
        nr_right_of_star = length - star_index - 1
        right_index = -nr_right_of_star if nr_right_of_star != 0 else len(items)
        left = items[:star_index]
        star = items[star_index:right_index]
        right = items[right_index:]
        return left, star, right

    def type_is_iterable(self, type: Type) -> bool:
        type = get_proper_type(type)
        if isinstance(type, CallableType) and type.is_type_obj():
            type = type.fallback
        return is_subtype(type, self.named_generic_type('typing.Iterable',
                                                        [AnyType(TypeOfAny.special_form)]))

    def check_multi_assignment_from_iterable(self, lvalues: List[Lvalue], rvalue_type: Type,
                                             context: Context,
                                             infer_lvalue_type: bool = True) -> None:
        rvalue_type = get_proper_type(rvalue_type)
        if self.type_is_iterable(rvalue_type) and isinstance(rvalue_type, Instance):
            item_type = self.iterable_item_type(rvalue_type)
            for lv in lvalues:
                if isinstance(lv, StarExpr):
                    items_type = self.named_generic_type('builtins.list', [item_type])
                    self.check_assignment(lv.expr, self.temp_node(items_type, context),
                                          infer_lvalue_type)
                else:
                    self.check_assignment(lv, self.temp_node(item_type, context),
                                          infer_lvalue_type)
        else:
            self.msg.type_not_iterable(rvalue_type, context)

    def check_lvalue(self, lvalue: Lvalue) -> Tuple[Optional[Type],
                                                    Optional[IndexExpr],
                                                    Optional[Var]]:
        lvalue_type = None
        index_lvalue = None
        inferred = None

        if self.is_definition(lvalue) and (
            not isinstance(lvalue, NameExpr) or isinstance(lvalue.node, Var)
        ):
            if isinstance(lvalue, NameExpr):
                inferred = cast(Var, lvalue.node)
                assert isinstance(inferred, Var)
            else:
                assert isinstance(lvalue, MemberExpr)
                self.expr_checker.accept(lvalue.expr)
                inferred = lvalue.def_var
        elif isinstance(lvalue, IndexExpr):
            index_lvalue = lvalue
        elif isinstance(lvalue, MemberExpr):
            lvalue_type = self.expr_checker.analyze_ordinary_member_access(lvalue, True)
            self.store_type(lvalue, lvalue_type)
        elif isinstance(lvalue, NameExpr):
            lvalue_type = self.expr_checker.analyze_ref_expr(lvalue, lvalue=True)
            self.store_type(lvalue, lvalue_type)
        elif isinstance(lvalue, TupleExpr) or isinstance(lvalue, ListExpr):
            types = [self.check_lvalue(sub_expr)[0] or
                     # This type will be used as a context for further inference of rvalue,
                     # we put Uninhabited if there is no information available from lvalue.
                     UninhabitedType() for sub_expr in lvalue.items]
            lvalue_type = TupleType(types, self.named_type('builtins.tuple'))
        elif isinstance(lvalue, StarExpr):
            typ, _, _ = self.check_lvalue(lvalue.expr)
            lvalue_type = StarType(typ) if typ else None
        else:
            lvalue_type = self.expr_checker.accept(lvalue)

        return lvalue_type, index_lvalue, inferred

    def is_definition(self, s: Lvalue) -> bool:
        if isinstance(s, NameExpr):
            if s.is_inferred_def:
                return True
            # If the node type is not defined, this must the first assignment
            # that we process => this is a definition, even though the semantic
            # analyzer did not recognize this as such. This can arise in code
            # that uses isinstance checks, if type checking of the primary
            # definition is skipped due to an always False type check.
            node = s.node
            if isinstance(node, Var):
                return node.type is None
        elif isinstance(s, MemberExpr):
            return s.is_inferred_def
        return False

    def infer_variable_type(self, name: Var, lvalue: Lvalue,
                            init_type: Type, context: Context) -> None:
        """Infer the type of initialized variables from initializer type."""
        init_type = get_proper_type(init_type)
        if isinstance(init_type, DeletedType):
            self.msg.deleted_as_rvalue(init_type, context)
        elif not is_valid_inferred_type(init_type) and not self.no_partial_types:
            # We cannot use the type of the initialization expression for full type
            # inference (it's not specific enough), but we might be able to give
            # partial type which will be made more specific later. A partial type
            # gets generated in assignment like 'x = []' where item type is not known.
            if not self.infer_partial_type(name, lvalue, init_type):
                self.msg.need_annotation_for_var(name, context, self.options.python_version)
                self.set_inference_error_fallback_type(name, lvalue, init_type)
        elif (isinstance(lvalue, MemberExpr) and self.inferred_attribute_types is not None
              and lvalue.def_var and lvalue.def_var in self.inferred_attribute_types
              and not is_same_type(self.inferred_attribute_types[lvalue.def_var], init_type)):
            # Multiple, inconsistent types inferred for an attribute.
            self.msg.need_annotation_for_var(name, context, self.options.python_version)
            name.type = AnyType(TypeOfAny.from_error)
        else:
            # Infer type of the target.

            # Make the type more general (strip away function names etc.).
            init_type = strip_type(init_type)

            self.set_inferred_type(name, lvalue, init_type)

    def infer_partial_type(self, name: Var, lvalue: Lvalue, init_type: Type) -> bool:
        init_type = get_proper_type(init_type)
        if isinstance(init_type, NoneType):
            partial_type = PartialType(None, name)
        elif isinstance(init_type, Instance):
            fullname = init_type.type.fullname
            is_ref = isinstance(lvalue, RefExpr)
            if (is_ref and
                    (fullname == 'builtins.list' or
                     fullname == 'builtins.set' or
                     fullname == 'builtins.dict' or
                     fullname == 'collections.OrderedDict') and
                    all(isinstance(t, (NoneType, UninhabitedType))
                        for t in get_proper_types(init_type.args))):
                partial_type = PartialType(init_type.type, name)
            elif is_ref and fullname == 'collections.defaultdict':
                arg0 = get_proper_type(init_type.args[0])
                arg1 = get_proper_type(init_type.args[1])
                if (isinstance(arg0, (NoneType, UninhabitedType)) and
                        self.is_valid_defaultdict_partial_value_type(arg1)):
                    arg1 = erase_type(arg1)
                    assert isinstance(arg1, Instance)
                    partial_type = PartialType(init_type.type, name, arg1)
                else:
                    return False
            else:
                return False
        else:
            return False
        self.set_inferred_type(name, lvalue, partial_type)
        self.partial_types[-1].map[name] = lvalue
        return True

    def is_valid_defaultdict_partial_value_type(self, t: ProperType) -> bool:
        """Check if t can be used as the basis for a partial defaultdict value type.

        Examples:

          * t is 'int' --> True
          * t is 'list[<nothing>]' --> True
          * t is 'dict[...]' --> False (only generic types with a single type
            argument supported)
        """
        if not isinstance(t, Instance):
            return False
        if len(t.args) == 0:
            return True
        if len(t.args) == 1:
            arg = get_proper_type(t.args[0])
            # TODO: This is too permissive -- we only allow TypeVarType since
            #       they leak in cases like defaultdict(list) due to a bug.
            #       This can result in incorrect types being inferred, but only
            #       in rare cases.
            if isinstance(arg, (TypeVarType, UninhabitedType, NoneType)):
                return True
        return False

    def set_inferred_type(self, var: Var, lvalue: Lvalue, type: Type) -> None:
        """Store inferred variable type.

        Store the type to both the variable node and the expression node that
        refers to the variable (lvalue). If var is None, do nothing.
        """
        if var and not self.current_node_deferred:
            var.type = type
            var.is_inferred = True
            if var not in self.var_decl_frames:
                # Used for the hack to improve optional type inference in conditionals
                self.var_decl_frames[var] = {frame.id for frame in self.binder.frames}
            if isinstance(lvalue, MemberExpr) and self.inferred_attribute_types is not None:
                # Store inferred attribute type so that we can check consistency afterwards.
                if lvalue.def_var is not None:
                    self.inferred_attribute_types[lvalue.def_var] = type
            self.store_type(lvalue, type)

    def set_inference_error_fallback_type(self, var: Var, lvalue: Lvalue, type: Type) -> None:
        """Store best known type for variable if type inference failed.

        If a program ignores error on type inference error, the variable should get some
        inferred type so that if can used later on in the program. Example:

          x = []  # type: ignore
          x.append(1)   # Should be ok!

        We implement this here by giving x a valid type (replacing inferred <nothing> with Any).
        """
        fallback = self.inference_error_fallback_type(type)
        self.set_inferred_type(var, lvalue, fallback)

    def inference_error_fallback_type(self, type: Type) -> Type:
        fallback = type.accept(SetNothingToAny())
        # Type variables may leak from inference, see https://github.com/python/mypy/issues/5738,
        # we therefore need to erase them.
        return erase_typevars(fallback)

    def check_simple_assignment(self, lvalue_type: Optional[Type], rvalue: Expression,
                                context: Context,
                                msg: str = message_registry.INCOMPATIBLE_TYPES_IN_ASSIGNMENT,
                                lvalue_name: str = 'variable',
                                rvalue_name: str = 'expression', *,
                                code: Optional[ErrorCode] = None) -> Type:
        if self.is_stub and isinstance(rvalue, EllipsisExpr):
            # '...' is always a valid initializer in a stub.
            return AnyType(TypeOfAny.special_form)
        else:
            lvalue_type = get_proper_type(lvalue_type)
            always_allow_any = lvalue_type is not None and not isinstance(lvalue_type, AnyType)
            rvalue_type = self.expr_checker.accept(rvalue, lvalue_type,
                                                   always_allow_any=always_allow_any)
            rvalue_type = get_proper_type(rvalue_type)
            if isinstance(rvalue_type, DeletedType):
                self.msg.deleted_as_rvalue(rvalue_type, context)
            if isinstance(lvalue_type, DeletedType):
                self.msg.deleted_as_lvalue(lvalue_type, context)
            elif lvalue_type:
                self.check_subtype(rvalue_type, lvalue_type, context, msg,
                                   '{} has type'.format(rvalue_name),
                                   '{} has type'.format(lvalue_name), code=code)
            return rvalue_type

    def check_member_assignment(self, instance_type: Type, attribute_type: Type,
                                rvalue: Expression, context: Context) -> Tuple[Type, Type, bool]:
        """Type member assignment.

        This defers to check_simple_assignment, unless the member expression
        is a descriptor, in which case this checks descriptor semantics as well.

        Return the inferred rvalue_type, inferred lvalue_type, and whether to use the binder
        for this assignment.

        Note: this method exists here and not in checkmember.py, because we need to take
        care about interaction between binder and __set__().
        """
        instance_type = get_proper_type(instance_type)
        attribute_type = get_proper_type(attribute_type)
        # Descriptors don't participate in class-attribute access
        if ((isinstance(instance_type, FunctionLike) and instance_type.is_type_obj()) or
                isinstance(instance_type, TypeType)):
            rvalue_type = self.check_simple_assignment(attribute_type, rvalue, context,
                                                       code=codes.ASSIGNMENT)
            return rvalue_type, attribute_type, True

        if not isinstance(attribute_type, Instance):
            # TODO: support __set__() for union types.
            rvalue_type = self.check_simple_assignment(attribute_type, rvalue, context,
                                                       code=codes.ASSIGNMENT)
            return rvalue_type, attribute_type, True

        mx = MemberContext(
            is_lvalue=False, is_super=False, is_operator=False,
            original_type=instance_type, context=context, self_type=None,
            msg=self.msg, chk=self,
        )
        get_type = analyze_descriptor_access(attribute_type, mx)
        if not attribute_type.type.has_readable_member('__set__'):
            # If there is no __set__, we type-check that the assigned value matches
            # the return type of __get__. This doesn't match the python semantics,
            # (which allow you to override the descriptor with any value), but preserves
            # the type of accessing the attribute (even after the override).
            rvalue_type = self.check_simple_assignment(get_type, rvalue, context,
                                                       code=codes.ASSIGNMENT)
            return rvalue_type, get_type, True

        dunder_set = attribute_type.type.get_method('__set__')
        if dunder_set is None:
            self.fail(message_registry.DESCRIPTOR_SET_NOT_CALLABLE.format(attribute_type), context)
            return AnyType(TypeOfAny.from_error), get_type, False

        bound_method = analyze_decorator_or_funcbase_access(
            defn=dunder_set, itype=attribute_type, info=attribute_type.type,
            self_type=attribute_type, name='__set__', mx=mx)
        typ = map_instance_to_supertype(attribute_type, dunder_set.info)
        dunder_set_type = expand_type_by_instance(bound_method, typ)

        callable_name = self.expr_checker.method_fullname(attribute_type, "__set__")
        dunder_set_type = self.expr_checker.transform_callee_type(
            callable_name, dunder_set_type,
            [TempNode(instance_type, context=context), rvalue],
            [nodes.ARG_POS, nodes.ARG_POS],
            context, object_type=attribute_type,
        )

        # For non-overloaded setters, the result should be type-checked like a regular assignment.
        # Hence, we first only try to infer the type by using the rvalue as type context.
        type_context = rvalue
        with self.msg.disable_errors():
            _, inferred_dunder_set_type = self.expr_checker.check_call(
                dunder_set_type,
                [TempNode(instance_type, context=context), type_context],
                [nodes.ARG_POS, nodes.ARG_POS],
                context, object_type=attribute_type,
                callable_name=callable_name)

        # And now we in fact type check the call, to show errors related to wrong arguments
        # count, etc., replacing the type context for non-overloaded setters only.
        inferred_dunder_set_type = get_proper_type(inferred_dunder_set_type)
        if isinstance(inferred_dunder_set_type, CallableType):
            type_context = TempNode(AnyType(TypeOfAny.special_form), context=context)
        self.expr_checker.check_call(
            dunder_set_type,
            [TempNode(instance_type, context=context), type_context],
            [nodes.ARG_POS, nodes.ARG_POS],
            context, object_type=attribute_type,
            callable_name=callable_name)

        # In the following cases, a message already will have been recorded in check_call.
        if ((not isinstance(inferred_dunder_set_type, CallableType)) or
                (len(inferred_dunder_set_type.arg_types) < 2)):
            return AnyType(TypeOfAny.from_error), get_type, False

        set_type = inferred_dunder_set_type.arg_types[1]
        # Special case: if the rvalue_type is a subtype of both '__get__' and '__set__' types,
        # and '__get__' type is narrower than '__set__', then we invoke the binder to narrow type
        # by this assignment. Technically, this is not safe, but in practice this is
        # what a user expects.
        rvalue_type = self.check_simple_assignment(set_type, rvalue, context,
                                                   code=codes.ASSIGNMENT)
        infer = is_subtype(rvalue_type, get_type) and is_subtype(get_type, set_type)
        return rvalue_type if infer else set_type, get_type, infer

    def check_indexed_assignment(self, lvalue: IndexExpr,
                                 rvalue: Expression, context: Context) -> None:
        """Type check indexed assignment base[index] = rvalue.

        The lvalue argument is the base[index] expression.
        """
        self.try_infer_partial_type_from_indexed_assignment(lvalue, rvalue)
        basetype = get_proper_type(self.expr_checker.accept(lvalue.base))
        method_type = self.expr_checker.analyze_external_member_access(
            '__setitem__', basetype, lvalue)

        lvalue.method_type = method_type
        self.expr_checker.check_method_call(
            '__setitem__', basetype, method_type, [lvalue.index, rvalue],
            [nodes.ARG_POS, nodes.ARG_POS], context)

    def try_infer_partial_type_from_indexed_assignment(
            self, lvalue: IndexExpr, rvalue: Expression) -> None:
        # TODO: Should we share some of this with try_infer_partial_type?
        var = None
        if isinstance(lvalue.base, RefExpr) and isinstance(lvalue.base.node, Var):
            var = lvalue.base.node
        elif isinstance(lvalue.base, MemberExpr):
            var = self.expr_checker.get_partial_self_var(lvalue.base)
        if isinstance(var, Var):
            if isinstance(var.type, PartialType):
                type_type = var.type.type
                if type_type is None:
                    return  # The partial type is None.
                partial_types = self.find_partial_types(var)
                if partial_types is None:
                    return
                typename = type_type.fullname
                if (typename == 'builtins.dict'
                        or typename == 'collections.OrderedDict'
                        or typename == 'collections.defaultdict'):
                    # TODO: Don't infer things twice.
                    key_type = self.expr_checker.accept(lvalue.index)
                    value_type = self.expr_checker.accept(rvalue)
                    if (is_valid_inferred_type(key_type) and
                            is_valid_inferred_type(value_type) and
                            not self.current_node_deferred and
                            not (typename == 'collections.defaultdict' and
                                 var.type.value_type is not None and
                                 not is_equivalent(value_type, var.type.value_type))):
                        var.type = self.named_generic_type(typename,
                                                           [key_type, value_type])
                        del partial_types[var]

    def visit_expression_stmt(self, s: ExpressionStmt) -> None:
        self.expr_checker.accept(s.expr, allow_none_return=True, always_allow_any=True)

    def visit_return_stmt(self, s: ReturnStmt) -> None:
        """Type check a return statement."""
        self.check_return_stmt(s)
        self.binder.unreachable()

    def check_return_stmt(self, s: ReturnStmt) -> None:
        defn = self.scope.top_function()
        if defn is not None:
            if defn.is_generator:
                return_type = self.get_generator_return_type(self.return_types[-1],
                                                             defn.is_coroutine)
            elif defn.is_coroutine:
                return_type = self.get_coroutine_return_type(self.return_types[-1])
            else:
                return_type = self.return_types[-1]
            return_type = get_proper_type(return_type)

            if isinstance(return_type, UninhabitedType):
                self.fail(message_registry.NO_RETURN_EXPECTED, s)
                return

            if s.expr:
                is_lambda = isinstance(self.scope.top_function(), LambdaExpr)
                declared_none_return = isinstance(return_type, NoneType)
                declared_any_return = isinstance(return_type, AnyType)

                # This controls whether or not we allow a function call that
                # returns None as the expression of this return statement.
                # E.g. `return f()` for some `f` that returns None.  We allow
                # this only if we're in a lambda or in a function that returns
                # `None` or `Any`.
                allow_none_func_call = is_lambda or declared_none_return or declared_any_return

                # Return with a value.
                typ = get_proper_type(self.expr_checker.accept(
                    s.expr, return_type, allow_none_return=allow_none_func_call))

                if defn.is_async_generator:
                    self.fail(message_registry.RETURN_IN_ASYNC_GENERATOR, s)
                    return
                # Returning a value of type Any is always fine.
                if isinstance(typ, AnyType):
                    # (Unless you asked to be warned in that case, and the
                    # function is not declared to return Any)
                    if (self.options.warn_return_any
                        and not self.current_node_deferred
                        and not is_proper_subtype(AnyType(TypeOfAny.special_form), return_type)
                        and not (defn.name in BINARY_MAGIC_METHODS and
                                 is_literal_not_implemented(s.expr))
                        and not (isinstance(return_type, Instance) and
                                 return_type.type.fullname == 'builtins.object')):
                        self.msg.incorrectly_returning_any(return_type, s)
                    return

                # Disallow return expressions in functions declared to return
                # None, subject to two exceptions below.
                if declared_none_return:
                    # Lambdas are allowed to have None returns.
                    # Functions returning a value of type None are allowed to have a None return.
                    if is_lambda or isinstance(typ, NoneType):
                        return
                    self.fail(message_registry.NO_RETURN_VALUE_EXPECTED, s)
                else:
                    self.check_subtype(
                        subtype_label='got',
                        subtype=typ,
                        supertype_label='expected',
                        supertype=return_type,
                        context=s.expr,
                        outer_context=s,
                        msg=message_registry.INCOMPATIBLE_RETURN_VALUE_TYPE,
                        code=codes.RETURN_VALUE)
            else:
                # Empty returns are valid in Generators with Any typed returns, but not in
                # coroutines.
                if (defn.is_generator and not defn.is_coroutine and
                        isinstance(return_type, AnyType)):
                    return

                if isinstance(return_type, (NoneType, AnyType)):
                    return

                if self.in_checked_function():
                    self.fail(message_registry.RETURN_VALUE_EXPECTED, s)

    def visit_if_stmt(self, s: IfStmt) -> None:
        """Type check an if statement."""
        # This frame records the knowledge from previous if/elif clauses not being taken.
        # Fall-through to the original frame is handled explicitly in each block.
        with self.binder.frame_context(can_skip=False, conditional_frame=True, fall_through=0):
            for e, b in zip(s.expr, s.body):
                t = get_proper_type(self.expr_checker.accept(e))

                if isinstance(t, DeletedType):
                    self.msg.deleted_as_rvalue(t, s)

                if_map, else_map = self.find_isinstance_check(e)

                # XXX Issue a warning if condition is always False?
                with self.binder.frame_context(can_skip=True, fall_through=2):
                    self.push_type_map(if_map)
                    self.accept(b)

                # XXX Issue a warning if condition is always True?
                self.push_type_map(else_map)

            with self.binder.frame_context(can_skip=False, fall_through=2):
                if s.else_body:
                    self.accept(s.else_body)

    def visit_while_stmt(self, s: WhileStmt) -> None:
        """Type check a while statement."""
        if_stmt = IfStmt([s.expr], [s.body], None)
        if_stmt.set_line(s.get_line(), s.get_column())
        self.accept_loop(if_stmt, s.else_body,
                         exit_condition=s.expr)

    def visit_operator_assignment_stmt(self,
                                       s: OperatorAssignmentStmt) -> None:
        """Type check an operator assignment statement, e.g. x += 1."""
        self.try_infer_partial_generic_type_from_assignment(s.lvalue, s.rvalue, s.op)
        if isinstance(s.lvalue, MemberExpr):
            # Special case, some additional errors may be given for
            # assignments to read-only or final attributes.
            lvalue_type = self.expr_checker.visit_member_expr(s.lvalue, True)
        else:
            lvalue_type = self.expr_checker.accept(s.lvalue)
        inplace, method = infer_operator_assignment_method(lvalue_type, s.op)
        if inplace:
            # There is __ifoo__, treat as x = x.__ifoo__(y)
            rvalue_type, method_type = self.expr_checker.check_op(
                method, lvalue_type, s.rvalue, s)
            if not is_subtype(rvalue_type, lvalue_type):
                self.msg.incompatible_operator_assignment(s.op, s)
        else:
            # There is no __ifoo__, treat as x = x <foo> y
            expr = OpExpr(s.op, s.lvalue, s.rvalue)
            expr.set_line(s)
            self.check_assignment(lvalue=s.lvalue, rvalue=expr,
                                  infer_lvalue_type=True, new_syntax=False)
        self.check_final(s)

    def visit_assert_stmt(self, s: AssertStmt) -> None:
        self.expr_checker.accept(s.expr)

        if isinstance(s.expr, TupleExpr) and len(s.expr.items) > 0:
            self.fail(message_registry.MALFORMED_ASSERT, s)

        # If this is asserting some isinstance check, bind that type in the following code
        true_map, else_map = self.find_isinstance_check(s.expr)
        if s.msg is not None:
            self.expr_checker.analyze_cond_branch(else_map, s.msg, None)
        self.push_type_map(true_map)

    def visit_raise_stmt(self, s: RaiseStmt) -> None:
        """Type check a raise statement."""
        if s.expr:
            self.type_check_raise(s.expr, s)
        if s.from_expr:
            self.type_check_raise(s.from_expr, s, True)
        self.binder.unreachable()

    def type_check_raise(self, e: Expression, s: RaiseStmt,
                         optional: bool = False) -> None:
        typ = get_proper_type(self.expr_checker.accept(e))
        if isinstance(typ, DeletedType):
            self.msg.deleted_as_rvalue(typ, e)
            return
        exc_type = self.named_type('builtins.BaseException')
        expected_type = UnionType([exc_type, TypeType(exc_type)])
        if optional:
            expected_type.items.append(NoneType())
        if self.options.python_version[0] == 2:
            # allow `raise type, value, traceback`
            # https://docs.python.org/2/reference/simple_stmts.html#the-raise-statement
            # TODO: Also check tuple item types.
            any_type = AnyType(TypeOfAny.implementation_artifact)
            tuple_type = self.named_type('builtins.tuple')
            expected_type.items.append(TupleType([any_type, any_type], tuple_type))
            expected_type.items.append(TupleType([any_type, any_type, any_type], tuple_type))
        self.check_subtype(typ, expected_type, s, message_registry.INVALID_EXCEPTION)

        if isinstance(typ, FunctionLike):
            # https://github.com/python/mypy/issues/11089
            self.expr_checker.check_call(typ, [], [], e)

    def visit_try_stmt(self, s: TryStmt) -> None:
        """Type check a try statement."""
        # Our enclosing frame will get the result if the try/except falls through.
        # This one gets all possible states after the try block exited abnormally
        # (by exception, return, break, etc.)
        with self.binder.frame_context(can_skip=False, fall_through=0):
            # Not only might the body of the try statement exit
            # abnormally, but so might an exception handler or else
            # clause. The finally clause runs in *all* cases, so we
            # need an outer try frame to catch all intermediate states
            # in case an exception is raised during an except or else
            # clause. As an optimization, only create the outer try
            # frame when there actually is a finally clause.
            self.visit_try_without_finally(s, try_frame=bool(s.finally_body))
            if s.finally_body:
                # First we check finally_body is type safe on all abnormal exit paths
                self.accept(s.finally_body)

        if s.finally_body:
            # Then we try again for the more restricted set of options
            # that can fall through. (Why do we need to check the
            # finally clause twice? Depending on whether the finally
            # clause was reached by the try clause falling off the end
            # or exiting abnormally, after completing the finally clause
            # either flow will continue to after the entire try statement
            # or the exception/return/etc. will be processed and control
            # flow will escape. We need to check that the finally clause
            # type checks in both contexts, but only the resulting types
            # from the latter context affect the type state in the code
            # that follows the try statement.)
            if not self.binder.is_unreachable():
                self.accept(s.finally_body)

    def visit_try_without_finally(self, s: TryStmt, try_frame: bool) -> None:
        """Type check a try statement, ignoring the finally block.

        On entry, the top frame should receive all flow that exits the
        try block abnormally (i.e., such that the else block does not
        execute), and its parent should receive all flow that exits
        the try block normally.
        """
        # This frame will run the else block if the try fell through.
        # In that case, control flow continues to the parent of what
        # was the top frame on entry.
        with self.binder.frame_context(can_skip=False, fall_through=2, try_frame=try_frame):
            # This frame receives exit via exception, and runs exception handlers
            with self.binder.frame_context(can_skip=False, conditional_frame=True, fall_through=2):
                # Finally, the body of the try statement
                with self.binder.frame_context(can_skip=False, fall_through=2, try_frame=True):
                    self.accept(s.body)
                for i in range(len(s.handlers)):
                    with self.binder.frame_context(can_skip=True, fall_through=4):
                        typ = s.types[i]
                        if typ:
                            t = self.check_except_handler_test(typ)
                            var = s.vars[i]
                            if var:
                                # To support local variables, we make this a definition line,
                                # causing assignment to set the variable's type.
                                var.is_inferred_def = True
                                # We also temporarily set current_node_deferred to False to
                                # make sure the inference happens.
                                # TODO: Use a better solution, e.g. a
                                # separate Var for each except block.
                                am_deferring = self.current_node_deferred
                                self.current_node_deferred = False
                                self.check_assignment(var, self.temp_node(t, var))
                                self.current_node_deferred = am_deferring
                        self.accept(s.handlers[i])
                        var = s.vars[i]
                        if var:
                            # Exception variables are deleted in python 3 but not python 2.
                            # But, since it's bad form in python 2 and the type checking
                            # wouldn't work very well, we delete it anyway.

                            # Unfortunately, this doesn't let us detect usage before the
                            # try/except block.
                            if self.options.python_version[0] >= 3:
                                source = var.name
                            else:
                                source = ('(exception variable "{}", which we do not '
                                          'accept outside except: blocks even in '
                                          'python 2)'.format(var.name))
                            if isinstance(var.node, Var):
                                var.node.type = DeletedType(source=source)
                            self.binder.cleanse(var)
            if s.else_body:
                self.accept(s.else_body)

    def check_except_handler_test(self, n: Expression) -> Type:
        """Type check an exception handler test clause."""
        typ = self.expr_checker.accept(n)

        all_types: List[Type] = []
        test_types = self.get_types_from_except_handler(typ, n)

        for ttype in get_proper_types(test_types):
            if isinstance(ttype, AnyType):
                all_types.append(ttype)
                continue

            if isinstance(ttype, FunctionLike):
                item = ttype.items[0]
                if not item.is_type_obj():
                    self.fail(message_registry.INVALID_EXCEPTION_TYPE, n)
                    return AnyType(TypeOfAny.from_error)
                exc_type = item.ret_type
            elif isinstance(ttype, TypeType):
                exc_type = ttype.item
            else:
                self.fail(message_registry.INVALID_EXCEPTION_TYPE, n)
                return AnyType(TypeOfAny.from_error)

            if not is_subtype(exc_type, self.named_type('builtins.BaseException')):
                self.fail(message_registry.INVALID_EXCEPTION_TYPE, n)
                return AnyType(TypeOfAny.from_error)

            all_types.append(exc_type)

        return make_simplified_union(all_types)

    def get_types_from_except_handler(self, typ: Type, n: Expression) -> List[Type]:
        """Helper for check_except_handler_test to retrieve handler types."""
        typ = get_proper_type(typ)
        if isinstance(typ, TupleType):
            return typ.items
        elif isinstance(typ, UnionType):
            return [
                union_typ
                for item in typ.relevant_items()
                for union_typ in self.get_types_from_except_handler(item, n)
            ]
        elif isinstance(typ, Instance) and is_named_instance(typ, 'builtins.tuple'):
            # variadic tuple
            return [typ.args[0]]
        else:
            return [typ]

    def visit_for_stmt(self, s: ForStmt) -> None:
        """Type check a for statement."""
        if s.is_async:
            iterator_type, item_type = self.analyze_async_iterable_item_type(s.expr)
        else:
            iterator_type, item_type = self.analyze_iterable_item_type(s.expr)
        s.inferred_item_type = item_type
        s.inferred_iterator_type = iterator_type
        self.analyze_index_variables(s.index, item_type, s.index_type is None, s)
        self.accept_loop(s.body, s.else_body)

    def analyze_async_iterable_item_type(self, expr: Expression) -> Tuple[Type, Type]:
        """Analyse async iterable expression and return iterator and iterator item types."""
        echk = self.expr_checker
        iterable = echk.accept(expr)
        iterator = echk.check_method_call_by_name('__aiter__', iterable, [], [], expr)[0]
        awaitable = echk.check_method_call_by_name('__anext__', iterator, [], [], expr)[0]
        item_type = echk.check_awaitable_expr(awaitable, expr,
                                              message_registry.INCOMPATIBLE_TYPES_IN_ASYNC_FOR)
        return iterator, item_type

    def analyze_iterable_item_type(self, expr: Expression) -> Tuple[Type, Type]:
        """Analyse iterable expression and return iterator and iterator item types."""
        echk = self.expr_checker
        iterable = get_proper_type(echk.accept(expr))
        iterator = echk.check_method_call_by_name('__iter__', iterable, [], [], expr)[0]

        if isinstance(iterable, TupleType):
            joined: Type = UninhabitedType()
            for item in iterable.items:
                joined = join_types(joined, item)
            return iterator, joined
        else:
            # Non-tuple iterable.
            if self.options.python_version[0] >= 3:
                nextmethod = '__next__'
            else:
                nextmethod = 'next'
            return iterator, echk.check_method_call_by_name(nextmethod, iterator, [], [], expr)[0]

    def analyze_container_item_type(self, typ: Type) -> Optional[Type]:
        """Check if a type is a nominal container of a union of such.

        Return the corresponding container item type.
        """
        typ = get_proper_type(typ)
        if isinstance(typ, UnionType):
            types: List[Type] = []
            for item in typ.items:
                c_type = self.analyze_container_item_type(item)
                if c_type:
                    types.append(c_type)
            return UnionType.make_union(types)
        if isinstance(typ, Instance) and typ.type.has_base('typing.Container'):
            supertype = self.named_type('typing.Container').type
            super_instance = map_instance_to_supertype(typ, supertype)
            assert len(super_instance.args) == 1
            return super_instance.args[0]
        if isinstance(typ, TupleType):
            return self.analyze_container_item_type(tuple_fallback(typ))
        return None

    def analyze_index_variables(self, index: Expression, item_type: Type,
                                infer_lvalue_type: bool, context: Context) -> None:
        """Type check or infer for loop or list comprehension index vars."""
        self.check_assignment(index, self.temp_node(item_type, context), infer_lvalue_type)

    def visit_del_stmt(self, s: DelStmt) -> None:
        if isinstance(s.expr, IndexExpr):
            e = s.expr
            m = MemberExpr(e.base, '__delitem__')
            m.line = s.line
            m.column = s.column
            c = CallExpr(m, [e.index], [nodes.ARG_POS], [None])
            c.line = s.line
            c.column = s.column
            self.expr_checker.accept(c, allow_none_return=True)
        else:
            s.expr.accept(self.expr_checker)
            for elt in flatten(s.expr):
                if isinstance(elt, NameExpr):
                    self.binder.assign_type(elt, DeletedType(source=elt.name),
                                            get_declaration(elt), False)

    def visit_decorator(self, e: Decorator) -> None:
        for d in e.decorators:
            if isinstance(d, RefExpr):
                if d.fullname == 'typing.no_type_check':
                    e.var.type = AnyType(TypeOfAny.special_form)
                    e.var.is_ready = True
                    return

        if self.recurse_into_functions:
            with self.tscope.function_scope(e.func):
                self.check_func_item(e.func, name=e.func.name)

        # Process decorators from the inside out to determine decorated signature, which
        # may be different from the declared signature.
        sig: Type = self.function_type(e.func)
        for d in reversed(e.decorators):
            if refers_to_fullname(d, 'typing.overload'):
                self.fail(message_registry.MULTIPLE_OVERLOADS_REQUIRED, e)
                continue
            dec = self.expr_checker.accept(d)
            temp = self.temp_node(sig, context=e)
            fullname = None
            if isinstance(d, RefExpr):
                fullname = d.fullname
            # if this is a expression like @b.a where b is an object, get the type of b
            # so we can pass it the method hook in the plugins
            object_type: Optional[Type] = None
            if fullname is None and isinstance(d, MemberExpr) and d.expr in self.type_map:
                object_type = self.type_map[d.expr]
                fullname = self.expr_checker.method_fullname(object_type, d.name)
            self.check_for_untyped_decorator(e.func, dec, d)
            sig, t2 = self.expr_checker.check_call(dec, [temp],
                                                   [nodes.ARG_POS], e,
                                                   callable_name=fullname,
                                                   object_type=object_type)
        self.check_untyped_after_decorator(sig, e.func)
        sig = set_callable_name(sig, e.func)
        e.var.type = sig
        e.var.is_ready = True
        if e.func.is_property:
            self.check_incompatible_property_override(e)
        if e.func.info and not e.func.is_dynamic():
            self.check_method_override(e)

        if e.func.info and e.func.name in ('__init__', '__new__'):
            if e.type and not isinstance(get_proper_type(e.type), (FunctionLike, AnyType)):
                self.fail(message_registry.BAD_CONSTRUCTOR_TYPE, e)

    def check_for_untyped_decorator(self,
                                    func: FuncDef,
                                    dec_type: Type,
                                    dec_expr: Expression) -> None:
        if (self.options.disallow_untyped_decorators and
                is_typed_callable(func.type) and
                is_untyped_decorator(dec_type)):
            self.msg.typed_function_untyped_decorator(func.name, dec_expr)

    def check_incompatible_property_override(self, e: Decorator) -> None:
        if not e.var.is_settable_property and e.func.info:
            name = e.func.name
            for base in e.func.info.mro[1:]:
                base_attr = base.names.get(name)
                if not base_attr:
                    continue
                if (isinstance(base_attr.node, OverloadedFuncDef) and
                        base_attr.node.is_property and
                        cast(Decorator,
                             base_attr.node.items[0]).var.is_settable_property):
                    self.fail(message_registry.READ_ONLY_PROPERTY_OVERRIDES_READ_WRITE, e)

    def visit_with_stmt(self, s: WithStmt) -> None:
        exceptions_maybe_suppressed = False
        for expr, target in zip(s.expr, s.target):
            if s.is_async:
                exit_ret_type = self.check_async_with_item(expr, target, s.unanalyzed_type is None)
            else:
                exit_ret_type = self.check_with_item(expr, target, s.unanalyzed_type is None)

            # Based on the return type, determine if this context manager 'swallows'
            # exceptions or not. We determine this using a heuristic based on the
            # return type of the __exit__ method -- see the discussion in
            # https://github.com/python/mypy/issues/7214 and the section about context managers
            # in https://github.com/python/typeshed/blob/master/CONTRIBUTING.md#conventions
            # for more details.

            exit_ret_type = get_proper_type(exit_ret_type)
            if is_literal_type(exit_ret_type, "builtins.bool", False):
                continue

            if (is_literal_type(exit_ret_type, "builtins.bool", True)
                    or (isinstance(exit_ret_type, Instance)
                        and exit_ret_type.type.fullname == 'builtins.bool'
                        and state.strict_optional)):
                # Note: if strict-optional is disabled, this bool instance
                # could actually be an Optional[bool].
                exceptions_maybe_suppressed = True

        if exceptions_maybe_suppressed:
            # Treat this 'with' block in the same way we'd treat a 'try: BODY; except: pass'
            # block. This means control flow can continue after the 'with' even if the 'with'
            # block immediately returns.
            with self.binder.frame_context(can_skip=True, try_frame=True):
                self.accept(s.body)
        else:
            self.accept(s.body)

    def check_untyped_after_decorator(self, typ: Type, func: FuncDef) -> None:
        if not self.options.disallow_any_decorated or self.is_stub:
            return

        if mypy.checkexpr.has_any_type(typ):
            self.msg.untyped_decorated_function(typ, func)

    def check_async_with_item(self, expr: Expression, target: Optional[Expression],
                              infer_lvalue_type: bool) -> Type:
        echk = self.expr_checker
        ctx = echk.accept(expr)
        obj = echk.check_method_call_by_name('__aenter__', ctx, [], [], expr)[0]
        obj = echk.check_awaitable_expr(
            obj, expr, message_registry.INCOMPATIBLE_TYPES_IN_ASYNC_WITH_AENTER)
        if target:
            self.check_assignment(target, self.temp_node(obj, expr), infer_lvalue_type)
        arg = self.temp_node(AnyType(TypeOfAny.special_form), expr)
        res, _ = echk.check_method_call_by_name(
            '__aexit__', ctx, [arg] * 3, [nodes.ARG_POS] * 3, expr)
        return echk.check_awaitable_expr(
            res, expr, message_registry.INCOMPATIBLE_TYPES_IN_ASYNC_WITH_AEXIT)

    def check_with_item(self, expr: Expression, target: Optional[Expression],
                        infer_lvalue_type: bool) -> Type:
        echk = self.expr_checker
        ctx = echk.accept(expr)
        obj = echk.check_method_call_by_name('__enter__', ctx, [], [], expr)[0]
        if target:
            self.check_assignment(target, self.temp_node(obj, expr), infer_lvalue_type)
        arg = self.temp_node(AnyType(TypeOfAny.special_form), expr)
        res, _ = echk.check_method_call_by_name(
            '__exit__', ctx, [arg] * 3, [nodes.ARG_POS] * 3, expr)
        return res

    def visit_print_stmt(self, s: PrintStmt) -> None:
        for arg in s.args:
            self.expr_checker.accept(arg)
        if s.target:
            target_type = get_proper_type(self.expr_checker.accept(s.target))
            if not isinstance(target_type, NoneType):
                write_type = self.expr_checker.analyze_external_member_access(
                    'write', target_type, s.target)
                required_type = CallableType(
                    arg_types=[self.named_type('builtins.str')],
                    arg_kinds=[ARG_POS],
                    arg_names=[None],
                    ret_type=AnyType(TypeOfAny.implementation_artifact),
                    fallback=self.named_type('builtins.function'),
                )
                # This has to be hard-coded, since it is a syntax pattern, not a function call.
                if not is_subtype(write_type, required_type):
                    self.fail(message_registry.PYTHON2_PRINT_FILE_TYPE.format(
                        write_type,
                        required_type,
                    ), s.target)

    def visit_break_stmt(self, s: BreakStmt) -> None:
        self.binder.handle_break()

    def visit_continue_stmt(self, s: ContinueStmt) -> None:
        self.binder.handle_continue()
        return None

    def make_fake_typeinfo(self,
                           curr_module_fullname: str,
                           class_gen_name: str,
                           class_short_name: str,
                           bases: List[Instance],
                           ) -> Tuple[ClassDef, TypeInfo]:
        # Build the fake ClassDef and TypeInfo together.
        # The ClassDef is full of lies and doesn't actually contain a body.
        # Use format_bare to generate a nice name for error messages.
        # We skip fully filling out a handful of TypeInfo fields because they
        # should be irrelevant for a generated type like this:
        # is_protocol, protocol_members, is_abstract
        cdef = ClassDef(class_short_name, Block([]))
        cdef.fullname = curr_module_fullname + '.' + class_gen_name
        info = TypeInfo(SymbolTable(), cdef, curr_module_fullname)
        cdef.info = info
        info.bases = bases
        calculate_mro(info)
        info.calculate_metaclass_type()
        return cdef, info

    def intersect_instances(self,
                            instances: Tuple[Instance, Instance],
                            ctx: Context,
                            ) -> Optional[Instance]:
        """Try creating an ad-hoc intersection of the given instances.

        Note that this function does *not* try and create a full-fledged
        intersection type. Instead, it returns an instance of a new ad-hoc
        subclass of the given instances.

        This is mainly useful when you need a way of representing some
        theoretical subclass of the instances the user may be trying to use
        the generated intersection can serve as a placeholder.

        This function will create a fresh subclass every time you call it,
        even if you pass in the exact same arguments. So this means calling
        `self.intersect_intersection([inst_1, inst_2], ctx)` twice will result
        in instances of two distinct subclasses of inst_1 and inst_2.

        This is by design: we want each ad-hoc intersection to be unique since
        they're supposed represent some other unknown subclass.

        Returns None if creating the subclass is impossible (e.g. due to
        MRO errors or incompatible signatures). If we do successfully create
        a subclass, its TypeInfo will automatically be added to the global scope.
        """
        curr_module = self.scope.stack[0]
        assert isinstance(curr_module, MypyFile)

        def _get_base_classes(instances_: Tuple[Instance, Instance]) -> List[Instance]:
            base_classes_ = []
            for inst in instances_:
                if inst.type.is_intersection:
                    expanded = inst.type.bases
                else:
                    expanded = [inst]

                for expanded_inst in expanded:
                    base_classes_.append(expanded_inst)
            return base_classes_

        def _make_fake_typeinfo_and_full_name(
                base_classes_: List[Instance],
                curr_module_: MypyFile,
        ) -> Tuple[TypeInfo, str]:
            names_list = pretty_seq([x.type.name for x in base_classes_], "and")
            short_name = '<subclass of {}>'.format(names_list)
            full_name_ = gen_unique_name(short_name, curr_module_.names)
            cdef, info_ = self.make_fake_typeinfo(
                curr_module_.fullname,
                full_name_,
                short_name,
                base_classes_,
            )
            return info_, full_name_

        old_msg = self.msg
        new_msg = old_msg.clean_copy()
        self.msg = new_msg
        base_classes = _get_base_classes(instances)
        # We use the pretty_names_list for error messages but can't
        # use it for the real name that goes into the symbol table
        # because it can have dots in it.
        pretty_names_list = pretty_seq(format_type_distinctly(*base_classes, bare=True), "and")
        try:
            info, full_name = _make_fake_typeinfo_and_full_name(base_classes, curr_module)
            self.check_multiple_inheritance(info)
            if new_msg.is_errors():
                # "class A(B, C)" unsafe, now check "class A(C, B)":
                new_msg = new_msg.clean_copy()
                self.msg = new_msg
                base_classes = _get_base_classes(instances[::-1])
                info, full_name = _make_fake_typeinfo_and_full_name(base_classes, curr_module)
                self.check_multiple_inheritance(info)
            info.is_intersection = True
        except MroError:
            if self.should_report_unreachable_issues():
                old_msg.impossible_intersection(
                    pretty_names_list, "inconsistent method resolution order", ctx)
            return None
        finally:
            self.msg = old_msg

        if new_msg.is_errors():
            if self.should_report_unreachable_issues():
                self.msg.impossible_intersection(
                    pretty_names_list, "incompatible method signatures", ctx)
            return None

        curr_module.names[full_name] = SymbolTableNode(GDEF, info)
        return Instance(info, [])

    def intersect_instance_callable(self, typ: Instance, callable_type: CallableType) -> Instance:
        """Creates a fake type that represents the intersection of an Instance and a CallableType.

        It operates by creating a bare-minimum dummy TypeInfo that
        subclasses type and adds a __call__ method matching callable_type.
        """

        # In order for this to work in incremental mode, the type we generate needs to
        # have a valid fullname and a corresponding entry in a symbol table. We generate
        # a unique name inside the symbol table of the current module.
        cur_module = cast(MypyFile, self.scope.stack[0])
        gen_name = gen_unique_name("<callable subtype of {}>".format(typ.type.name),
                                   cur_module.names)

        # Synthesize a fake TypeInfo
        short_name = format_type_bare(typ)
        cdef, info = self.make_fake_typeinfo(cur_module.fullname, gen_name, short_name, [typ])

        # Build up a fake FuncDef so we can populate the symbol table.
        func_def = FuncDef('__call__', [], Block([]), callable_type)
        func_def._fullname = cdef.fullname + '.__call__'
        func_def.info = info
        info.names['__call__'] = SymbolTableNode(MDEF, func_def)

        cur_module.names[gen_name] = SymbolTableNode(GDEF, info)

        return Instance(info, [])

    def make_fake_callable(self, typ: Instance) -> Instance:
        """Produce a new type that makes type Callable with a generic callable type."""

        fallback = self.named_type('builtins.function')
        callable_type = CallableType([AnyType(TypeOfAny.explicit),
                                      AnyType(TypeOfAny.explicit)],
                                     [nodes.ARG_STAR, nodes.ARG_STAR2],
                                     [None, None],
                                     ret_type=AnyType(TypeOfAny.explicit),
                                     fallback=fallback,
                                     is_ellipsis_args=True)

        return self.intersect_instance_callable(typ, callable_type)

    def partition_by_callable(self, typ: Type,
                              unsound_partition: bool) -> Tuple[List[Type], List[Type]]:
        """Partitions a type into callable subtypes and uncallable subtypes.

        Thus, given:
        `callables, uncallables = partition_by_callable(type)`

        If we assert `callable(type)` then `type` has type Union[*callables], and
        If we assert `not callable(type)` then `type` has type Union[*uncallables]

        If unsound_partition is set, assume that anything that is not
        clearly callable is in fact not callable. Otherwise we generate a
        new subtype that *is* callable.

        Guaranteed to not return [], [].
        """
        typ = get_proper_type(typ)

        if isinstance(typ, FunctionLike) or isinstance(typ, TypeType):
            return [typ], []

        if isinstance(typ, AnyType):
            return [typ], [typ]

        if isinstance(typ, NoneType):
            return [], [typ]

        if isinstance(typ, UnionType):
            callables = []
            uncallables = []
            for subtype in typ.items:
                # Use unsound_partition when handling unions in order to
                # allow the expected type discrimination.
                subcallables, subuncallables = self.partition_by_callable(subtype,
                                                                          unsound_partition=True)
                callables.extend(subcallables)
                uncallables.extend(subuncallables)
            return callables, uncallables

        if isinstance(typ, TypeVarType):
            # We could do better probably?
            # Refine the the type variable's bound as our type in the case that
            # callable() is true. This unfortunately loses the information that
            # the type is a type variable in that branch.
            # This matches what is done for isinstance, but it may be possible to
            # do better.
            # If it is possible for the false branch to execute, return the original
            # type to avoid losing type information.
            callables, uncallables = self.partition_by_callable(erase_to_union_or_bound(typ),
                                                                unsound_partition)
            uncallables = [typ] if len(uncallables) else []
            return callables, uncallables

        # A TupleType is callable if its fallback is, but needs special handling
        # when we dummy up a new type.
        ityp = typ
        if isinstance(typ, TupleType):
            ityp = tuple_fallback(typ)

        if isinstance(ityp, Instance):
            method = ityp.type.get_method('__call__')
            if method and method.type:
                callables, uncallables = self.partition_by_callable(method.type,
                                                                    unsound_partition=False)
                if len(callables) and not len(uncallables):
                    # Only consider the type callable if its __call__ method is
                    # definitely callable.
                    return [typ], []

            if not unsound_partition:
                fake = self.make_fake_callable(ityp)
                if isinstance(typ, TupleType):
                    fake.type.tuple_type = TupleType(typ.items, fake)
                    return [fake.type.tuple_type], [typ]
                return [fake], [typ]

        if unsound_partition:
            return [], [typ]
        else:
            # We don't know how properly make the type callable.
            return [typ], [typ]

    def conditional_callable_type_map(self, expr: Expression,
                                      current_type: Optional[Type],
                                      ) -> Tuple[TypeMap, TypeMap]:
        """Takes in an expression and the current type of the expression.

        Returns a 2-tuple: The first element is a map from the expression to
        the restricted type if it were callable. The second element is a
        map from the expression to the type it would hold if it weren't
        callable.
        """
        if not current_type:
            return {}, {}

        if isinstance(get_proper_type(current_type), AnyType):
            return {}, {}

        callables, uncallables = self.partition_by_callable(current_type,
                                                            unsound_partition=False)

        if len(callables) and len(uncallables):
            callable_map = {expr: UnionType.make_union(callables)} if len(callables) else None
            uncallable_map = {
                expr: UnionType.make_union(uncallables)} if len(uncallables) else None
            return callable_map, uncallable_map

        elif len(callables):
            return {}, None

        return None, {}

    def _is_truthy_type(self, t: ProperType) -> bool:
        return (
            (
                isinstance(t, Instance) and
                bool(t.type) and
                not t.type.has_readable_member('__bool__') and
                not t.type.has_readable_member('__len__')
            )
            or isinstance(t, FunctionLike)
            or (
                isinstance(t, UnionType) and
                all(self._is_truthy_type(t) for t in get_proper_types(t.items))
            )
        )

    def _check_for_truthy_type(self, t: Type, expr: Expression) -> None:
        if not state.strict_optional:
            return  # if everything can be None, all bets are off

        t = get_proper_type(t)
        if not self._is_truthy_type(t):
            return

        def format_expr_type() -> str:
            typ = format_type(t)
            if isinstance(expr, MemberExpr):
                return f'Member "{expr.name}" has type {typ}'
            elif isinstance(expr, RefExpr) and expr.fullname:
                return f'"{expr.fullname}" has type {typ}'
            elif isinstance(expr, CallExpr):
                if isinstance(expr.callee, MemberExpr):
                    return f'"{expr.callee.name}" returns {typ}'
                elif isinstance(expr.callee, RefExpr) and expr.callee.fullname:
                    return f'"{expr.callee.fullname}" returns {typ}'
                return f'Call returns {typ}'
            else:
                return f'Expression has type {typ}'

        if isinstance(t, FunctionLike):
            self.fail(message_registry.FUNCTION_ALWAYS_TRUE.format(format_type(t)), expr)
        elif isinstance(t, UnionType):
            self.fail(message_registry.TYPE_ALWAYS_TRUE_UNIONTYPE.format(format_expr_type()),
                      expr)
        else:
            self.fail(message_registry.TYPE_ALWAYS_TRUE.format(format_expr_type()), expr)

    def find_type_equals_check(self, node: ComparisonExpr, expr_indices: List[int]
                               ) -> Tuple[TypeMap, TypeMap]:
        """Narrow types based on any checks of the type ``type(x) == T``

        Args:
            node: The node that might contain the comparison
            expr_indices: The list of indices of expressions in ``node`` that are being
                compared
        """
        type_map = self.type_map

        def is_type_call(expr: CallExpr) -> bool:
            """Is expr a call to type with one argument?"""
            return (refers_to_fullname(expr.callee, 'builtins.type')
                    and len(expr.args) == 1)

        # exprs that are being passed into type
        exprs_in_type_calls: List[Expression] = []
        # type that is being compared to type(expr)
        type_being_compared: Optional[List[TypeRange]] = None
        # whether the type being compared to is final
        is_final = False

        for index in expr_indices:
            expr = node.operands[index]

            if isinstance(expr, CallExpr) and is_type_call(expr):
                exprs_in_type_calls.append(expr.args[0])
            else:
                current_type = get_isinstance_type(expr, type_map)
                if current_type is None:
                    continue
                if type_being_compared is not None:
                    # It doesn't really make sense to have several types being
                    # compared to the output of type (like type(x) == int == str)
                    # because whether that's true is solely dependent on what the
                    # types being compared are, so we don't try to narrow types any
                    # further because we can't really get any information about the
                    # type of x from that check
                    return {}, {}
                else:
                    if isinstance(expr, RefExpr) and isinstance(expr.node, TypeInfo):
                        is_final = expr.node.is_final
                    type_being_compared = current_type

        if not exprs_in_type_calls:
            return {}, {}

        if_maps: List[TypeMap] = []
        else_maps: List[TypeMap] = []
        for expr in exprs_in_type_calls:
            current_if_map, current_else_map = self.conditional_type_map_with_intersection(
                expr,
                type_map[expr],
                type_being_compared
            )
            if_maps.append(current_if_map)
            else_maps.append(current_else_map)

        def combine_maps(list_maps: List[TypeMap]) -> TypeMap:
            """Combine all typemaps in list_maps into one typemap"""
            result_map = {}
            for d in list_maps:
                if d is not None:
                    result_map.update(d)
            return result_map

        if_map = combine_maps(if_maps)
        # type(x) == T is only true when x has the same type as T, meaning
        # that it can be false if x is an instance of a subclass of T. That means
        # we can't do any narrowing in the else case unless T is final, in which
        # case T can't be subclassed
        if is_final:
            else_map = combine_maps(else_maps)
        else:
            else_map = {}
        return if_map, else_map

    def find_isinstance_check(self, node: Expression
                              ) -> Tuple[TypeMap, TypeMap]:
        """Find any isinstance checks (within a chain of ands).  Includes
        implicit and explicit checks for None and calls to callable.
        Also includes TypeGuard functions.

        Return value is a map of variables to their types if the condition
        is true and a map of variables to their types if the condition is false.

        If either of the values in the tuple is None, then that particular
        branch can never occur.

        May return {}, {}.
        Can return None, None in situations involving NoReturn.
        """
        if_map, else_map = self.find_isinstance_check_helper(node)
        new_if_map = self.propagate_up_typemap_info(self.type_map, if_map)
        new_else_map = self.propagate_up_typemap_info(self.type_map, else_map)
        return new_if_map, new_else_map

    def find_isinstance_check_helper(self, node: Expression) -> Tuple[TypeMap, TypeMap]:
        type_map = self.type_map
        if is_true_literal(node):
            return {}, None
        if is_false_literal(node):
            return None, {}

        if isinstance(node, CallExpr) and len(node.args) != 0:
            expr = collapse_walrus(node.args[0])
            if refers_to_fullname(node.callee, 'builtins.isinstance'):
                if len(node.args) != 2:  # the error will be reported elsewhere
                    return {}, {}
                if literal(expr) == LITERAL_TYPE:
                    return self.conditional_type_map_with_intersection(
                        expr,
                        type_map[expr],
                        get_isinstance_type(node.args[1], type_map),
                    )
            elif refers_to_fullname(node.callee, 'builtins.issubclass'):
                if len(node.args) != 2:  # the error will be reported elsewhere
                    return {}, {}
                if literal(expr) == LITERAL_TYPE:
                    return self.infer_issubclass_maps(node, expr, type_map)
            elif refers_to_fullname(node.callee, 'builtins.callable'):
                if len(node.args) != 1:  # the error will be reported elsewhere
                    return {}, {}
                if literal(expr) == LITERAL_TYPE:
                    vartype = type_map[expr]
                    return self.conditional_callable_type_map(expr, vartype)
            elif isinstance(node.callee, RefExpr):
                if node.callee.type_guard is not None:
                    # TODO: Follow keyword args or *args, **kwargs
                    if node.arg_kinds[0] != nodes.ARG_POS:
                        self.fail(message_registry.TYPE_GUARD_POS_ARG_REQUIRED, node)
                        return {}, {}
                    if literal(expr) == LITERAL_TYPE:
                        # Note: we wrap the target type, so that we can special case later.
                        # Namely, for isinstance() we use a normal meet, while TypeGuard is
                        # considered "always right" (i.e. even if the types are not overlapping).
                        # Also note that a care must be taken to unwrap this back at read places
                        # where we use this to narrow down declared type.
                        return {expr: TypeGuardedType(node.callee.type_guard)}, {}
        elif isinstance(node, ComparisonExpr):
            # Step 1: Obtain the types of each operand and whether or not we can
            # narrow their types. (For example, we shouldn't try narrowing the
            # types of literal string or enum expressions).

            operands = [collapse_walrus(x) for x in node.operands]
            operand_types = []
            narrowable_operand_index_to_hash = {}
            for i, expr in enumerate(operands):
                if expr not in type_map:
                    return {}, {}
                expr_type = type_map[expr]
                operand_types.append(expr_type)

                if (literal(expr) == LITERAL_TYPE
                        and not is_literal_none(expr)
                        and not is_literal_enum(type_map, expr)):
                    h = literal_hash(expr)
                    if h is not None:
                        narrowable_operand_index_to_hash[i] = h

            # Step 2: Group operands chained by either the 'is' or '==' operands
            # together. For all other operands, we keep them in groups of size 2.
            # So the expression:
            #
            #   x0 == x1 == x2 < x3 < x4 is x5 is x6 is not x7 is not x8
            #
            # ...is converted into the simplified operator list:
            #
            #  [("==", [0, 1, 2]), ("<", [2, 3]), ("<", [3, 4]),
            #   ("is", [4, 5, 6]), ("is not", [6, 7]), ("is not", [7, 8])]
            #
            # We group identity/equality expressions so we can propagate information
            # we discover about one operand across the entire chain. We don't bother
            # handling 'is not' and '!=' chains in a special way: those are very rare
            # in practice.

            simplified_operator_list = group_comparison_operands(
                node.pairwise(),
                narrowable_operand_index_to_hash,
                {'==', 'is'},
            )

            # Step 3: Analyze each group and infer more precise type maps for each
            # assignable operand, if possible. We combine these type maps together
            # in the final step.

            partial_type_maps = []
            for operator, expr_indices in simplified_operator_list:
                if operator in {'is', 'is not', '==', '!='}:
                    # is_valid_target:
                    #   Controls which types we're allowed to narrow exprs to. Note that
                    #   we cannot use 'is_literal_type_like' in both cases since doing
                    #   'x = 10000 + 1; x is 10001' is not always True in all Python
                    #   implementations.
                    #
                    # coerce_only_in_literal_context:
                    #   If true, coerce types into literal types only if one or more of
                    #   the provided exprs contains an explicit Literal type. This could
                    #   technically be set to any arbitrary value, but it seems being liberal
                    #   with narrowing when using 'is' and conservative when using '==' seems
                    #   to break the least amount of real-world code.
                    #
                    # should_narrow_by_identity:
                    #   Set to 'false' only if the user defines custom __eq__ or __ne__ methods
                    #   that could cause identity-based narrowing to produce invalid results.
                    if operator in {'is', 'is not'}:
                        is_valid_target: Callable[[Type], bool] = is_singleton_type
                        coerce_only_in_literal_context = False
                        should_narrow_by_identity = True
                    else:
                        def is_exactly_literal_type(t: Type) -> bool:
                            return isinstance(get_proper_type(t), LiteralType)

                        def has_no_custom_eq_checks(t: Type) -> bool:
                            return (not custom_special_method(t, '__eq__', check_all=False)
                                    and not custom_special_method(t, '__ne__', check_all=False))

                        is_valid_target = is_exactly_literal_type
                        coerce_only_in_literal_context = True

                        expr_types = [operand_types[i] for i in expr_indices]
                        should_narrow_by_identity = all(map(has_no_custom_eq_checks, expr_types))

                    if_map: TypeMap = {}
                    else_map: TypeMap = {}
                    if should_narrow_by_identity:
                        if_map, else_map = self.refine_identity_comparison_expression(
                            operands,
                            operand_types,
                            expr_indices,
                            narrowable_operand_index_to_hash.keys(),
                            is_valid_target,
                            coerce_only_in_literal_context,
                        )

                    # Strictly speaking, we should also skip this check if the objects in the expr
                    # chain have custom __eq__ or __ne__ methods. But we (maybe optimistically)
                    # assume nobody would actually create a custom objects that considers itself
                    # equal to None.
                    if if_map == {} and else_map == {}:
                        if_map, else_map = self.refine_away_none_in_comparison(
                            operands,
                            operand_types,
                            expr_indices,
                            narrowable_operand_index_to_hash.keys(),
                        )

                    # If we haven't been able to narrow types yet, we might be dealing with a
                    # explicit type(x) == some_type check
                    if if_map == {} and else_map == {}:
                        if_map, else_map = self.find_type_equals_check(node, expr_indices)
                elif operator in {'in', 'not in'}:
                    assert len(expr_indices) == 2
                    left_index, right_index = expr_indices
                    if left_index not in narrowable_operand_index_to_hash:
                        continue

                    item_type = operand_types[left_index]
                    collection_type = operand_types[right_index]

                    # We only try and narrow away 'None' for now
                    if not is_optional(item_type):
                        continue

                    collection_item_type = get_proper_type(builtin_item_type(collection_type))
                    if collection_item_type is None or is_optional(collection_item_type):
                        continue
                    if (isinstance(collection_item_type, Instance)
                            and collection_item_type.type.fullname == 'builtins.object'):
                        continue
                    if is_overlapping_erased_types(item_type, collection_item_type):
                        if_map, else_map = {operands[left_index]: remove_optional(item_type)}, {}
                    else:
                        continue
                else:
                    if_map = {}
                    else_map = {}

                if operator in {'is not', '!=', 'not in'}:
                    if_map, else_map = else_map, if_map

                partial_type_maps.append((if_map, else_map))

            return reduce_conditional_maps(partial_type_maps)
        elif isinstance(node, AssignmentExpr):
            if_map = {}
            else_map = {}

            if_assignment_map, else_assignment_map = self.find_isinstance_check(node.target)

            if if_assignment_map is not None:
                if_map.update(if_assignment_map)
            if else_assignment_map is not None:
                else_map.update(else_assignment_map)

            if_condition_map, else_condition_map = self.find_isinstance_check(node.value)

            if if_condition_map is not None:
                if_map.update(if_condition_map)
            if else_condition_map is not None:
                else_map.update(else_condition_map)

            return (
                (None if if_assignment_map is None or if_condition_map is None else if_map),
                (None if else_assignment_map is None or else_condition_map is None else else_map),
            )
        elif isinstance(node, OpExpr) and node.op == 'and':
            left_if_vars, left_else_vars = self.find_isinstance_check(node.left)
            right_if_vars, right_else_vars = self.find_isinstance_check(node.right)

            # (e1 and e2) is true if both e1 and e2 are true,
            # and false if at least one of e1 and e2 is false.
            return (and_conditional_maps(left_if_vars, right_if_vars),
                    or_conditional_maps(left_else_vars, right_else_vars))
        elif isinstance(node, OpExpr) and node.op == 'or':
            left_if_vars, left_else_vars = self.find_isinstance_check(node.left)
            right_if_vars, right_else_vars = self.find_isinstance_check(node.right)

            # (e1 or e2) is true if at least one of e1 or e2 is true,
            # and false if both e1 and e2 are false.
            return (or_conditional_maps(left_if_vars, right_if_vars),
                    and_conditional_maps(left_else_vars, right_else_vars))
        elif isinstance(node, UnaryExpr) and node.op == 'not':
            left, right = self.find_isinstance_check(node.expr)
            return right, left

        # Restrict the type of the variable to True-ish/False-ish in the if and else branches
        # respectively
        original_vartype = type_map[node]
        self._check_for_truthy_type(original_vartype, node)
        vartype = try_expanding_sum_type_to_union(original_vartype, "builtins.bool")

        if_type = true_only(vartype)
        else_type = false_only(vartype)
        if_map = (
            {node: if_type}
            if not isinstance(if_type, UninhabitedType)
            else None
        )
        else_map = (
            {node: else_type}
            if not isinstance(else_type, UninhabitedType)
            else None
        )
        return if_map, else_map

    def propagate_up_typemap_info(self,
                                  existing_types: Mapping[Expression, Type],
                                  new_types: TypeMap) -> TypeMap:
        """Attempts refining parent expressions of any MemberExpr or IndexExprs in new_types.

        Specifically, this function accepts two mappings of expression to original types:
        the original mapping (existing_types), and a new mapping (new_types) intended to
        update the original.

        This function iterates through new_types and attempts to use the information to try
        refining any parent types that happen to be unions.

        For example, suppose there are two types "A = Tuple[int, int]" and "B = Tuple[str, str]".
        Next, suppose that 'new_types' specifies the expression 'foo[0]' has a refined type
        of 'int' and that 'foo' was previously deduced to be of type Union[A, B].

        Then, this function will observe that since A[0] is an int and B[0] is not, the type of
        'foo' can be further refined from Union[A, B] into just B.

        We perform this kind of "parent narrowing" for member lookup expressions and indexing
        expressions into tuples, namedtuples, and typeddicts. We repeat this narrowing
        recursively if the parent is also a "lookup expression". So for example, if we have
        the expression "foo['bar'].baz[0]", we'd potentially end up refining types for the
        expressions "foo", "foo['bar']", and "foo['bar'].baz".

        We return the newly refined map. This map is guaranteed to be a superset of 'new_types'.
        """
        if new_types is None:
            return None
        output_map = {}
        for expr, expr_type in new_types.items():
            # The original inferred type should always be present in the output map, of course
            output_map[expr] = expr_type

            # Next, try using this information to refine the parent types, if applicable.
            new_mapping = self.refine_parent_types(existing_types, expr, expr_type)
            for parent_expr, proposed_parent_type in new_mapping.items():
                # We don't try inferring anything if we've already inferred something for
                # the parent expression.
                # TODO: Consider picking the narrower type instead of always discarding this?
                if parent_expr in new_types:
                    continue
                output_map[parent_expr] = proposed_parent_type
        return output_map

    def refine_parent_types(self,
                            existing_types: Mapping[Expression, Type],
                            expr: Expression,
                            expr_type: Type) -> Mapping[Expression, Type]:
        """Checks if the given expr is a 'lookup operation' into a union and iteratively refines
        the parent types based on the 'expr_type'.

        For example, if 'expr' is an expression like 'a.b.c.d', we'll potentially return refined
        types for expressions 'a', 'a.b', and 'a.b.c'.

        For more details about what a 'lookup operation' is and how we use the expr_type to refine
        the parent types of lookup_expr, see the docstring in 'propagate_up_typemap_info'.
        """
        output: Dict[Expression, Type] = {}

        # Note: parent_expr and parent_type are progressively refined as we crawl up the
        # parent lookup chain.
        while True:
            # First, check if this expression is one that's attempting to
            # "lookup" some key in the parent type. If so, save the parent type
            # and create function that will try replaying the same lookup
            # operation against arbitrary types.
            if isinstance(expr, MemberExpr):
                parent_expr = expr.expr
                parent_type = existing_types.get(parent_expr)
                member_name = expr.name

                def replay_lookup(new_parent_type: ProperType) -> Optional[Type]:
                    msg_copy = self.msg.clean_copy()
                    member_type = analyze_member_access(
                        name=member_name,
                        typ=new_parent_type,
                        context=parent_expr,
                        is_lvalue=False,
                        is_super=False,
                        is_operator=False,
                        msg=msg_copy,
                        original_type=new_parent_type,
                        chk=self,
                        in_literal_context=False,
                    )
                    if msg_copy.is_errors():
                        return None
                    else:
                        return member_type
            elif isinstance(expr, IndexExpr):
                parent_expr = expr.base
                parent_type = existing_types.get(parent_expr)

                index_type = existing_types.get(expr.index)
                if index_type is None:
                    return output

                str_literals = try_getting_str_literals_from_type(index_type)
                if str_literals is not None:
                    # Refactoring these two indexing replay functions is surprisingly
                    # tricky -- see https://github.com/python/mypy/pull/7917, which
                    # was blocked by https://github.com/mypyc/mypyc/issues/586
                    def replay_lookup(new_parent_type: ProperType) -> Optional[Type]:
                        if not isinstance(new_parent_type, TypedDictType):
                            return None
                        try:
                            assert str_literals is not None
                            member_types = [new_parent_type.items[key] for key in str_literals]
                        except KeyError:
                            return None
                        return make_simplified_union(member_types)
                else:
                    int_literals = try_getting_int_literals_from_type(index_type)
                    if int_literals is not None:
                        def replay_lookup(new_parent_type: ProperType) -> Optional[Type]:
                            if not isinstance(new_parent_type, TupleType):
                                return None
                            try:
                                assert int_literals is not None
                                member_types = [new_parent_type.items[key] for key in int_literals]
                            except IndexError:
                                return None
                            return make_simplified_union(member_types)
                    else:
                        return output
            else:
                return output

            # If we somehow didn't previously derive the parent type, abort completely
            # with what we have so far: something went wrong at an earlier stage.
            if parent_type is None:
                return output

            # We currently only try refining the parent type if it's a Union.
            # If not, there's no point in trying to refine any further parents
            # since we have no further information we can use to refine the lookup
            # chain, so we end early as an optimization.
            parent_type = get_proper_type(parent_type)
            if not isinstance(parent_type, UnionType):
                return output

            # Take each element in the parent union and replay the original lookup procedure
            # to figure out which parents are compatible.
            new_parent_types = []
            for item in union_items(parent_type):
                member_type = replay_lookup(item)
                if member_type is None:
                    # We were unable to obtain the member type. So, we give up on refining this
                    # parent type entirely and abort.
                    return output

                if is_overlapping_types(member_type, expr_type):
                    new_parent_types.append(item)

            # If none of the parent types overlap (if we derived an empty union), something
            # went wrong. We should never hit this case, but deriving the uninhabited type or
            # reporting an error both seem unhelpful. So we abort.
            if not new_parent_types:
                return output

            expr = parent_expr
            expr_type = output[parent_expr] = make_simplified_union(new_parent_types)

    def refine_identity_comparison_expression(self,
                                              operands: List[Expression],
                                              operand_types: List[Type],
                                              chain_indices: List[int],
                                              narrowable_operand_indices: AbstractSet[int],
                                              is_valid_target: Callable[[ProperType], bool],
                                              coerce_only_in_literal_context: bool,
                                              ) -> Tuple[TypeMap, TypeMap]:
        """Produce conditional type maps refining expressions by an identity/equality comparison.

        The 'operands' and 'operand_types' lists should be the full list of operands used
        in the overall comparison expression. The 'chain_indices' list is the list of indices
        actually used within this identity comparison chain.

        So if we have the expression:

            a <= b is c is d <= e

        ...then 'operands' and 'operand_types' would be lists of length 5 and 'chain_indices'
        would be the list [1, 2, 3].

        The 'narrowable_operand_indices' parameter is the set of all indices we are allowed
        to refine the types of: that is, all operands that will potentially be a part of
        the output TypeMaps.

        Although this function could theoretically try setting the types of the operands
        in the chains to the meet, doing that causes too many issues in real-world code.
        Instead, we use 'is_valid_target' to identify which of the given chain types
        we could plausibly use as the refined type for the expressions in the chain.

        Similarly, 'coerce_only_in_literal_context' controls whether we should try coercing
        expressions in the chain to a Literal type. Performing this coercion is sometimes
        too aggressive of a narrowing, depending on context.
        """
        should_coerce = True
        if coerce_only_in_literal_context:
            should_coerce = any(is_literal_type_like(operand_types[i]) for i in chain_indices)

        target: Optional[Type] = None
        possible_target_indices = []
        for i in chain_indices:
            expr_type = operand_types[i]
            if should_coerce:
                expr_type = coerce_to_literal(expr_type)
            if not is_valid_target(get_proper_type(expr_type)):
                continue
            if target and not is_same_type(target, expr_type):
                # We have multiple disjoint target types. So the 'if' branch
                # must be unreachable.
                return None, {}
            target = expr_type
            possible_target_indices.append(i)

        # There's nothing we can currently infer if none of the operands are valid targets,
        # so we end early and infer nothing.
        if target is None:
            return {}, {}

        # If possible, use an unassignable expression as the target.
        # We skip refining the type of the target below, so ideally we'd
        # want to pick an expression we were going to skip anyways.
        singleton_index = -1
        for i in possible_target_indices:
            if i not in narrowable_operand_indices:
                singleton_index = i

        # But if none of the possible singletons are unassignable ones, we give up
        # and arbitrarily pick the last item, mostly because other parts of the
        # type narrowing logic bias towards picking the rightmost item and it'd be
        # nice to stay consistent.
        #
        # That said, it shouldn't matter which index we pick. For example, suppose we
        # have this if statement, where 'x' and 'y' both have singleton types:
        #
        #     if x is y:
        #         reveal_type(x)
        #         reveal_type(y)
        #     else:
        #         reveal_type(x)
        #         reveal_type(y)
        #
        # At this point, 'x' and 'y' *must* have the same singleton type: we would have
        # ended early in the first for-loop in this function if they weren't.
        #
        # So, we should always get the same result in the 'if' case no matter which
        # index we pick. And while we do end up getting different results in the 'else'
        # case depending on the index (e.g. if we pick 'y', then its type stays the same
        # while 'x' is narrowed to '<uninhabited>'), this distinction is also moot: mypy
        # currently will just mark the whole branch as unreachable if either operand is
        # narrowed to <uninhabited>.
        if singleton_index == -1:
            singleton_index = possible_target_indices[-1]

        sum_type_name = None
        target = get_proper_type(target)
        if (isinstance(target, LiteralType) and
                (target.is_enum_literal() or isinstance(target.value, bool))):
            sum_type_name = target.fallback.type.fullname

        target_type = [TypeRange(target, is_upper_bound=False)]

        partial_type_maps = []
        for i in chain_indices:
            # If we try refining a type against itself, conditional_type_map
            # will end up assuming that the 'else' branch is unreachable. This is
            # typically not what we want: generally the user will intend for the
            # target type to be some fixed 'sentinel' value and will want to refine
            # the other exprs against this one instead.
            if i == singleton_index:
                continue

            # Naturally, we can't refine operands which are not permitted to be refined.
            if i not in narrowable_operand_indices:
                continue

            expr = operands[i]
            expr_type = coerce_to_literal(operand_types[i])

            if sum_type_name is not None:
                expr_type = try_expanding_sum_type_to_union(expr_type, sum_type_name)

            # We intentionally use 'conditional_type_map' directly here instead of
            # 'self.conditional_type_map_with_intersection': we only compute ad-hoc
            # intersections when working with pure instances.
            partial_type_maps.append(conditional_type_map(expr, expr_type, target_type))

        return reduce_conditional_maps(partial_type_maps)

    def refine_away_none_in_comparison(self,
                                       operands: List[Expression],
                                       operand_types: List[Type],
                                       chain_indices: List[int],
                                       narrowable_operand_indices: AbstractSet[int],
                                       ) -> Tuple[TypeMap, TypeMap]:
        """Produces conditional type maps refining away None in an identity/equality chain.

        For more details about what the different arguments mean, see the
        docstring of 'refine_identity_comparison_expression' up above.
        """
        non_optional_types = []
        for i in chain_indices:
            typ = operand_types[i]
            if not is_optional(typ):
                non_optional_types.append(typ)

        # Make sure we have a mixture of optional and non-optional types.
        if len(non_optional_types) == 0 or len(non_optional_types) == len(chain_indices):
            return {}, {}

        if_map = {}
        for i in narrowable_operand_indices:
            expr_type = operand_types[i]
            if not is_optional(expr_type):
                continue
            if any(is_overlapping_erased_types(expr_type, t) for t in non_optional_types):
                if_map[operands[i]] = remove_optional(expr_type)

        return if_map, {}

    #
    # Helpers
    #

    def check_subtype(self,
                      subtype: Type,
                      supertype: Type,
                      context: Context,
                      msg: Union[str, ErrorMessage] = message_registry.INCOMPATIBLE_TYPES,
                      subtype_label: Optional[str] = None,
                      supertype_label: Optional[str] = None,
                      *,
                      code: Optional[ErrorCode] = None,
                      outer_context: Optional[Context] = None) -> bool:
        """Generate an error if the subtype is not compatible with supertype."""
        if is_subtype(subtype, supertype):
            return True

        if isinstance(msg, ErrorMessage):
            msg_text = msg.value
            code = msg.code
        else:
            msg_text = msg
        subtype = get_proper_type(subtype)
        supertype = get_proper_type(supertype)
        if self.msg.try_report_long_tuple_assignment_error(subtype, supertype, context, msg_text,
                                       subtype_label, supertype_label, code=code):
            return False
        if self.should_suppress_optional_error([subtype]):
            return False
        extra_info: List[str] = []
        note_msg = ''
        notes: List[str] = []
        if subtype_label is not None or supertype_label is not None:
            subtype_str, supertype_str = format_type_distinctly(subtype, supertype)
            if subtype_label is not None:
                extra_info.append(subtype_label + ' ' + subtype_str)
            if supertype_label is not None:
                extra_info.append(supertype_label + ' ' + supertype_str)
            note_msg = make_inferred_type_note(outer_context or context, subtype,
                                               supertype, supertype_str)
            if isinstance(subtype, Instance) and isinstance(supertype, Instance):
                notes = append_invariance_notes([], subtype, supertype)
        if extra_info:
            msg_text += ' (' + ', '.join(extra_info) + ')'

        self.fail(ErrorMessage(msg_text, code=code), context)
        for note in notes:
            self.msg.note(note, context, code=code)
        if note_msg:
            self.note(note_msg, context, code=code)
        if (isinstance(supertype, Instance) and supertype.type.is_protocol and
                isinstance(subtype, (Instance, TupleType, TypedDictType))):
            self.msg.report_protocol_problems(subtype, supertype, context, code=code)
        if isinstance(supertype, CallableType) and isinstance(subtype, Instance):
            call = find_member('__call__', subtype, subtype, is_operator=True)
            if call:
                self.msg.note_call(subtype, call, context, code=code)
        if isinstance(subtype, (CallableType, Overloaded)) and isinstance(supertype, Instance):
            if supertype.type.is_protocol and supertype.type.protocol_members == ['__call__']:
                call = find_member('__call__', supertype, subtype, is_operator=True)
                assert call is not None
                self.msg.note_call(supertype, call, context, code=code)
        return False

    def contains_none(self, t: Type) -> bool:
        t = get_proper_type(t)
        return (
            isinstance(t, NoneType) or
            (isinstance(t, UnionType) and any(self.contains_none(ut) for ut in t.items)) or
            (isinstance(t, TupleType) and any(self.contains_none(tt) for tt in t.items)) or
            (isinstance(t, Instance) and bool(t.args)
             and any(self.contains_none(it) for it in t.args))
        )

    def should_suppress_optional_error(self, related_types: List[Type]) -> bool:
        return self.suppress_none_errors and any(self.contains_none(t) for t in related_types)

    def named_type(self, name: str) -> Instance:
        """Return an instance type with given name and implicit Any type args.

        For example, named_type('builtins.object') produces the 'object' type.
        """
        # Assume that the name refers to a type.
        sym = self.lookup_qualified(name)
        node = sym.node
        if isinstance(node, TypeAlias):
            assert isinstance(node.target, Instance)  # type: ignore
            node = node.target.type
        assert isinstance(node, TypeInfo)
        any_type = AnyType(TypeOfAny.from_omitted_generics)
        return Instance(node, [any_type] * len(node.defn.type_vars))

    def named_generic_type(self, name: str, args: List[Type]) -> Instance:
        """Return an instance with the given name and type arguments.

        Assume that the number of arguments is correct.  Assume that
        the name refers to a compatible generic type.
        """
        info = self.lookup_typeinfo(name)
        args = [remove_instance_last_known_values(arg) for arg in args]
        # TODO: assert len(args) == len(info.defn.type_vars)
        return Instance(info, args)

    def lookup_typeinfo(self, fullname: str) -> TypeInfo:
        # Assume that the name refers to a class.
        sym = self.lookup_qualified(fullname)
        node = sym.node
        assert isinstance(node, TypeInfo)
        return node

    def type_type(self) -> Instance:
        """Return instance type 'type'."""
        return self.named_type('builtins.type')

    def str_type(self) -> Instance:
        """Return instance type 'str'."""
        return self.named_type('builtins.str')

    def store_type(self, node: Expression, typ: Type) -> None:
        """Store the type of a node in the type map."""
        self.type_map[node] = typ

    def in_checked_function(self) -> bool:
        """Should we type-check the current function?

        - Yes if --check-untyped-defs is set.
        - Yes outside functions.
        - Yes in annotated functions.
        - No otherwise.
        """
        return (self.options.check_untyped_defs
                or not self.dynamic_funcs
                or not self.dynamic_funcs[-1])

    def lookup(self, name: str) -> SymbolTableNode:
        """Look up a definition from the symbol table with the given name.
        """
        if name in self.globals:
            return self.globals[name]
        else:
            b = self.globals.get('__builtins__', None)
            if b:
                table = cast(MypyFile, b.node).names
                if name in table:
                    return table[name]
            raise KeyError('Failed lookup: {}'.format(name))

    def lookup_qualified(self, name: str) -> SymbolTableNode:
        if '.' not in name:
            return self.lookup(name)
        else:
            parts = name.split('.')
            n = self.modules[parts[0]]
            for i in range(1, len(parts) - 1):
                sym = n.names.get(parts[i])
                assert sym is not None, "Internal error: attempted lookup of unknown name"
                n = cast(MypyFile, sym.node)
            last = parts[-1]
            if last in n.names:
                return n.names[last]
            elif len(parts) == 2 and parts[0] == 'builtins':
                fullname = 'builtins.' + last
                if fullname in SUGGESTED_TEST_FIXTURES:
                    suggestion = ", e.g. add '[builtins fixtures/{}]' to your test".format(
                        SUGGESTED_TEST_FIXTURES[fullname])
                else:
                    suggestion = ''
                raise KeyError("Could not find builtin symbol '{}' (If you are running a "
                               "test case, use a fixture that "
                               "defines this symbol{})".format(last, suggestion))
            else:
                msg = "Failed qualified lookup: '{}' (fullname = '{}')."
                raise KeyError(msg.format(last, name))

    @contextmanager
    def enter_partial_types(self, *, is_function: bool = False,
                            is_class: bool = False) -> Iterator[None]:
        """Enter a new scope for collecting partial types.

        Also report errors for (some) variables which still have partial
        types, i.e. we couldn't infer a complete type.
        """
        is_local = (self.partial_types and self.partial_types[-1].is_local) or is_function
        self.partial_types.append(PartialTypeScope({}, is_function, is_local))
        yield

        # Don't complain about not being able to infer partials if it is
        # at the toplevel (with allow_untyped_globals) or if it is in an
        # untyped function being checked with check_untyped_defs.
        permissive = (self.options.allow_untyped_globals and not is_local) or (
            self.options.check_untyped_defs
            and self.dynamic_funcs
            and self.dynamic_funcs[-1]
        )

        partial_types, _, _ = self.partial_types.pop()
        if not self.current_node_deferred:
            for var, context in partial_types.items():
                # If we require local partial types, there are a few exceptions where
                # we fall back to inferring just "None" as the type from a None initializer:
                #
                # 1. If all happens within a single function this is acceptable, since only
                #    the topmost function is a separate target in fine-grained incremental mode.
                #    We primarily want to avoid "splitting" partial types across targets.
                #
                # 2. A None initializer in the class body if the attribute is defined in a base
                #    class is fine, since the attribute is already defined and it's currently okay
                #    to vary the type of an attribute covariantly. The None type will still be
                #    checked for compatibility with base classes elsewhere. Without this exception
                #    mypy could require an annotation for an attribute that already has been
                #    declared in a base class, which would be bad.
                allow_none = (not self.options.local_partial_types
                              or is_function
                              or (is_class and self.is_defined_in_base_class(var)))
                if (allow_none
                        and isinstance(var.type, PartialType)
                        and var.type.type is None
                        and not permissive):
                    var.type = NoneType()
                else:
                    if var not in self.partial_reported and not permissive:
                        self.msg.need_annotation_for_var(var, context, self.options.python_version)
                        self.partial_reported.add(var)
                    if var.type:
                        var.type = self.fixup_partial_type(var.type)

    def handle_partial_var_type(
            self, typ: PartialType, is_lvalue: bool, node: Var, context: Context) -> Type:
        """Handle a reference to a partial type through a var.

        (Used by checkexpr and checkmember.)
        """
        in_scope, is_local, partial_types = self.find_partial_types_in_all_scopes(node)
        if typ.type is None and in_scope:
            # 'None' partial type. It has a well-defined type. In an lvalue context
            # we want to preserve the knowledge of it being a partial type.
            if not is_lvalue:
                return NoneType()
            else:
                return typ
        else:
            if partial_types is not None and not self.current_node_deferred:
                if in_scope:
                    context = partial_types[node]
                    if is_local or not self.options.allow_untyped_globals:
                        self.msg.need_annotation_for_var(node, context,
                                                         self.options.python_version)
                        self.partial_reported.add(node)
                else:
                    # Defer the node -- we might get a better type in the outer scope
                    self.handle_cannot_determine_type(node.name, context)
            return self.fixup_partial_type(typ)

    def fixup_partial_type(self, typ: Type) -> Type:
        """Convert a partial type that we couldn't resolve into something concrete.

        This means, for None we make it Optional[Any], and for anything else we
        fill in all of the type arguments with Any.
        """
        if not isinstance(typ, PartialType):
            return typ
        if typ.type is None:
            return UnionType.make_union([AnyType(TypeOfAny.unannotated), NoneType()])
        else:
            return Instance(
                typ.type,
                [AnyType(TypeOfAny.unannotated)] * len(typ.type.type_vars))

    def is_defined_in_base_class(self, var: Var) -> bool:
        if var.info:
            for base in var.info.mro[1:]:
                if base.get(var.name) is not None:
                    return True
            if var.info.fallback_to_any:
                return True
        return False

    def find_partial_types(self, var: Var) -> Optional[Dict[Var, Context]]:
        """Look for an active partial type scope containing variable.

        A scope is active if assignments in the current context can refine a partial
        type originally defined in the scope. This is affected by the local_partial_types
        configuration option.
        """
        in_scope, _, partial_types = self.find_partial_types_in_all_scopes(var)
        if in_scope:
            return partial_types
        return None

    def find_partial_types_in_all_scopes(
            self, var: Var) -> Tuple[bool, bool, Optional[Dict[Var, Context]]]:
        """Look for partial type scope containing variable.

        Return tuple (is the scope active, is the scope a local scope, scope).
        """
        for scope in reversed(self.partial_types):
            if var in scope.map:
                # All scopes within the outermost function are active. Scopes out of
                # the outermost function are inactive to allow local reasoning (important
                # for fine-grained incremental mode).
                disallow_other_scopes = self.options.local_partial_types

                if isinstance(var.type, PartialType) and var.type.type is not None and var.info:
                    # This is an ugly hack to make partial generic self attributes behave
                    # as if --local-partial-types is always on (because it used to be like this).
                    disallow_other_scopes = True

                scope_active = (not disallow_other_scopes
                                or scope.is_local == self.partial_types[-1].is_local)
                return scope_active, scope.is_local, scope.map
        return False, False, None

    def temp_node(self, t: Type, context: Optional[Context] = None) -> TempNode:
        """Create a temporary node with the given, fixed type."""
        return TempNode(t, context=context)

    def fail(self, msg: Union[str, ErrorMessage], context: Context, *,
             code: Optional[ErrorCode] = None) -> None:
        """Produce an error message."""
        if isinstance(msg, ErrorMessage):
            self.msg.fail(msg.value, context, code=msg.code)
            return
        self.msg.fail(msg, context, code=code)

    def note(self,
             msg: str,
             context: Context,
             offset: int = 0,
             *,
             code: Optional[ErrorCode] = None) -> None:
        """Produce a note."""
        self.msg.note(msg, context, offset=offset, code=code)

    def iterable_item_type(self, instance: Instance) -> Type:
        iterable = map_instance_to_supertype(
            instance,
            self.lookup_typeinfo('typing.Iterable'))
        item_type = iterable.args[0]
        if not isinstance(get_proper_type(item_type), AnyType):
            # This relies on 'map_instance_to_supertype' returning 'Iterable[Any]'
            # in case there is no explicit base class.
            return item_type
        # Try also structural typing.
        iter_type = get_proper_type(find_member('__iter__', instance, instance, is_operator=True))
        if iter_type and isinstance(iter_type, CallableType):
            ret_type = get_proper_type(iter_type.ret_type)
            if isinstance(ret_type, Instance):
                iterator = map_instance_to_supertype(ret_type,
                                                     self.lookup_typeinfo('typing.Iterator'))
                item_type = iterator.args[0]
        return item_type

    def function_type(self, func: FuncBase) -> FunctionLike:
        return function_type(func, self.named_type('builtins.function'))

    def push_type_map(self, type_map: 'TypeMap') -> None:
        if type_map is None:
            self.binder.unreachable()
        else:
            for expr, type in type_map.items():
                self.binder.put(expr, type)

    def infer_issubclass_maps(self, node: CallExpr,
                              expr: Expression,
                              type_map: Dict[Expression, Type]
                              ) -> Tuple[TypeMap, TypeMap]:
        """Infer type restrictions for an expression in issubclass call."""
        vartype = type_map[expr]
        type = get_isinstance_type(node.args[1], type_map)
        if isinstance(vartype, TypeVarType):
            vartype = vartype.upper_bound
        vartype = get_proper_type(vartype)
        if isinstance(vartype, UnionType):
            union_list = []
            for t in get_proper_types(vartype.items):
                if isinstance(t, TypeType):
                    union_list.append(t.item)
                else:
                    # This is an error that should be reported earlier
                    # if we reach here, we refuse to do any type inference.
                    return {}, {}
            vartype = UnionType(union_list)
        elif isinstance(vartype, TypeType):
            vartype = vartype.item
        elif (isinstance(vartype, Instance) and
                vartype.type.fullname == 'builtins.type'):
            vartype = self.named_type('builtins.object')
        else:
            # Any other object whose type we don't know precisely
            # for example, Any or a custom metaclass.
            return {}, {}  # unknown type
        yes_map, no_map = self.conditional_type_map_with_intersection(expr, vartype, type)
        yes_map, no_map = map(convert_to_typetype, (yes_map, no_map))
        return yes_map, no_map

    def conditional_type_map_with_intersection(self,
                                               expr: Expression,
                                               expr_type: Type,
                                               type_ranges: Optional[List[TypeRange]],
                                               ) -> Tuple[TypeMap, TypeMap]:
        # For some reason, doing "yes_map, no_map = conditional_type_maps(...)"
        # doesn't work: mypyc will decide that 'yes_map' is of type None if we try.
        initial_maps = conditional_type_map(expr, expr_type, type_ranges)
        yes_map: TypeMap = initial_maps[0]
        no_map: TypeMap = initial_maps[1]

        if yes_map is not None or type_ranges is None:
            return yes_map, no_map

        # If conditions_type_map was unable to successfully narrow the expr_type
        # using the type_ranges and concluded if-branch is unreachable, we try
        # computing it again using a different algorithm that tries to generate
        # an ad-hoc intersection between the expr_type and the type_ranges.
        expr_type = get_proper_type(expr_type)
        if isinstance(expr_type, UnionType):
            possible_expr_types = get_proper_types(expr_type.relevant_items())
        else:
            possible_expr_types = [expr_type]

        possible_target_types = []
        for tr in type_ranges:
            item = get_proper_type(tr.item)
            if not isinstance(item, Instance) or tr.is_upper_bound:
                return yes_map, no_map
            possible_target_types.append(item)

        out = []
        for v in possible_expr_types:
            if not isinstance(v, Instance):
                return yes_map, no_map
            for t in possible_target_types:
                intersection = self.intersect_instances((v, t), expr)
                if intersection is None:
                    continue
                out.append(intersection)
        if len(out) == 0:
            return None, {}
        new_yes_type = make_simplified_union(out)
        return {expr: new_yes_type}, {}

    def is_writable_attribute(self, node: Node) -> bool:
        """Check if an attribute is writable"""
        if isinstance(node, Var):
            return True
        elif isinstance(node, OverloadedFuncDef) and node.is_property:
            first_item = cast(Decorator, node.items[0])
            return first_item.var.is_settable_property
        else:
            return False


def conditional_type_map(expr: Expression,
                         current_type: Optional[Type],
                         proposed_type_ranges: Optional[List[TypeRange]],
                         ) -> Tuple[TypeMap, TypeMap]:
    """Takes in an expression, the current type of the expression, and a
    proposed type of that expression.

    Returns a 2-tuple: The first element is a map from the expression to
    the proposed type, if the expression can be the proposed type.  The
    second element is a map from the expression to the type it would hold
    if it was not the proposed type, if any. None means bot, {} means top"""
    if proposed_type_ranges:
        proposed_items = [type_range.item for type_range in proposed_type_ranges]
        proposed_type = make_simplified_union(proposed_items)
        if current_type:
            if isinstance(proposed_type, AnyType):
                # We don't really know much about the proposed type, so we shouldn't
                # attempt to narrow anything. Instead, we broaden the expr to Any to
                # avoid false positives
                return {expr: proposed_type}, {}
            elif (not any(type_range.is_upper_bound for type_range in proposed_type_ranges)
               and is_proper_subtype(current_type, proposed_type)):
                # Expression is always of one of the types in proposed_type_ranges
                return {}, None
            elif not is_overlapping_types(current_type, proposed_type,
                                          prohibit_none_typevar_overlap=True):
                # Expression is never of any type in proposed_type_ranges
                return None, {}
            else:
                # we can only restrict when the type is precise, not bounded
                proposed_precise_type = UnionType.make_union([
                    type_range.item
                    for type_range in proposed_type_ranges
                    if not type_range.is_upper_bound
                ])
                remaining_type = restrict_subtype_away(current_type, proposed_precise_type)
                return {expr: proposed_type}, {expr: remaining_type}
        else:
            return {expr: proposed_type}, {}
    else:
        # An isinstance check, but we don't understand the type
        return {}, {}


def gen_unique_name(base: str, table: SymbolTable) -> str:
    """Generate a name that does not appear in table by appending numbers to base."""
    if base not in table:
        return base
    i = 1
    while base + str(i) in table:
        i += 1
    return base + str(i)


def is_true_literal(n: Expression) -> bool:
    """Returns true if this expression is the 'True' literal/keyword."""
    return (refers_to_fullname(n, 'builtins.True')
            or isinstance(n, IntExpr) and n.value != 0)


def is_false_literal(n: Expression) -> bool:
    """Returns true if this expression is the 'False' literal/keyword."""
    return (refers_to_fullname(n, 'builtins.False')
            or isinstance(n, IntExpr) and n.value == 0)


def is_literal_enum(type_map: Mapping[Expression, Type], n: Expression) -> bool:
    """Returns true if this expression (with the given type context) is an Enum literal.

    For example, if we had an enum:

        class Foo(Enum):
            A = 1
            B = 2

    ...and if the expression 'Foo' referred to that enum within the current type context,
    then the expression 'Foo.A' would be a literal enum. However, if we did 'a = Foo.A',
    then the variable 'a' would *not* be a literal enum.

    We occasionally special-case expressions like 'Foo.A' and treat them as a single primitive
    unit for the same reasons we sometimes treat 'True', 'False', or 'None' as a single
    primitive unit.
    """
    if not isinstance(n, MemberExpr) or not isinstance(n.expr, NameExpr):
        return False

    parent_type = type_map.get(n.expr)
    member_type = type_map.get(n)
    if member_type is None or parent_type is None:
        return False

    parent_type = get_proper_type(parent_type)
    member_type = get_proper_type(coerce_to_literal(member_type))
    if not isinstance(parent_type, FunctionLike) or not isinstance(member_type, LiteralType):
        return False

    if not parent_type.is_type_obj():
        return False

    return member_type.is_enum_literal() and member_type.fallback.type == parent_type.type_object()


def is_literal_none(n: Expression) -> bool:
    """Returns true if this expression is the 'None' literal/keyword."""
    return isinstance(n, NameExpr) and n.fullname == 'builtins.None'


def is_literal_not_implemented(n: Expression) -> bool:
    return isinstance(n, NameExpr) and n.fullname == 'builtins.NotImplemented'


def builtin_item_type(tp: Type) -> Optional[Type]:
    """Get the item type of a builtin container.

    If 'tp' is not one of the built containers (these includes NamedTuple and TypedDict)
    or if the container is not parameterized (like List or List[Any])
    return None. This function is used to narrow optional types in situations like this:

        x: Optional[int]
        if x in (1, 2, 3):
            x + 42  # OK

    Note: this is only OK for built-in containers, where we know the behavior
    of __contains__.
    """
    tp = get_proper_type(tp)

    if isinstance(tp, Instance):
        if tp.type.fullname in [
            'builtins.list', 'builtins.tuple', 'builtins.dict',
            'builtins.set', 'builtins.frozenset',
        ]:
            if not tp.args:
                # TODO: fix tuple in lib-stub/builtins.pyi (it should be generic).
                return None
            if not isinstance(get_proper_type(tp.args[0]), AnyType):
                return tp.args[0]
    elif isinstance(tp, TupleType) and all(not isinstance(it, AnyType)
                                           for it in get_proper_types(tp.items)):
        return make_simplified_union(tp.items)  # this type is not externally visible
    elif isinstance(tp, TypedDictType):
        # TypedDict always has non-optional string keys. Find the key type from the Mapping
        # base class.
        for base in tp.fallback.type.mro:
            if base.fullname == 'typing.Mapping':
                return map_instance_to_supertype(tp.fallback, base).args[0]
        assert False, 'No Mapping base class found for TypedDict fallback'
    return None


def and_conditional_maps(m1: TypeMap, m2: TypeMap) -> TypeMap:
    """Calculate what information we can learn from the truth of (e1 and e2)
    in terms of the information that we can learn from the truth of e1 and
    the truth of e2.
    """

    if m1 is None or m2 is None:
        # One of the conditions can never be true.
        return None
    # Both conditions can be true; combine the information. Anything
    # we learn from either conditions's truth is valid. If the same
    # expression's type is refined by both conditions, we somewhat
    # arbitrarily give precedence to m2. (In the future, we could use
    # an intersection type.)
    result = m2.copy()
    m2_keys = set(literal_hash(n2) for n2 in m2)
    for n1 in m1:
        if literal_hash(n1) not in m2_keys:
            result[n1] = m1[n1]
    return result


def or_conditional_maps(m1: TypeMap, m2: TypeMap) -> TypeMap:
    """Calculate what information we can learn from the truth of (e1 or e2)
    in terms of the information that we can learn from the truth of e1 and
    the truth of e2.
    """

    if m1 is None:
        return m2
    if m2 is None:
        return m1
    # Both conditions can be true. Combine information about
    # expressions whose type is refined by both conditions. (We do not
    # learn anything about expressions whose type is refined by only
    # one condition.)
    result: Dict[Expression, Type] = {}
    for n1 in m1:
        for n2 in m2:
            if literal_hash(n1) == literal_hash(n2):
                result[n1] = make_simplified_union([m1[n1], m2[n2]])
    return result


def reduce_conditional_maps(type_maps: List[Tuple[TypeMap, TypeMap]],
                            ) -> Tuple[TypeMap, TypeMap]:
    """Reduces a list containing pairs of if/else TypeMaps into a single pair.

    We "and" together all of the if TypeMaps and "or" together the else TypeMaps. So
    for example, if we had the input:

        [
            ({x: TypeIfX, shared: TypeIfShared1}, {x: TypeElseX, shared: TypeElseShared1}),
            ({y: TypeIfY, shared: TypeIfShared2}, {y: TypeElseY, shared: TypeElseShared2}),
        ]

    ...we'd return the output:

        (
            {x: TypeIfX,   y: TypeIfY,   shared: PseudoIntersection[TypeIfShared1, TypeIfShared2]},
            {shared: Union[TypeElseShared1, TypeElseShared2]},
        )

    ...where "PseudoIntersection[X, Y] == Y" because mypy actually doesn't understand intersections
    yet, so we settle for just arbitrarily picking the right expr's type.

    We only retain the shared expression in the 'else' case because we don't actually know
    whether x was refined or y was refined -- only just that one of the two was refined.
    """
    if len(type_maps) == 0:
        return {}, {}
    elif len(type_maps) == 1:
        return type_maps[0]
    else:
        final_if_map, final_else_map = type_maps[0]
        for if_map, else_map in type_maps[1:]:
            final_if_map = and_conditional_maps(final_if_map, if_map)
            final_else_map = or_conditional_maps(final_else_map, else_map)

        return final_if_map, final_else_map


def convert_to_typetype(type_map: TypeMap) -> TypeMap:
    converted_type_map: Dict[Expression, Type] = {}
    if type_map is None:
        return None
    for expr, typ in type_map.items():
        t = typ
        if isinstance(t, TypeVarType):
            t = t.upper_bound
        # TODO: should we only allow unions of instances as per PEP 484?
        if not isinstance(get_proper_type(t), (UnionType, Instance)):
            # unknown type; error was likely reported earlier
            return {}
        converted_type_map[expr] = TypeType.make_normalized(typ)
    return converted_type_map


def flatten(t: Expression) -> List[Expression]:
    """Flatten a nested sequence of tuples/lists into one list of nodes."""
    if isinstance(t, TupleExpr) or isinstance(t, ListExpr):
        return [b for a in t.items for b in flatten(a)]
    elif isinstance(t, StarExpr):
        return flatten(t.expr)
    else:
        return [t]


def flatten_types(t: Type) -> List[Type]:
    """Flatten a nested sequence of tuples into one list of nodes."""
    t = get_proper_type(t)
    if isinstance(t, TupleType):
        return [b for a in t.items for b in flatten_types(a)]
    else:
        return [t]


def get_isinstance_type(expr: Expression,
                        type_map: Dict[Expression, Type]) -> Optional[List[TypeRange]]:
    if isinstance(expr, OpExpr) and expr.op == '|':
        left = get_isinstance_type(expr.left, type_map)
        right = get_isinstance_type(expr.right, type_map)
        if left is None or right is None:
            return None
        return left + right
    all_types = get_proper_types(flatten_types(type_map[expr]))
    types: List[TypeRange] = []
    for typ in all_types:
        if isinstance(typ, FunctionLike) and typ.is_type_obj():
            # Type variables may be present -- erase them, which is the best
            # we can do (outside disallowing them here).
            erased_type = erase_typevars(typ.items[0].ret_type)
            types.append(TypeRange(erased_type, is_upper_bound=False))
        elif isinstance(typ, TypeType):
            # Type[A] means "any type that is a subtype of A" rather than "precisely type A"
            # we indicate this by setting is_upper_bound flag
            types.append(TypeRange(typ.item, is_upper_bound=True))
        elif isinstance(typ, Instance) and typ.type.fullname == 'builtins.type':
            object_type = Instance(typ.type.mro[-1], [])
            types.append(TypeRange(object_type, is_upper_bound=True))
        elif isinstance(typ, AnyType):
            types.append(TypeRange(typ, is_upper_bound=False))
        else:  # we didn't see an actual type, but rather a variable whose value is unknown to us
            return None
    if not types:
        # this can happen if someone has empty tuple as 2nd argument to isinstance
        # strictly speaking, we should return UninhabitedType but for simplicity we will simply
        # refuse to do any type inference for now
        return None
    return types


def expand_func(defn: FuncItem, map: Dict[TypeVarId, Type]) -> FuncItem:
    visitor = TypeTransformVisitor(map)
    ret = defn.accept(visitor)
    assert isinstance(ret, FuncItem)
    return ret


class TypeTransformVisitor(TransformVisitor):
    def __init__(self, map: Dict[TypeVarId, Type]) -> None:
        super().__init__()
        self.map = map

    def type(self, type: Type) -> Type:
        return expand_type(type, self.map)


def are_argument_counts_overlapping(t: CallableType, s: CallableType) -> bool:
    """Can a single call match both t and s, based just on positional argument counts?
    """
    min_args = max(t.min_args, s.min_args)
    max_args = min(t.max_possible_positional_args(), s.max_possible_positional_args())
    return min_args <= max_args


def is_unsafe_overlapping_overload_signatures(signature: CallableType,
                                              other: CallableType) -> bool:
    """Check if two overloaded signatures are unsafely overlapping or partially overlapping.

    We consider two functions 's' and 't' to be unsafely overlapping if both
    of the following are true:

    1.  s's parameters are all more precise or partially overlapping with t's
    2.  s's return type is NOT a subtype of t's.

    Assumes that 'signature' appears earlier in the list of overload
    alternatives then 'other' and that their argument counts are overlapping.
    """
    # Try detaching callables from the containing class so that all TypeVars
    # are treated as being free.
    #
    # This lets us identify cases where the two signatures use completely
    # incompatible types -- e.g. see the testOverloadingInferUnionReturnWithMixedTypevars
    # test case.
    signature = detach_callable(signature)
    other = detach_callable(other)

    # Note: We repeat this check twice in both directions due to a slight
    # asymmetry in 'is_callable_compatible'. When checking for partial overlaps,
    # we attempt to unify 'signature' and 'other' both against each other.
    #
    # If 'signature' cannot be unified with 'other', we end early. However,
    # if 'other' cannot be modified with 'signature', the function continues
    # using the older version of 'other'.
    #
    # This discrepancy is unfortunately difficult to get rid of, so we repeat the
    # checks twice in both directions for now.
    return (is_callable_compatible(signature, other,
                                  is_compat=is_overlapping_types_no_promote,
                                  is_compat_return=lambda l, r: not is_subtype_no_promote(l, r),
                                  ignore_return=False,
                                  check_args_covariantly=True,
                                  allow_partial_overlap=True) or
            is_callable_compatible(other, signature,
                                   is_compat=is_overlapping_types_no_promote,
                                   is_compat_return=lambda l, r: not is_subtype_no_promote(r, l),
                                   ignore_return=False,
                                   check_args_covariantly=False,
                                   allow_partial_overlap=True))


def detach_callable(typ: CallableType) -> CallableType:
    """Ensures that the callable's type variables are 'detached' and independent of the context.

    A callable normally keeps track of the type variables it uses within its 'variables' field.
    However, if the callable is from a method and that method is using a class type variable,
    the callable will not keep track of that type variable since it belongs to the class.

    This function will traverse the callable and find all used type vars and add them to the
    variables field if it isn't already present.

    The caller can then unify on all type variables whether or not the callable is originally
    from a class or not."""
    type_list = typ.arg_types + [typ.ret_type]

    appear_map: Dict[str, List[int]] = {}
    for i, inner_type in enumerate(type_list):
        typevars_available = get_type_vars(inner_type)
        for var in typevars_available:
            if var.fullname not in appear_map:
                appear_map[var.fullname] = []
            appear_map[var.fullname].append(i)

    used_type_var_names = set()
    for var_name, appearances in appear_map.items():
        used_type_var_names.add(var_name)

    all_type_vars = get_type_vars(typ)
    new_variables = []
    for var in set(all_type_vars):
        if var.fullname not in used_type_var_names:
            continue
        new_variables.append(TypeVarType(
            name=var.name,
            fullname=var.fullname,
            id=var.id,
            values=var.values,
            upper_bound=var.upper_bound,
            variance=var.variance,
        ))
    out = typ.copy_modified(
        variables=new_variables,
        arg_types=type_list[:-1],
        ret_type=type_list[-1],
    )
    return out


def overload_can_never_match(signature: CallableType, other: CallableType) -> bool:
    """Check if the 'other' method can never be matched due to 'signature'.

    This can happen if signature's parameters are all strictly broader then
    other's parameters.

    Assumes that both signatures have overlapping argument counts.
    """
    # The extra erasure is needed to prevent spurious errors
    # in situations where an `Any` overload is used as a fallback
    # for an overload with type variables. The spurious error appears
    # because the type variables turn into `Any` during unification in
    # the below subtype check and (surprisingly?) `is_proper_subtype(Any, Any)`
    # returns `True`.
    # TODO: find a cleaner solution instead of this ad-hoc erasure.
    exp_signature = expand_type(signature, {tvar.id: erase_def_to_union_or_bound(tvar)
                                for tvar in signature.variables})
    assert isinstance(exp_signature, ProperType)
    assert isinstance(exp_signature, CallableType)
    return is_callable_compatible(exp_signature, other,
                                  is_compat=is_more_precise,
                                  ignore_return=True)


def is_more_general_arg_prefix(t: FunctionLike, s: FunctionLike) -> bool:
    """Does t have wider arguments than s?"""
    # TODO should an overload with additional items be allowed to be more
    #      general than one with fewer items (or just one item)?
    if isinstance(t, CallableType):
        if isinstance(s, CallableType):
            return is_callable_compatible(t, s,
                                          is_compat=is_proper_subtype,
                                          ignore_return=True)
    elif isinstance(t, FunctionLike):
        if isinstance(s, FunctionLike):
            if len(t.items) == len(s.items):
                return all(is_same_arg_prefix(items, itemt)
                           for items, itemt in zip(t.items, s.items))
    return False


def is_same_arg_prefix(t: CallableType, s: CallableType) -> bool:
    return is_callable_compatible(t, s,
                                  is_compat=is_same_type,
                                  ignore_return=True,
                                  check_args_covariantly=True,
                                  ignore_pos_arg_names=True)


def infer_operator_assignment_method(typ: Type, operator: str) -> Tuple[bool, str]:
    """Determine if operator assignment on given value type is in-place, and the method name.

    For example, if operator is '+', return (True, '__iadd__') or (False, '__add__')
    depending on which method is supported by the type.
    """
    typ = get_proper_type(typ)
    method = operators.op_methods[operator]
    if isinstance(typ, Instance):
        if operator in operators.ops_with_inplace_method:
            inplace_method = '__i' + method[2:]
            if typ.type.has_readable_member(inplace_method):
                return True, inplace_method
    return False, method


def is_valid_inferred_type(typ: Type) -> bool:
    """Is an inferred type valid?

    Examples of invalid types include the None type or List[<uninhabited>].

    When not doing strict Optional checking, all types containing None are
    invalid.  When doing strict Optional checking, only None and types that are
    incompletely defined (i.e. contain UninhabitedType) are invalid.
    """
    if isinstance(get_proper_type(typ), (NoneType, UninhabitedType)):
        # With strict Optional checking, we *may* eventually infer NoneType when
        # the initializer is None, but we only do that if we can't infer a
        # specific Optional type.  This resolution happens in
        # leave_partial_types when we pop a partial types scope.
        return False
    return not typ.accept(NothingSeeker())


class NothingSeeker(TypeQuery[bool]):
    """Find any <nothing> types resulting from failed (ambiguous) type inference."""

    def __init__(self) -> None:
        super().__init__(any)

    def visit_uninhabited_type(self, t: UninhabitedType) -> bool:
        return t.ambiguous


class SetNothingToAny(TypeTranslator):
    """Replace all ambiguous <nothing> types with Any (to avoid spurious extra errors)."""

    def visit_uninhabited_type(self, t: UninhabitedType) -> Type:
        if t.ambiguous:
            return AnyType(TypeOfAny.from_error)
        return t

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        # Target of the alias cannot by an ambiguous <nothing>, so we just
        # replace the arguments.
        return t.copy_modified(args=[a.accept(self) for a in t.args])


def is_node_static(node: Optional[Node]) -> Optional[bool]:
    """Find out if a node describes a static function method."""

    if isinstance(node, FuncDef):
        return node.is_static

    if isinstance(node, Var):
        return node.is_staticmethod

    return None


class CheckerScope:
    # We keep two stacks combined, to maintain the relative order
    stack: List[Union[TypeInfo, FuncItem, MypyFile]]

    def __init__(self, module: MypyFile) -> None:
        self.stack = [module]

    def top_function(self) -> Optional[FuncItem]:
        for e in reversed(self.stack):
            if isinstance(e, FuncItem):
                return e
        return None

    def top_non_lambda_function(self) -> Optional[FuncItem]:
        for e in reversed(self.stack):
            if isinstance(e, FuncItem) and not isinstance(e, LambdaExpr):
                return e
        return None

    def active_class(self) -> Optional[TypeInfo]:
        if isinstance(self.stack[-1], TypeInfo):
            return self.stack[-1]
        return None

    def enclosing_class(self) -> Optional[TypeInfo]:
        """Is there a class *directly* enclosing this function?"""
        top = self.top_function()
        assert top, "This method must be called from inside a function"
        index = self.stack.index(top)
        assert index, "CheckerScope stack must always start with a module"
        enclosing = self.stack[index - 1]
        if isinstance(enclosing, TypeInfo):
            return enclosing
        return None

    def active_self_type(self) -> Optional[Union[Instance, TupleType]]:
        """An instance or tuple type representing the current class.

        This returns None unless we are in class body or in a method.
        In particular, inside a function nested in method this returns None.
        """
        info = self.active_class()
        if not info and self.top_function():
            info = self.enclosing_class()
        if info:
            return fill_typevars(info)
        return None

    @contextmanager
    def push_function(self, item: FuncItem) -> Iterator[None]:
        self.stack.append(item)
        yield
        self.stack.pop()

    @contextmanager
    def push_class(self, info: TypeInfo) -> Iterator[None]:
        self.stack.append(info)
        yield
        self.stack.pop()


TKey = TypeVar('TKey')
TValue = TypeVar('TValue')


class DisjointDict(Generic[TKey, TValue]):
    """An variation of the union-find algorithm/data structure where instead of keeping
    track of just disjoint sets, we keep track of disjoint dicts -- keep track of multiple
    Set[Key] -> Set[Value] mappings, where each mapping's keys are guaranteed to be disjoint.

    This data structure is currently used exclusively by 'group_comparison_operands' below
    to merge chains of '==' and 'is' comparisons when two or more chains use the same expression
    in best-case O(n), where n is the number of operands.

    Specifically, the `add_mapping()` function and `items()` functions will take on average
    O(k + v) and O(n) respectively, where k and v are the number of keys and values we're adding
    for a given chain. Note that k <= n and v <= n.

    We hit these average/best-case scenarios for most user code: e.g. when the user has just
    a single chain like 'a == b == c == d == ...' or multiple disjoint chains like
    'a==b < c==d < e==f < ...'. (Note that a naive iterative merging would be O(n^2) for
    the latter case).

    In comparison, this data structure will make 'group_comparison_operands' have a worst-case
    runtime of O(n*log(n)): 'add_mapping()' and 'items()' are worst-case O(k*log(n) + v) and
    O(k*log(n)) respectively. This happens only in the rare case where the user keeps repeatedly
    making disjoint mappings before merging them in a way that persistently dodges the path
    compression optimization in '_lookup_root_id', which would end up constructing a single
    tree of height log_2(n). This makes root lookups no longer amoritized constant time when we
    finally call 'items()'.
    """
    def __init__(self) -> None:
        # Each key maps to a unique ID
        self._key_to_id: Dict[TKey, int] = {}

        # Each id points to the parent id, forming a forest of upwards-pointing trees. If the
        # current id already is the root, it points to itself. We gradually flatten these trees
        # as we perform root lookups: eventually all nodes point directly to its root.
        self._id_to_parent_id: Dict[int, int] = {}

        # Each root id in turn maps to the set of values.
        self._root_id_to_values: Dict[int, Set[TValue]] = {}

    def add_mapping(self, keys: Set[TKey], values: Set[TValue]) -> None:
        """Adds a 'Set[TKey] -> Set[TValue]' mapping. If there already exists a mapping
        containing one or more of the given keys, we merge the input mapping with the old one.

        Note that the given set of keys must be non-empty -- otherwise, nothing happens.
        """
        if len(keys) == 0:
            return

        subtree_roots = [self._lookup_or_make_root_id(key) for key in keys]
        new_root = subtree_roots[0]

        root_values = self._root_id_to_values[new_root]
        root_values.update(values)
        for subtree_root in subtree_roots[1:]:
            if subtree_root == new_root or subtree_root not in self._root_id_to_values:
                continue
            self._id_to_parent_id[subtree_root] = new_root
            root_values.update(self._root_id_to_values.pop(subtree_root))

    def items(self) -> List[Tuple[Set[TKey], Set[TValue]]]:
        """Returns all disjoint mappings in key-value pairs."""
        root_id_to_keys: Dict[int, Set[TKey]] = {}
        for key in self._key_to_id:
            root_id = self._lookup_root_id(key)
            if root_id not in root_id_to_keys:
                root_id_to_keys[root_id] = set()
            root_id_to_keys[root_id].add(key)

        output = []
        for root_id, keys in root_id_to_keys.items():
            output.append((keys, self._root_id_to_values[root_id]))

        return output

    def _lookup_or_make_root_id(self, key: TKey) -> int:
        if key in self._key_to_id:
            return self._lookup_root_id(key)
        else:
            new_id = len(self._key_to_id)
            self._key_to_id[key] = new_id
            self._id_to_parent_id[new_id] = new_id
            self._root_id_to_values[new_id] = set()
            return new_id

    def _lookup_root_id(self, key: TKey) -> int:
        i = self._key_to_id[key]
        while i != self._id_to_parent_id[i]:
            # Optimization: make keys directly point to their grandparents to speed up
            # future traversals. This prevents degenerate trees of height n from forming.
            new_parent = self._id_to_parent_id[self._id_to_parent_id[i]]
            self._id_to_parent_id[i] = new_parent
            i = new_parent
        return i


def group_comparison_operands(pairwise_comparisons: Iterable[Tuple[str, Expression, Expression]],
                              operand_to_literal_hash: Mapping[int, Key],
                              operators_to_group: Set[str],
                              ) -> List[Tuple[str, List[int]]]:
    """Group a series of comparison operands together chained by any operand
    in the 'operators_to_group' set. All other pairwise operands are kept in
    groups of size 2.

    For example, suppose we have the input comparison expression:

        x0 == x1 == x2 < x3 < x4 is x5 is x6 is not x7 is not x8

    If we get these expressions in a pairwise way (e.g. by calling ComparisionExpr's
    'pairwise()' method), we get the following as input:

        [('==', x0, x1), ('==', x1, x2), ('<', x2, x3), ('<', x3, x4),
         ('is', x4, x5), ('is', x5, x6), ('is not', x6, x7), ('is not', x7, x8)]

    If `operators_to_group` is the set {'==', 'is'}, this function will produce
    the following "simplified operator list":

       [("==", [0, 1, 2]), ("<", [2, 3]), ("<", [3, 4]),
        ("is", [4, 5, 6]), ("is not", [6, 7]), ("is not", [7, 8])]

    Note that (a) we yield *indices* to the operands rather then the operand
    expressions themselves and that (b) operands used in a consecutive chain
    of '==' or 'is' are grouped together.

    If two of these chains happen to contain operands with the same underlying
    literal hash (e.g. are assignable and correspond to the same expression),
    we combine those chains together. For example, if we had:

        same == x < y == same

    ...and if 'operand_to_literal_hash' contained the same values for the indices
    0 and 3, we'd produce the following output:

        [("==", [0, 1, 2, 3]), ("<", [1, 2])]

    But if the 'operand_to_literal_hash' did *not* contain an entry, we'd instead
    default to returning:

        [("==", [0, 1]), ("<", [1, 2]), ("==", [2, 3])]

    This function is currently only used to assist with type-narrowing refinements
    and is extracted out to a helper function so we can unit test it.
    """
    groups: Dict[str, DisjointDict[Key, int]] = {op: DisjointDict() for op in operators_to_group}

    simplified_operator_list: List[Tuple[str, List[int]]] = []
    last_operator: Optional[str] = None
    current_indices: Set[int] = set()
    current_hashes: Set[Key] = set()
    for i, (operator, left_expr, right_expr) in enumerate(pairwise_comparisons):
        if last_operator is None:
            last_operator = operator

        if current_indices and (operator != last_operator or operator not in operators_to_group):
            # If some of the operands in the chain are assignable, defer adding it: we might
            # end up needing to merge it with other chains that appear later.
            if len(current_hashes) == 0:
                simplified_operator_list.append((last_operator, sorted(current_indices)))
            else:
                groups[last_operator].add_mapping(current_hashes, current_indices)
            last_operator = operator
            current_indices = set()
            current_hashes = set()

        # Note: 'i' corresponds to the left operand index, so 'i + 1' is the
        # right operand.
        current_indices.add(i)
        current_indices.add(i + 1)

        # We only ever want to combine operands/combine chains for these operators
        if operator in operators_to_group:
            left_hash = operand_to_literal_hash.get(i)
            if left_hash is not None:
                current_hashes.add(left_hash)
            right_hash = operand_to_literal_hash.get(i + 1)
            if right_hash is not None:
                current_hashes.add(right_hash)

    if last_operator is not None:
        if len(current_hashes) == 0:
            simplified_operator_list.append((last_operator, sorted(current_indices)))
        else:
            groups[last_operator].add_mapping(current_hashes, current_indices)

    # Now that we know which chains happen to contain the same underlying expressions
    # and can be merged together, add in this info back to the output.
    for operator, disjoint_dict in groups.items():
        for keys, indices in disjoint_dict.items():
            simplified_operator_list.append((operator, sorted(indices)))

    # For stability, reorder list by the first operand index to appear
    simplified_operator_list.sort(key=lambda item: item[1][0])
    return simplified_operator_list


def is_typed_callable(c: Optional[Type]) -> bool:
    c = get_proper_type(c)
    if not c or not isinstance(c, CallableType):
        return False
    return not all(isinstance(t, AnyType) and t.type_of_any == TypeOfAny.unannotated
                   for t in get_proper_types(c.arg_types + [c.ret_type]))


def is_untyped_decorator(typ: Optional[Type]) -> bool:
    typ = get_proper_type(typ)
    if not typ:
        return True
    elif isinstance(typ, CallableType):
        return not is_typed_callable(typ)
    elif isinstance(typ, Instance):
        method = typ.type.get_method('__call__')
        if method:
            if isinstance(method, Decorator):
                return (
                    is_untyped_decorator(method.func.type)
                    or is_untyped_decorator(method.var.type)
                )

            if isinstance(method.type, Overloaded):
                return any(is_untyped_decorator(item) for item in method.type.items)
            else:
                return not is_typed_callable(method.type)
        else:
            return False
    elif isinstance(typ, Overloaded):
        return any(is_untyped_decorator(item) for item in typ.items)
    return True


def is_static(func: Union[FuncBase, Decorator]) -> bool:
    if isinstance(func, Decorator):
        return is_static(func.func)
    elif isinstance(func, FuncBase):
        return func.is_static
    assert False, "Unexpected func type: {}".format(type(func))


def is_subtype_no_promote(left: Type, right: Type) -> bool:
    return is_subtype(left, right, ignore_promotions=True)


def is_overlapping_types_no_promote(left: Type, right: Type) -> bool:
    return is_overlapping_types(left, right, ignore_promotions=True)


def is_private(node_name: str) -> bool:
    """Check if node is private to class definition."""
    return node_name.startswith('__') and not node_name.endswith('__')


def has_bool_item(typ: ProperType) -> bool:
    """Return True if type is 'bool' or a union with a 'bool' item."""
    if is_named_instance(typ, 'builtins.bool'):
        return True
    if isinstance(typ, UnionType):
        return any(is_named_instance(item, 'builtins.bool')
                   for item in typ.items)
    return False


def collapse_walrus(e: Expression) -> Expression:
    """If an expression is an AssignmentExpr, pull out the assignment target.

    We don't make any attempt to pull out all the targets in code like `x := (y := z)`.
    We could support narrowing those if that sort of code turns out to be common.
    """
    if isinstance(e, AssignmentExpr):
        return e.target
    return e
