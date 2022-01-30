"""Expression type checker. This file is conceptually part of TypeChecker."""

from mypy.backports import OrderedDict, nullcontext
from contextlib import contextmanager
import itertools
from typing import (
    Any, cast, Dict, Set, List, Tuple, Callable, Union, Optional, Sequence, Iterator
)
from typing_extensions import ClassVar, Final, overload, TypeAlias as _TypeAlias

from mypy.errors import report_internal_error
from mypy.typeanal import (
    has_any_from_unimported_type, check_for_explicit_any, set_any_tvars, expand_type_alias,
    make_optional_type,
)
from mypy.semanal_enum import ENUM_BASES
from mypy.types import (
    Type, AnyType, CallableType, Overloaded, NoneType, TypeVarType,
    TupleType, TypedDictType, Instance, ErasedType, UnionType,
    PartialType, DeletedType, UninhabitedType, TypeType, TypeOfAny, LiteralType, LiteralValue,
    is_named_instance, FunctionLike, ParamSpecType, ParamSpecFlavor,
    StarType, is_optional, remove_optional, is_generic_instance, get_proper_type, ProperType,
    get_proper_types, flatten_nested_unions
)
from mypy.nodes import (
    NameExpr, RefExpr, Var, FuncDef, OverloadedFuncDef, TypeInfo, CallExpr,
    MemberExpr, IntExpr, StrExpr, BytesExpr, UnicodeExpr, FloatExpr,
    OpExpr, UnaryExpr, IndexExpr, CastExpr, RevealExpr, TypeApplication, ListExpr,
    TupleExpr, DictExpr, LambdaExpr, SuperExpr, SliceExpr, Context, Expression,
    ListComprehension, GeneratorExpr, SetExpr, MypyFile, Decorator,
    ConditionalExpr, ComparisonExpr, TempNode, SetComprehension, AssignmentExpr,
    DictionaryComprehension, ComplexExpr, EllipsisExpr, StarExpr, AwaitExpr, YieldExpr,
    YieldFromExpr, TypedDictExpr, PromoteExpr, NewTypeExpr, NamedTupleExpr, TypeVarExpr,
    TypeAliasExpr, BackquoteExpr, EnumCallExpr, TypeAlias, SymbolNode, PlaceholderNode,
    ParamSpecExpr,
    ArgKind, ARG_POS, ARG_NAMED, ARG_STAR, ARG_STAR2, LITERAL_TYPE, REVEAL_TYPE,
)
from mypy.literals import literal
from mypy import nodes
from mypy import operators
import mypy.checker
from mypy import types
from mypy.sametypes import is_same_type
from mypy.erasetype import replace_meta_vars, erase_type, remove_instance_last_known_values
from mypy.maptype import map_instance_to_supertype
from mypy.messages import MessageBuilder
from mypy import message_registry
from mypy.infer import (
    ArgumentInferContext, infer_type_arguments, infer_function_type_arguments,
)
from mypy import join
from mypy.meet import narrow_declared_type, is_overlapping_types
from mypy.subtypes import is_subtype, is_proper_subtype, is_equivalent, non_method_protocol_members
from mypy import applytype
from mypy import erasetype
from mypy.checkmember import analyze_member_access, type_object_type
from mypy.argmap import ArgTypeExpander, map_actuals_to_formals, map_formals_to_actuals
from mypy.checkstrformat import StringFormatterChecker
from mypy.expandtype import expand_type, expand_type_by_instance, freshen_function_type_vars
from mypy.util import split_module_names
from mypy.typevars import fill_typevars
from mypy.visitor import ExpressionVisitor
from mypy.plugin import (
    Plugin,
    MethodContext, MethodSigContext,
    FunctionContext, FunctionSigContext,
)
from mypy.typeops import (
    try_expanding_sum_type_to_union, tuple_fallback, make_simplified_union,
    true_only, false_only, erase_to_union_or_bound, function_type,
    callable_type, try_getting_str_literals, custom_special_method,
    is_literal_type_like,
)
from mypy.message_registry import ErrorMessage
import mypy.errorcodes as codes

# Type of callback user for checking individual function arguments. See
# check_args() below for details.
ArgChecker: _TypeAlias = Callable[[
        Type,
        Type,
        ArgKind,
        Type,
        int,
        int,
        CallableType,
        Optional[Type],
        Context,
        Context,
        MessageBuilder,
    ],
    None,
]

# Maximum nesting level for math union in overloads, setting this to large values
# may cause performance issues. The reason is that although union math algorithm we use
# nicely captures most corner cases, its worst case complexity is exponential,
# see https://github.com/python/mypy/pull/5255#discussion_r196896335 for discussion.
MAX_UNIONS: Final = 5


# Types considered safe for comparisons with --strict-equality due to known behaviour of __eq__.
# NOTE: All these types are subtypes of AbstractSet.
OVERLAPPING_TYPES_ALLOWLIST: Final = [
    "builtins.set",
    "builtins.frozenset",
    "typing.KeysView",
    "typing.ItemsView",
    "builtins._dict_keys",
    "builtins._dict_items",
    "_collections_abc.dict_keys",
    "_collections_abc.dict_items",
]


class TooManyUnions(Exception):
    """Indicates that we need to stop splitting unions in an attempt
    to match an overload in order to save performance.
    """


def extract_refexpr_names(expr: RefExpr) -> Set[str]:
    """Recursively extracts all module references from a reference expression.

    Note that currently, the only two subclasses of RefExpr are NameExpr and
    MemberExpr."""
    output: Set[str] = set()
    while isinstance(expr.node, MypyFile) or expr.fullname is not None:
        if isinstance(expr.node, MypyFile) and expr.fullname is not None:
            # If it's None, something's wrong (perhaps due to an
            # import cycle or a suppressed error).  For now we just
            # skip it.
            output.add(expr.fullname)

        if isinstance(expr, NameExpr):
            is_suppressed_import = isinstance(expr.node, Var) and expr.node.is_suppressed_import
            if isinstance(expr.node, TypeInfo):
                # Reference to a class or a nested class
                output.update(split_module_names(expr.node.module_name))
            elif expr.fullname is not None and '.' in expr.fullname and not is_suppressed_import:
                # Everything else (that is not a silenced import within a class)
                output.add(expr.fullname.rsplit('.', 1)[0])
            break
        elif isinstance(expr, MemberExpr):
            if isinstance(expr.expr, RefExpr):
                expr = expr.expr
            else:
                break
        else:
            raise AssertionError("Unknown RefExpr subclass: {}".format(type(expr)))
    return output


class Finished(Exception):
    """Raised if we can terminate overload argument check early (no match)."""


class ExpressionChecker(ExpressionVisitor[Type]):
    """Expression type checker.

    This class works closely together with checker.TypeChecker.
    """

    # Some services are provided by a TypeChecker instance.
    chk: "mypy.checker.TypeChecker"
    # This is shared with TypeChecker, but stored also here for convenience.
    msg: MessageBuilder
    # Type context for type inference
    type_context: List[Optional[Type]]

    strfrm_checker: StringFormatterChecker
    plugin: Plugin

    def __init__(self,
                 chk: 'mypy.checker.TypeChecker',
                 msg: MessageBuilder,
                 plugin: Plugin) -> None:
        """Construct an expression type checker."""
        self.chk = chk
        self.msg = msg
        self.plugin = plugin
        self.type_context = [None]

        # Temporary overrides for expression types. This is currently
        # used by the union math in overloads.
        # TODO: refactor this to use a pattern similar to one in
        # multiassign_from_union, or maybe even combine the two?
        self.type_overrides: Dict[Expression, Type] = {}
        self.strfrm_checker = StringFormatterChecker(self, self.chk, self.msg)

    def visit_name_expr(self, e: NameExpr) -> Type:
        """Type check a name expression.

        It can be of any kind: local, member or global.
        """
        self.chk.module_refs.update(extract_refexpr_names(e))
        result = self.analyze_ref_expr(e)
        return self.narrow_type_from_binder(e, result)

    def analyze_ref_expr(self, e: RefExpr, lvalue: bool = False) -> Type:
        result: Optional[Type] = None
        node = e.node

        if isinstance(e, NameExpr) and e.is_special_form:
            # A special form definition, nothing to check here.
            return AnyType(TypeOfAny.special_form)

        if isinstance(node, Var):
            # Variable reference.
            result = self.analyze_var_ref(node, e)
            if isinstance(result, PartialType):
                result = self.chk.handle_partial_var_type(result, lvalue, node, e)
        elif isinstance(node, FuncDef):
            # Reference to a global function.
            result = function_type(node, self.named_type('builtins.function'))
        elif isinstance(node, OverloadedFuncDef) and node.type is not None:
            # node.type is None when there are multiple definitions of a function
            # and it's decorated by something that is not typing.overload
            # TODO: use a dummy Overloaded instead of AnyType in this case
            # like we do in mypy.types.function_type()?
            result = node.type
        elif isinstance(node, TypeInfo):
            # Reference to a type object.
            result = type_object_type(node, self.named_type)
            if (isinstance(result, CallableType) and
                    isinstance(result.ret_type, Instance)):  # type: ignore
                # We need to set correct line and column
                # TODO: always do this in type_object_type by passing the original context
                result.ret_type.line = e.line
                result.ret_type.column = e.column
            if isinstance(get_proper_type(self.type_context[-1]), TypeType):
                # This is the type in a Type[] expression, so substitute type
                # variables with Any.
                result = erasetype.erase_typevars(result)
        elif isinstance(node, MypyFile):
            # Reference to a module object.
            try:
                result = self.named_type('types.ModuleType')
            except KeyError:
                # In test cases might 'types' may not be available.
                # Fall back to a dummy 'object' type instead to
                # avoid a crash.
                result = self.named_type('builtins.object')
        elif isinstance(node, Decorator):
            result = self.analyze_var_ref(node.var, e)
        elif isinstance(node, TypeAlias):
            # Something that refers to a type alias appears in runtime context.
            # Note that we suppress bogus errors for alias redefinitions,
            # they are already reported in semanal.py.
            result = self.alias_type_in_runtime_context(node, node.no_args, e,
                                                        alias_definition=e.is_alias_rvalue
                                                        or lvalue)
        elif isinstance(node, (TypeVarExpr, ParamSpecExpr)):
            result = self.object_type()
        else:
            if isinstance(node, PlaceholderNode):
                assert False, 'PlaceholderNode %r leaked to checker' % node.fullname
            # Unknown reference; use any type implicitly to avoid
            # generating extra type errors.
            result = AnyType(TypeOfAny.from_error)
        assert result is not None
        return result

    def analyze_var_ref(self, var: Var, context: Context) -> Type:
        if var.type:
            var_type = get_proper_type(var.type)
            if isinstance(var_type, Instance):
                if self.is_literal_context() and var_type.last_known_value is not None:
                    return var_type.last_known_value
                if var.name in {'True', 'False'}:
                    return self.infer_literal_expr_type(var.name == 'True', 'builtins.bool')
            return var.type
        else:
            if not var.is_ready and self.chk.in_checked_function():
                self.chk.handle_cannot_determine_type(var.name, context)
            # Implicit 'Any' type.
            return AnyType(TypeOfAny.special_form)

    def visit_call_expr(self, e: CallExpr, allow_none_return: bool = False) -> Type:
        """Type check a call expression."""
        if e.analyzed:
            if isinstance(e.analyzed, NamedTupleExpr) and not e.analyzed.is_typed:
                # Type check the arguments, but ignore the results. This relies
                # on the typeshed stubs to type check the arguments.
                self.visit_call_expr_inner(e)
            # It's really a special form that only looks like a call.
            return self.accept(e.analyzed, self.type_context[-1])
        return self.visit_call_expr_inner(e, allow_none_return=allow_none_return)

    def visit_call_expr_inner(self, e: CallExpr, allow_none_return: bool = False) -> Type:
        if isinstance(e.callee, RefExpr) and isinstance(e.callee.node, TypeInfo) and \
                e.callee.node.typeddict_type is not None:
            # Use named fallback for better error messages.
            typeddict_type = e.callee.node.typeddict_type.copy_modified(
                fallback=Instance(e.callee.node, []))
            return self.check_typeddict_call(typeddict_type, e.arg_kinds, e.arg_names, e.args, e)
        if (isinstance(e.callee, NameExpr) and e.callee.name in ('isinstance', 'issubclass')
                and len(e.args) == 2):
            for typ in mypy.checker.flatten(e.args[1]):
                node = None
                if isinstance(typ, NameExpr):
                    try:
                        node = self.chk.lookup_qualified(typ.name)
                    except KeyError:
                        # Undefined names should already be reported in semantic analysis.
                        pass
                if is_expr_literal_type(typ):
                    self.msg.cannot_use_function_with_type(e.callee.name, "Literal", e)
                    continue
                if (node and isinstance(node.node, TypeAlias)
                        and isinstance(get_proper_type(node.node.target), AnyType)):
                    self.msg.cannot_use_function_with_type(e.callee.name, "Any", e)
                    continue
                if ((isinstance(typ, IndexExpr)
                        and isinstance(typ.analyzed, (TypeApplication, TypeAliasExpr)))
                        or (isinstance(typ, NameExpr) and node and
                            isinstance(node.node, TypeAlias) and not node.node.no_args)):
                    self.msg.type_arguments_not_allowed(e)
                if isinstance(typ, RefExpr) and isinstance(typ.node, TypeInfo):
                    if typ.node.typeddict_type:
                        self.msg.cannot_use_function_with_type(e.callee.name, "TypedDict", e)
                    elif typ.node.is_newtype:
                        self.msg.cannot_use_function_with_type(e.callee.name, "NewType", e)
        self.try_infer_partial_type(e)
        type_context = None
        if isinstance(e.callee, LambdaExpr):
            formal_to_actual = map_actuals_to_formals(
                e.arg_kinds, e.arg_names,
                e.callee.arg_kinds, e.callee.arg_names,
                lambda i: self.accept(e.args[i]))

            arg_types = [join.join_type_list([self.accept(e.args[j]) for j in formal_to_actual[i]])
                         for i in range(len(e.callee.arg_kinds))]
            type_context = CallableType(arg_types, e.callee.arg_kinds, e.callee.arg_names,
                                        ret_type=self.object_type(),
                                        fallback=self.named_type('builtins.function'))
        callee_type = get_proper_type(self.accept(e.callee, type_context, always_allow_any=True))
        if (isinstance(e.callee, RefExpr)
                and isinstance(callee_type, CallableType)
                and callee_type.type_guard is not None):
            # Cache it for find_isinstance_check()
            e.callee.type_guard = callee_type.type_guard
        if (self.chk.options.disallow_untyped_calls and
                self.chk.in_checked_function() and
                isinstance(callee_type, CallableType)
                and callee_type.implicit):
            self.msg.untyped_function_call(callee_type, e)

        # Figure out the full name of the callee for plugin lookup.
        object_type = None
        member = None
        fullname = None
        if isinstance(e.callee, RefExpr):
            # There are two special cases where plugins might act:
            # * A "static" reference/alias to a class or function;
            #   get_function_hook() will be invoked for these.
            fullname = e.callee.fullname
            if isinstance(e.callee.node, TypeAlias):
                target = get_proper_type(e.callee.node.target)
                if isinstance(target, Instance):
                    fullname = target.type.fullname
            # * Call to a method on object that has a full name (see
            #   method_fullname() for details on supported objects);
            #   get_method_hook() and get_method_signature_hook() will
            #   be invoked for these.
            if (fullname is None
                    and isinstance(e.callee, MemberExpr)
                    and e.callee.expr in self.chk.type_map):
                member = e.callee.name
                object_type = self.chk.type_map[e.callee.expr]
        ret_type = self.check_call_expr_with_callee_type(callee_type, e, fullname,
                                                         object_type, member)
        if isinstance(e.callee, RefExpr) and len(e.args) == 2:
            if e.callee.fullname in ('builtins.isinstance', 'builtins.issubclass'):
                self.check_runtime_protocol_test(e)
            if e.callee.fullname == 'builtins.issubclass':
                self.check_protocol_issubclass(e)
        if isinstance(e.callee, MemberExpr) and e.callee.name == 'format':
            self.check_str_format_call(e)
        ret_type = get_proper_type(ret_type)
        if isinstance(ret_type, UninhabitedType) and not ret_type.ambiguous:
            self.chk.binder.unreachable()
        # Warn on calls to functions that always return None. The check
        # of ret_type is both a common-case optimization and prevents reporting
        # the error in dynamic functions (where it will be Any).
        if (not allow_none_return and isinstance(ret_type, NoneType)
                and self.always_returns_none(e.callee)):
            self.chk.msg.does_not_return_value(callee_type, e)
            return AnyType(TypeOfAny.from_error)
        return ret_type

    def check_str_format_call(self, e: CallExpr) -> None:
        """More precise type checking for str.format() calls on literals."""
        assert isinstance(e.callee, MemberExpr)
        format_value = None
        if isinstance(e.callee.expr, (StrExpr, UnicodeExpr)):
            format_value = e.callee.expr.value
        elif e.callee.expr in self.chk.type_map:
            base_typ = try_getting_literal(self.chk.type_map[e.callee.expr])
            if isinstance(base_typ, LiteralType) and isinstance(base_typ.value, str):
                format_value = base_typ.value
        if format_value is not None:
            self.strfrm_checker.check_str_format_call(e, format_value)

    def method_fullname(self, object_type: Type, method_name: str) -> Optional[str]:
        """Convert a method name to a fully qualified name, based on the type of the object that
        it is invoked on. Return `None` if the name of `object_type` cannot be determined.
        """
        object_type = get_proper_type(object_type)

        if isinstance(object_type, CallableType) and object_type.is_type_obj():
            # For class method calls, object_type is a callable representing the class object.
            # We "unwrap" it to a regular type, as the class/instance method difference doesn't
            # affect the fully qualified name.
            object_type = get_proper_type(object_type.ret_type)
        elif isinstance(object_type, TypeType):
            object_type = object_type.item

        type_name = None
        if isinstance(object_type, Instance):
            type_name = object_type.type.fullname
        elif isinstance(object_type, (TypedDictType, LiteralType)):
            info = object_type.fallback.type.get_containing_type_info(method_name)
            type_name = info.fullname if info is not None else None
        elif isinstance(object_type, TupleType):
            type_name = tuple_fallback(object_type).type.fullname

        if type_name is not None:
            return '{}.{}'.format(type_name, method_name)
        else:
            return None

    def always_returns_none(self, node: Expression) -> bool:
        """Check if `node` refers to something explicitly annotated as only returning None."""
        if isinstance(node, RefExpr):
            if self.defn_returns_none(node.node):
                return True
        if isinstance(node, MemberExpr) and node.node is None:  # instance or class attribute
            typ = get_proper_type(self.chk.type_map.get(node.expr))
            if isinstance(typ, Instance):
                info = typ.type
            elif isinstance(typ, CallableType) and typ.is_type_obj():
                ret_type = get_proper_type(typ.ret_type)
                if isinstance(ret_type, Instance):
                    info = ret_type.type
                else:
                    return False
            else:
                return False
            sym = info.get(node.name)
            if sym and self.defn_returns_none(sym.node):
                return True
        return False

    def defn_returns_none(self, defn: Optional[SymbolNode]) -> bool:
        """Check if `defn` can _only_ return None."""
        if isinstance(defn, FuncDef):
            return (isinstance(defn.type, CallableType) and
                    isinstance(get_proper_type(defn.type.ret_type), NoneType))
        if isinstance(defn, OverloadedFuncDef):
            return all(self.defn_returns_none(item) for item in defn.items)
        if isinstance(defn, Var):
            typ = get_proper_type(defn.type)
            if (not defn.is_inferred and isinstance(typ, CallableType) and
                    isinstance(get_proper_type(typ.ret_type), NoneType)):
                return True
            if isinstance(typ, Instance):
                sym = typ.type.get('__call__')
                if sym and self.defn_returns_none(sym.node):
                    return True
        return False

    def check_runtime_protocol_test(self, e: CallExpr) -> None:
        for expr in mypy.checker.flatten(e.args[1]):
            tp = get_proper_type(self.chk.type_map[expr])
            if (isinstance(tp, CallableType) and tp.is_type_obj() and
                    tp.type_object().is_protocol and
                    not tp.type_object().runtime_protocol):
                self.chk.fail(message_registry.RUNTIME_PROTOCOL_EXPECTED, e)

    def check_protocol_issubclass(self, e: CallExpr) -> None:
        for expr in mypy.checker.flatten(e.args[1]):
            tp = get_proper_type(self.chk.type_map[expr])
            if (isinstance(tp, CallableType) and tp.is_type_obj() and
                    tp.type_object().is_protocol):
                attr_members = non_method_protocol_members(tp.type_object())
                if attr_members:
                    self.chk.msg.report_non_method_protocol(tp.type_object(),
                                                            attr_members, e)

    def check_typeddict_call(self, callee: TypedDictType,
                             arg_kinds: List[ArgKind],
                             arg_names: Sequence[Optional[str]],
                             args: List[Expression],
                             context: Context) -> Type:
        if len(args) >= 1 and all([ak == ARG_NAMED for ak in arg_kinds]):
            # ex: Point(x=42, y=1337)
            assert all(arg_name is not None for arg_name in arg_names)
            item_names = cast(List[str], arg_names)
            item_args = args
            return self.check_typeddict_call_with_kwargs(
                callee, OrderedDict(zip(item_names, item_args)), context)

        if len(args) == 1 and arg_kinds[0] == ARG_POS:
            unique_arg = args[0]
            if isinstance(unique_arg, DictExpr):
                # ex: Point({'x': 42, 'y': 1337})
                return self.check_typeddict_call_with_dict(callee, unique_arg, context)
            if isinstance(unique_arg, CallExpr) and isinstance(unique_arg.analyzed, DictExpr):
                # ex: Point(dict(x=42, y=1337))
                return self.check_typeddict_call_with_dict(callee, unique_arg.analyzed, context)

        if len(args) == 0:
            # ex: EmptyDict()
            return self.check_typeddict_call_with_kwargs(
                callee, OrderedDict(), context)

        self.chk.fail(message_registry.INVALID_TYPEDDICT_ARGS, context)
        return AnyType(TypeOfAny.from_error)

    def validate_typeddict_kwargs(
            self, kwargs: DictExpr) -> 'Optional[OrderedDict[str, Expression]]':
        item_args = [item[1] for item in kwargs.items]

        item_names = []  # List[str]
        for item_name_expr, item_arg in kwargs.items:
            literal_value = None
            if item_name_expr:
                key_type = self.accept(item_name_expr)
                values = try_getting_str_literals(item_name_expr, key_type)
                if values and len(values) == 1:
                    literal_value = values[0]
            if literal_value is None:
                key_context = item_name_expr or item_arg
                self.chk.fail(message_registry.TYPEDDICT_KEY_MUST_BE_STRING_LITERAL,
                              key_context)
                return None
            else:
                item_names.append(literal_value)
        return OrderedDict(zip(item_names, item_args))

    def match_typeddict_call_with_dict(self, callee: TypedDictType,
                                       kwargs: DictExpr,
                                       context: Context) -> bool:
        validated_kwargs = self.validate_typeddict_kwargs(kwargs=kwargs)
        if validated_kwargs is not None:
            return (callee.required_keys <= set(validated_kwargs.keys())
                <= set(callee.items.keys()))
        else:
            return False

    def check_typeddict_call_with_dict(self, callee: TypedDictType,
                                       kwargs: DictExpr,
                                       context: Context) -> Type:
        validated_kwargs = self.validate_typeddict_kwargs(kwargs=kwargs)
        if validated_kwargs is not None:
            return self.check_typeddict_call_with_kwargs(
                callee,
                kwargs=validated_kwargs,
                context=context)
        else:
            return AnyType(TypeOfAny.from_error)

    def check_typeddict_call_with_kwargs(self, callee: TypedDictType,
                                         kwargs: 'OrderedDict[str, Expression]',
                                         context: Context) -> Type:
        if not (callee.required_keys <= set(kwargs.keys()) <= set(callee.items.keys())):
            expected_keys = [key for key in callee.items.keys()
                             if key in callee.required_keys or key in kwargs.keys()]
            actual_keys = kwargs.keys()
            self.msg.unexpected_typeddict_keys(
                callee,
                expected_keys=expected_keys,
                actual_keys=list(actual_keys),
                context=context)
            return AnyType(TypeOfAny.from_error)

        for (item_name, item_expected_type) in callee.items.items():
            if item_name in kwargs:
                item_value = kwargs[item_name]
                self.chk.check_simple_assignment(
                    lvalue_type=item_expected_type, rvalue=item_value, context=item_value,
                    msg=message_registry.INCOMPATIBLE_TYPES,
                    lvalue_name='TypedDict item "{}"'.format(item_name),
                    rvalue_name='expression',
                    code=codes.TYPEDDICT_ITEM)

        return callee

    def get_partial_self_var(self, expr: MemberExpr) -> Optional[Var]:
        """Get variable node for a partial self attribute.

        If the expression is not a self attribute, or attribute is not variable,
        or variable is not partial, return None.
        """
        if not (isinstance(expr.expr, NameExpr) and
                isinstance(expr.expr.node, Var) and expr.expr.node.is_self):
            # Not a self.attr expression.
            return None
        info = self.chk.scope.enclosing_class()
        if not info or expr.name not in info.names:
            # Don't mess with partial types in superclasses.
            return None
        sym = info.names[expr.name]
        if isinstance(sym.node, Var) and isinstance(sym.node.type, PartialType):
            return sym.node
        return None

    # Types and methods that can be used to infer partial types.
    item_args: ClassVar[Dict[str, List[str]]] = {
        "builtins.list": ["append"],
        "builtins.set": ["add", "discard"],
    }
    container_args: ClassVar[Dict[str, Dict[str, List[str]]]] = {
        "builtins.list": {"extend": ["builtins.list"]},
        "builtins.dict": {"update": ["builtins.dict"]},
        "collections.OrderedDict": {"update": ["builtins.dict"]},
        "builtins.set": {"update": ["builtins.set", "builtins.list"]},
    }

    def try_infer_partial_type(self, e: CallExpr) -> None:
        """Try to make partial type precise from a call."""
        if not isinstance(e.callee, MemberExpr):
            return
        callee = e.callee
        if isinstance(callee.expr, RefExpr):
            # Call a method with a RefExpr callee, such as 'x.method(...)'.
            ret = self.get_partial_var(callee.expr)
            if ret is None:
                return
            var, partial_types = ret
            typ = self.try_infer_partial_value_type_from_call(e, callee.name, var)
            if typ is not None:
                var.type = typ
                del partial_types[var]
        elif isinstance(callee.expr, IndexExpr) and isinstance(callee.expr.base, RefExpr):
            # Call 'x[y].method(...)'; may infer type of 'x' if it's a partial defaultdict.
            if callee.expr.analyzed is not None:
                return  # A special form
            base = callee.expr.base
            index = callee.expr.index
            ret = self.get_partial_var(base)
            if ret is None:
                return
            var, partial_types = ret
            partial_type = get_partial_instance_type(var.type)
            if partial_type is None or partial_type.value_type is None:
                return
            value_type = self.try_infer_partial_value_type_from_call(e, callee.name, var)
            if value_type is not None:
                # Infer key type.
                key_type = self.accept(index)
                if mypy.checker.is_valid_inferred_type(key_type):
                    # Store inferred partial type.
                    assert partial_type.type is not None
                    typename = partial_type.type.fullname
                    var.type = self.chk.named_generic_type(typename,
                                                           [key_type, value_type])
                    del partial_types[var]

    def get_partial_var(self, ref: RefExpr) -> Optional[Tuple[Var, Dict[Var, Context]]]:
        var = ref.node
        if var is None and isinstance(ref, MemberExpr):
            var = self.get_partial_self_var(ref)
        if not isinstance(var, Var):
            return None
        partial_types = self.chk.find_partial_types(var)
        if partial_types is None:
            return None
        return var, partial_types

    def try_infer_partial_value_type_from_call(
            self,
            e: CallExpr,
            methodname: str,
            var: Var) -> Optional[Instance]:
        """Try to make partial type precise from a call such as 'x.append(y)'."""
        if self.chk.current_node_deferred:
            return None
        partial_type = get_partial_instance_type(var.type)
        if partial_type is None:
            return None
        if partial_type.value_type:
            typename = partial_type.value_type.type.fullname
        else:
            assert partial_type.type is not None
            typename = partial_type.type.fullname
        # Sometimes we can infer a full type for a partial List, Dict or Set type.
        # TODO: Don't infer argument expression twice.
        if (typename in self.item_args and methodname in self.item_args[typename]
                and e.arg_kinds == [ARG_POS]):
            item_type = self.accept(e.args[0])
            if mypy.checker.is_valid_inferred_type(item_type):
                return self.chk.named_generic_type(typename, [item_type])
        elif (typename in self.container_args
              and methodname in self.container_args[typename]
              and e.arg_kinds == [ARG_POS]):
            arg_type = get_proper_type(self.accept(e.args[0]))
            if isinstance(arg_type, Instance):
                arg_typename = arg_type.type.fullname
                if arg_typename in self.container_args[typename][methodname]:
                    if all(mypy.checker.is_valid_inferred_type(item_type)
                           for item_type in arg_type.args):
                        return self.chk.named_generic_type(typename,
                                                           list(arg_type.args))
            elif isinstance(arg_type, AnyType):
                return self.chk.named_type(typename)

        return None

    def apply_function_plugin(self,
                              callee: CallableType,
                              arg_kinds: List[ArgKind],
                              arg_types: List[Type],
                              arg_names: Optional[Sequence[Optional[str]]],
                              formal_to_actual: List[List[int]],
                              args: List[Expression],
                              fullname: str,
                              object_type: Optional[Type],
                              context: Context) -> Type:
        """Use special case logic to infer the return type of a specific named function/method.

        Caller must ensure that a plugin hook exists. There are two different cases:

        - If object_type is None, the caller must ensure that a function hook exists
          for fullname.
        - If object_type is not None, the caller must ensure that a method hook exists
          for fullname.

        Return the inferred return type.
        """
        num_formals = len(callee.arg_types)
        formal_arg_types: List[List[Type]] = [[] for _ in range(num_formals)]
        formal_arg_exprs: List[List[Expression]] = [[] for _ in range(num_formals)]
        formal_arg_names: List[List[Optional[str]]] = [[] for _ in range(num_formals)]
        formal_arg_kinds: List[List[ArgKind]] = [[] for _ in range(num_formals)]
        for formal, actuals in enumerate(formal_to_actual):
            for actual in actuals:
                formal_arg_types[formal].append(arg_types[actual])
                formal_arg_exprs[formal].append(args[actual])
                if arg_names:
                    formal_arg_names[formal].append(arg_names[actual])
                formal_arg_kinds[formal].append(arg_kinds[actual])

        if object_type is None:
            # Apply function plugin
            callback = self.plugin.get_function_hook(fullname)
            assert callback is not None  # Assume that caller ensures this
            return callback(
                FunctionContext(formal_arg_types, formal_arg_kinds,
                                callee.arg_names, formal_arg_names,
                                callee.ret_type, formal_arg_exprs, context, self.chk))
        else:
            # Apply method plugin
            method_callback = self.plugin.get_method_hook(fullname)
            assert method_callback is not None  # Assume that caller ensures this
            object_type = get_proper_type(object_type)
            return method_callback(
                MethodContext(object_type, formal_arg_types, formal_arg_kinds,
                              callee.arg_names, formal_arg_names,
                              callee.ret_type, formal_arg_exprs, context, self.chk))

    def apply_signature_hook(
            self, callee: FunctionLike, args: List[Expression],
            arg_kinds: List[ArgKind],
            arg_names: Optional[Sequence[Optional[str]]],
            hook: Callable[
                [List[List[Expression]], CallableType],
                FunctionLike,
            ]) -> FunctionLike:
        """Helper to apply a signature hook for either a function or method"""
        if isinstance(callee, CallableType):
            num_formals = len(callee.arg_kinds)
            formal_to_actual = map_actuals_to_formals(
                arg_kinds, arg_names,
                callee.arg_kinds, callee.arg_names,
                lambda i: self.accept(args[i]))
            formal_arg_exprs: List[List[Expression]] = [[] for _ in range(num_formals)]
            for formal, actuals in enumerate(formal_to_actual):
                for actual in actuals:
                    formal_arg_exprs[formal].append(args[actual])
            return hook(formal_arg_exprs, callee)
        else:
            assert isinstance(callee, Overloaded)
            items = []
            for item in callee.items:
                adjusted = self.apply_signature_hook(
                    item, args, arg_kinds, arg_names, hook)
                assert isinstance(adjusted, CallableType)
                items.append(adjusted)
            return Overloaded(items)

    def apply_function_signature_hook(
            self, callee: FunctionLike, args: List[Expression],
            arg_kinds: List[ArgKind], context: Context,
            arg_names: Optional[Sequence[Optional[str]]],
            signature_hook: Callable[[FunctionSigContext], FunctionLike]) -> FunctionLike:
        """Apply a plugin hook that may infer a more precise signature for a function."""
        return self.apply_signature_hook(
            callee, args, arg_kinds, arg_names,
            (lambda args, sig:
             signature_hook(FunctionSigContext(args, sig, context, self.chk))))

    def apply_method_signature_hook(
            self, callee: FunctionLike, args: List[Expression],
            arg_kinds: List[ArgKind], context: Context,
            arg_names: Optional[Sequence[Optional[str]]], object_type: Type,
            signature_hook: Callable[[MethodSigContext], FunctionLike]) -> FunctionLike:
        """Apply a plugin hook that may infer a more precise signature for a method."""
        pobject_type = get_proper_type(object_type)
        return self.apply_signature_hook(
            callee, args, arg_kinds, arg_names,
            (lambda args, sig:
             signature_hook(MethodSigContext(pobject_type, args, sig, context, self.chk))))

    def transform_callee_type(
            self, callable_name: Optional[str], callee: Type, args: List[Expression],
            arg_kinds: List[ArgKind], context: Context,
            arg_names: Optional[Sequence[Optional[str]]] = None,
            object_type: Optional[Type] = None) -> Type:
        """Attempt to determine a more accurate signature for a method call.

        This is done by looking up and applying a method signature hook (if one exists for the
        given method name).

        If no matching method signature hook is found, callee is returned unmodified. The same
        happens if the arguments refer to a non-method callable (this is allowed so that the code
        calling transform_callee_type needs to perform fewer boilerplate checks).

        Note: this method is *not* called automatically as part of check_call, because in some
        cases check_call is called multiple times while checking a single call (for example when
        dealing with overloads). Instead, this method needs to be called explicitly
        (if appropriate) before the signature is passed to check_call.
        """
        callee = get_proper_type(callee)
        if callable_name is not None and isinstance(callee, FunctionLike):
            if object_type is not None:
                method_sig_hook = self.plugin.get_method_signature_hook(callable_name)
                if method_sig_hook:
                    return self.apply_method_signature_hook(
                        callee, args, arg_kinds, context, arg_names, object_type, method_sig_hook)
            else:
                function_sig_hook = self.plugin.get_function_signature_hook(callable_name)
                if function_sig_hook:
                    return self.apply_function_signature_hook(
                        callee, args, arg_kinds, context, arg_names, function_sig_hook)

        return callee

    def check_call_expr_with_callee_type(self,
                                         callee_type: Type,
                                         e: CallExpr,
                                         callable_name: Optional[str],
                                         object_type: Optional[Type],
                                         member: Optional[str] = None) -> Type:
        """Type check call expression.

        The callee_type should be used as the type of callee expression. In particular,
        in case of a union type this can be a particular item of the union, so that we can
        apply plugin hooks to each item.

        The 'member', 'callable_name' and 'object_type' are only used to call plugin hooks.
        If 'callable_name' is None but 'member' is not None (member call), try constructing
        'callable_name' using 'object_type' (the base type on which the method is called),
        for example 'typing.Mapping.get'.
        """
        if callable_name is None and member is not None:
            assert object_type is not None
            callable_name = self.method_fullname(object_type, member)
        object_type = get_proper_type(object_type)
        if callable_name:
            # Try to refine the call signature using plugin hooks before checking the call.
            callee_type = self.transform_callee_type(
                callable_name, callee_type, e.args, e.arg_kinds, e, e.arg_names, object_type)
        # Unions are special-cased to allow plugins to act on each item in the union.
        elif member is not None and isinstance(object_type, UnionType):
            return self.check_union_call_expr(e, object_type, member)
        return self.check_call(callee_type, e.args, e.arg_kinds, e,
                               e.arg_names, callable_node=e.callee,
                               callable_name=callable_name,
                               object_type=object_type)[0]

    def check_union_call_expr(self, e: CallExpr, object_type: UnionType, member: str) -> Type:
        """"Type check calling a member expression where the base type is a union."""
        res: List[Type] = []
        for typ in object_type.relevant_items():
            # Member access errors are already reported when visiting the member expression.
            with self.msg.disable_errors():
                item = analyze_member_access(member, typ, e, False, False, False,
                                             self.msg, original_type=object_type, chk=self.chk,
                                             in_literal_context=self.is_literal_context(),
                                             self_type=typ)
            narrowed = self.narrow_type_from_binder(e.callee, item, skip_non_overlapping=True)
            if narrowed is None:
                continue
            callable_name = self.method_fullname(typ, member)
            item_object_type = typ if callable_name else None
            res.append(self.check_call_expr_with_callee_type(narrowed, e, callable_name,
                                                             item_object_type))
        return make_simplified_union(res)

    def check_call(self,
                   callee: Type,
                   args: List[Expression],
                   arg_kinds: List[ArgKind],
                   context: Context,
                   arg_names: Optional[Sequence[Optional[str]]] = None,
                   callable_node: Optional[Expression] = None,
                   arg_messages: Optional[MessageBuilder] = None,
                   callable_name: Optional[str] = None,
                   object_type: Optional[Type] = None) -> Tuple[Type, Type]:
        """Type check a call.

        Also infer type arguments if the callee is a generic function.

        Return (result type, inferred callee type).

        Arguments:
            callee: type of the called value
            args: actual argument expressions
            arg_kinds: contains nodes.ARG_* constant for each argument in args
                 describing whether the argument is positional, *arg, etc.
            context: current expression context, used for inference.
            arg_names: names of arguments (optional)
            callable_node: associate the inferred callable type to this node,
                if specified
            arg_messages: utility for generating messages, can be swapped to suppress errors,
                by default uses 'self.msg' to show errors
            callable_name: Fully-qualified name of the function/method to call,
                or None if unavailable (examples: 'builtins.open', 'typing.Mapping.get')
            object_type: If callable_name refers to a method, the type of the object
                on which the method is being called
        """
        arg_messages = arg_messages or self.msg
        callee = get_proper_type(callee)

        if isinstance(callee, CallableType):
            return self.check_callable_call(callee, args, arg_kinds, context, arg_names,
                                            callable_node, arg_messages, callable_name,
                                            object_type)
        elif isinstance(callee, Overloaded):
            return self.check_overload_call(callee, args, arg_kinds, arg_names, callable_name,
                                            object_type, context, arg_messages)
        elif isinstance(callee, AnyType) or not self.chk.in_checked_function():
            return self.check_any_type_call(args, callee)
        elif isinstance(callee, UnionType):
            return self.check_union_call(callee, args, arg_kinds, arg_names, context, arg_messages)
        elif isinstance(callee, Instance):
            call_function = analyze_member_access('__call__', callee, context, is_lvalue=False,
                                                  is_super=False, is_operator=True, msg=self.msg,
                                                  original_type=callee, chk=self.chk,
                                                  in_literal_context=self.is_literal_context())
            callable_name = callee.type.fullname + ".__call__"
            # Apply method signature hook, if one exists
            call_function = self.transform_callee_type(
                callable_name, call_function, args, arg_kinds, context, arg_names, callee)
            result = self.check_call(call_function, args, arg_kinds, context, arg_names,
                                     callable_node, arg_messages, callable_name, callee)
            if callable_node:
                # check_call() stored "call_function" as the type, which is incorrect.
                # Override the type.
                self.chk.store_type(callable_node, callee)
            return result
        elif isinstance(callee, TypeVarType):
            return self.check_call(callee.upper_bound, args, arg_kinds, context, arg_names,
                                   callable_node, arg_messages)
        elif isinstance(callee, TypeType):
            item = self.analyze_type_type_callee(callee.item, context)
            return self.check_call(item, args, arg_kinds, context, arg_names,
                                   callable_node, arg_messages)
        elif isinstance(callee, TupleType):
            return self.check_call(tuple_fallback(callee), args, arg_kinds, context,
                                   arg_names, callable_node, arg_messages, callable_name,
                                   object_type)
        else:
            return self.msg.not_callable(callee, context), AnyType(TypeOfAny.from_error)

    def check_callable_call(self,
                            callee: CallableType,
                            args: List[Expression],
                            arg_kinds: List[ArgKind],
                            context: Context,
                            arg_names: Optional[Sequence[Optional[str]]],
                            callable_node: Optional[Expression],
                            arg_messages: MessageBuilder,
                            callable_name: Optional[str],
                            object_type: Optional[Type]) -> Tuple[Type, Type]:
        """Type check a call that targets a callable value.

        See the docstring of check_call for more information.
        """
        if callable_name is None and callee.name:
            callable_name = callee.name
        ret_type = get_proper_type(callee.ret_type)
        if callee.is_type_obj() and isinstance(ret_type, Instance):
            callable_name = ret_type.type.fullname
        if isinstance(callable_node, RefExpr) and callable_node.fullname in ENUM_BASES:
            # An Enum() call that failed SemanticAnalyzerPass2.check_enum_call().
            return callee.ret_type, callee

        if (callee.is_type_obj() and callee.type_object().is_abstract
                # Exception for Type[...]
                and not callee.from_type_type
                and not callee.type_object().fallback_to_any):
            type = callee.type_object()
            self.msg.cannot_instantiate_abstract_class(
                callee.type_object().name, type.abstract_attributes,
                context)
        elif (callee.is_type_obj() and callee.type_object().is_protocol
              # Exception for Type[...]
              and not callee.from_type_type):
            self.chk.fail(message_registry.CANNOT_INSTANTIATE_PROTOCOL
                          .format(callee.type_object().name), context)

        formal_to_actual = map_actuals_to_formals(
            arg_kinds, arg_names,
            callee.arg_kinds, callee.arg_names,
            lambda i: self.accept(args[i]))

        if callee.is_generic():
            need_refresh = any(isinstance(v, ParamSpecType) for v in callee.variables)
            callee = freshen_function_type_vars(callee)
            callee = self.infer_function_type_arguments_using_context(
                callee, context)
            callee = self.infer_function_type_arguments(
                callee, args, arg_kinds, formal_to_actual, context)
            if need_refresh:
                # Argument kinds etc. may have changed due to
                # ParamSpec variables being replaced with an arbitrary
                # number of arguments; recalculate actual-to-formal map
                formal_to_actual = map_actuals_to_formals(
                    arg_kinds, arg_names,
                    callee.arg_kinds, callee.arg_names,
                    lambda i: self.accept(args[i]))

        param_spec = callee.param_spec()
        if param_spec is not None and arg_kinds == [ARG_STAR, ARG_STAR2]:
            arg1 = self.accept(args[0])
            arg2 = self.accept(args[1])
            if (isinstance(arg1, ParamSpecType)
                    and isinstance(arg2, ParamSpecType)
                    and arg1.flavor == ParamSpecFlavor.ARGS
                    and arg2.flavor == ParamSpecFlavor.KWARGS
                    and arg1.id == arg2.id == param_spec.id):
                return callee.ret_type, callee

        arg_types = self.infer_arg_types_in_context(
            callee, args, arg_kinds, formal_to_actual)

        self.check_argument_count(callee, arg_types, arg_kinds,
                                  arg_names, formal_to_actual, context, self.msg)

        self.check_argument_types(arg_types, arg_kinds, args, callee, formal_to_actual, context,
                                  messages=arg_messages, object_type=object_type)

        if (callee.is_type_obj() and (len(arg_types) == 1)
                and is_equivalent(callee.ret_type, self.named_type('builtins.type'))):
            callee = callee.copy_modified(ret_type=TypeType.make_normalized(arg_types[0]))

        if callable_node:
            # Store the inferred callable type.
            self.chk.store_type(callable_node, callee)

        if (callable_name
                and ((object_type is None and self.plugin.get_function_hook(callable_name))
                     or (object_type is not None
                         and self.plugin.get_method_hook(callable_name)))):
            new_ret_type = self.apply_function_plugin(
                callee, arg_kinds, arg_types, arg_names, formal_to_actual, args,
                callable_name, object_type, context)
            callee = callee.copy_modified(ret_type=new_ret_type)
        return callee.ret_type, callee

    def analyze_type_type_callee(self, item: ProperType, context: Context) -> Type:
        """Analyze the callee X in X(...) where X is Type[item].

        Return a Y that we can pass to check_call(Y, ...).
        """
        if isinstance(item, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=item)
        if isinstance(item, Instance):
            res = type_object_type(item.type, self.named_type)
            if isinstance(res, CallableType):
                res = res.copy_modified(from_type_type=True)
            expanded = get_proper_type(expand_type_by_instance(res, item))
            if isinstance(expanded, CallableType):
                # Callee of the form Type[...] should never be generic, only
                # proper class objects can be.
                expanded = expanded.copy_modified(variables=[])
            return expanded
        if isinstance(item, UnionType):
            return UnionType([self.analyze_type_type_callee(get_proper_type(tp), context)
                              for tp in item.relevant_items()], item.line)
        if isinstance(item, TypeVarType):
            # Pretend we're calling the typevar's upper bound,
            # i.e. its constructor (a poor approximation for reality,
            # but better than AnyType...), but replace the return type
            # with typevar.
            callee = self.analyze_type_type_callee(get_proper_type(item.upper_bound), context)
            callee = get_proper_type(callee)
            if isinstance(callee, CallableType):
                callee = callee.copy_modified(ret_type=item)
            elif isinstance(callee, Overloaded):
                callee = Overloaded([c.copy_modified(ret_type=item)
                                     for c in callee.items])
            return callee
        # We support Type of namedtuples but not of tuples in general
        if (isinstance(item, TupleType)
                and tuple_fallback(item).type.fullname != 'builtins.tuple'):
            return self.analyze_type_type_callee(tuple_fallback(item), context)

        self.msg.unsupported_type_type(item, context)
        return AnyType(TypeOfAny.from_error)

    def infer_arg_types_in_empty_context(self, args: List[Expression]) -> List[Type]:
        """Infer argument expression types in an empty context.

        In short, we basically recurse on each argument without considering
        in what context the argument was called.
        """
        res: List[Type] = []

        for arg in args:
            arg_type = self.accept(arg)
            if has_erased_component(arg_type):
                res.append(NoneType())
            else:
                res.append(arg_type)
        return res

    def infer_arg_types_in_context(
            self, callee: CallableType, args: List[Expression], arg_kinds: List[ArgKind],
            formal_to_actual: List[List[int]]) -> List[Type]:
        """Infer argument expression types using a callable type as context.

        For example, if callee argument 2 has type List[int], infer the
        argument expression with List[int] type context.

        Returns the inferred types of *actual arguments*.
        """
        res: List[Optional[Type]] = [None] * len(args)

        for i, actuals in enumerate(formal_to_actual):
            for ai in actuals:
                if not arg_kinds[ai].is_star():
                    res[ai] = self.accept(args[ai], callee.arg_types[i])

        # Fill in the rest of the argument types.
        for i, t in enumerate(res):
            if not t:
                res[i] = self.accept(args[i])
        assert all(tp is not None for tp in res)
        return cast(List[Type], res)

    def infer_function_type_arguments_using_context(
            self, callable: CallableType, error_context: Context) -> CallableType:
        """Unify callable return type to type context to infer type vars.

        For example, if the return type is set[t] where 't' is a type variable
        of callable, and if the context is set[int], return callable modified
        by substituting 't' with 'int'.
        """
        ctx = self.type_context[-1]
        if not ctx:
            return callable
        # The return type may have references to type metavariables that
        # we are inferring right now. We must consider them as indeterminate
        # and they are not potential results; thus we replace them with the
        # special ErasedType type. On the other hand, class type variables are
        # valid results.
        erased_ctx = replace_meta_vars(ctx, ErasedType())
        ret_type = callable.ret_type
        if is_optional(ret_type) and is_optional(ctx):
            # If both the context and the return type are optional, unwrap the optional,
            # since in 99% cases this is what a user expects. In other words, we replace
            #     Optional[T] <: Optional[int]
            # with
            #     T <: int
            # while the former would infer T <: Optional[int].
            ret_type = remove_optional(ret_type)
            erased_ctx = remove_optional(erased_ctx)
            #
            # TODO: Instead of this hack and the one below, we need to use outer and
            # inner contexts at the same time. This is however not easy because of two
            # reasons:
            #   * We need to support constraints like [1 <: 2, 2 <: X], i.e. with variables
            #     on both sides. (This is not too hard.)
            #   * We need to update all the inference "infrastructure", so that all
            #     variables in an expression are inferred at the same time.
            #     (And this is hard, also we need to be careful with lambdas that require
            #     two passes.)
        if isinstance(ret_type, TypeVarType):
            # Another special case: the return type is a type variable. If it's unrestricted,
            # we could infer a too general type for the type variable if we use context,
            # and this could result in confusing and spurious type errors elsewhere.
            #
            # So we give up and just use function arguments for type inference, with just two
            # exceptions:
            #
            # 1. If the context is a generic instance type, actually use it as context, as
            #    this *seems* to usually be the reasonable thing to do.
            #
            #    See also github issues #462 and #360.
            #
            # 2. If the context is some literal type, we want to "propagate" that information
            #    down so that we infer a more precise type for literal expressions. For example,
            #    the expression `3` normally has an inferred type of `builtins.int`: but if it's
            #    in a literal context like below, we want it to infer `Literal[3]` instead.
            #
            #        def expects_literal(x: Literal[3]) -> None: pass
            #        def identity(x: T) -> T: return x
            #
            #        expects_literal(identity(3))  # Should type-check
            if not is_generic_instance(ctx) and not is_literal_type_like(ctx):
                return callable.copy_modified()
        args = infer_type_arguments(callable.type_var_ids(), ret_type, erased_ctx)
        # Only substitute non-Uninhabited and non-erased types.
        new_args: List[Optional[Type]] = []
        for arg in args:
            if has_uninhabited_component(arg) or has_erased_component(arg):
                new_args.append(None)
            else:
                new_args.append(arg)
        # Don't show errors after we have only used the outer context for inference.
        # We will use argument context to infer more variables.
        return self.apply_generic_arguments(callable, new_args, error_context,
                                            skip_unsatisfied=True)

    def infer_function_type_arguments(self, callee_type: CallableType,
                                      args: List[Expression],
                                      arg_kinds: List[ArgKind],
                                      formal_to_actual: List[List[int]],
                                      context: Context) -> CallableType:
        """Infer the type arguments for a generic callee type.

        Infer based on the types of arguments.

        Return a derived callable type that has the arguments applied.
        """
        if self.chk.in_checked_function():
            # Disable type errors during type inference. There may be errors
            # due to partial available context information at this time, but
            # these errors can be safely ignored as the arguments will be
            # inferred again later.
            with self.msg.disable_errors():
                arg_types = self.infer_arg_types_in_context(
                    callee_type, args, arg_kinds, formal_to_actual)

            arg_pass_nums = self.get_arg_infer_passes(
                callee_type.arg_types, formal_to_actual, len(args))

            pass1_args: List[Optional[Type]] = []
            for i, arg in enumerate(arg_types):
                if arg_pass_nums[i] > 1:
                    pass1_args.append(None)
                else:
                    pass1_args.append(arg)

            inferred_args = infer_function_type_arguments(
                callee_type, pass1_args, arg_kinds, formal_to_actual,
                context=self.argument_infer_context(),
                strict=self.chk.in_checked_function())

            if 2 in arg_pass_nums:
                # Second pass of type inference.
                (callee_type,
                 inferred_args) = self.infer_function_type_arguments_pass2(
                    callee_type, args, arg_kinds, formal_to_actual,
                    inferred_args, context)

            if callee_type.special_sig == 'dict' and len(inferred_args) == 2 and (
                    ARG_NAMED in arg_kinds or ARG_STAR2 in arg_kinds):
                # HACK: Infer str key type for dict(...) with keyword args. The type system
                #       can't represent this so we special case it, as this is a pretty common
                #       thing. This doesn't quite work with all possible subclasses of dict
                #       if they shuffle type variables around, as we assume that there is a 1-1
                #       correspondence with dict type variables. This is a marginal issue and
                #       a little tricky to fix so it's left unfixed for now.
                first_arg = get_proper_type(inferred_args[0])
                if isinstance(first_arg, (NoneType, UninhabitedType)):
                    inferred_args[0] = self.named_type('builtins.str')
                elif not first_arg or not is_subtype(self.named_type('builtins.str'), first_arg):
                    self.chk.fail(message_registry.KEYWORD_ARGUMENT_REQUIRES_STR_KEY_TYPE,
                                  context)
        else:
            # In dynamically typed functions use implicit 'Any' types for
            # type variables.
            inferred_args = [AnyType(TypeOfAny.unannotated)] * len(callee_type.variables)
        return self.apply_inferred_arguments(callee_type, inferred_args,
                                             context)

    def infer_function_type_arguments_pass2(
            self, callee_type: CallableType,
            args: List[Expression],
            arg_kinds: List[ArgKind],
            formal_to_actual: List[List[int]],
            old_inferred_args: Sequence[Optional[Type]],
            context: Context) -> Tuple[CallableType, List[Optional[Type]]]:
        """Perform second pass of generic function type argument inference.

        The second pass is needed for arguments with types such as Callable[[T], S],
        where both T and S are type variables, when the actual argument is a
        lambda with inferred types.  The idea is to infer the type variable T
        in the first pass (based on the types of other arguments).  This lets
        us infer the argument and return type of the lambda expression and
        thus also the type variable S in this second pass.

        Return (the callee with type vars applied, inferred actual arg types).
        """
        # None or erased types in inferred types mean that there was not enough
        # information to infer the argument. Replace them with None values so
        # that they are not applied yet below.
        inferred_args = list(old_inferred_args)
        for i, arg in enumerate(get_proper_types(inferred_args)):
            if isinstance(arg, (NoneType, UninhabitedType)) or has_erased_component(arg):
                inferred_args[i] = None
        callee_type = self.apply_generic_arguments(callee_type, inferred_args, context)

        arg_types = self.infer_arg_types_in_context(
            callee_type, args, arg_kinds, formal_to_actual)

        inferred_args = infer_function_type_arguments(
            callee_type, arg_types, arg_kinds, formal_to_actual,
            context=self.argument_infer_context(),
        )

        return callee_type, inferred_args

    def argument_infer_context(self) -> ArgumentInferContext:
        return ArgumentInferContext(
            self.chk.named_type('typing.Mapping'),
            self.chk.named_type('typing.Iterable'),
        )

    def get_arg_infer_passes(self, arg_types: List[Type],
                             formal_to_actual: List[List[int]],
                             num_actuals: int) -> List[int]:
        """Return pass numbers for args for two-pass argument type inference.

        For each actual, the pass number is either 1 (first pass) or 2 (second
        pass).

        Two-pass argument type inference primarily lets us infer types of
        lambdas more effectively.
        """
        res = [1] * num_actuals
        for i, arg in enumerate(arg_types):
            if arg.accept(ArgInferSecondPassQuery()):
                for j in formal_to_actual[i]:
                    res[j] = 2
        return res

    def apply_inferred_arguments(self, callee_type: CallableType,
                                 inferred_args: Sequence[Optional[Type]],
                                 context: Context) -> CallableType:
        """Apply inferred values of type arguments to a generic function.

        Inferred_args contains the values of function type arguments.
        """
        # Report error if some of the variables could not be solved. In that
        # case assume that all variables have type Any to avoid extra
        # bogus error messages.
        for i, inferred_type in enumerate(inferred_args):
            if not inferred_type or has_erased_component(inferred_type):
                # Could not infer a non-trivial type for a type variable.
                self.msg.could_not_infer_type_arguments(
                    callee_type, i + 1, context)
                inferred_args = [AnyType(TypeOfAny.from_error)] * len(inferred_args)
        # Apply the inferred types to the function type. In this case the
        # return type must be CallableType, since we give the right number of type
        # arguments.
        return self.apply_generic_arguments(callee_type, inferred_args, context)

    def check_argument_count(self,
                             callee: CallableType,
                             actual_types: List[Type],
                             actual_kinds: List[ArgKind],
                             actual_names: Optional[Sequence[Optional[str]]],
                             formal_to_actual: List[List[int]],
                             context: Optional[Context],
                             messages: Optional[MessageBuilder]) -> bool:
        """Check that there is a value for all required arguments to a function.

        Also check that there are no duplicate values for arguments. Report found errors
        using 'messages' if it's not None. If 'messages' is given, 'context' must also be given.

        Return False if there were any errors. Otherwise return True
        """
        if messages:
            assert context, "Internal error: messages given without context"
        elif context is None:
            # Avoid "is None" checks
            context = TempNode(AnyType(TypeOfAny.special_form))

        # TODO(jukka): We could return as soon as we find an error if messages is None.

        # Collect list of all actual arguments matched to formal arguments.
        all_actuals: List[int] = []
        for actuals in formal_to_actual:
            all_actuals.extend(actuals)

        ok, is_unexpected_arg_error = self.check_for_extra_actual_arguments(
            callee, actual_types, actual_kinds, actual_names, all_actuals, context, messages)

        # Check for too many or few values for formals.
        for i, kind in enumerate(callee.arg_kinds):
            if kind.is_required() and not formal_to_actual[i] and not is_unexpected_arg_error:
                # No actual for a mandatory formal
                if messages:
                    if kind.is_positional():
                        messages.too_few_arguments(callee, context, actual_names)
                    else:
                        argname = callee.arg_names[i] or "?"
                        messages.missing_named_argument(callee, context, argname)
                ok = False
            elif not kind.is_star() and is_duplicate_mapping(
                    formal_to_actual[i], actual_types, actual_kinds):
                if (self.chk.in_checked_function() or
                        isinstance(get_proper_type(actual_types[formal_to_actual[i][0]]),
                                   TupleType)):
                    if messages:
                        messages.duplicate_argument_value(callee, i, context)
                    ok = False
            elif (kind.is_named() and formal_to_actual[i] and
                  actual_kinds[formal_to_actual[i][0]] not in [nodes.ARG_NAMED, nodes.ARG_STAR2]):
                # Positional argument when expecting a keyword argument.
                if messages:
                    messages.too_many_positional_arguments(callee, context)
                ok = False
        return ok

    def check_for_extra_actual_arguments(self,
                                         callee: CallableType,
                                         actual_types: List[Type],
                                         actual_kinds: List[ArgKind],
                                         actual_names: Optional[Sequence[Optional[str]]],
                                         all_actuals: List[int],
                                         context: Context,
                                         messages: Optional[MessageBuilder]) -> Tuple[bool, bool]:
        """Check for extra actual arguments.

        Return tuple (was everything ok,
                      was there an extra keyword argument error [used to avoid duplicate errors]).
        """

        is_unexpected_arg_error = False  # Keep track of errors to avoid duplicate errors
        ok = True  # False if we've found any error

        for i, kind in enumerate(actual_kinds):
            if (i not in all_actuals and
                    # We accept the other iterables than tuple (including Any)
                    # as star arguments because they could be empty, resulting no arguments.
                    (kind != nodes.ARG_STAR or is_non_empty_tuple(actual_types[i])) and
                    # Accept all types for double-starred arguments, because they could be empty
                    # dictionaries and we can't tell it from their types
                    kind != nodes.ARG_STAR2):
                # Extra actual: not matched by a formal argument.
                ok = False
                if kind != nodes.ARG_NAMED:
                    if messages:
                        messages.too_many_arguments(callee, context)
                else:
                    if messages:
                        assert actual_names, "Internal error: named kinds without names given"
                        act_name = actual_names[i]
                        assert act_name is not None
                        act_type = actual_types[i]
                        messages.unexpected_keyword_argument(callee, act_name, act_type, context)
                    is_unexpected_arg_error = True
            elif ((kind == nodes.ARG_STAR and nodes.ARG_STAR not in callee.arg_kinds)
                  or kind == nodes.ARG_STAR2):
                actual_type = get_proper_type(actual_types[i])
                if isinstance(actual_type, (TupleType, TypedDictType)):
                    if all_actuals.count(i) < len(actual_type.items):
                        # Too many tuple/dict items as some did not match.
                        if messages:
                            if (kind != nodes.ARG_STAR2
                                    or not isinstance(actual_type, TypedDictType)):
                                messages.too_many_arguments(callee, context)
                            else:
                                messages.too_many_arguments_from_typed_dict(callee, actual_type,
                                                                            context)
                                is_unexpected_arg_error = True
                        ok = False
                # *args/**kwargs can be applied even if the function takes a fixed
                # number of positional arguments. This may succeed at runtime.

        return ok, is_unexpected_arg_error

    def check_argument_types(self,
                             arg_types: List[Type],
                             arg_kinds: List[ArgKind],
                             args: List[Expression],
                             callee: CallableType,
                             formal_to_actual: List[List[int]],
                             context: Context,
                             messages: Optional[MessageBuilder] = None,
                             check_arg: Optional[ArgChecker] = None,
                             object_type: Optional[Type] = None) -> None:
        """Check argument types against a callable type.

        Report errors if the argument types are not compatible.

        The check_call docstring describes some of the arguments.
        """
        messages = messages or self.msg
        check_arg = check_arg or self.check_arg
        # Keep track of consumed tuple *arg items.
        mapper = ArgTypeExpander(self.argument_infer_context())
        for i, actuals in enumerate(formal_to_actual):
            for actual in actuals:
                actual_type = arg_types[actual]
                if actual_type is None:
                    continue  # Some kind of error was already reported.
                actual_kind = arg_kinds[actual]
                # Check that a *arg is valid as varargs.
                if (actual_kind == nodes.ARG_STAR and
                        not self.is_valid_var_arg(actual_type)):
                    messages.invalid_var_arg(actual_type, context)
                if (actual_kind == nodes.ARG_STAR2 and
                        not self.is_valid_keyword_var_arg(actual_type)):
                    is_mapping = is_subtype(actual_type, self.chk.named_type('typing.Mapping'))
                    messages.invalid_keyword_var_arg(actual_type, is_mapping, context)
                expanded_actual = mapper.expand_actual_type(
                    actual_type, actual_kind,
                    callee.arg_names[i], callee.arg_kinds[i])
                check_arg(expanded_actual, actual_type, arg_kinds[actual],
                          callee.arg_types[i],
                          actual + 1, i + 1, callee, object_type, args[actual], context, messages)

    def check_arg(self,
                  caller_type: Type,
                  original_caller_type: Type,
                  caller_kind: ArgKind,
                  callee_type: Type,
                  n: int,
                  m: int,
                  callee: CallableType,
                  object_type: Optional[Type],
                  context: Context,
                  outer_context: Context,
                  messages: MessageBuilder) -> None:
        """Check the type of a single argument in a call."""
        caller_type = get_proper_type(caller_type)
        original_caller_type = get_proper_type(original_caller_type)
        callee_type = get_proper_type(callee_type)

        if isinstance(caller_type, DeletedType):
            messages.deleted_as_rvalue(caller_type, context)
        # Only non-abstract non-protocol class can be given where Type[...] is expected...
        elif (isinstance(caller_type, CallableType) and isinstance(callee_type, TypeType) and
              caller_type.is_type_obj() and
              (caller_type.type_object().is_abstract or caller_type.type_object().is_protocol) and
              isinstance(callee_type.item, Instance) and
              (callee_type.item.type.is_abstract or callee_type.item.type.is_protocol)):
            self.msg.concrete_only_call(callee_type, context)
        elif not is_subtype(caller_type, callee_type):
            if self.chk.should_suppress_optional_error([caller_type, callee_type]):
                return
            code = messages.incompatible_argument(n,
                                                  m,
                                                  callee,
                                                  original_caller_type,
                                                  caller_kind,
                                                  object_type=object_type,
                                                  context=context,
                                                  outer_context=outer_context)
            messages.incompatible_argument_note(original_caller_type, callee_type, context,
                                                code=code)

    def check_overload_call(self,
                            callee: Overloaded,
                            args: List[Expression],
                            arg_kinds: List[ArgKind],
                            arg_names: Optional[Sequence[Optional[str]]],
                            callable_name: Optional[str],
                            object_type: Optional[Type],
                            context: Context,
                            arg_messages: MessageBuilder) -> Tuple[Type, Type]:
        """Checks a call to an overloaded function."""
        arg_types = self.infer_arg_types_in_empty_context(args)
        # Step 1: Filter call targets to remove ones where the argument counts don't match
        plausible_targets = self.plausible_overload_call_targets(arg_types, arg_kinds,
                                                                 arg_names, callee)

        # Step 2: If the arguments contain a union, we try performing union math first,
        #         instead of picking the first matching overload.
        #         This is because picking the first overload often ends up being too greedy:
        #         for example, when we have a fallback alternative that accepts an unrestricted
        #         typevar. See https://github.com/python/mypy/issues/4063 for related discussion.
        erased_targets: Optional[List[CallableType]] = None
        unioned_result: Optional[Tuple[Type, Type]] = None
        union_interrupted = False  # did we try all union combinations?
        if any(self.real_union(arg) for arg in arg_types):
            unioned_errors = arg_messages.clean_copy()
            try:
                unioned_return = self.union_overload_result(plausible_targets, args,
                                                            arg_types, arg_kinds, arg_names,
                                                            callable_name, object_type,
                                                            context,
                                                            arg_messages=unioned_errors)
            except TooManyUnions:
                union_interrupted = True
            else:
                # Record if we succeeded. Next we need to see if maybe normal procedure
                # gives a narrower type.
                if unioned_return:
                    # TODO: fix signature of zip() in typeshed.
                    returns, inferred_types = cast(Any, zip)(*unioned_return)
                    # Note that we use `combine_function_signatures` instead of just returning
                    # a union of inferred callables because for example a call
                    # Union[int -> int, str -> str](Union[int, str]) is invalid and
                    # we don't want to introduce internal inconsistencies.
                    unioned_result = (make_simplified_union(list(returns),
                                                            context.line,
                                                            context.column),
                                      self.combine_function_signatures(inferred_types))

        # Step 3: We try checking each branch one-by-one.
        inferred_result = self.infer_overload_return_type(plausible_targets, args, arg_types,
                                                          arg_kinds, arg_names, callable_name,
                                                          object_type, context, arg_messages)
        # If any of checks succeed, stop early.
        if inferred_result is not None and unioned_result is not None:
            # Both unioned and direct checks succeeded, choose the more precise type.
            if (is_subtype(inferred_result[0], unioned_result[0]) and
                    not isinstance(get_proper_type(inferred_result[0]), AnyType)):
                return inferred_result
            return unioned_result
        elif unioned_result is not None:
            return unioned_result
        elif inferred_result is not None:
            return inferred_result

        # Step 4: Failure. At this point, we know there is no match. We fall back to trying
        #         to find a somewhat plausible overload target using the erased types
        #         so we can produce a nice error message.
        #
        #         For example, suppose the user passes a value of type 'List[str]' into an
        #         overload with signatures f(x: int) -> int and f(x: List[int]) -> List[int].
        #
        #         Neither alternative matches, but we can guess the user probably wants the
        #         second one.
        erased_targets = self.overload_erased_call_targets(plausible_targets, arg_types,
                                                           arg_kinds, arg_names, args, context)

        # Step 5: We try and infer a second-best alternative if possible. If not, fall back
        #         to using 'Any'.
        if len(erased_targets) > 0:
            # Pick the first plausible erased target as the fallback
            # TODO: Adjust the error message here to make it clear there was no match.
            #       In order to do this, we need to find a clean way of associating
            #       a note with whatever error message 'self.check_call' will generate.
            #       In particular, the note's line and column numbers need to be the same
            #       as the error's.
            target: Type = erased_targets[0]
        else:
            # There was no plausible match: give up
            target = AnyType(TypeOfAny.from_error)

            if not self.chk.should_suppress_optional_error(arg_types):
                if not is_operator_method(callable_name):
                    code = None
                else:
                    code = codes.OPERATOR
                arg_messages.no_variant_matches_arguments(
                    callee, arg_types, context, code=code)

        result = self.check_call(target, args, arg_kinds, context, arg_names,
                                 arg_messages=arg_messages,
                                 callable_name=callable_name,
                                 object_type=object_type)
        if union_interrupted:
            self.chk.fail(message_registry.TOO_MANY_UNION_COMBINATIONS, context)
        return result

    def plausible_overload_call_targets(self,
                                        arg_types: List[Type],
                                        arg_kinds: List[ArgKind],
                                        arg_names: Optional[Sequence[Optional[str]]],
                                        overload: Overloaded) -> List[CallableType]:
        """Returns all overload call targets that having matching argument counts.

        If the given args contains a star-arg (*arg or **kwarg argument), this method
        will ensure all star-arg overloads appear at the start of the list, instead
        of their usual location.

        The only exception is if the starred argument is something like a Tuple or a
        NamedTuple, which has a definitive "shape". If so, we don't move the corresponding
        alternative to the front since we can infer a more precise match using the original
        order."""

        def has_shape(typ: Type) -> bool:
            typ = get_proper_type(typ)
            return (isinstance(typ, TupleType) or isinstance(typ, TypedDictType)
                    or (isinstance(typ, Instance) and typ.type.is_named_tuple))

        matches: List[CallableType] = []
        star_matches: List[CallableType] = []

        args_have_var_arg = False
        args_have_kw_arg = False
        for kind, typ in zip(arg_kinds, arg_types):
            if kind == ARG_STAR and not has_shape(typ):
                args_have_var_arg = True
            if kind == ARG_STAR2 and not has_shape(typ):
                args_have_kw_arg = True

        for typ in overload.items:
            formal_to_actual = map_actuals_to_formals(arg_kinds, arg_names,
                                                      typ.arg_kinds, typ.arg_names,
                                                      lambda i: arg_types[i])

            if self.check_argument_count(typ, arg_types, arg_kinds, arg_names,
                                         formal_to_actual, None, None):
                if args_have_var_arg and typ.is_var_arg:
                    star_matches.append(typ)
                elif args_have_kw_arg and typ.is_kw_arg:
                    star_matches.append(typ)
                else:
                    matches.append(typ)

        return star_matches + matches

    def infer_overload_return_type(self,
                                   plausible_targets: List[CallableType],
                                   args: List[Expression],
                                   arg_types: List[Type],
                                   arg_kinds: List[ArgKind],
                                   arg_names: Optional[Sequence[Optional[str]]],
                                   callable_name: Optional[str],
                                   object_type: Optional[Type],
                                   context: Context,
                                   arg_messages: Optional[MessageBuilder] = None,
                                   ) -> Optional[Tuple[Type, Type]]:
        """Attempts to find the first matching callable from the given list.

        If a match is found, returns a tuple containing the result type and the inferred
        callee type. (This tuple is meant to be eventually returned by check_call.)
        If multiple targets match due to ambiguous Any parameters, returns (AnyType, AnyType).
        If no targets match, returns None.

        Assumes all of the given targets have argument counts compatible with the caller.
        """

        arg_messages = self.msg if arg_messages is None else arg_messages
        matches: List[CallableType] = []
        return_types: List[Type] = []
        inferred_types: List[Type] = []
        args_contain_any = any(map(has_any_type, arg_types))

        for typ in plausible_targets:
            overload_messages = self.msg.clean_copy()
            prev_messages = self.msg
            assert self.msg is self.chk.msg
            self.msg = overload_messages
            self.chk.msg = overload_messages
            try:
                # Passing `overload_messages` as the `arg_messages` parameter doesn't
                # seem to reliably catch all possible errors.
                # TODO: Figure out why
                ret_type, infer_type = self.check_call(
                    callee=typ,
                    args=args,
                    arg_kinds=arg_kinds,
                    arg_names=arg_names,
                    context=context,
                    arg_messages=overload_messages,
                    callable_name=callable_name,
                    object_type=object_type)
            finally:
                self.chk.msg = prev_messages
                self.msg = prev_messages

            is_match = not overload_messages.is_errors()
            if is_match:
                # Return early if possible; otherwise record info so we can
                # check for ambiguity due to 'Any' below.
                if not args_contain_any:
                    return ret_type, infer_type
                matches.append(typ)
                return_types.append(ret_type)
                inferred_types.append(infer_type)

        if len(matches) == 0:
            # No match was found
            return None
        elif any_causes_overload_ambiguity(matches, return_types, arg_types, arg_kinds, arg_names):
            # An argument of type or containing the type 'Any' caused ambiguity.
            # We try returning a precise type if we can. If not, we give up and just return 'Any'.
            if all_same_types(return_types):
                return return_types[0], inferred_types[0]
            elif all_same_types([erase_type(typ) for typ in return_types]):
                return erase_type(return_types[0]), erase_type(inferred_types[0])
            else:
                return self.check_call(callee=AnyType(TypeOfAny.special_form),
                                       args=args,
                                       arg_kinds=arg_kinds,
                                       arg_names=arg_names,
                                       context=context,
                                       arg_messages=arg_messages,
                                       callable_name=callable_name,
                                       object_type=object_type)
        else:
            # Success! No ambiguity; return the first match.
            return return_types[0], inferred_types[0]

    def overload_erased_call_targets(self,
                                     plausible_targets: List[CallableType],
                                     arg_types: List[Type],
                                     arg_kinds: List[ArgKind],
                                     arg_names: Optional[Sequence[Optional[str]]],
                                     args: List[Expression],
                                     context: Context) -> List[CallableType]:
        """Returns a list of all targets that match the caller after erasing types.

        Assumes all of the given targets have argument counts compatible with the caller.
        """
        matches: List[CallableType] = []
        for typ in plausible_targets:
            if self.erased_signature_similarity(arg_types, arg_kinds, arg_names, args, typ,
                                                context):
                matches.append(typ)
        return matches

    def union_overload_result(self,
                              plausible_targets: List[CallableType],
                              args: List[Expression],
                              arg_types: List[Type],
                              arg_kinds: List[ArgKind],
                              arg_names: Optional[Sequence[Optional[str]]],
                              callable_name: Optional[str],
                              object_type: Optional[Type],
                              context: Context,
                              arg_messages: Optional[MessageBuilder] = None,
                              level: int = 0
                              ) -> Optional[List[Tuple[Type, Type]]]:
        """Accepts a list of overload signatures and attempts to match calls by destructuring
        the first union.

        Return a list of (<return type>, <inferred variant type>) if call succeeds for every
        item of the desctructured union. Returns None if there is no match.
        """
        # Step 1: If we are already too deep, then stop immediately. Otherwise mypy might
        # hang for long time because of a weird overload call. The caller will get
        # the exception and generate an appropriate note message, if needed.
        if level >= MAX_UNIONS:
            raise TooManyUnions

        # Step 2: Find position of the first union in arguments. Return the normal inferred
        # type if no more unions left.
        for idx, typ in enumerate(arg_types):
            if self.real_union(typ):
                break
        else:
            # No unions in args, just fall back to normal inference
            with self.type_overrides_set(args, arg_types):
                res = self.infer_overload_return_type(plausible_targets, args, arg_types,
                                                      arg_kinds, arg_names, callable_name,
                                                      object_type, context, arg_messages)
            if res is not None:
                return [res]
            return None

        # Step 3: Try a direct match before splitting to avoid unnecessary union splits
        # and save performance.
        with self.type_overrides_set(args, arg_types):
            direct = self.infer_overload_return_type(plausible_targets, args, arg_types,
                                                     arg_kinds, arg_names, callable_name,
                                                     object_type, context, arg_messages)
        if direct is not None and not isinstance(get_proper_type(direct[0]),
                                                 (UnionType, AnyType)):
            # We only return non-unions soon, to avoid greedy match.
            return [direct]

        # Step 4: Split the first remaining union type in arguments into items and
        # try to match each item individually (recursive).
        first_union = get_proper_type(arg_types[idx])
        assert isinstance(first_union, UnionType)
        res_items = []
        for item in first_union.relevant_items():
            new_arg_types = arg_types.copy()
            new_arg_types[idx] = item
            sub_result = self.union_overload_result(plausible_targets, args, new_arg_types,
                                                    arg_kinds, arg_names, callable_name,
                                                    object_type, context, arg_messages,
                                                    level + 1)
            if sub_result is not None:
                res_items.extend(sub_result)
            else:
                # Some item doesn't match, return soon.
                return None

        # Step 5: If splitting succeeded, then filter out duplicate items before returning.
        seen: Set[Tuple[Type, Type]] = set()
        result = []
        for pair in res_items:
            if pair not in seen:
                seen.add(pair)
                result.append(pair)
        return result

    def real_union(self, typ: Type) -> bool:
        typ = get_proper_type(typ)
        return isinstance(typ, UnionType) and len(typ.relevant_items()) > 1

    @contextmanager
    def type_overrides_set(self, exprs: Sequence[Expression],
                           overrides: Sequence[Type]) -> Iterator[None]:
        """Set _temporary_ type overrides for given expressions."""
        assert len(exprs) == len(overrides)
        for expr, typ in zip(exprs, overrides):
            self.type_overrides[expr] = typ
        try:
            yield
        finally:
            for expr in exprs:
                del self.type_overrides[expr]

    def combine_function_signatures(self, types: Sequence[Type]) -> Union[AnyType, CallableType]:
        """Accepts a list of function signatures and attempts to combine them together into a
        new CallableType consisting of the union of all of the given arguments and return types.

        If there is at least one non-callable type, return Any (this can happen if there is
        an ambiguity because of Any in arguments).
        """
        assert types, "Trying to merge no callables"
        types = get_proper_types(types)
        if not all(isinstance(c, CallableType) for c in types):
            return AnyType(TypeOfAny.special_form)
        callables = cast(Sequence[CallableType], types)
        if len(callables) == 1:
            return callables[0]

        # Note: we are assuming here that if a user uses some TypeVar 'T' in
        # two different functions, they meant for that TypeVar to mean the
        # same thing.
        #
        # This function will make sure that all instances of that TypeVar 'T'
        # refer to the same underlying TypeVarType objects to simplify the union-ing
        # logic below.
        #
        # (If the user did *not* mean for 'T' to be consistently bound to the
        # same type in their overloads, well, their code is probably too
        # confusing and ought to be re-written anyways.)
        callables, variables = merge_typevars_in_callables_by_name(callables)

        new_args: List[List[Type]] = [[] for _ in range(len(callables[0].arg_types))]
        new_kinds = list(callables[0].arg_kinds)
        new_returns: List[Type] = []

        too_complex = False
        for target in callables:
            # We fall back to Callable[..., Union[<returns>]] if the functions do not have
            # the exact same signature. The only exception is if one arg is optional and
            # the other is positional: in that case, we continue unioning (and expect a
            # positional arg).
            # TODO: Enhance the merging logic to handle a wider variety of signatures.
            if len(new_kinds) != len(target.arg_kinds):
                too_complex = True
                break
            for i, (new_kind, target_kind) in enumerate(zip(new_kinds, target.arg_kinds)):
                if new_kind == target_kind:
                    continue
                elif new_kind.is_positional() and target_kind.is_positional():
                    new_kinds[i] = ARG_POS
                else:
                    too_complex = True
                    break

            if too_complex:
                break  # outer loop

            for i, arg in enumerate(target.arg_types):
                new_args[i].append(arg)
            new_returns.append(target.ret_type)

        union_return = make_simplified_union(new_returns)
        if too_complex:
            any = AnyType(TypeOfAny.special_form)
            return callables[0].copy_modified(
                arg_types=[any, any],
                arg_kinds=[ARG_STAR, ARG_STAR2],
                arg_names=[None, None],
                ret_type=union_return,
                variables=variables,
                implicit=True)

        final_args = []
        for args_list in new_args:
            new_type = make_simplified_union(args_list)
            final_args.append(new_type)

        return callables[0].copy_modified(
            arg_types=final_args,
            arg_kinds=new_kinds,
            ret_type=union_return,
            variables=variables,
            implicit=True)

    def erased_signature_similarity(self,
                                    arg_types: List[Type],
                                    arg_kinds: List[ArgKind],
                                    arg_names: Optional[Sequence[Optional[str]]],
                                    args: List[Expression],
                                    callee: CallableType,
                                    context: Context) -> bool:
        """Determine whether arguments could match the signature at runtime, after
        erasing types."""
        formal_to_actual = map_actuals_to_formals(arg_kinds,
                                                  arg_names,
                                                  callee.arg_kinds,
                                                  callee.arg_names,
                                                  lambda i: arg_types[i])

        if not self.check_argument_count(callee, arg_types, arg_kinds, arg_names,
                                         formal_to_actual, None, None):
            # Too few or many arguments -> no match.
            return False

        def check_arg(caller_type: Type,
                      original_ccaller_type: Type,
                      caller_kind: ArgKind,
                      callee_type: Type,
                      n: int,
                      m: int,
                      callee: CallableType,
                      object_type: Optional[Type],
                      context: Context,
                      outer_context: Context,
                      messages: MessageBuilder) -> None:
            if not arg_approximate_similarity(caller_type, callee_type):
                # No match -- exit early since none of the remaining work can change
                # the result.
                raise Finished

        try:
            self.check_argument_types(arg_types, arg_kinds, args, callee,
                                      formal_to_actual, context=context, check_arg=check_arg)
            return True
        except Finished:
            return False

    def apply_generic_arguments(self, callable: CallableType, types: Sequence[Optional[Type]],
                                context: Context, skip_unsatisfied: bool = False) -> CallableType:
        """Simple wrapper around mypy.applytype.apply_generic_arguments."""
        return applytype.apply_generic_arguments(callable, types,
                                                 self.msg.incompatible_typevar_value, context,
                                                 skip_unsatisfied=skip_unsatisfied)

    def check_any_type_call(self, args: List[Expression], callee: Type) -> Tuple[Type, Type]:
        self.infer_arg_types_in_empty_context(args)
        callee = get_proper_type(callee)
        if isinstance(callee, AnyType):
            return (AnyType(TypeOfAny.from_another_any, source_any=callee),
                    AnyType(TypeOfAny.from_another_any, source_any=callee))
        else:
            return AnyType(TypeOfAny.special_form), AnyType(TypeOfAny.special_form)

    def check_union_call(self,
                         callee: UnionType,
                         args: List[Expression],
                         arg_kinds: List[ArgKind],
                         arg_names: Optional[Sequence[Optional[str]]],
                         context: Context,
                         arg_messages: MessageBuilder) -> Tuple[Type, Type]:
        self.msg.disable_type_names += 1
        results = [self.check_call(subtype, args, arg_kinds, context, arg_names,
                                   arg_messages=arg_messages)
                   for subtype in callee.relevant_items()]
        self.msg.disable_type_names -= 1
        return (make_simplified_union([res[0] for res in results]),
                callee)

    def visit_member_expr(self, e: MemberExpr, is_lvalue: bool = False) -> Type:
        """Visit member expression (of form e.id)."""
        self.chk.module_refs.update(extract_refexpr_names(e))
        result = self.analyze_ordinary_member_access(e, is_lvalue)
        return self.narrow_type_from_binder(e, result)

    def analyze_ordinary_member_access(self, e: MemberExpr,
                                       is_lvalue: bool) -> Type:
        """Analyse member expression or member lvalue."""
        if e.kind is not None:
            # This is a reference to a module attribute.
            return self.analyze_ref_expr(e)
        else:
            # This is a reference to a non-module attribute.
            original_type = self.accept(e.expr)
            base = e.expr
            module_symbol_table = None

            if isinstance(base, RefExpr) and isinstance(base.node, MypyFile):
                module_symbol_table = base.node.names

            member_type = analyze_member_access(
                e.name, original_type, e, is_lvalue, False, False,
                self.msg, original_type=original_type, chk=self.chk,
                in_literal_context=self.is_literal_context(),
                module_symbol_table=module_symbol_table)

            return member_type

    def analyze_external_member_access(self, member: str, base_type: Type,
                                       context: Context) -> Type:
        """Analyse member access that is external, i.e. it cannot
        refer to private definitions. Return the result type.
        """
        # TODO remove; no private definitions in mypy
        return analyze_member_access(member, base_type, context, False, False, False,
                                     self.msg, original_type=base_type, chk=self.chk,
                                     in_literal_context=self.is_literal_context())

    def is_literal_context(self) -> bool:
        return is_literal_type_like(self.type_context[-1])

    def infer_literal_expr_type(self, value: LiteralValue, fallback_name: str) -> Type:
        """Analyzes the given literal expression and determines if we should be
        inferring an Instance type, a Literal[...] type, or an Instance that
        remembers the original literal. We...

        1. ...Infer a normal Instance in most circumstances.

        2. ...Infer a Literal[...] if we're in a literal context. For example, if we
           were analyzing the "3" in "foo(3)" where "foo" has a signature of
           "def foo(Literal[3]) -> None", we'd want to infer that the "3" has a
           type of Literal[3] instead of Instance.

        3. ...Infer an Instance that remembers the original Literal if we're declaring
           a Final variable with an inferred type -- for example, "bar" in "bar: Final = 3"
           would be assigned an Instance that remembers it originated from a '3'. See
           the comments in Instance's constructor for more details.
        """
        typ = self.named_type(fallback_name)
        if self.is_literal_context():
            return LiteralType(value=value, fallback=typ)
        else:
            return typ.copy_modified(last_known_value=LiteralType(
                value=value,
                fallback=typ,
                line=typ.line,
                column=typ.column,
            ))

    def concat_tuples(self, left: TupleType, right: TupleType) -> TupleType:
        """Concatenate two fixed length tuples."""
        return TupleType(items=left.items + right.items,
                         fallback=self.named_type('builtins.tuple'))

    def visit_int_expr(self, e: IntExpr) -> Type:
        """Type check an integer literal (trivial)."""
        return self.infer_literal_expr_type(e.value, 'builtins.int')

    def visit_str_expr(self, e: StrExpr) -> Type:
        """Type check a string literal (trivial)."""
        return self.infer_literal_expr_type(e.value, 'builtins.str')

    def visit_bytes_expr(self, e: BytesExpr) -> Type:
        """Type check a bytes literal (trivial)."""
        return self.infer_literal_expr_type(e.value, 'builtins.bytes')

    def visit_unicode_expr(self, e: UnicodeExpr) -> Type:
        """Type check a unicode literal (trivial)."""
        return self.infer_literal_expr_type(e.value, 'builtins.unicode')

    def visit_float_expr(self, e: FloatExpr) -> Type:
        """Type check a float literal (trivial)."""
        return self.named_type('builtins.float')

    def visit_complex_expr(self, e: ComplexExpr) -> Type:
        """Type check a complex literal."""
        return self.named_type('builtins.complex')

    def visit_ellipsis(self, e: EllipsisExpr) -> Type:
        """Type check '...'."""
        if self.chk.options.python_version[0] >= 3:
            return self.named_type('builtins.ellipsis')
        else:
            # '...' is not valid in normal Python 2 code, but it can
            # be used in stubs.  The parser makes sure that we only
            # get this far if we are in a stub, and we can safely
            # return 'object' as ellipsis is special cased elsewhere.
            # The builtins.ellipsis type does not exist in Python 2.
            return self.named_type('builtins.object')

    def visit_op_expr(self, e: OpExpr) -> Type:
        """Type check a binary operator expression."""
        if e.op == 'and' or e.op == 'or':
            return self.check_boolean_op(e, e)
        if e.op == '*' and isinstance(e.left, ListExpr):
            # Expressions of form [...] * e get special type inference.
            return self.check_list_multiply(e)
        if e.op == '%':
            pyversion = self.chk.options.python_version
            if pyversion[0] == 3:
                if isinstance(e.left, BytesExpr) and pyversion[1] >= 5:
                    return self.strfrm_checker.check_str_interpolation(e.left, e.right)
                if isinstance(e.left, StrExpr):
                    return self.strfrm_checker.check_str_interpolation(e.left, e.right)
            elif pyversion[0] <= 2:
                if isinstance(e.left, (StrExpr, BytesExpr, UnicodeExpr)):
                    return self.strfrm_checker.check_str_interpolation(e.left, e.right)
        left_type = self.accept(e.left)

        proper_left_type = get_proper_type(left_type)
        if isinstance(proper_left_type, TupleType) and e.op == '+':
            left_add_method = proper_left_type.partial_fallback.type.get('__add__')
            if left_add_method and left_add_method.fullname == 'builtins.tuple.__add__':
                proper_right_type = get_proper_type(self.accept(e.right))
                if isinstance(proper_right_type, TupleType):
                    right_radd_method = proper_right_type.partial_fallback.type.get('__radd__')
                    if right_radd_method is None:
                        return self.concat_tuples(proper_left_type, proper_right_type)

        if e.op in operators.op_methods:
            method = self.get_operator_method(e.op)
            result, method_type = self.check_op(method, left_type, e.right, e,
                                                allow_reverse=True)
            e.method_type = method_type
            return result
        else:
            raise RuntimeError('Unknown operator {}'.format(e.op))

    def visit_comparison_expr(self, e: ComparisonExpr) -> Type:
        """Type check a comparison expression.

        Comparison expressions are type checked consecutive-pair-wise
        That is, 'a < b > c == d' is check as 'a < b and b > c and c == d'
        """
        result: Optional[Type] = None
        sub_result: Optional[Type] = None

        # Check each consecutive operand pair and their operator
        for left, right, operator in zip(e.operands, e.operands[1:], e.operators):
            left_type = self.accept(left)

            method_type: Optional[mypy.types.Type] = None

            if operator == 'in' or operator == 'not in':
                # If the right operand has partial type, look it up without triggering
                # a "Need type annotation ..." message, as it would be noise.
                right_type = self.find_partial_type_ref_fast_path(right)
                if right_type is None:
                    right_type = self.accept(right)  # Validate the right operand

                # Keep track of whether we get type check errors (these won't be reported, they
                # are just to verify whether something is valid typing wise).
                local_errors = self.msg.clean_copy()
                _, method_type = self.check_method_call_by_name(
                    method='__contains__',
                    base_type=right_type,
                    args=[left],
                    arg_kinds=[ARG_POS],
                    context=e,
                    local_errors=local_errors,
                )
                sub_result = self.bool_type()
                # Container item type for strict type overlap checks. Note: we need to only
                # check for nominal type, because a usual "Unsupported operands for in"
                # will be reported for types incompatible with __contains__().
                # See testCustomContainsCheckStrictEquality for an example.
                cont_type = self.chk.analyze_container_item_type(right_type)
                if isinstance(right_type, PartialType):
                    # We don't really know if this is an error or not, so just shut up.
                    pass
                elif (local_errors.is_errors() and
                    # is_valid_var_arg is True for any Iterable
                        self.is_valid_var_arg(right_type)):
                    _, itertype = self.chk.analyze_iterable_item_type(right)
                    method_type = CallableType(
                        [left_type],
                        [nodes.ARG_POS],
                        [None],
                        self.bool_type(),
                        self.named_type('builtins.function'))
                    if not is_subtype(left_type, itertype):
                        self.msg.unsupported_operand_types('in', left_type, right_type, e)
                # Only show dangerous overlap if there are no other errors.
                elif (not local_errors.is_errors() and cont_type and
                        self.dangerous_comparison(left_type, cont_type,
                                                  original_container=right_type)):
                    self.msg.dangerous_comparison(left_type, cont_type, 'container', e)
                else:
                    self.msg.add_errors(local_errors)
            elif operator in operators.op_methods:
                method = self.get_operator_method(operator)
                err_count = self.msg.errors.total_errors()
                sub_result, method_type = self.check_op(method, left_type, right, e,
                                                        allow_reverse=True)
                # Only show dangerous overlap if there are no other errors. See
                # testCustomEqCheckStrictEquality for an example.
                if self.msg.errors.total_errors() == err_count and operator in ('==', '!='):
                    right_type = self.accept(right)
                    # We suppress the error if there is a custom __eq__() method on either
                    # side. User defined (or even standard library) classes can define this
                    # to return True for comparisons between non-overlapping types.
                    if (not custom_special_method(left_type, '__eq__') and
                            not custom_special_method(right_type, '__eq__')):
                        # Also flag non-overlapping literals in situations like:
                        #    x: Literal['a', 'b']
                        #    if x == 'c':
                        #        ...
                        left_type = try_getting_literal(left_type)
                        right_type = try_getting_literal(right_type)
                        if self.dangerous_comparison(left_type, right_type):
                            self.msg.dangerous_comparison(left_type, right_type, 'equality', e)

            elif operator == 'is' or operator == 'is not':
                right_type = self.accept(right)  # validate the right operand
                sub_result = self.bool_type()
                left_type = try_getting_literal(left_type)
                right_type = try_getting_literal(right_type)
                if self.dangerous_comparison(left_type, right_type):
                    self.msg.dangerous_comparison(left_type, right_type, 'identity', e)
                method_type = None
            else:
                raise RuntimeError('Unknown comparison operator {}'.format(operator))

            e.method_types.append(method_type)

            #  Determine type of boolean-and of result and sub_result
            if result is None:
                result = sub_result
            else:
                result = join.join_types(result, sub_result)

        assert result is not None
        return result

    def find_partial_type_ref_fast_path(self, expr: Expression) -> Optional[Type]:
        """If expression has a partial generic type, return it without additional checks.

        In particular, this does not generate an error about a missing annotation.

        Otherwise, return None.
        """
        if not isinstance(expr, RefExpr):
            return None
        if isinstance(expr.node, Var):
            result = self.analyze_var_ref(expr.node, expr)
            if isinstance(result, PartialType) and result.type is not None:
                self.chk.store_type(expr, self.chk.fixup_partial_type(result))
                return result
        return None

    def dangerous_comparison(self, left: Type, right: Type,
                             original_container: Optional[Type] = None) -> bool:
        """Check for dangerous non-overlapping comparisons like 42 == 'no'.

        The original_container is the original container type for 'in' checks
        (and None for equality checks).

        Rules:
            * X and None are overlapping even in strict-optional mode. This is to allow
            'assert x is not None' for x defined as 'x = None  # type: str' in class body
            (otherwise mypy itself would have couple dozen errors because of this).
            * Optional[X] and Optional[Y] are non-overlapping if X and Y are
            non-overlapping, although technically None is overlap, it is most
            likely an error.
            * Any overlaps with everything, i.e. always safe.
            * Special case: b'abc' in b'cde' is safe.
        """
        if not self.chk.options.strict_equality:
            return False

        left, right = get_proper_types((left, right))

        if self.chk.binder.is_unreachable_warning_suppressed():
            # We are inside a function that contains type variables with value restrictions in
            # its signature. In this case we just suppress all strict-equality checks to avoid
            # false positives for code like:
            #
            #     T = TypeVar('T', str, int)
            #     def f(x: T) -> T:
            #         if x == 0:
            #             ...
            #         return x
            #
            # TODO: find a way of disabling the check only for types resulted from the expansion.
            return False
        if isinstance(left, NoneType) or isinstance(right, NoneType):
            return False
        if isinstance(left, UnionType) and isinstance(right, UnionType):
            left = remove_optional(left)
            right = remove_optional(right)
            left, right = get_proper_types((left, right))
        py2 = self.chk.options.python_version < (3, 0)
        if (original_container and has_bytes_component(original_container, py2) and
                has_bytes_component(left, py2)):
            # We need to special case bytes and bytearray, because 97 in b'abc', b'a' in b'abc',
            # b'a' in bytearray(b'abc') etc. all return True (and we want to show the error only
            # if the check can _never_ be True).
            return False
        if isinstance(left, Instance) and isinstance(right, Instance):
            # Special case some builtin implementations of AbstractSet.
            if (left.type.fullname in OVERLAPPING_TYPES_ALLOWLIST and
                    right.type.fullname in OVERLAPPING_TYPES_ALLOWLIST):
                abstract_set = self.chk.lookup_typeinfo('typing.AbstractSet')
                left = map_instance_to_supertype(left, abstract_set)
                right = map_instance_to_supertype(right, abstract_set)
                return not is_overlapping_types(left.args[0], right.args[0])
        if isinstance(left, LiteralType) and isinstance(right, LiteralType):
            if isinstance(left.value, bool) and isinstance(right.value, bool):
                # Comparing different booleans is not dangerous.
                return False
        return not is_overlapping_types(left, right, ignore_promotions=False)

    def get_operator_method(self, op: str) -> str:
        if op == '/' and self.chk.options.python_version[0] == 2:
            # TODO also check for "from __future__ import division"
            return '__div__'
        else:
            return operators.op_methods[op]

    def check_method_call_by_name(self,
                                  method: str,
                                  base_type: Type,
                                  args: List[Expression],
                                  arg_kinds: List[ArgKind],
                                  context: Context,
                                  local_errors: Optional[MessageBuilder] = None,
                                  original_type: Optional[Type] = None
                                  ) -> Tuple[Type, Type]:
        """Type check a call to a named method on an object.

        Return tuple (result type, inferred method type). The 'original_type'
        is used for error messages.
        """
        local_errors = local_errors or self.msg
        original_type = original_type or base_type
        # Unions are special-cased to allow plugins to act on each element of the union.
        base_type = get_proper_type(base_type)
        if isinstance(base_type, UnionType):
            return self.check_union_method_call_by_name(method, base_type,
                                                        args, arg_kinds,
                                                        context, local_errors, original_type)

        method_type = analyze_member_access(method, base_type, context, False, False, True,
                                            local_errors, original_type=original_type,
                                            chk=self.chk,
                                            in_literal_context=self.is_literal_context())
        return self.check_method_call(
            method, base_type, method_type, args, arg_kinds, context, local_errors)

    def check_union_method_call_by_name(self,
                                        method: str,
                                        base_type: UnionType,
                                        args: List[Expression],
                                        arg_kinds: List[ArgKind],
                                        context: Context,
                                        local_errors: MessageBuilder,
                                        original_type: Optional[Type] = None
                                        ) -> Tuple[Type, Type]:
        """Type check a call to a named method on an object with union type.

        This essentially checks the call using check_method_call_by_name() for each
        union item and unions the result. We do this to allow plugins to act on
        individual union items.
        """
        res: List[Type] = []
        meth_res: List[Type] = []
        for typ in base_type.relevant_items():
            # Format error messages consistently with
            # mypy.checkmember.analyze_union_member_access().
            local_errors.disable_type_names += 1
            item, meth_item = self.check_method_call_by_name(method, typ, args, arg_kinds,
                                                             context, local_errors,
                                                             original_type)
            local_errors.disable_type_names -= 1
            res.append(item)
            meth_res.append(meth_item)
        return make_simplified_union(res), make_simplified_union(meth_res)

    def check_method_call(self,
                          method_name: str,
                          base_type: Type,
                          method_type: Type,
                          args: List[Expression],
                          arg_kinds: List[ArgKind],
                          context: Context,
                          local_errors: Optional[MessageBuilder] = None) -> Tuple[Type, Type]:
        """Type check a call to a method with the given name and type on an object.

        Return tuple (result type, inferred method type).
        """
        callable_name = self.method_fullname(base_type, method_name)
        object_type = base_type if callable_name is not None else None

        # Try to refine the method signature using plugin hooks before checking the call.
        method_type = self.transform_callee_type(
            callable_name, method_type, args, arg_kinds, context, object_type=object_type)

        return self.check_call(method_type, args, arg_kinds,
                               context, arg_messages=local_errors,
                               callable_name=callable_name, object_type=base_type)

    def check_op_reversible(self,
                            op_name: str,
                            left_type: Type,
                            left_expr: Expression,
                            right_type: Type,
                            right_expr: Expression,
                            context: Context,
                            msg: MessageBuilder) -> Tuple[Type, Type]:
        def lookup_operator(op_name: str, base_type: Type) -> Optional[Type]:
            """Looks up the given operator and returns the corresponding type,
            if it exists."""

            # This check is an important performance optimization,
            # even though it is mostly a subset of
            # analyze_member_access.
            # TODO: Find a way to remove this call without performance implications.
            if not self.has_member(base_type, op_name):
                return None

            local_errors = msg.clean_copy()

            member = analyze_member_access(
                name=op_name,
                typ=base_type,
                is_lvalue=False,
                is_super=False,
                is_operator=True,
                original_type=base_type,
                context=context,
                msg=local_errors,
                chk=self.chk,
                in_literal_context=self.is_literal_context()
            )
            if local_errors.is_errors():
                return None
            else:
                return member

        def lookup_definer(typ: Instance, attr_name: str) -> Optional[str]:
            """Returns the name of the class that contains the actual definition of attr_name.

            So if class A defines foo and class B subclasses A, running
            'get_class_defined_in(B, "foo")` would return the full name of A.

            However, if B were to override and redefine foo, that method call would
            return the full name of B instead.

            If the attr name is not present in the given class or its MRO, returns None.
            """
            for cls in typ.type.mro:
                if cls.names.get(attr_name):
                    return cls.fullname
            return None

        left_type = get_proper_type(left_type)
        right_type = get_proper_type(right_type)

        # If either the LHS or the RHS are Any, we can't really concluding anything
        # about the operation since the Any type may or may not define an
        # __op__ or __rop__ method. So, we punt and return Any instead.

        if isinstance(left_type, AnyType):
            any_type = AnyType(TypeOfAny.from_another_any, source_any=left_type)
            return any_type, any_type
        if isinstance(right_type, AnyType):
            any_type = AnyType(TypeOfAny.from_another_any, source_any=right_type)
            return any_type, any_type

        # STEP 1:
        # We start by getting the __op__ and __rop__ methods, if they exist.

        rev_op_name = self.get_reverse_op_method(op_name)

        left_op = lookup_operator(op_name, left_type)
        right_op = lookup_operator(rev_op_name, right_type)

        # STEP 2a:
        # We figure out in which order Python will call the operator methods. As it
        # turns out, it's not as simple as just trying to call __op__ first and
        # __rop__ second.
        #
        # We store the determined order inside the 'variants_raw' variable,
        # which records tuples containing the method, base type, and the argument.

        bias_right = is_proper_subtype(right_type, left_type)
        if op_name in operators.op_methods_that_shortcut and is_same_type(left_type, right_type):
            # When we do "A() + A()", for example, Python will only call the __add__ method,
            # never the __radd__ method.
            #
            # This is the case even if the __add__ method is completely missing and the __radd__
            # method is defined.

            variants_raw = [
                (left_op, left_type, right_expr)
            ]
        elif (is_subtype(right_type, left_type)
                and isinstance(left_type, Instance)
                and isinstance(right_type, Instance)
                and lookup_definer(left_type, op_name) != lookup_definer(right_type, rev_op_name)):
            # When we do "A() + B()" where B is a subclass of B, we'll actually try calling
            # B's __radd__ method first, but ONLY if B explicitly defines or overrides the
            # __radd__ method.
            #
            # This mechanism lets subclasses "refine" the expected outcome of the operation, even
            # if they're located on the RHS.

            variants_raw = [
                (right_op, right_type, left_expr),
                (left_op, left_type, right_expr),
            ]
        else:
            # In all other cases, we do the usual thing and call __add__ first and
            # __radd__ second when doing "A() + B()".

            variants_raw = [
                (left_op, left_type, right_expr),
                (right_op, right_type, left_expr),
            ]

        # STEP 2b:
        # When running Python 2, we might also try calling the __cmp__ method.

        is_python_2 = self.chk.options.python_version[0] == 2
        if is_python_2 and op_name in operators.ops_falling_back_to_cmp:
            cmp_method = operators.comparison_fallback_method
            left_cmp_op = lookup_operator(cmp_method, left_type)
            right_cmp_op = lookup_operator(cmp_method, right_type)

            if bias_right:
                variants_raw.append((right_cmp_op, right_type, left_expr))
                variants_raw.append((left_cmp_op, left_type, right_expr))
            else:
                variants_raw.append((left_cmp_op, left_type, right_expr))
                variants_raw.append((right_cmp_op, right_type, left_expr))

        # STEP 3:
        # We now filter out all non-existent operators. The 'variants' list contains
        # all operator methods that are actually present, in the order that Python
        # attempts to invoke them.

        variants = [(op, obj, arg) for (op, obj, arg) in variants_raw if op is not None]

        # STEP 4:
        # We now try invoking each one. If an operation succeeds, end early and return
        # the corresponding result. Otherwise, return the result and errors associated
        # with the first entry.

        errors = []
        results = []
        for method, obj, arg in variants:
            local_errors = msg.clean_copy()
            result = self.check_method_call(
                op_name, obj, method, [arg], [ARG_POS], context, local_errors)
            if local_errors.is_errors():
                errors.append(local_errors)
                results.append(result)
            else:
                return result

        # We finish invoking above operators and no early return happens. Therefore,
        # we check if either the LHS or the RHS is Instance and fallbacks to Any,
        # if so, we also return Any
        if ((isinstance(left_type, Instance) and left_type.type.fallback_to_any) or
                (isinstance(right_type, Instance) and right_type.type.fallback_to_any)):
            any_type = AnyType(TypeOfAny.special_form)
            return any_type, any_type

        # STEP 4b:
        # Sometimes, the variants list is empty. In that case, we fall-back to attempting to
        # call the __op__ method (even though it's missing).

        if not variants:
            local_errors = msg.clean_copy()
            result = self.check_method_call_by_name(
                op_name, left_type, [right_expr], [ARG_POS], context, local_errors)

            if local_errors.is_errors():
                errors.append(local_errors)
                results.append(result)
            else:
                # In theory, we should never enter this case, but it seems
                # we sometimes do, when dealing with Type[...]? E.g. see
                # check-classes.testTypeTypeComparisonWorks.
                #
                # This is probably related to the TODO in lookup_operator(...)
                # up above.
                #
                # TODO: Remove this extra case
                return result

        msg.add_errors(errors[0])
        if len(results) == 1:
            return results[0]
        else:
            error_any = AnyType(TypeOfAny.from_error)
            result = error_any, error_any
            return result

    def check_op(self, method: str, base_type: Type,
                 arg: Expression, context: Context,
                 allow_reverse: bool = False) -> Tuple[Type, Type]:
        """Type check a binary operation which maps to a method call.

        Return tuple (result type, inferred operator method type).
        """

        if allow_reverse:
            left_variants = [base_type]
            base_type = get_proper_type(base_type)
            if isinstance(base_type, UnionType):
                left_variants = [item for item in
                                 flatten_nested_unions(base_type.relevant_items(),
                                                       handle_type_alias_type=True)]
            right_type = self.accept(arg)

            # Step 1: We first try leaving the right arguments alone and destructure
            # just the left ones. (Mypy can sometimes perform some more precise inference
            # if we leave the right operands a union -- see testOperatorWithEmptyListAndSum.)
            msg = self.msg.clean_copy()
            all_results = []
            all_inferred = []

            for left_possible_type in left_variants:
                result, inferred = self.check_op_reversible(
                    op_name=method,
                    left_type=left_possible_type,
                    left_expr=TempNode(left_possible_type, context=context),
                    right_type=right_type,
                    right_expr=arg,
                    context=context,
                    msg=msg)
                all_results.append(result)
                all_inferred.append(inferred)

            if not msg.is_errors():
                results_final = make_simplified_union(all_results)
                inferred_final = make_simplified_union(all_inferred)
                return results_final, inferred_final

            # Step 2: If that fails, we try again but also destructure the right argument.
            # This is also necessary to make certain edge cases work -- see
            # testOperatorDoubleUnionInterwovenUnionAdd, for example.

            # Note: We want to pass in the original 'arg' for 'left_expr' and 'right_expr'
            # whenever possible so that plugins and similar things can introspect on the original
            # node if possible.
            #
            # We don't do the same for the base expression because it could lead to weird
            # type inference errors -- e.g. see 'testOperatorDoubleUnionSum'.
            # TODO: Can we use `type_overrides_set()` here?
            right_variants = [(right_type, arg)]
            right_type = get_proper_type(right_type)
            if isinstance(right_type, UnionType):
                right_variants = [
                    (item, TempNode(item, context=context))
                    for item in flatten_nested_unions(right_type.relevant_items(),
                                                      handle_type_alias_type=True)
                ]

            msg = self.msg.clean_copy()
            all_results = []
            all_inferred = []

            for left_possible_type in left_variants:
                for right_possible_type, right_expr in right_variants:
                    result, inferred = self.check_op_reversible(
                        op_name=method,
                        left_type=left_possible_type,
                        left_expr=TempNode(left_possible_type, context=context),
                        right_type=right_possible_type,
                        right_expr=right_expr,
                        context=context,
                        msg=msg)
                    all_results.append(result)
                    all_inferred.append(inferred)

            if msg.is_errors():
                self.msg.add_errors(msg)
                # Point any notes to the same location as an existing message.
                recent_context = msg.most_recent_context()
                if len(left_variants) >= 2 and len(right_variants) >= 2:
                    self.msg.warn_both_operands_are_from_unions(recent_context)
                elif len(left_variants) >= 2:
                    self.msg.warn_operand_was_from_union(
                        "Left", base_type, context=recent_context)
                elif len(right_variants) >= 2:
                    self.msg.warn_operand_was_from_union(
                        "Right", right_type, context=recent_context)

            # See the comment in 'check_overload_call' for more details on why
            # we call 'combine_function_signature' instead of just unioning the inferred
            # callable types.
            results_final = make_simplified_union(all_results)
            inferred_final = self.combine_function_signatures(all_inferred)
            return results_final, inferred_final
        else:
            return self.check_method_call_by_name(
                method=method,
                base_type=base_type,
                args=[arg],
                arg_kinds=[ARG_POS],
                context=context,
                local_errors=self.msg,
            )

    def get_reverse_op_method(self, method: str) -> str:
        if method == '__div__' and self.chk.options.python_version[0] == 2:
            return '__rdiv__'
        else:
            return operators.reverse_op_methods[method]

    def check_boolean_op(self, e: OpExpr, context: Context) -> Type:
        """Type check a boolean operation ('and' or 'or')."""

        # A boolean operation can evaluate to either of the operands.

        # We use the current type context to guide the type inference of of
        # the left operand. We also use the left operand type to guide the type
        # inference of the right operand so that expressions such as
        # '[1] or []' are inferred correctly.
        ctx = self.type_context[-1]
        left_type = self.accept(e.left, ctx)
        expanded_left_type = try_expanding_sum_type_to_union(
            self.accept(e.left, ctx), "builtins.bool"
        )

        assert e.op in ('and', 'or')  # Checked by visit_op_expr

        if e.right_always:
            left_map, right_map = None, {}  # type: mypy.checker.TypeMap, mypy.checker.TypeMap
        elif e.right_unreachable:
            left_map, right_map = {}, None
        elif e.op == 'and':
            right_map, left_map = self.chk.find_isinstance_check(e.left)
        elif e.op == 'or':
            left_map, right_map = self.chk.find_isinstance_check(e.left)

        # If left_map is None then we know mypy considers the left expression
        # to be redundant.
        if (
            codes.REDUNDANT_EXPR in self.chk.options.enabled_error_codes
            and left_map is None
            # don't report an error if it's intentional
            and not e.right_always
        ):
            self.msg.redundant_left_operand(e.op, e.left)

        if (
            self.chk.should_report_unreachable_issues()
            and right_map is None
            # don't report an error if it's intentional
            and not e.right_unreachable
        ):
            self.msg.unreachable_right_operand(e.op, e.right)

        # If right_map is None then we know mypy considers the right branch
        # to be unreachable and therefore any errors found in the right branch
        # should be suppressed.
        with (self.msg.disable_errors() if right_map is None else nullcontext()):
            right_type = self.analyze_cond_branch(right_map, e.right, expanded_left_type)

        if left_map is None and right_map is None:
            return UninhabitedType()

        if right_map is None:
            # The boolean expression is statically known to be the left value
            assert left_map is not None
            return left_type
        if left_map is None:
            # The boolean expression is statically known to be the right value
            assert right_map is not None
            return right_type

        if e.op == 'and':
            restricted_left_type = false_only(expanded_left_type)
            result_is_left = not expanded_left_type.can_be_true
        elif e.op == 'or':
            restricted_left_type = true_only(expanded_left_type)
            result_is_left = not expanded_left_type.can_be_false

        if isinstance(restricted_left_type, UninhabitedType):
            # The left operand can never be the result
            return right_type
        elif result_is_left:
            # The left operand is always the result
            return left_type
        else:
            return make_simplified_union([restricted_left_type, right_type])

    def check_list_multiply(self, e: OpExpr) -> Type:
        """Type check an expression of form '[...] * e'.

        Type inference is special-cased for this common construct.
        """
        right_type = self.accept(e.right)
        if is_subtype(right_type, self.named_type('builtins.int')):
            # Special case: [...] * <int value>. Use the type context of the
            # OpExpr, since the multiplication does not affect the type.
            left_type = self.accept(e.left, type_context=self.type_context[-1])
        else:
            left_type = self.accept(e.left)
        result, method_type = self.check_op('__mul__', left_type, e.right, e)
        e.method_type = method_type
        return result

    def visit_assignment_expr(self, e: AssignmentExpr) -> Type:
        value = self.accept(e.value)
        self.chk.check_assignment(e.target, e.value)
        self.chk.check_final(e)
        self.chk.store_type(e.target, value)
        self.find_partial_type_ref_fast_path(e.target)
        return value

    def visit_unary_expr(self, e: UnaryExpr) -> Type:
        """Type check an unary operation ('not', '-', '+' or '~')."""
        operand_type = self.accept(e.expr)
        op = e.op
        if op == 'not':
            result: Type = self.bool_type()
        else:
            method = operators.unary_op_methods[op]
            result, method_type = self.check_method_call_by_name(method, operand_type, [], [], e)
            e.method_type = method_type
        return result

    def visit_index_expr(self, e: IndexExpr) -> Type:
        """Type check an index expression (base[index]).

        It may also represent type application.
        """
        result = self.visit_index_expr_helper(e)
        result = get_proper_type(self.narrow_type_from_binder(e, result))
        if (self.is_literal_context() and isinstance(result, Instance)
                and result.last_known_value is not None):
            result = result.last_known_value
        return result

    def visit_index_expr_helper(self, e: IndexExpr) -> Type:
        if e.analyzed:
            # It's actually a type application.
            return self.accept(e.analyzed)
        left_type = self.accept(e.base)
        return self.visit_index_with_type(left_type, e)

    def visit_index_with_type(self, left_type: Type, e: IndexExpr,
                              original_type: Optional[ProperType] = None) -> Type:
        """Analyze type of an index expression for a given type of base expression.

        The 'original_type' is used for error messages (currently used for union types).
        """
        index = e.index
        left_type = get_proper_type(left_type)

        # Visit the index, just to make sure we have a type for it available
        self.accept(index)

        if isinstance(left_type, UnionType):
            original_type = original_type or left_type
            # Don't combine literal types, since we may need them for type narrowing.
            return make_simplified_union([self.visit_index_with_type(typ, e,
                                                                     original_type)
                                          for typ in left_type.relevant_items()],
                                         contract_literals=False)
        elif isinstance(left_type, TupleType) and self.chk.in_checked_function():
            # Special case for tuples. They return a more specific type when
            # indexed by an integer literal.
            if isinstance(index, SliceExpr):
                return self.visit_tuple_slice_helper(left_type, index)

            ns = self.try_getting_int_literals(index)
            if ns is not None:
                out = []
                for n in ns:
                    if n < 0:
                        n += len(left_type.items)
                    if 0 <= n < len(left_type.items):
                        out.append(left_type.items[n])
                    else:
                        self.chk.fail(message_registry.TUPLE_INDEX_OUT_OF_RANGE, e)
                        return AnyType(TypeOfAny.from_error)
                return make_simplified_union(out)
            else:
                return self.nonliteral_tuple_index_helper(left_type, index)
        elif isinstance(left_type, TypedDictType):
            return self.visit_typeddict_index_expr(left_type, e.index)
        elif (isinstance(left_type, CallableType)
              and left_type.is_type_obj() and left_type.type_object().is_enum):
            return self.visit_enum_index_expr(left_type.type_object(), e.index, e)
        elif (isinstance(left_type, TypeVarType)
              and not self.has_member(left_type.upper_bound, "__getitem__")):
            return self.visit_index_with_type(left_type.upper_bound, e, original_type)
        else:
            result, method_type = self.check_method_call_by_name(
                '__getitem__', left_type, [e.index], [ARG_POS], e,
                original_type=original_type)
            e.method_type = method_type
            return result

    def visit_tuple_slice_helper(self, left_type: TupleType, slic: SliceExpr) -> Type:
        begin: Sequence[Optional[int]] = [None]
        end: Sequence[Optional[int]] = [None]
        stride: Sequence[Optional[int]] = [None]

        if slic.begin_index:
            begin_raw = self.try_getting_int_literals(slic.begin_index)
            if begin_raw is None:
                return self.nonliteral_tuple_index_helper(left_type, slic)
            begin = begin_raw

        if slic.end_index:
            end_raw = self.try_getting_int_literals(slic.end_index)
            if end_raw is None:
                return self.nonliteral_tuple_index_helper(left_type, slic)
            end = end_raw

        if slic.stride:
            stride_raw = self.try_getting_int_literals(slic.stride)
            if stride_raw is None:
                return self.nonliteral_tuple_index_helper(left_type, slic)
            stride = stride_raw

        items: List[Type] = []
        for b, e, s in itertools.product(begin, end, stride):
            items.append(left_type.slice(b, e, s))
        return make_simplified_union(items)

    def try_getting_int_literals(self, index: Expression) -> Optional[List[int]]:
        """If the given expression or type corresponds to an int literal
        or a union of int literals, returns a list of the underlying ints.
        Otherwise, returns None.

        Specifically, this function is guaranteed to return a list with
        one or more ints if one one the following is true:

        1. 'expr' is a IntExpr or a UnaryExpr backed by an IntExpr
        2. 'typ' is a LiteralType containing an int
        3. 'typ' is a UnionType containing only LiteralType of ints
        """
        if isinstance(index, IntExpr):
            return [index.value]
        elif isinstance(index, UnaryExpr):
            if index.op == '-':
                operand = index.expr
                if isinstance(operand, IntExpr):
                    return [-1 * operand.value]
        typ = get_proper_type(self.accept(index))
        if isinstance(typ, Instance) and typ.last_known_value is not None:
            typ = typ.last_known_value
        if isinstance(typ, LiteralType) and isinstance(typ.value, int):
            return [typ.value]
        if isinstance(typ, UnionType):
            out = []
            for item in get_proper_types(typ.items):
                if isinstance(item, LiteralType) and isinstance(item.value, int):
                    out.append(item.value)
                else:
                    return None
            return out
        return None

    def nonliteral_tuple_index_helper(self, left_type: TupleType, index: Expression) -> Type:
        index_type = self.accept(index)
        expected_type = UnionType.make_union([self.named_type('builtins.int'),
                                              self.named_type('builtins.slice')])
        if not self.chk.check_subtype(index_type, expected_type, index,
                                      message_registry.INVALID_TUPLE_INDEX_TYPE,
                                      'actual type', 'expected type'):
            return AnyType(TypeOfAny.from_error)
        else:
            union = make_simplified_union(left_type.items)
            if isinstance(index, SliceExpr):
                return self.chk.named_generic_type('builtins.tuple', [union])
            else:
                return union

    def visit_typeddict_index_expr(self, td_type: TypedDictType, index: Expression) -> Type:
        if isinstance(index, (StrExpr, UnicodeExpr)):
            key_names = [index.value]
        else:
            typ = get_proper_type(self.accept(index))
            if isinstance(typ, UnionType):
                key_types: List[Type] = list(typ.items)
            else:
                key_types = [typ]

            key_names = []
            for key_type in get_proper_types(key_types):
                if isinstance(key_type, Instance) and key_type.last_known_value is not None:
                    key_type = key_type.last_known_value

                if (isinstance(key_type, LiteralType)
                        and isinstance(key_type.value, str)
                        and key_type.fallback.type.fullname != 'builtins.bytes'):
                    key_names.append(key_type.value)
                else:
                    self.msg.typeddict_key_must_be_string_literal(td_type, index)
                    return AnyType(TypeOfAny.from_error)

        value_types = []
        for key_name in key_names:
            value_type = td_type.items.get(key_name)
            if value_type is None:
                self.msg.typeddict_key_not_found(td_type, key_name, index)
                return AnyType(TypeOfAny.from_error)
            else:
                value_types.append(value_type)
        return make_simplified_union(value_types)

    def visit_enum_index_expr(
        self, enum_type: TypeInfo, index: Expression, context: Context
    ) -> Type:
        string_type: Type = self.named_type("builtins.str")
        if self.chk.options.python_version[0] < 3:
            string_type = UnionType.make_union([string_type,
                                                self.named_type('builtins.unicode')])
        self.chk.check_subtype(self.accept(index), string_type, context,
                               "Enum index should be a string", "actual index type")
        return Instance(enum_type, [])

    def visit_cast_expr(self, expr: CastExpr) -> Type:
        """Type check a cast expression."""
        source_type = self.accept(expr.expr, type_context=AnyType(TypeOfAny.special_form),
                                  allow_none_return=True, always_allow_any=True)
        target_type = expr.type
        options = self.chk.options
        if (options.warn_redundant_casts and not isinstance(get_proper_type(target_type), AnyType)
                and is_same_type(source_type, target_type)):
            self.msg.redundant_cast(target_type, expr)
        if options.disallow_any_unimported and has_any_from_unimported_type(target_type):
            self.msg.unimported_type_becomes_any("Target type of cast", target_type, expr)
        check_for_explicit_any(target_type, self.chk.options, self.chk.is_typeshed_stub, self.msg,
                               context=expr)
        return target_type

    def visit_reveal_expr(self, expr: RevealExpr) -> Type:
        """Type check a reveal_type expression."""
        if expr.kind == REVEAL_TYPE:
            assert expr.expr is not None
            revealed_type = self.accept(expr.expr, type_context=self.type_context[-1],
                                        allow_none_return=True)
            if not self.chk.current_node_deferred:
                self.msg.reveal_type(revealed_type, expr.expr)
                if not self.chk.in_checked_function():
                    self.msg.note("'reveal_type' always outputs 'Any' in unchecked functions",
                                  expr.expr)
            return revealed_type
        else:
            # REVEAL_LOCALS
            if not self.chk.current_node_deferred:
                # the RevealExpr contains a local_nodes attribute,
                # calculated at semantic analysis time. Use it to pull out the
                # corresponding subset of variables in self.chk.type_map
                names_to_types = {
                    var_node.name: var_node.type for var_node in expr.local_nodes
                } if expr.local_nodes is not None else {}

                self.msg.reveal_locals(names_to_types, expr)
            return NoneType()

    def visit_type_application(self, tapp: TypeApplication) -> Type:
        """Type check a type application (expr[type, ...]).

        There are two different options here, depending on whether expr refers
        to a type alias or directly to a generic class. In the first case we need
        to use a dedicated function typeanal.expand_type_aliases. This
        is due to the fact that currently type aliases machinery uses
        unbound type variables, while normal generics use bound ones;
        see TypeAlias docstring for more details.
        """
        if isinstance(tapp.expr, RefExpr) and isinstance(tapp.expr.node, TypeAlias):
            # Subscription of a (generic) alias in runtime context, expand the alias.
            item = expand_type_alias(tapp.expr.node, tapp.types, self.chk.fail,
                                     tapp.expr.node.no_args, tapp)
            item = get_proper_type(item)
            if isinstance(item, Instance):
                tp = type_object_type(item.type, self.named_type)
                return self.apply_type_arguments_to_callable(tp, item.args, tapp)
            else:
                self.chk.fail(message_registry.ONLY_CLASS_APPLICATION, tapp)
                return AnyType(TypeOfAny.from_error)
        # Type application of a normal generic class in runtime context.
        # This is typically used as `x = G[int]()`.
        tp = get_proper_type(self.accept(tapp.expr))
        if isinstance(tp, (CallableType, Overloaded)):
            if not tp.is_type_obj():
                self.chk.fail(message_registry.ONLY_CLASS_APPLICATION, tapp)
            return self.apply_type_arguments_to_callable(tp, tapp.types, tapp)
        if isinstance(tp, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=tp)
        return AnyType(TypeOfAny.special_form)

    def visit_type_alias_expr(self, alias: TypeAliasExpr) -> Type:
        """Right hand side of a type alias definition.

        It has the same type as if the alias itself was used in a runtime context.
        For example, here:

            A = reveal_type(List[T])
            reveal_type(A)

        both `reveal_type` instances will reveal the same type `def (...) -> builtins.list[Any]`.
        Note that type variables are implicitly substituted with `Any`.
        """
        return self.alias_type_in_runtime_context(alias.node, alias.no_args,
                                                  alias, alias_definition=True)

    def alias_type_in_runtime_context(self, alias: TypeAlias,
                                      no_args: bool, ctx: Context,
                                      *,
                                      alias_definition: bool = False) -> Type:
        """Get type of a type alias (could be generic) in a runtime expression.

        Note that this function can be called only if the alias appears _not_
        as a target of type application, which is treated separately in the
        visit_type_application method. Some examples where this method is called are
        casts and instantiation:

            class LongName(Generic[T]): ...
            A = LongName[int]

            x = A()
            y = cast(A, ...)
        """
        if isinstance(alias.target, Instance) and alias.target.invalid:  # type: ignore
            # An invalid alias, error already has been reported
            return AnyType(TypeOfAny.from_error)
        # If this is a generic alias, we set all variables to `Any`.
        # For example:
        #     A = List[Tuple[T, T]]
        #     x = A() <- same as List[Tuple[Any, Any]], see PEP 484.
        item = get_proper_type(set_any_tvars(alias, ctx.line, ctx.column))
        if isinstance(item, Instance):
            # Normally we get a callable type (or overloaded) with .is_type_obj() true
            # representing the class's constructor
            tp = type_object_type(item.type, self.named_type)
            if no_args:
                return tp
            return self.apply_type_arguments_to_callable(tp, item.args, ctx)
        elif (isinstance(item, TupleType) and
              # Tuple[str, int]() fails at runtime, only named tuples and subclasses work.
              tuple_fallback(item).type.fullname != 'builtins.tuple'):
            return type_object_type(tuple_fallback(item).type, self.named_type)
        elif isinstance(item, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=item)
        else:
            if alias_definition:
                return AnyType(TypeOfAny.special_form)
            # This type is invalid in most runtime contexts, give it an 'object' type.
            return self.named_type('builtins.object')

    def apply_type_arguments_to_callable(
        self, tp: Type, args: Sequence[Type], ctx: Context
    ) -> Type:
        """Apply type arguments to a generic callable type coming from a type object.

        This will first perform type arguments count checks, report the
        error as needed, and return the correct kind of Any. As a special
        case this returns Any for non-callable types, because if type object type
        is not callable, then an error should be already reported.
        """
        tp = get_proper_type(tp)

        if isinstance(tp, CallableType):
            if len(tp.variables) != len(args):
                self.msg.incompatible_type_application(len(tp.variables),
                                                       len(args), ctx)
                return AnyType(TypeOfAny.from_error)
            return self.apply_generic_arguments(tp, args, ctx)
        if isinstance(tp, Overloaded):
            for it in tp.items:
                if len(it.variables) != len(args):
                    self.msg.incompatible_type_application(len(it.variables),
                                                           len(args), ctx)
                    return AnyType(TypeOfAny.from_error)
            return Overloaded([self.apply_generic_arguments(it, args, ctx)
                               for it in tp.items])
        return AnyType(TypeOfAny.special_form)

    def visit_list_expr(self, e: ListExpr) -> Type:
        """Type check a list expression [...]."""
        return self.check_lst_expr(e.items, 'builtins.list', '<list>', e)

    def visit_set_expr(self, e: SetExpr) -> Type:
        return self.check_lst_expr(e.items, 'builtins.set', '<set>', e)

    def fast_container_type(
            self, items: List[Expression], container_fullname: str
    ) -> Optional[Type]:
        """
        Fast path to determine the type of a list or set literal,
        based on the list of entries. This mostly impacts large
        module-level constant definitions.

        Limitations:
         - no active type context
         - no star expressions
         - the joined type of all entries must be an Instance type
        """
        ctx = self.type_context[-1]
        if ctx:
            return None
        values: List[Type] = []
        for item in items:
            if isinstance(item, StarExpr):
                # fallback to slow path
                return None
            values.append(self.accept(item))
        vt = join.join_type_list(values)
        if not isinstance(vt, Instance):
            return None
        # TODO: update tests instead?
        vt.erased = True
        return self.chk.named_generic_type(container_fullname, [vt])

    def check_lst_expr(self, items: List[Expression], fullname: str,
                       tag: str, context: Context) -> Type:
        # fast path
        t = self.fast_container_type(items, fullname)
        if t:
            return t

        # Translate into type checking a generic function call.
        # Used for list and set expressions, as well as for tuples
        # containing star expressions that don't refer to a
        # Tuple. (Note: "lst" stands for list-set-tuple. :-)
        tv = TypeVarType('T', 'T', -1, [], self.object_type())
        constructor = CallableType(
            [tv],
            [nodes.ARG_STAR],
            [None],
            self.chk.named_generic_type(fullname, [tv]),
            self.named_type('builtins.function'),
            name=tag,
            variables=[tv])
        out = self.check_call(constructor,
                              [(i.expr if isinstance(i, StarExpr) else i)
                               for i in items],
                              [(nodes.ARG_STAR if isinstance(i, StarExpr) else nodes.ARG_POS)
                               for i in items],
                              context)[0]
        return remove_instance_last_known_values(out)

    def visit_tuple_expr(self, e: TupleExpr) -> Type:
        """Type check a tuple expression."""
        # Try to determine type context for type inference.
        type_context = get_proper_type(self.type_context[-1])
        type_context_items = None
        if isinstance(type_context, UnionType):
            tuples_in_context = [t for t in get_proper_types(type_context.items)
                                 if (isinstance(t, TupleType) and len(t.items) == len(e.items)) or
                                 is_named_instance(t, 'builtins.tuple')]
            if len(tuples_in_context) == 1:
                type_context = tuples_in_context[0]
            else:
                # There are either no relevant tuples in the Union, or there is
                # more than one.  Either way, we can't decide on a context.
                pass

        if isinstance(type_context, TupleType):
            type_context_items = type_context.items
        elif type_context and is_named_instance(type_context, 'builtins.tuple'):
            assert isinstance(type_context, Instance)
            if type_context.args:
                type_context_items = [type_context.args[0]] * len(e.items)
        # NOTE: it's possible for the context to have a different
        # number of items than e.  In that case we use those context
        # items that match a position in e, and we'll worry about type
        # mismatches later.

        # Infer item types.  Give up if there's a star expression
        # that's not a Tuple.
        items: List[Type] = []
        j = 0  # Index into type_context_items; irrelevant if type_context_items is none
        for i in range(len(e.items)):
            item = e.items[i]
            if isinstance(item, StarExpr):
                # Special handling for star expressions.
                # TODO: If there's a context, and item.expr is a
                # TupleExpr, flatten it, so we can benefit from the
                # context?  Counterargument: Why would anyone write
                # (1, *(2, 3)) instead of (1, 2, 3) except in a test?
                tt = self.accept(item.expr)
                tt = get_proper_type(tt)
                if isinstance(tt, TupleType):
                    items.extend(tt.items)
                    j += len(tt.items)
                else:
                    # A star expression that's not a Tuple.
                    # Treat the whole thing as a variable-length tuple.
                    return self.check_lst_expr(e.items, 'builtins.tuple', '<tuple>', e)
            else:
                if not type_context_items or j >= len(type_context_items):
                    tt = self.accept(item)
                else:
                    tt = self.accept(item, type_context_items[j])
                    j += 1
                items.append(tt)
        # This is a partial fallback item type. A precise type will be calculated on demand.
        fallback_item = AnyType(TypeOfAny.special_form)
        return TupleType(items, self.chk.named_generic_type('builtins.tuple', [fallback_item]))

    def fast_dict_type(self, e: DictExpr) -> Optional[Type]:
        """
        Fast path to determine the type of a dict literal,
        based on the list of entries. This mostly impacts large
        module-level constant definitions.

        Limitations:
         - no active type context
         - only supported star expressions are other dict instances
         - the joined types of all keys and values must be Instance types
        """
        ctx = self.type_context[-1]
        if ctx:
            return None
        keys: List[Type] = []
        values: List[Type] = []
        stargs: Optional[Tuple[Type, Type]] = None
        for key, value in e.items:
            if key is None:
                st = get_proper_type(self.accept(value))
                if (
                        isinstance(st, Instance)
                        and st.type.fullname == 'builtins.dict'
                        and len(st.args) == 2
                ):
                    stargs = (st.args[0], st.args[1])
                else:
                    return None
            else:
                keys.append(self.accept(key))
                values.append(self.accept(value))
        kt = join.join_type_list(keys)
        vt = join.join_type_list(values)
        if not (isinstance(kt, Instance) and isinstance(vt, Instance)):
            return None
        if stargs and (stargs[0] != kt or stargs[1] != vt):
            return None
        # TODO: update tests instead?
        kt.erased = True
        vt.erased = True
        return self.chk.named_generic_type('builtins.dict', [kt, vt])

    def visit_dict_expr(self, e: DictExpr) -> Type:
        """Type check a dict expression.

        Translate it into a call to dict(), with provisions for **expr.
        """
        # if the dict literal doesn't match TypedDict, check_typeddict_call_with_dict reports
        # an error, but returns the TypedDict type that matches the literal it found
        # that would cause a second error when that TypedDict type is returned upstream
        # to avoid the second error, we always return TypedDict type that was requested
        typeddict_context = self.find_typeddict_context(self.type_context[-1], e)
        if typeddict_context:
            self.check_typeddict_call_with_dict(
                callee=typeddict_context,
                kwargs=e,
                context=e
            )
            return typeddict_context.copy_modified()

        # fast path attempt
        dt = self.fast_dict_type(e)
        if dt:
            return dt

        # Collect function arguments, watching out for **expr.
        args: List[Expression] = []  # Regular "key: value"
        stargs: List[Expression] = []  # For "**expr"
        for key, value in e.items:
            if key is None:
                stargs.append(value)
            else:
                tup = TupleExpr([key, value])
                if key.line >= 0:
                    tup.line = key.line
                    tup.column = key.column
                else:
                    tup.line = value.line
                    tup.column = value.column
                args.append(tup)
        # Define type variables (used in constructors below).
        kt = TypeVarType('KT', 'KT', -1, [], self.object_type())
        vt = TypeVarType('VT', 'VT', -2, [], self.object_type())
        rv = None
        # Call dict(*args), unless it's empty and stargs is not.
        if args or not stargs:
            # The callable type represents a function like this:
            #
            #   def <unnamed>(*v: Tuple[kt, vt]) -> Dict[kt, vt]: ...
            constructor = CallableType(
                [TupleType([kt, vt], self.named_type('builtins.tuple'))],
                [nodes.ARG_STAR],
                [None],
                self.chk.named_generic_type('builtins.dict', [kt, vt]),
                self.named_type('builtins.function'),
                name='<dict>',
                variables=[kt, vt])
            rv = self.check_call(constructor, args, [nodes.ARG_POS] * len(args), e)[0]
        else:
            # dict(...) will be called below.
            pass
        # Call rv.update(arg) for each arg in **stargs,
        # except if rv isn't set yet, then set rv = dict(arg).
        if stargs:
            for arg in stargs:
                if rv is None:
                    constructor = CallableType(
                        [self.chk.named_generic_type('typing.Mapping', [kt, vt])],
                        [nodes.ARG_POS],
                        [None],
                        self.chk.named_generic_type('builtins.dict', [kt, vt]),
                        self.named_type('builtins.function'),
                        name='<list>',
                        variables=[kt, vt])
                    rv = self.check_call(constructor, [arg], [nodes.ARG_POS], arg)[0]
                else:
                    self.check_method_call_by_name('update', rv, [arg], [nodes.ARG_POS], arg)
        assert rv is not None
        return rv

    def find_typeddict_context(self, context: Optional[Type],
                               dict_expr: DictExpr) -> Optional[TypedDictType]:
        context = get_proper_type(context)
        if isinstance(context, TypedDictType):
            return context
        elif isinstance(context, UnionType):
            items = []
            for item in context.items:
                item_context = self.find_typeddict_context(item, dict_expr)
                if (item_context is not None
                        and self.match_typeddict_call_with_dict(
                            item_context, dict_expr, dict_expr)):
                    items.append(item_context)
            if len(items) == 1:
                # Only one union item is valid TypedDict for the given dict_expr, so use the
                # context as it's unambiguous.
                return items[0]
            if len(items) > 1:
                self.msg.typeddict_context_ambiguous(items, dict_expr)
        # No TypedDict type in context.
        return None

    def visit_lambda_expr(self, e: LambdaExpr) -> Type:
        """Type check lambda expression."""
        self.chk.check_default_args(e, body_is_trivial=False)
        inferred_type, type_override = self.infer_lambda_type_using_context(e)
        if not inferred_type:
            self.chk.return_types.append(AnyType(TypeOfAny.special_form))
            # Type check everything in the body except for the final return
            # statement (it can contain tuple unpacking before return).
            with self.chk.scope.push_function(e):
                # Lambdas can have more than one element in body,
                # when we add "fictional" AssigmentStatement nodes, like in:
                # `lambda (a, b): a`
                for stmt in e.body.body[:-1]:
                    stmt.accept(self.chk)
                # Only type check the return expression, not the return statement.
                # This is important as otherwise the following statements would be
                # considered unreachable. There's no useful type context.
                ret_type = self.accept(e.expr(), allow_none_return=True)
            fallback = self.named_type('builtins.function')
            self.chk.return_types.pop()
            return callable_type(e, fallback, ret_type)
        else:
            # Type context available.
            self.chk.return_types.append(inferred_type.ret_type)
            self.chk.check_func_item(e, type_override=type_override)
            if e.expr() not in self.chk.type_map:
                # TODO: return expression must be accepted before exiting function scope.
                self.accept(e.expr(), allow_none_return=True)
            ret_type = self.chk.type_map[e.expr()]
            self.chk.return_types.pop()
            return replace_callable_return_type(inferred_type, ret_type)

    def infer_lambda_type_using_context(self, e: LambdaExpr) -> Tuple[Optional[CallableType],
                                                                    Optional[CallableType]]:
        """Try to infer lambda expression type using context.

        Return None if could not infer type.
        The second item in the return type is the type_override parameter for check_func_item.
        """
        # TODO also accept 'Any' context
        ctx = get_proper_type(self.type_context[-1])

        if isinstance(ctx, UnionType):
            callables = [t for t in get_proper_types(ctx.relevant_items())
                         if isinstance(t, CallableType)]
            if len(callables) == 1:
                ctx = callables[0]

        if not ctx or not isinstance(ctx, CallableType):
            return None, None

        # The context may have function type variables in it. We replace them
        # since these are the type variables we are ultimately trying to infer;
        # they must be considered as indeterminate. We use ErasedType since it
        # does not affect type inference results (it is for purposes like this
        # only).
        callable_ctx = get_proper_type(replace_meta_vars(ctx, ErasedType()))
        assert isinstance(callable_ctx, CallableType)

        if callable_ctx.type_guard is not None:
            # Lambda's return type cannot be treated as a `TypeGuard`,
            # because it is implicit. And `TypeGuard`s must be explicit.
            # See https://github.com/python/mypy/issues/9927
            return None, None

        arg_kinds = [arg.kind for arg in e.arguments]

        if callable_ctx.is_ellipsis_args or ctx.param_spec() is not None:
            # Fill in Any arguments to match the arguments of the lambda.
            callable_ctx = callable_ctx.copy_modified(
                is_ellipsis_args=False,
                arg_types=[AnyType(TypeOfAny.special_form)] * len(arg_kinds),
                arg_kinds=arg_kinds,
                arg_names=e.arg_names[:],
            )

        if ARG_STAR in arg_kinds or ARG_STAR2 in arg_kinds:
            # TODO treat this case appropriately
            return callable_ctx, None

        if callable_ctx.arg_kinds != arg_kinds:
            # Incompatible context; cannot use it to infer types.
            self.chk.fail(message_registry.CANNOT_INFER_LAMBDA_TYPE, e)
            return None, None

        return callable_ctx, callable_ctx

    def visit_super_expr(self, e: SuperExpr) -> Type:
        """Type check a super expression (non-lvalue)."""

        # We have an expression like super(T, var).member

        # First compute the types of T and var
        types = self._super_arg_types(e)
        if isinstance(types, tuple):
            type_type, instance_type = types
        else:
            return types

        # Now get the MRO
        type_info = type_info_from_type(type_type)
        if type_info is None:
            self.chk.fail(message_registry.UNSUPPORTED_ARG_1_FOR_SUPER, e)
            return AnyType(TypeOfAny.from_error)

        instance_info = type_info_from_type(instance_type)
        if instance_info is None:
            self.chk.fail(message_registry.UNSUPPORTED_ARG_2_FOR_SUPER, e)
            return AnyType(TypeOfAny.from_error)

        mro = instance_info.mro

        # The base is the first MRO entry *after* type_info that has a member
        # with the right name
        try:
            index = mro.index(type_info)
        except ValueError:
            self.chk.fail(message_registry.SUPER_ARG_2_NOT_INSTANCE_OF_ARG_1, e)
            return AnyType(TypeOfAny.from_error)

        if len(mro) == index + 1:
            self.chk.fail(message_registry.TARGET_CLASS_HAS_NO_BASE_CLASS, e)
            return AnyType(TypeOfAny.from_error)

        for base in mro[index+1:]:
            if e.name in base.names or base == mro[-1]:
                if e.info and e.info.fallback_to_any and base == mro[-1]:
                    # There's an undefined base class, and we're at the end of the
                    # chain.  That's not an error.
                    return AnyType(TypeOfAny.special_form)

                return analyze_member_access(name=e.name,
                                             typ=instance_type,
                                             is_lvalue=False,
                                             is_super=True,
                                             is_operator=False,
                                             original_type=instance_type,
                                             override_info=base,
                                             context=e,
                                             msg=self.msg,
                                             chk=self.chk,
                                             in_literal_context=self.is_literal_context())

        assert False, 'unreachable'

    def _super_arg_types(self, e: SuperExpr) -> Union[Type, Tuple[Type, Type]]:
        """
        Computes the types of the type and instance expressions in super(T, instance), or the
        implicit ones for zero-argument super() expressions.  Returns a single type for the whole
        super expression when possible (for errors, anys), otherwise the pair of computed types.
        """

        if not self.chk.in_checked_function():
            return AnyType(TypeOfAny.unannotated)
        elif len(e.call.args) == 0:
            if self.chk.options.python_version[0] == 2:
                self.chk.fail(message_registry.TOO_FEW_ARGS_FOR_SUPER, e)
                return AnyType(TypeOfAny.from_error)
            elif not e.info:
                # This has already been reported by the semantic analyzer.
                return AnyType(TypeOfAny.from_error)
            elif self.chk.scope.active_class():
                self.chk.fail(message_registry.SUPER_OUTSIDE_OF_METHOD_NOT_SUPPORTED, e)
                return AnyType(TypeOfAny.from_error)

            # Zero-argument super() is like super(<current class>, <self>)
            current_type = fill_typevars(e.info)
            type_type: ProperType = TypeType(current_type)

            # Use the type of the self argument, in case it was annotated
            method = self.chk.scope.top_function()
            assert method is not None
            if method.arguments:
                instance_type: Type = method.arguments[0].variable.type or current_type
            else:
                self.chk.fail(message_registry.SUPER_ENCLOSING_POSITIONAL_ARGS_REQUIRED, e)
                return AnyType(TypeOfAny.from_error)
        elif ARG_STAR in e.call.arg_kinds:
            self.chk.fail(message_registry.SUPER_VARARGS_NOT_SUPPORTED, e)
            return AnyType(TypeOfAny.from_error)
        elif set(e.call.arg_kinds) != {ARG_POS}:
            self.chk.fail(message_registry.SUPER_POSITIONAL_ARGS_REQUIRED, e)
            return AnyType(TypeOfAny.from_error)
        elif len(e.call.args) == 1:
            self.chk.fail(message_registry.SUPER_WITH_SINGLE_ARG_NOT_SUPPORTED, e)
            return AnyType(TypeOfAny.from_error)
        elif len(e.call.args) == 2:
            type_type = get_proper_type(self.accept(e.call.args[0]))
            instance_type = self.accept(e.call.args[1])
        else:
            self.chk.fail(message_registry.TOO_MANY_ARGS_FOR_SUPER, e)
            return AnyType(TypeOfAny.from_error)

        # Imprecisely assume that the type is the current class
        if isinstance(type_type, AnyType):
            if e.info:
                type_type = TypeType(fill_typevars(e.info))
            else:
                return AnyType(TypeOfAny.from_another_any, source_any=type_type)
        elif isinstance(type_type, TypeType):
            type_item = type_type.item
            if isinstance(type_item, AnyType):
                if e.info:
                    type_type = TypeType(fill_typevars(e.info))
                else:
                    return AnyType(TypeOfAny.from_another_any, source_any=type_item)

        if (not isinstance(type_type, TypeType)
                and not (isinstance(type_type, FunctionLike) and type_type.is_type_obj())):
            self.msg.first_argument_for_super_must_be_type(type_type, e)
            return AnyType(TypeOfAny.from_error)

        # Imprecisely assume that the instance is of the current class
        instance_type = get_proper_type(instance_type)
        if isinstance(instance_type, AnyType):
            if e.info:
                instance_type = fill_typevars(e.info)
            else:
                return AnyType(TypeOfAny.from_another_any, source_any=instance_type)
        elif isinstance(instance_type, TypeType):
            instance_item = instance_type.item
            if isinstance(instance_item, AnyType):
                if e.info:
                    instance_type = TypeType(fill_typevars(e.info))
                else:
                    return AnyType(TypeOfAny.from_another_any, source_any=instance_item)

        return type_type, instance_type

    def visit_slice_expr(self, e: SliceExpr) -> Type:
        expected = make_optional_type(self.named_type('builtins.int'))
        for index in [e.begin_index, e.end_index, e.stride]:
            if index:
                t = self.accept(index)
                self.chk.check_subtype(t, expected,
                                       index, message_registry.INVALID_SLICE_INDEX)
        return self.named_type('builtins.slice')

    def visit_list_comprehension(self, e: ListComprehension) -> Type:
        return self.check_generator_or_comprehension(
            e.generator, 'builtins.list', '<list-comprehension>')

    def visit_set_comprehension(self, e: SetComprehension) -> Type:
        return self.check_generator_or_comprehension(
            e.generator, 'builtins.set', '<set-comprehension>')

    def visit_generator_expr(self, e: GeneratorExpr) -> Type:
        # If any of the comprehensions use async for, the expression will return an async generator
        # object
        if any(e.is_async):
            typ = 'typing.AsyncGenerator'
            # received type is always None in async generator expressions
            additional_args: List[Type] = [NoneType()]
        else:
            typ = 'typing.Generator'
            # received type and returned type are None
            additional_args = [NoneType(), NoneType()]
        return self.check_generator_or_comprehension(e, typ, '<generator>',
                                                     additional_args=additional_args)

    def check_generator_or_comprehension(self, gen: GeneratorExpr,
                                         type_name: str,
                                         id_for_messages: str,
                                         additional_args: Optional[List[Type]] = None) -> Type:
        """Type check a generator expression or a list comprehension."""
        additional_args = additional_args or []
        with self.chk.binder.frame_context(can_skip=True, fall_through=0):
            self.check_for_comp(gen)

            # Infer the type of the list comprehension by using a synthetic generic
            # callable type.
            tv = TypeVarType('T', 'T', -1, [], self.object_type())
            tv_list: List[Type] = [tv]
            constructor = CallableType(
                tv_list,
                [nodes.ARG_POS],
                [None],
                self.chk.named_generic_type(type_name, tv_list + additional_args),
                self.chk.named_type('builtins.function'),
                name=id_for_messages,
                variables=[tv])
            return self.check_call(constructor,
                                [gen.left_expr], [nodes.ARG_POS], gen)[0]

    def visit_dictionary_comprehension(self, e: DictionaryComprehension) -> Type:
        """Type check a dictionary comprehension."""
        with self.chk.binder.frame_context(can_skip=True, fall_through=0):
            self.check_for_comp(e)

            # Infer the type of the list comprehension by using a synthetic generic
            # callable type.
            ktdef = TypeVarType('KT', 'KT', -1, [], self.object_type())
            vtdef = TypeVarType('VT', 'VT', -2, [], self.object_type())
            constructor = CallableType(
                [ktdef, vtdef],
                [nodes.ARG_POS, nodes.ARG_POS],
                [None, None],
                self.chk.named_generic_type('builtins.dict', [ktdef, vtdef]),
                self.chk.named_type('builtins.function'),
                name='<dictionary-comprehension>',
                variables=[ktdef, vtdef])
            return self.check_call(constructor,
                                   [e.key, e.value], [nodes.ARG_POS, nodes.ARG_POS], e)[0]

    def check_for_comp(self, e: Union[GeneratorExpr, DictionaryComprehension]) -> None:
        """Check the for_comp part of comprehensions. That is the part from 'for':
        ... for x in y if z

        Note: This adds the type information derived from the condlists to the current binder.
        """
        for index, sequence, conditions, is_async in zip(e.indices, e.sequences,
                                                         e.condlists, e.is_async):
            if is_async:
                _, sequence_type = self.chk.analyze_async_iterable_item_type(sequence)
            else:
                _, sequence_type = self.chk.analyze_iterable_item_type(sequence)
            self.chk.analyze_index_variables(index, sequence_type, True, e)
            for condition in conditions:
                self.accept(condition)

                # values are only part of the comprehension when all conditions are true
                true_map, false_map = self.chk.find_isinstance_check(condition)

                if true_map:
                    self.chk.push_type_map(true_map)

                if codes.REDUNDANT_EXPR in self.chk.options.enabled_error_codes:
                    if true_map is None:
                        self.msg.redundant_condition_in_comprehension(False, condition)
                    elif false_map is None:
                        self.msg.redundant_condition_in_comprehension(True, condition)

    def visit_conditional_expr(self, e: ConditionalExpr, allow_none_return: bool = False) -> Type:
        self.accept(e.cond)
        ctx = self.type_context[-1]

        # Gain type information from isinstance if it is there
        # but only for the current expression
        if_map, else_map = self.chk.find_isinstance_check(e.cond)
        if codes.REDUNDANT_EXPR in self.chk.options.enabled_error_codes:
            if if_map is None:
                self.msg.redundant_condition_in_if(False, e.cond)
            elif else_map is None:
                self.msg.redundant_condition_in_if(True, e.cond)

        if_type = self.analyze_cond_branch(if_map, e.if_expr, context=ctx,
                                           allow_none_return=allow_none_return)

        # Analyze the right branch using full type context and store the type
        full_context_else_type = self.analyze_cond_branch(else_map, e.else_expr, context=ctx,
                                                          allow_none_return=allow_none_return)
        if not mypy.checker.is_valid_inferred_type(if_type):
            # Analyze the right branch disregarding the left branch.
            else_type = full_context_else_type

            # If it would make a difference, re-analyze the left
            # branch using the right branch's type as context.
            if ctx is None or not is_equivalent(else_type, ctx):
                # TODO: If it's possible that the previous analysis of
                # the left branch produced errors that are avoided
                # using this context, suppress those errors.
                if_type = self.analyze_cond_branch(if_map, e.if_expr, context=else_type,
                                                   allow_none_return=allow_none_return)

        else:
            # Analyze the right branch in the context of the left
            # branch's type.
            else_type = self.analyze_cond_branch(else_map, e.else_expr, context=if_type,
                                                 allow_none_return=allow_none_return)

        # Only create a union type if the type context is a union, to be mostly
        # compatible with older mypy versions where we always did a join.
        #
        # TODO: Always create a union or at least in more cases?
        if isinstance(get_proper_type(self.type_context[-1]), UnionType):
            res = make_simplified_union([if_type, full_context_else_type])
        else:
            res = join.join_types(if_type, else_type)

        return res

    def analyze_cond_branch(self, map: Optional[Dict[Expression, Type]],
                            node: Expression, context: Optional[Type],
                            allow_none_return: bool = False) -> Type:
        with self.chk.binder.frame_context(can_skip=True, fall_through=0):
            if map is None:
                # We still need to type check node, in case we want to
                # process it for isinstance checks later
                self.accept(node, type_context=context, allow_none_return=allow_none_return)
                return UninhabitedType()
            self.chk.push_type_map(map)
            return self.accept(node, type_context=context, allow_none_return=allow_none_return)

    def visit_backquote_expr(self, e: BackquoteExpr) -> Type:
        self.accept(e.expr)
        return self.named_type('builtins.str')

    #
    # Helpers
    #

    def accept(self,
               node: Expression,
               type_context: Optional[Type] = None,
               allow_none_return: bool = False,
               always_allow_any: bool = False,
               ) -> Type:
        """Type check a node in the given type context.  If allow_none_return
        is True and this expression is a call, allow it to return None.  This
        applies only to this expression and not any subexpressions.
        """
        if node in self.type_overrides:
            return self.type_overrides[node]
        self.type_context.append(type_context)
        try:
            if allow_none_return and isinstance(node, CallExpr):
                typ = self.visit_call_expr(node, allow_none_return=True)
            elif allow_none_return and isinstance(node, YieldFromExpr):
                typ = self.visit_yield_from_expr(node, allow_none_return=True)
            elif allow_none_return and isinstance(node, ConditionalExpr):
                typ = self.visit_conditional_expr(node, allow_none_return=True)
            else:
                typ = node.accept(self)
        except Exception as err:
            report_internal_error(err, self.chk.errors.file,
                                  node.line, self.chk.errors, self.chk.options)

        self.type_context.pop()
        assert typ is not None
        self.chk.store_type(node, typ)

        if (self.chk.options.disallow_any_expr and
                not always_allow_any and
                not self.chk.is_stub and
                self.chk.in_checked_function() and
                has_any_type(typ) and not self.chk.current_node_deferred):
            self.msg.disallowed_any_type(typ, node)

        if not self.chk.in_checked_function() or self.chk.current_node_deferred:
            return AnyType(TypeOfAny.unannotated)
        else:
            return typ

    def named_type(self, name: str) -> Instance:
        """Return an instance type with type given by the name and no type
        arguments. Alias for TypeChecker.named_type.
        """
        return self.chk.named_type(name)

    def is_valid_var_arg(self, typ: Type) -> bool:
        """Is a type valid as a *args argument?"""
        typ = get_proper_type(typ)
        return (isinstance(typ, TupleType) or
                is_subtype(typ, self.chk.named_generic_type('typing.Iterable',
                                                            [AnyType(TypeOfAny.special_form)])) or
                isinstance(typ, AnyType) or
                isinstance(typ, ParamSpecType))

    def is_valid_keyword_var_arg(self, typ: Type) -> bool:
        """Is a type valid as a **kwargs argument?"""
        ret = (
                is_subtype(typ, self.chk.named_generic_type('typing.Mapping',
                    [self.named_type('builtins.str'), AnyType(TypeOfAny.special_form)])) or
                is_subtype(typ, self.chk.named_generic_type('typing.Mapping',
                    [UninhabitedType(), UninhabitedType()])) or
                isinstance(typ, ParamSpecType)
        )
        if self.chk.options.python_version[0] < 3:
            ret = ret or is_subtype(typ, self.chk.named_generic_type('typing.Mapping',
                [self.named_type('builtins.unicode'), AnyType(TypeOfAny.special_form)]))
        return ret

    def has_member(self, typ: Type, member: str) -> bool:
        """Does type have member with the given name?"""
        # TODO: refactor this to use checkmember.analyze_member_access, otherwise
        # these two should be carefully kept in sync.
        # This is much faster than analyze_member_access, though, and so using
        # it first as a filter is important for performance.
        typ = get_proper_type(typ)

        if isinstance(typ, TypeVarType):
            typ = get_proper_type(typ.upper_bound)
        if isinstance(typ, TupleType):
            typ = tuple_fallback(typ)
        if isinstance(typ, LiteralType):
            typ = typ.fallback
        if isinstance(typ, Instance):
            return typ.type.has_readable_member(member)
        if isinstance(typ, CallableType) and typ.is_type_obj():
            return typ.fallback.type.has_readable_member(member)
        elif isinstance(typ, AnyType):
            return True
        elif isinstance(typ, UnionType):
            result = all(self.has_member(x, member) for x in typ.relevant_items())
            return result
        elif isinstance(typ, TypeType):
            # Type[Union[X, ...]] is always normalized to Union[Type[X], ...],
            # so we don't need to care about unions here.
            item = typ.item
            if isinstance(item, TypeVarType):
                item = get_proper_type(item.upper_bound)
            if isinstance(item, TupleType):
                item = tuple_fallback(item)
            if isinstance(item, Instance) and item.type.metaclass_type is not None:
                return self.has_member(item.type.metaclass_type, member)
            if isinstance(item, AnyType):
                return True
            return False
        else:
            return False

    def not_ready_callback(self, name: str, context: Context) -> None:
        """Called when we can't infer the type of a variable because it's not ready yet.

        Either defer type checking of the enclosing function to the next
        pass or report an error.
        """
        self.chk.handle_cannot_determine_type(name, context)

    def visit_yield_expr(self, e: YieldExpr) -> Type:
        return_type = self.chk.return_types[-1]
        expected_item_type = self.chk.get_generator_yield_type(return_type, False)
        if e.expr is None:
            if (not isinstance(get_proper_type(expected_item_type), (NoneType, AnyType))
                    and self.chk.in_checked_function()):
                self.chk.fail(message_registry.YIELD_VALUE_EXPECTED, e)
        else:
            actual_item_type = self.accept(e.expr, expected_item_type)
            self.chk.check_subtype(actual_item_type, expected_item_type, e,
                                   message_registry.INCOMPATIBLE_TYPES_IN_YIELD,
                                   'actual type', 'expected type')
        return self.chk.get_generator_receive_type(return_type, False)

    def visit_await_expr(self, e: AwaitExpr) -> Type:
        expected_type = self.type_context[-1]
        if expected_type is not None:
            expected_type = self.chk.named_generic_type('typing.Awaitable', [expected_type])
        actual_type = get_proper_type(self.accept(e.expr, expected_type))
        if isinstance(actual_type, AnyType):
            return AnyType(TypeOfAny.from_another_any, source_any=actual_type)
        return self.check_awaitable_expr(actual_type, e,
                                         message_registry.INCOMPATIBLE_TYPES_IN_AWAIT)

    def check_awaitable_expr(self, t: Type, ctx: Context, msg: Union[str, ErrorMessage]) -> Type:
        """Check the argument to `await` and extract the type of value.

        Also used by `async for` and `async with`.
        """
        if not self.chk.check_subtype(t, self.named_type('typing.Awaitable'), ctx,
                                      msg, 'actual type', 'expected type'):
            return AnyType(TypeOfAny.special_form)
        else:
            generator = self.check_method_call_by_name('__await__', t, [], [], ctx)[0]
            ret_type = self.chk.get_generator_return_type(generator, False)
            ret_type = get_proper_type(ret_type)
            if isinstance(ret_type, UninhabitedType) and not ret_type.ambiguous:
                self.chk.binder.unreachable()
            return ret_type

    def visit_yield_from_expr(self, e: YieldFromExpr, allow_none_return: bool = False) -> Type:
        # NOTE: Whether `yield from` accepts an `async def` decorated
        # with `@types.coroutine` (or `@asyncio.coroutine`) depends on
        # whether the generator containing the `yield from` is itself
        # thus decorated.  But it accepts a generator regardless of
        # how it's decorated.
        return_type = self.chk.return_types[-1]
        # TODO: What should the context for the sub-expression be?
        # If the containing function has type Generator[X, Y, ...],
        # the context should be Generator[X, Y, T], where T is the
        # context of the 'yield from' itself (but it isn't known).
        subexpr_type = get_proper_type(self.accept(e.expr))

        # Check that the expr is an instance of Iterable and get the type of the iterator produced
        # by __iter__.
        if isinstance(subexpr_type, AnyType):
            iter_type: Type = AnyType(TypeOfAny.from_another_any, source_any=subexpr_type)
        elif self.chk.type_is_iterable(subexpr_type):
            if is_async_def(subexpr_type) and not has_coroutine_decorator(return_type):
                self.chk.msg.yield_from_invalid_operand_type(subexpr_type, e)

            any_type = AnyType(TypeOfAny.special_form)
            generic_generator_type = self.chk.named_generic_type('typing.Generator',
                                                                 [any_type, any_type, any_type])
            iter_type, _ = self.check_method_call_by_name(
                '__iter__', subexpr_type, [], [], context=generic_generator_type)
        else:
            if not (is_async_def(subexpr_type) and has_coroutine_decorator(return_type)):
                self.chk.msg.yield_from_invalid_operand_type(subexpr_type, e)
                iter_type = AnyType(TypeOfAny.from_error)
            else:
                iter_type = self.check_awaitable_expr(
                    subexpr_type, e, message_registry.INCOMPATIBLE_TYPES_IN_YIELD_FROM)

        # Check that the iterator's item type matches the type yielded by the Generator function
        # containing this `yield from` expression.
        expected_item_type = self.chk.get_generator_yield_type(return_type, False)
        actual_item_type = self.chk.get_generator_yield_type(iter_type, False)

        self.chk.check_subtype(actual_item_type, expected_item_type, e,
                               message_registry.INCOMPATIBLE_TYPES_IN_YIELD_FROM,
                               'actual type', 'expected type')

        # Determine the type of the entire yield from expression.
        iter_type = get_proper_type(iter_type)
        if (isinstance(iter_type, Instance) and
                iter_type.type.fullname == 'typing.Generator'):
            expr_type = self.chk.get_generator_return_type(iter_type, False)
        else:
            # Non-Generators don't return anything from `yield from` expressions.
            # However special-case Any (which might be produced by an error).
            actual_item_type = get_proper_type(actual_item_type)
            if isinstance(actual_item_type, AnyType):
                expr_type = AnyType(TypeOfAny.from_another_any, source_any=actual_item_type)
            else:
                # Treat `Iterator[X]` as a shorthand for `Generator[X, None, Any]`.
                expr_type = NoneType()

        if not allow_none_return and isinstance(get_proper_type(expr_type), NoneType):
            self.chk.msg.does_not_return_value(None, e)
        return expr_type

    def visit_temp_node(self, e: TempNode) -> Type:
        return e.type

    def visit_type_var_expr(self, e: TypeVarExpr) -> Type:
        return AnyType(TypeOfAny.special_form)

    def visit_paramspec_expr(self, e: ParamSpecExpr) -> Type:
        return AnyType(TypeOfAny.special_form)

    def visit_newtype_expr(self, e: NewTypeExpr) -> Type:
        return AnyType(TypeOfAny.special_form)

    def visit_namedtuple_expr(self, e: NamedTupleExpr) -> Type:
        tuple_type = e.info.tuple_type
        if tuple_type:
            if (self.chk.options.disallow_any_unimported and
                    has_any_from_unimported_type(tuple_type)):
                self.msg.unimported_type_becomes_any("NamedTuple type", tuple_type, e)
            check_for_explicit_any(tuple_type, self.chk.options, self.chk.is_typeshed_stub,
                                   self.msg, context=e)
        return AnyType(TypeOfAny.special_form)

    def visit_enum_call_expr(self, e: EnumCallExpr) -> Type:
        for name, value in zip(e.items, e.values):
            if value is not None:
                typ = self.accept(value)
                if not isinstance(get_proper_type(typ), AnyType):
                    var = e.info.names[name].node
                    if isinstance(var, Var):
                        # Inline TypeChecker.set_inferred_type(),
                        # without the lvalue.  (This doesn't really do
                        # much, since the value attribute is defined
                        # to have type Any in the typeshed stub.)
                        var.type = typ
                        var.is_inferred = True
        return AnyType(TypeOfAny.special_form)

    def visit_typeddict_expr(self, e: TypedDictExpr) -> Type:
        return AnyType(TypeOfAny.special_form)

    def visit__promote_expr(self, e: PromoteExpr) -> Type:
        return e.type

    def visit_star_expr(self, e: StarExpr) -> StarType:
        return StarType(self.accept(e.expr))

    def object_type(self) -> Instance:
        """Return instance type 'object'."""
        return self.named_type('builtins.object')

    def bool_type(self) -> Instance:
        """Return instance type 'bool'."""
        return self.named_type('builtins.bool')

    @overload
    def narrow_type_from_binder(self, expr: Expression, known_type: Type) -> Type: ...

    @overload
    def narrow_type_from_binder(self, expr: Expression, known_type: Type,
                                skip_non_overlapping: bool) -> Optional[Type]: ...

    def narrow_type_from_binder(self, expr: Expression, known_type: Type,
                                skip_non_overlapping: bool = False) -> Optional[Type]:
        """Narrow down a known type of expression using information in conditional type binder.

        If 'skip_non_overlapping' is True, return None if the type and restriction are
        non-overlapping.
        """
        if literal(expr) >= LITERAL_TYPE:
            restriction = self.chk.binder.get(expr)
            # If the current node is deferred, some variables may get Any types that they
            # otherwise wouldn't have. We don't want to narrow down these since it may
            # produce invalid inferred Optional[Any] types, at least.
            if restriction and not (isinstance(get_proper_type(known_type), AnyType)
                                    and self.chk.current_node_deferred):
                # Note: this call should match the one in narrow_declared_type().
                if (skip_non_overlapping and
                        not is_overlapping_types(known_type, restriction,
                                                 prohibit_none_typevar_overlap=True)):
                    return None
                return narrow_declared_type(known_type, restriction)
        return known_type


def has_any_type(t: Type, ignore_in_type_obj: bool = False) -> bool:
    """Whether t contains an Any type"""
    return t.accept(HasAnyType(ignore_in_type_obj))


class HasAnyType(types.TypeQuery[bool]):
    def __init__(self, ignore_in_type_obj: bool) -> None:
        super().__init__(any)
        self.ignore_in_type_obj = ignore_in_type_obj

    def visit_any(self, t: AnyType) -> bool:
        return t.type_of_any != TypeOfAny.special_form  # special forms are not real Any types

    def visit_callable_type(self, t: CallableType) -> bool:
        if self.ignore_in_type_obj and t.is_type_obj():
            return False
        return super().visit_callable_type(t)


def has_coroutine_decorator(t: Type) -> bool:
    """Whether t came from a function decorated with `@coroutine`."""
    t = get_proper_type(t)
    return isinstance(t, Instance) and t.type.fullname == 'typing.AwaitableGenerator'


def is_async_def(t: Type) -> bool:
    """Whether t came from a function defined using `async def`."""
    # In check_func_def(), when we see a function decorated with
    # `@typing.coroutine` or `@async.coroutine`, we change the
    # return type to typing.AwaitableGenerator[...], so that its
    # type is compatible with either Generator or Awaitable.
    # But for the check here we need to know whether the original
    # function (before decoration) was an `async def`.  The
    # AwaitableGenerator type conveniently preserves the original
    # type as its 4th parameter (3rd when using 0-origin indexing
    # :-), so that we can recover that information here.
    # (We really need to see whether the original, undecorated
    # function was an `async def`, which is orthogonal to its
    # decorations.)
    t = get_proper_type(t)
    if (isinstance(t, Instance)
            and t.type.fullname == 'typing.AwaitableGenerator'
            and len(t.args) >= 4):
        t = get_proper_type(t.args[3])
    return isinstance(t, Instance) and t.type.fullname == 'typing.Coroutine'


def is_non_empty_tuple(t: Type) -> bool:
    t = get_proper_type(t)
    return isinstance(t, TupleType) and bool(t.items)


def is_duplicate_mapping(mapping: List[int],
                         actual_types: List[Type],
                         actual_kinds: List[ArgKind]) -> bool:
    return (
        len(mapping) > 1
        # Multiple actuals can map to the same formal if they both come from
        # varargs (*args and **kwargs); in this case at runtime it is possible
        # that here are no duplicates. We need to allow this, as the convention
        # f(..., *args, **kwargs) is common enough.
        and not (len(mapping) == 2
                 and actual_kinds[mapping[0]] == nodes.ARG_STAR
                 and actual_kinds[mapping[1]] == nodes.ARG_STAR2)
        # Multiple actuals can map to the same formal if there are multiple
        # **kwargs which cannot be mapped with certainty (non-TypedDict
        # **kwargs).
        and not all(actual_kinds[m] == nodes.ARG_STAR2 and
                    not isinstance(get_proper_type(actual_types[m]), TypedDictType)
                    for m in mapping)
    )


def replace_callable_return_type(c: CallableType, new_ret_type: Type) -> CallableType:
    """Return a copy of a callable type with a different return type."""
    return c.copy_modified(ret_type=new_ret_type)


class ArgInferSecondPassQuery(types.TypeQuery[bool]):
    """Query whether an argument type should be inferred in the second pass.

    The result is True if the type has a type variable in a callable return
    type anywhere. For example, the result for Callable[[], T] is True if t is
    a type variable.
    """
    def __init__(self) -> None:
        super().__init__(any)

    def visit_callable_type(self, t: CallableType) -> bool:
        return self.query_types(t.arg_types) or t.accept(HasTypeVarQuery())


class HasTypeVarQuery(types.TypeQuery[bool]):
    """Visitor for querying whether a type has a type variable component."""
    def __init__(self) -> None:
        super().__init__(any)

    def visit_type_var(self, t: TypeVarType) -> bool:
        return True


def has_erased_component(t: Optional[Type]) -> bool:
    return t is not None and t.accept(HasErasedComponentsQuery())


class HasErasedComponentsQuery(types.TypeQuery[bool]):
    """Visitor for querying whether a type has an erased component."""
    def __init__(self) -> None:
        super().__init__(any)

    def visit_erased_type(self, t: ErasedType) -> bool:
        return True


def has_uninhabited_component(t: Optional[Type]) -> bool:
    return t is not None and t.accept(HasUninhabitedComponentsQuery())


class HasUninhabitedComponentsQuery(types.TypeQuery[bool]):
    """Visitor for querying whether a type has an UninhabitedType component."""
    def __init__(self) -> None:
        super().__init__(any)

    def visit_uninhabited_type(self, t: UninhabitedType) -> bool:
        return True


def arg_approximate_similarity(actual: Type, formal: Type) -> bool:
    """Return if caller argument (actual) is roughly compatible with signature arg (formal).

    This function is deliberately loose and will report two types are similar
    as long as their "shapes" are plausibly the same.

    This is useful when we're doing error reporting: for example, if we're trying
    to select an overload alternative and there's no exact match, we can use
    this function to help us identify which alternative the user might have
    *meant* to match.
    """
    actual = get_proper_type(actual)
    formal = get_proper_type(formal)

    # Erase typevars: we'll consider them all to have the same "shape".
    if isinstance(actual, TypeVarType):
        actual = erase_to_union_or_bound(actual)
    if isinstance(formal, TypeVarType):
        formal = erase_to_union_or_bound(formal)

    # Callable or Type[...]-ish types
    def is_typetype_like(typ: ProperType) -> bool:
        return (isinstance(typ, TypeType)
                or (isinstance(typ, FunctionLike) and typ.is_type_obj())
                or (isinstance(typ, Instance) and typ.type.fullname == "builtins.type"))

    if isinstance(formal, CallableType):
        if isinstance(actual, (CallableType, Overloaded, TypeType)):
            return True
    if is_typetype_like(actual) and is_typetype_like(formal):
        return True

    # Unions
    if isinstance(actual, UnionType):
        return any(arg_approximate_similarity(item, formal) for item in actual.relevant_items())
    if isinstance(formal, UnionType):
        return any(arg_approximate_similarity(actual, item) for item in formal.relevant_items())

    # TypedDicts
    if isinstance(actual, TypedDictType):
        if isinstance(formal, TypedDictType):
            return True
        return arg_approximate_similarity(actual.fallback, formal)

    # Instances
    # For instances, we mostly defer to the existing is_subtype check.
    if isinstance(formal, Instance):
        if isinstance(actual, CallableType):
            actual = actual.fallback
        if isinstance(actual, Overloaded):
            actual = actual.items[0].fallback
        if isinstance(actual, TupleType):
            actual = tuple_fallback(actual)
        if isinstance(actual, Instance) and formal.type in actual.type.mro:
            # Try performing a quick check as an optimization
            return True

    # Fall back to a standard subtype check for the remaining kinds of type.
    return is_subtype(erasetype.erase_type(actual), erasetype.erase_type(formal))


def any_causes_overload_ambiguity(items: List[CallableType],
                                  return_types: List[Type],
                                  arg_types: List[Type],
                                  arg_kinds: List[ArgKind],
                                  arg_names: Optional[Sequence[Optional[str]]]) -> bool:
    """May an argument containing 'Any' cause ambiguous result type on call to overloaded function?

    Note that this sometimes returns True even if there is no ambiguity, since a correct
    implementation would be complex (and the call would be imprecisely typed due to Any
    types anyway).

    Args:
        items: Overload items matching the actual arguments
        arg_types: Actual argument types
        arg_kinds: Actual argument kinds
        arg_names: Actual argument names
    """
    if all_same_types(return_types):
        return False

    actual_to_formal = [
        map_formals_to_actuals(
            arg_kinds, arg_names, item.arg_kinds, item.arg_names, lambda i: arg_types[i])
        for item in items
    ]

    for arg_idx, arg_type in enumerate(arg_types):
        # We ignore Anys in type object callables as ambiguity
        # creators, since that can lead to falsely claiming ambiguity
        # for overloads between Type and Callable.
        if has_any_type(arg_type, ignore_in_type_obj=True):
            matching_formals_unfiltered = [(item_idx, lookup[arg_idx])
                                           for item_idx, lookup in enumerate(actual_to_formal)
                                           if lookup[arg_idx]]

            matching_returns = []
            matching_formals = []
            for item_idx, formals in matching_formals_unfiltered:
                matched_callable = items[item_idx]
                matching_returns.append(matched_callable.ret_type)

                # Note: if an actual maps to multiple formals of differing types within
                # a single callable, then we know at least one of those formals must be
                # a different type then the formal(s) in some other callable.
                # So it's safe to just append everything to the same list.
                for formal in formals:
                    matching_formals.append(matched_callable.arg_types[formal])
            if not all_same_types(matching_formals) and not all_same_types(matching_returns):
                # Any maps to multiple different types, and the return types of these items differ.
                return True
    return False


def all_same_types(types: List[Type]) -> bool:
    if len(types) == 0:
        return True
    return all(is_same_type(t, types[0]) for t in types[1:])


def merge_typevars_in_callables_by_name(
        callables: Sequence[CallableType]) -> Tuple[List[CallableType], List[TypeVarType]]:
    """Takes all the typevars present in the callables and 'combines' the ones with the same name.

    For example, suppose we have two callables with signatures "f(x: T, y: S) -> T" and
    "f(x: List[Tuple[T, S]]) -> Tuple[T, S]". Both callables use typevars named "T" and
    "S", but we treat them as distinct, unrelated typevars. (E.g. they could both have
    distinct ids.)

    If we pass in both callables into this function, it returns a list containing two
    new callables that are identical in signature, but use the same underlying TypeVarType
    for T and S.

    This is useful if we want to take the output lists and "merge" them into one callable
    in some way -- for example, when unioning together overloads.

    Returns both the new list of callables and a list of all distinct TypeVarType objects used.
    """
    output: List[CallableType] = []
    unique_typevars: Dict[str, TypeVarType] = {}
    variables: List[TypeVarType] = []

    for target in callables:
        if target.is_generic():
            target = freshen_function_type_vars(target)

            rename = {}  # Dict[TypeVarId, TypeVar]
            for tv in target.variables:
                name = tv.fullname
                if name not in unique_typevars:
                    # TODO(PEP612): fix for ParamSpecType
                    if isinstance(tv, ParamSpecType):
                        continue
                    assert isinstance(tv, TypeVarType)
                    unique_typevars[name] = tv
                    variables.append(tv)
                rename[tv.id] = unique_typevars[name]

            target = cast(CallableType, expand_type(target, rename))
        output.append(target)

    return output, variables


def try_getting_literal(typ: Type) -> ProperType:
    """If possible, get a more precise literal type for a given type."""
    typ = get_proper_type(typ)
    if isinstance(typ, Instance) and typ.last_known_value is not None:
        return typ.last_known_value
    return typ


def is_expr_literal_type(node: Expression) -> bool:
    """Returns 'true' if the given node is a Literal"""
    valid = ('typing.Literal', 'typing_extensions.Literal')
    if isinstance(node, IndexExpr):
        base = node.base
        return isinstance(base, RefExpr) and base.fullname in valid
    if isinstance(node, NameExpr):
        underlying = node.node
        return isinstance(underlying, TypeAlias) and isinstance(get_proper_type(underlying.target),
                                                                LiteralType)
    return False


def has_bytes_component(typ: Type, py2: bool = False) -> bool:
    """Is this one of builtin byte types, or a union that contains it?"""
    typ = get_proper_type(typ)
    if py2:
        byte_types = {'builtins.str', 'builtins.bytearray'}
    else:
        byte_types = {'builtins.bytes', 'builtins.bytearray'}
    if isinstance(typ, UnionType):
        return any(has_bytes_component(t) for t in typ.items)
    if isinstance(typ, Instance) and typ.type.fullname in byte_types:
        return True
    return False


def type_info_from_type(typ: Type) -> Optional[TypeInfo]:
    """Gets the TypeInfo for a type, indirecting through things like type variables and tuples."""
    typ = get_proper_type(typ)
    if isinstance(typ, FunctionLike) and typ.is_type_obj():
        return typ.type_object()
    if isinstance(typ, TypeType):
        typ = typ.item
    if isinstance(typ, TypeVarType):
        typ = get_proper_type(typ.upper_bound)
    if isinstance(typ, TupleType):
        typ = tuple_fallback(typ)
    if isinstance(typ, Instance):
        return typ.type

    # A complicated type. Too tricky, give up.
    # TODO: Do something more clever here.
    return None


def is_operator_method(fullname: Optional[str]) -> bool:
    if fullname is None:
        return False
    short_name = fullname.split('.')[-1]
    return (
        short_name in operators.op_methods.values() or
        short_name in operators.reverse_op_methods.values() or
        short_name in operators.unary_op_methods.values())


def get_partial_instance_type(t: Optional[Type]) -> Optional[PartialType]:
    if t is None or not isinstance(t, PartialType) or t.type is None:
        return None
    return t
