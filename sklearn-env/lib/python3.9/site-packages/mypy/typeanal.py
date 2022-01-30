"""Semantic analysis of types"""

import itertools
from itertools import chain
from contextlib import contextmanager
from mypy.backports import OrderedDict

from typing import Callable, List, Optional, Set, Tuple, Iterator, TypeVar, Iterable, Sequence
from typing_extensions import Final
from mypy_extensions import DefaultNamedArg

from mypy.messages import MessageBuilder, quote_type_string, format_type_bare
from mypy.options import Options
from mypy.types import (
    Type, UnboundType, TupleType, TypedDictType, UnionType, Instance, AnyType,
    CallableType, NoneType, ErasedType, DeletedType, TypeList, TypeVarType, SyntheticTypeVisitor,
    StarType, PartialType, EllipsisType, UninhabitedType, TypeType, CallableArgument,
    TypeQuery, union_items, TypeOfAny, LiteralType, RawExpressionType,
    PlaceholderType, Overloaded, get_proper_type, TypeAliasType, RequiredType,
    TypeVarLikeType, ParamSpecType, ParamSpecFlavor, callable_with_ellipsis
)

from mypy.nodes import (
    TypeInfo, Context, SymbolTableNode, Var, Expression,
    get_nongen_builtins, check_arg_names, check_arg_kinds, ArgKind, ARG_POS, ARG_NAMED,
    ARG_OPT, ARG_NAMED_OPT, ARG_STAR, ARG_STAR2, TypeVarExpr, TypeVarLikeExpr, ParamSpecExpr,
    TypeAlias, PlaceholderNode, SYMBOL_FUNCBASE_TYPES, Decorator, MypyFile
)
from mypy.typetraverser import TypeTraverserVisitor
from mypy.tvar_scope import TypeVarLikeScope
from mypy.exprtotype import expr_to_unanalyzed_type, TypeTranslationError
from mypy.plugin import Plugin, TypeAnalyzerPluginInterface, AnalyzeTypeContext
from mypy.semanal_shared import SemanticAnalyzerCoreInterface
from mypy.errorcodes import ErrorCode
from mypy import nodes, message_registry, errorcodes as codes

T = TypeVar('T')

type_constructors: Final = {
    'typing.Callable',
    'typing.Optional',
    'typing.Tuple',
    'typing.Type',
    'typing.Union',
    'typing.Literal',
    'typing_extensions.Literal',
    'typing.Annotated',
    'typing_extensions.Annotated',
}

ARG_KINDS_BY_CONSTRUCTOR: Final = {
    'mypy_extensions.Arg': ARG_POS,
    'mypy_extensions.DefaultArg': ARG_OPT,
    'mypy_extensions.NamedArg': ARG_NAMED,
    'mypy_extensions.DefaultNamedArg': ARG_NAMED_OPT,
    'mypy_extensions.VarArg': ARG_STAR,
    'mypy_extensions.KwArg': ARG_STAR2,
}

GENERIC_STUB_NOT_AT_RUNTIME_TYPES: Final = {
    'queue.Queue',
    'builtins._PathLike',
    'asyncio.futures.Future',
}


def analyze_type_alias(node: Expression,
                       api: SemanticAnalyzerCoreInterface,
                       tvar_scope: TypeVarLikeScope,
                       plugin: Plugin,
                       options: Options,
                       is_typeshed_stub: bool,
                       allow_placeholder: bool = False,
                       in_dynamic_func: bool = False,
                       global_scope: bool = True) -> Optional[Tuple[Type, Set[str]]]:
    """Analyze r.h.s. of a (potential) type alias definition.

    If `node` is valid as a type alias rvalue, return the resulting type and a set of
    full names of type aliases it depends on (directly or indirectly).
    Return None otherwise. 'node' must have been semantically analyzed.
    """
    try:
        type = expr_to_unanalyzed_type(node, options, api.is_stub_file)
    except TypeTranslationError:
        api.fail('Invalid type alias: expression is not a valid type', node)
        return None
    analyzer = TypeAnalyser(api, tvar_scope, plugin, options, is_typeshed_stub,
                            defining_alias=True,
                            allow_placeholder=allow_placeholder)
    analyzer.in_dynamic_func = in_dynamic_func
    analyzer.global_scope = global_scope
    res = type.accept(analyzer)
    return res, analyzer.aliases_used


def no_subscript_builtin_alias(name: str, propose_alt: bool = True) -> str:
    msg = '"{}" is not subscriptable'.format(name.split('.')[-1])
    # This should never be called if the python_version is 3.9 or newer
    nongen_builtins = get_nongen_builtins((3, 8))
    replacement = nongen_builtins[name]
    if replacement and propose_alt:
        msg += ', use "{}" instead'.format(replacement)
    return msg


class TypeAnalyser(SyntheticTypeVisitor[Type], TypeAnalyzerPluginInterface):
    """Semantic analyzer for types.

    Converts unbound types into bound types. This is a no-op for already
    bound types.

    If an incomplete reference is encountered, this does a defer. The
    caller never needs to defer.
    """

    # Is this called from an untyped function definition?
    in_dynamic_func: bool = False
    # Is this called from global scope?
    global_scope: bool = True

    def __init__(self,
                 api: SemanticAnalyzerCoreInterface,
                 tvar_scope: TypeVarLikeScope,
                 plugin: Plugin,
                 options: Options,
                 is_typeshed_stub: bool, *,
                 defining_alias: bool = False,
                 allow_tuple_literal: bool = False,
                 allow_unbound_tvars: bool = False,
                 allow_placeholder: bool = False,
                 allow_required: bool = False,
                 report_invalid_types: bool = True) -> None:
        self.api = api
        self.lookup_qualified = api.lookup_qualified
        self.lookup_fqn_func = api.lookup_fully_qualified
        self.fail_func = api.fail
        self.note_func = api.note
        self.tvar_scope = tvar_scope
        # Are we analysing a type alias definition rvalue?
        self.defining_alias = defining_alias
        self.allow_tuple_literal = allow_tuple_literal
        # Positive if we are analyzing arguments of another (outer) type
        self.nesting_level = 0
        # Should we allow new type syntax when targeting older Python versions
        # like 'list[int]' or 'X | Y' (allowed in stubs and with `__future__` import)?
        self.always_allow_new_syntax = (
            self.api.is_stub_file
            or self.api.is_future_flag_set('annotations')
        )
        # Should we accept unbound type variables (always OK in aliases)?
        self.allow_unbound_tvars = allow_unbound_tvars or defining_alias
        # If false, record incomplete ref if we generate PlaceholderType.
        self.allow_placeholder = allow_placeholder
        # Are we in a context where Required[] is allowed?
        self.allow_required = allow_required
        # Should we report an error whenever we encounter a RawExpressionType outside
        # of a Literal context: e.g. whenever we encounter an invalid type? Normally,
        # we want to report an error, but the caller may want to do more specialized
        # error handling.
        self.report_invalid_types = report_invalid_types
        self.plugin = plugin
        self.options = options
        self.is_typeshed_stub = is_typeshed_stub
        # Names of type aliases encountered while analysing a type will be collected here.
        self.aliases_used: Set[str] = set()

    def visit_unbound_type(self, t: UnboundType, defining_literal: bool = False) -> Type:
        typ = self.visit_unbound_type_nonoptional(t, defining_literal)
        if t.optional:
            # We don't need to worry about double-wrapping Optionals or
            # wrapping Anys: Union simplification will take care of that.
            return make_optional_type(typ)
        return typ

    def visit_unbound_type_nonoptional(self, t: UnboundType, defining_literal: bool) -> Type:
        sym = self.lookup_qualified(t.name, t)
        if sym is not None:
            node = sym.node
            if isinstance(node, PlaceholderNode):
                if node.becomes_typeinfo:
                    # Reference to placeholder type.
                    if self.api.final_iteration:
                        self.cannot_resolve_type(t)
                        return AnyType(TypeOfAny.from_error)
                    elif self.allow_placeholder:
                        self.api.defer()
                    else:
                        self.api.record_incomplete_ref()
                    return PlaceholderType(node.fullname, self.anal_array(t.args), t.line)
                else:
                    if self.api.final_iteration:
                        self.cannot_resolve_type(t)
                        return AnyType(TypeOfAny.from_error)
                    else:
                        # Reference to an unknown placeholder node.
                        self.api.record_incomplete_ref()
                        return AnyType(TypeOfAny.special_form)
            if node is None:
                self.fail('Internal error (node is None, kind={})'.format(sym.kind), t)
                return AnyType(TypeOfAny.special_form)
            fullname = node.fullname
            hook = self.plugin.get_type_analyze_hook(fullname)
            if hook is not None:
                return hook(AnalyzeTypeContext(t, t, self))
            if (fullname in get_nongen_builtins(self.options.python_version)
                    and t.args
                    and not self.always_allow_new_syntax):
                self.fail(no_subscript_builtin_alias(fullname,
                                                     propose_alt=not self.defining_alias), t)
            tvar_def = self.tvar_scope.get_binding(sym)
            if isinstance(sym.node, ParamSpecExpr):
                if tvar_def is None:
                    self.fail('ParamSpec "{}" is unbound'.format(t.name), t)
                    return AnyType(TypeOfAny.from_error)
                assert isinstance(tvar_def, ParamSpecType)
                if len(t.args) > 0:
                    self.fail('ParamSpec "{}" used with arguments'.format(t.name), t)
                # Change the line number
                return ParamSpecType(
                    tvar_def.name, tvar_def.fullname, tvar_def.id, tvar_def.flavor,
                    tvar_def.upper_bound, line=t.line, column=t.column,
                )
            if isinstance(sym.node, TypeVarExpr) and tvar_def is not None and self.defining_alias:
                self.fail('Can\'t use bound type variable "{}"'
                          ' to define generic alias'.format(t.name), t)
                return AnyType(TypeOfAny.from_error)
            if isinstance(sym.node, TypeVarExpr) and tvar_def is not None:
                assert isinstance(tvar_def, TypeVarType)
                if len(t.args) > 0:
                    self.fail('Type variable "{}" used with arguments'.format(t.name), t)
                # Change the line number
                return TypeVarType(
                    tvar_def.name, tvar_def.fullname, tvar_def.id, tvar_def.values,
                    tvar_def.upper_bound, tvar_def.variance, line=t.line, column=t.column,
                )
            special = self.try_analyze_special_unbound_type(t, fullname)
            if special is not None:
                return special
            if isinstance(node, TypeAlias):
                self.aliases_used.add(fullname)
                an_args = self.anal_array(t.args)
                disallow_any = self.options.disallow_any_generics and not self.is_typeshed_stub
                res = expand_type_alias(node, an_args, self.fail, node.no_args, t,
                                        unexpanded_type=t,
                                        disallow_any=disallow_any)
                # The only case where expand_type_alias() can return an incorrect instance is
                # when it is top-level instance, so no need to recurse.
                if (isinstance(res, Instance) and  # type: ignore[misc]
                        len(res.args) != len(res.type.type_vars) and
                        not self.defining_alias):
                    fix_instance(
                        res,
                        self.fail,
                        self.note,
                        disallow_any=disallow_any,
                        python_version=self.options.python_version,
                        use_generic_error=True,
                        unexpanded_type=t)
                if node.eager:
                    # TODO: Generate error if recursive (once we have recursive types)
                    res = get_proper_type(res)
                return res
            elif isinstance(node, TypeInfo):
                return self.analyze_type_with_type_info(node, t.args, t)
            elif node.fullname in ("typing_extensions.TypeAlias", "typing.TypeAlias"):
                return AnyType(TypeOfAny.special_form)
            else:
                return self.analyze_unbound_type_without_type_info(t, sym, defining_literal)
        else:  # sym is None
            return AnyType(TypeOfAny.special_form)

    def cannot_resolve_type(self, t: UnboundType) -> None:
        # TODO: Move error message generation to messages.py. We'd first
        #       need access to MessageBuilder here. Also move the similar
        #       message generation logic in semanal.py.
        self.api.fail(
            'Cannot resolve name "{}" (possible cyclic definition)'.format(t.name),
            t)

    def try_analyze_special_unbound_type(self, t: UnboundType, fullname: str) -> Optional[Type]:
        """Bind special type that is recognized through magic name such as 'typing.Any'.

        Return the bound type if successful, and return None if the type is a normal type.
        """
        if fullname == 'builtins.None':
            return NoneType()
        elif fullname == 'typing.Any' or fullname == 'builtins.Any':
            return AnyType(TypeOfAny.explicit)
        elif fullname in ('typing.Final', 'typing_extensions.Final'):
            self.fail("Final can be only used as an outermost qualifier"
                      " in a variable annotation", t)
            return AnyType(TypeOfAny.from_error)
        elif (fullname == 'typing.Tuple' or
             (fullname == 'builtins.tuple'
                and (self.always_allow_new_syntax or self.options.python_version >= (3, 9)))):
            # Tuple is special because it is involved in builtin import cycle
            # and may be not ready when used.
            sym = self.api.lookup_fully_qualified_or_none('builtins.tuple')
            if not sym or isinstance(sym.node, PlaceholderNode):
                if self.api.is_incomplete_namespace('builtins'):
                    self.api.record_incomplete_ref()
                else:
                    self.fail('Name "tuple" is not defined', t)
                return AnyType(TypeOfAny.special_form)
            if len(t.args) == 0 and not t.empty_tuple_index:
                # Bare 'Tuple' is same as 'tuple'
                any_type = self.get_omitted_any(t)
                return self.named_type('builtins.tuple', [any_type],
                                       line=t.line, column=t.column)
            if len(t.args) == 2 and isinstance(t.args[1], EllipsisType):
                # Tuple[T, ...] (uniform, variable-length tuple)
                instance = self.named_type('builtins.tuple', [self.anal_type(t.args[0])])
                instance.line = t.line
                return instance
            return self.tuple_type(self.anal_array(t.args))
        elif fullname == 'typing.Union':
            items = self.anal_array(t.args)
            return UnionType.make_union(items)
        elif fullname == 'typing.Optional':
            if len(t.args) != 1:
                self.fail('Optional[...] must have exactly one type argument', t)
                return AnyType(TypeOfAny.from_error)
            item = self.anal_type(t.args[0])
            return make_optional_type(item)
        elif fullname == 'typing.Callable':
            return self.analyze_callable_type(t)
        elif (fullname == 'typing.Type' or
             (fullname == 'builtins.type'
                and (self.always_allow_new_syntax or self.options.python_version >= (3, 9)))):
            if len(t.args) == 0:
                if fullname == 'typing.Type':
                    any_type = self.get_omitted_any(t)
                    return TypeType(any_type, line=t.line, column=t.column)
                else:
                    # To prevent assignment of 'builtins.type' inferred as 'builtins.object'
                    # See https://github.com/python/mypy/issues/9476 for more information
                    return None
            if len(t.args) != 1:
                type_str = 'Type[...]' if fullname == 'typing.Type' else 'type[...]'
                self.fail(type_str + ' must have exactly one type argument', t)
            item = self.anal_type(t.args[0])
            return TypeType.make_normalized(item, line=t.line)
        elif fullname == 'typing.ClassVar':
            if self.nesting_level > 0:
                self.fail('Invalid type: ClassVar nested inside other type', t)
            if len(t.args) == 0:
                return AnyType(TypeOfAny.from_omitted_generics, line=t.line, column=t.column)
            if len(t.args) != 1:
                self.fail('ClassVar[...] must have at most one type argument', t)
                return AnyType(TypeOfAny.from_error)
            return self.anal_type(t.args[0])
        elif fullname in ('mypy_extensions.NoReturn', 'typing.NoReturn'):
            return UninhabitedType(is_noreturn=True)
        elif fullname in ('typing_extensions.Literal', 'typing.Literal'):
            return self.analyze_literal_type(t)
        elif fullname in ('typing_extensions.Annotated', 'typing.Annotated'):
            if len(t.args) < 2:
                self.fail("Annotated[...] must have exactly one type argument"
                          " and at least one annotation", t)
                return AnyType(TypeOfAny.from_error)
            return self.anal_type(t.args[0])
        elif fullname in ('typing_extensions.Required', 'typing.Required'):
            if not self.allow_required:
                self.fail("Required[] can be only used in a TypedDict definition", t)
                return AnyType(TypeOfAny.from_error)
            if len(t.args) != 1:
                self.fail("Required[] must have exactly one type argument", t)
                return AnyType(TypeOfAny.from_error)
            return RequiredType(self.anal_type(t.args[0]), required=True)
        elif fullname in ('typing_extensions.NotRequired', 'typing.NotRequired'):
            if not self.allow_required:
                self.fail("NotRequired[] can be only used in a TypedDict definition", t)
                return AnyType(TypeOfAny.from_error)
            if len(t.args) != 1:
                self.fail("NotRequired[] must have exactly one type argument", t)
                return AnyType(TypeOfAny.from_error)
            return RequiredType(self.anal_type(t.args[0]), required=False)
        elif self.anal_type_guard_arg(t, fullname) is not None:
            # In most contexts, TypeGuard[...] acts as an alias for bool (ignoring its args)
            return self.named_type('builtins.bool')
        return None

    def get_omitted_any(self, typ: Type, fullname: Optional[str] = None) -> AnyType:
        disallow_any = not self.is_typeshed_stub and self.options.disallow_any_generics
        return get_omitted_any(disallow_any, self.fail, self.note, typ,
                               self.options.python_version, fullname)

    def analyze_type_with_type_info(
            self, info: TypeInfo, args: Sequence[Type], ctx: Context) -> Type:
        """Bind unbound type when were able to find target TypeInfo.

        This handles simple cases like 'int', 'modname.UserClass[str]', etc.
        """

        if len(args) > 0 and info.fullname == 'builtins.tuple':
            fallback = Instance(info, [AnyType(TypeOfAny.special_form)], ctx.line)
            return TupleType(self.anal_array(args), fallback, ctx.line)
        # Analyze arguments and (usually) construct Instance type. The
        # number of type arguments and their values are
        # checked only later, since we do not always know the
        # valid count at this point. Thus we may construct an
        # Instance with an invalid number of type arguments.
        instance = Instance(info, self.anal_array(args, allow_param_spec=True),
                            ctx.line, ctx.column)
        # Check type argument count.
        if len(instance.args) != len(info.type_vars) and not self.defining_alias:
            fix_instance(instance, self.fail, self.note,
                         disallow_any=self.options.disallow_any_generics and
                         not self.is_typeshed_stub,
                         python_version=self.options.python_version)

        tup = info.tuple_type
        if tup is not None:
            # The class has a Tuple[...] base class so it will be
            # represented as a tuple type.
            if args:
                self.fail('Generic tuple types not supported', ctx)
                return AnyType(TypeOfAny.from_error)
            return tup.copy_modified(items=self.anal_array(tup.items),
                                     fallback=instance)
        td = info.typeddict_type
        if td is not None:
            # The class has a TypedDict[...] base class so it will be
            # represented as a typeddict type.
            if args:
                self.fail('Generic TypedDict types not supported', ctx)
                return AnyType(TypeOfAny.from_error)
            # Create a named TypedDictType
            return td.copy_modified(item_types=self.anal_array(list(td.items.values())),
                                    fallback=instance)
        return instance

    def analyze_unbound_type_without_type_info(self, t: UnboundType, sym: SymbolTableNode,
                                               defining_literal: bool) -> Type:
        """Figure out what an unbound type that doesn't refer to a TypeInfo node means.

        This is something unusual. We try our best to find out what it is.
        """
        name = sym.fullname
        if name is None:
            assert sym.node is not None
            name = sym.node.name
        # Option 1:
        # Something with an Any type -- make it an alias for Any in a type
        # context. This is slightly problematic as it allows using the type 'Any'
        # as a base class -- however, this will fail soon at runtime so the problem
        # is pretty minor.
        if isinstance(sym.node, Var):
            typ = get_proper_type(sym.node.type)
            if isinstance(typ, AnyType):
                return AnyType(TypeOfAny.from_unimported_type,
                               missing_import_name=typ.missing_import_name)
        # Option 2:
        # Unbound type variable. Currently these may be still valid,
        # for example when defining a generic type alias.
        unbound_tvar = (isinstance(sym.node, TypeVarExpr) and
                        self.tvar_scope.get_binding(sym) is None)
        if self.allow_unbound_tvars and unbound_tvar:
            return t

        # Option 3:
        # Enum value. Note: we only want to return a LiteralType when
        # we're using this enum value specifically within context of
        # a "Literal[...]" type. So, if `defining_literal` is not set,
        # we bail out early with an error.
        #
        # If, in the distant future, we decide to permit things like
        # `def foo(x: Color.RED) -> None: ...`, we can remove that
        # check entirely.
        if isinstance(sym.node, Var) and sym.node.info and sym.node.info.is_enum:
            value = sym.node.name
            base_enum_short_name = sym.node.info.name
            if not defining_literal:
                msg = message_registry.INVALID_TYPE_RAW_ENUM_VALUE.format(
                    base_enum_short_name, value)
                self.fail(msg, t)
                return AnyType(TypeOfAny.from_error)
            return LiteralType(
                value=value,
                fallback=Instance(sym.node.info, [], line=t.line, column=t.column),
                line=t.line,
                column=t.column,
            )

        # None of the above options worked. We parse the args (if there are any)
        # to make sure there are no remaining semanal-only types, then give up.
        t = t.copy_modified(args=self.anal_array(t.args))
        # TODO: Move this message building logic to messages.py.
        notes: List[str] = []
        if isinstance(sym.node, Var):
            notes.append('See https://mypy.readthedocs.io/en/'
                         'stable/common_issues.html#variables-vs-type-aliases')
            message = 'Variable "{}" is not valid as a type'
        elif isinstance(sym.node, (SYMBOL_FUNCBASE_TYPES, Decorator)):
            message = 'Function "{}" is not valid as a type'
            notes.append('Perhaps you need "Callable[...]" or a callback protocol?')
        elif isinstance(sym.node, MypyFile):
            # TODO: suggest a protocol when supported.
            message = 'Module "{}" is not valid as a type'
        elif unbound_tvar:
            message = 'Type variable "{}" is unbound'
            short = name.split('.')[-1]
            notes.append(('(Hint: Use "Generic[{}]" or "Protocol[{}]" base class'
                          ' to bind "{}" inside a class)').format(short, short, short))
            notes.append('(Hint: Use "{}" in function signature to bind "{}"'
                         ' inside a function)'.format(short, short))
        else:
            message = 'Cannot interpret reference "{}" as a type'
        self.fail(message.format(name), t, code=codes.VALID_TYPE)
        for note in notes:
            self.note(note, t, code=codes.VALID_TYPE)

        # TODO: Would it be better to always return Any instead of UnboundType
        # in case of an error? On one hand, UnboundType has a name so error messages
        # are more detailed, on the other hand, some of them may be bogus,
        # see https://github.com/python/mypy/issues/4987.
        return t

    def visit_any(self, t: AnyType) -> Type:
        return t

    def visit_none_type(self, t: NoneType) -> Type:
        return t

    def visit_uninhabited_type(self, t: UninhabitedType) -> Type:
        return t

    def visit_erased_type(self, t: ErasedType) -> Type:
        # This type should exist only temporarily during type inference
        assert False, "Internal error: Unexpected erased type"

    def visit_deleted_type(self, t: DeletedType) -> Type:
        return t

    def visit_type_list(self, t: TypeList) -> Type:
        self.fail('Bracketed expression "[...]" is not valid as a type', t)
        self.note('Did you mean "List[...]"?', t)
        return AnyType(TypeOfAny.from_error)

    def visit_callable_argument(self, t: CallableArgument) -> Type:
        self.fail('Invalid type', t)
        return AnyType(TypeOfAny.from_error)

    def visit_instance(self, t: Instance) -> Type:
        return t

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        # TODO: should we do something here?
        return t

    def visit_type_var(self, t: TypeVarType) -> Type:
        return t

    def visit_param_spec(self, t: ParamSpecType) -> Type:
        return t

    def visit_callable_type(self, t: CallableType, nested: bool = True) -> Type:
        # Every Callable can bind its own type variables, if they're not in the outer scope
        with self.tvar_scope_frame():
            if self.defining_alias:
                variables = t.variables
            else:
                variables = self.bind_function_type_variables(t, t)
            special = self.anal_type_guard(t.ret_type)
            arg_kinds = t.arg_kinds
            if len(arg_kinds) >= 2 and arg_kinds[-2] == ARG_STAR and arg_kinds[-1] == ARG_STAR2:
                arg_types = self.anal_array(t.arg_types[:-2], nested=nested) + [
                    self.anal_star_arg_type(t.arg_types[-2], ARG_STAR, nested=nested),
                    self.anal_star_arg_type(t.arg_types[-1], ARG_STAR2, nested=nested),
                ]
            else:
                arg_types = self.anal_array(t.arg_types, nested=nested)
            ret = t.copy_modified(arg_types=arg_types,
                                  ret_type=self.anal_type(t.ret_type, nested=nested),
                                  # If the fallback isn't filled in yet,
                                  # its type will be the falsey FakeInfo
                                  fallback=(t.fallback if t.fallback.type
                                            else self.named_type('builtins.function')),
                                  variables=self.anal_var_defs(variables),
                                  type_guard=special,
                                  )
        return ret

    def anal_type_guard(self, t: Type) -> Optional[Type]:
        if isinstance(t, UnboundType):
            sym = self.lookup_qualified(t.name, t)
            if sym is not None and sym.node is not None:
                return self.anal_type_guard_arg(t, sym.node.fullname)
        # TODO: What if it's an Instance? Then use t.type.fullname?
        return None

    def anal_type_guard_arg(self, t: UnboundType, fullname: str) -> Optional[Type]:
        if fullname in ('typing_extensions.TypeGuard', 'typing.TypeGuard'):
            if len(t.args) != 1:
                self.fail("TypeGuard must have exactly one type argument", t)
                return AnyType(TypeOfAny.from_error)
            return self.anal_type(t.args[0])
        return None

    def anal_star_arg_type(self, t: Type, kind: ArgKind, nested: bool) -> Type:
        """Analyze signature argument type for *args and **kwargs argument."""
        # TODO: Check that suffix and kind match
        if isinstance(t, UnboundType) and t.name and '.' in t.name and not t.args:
            components = t.name.split('.')
            sym = self.lookup_qualified('.'.join(components[:-1]), t)
            if sym is not None and isinstance(sym.node, ParamSpecExpr):
                tvar_def = self.tvar_scope.get_binding(sym)
                if isinstance(tvar_def, ParamSpecType):
                    if kind == ARG_STAR:
                        flavor = ParamSpecFlavor.ARGS
                    elif kind == ARG_STAR2:
                        flavor = ParamSpecFlavor.KWARGS
                    else:
                        assert False, kind
                    return ParamSpecType(tvar_def.name, tvar_def.fullname, tvar_def.id, flavor,
                                         upper_bound=self.named_type('builtins.object'),
                                         line=t.line, column=t.column)
        return self.anal_type(t, nested=nested)

    def visit_overloaded(self, t: Overloaded) -> Type:
        # Overloaded types are manually constructed in semanal.py by analyzing the
        # AST and combining together the Callable types this visitor converts.
        #
        # So if we're ever asked to reanalyze an Overloaded type, we know it's
        # fine to just return it as-is.
        return t

    def visit_tuple_type(self, t: TupleType) -> Type:
        # Types such as (t1, t2, ...) only allowed in assignment statements. They'll
        # generate errors elsewhere, and Tuple[t1, t2, ...] must be used instead.
        if t.implicit and not self.allow_tuple_literal:
            self.fail('Syntax error in type annotation', t, code=codes.SYNTAX)
            if len(t.items) == 0:
                self.note('Suggestion: Use Tuple[()] instead of () for an empty tuple, or '
                'None for a function without a return value', t, code=codes.SYNTAX)
            elif len(t.items) == 1:
                self.note('Suggestion: Is there a spurious trailing comma?', t, code=codes.SYNTAX)
            else:
                self.note('Suggestion: Use Tuple[T1, ..., Tn] instead of (T1, ..., Tn)', t,
                          code=codes.SYNTAX)
            return AnyType(TypeOfAny.from_error)
        star_count = sum(1 for item in t.items if isinstance(item, StarType))
        if star_count > 1:
            self.fail('At most one star type allowed in a tuple', t)
            if t.implicit:
                return TupleType([AnyType(TypeOfAny.from_error) for _ in t.items],
                                 self.named_type('builtins.tuple'),
                                 t.line)
            else:
                return AnyType(TypeOfAny.from_error)
        any_type = AnyType(TypeOfAny.special_form)
        # If the fallback isn't filled in yet, its type will be the falsey FakeInfo
        fallback = (t.partial_fallback if t.partial_fallback.type
                    else self.named_type('builtins.tuple', [any_type]))
        return TupleType(self.anal_array(t.items), fallback, t.line)

    def visit_typeddict_type(self, t: TypedDictType) -> Type:
        items = OrderedDict([
            (item_name, self.anal_type(item_type))
            for (item_name, item_type) in t.items.items()
        ])
        return TypedDictType(items, set(t.required_keys), t.fallback)

    def visit_raw_expression_type(self, t: RawExpressionType) -> Type:
        # We should never see a bare Literal. We synthesize these raw literals
        # in the earlier stages of semantic analysis, but those
        # "fake literals" should always be wrapped in an UnboundType
        # corresponding to 'Literal'.
        #
        # Note: if at some point in the distant future, we decide to
        # make signatures like "foo(x: 20) -> None" legal, we can change
        # this method so it generates and returns an actual LiteralType
        # instead.

        if self.report_invalid_types:
            if t.base_type_name in ('builtins.int', 'builtins.bool'):
                # The only time it makes sense to use an int or bool is inside of
                # a literal type.
                msg = "Invalid type: try using Literal[{}] instead?".format(repr(t.literal_value))
            elif t.base_type_name in ('builtins.float', 'builtins.complex'):
                # We special-case warnings for floats and complex numbers.
                msg = "Invalid type: {} literals cannot be used as a type".format(t.simple_name())
            else:
                # And in all other cases, we default to a generic error message.
                # Note: the reason why we use a generic error message for strings
                # but not ints or bools is because whenever we see an out-of-place
                # string, it's unclear if the user meant to construct a literal type
                # or just misspelled a regular type. So we avoid guessing.
                msg = 'Invalid type comment or annotation'

            self.fail(msg, t, code=codes.VALID_TYPE)
            if t.note is not None:
                self.note(t.note, t, code=codes.VALID_TYPE)

        return AnyType(TypeOfAny.from_error, line=t.line, column=t.column)

    def visit_literal_type(self, t: LiteralType) -> Type:
        return t

    def visit_star_type(self, t: StarType) -> Type:
        return StarType(self.anal_type(t.type), t.line)

    def visit_union_type(self, t: UnionType) -> Type:
        if (t.uses_pep604_syntax is True
                and t.is_evaluated is True
                and not self.always_allow_new_syntax
                and not self.options.python_version >= (3, 10)):
            self.fail("X | Y syntax for unions requires Python 3.10", t)
        return UnionType(self.anal_array(t.items), t.line)

    def visit_partial_type(self, t: PartialType) -> Type:
        assert False, "Internal error: Unexpected partial type"

    def visit_ellipsis_type(self, t: EllipsisType) -> Type:
        self.fail('Unexpected "..."', t)
        return AnyType(TypeOfAny.from_error)

    def visit_type_type(self, t: TypeType) -> Type:
        return TypeType.make_normalized(self.anal_type(t.item), line=t.line)

    def visit_placeholder_type(self, t: PlaceholderType) -> Type:
        n = None if t.fullname is None else self.api.lookup_fully_qualified(t.fullname)
        if not n or isinstance(n.node, PlaceholderNode):
            self.api.defer()  # Still incomplete
            return t
        else:
            # TODO: Handle non-TypeInfo
            assert isinstance(n.node, TypeInfo)
            return self.analyze_type_with_type_info(n.node, t.args, t)

    def analyze_callable_args_for_paramspec(
        self,
        callable_args: Type,
        ret_type: Type,
        fallback: Instance,
    ) -> Optional[CallableType]:
        """Construct a 'Callable[P, RET]', where P is ParamSpec, return None if we cannot."""
        if not isinstance(callable_args, UnboundType):
            return None
        sym = self.lookup_qualified(callable_args.name, callable_args)
        if sym is None:
            return None
        tvar_def = self.tvar_scope.get_binding(sym)
        if not isinstance(tvar_def, ParamSpecType):
            return None

        # TODO: Use tuple[...] or Mapping[..] instead?
        obj = self.named_type('builtins.object')
        return CallableType(
            [ParamSpecType(tvar_def.name, tvar_def.fullname, tvar_def.id, ParamSpecFlavor.ARGS,
                           upper_bound=obj),
             ParamSpecType(tvar_def.name, tvar_def.fullname, tvar_def.id, ParamSpecFlavor.KWARGS,
                           upper_bound=obj)],
            [nodes.ARG_STAR, nodes.ARG_STAR2],
            [None, None],
            ret_type=ret_type,
            fallback=fallback,
        )

    def analyze_callable_type(self, t: UnboundType) -> Type:
        fallback = self.named_type('builtins.function')
        if len(t.args) == 0:
            # Callable (bare). Treat as Callable[..., Any].
            any_type = self.get_omitted_any(t)
            ret = callable_with_ellipsis(any_type, any_type, fallback)
        elif len(t.args) == 2:
            callable_args = t.args[0]
            ret_type = t.args[1]
            if isinstance(callable_args, TypeList):
                # Callable[[ARG, ...], RET] (ordinary callable type)
                analyzed_args = self.analyze_callable_args(callable_args)
                if analyzed_args is None:
                    return AnyType(TypeOfAny.from_error)
                args, kinds, names = analyzed_args
                ret = CallableType(args,
                                   kinds,
                                   names,
                                   ret_type=ret_type,
                                   fallback=fallback)
            elif isinstance(callable_args, EllipsisType):
                # Callable[..., RET] (with literal ellipsis; accept arbitrary arguments)
                ret = callable_with_ellipsis(AnyType(TypeOfAny.explicit),
                                             ret_type=ret_type,
                                             fallback=fallback)
            else:
                # Callable[P, RET] (where P is ParamSpec)
                maybe_ret = self.analyze_callable_args_for_paramspec(
                    callable_args,
                    ret_type,
                    fallback
                )
                if maybe_ret is None:
                    # Callable[?, RET] (where ? is something invalid)
                    # TODO(PEP612): change error to mention paramspec, once we actually have some
                    # support for it
                    self.fail('The first argument to Callable must be a list of types or "..."', t)
                    return AnyType(TypeOfAny.from_error)
                ret = maybe_ret
        else:
            self.fail('Please use "Callable[[<parameters>], <return type>]" or "Callable"', t)
            return AnyType(TypeOfAny.from_error)
        assert isinstance(ret, CallableType)
        return ret.accept(self)

    def analyze_callable_args(self, arglist: TypeList) -> Optional[Tuple[List[Type],
                                                                         List[ArgKind],
                                                                         List[Optional[str]]]]:
        args: List[Type] = []
        kinds: List[ArgKind] = []
        names: List[Optional[str]] = []
        for arg in arglist.items:
            if isinstance(arg, CallableArgument):
                args.append(arg.typ)
                names.append(arg.name)
                if arg.constructor is None:
                    return None
                found = self.lookup_qualified(arg.constructor, arg)
                if found is None:
                    # Looking it up already put an error message in
                    return None
                elif found.fullname not in ARG_KINDS_BY_CONSTRUCTOR:
                    self.fail('Invalid argument constructor "{}"'.format(
                        found.fullname), arg)
                    return None
                else:
                    assert found.fullname is not None
                    kind = ARG_KINDS_BY_CONSTRUCTOR[found.fullname]
                    kinds.append(kind)
                    if arg.name is not None and kind.is_star():
                        self.fail("{} arguments should not have names".format(
                            arg.constructor), arg)
                        return None
            else:
                args.append(arg)
                kinds.append(ARG_POS)
                names.append(None)
        # Note that arglist below is only used for error context.
        check_arg_names(names, [arglist] * len(args), self.fail, "Callable")
        check_arg_kinds(kinds, [arglist] * len(args), self.fail)
        return args, kinds, names

    def analyze_literal_type(self, t: UnboundType) -> Type:
        if len(t.args) == 0:
            self.fail('Literal[...] must have at least one parameter', t)
            return AnyType(TypeOfAny.from_error)

        output: List[Type] = []
        for i, arg in enumerate(t.args):
            analyzed_types = self.analyze_literal_param(i + 1, arg, t)
            if analyzed_types is None:
                return AnyType(TypeOfAny.from_error)
            else:
                output.extend(analyzed_types)
        return UnionType.make_union(output, line=t.line)

    def analyze_literal_param(self, idx: int, arg: Type, ctx: Context) -> Optional[List[Type]]:
        # This UnboundType was originally defined as a string.
        if isinstance(arg, UnboundType) and arg.original_str_expr is not None:
            assert arg.original_str_fallback is not None
            return [LiteralType(
                value=arg.original_str_expr,
                fallback=self.named_type_with_normalized_str(arg.original_str_fallback),
                line=arg.line,
                column=arg.column,
            )]

        # If arg is an UnboundType that was *not* originally defined as
        # a string, try expanding it in case it's a type alias or something.
        if isinstance(arg, UnboundType):
            self.nesting_level += 1
            try:
                arg = self.visit_unbound_type(arg, defining_literal=True)
            finally:
                self.nesting_level -= 1

        # Literal[...] cannot contain Any. Give up and add an error message
        # (if we haven't already).
        arg = get_proper_type(arg)
        if isinstance(arg, AnyType):
            # Note: We can encounter Literals containing 'Any' under three circumstances:
            #
            # 1. If the user attempts use an explicit Any as a parameter
            # 2. If the user is trying to use an enum value imported from a module with
            #    no type hints, giving it an implicit type of 'Any'
            # 3. If there's some other underlying problem with the parameter.
            #
            # We report an error in only the first two cases. In the third case, we assume
            # some other region of the code has already reported a more relevant error.
            #
            # TODO: Once we start adding support for enums, make sure we report a custom
            # error for case 2 as well.
            if arg.type_of_any not in (TypeOfAny.from_error, TypeOfAny.special_form):
                self.fail('Parameter {} of Literal[...] cannot be of type "Any"'.format(idx), ctx)
            return None
        elif isinstance(arg, RawExpressionType):
            # A raw literal. Convert it directly into a literal if we can.
            if arg.literal_value is None:
                name = arg.simple_name()
                if name in ('float', 'complex'):
                    msg = 'Parameter {} of Literal[...] cannot be of type "{}"'.format(idx, name)
                else:
                    msg = 'Invalid type: Literal[...] cannot contain arbitrary expressions'
                self.fail(msg, ctx)
                # Note: we deliberately ignore arg.note here: the extra info might normally be
                # helpful, but it generally won't make sense in the context of a Literal[...].
                return None

            # Remap bytes and unicode into the appropriate type for the correct Python version
            fallback = self.named_type_with_normalized_str(arg.base_type_name)
            assert isinstance(fallback, Instance)
            return [LiteralType(arg.literal_value, fallback, line=arg.line, column=arg.column)]
        elif isinstance(arg, (NoneType, LiteralType)):
            # Types that we can just add directly to the literal/potential union of literals.
            return [arg]
        elif isinstance(arg, Instance) and arg.last_known_value is not None:
            # Types generated from declarations like "var: Final = 4".
            return [arg.last_known_value]
        elif isinstance(arg, UnionType):
            out = []
            for union_arg in arg.items:
                union_result = self.analyze_literal_param(idx, union_arg, ctx)
                if union_result is None:
                    return None
                out.extend(union_result)
            return out
        else:
            self.fail('Parameter {} of Literal[...] is invalid'.format(idx), ctx)
            return None

    def analyze_type(self, t: Type) -> Type:
        return t.accept(self)

    def fail(self, msg: str, ctx: Context, *, code: Optional[ErrorCode] = None) -> None:
        self.fail_func(msg, ctx, code=code)

    def note(self, msg: str, ctx: Context, *, code: Optional[ErrorCode] = None) -> None:
        self.note_func(msg, ctx, code=code)

    @contextmanager
    def tvar_scope_frame(self) -> Iterator[None]:
        old_scope = self.tvar_scope
        self.tvar_scope = self.tvar_scope.method_frame()
        yield
        self.tvar_scope = old_scope

    def infer_type_variables(self,
                             type: CallableType) -> List[Tuple[str, TypeVarLikeExpr]]:
        """Return list of unique type variables referred to in a callable."""
        names: List[str] = []
        tvars: List[TypeVarLikeExpr] = []
        for arg in type.arg_types:
            for name, tvar_expr in arg.accept(
                TypeVarLikeQuery(self.lookup_qualified, self.tvar_scope)
            ):
                if name not in names:
                    names.append(name)
                    tvars.append(tvar_expr)
        # When finding type variables in the return type of a function, don't
        # look inside Callable types.  Type variables only appearing in
        # functions in the return type belong to those functions, not the
        # function we're currently analyzing.
        for name, tvar_expr in type.ret_type.accept(
            TypeVarLikeQuery(self.lookup_qualified, self.tvar_scope, include_callables=False)
        ):
            if name not in names:
                names.append(name)
                tvars.append(tvar_expr)
        return list(zip(names, tvars))

    def bind_function_type_variables(
        self, fun_type: CallableType, defn: Context
    ) -> Sequence[TypeVarLikeType]:
        """Find the type variables of the function type and bind them in our tvar_scope"""
        if fun_type.variables:
            for var in fun_type.variables:
                var_node = self.lookup_qualified(var.name, defn)
                assert var_node, "Binding for function type variable not found within function"
                var_expr = var_node.node
                assert isinstance(var_expr, TypeVarLikeExpr)
                self.tvar_scope.bind_new(var.name, var_expr)
            return fun_type.variables
        typevars = self.infer_type_variables(fun_type)
        # Do not define a new type variable if already defined in scope.
        typevars = [(name, tvar) for name, tvar in typevars
                    if not self.is_defined_type_var(name, defn)]
        defs: List[TypeVarLikeType] = []
        for name, tvar in typevars:
            if not self.tvar_scope.allow_binding(tvar.fullname):
                self.fail('Type variable "{}" is bound by an outer class'.format(name), defn)
            self.tvar_scope.bind_new(name, tvar)
            binding = self.tvar_scope.get_binding(tvar.fullname)
            assert binding is not None
            defs.append(binding)

        return defs

    def is_defined_type_var(self, tvar: str, context: Context) -> bool:
        tvar_node = self.lookup_qualified(tvar, context)
        if not tvar_node:
            return False
        return self.tvar_scope.get_binding(tvar_node) is not None

    def anal_array(self,
                   a: Iterable[Type],
                   nested: bool = True, *,
                   allow_param_spec: bool = False) -> List[Type]:
        res: List[Type] = []
        for t in a:
            res.append(self.anal_type(t, nested, allow_param_spec=allow_param_spec))
        return res

    def anal_type(self, t: Type, nested: bool = True, *, allow_param_spec: bool = False) -> Type:
        if nested:
            self.nesting_level += 1
        old_allow_required = self.allow_required
        self.allow_required = False
        try:
            analyzed = t.accept(self)
        finally:
            if nested:
                self.nesting_level -= 1
            self.allow_required = old_allow_required
        if (not allow_param_spec
                and isinstance(analyzed, ParamSpecType)
                and analyzed.flavor == ParamSpecFlavor.BARE):
            self.fail('Invalid location for ParamSpec "{}"'.format(analyzed.name), t)
            self.note(
                'You can use ParamSpec as the first argument to Callable, e.g., '
                "'Callable[{}, int]'".format(analyzed.name),
                t
            )
        return analyzed

    def anal_var_def(self, var_def: TypeVarLikeType) -> TypeVarLikeType:
        if isinstance(var_def, TypeVarType):
            return TypeVarType(
                var_def.name,
                var_def.fullname,
                var_def.id.raw_id,
                self.anal_array(var_def.values),
                var_def.upper_bound.accept(self),
                var_def.variance,
                var_def.line
            )
        else:
            return var_def

    def anal_var_defs(self, var_defs: Sequence[TypeVarLikeType]) -> List[TypeVarLikeType]:
        return [self.anal_var_def(vd) for vd in var_defs]

    def named_type_with_normalized_str(self, fully_qualified_name: str) -> Instance:
        """Does almost the same thing as `named_type`, except that we immediately
        unalias `builtins.bytes` and `builtins.unicode` to `builtins.str` as appropriate.
        """
        python_version = self.options.python_version
        if python_version[0] == 2 and fully_qualified_name == 'builtins.bytes':
            fully_qualified_name = 'builtins.str'
        if python_version[0] >= 3 and fully_qualified_name == 'builtins.unicode':
            fully_qualified_name = 'builtins.str'
        return self.named_type(fully_qualified_name)

    def named_type(self, fully_qualified_name: str,
                   args: Optional[List[Type]] = None,
                   line: int = -1,
                   column: int = -1) -> Instance:
        node = self.lookup_fqn_func(fully_qualified_name)
        assert isinstance(node.node, TypeInfo)
        any_type = AnyType(TypeOfAny.special_form)
        return Instance(node.node, args or [any_type] * len(node.node.defn.type_vars),
                        line=line, column=column)

    def tuple_type(self, items: List[Type]) -> TupleType:
        any_type = AnyType(TypeOfAny.special_form)
        return TupleType(items, fallback=self.named_type('builtins.tuple', [any_type]))


TypeVarLikeList = List[Tuple[str, TypeVarLikeExpr]]

# Mypyc doesn't support callback protocols yet.
MsgCallback = Callable[[str, Context, DefaultNamedArg(Optional[ErrorCode], 'code')], None]


def get_omitted_any(disallow_any: bool, fail: MsgCallback, note: MsgCallback,
                    orig_type: Type, python_version: Tuple[int, int],
                    fullname: Optional[str] = None,
                    unexpanded_type: Optional[Type] = None) -> AnyType:
    if disallow_any:
        nongen_builtins = get_nongen_builtins(python_version)
        if fullname in nongen_builtins:
            typ = orig_type
            # We use a dedicated error message for builtin generics (as the most common case).
            alternative = nongen_builtins[fullname]
            fail(message_registry.IMPLICIT_GENERIC_ANY_BUILTIN.format(alternative), typ,
                 code=codes.TYPE_ARG)
        else:
            typ = unexpanded_type or orig_type
            type_str = typ.name if isinstance(typ, UnboundType) else format_type_bare(typ)

            fail(
                message_registry.BARE_GENERIC.format(quote_type_string(type_str)),
                typ,
                code=codes.TYPE_ARG)
            base_type = get_proper_type(orig_type)
            base_fullname = (
                base_type.type.fullname if isinstance(base_type, Instance) else fullname
            )
            # Ideally, we'd check whether the type is quoted or `from __future__ annotations`
            # is set before issuing this note
            if python_version < (3, 9) and base_fullname in GENERIC_STUB_NOT_AT_RUNTIME_TYPES:
                # Recommend `from __future__ import annotations` or to put type in quotes
                # (string literal escaping) for classes not generic at runtime
                note(
                    "Subscripting classes that are not generic at runtime may require "
                    "escaping, see https://mypy.readthedocs.io/en/stable/runtime_troubles.html"
                    "#not-generic-runtime",
                    typ,
                    code=codes.TYPE_ARG)

        any_type = AnyType(TypeOfAny.from_error, line=typ.line, column=typ.column)
    else:
        any_type = AnyType(
            TypeOfAny.from_omitted_generics, line=orig_type.line, column=orig_type.column
        )
    return any_type


def fix_instance(t: Instance, fail: MsgCallback, note: MsgCallback,
                 disallow_any: bool, python_version: Tuple[int, int],
                 use_generic_error: bool = False,
                 unexpanded_type: Optional[Type] = None,) -> None:
    """Fix a malformed instance by replacing all type arguments with Any.

    Also emit a suitable error if this is not due to implicit Any's.
    """
    if len(t.args) == 0:
        if use_generic_error:
            fullname: Optional[str] = None
        else:
            fullname = t.type.fullname
        any_type = get_omitted_any(disallow_any, fail, note, t, python_version, fullname,
                                   unexpanded_type)
        t.args = (any_type,) * len(t.type.type_vars)
        return
    # Invalid number of type parameters.
    n = len(t.type.type_vars)
    s = '{} type arguments'.format(n)
    if n == 0:
        s = 'no type arguments'
    elif n == 1:
        s = '1 type argument'
    act = str(len(t.args))
    if act == '0':
        act = 'none'
    fail('"{}" expects {}, but {} given'.format(
        t.type.name, s, act), t, code=codes.TYPE_ARG)
    # Construct the correct number of type arguments, as
    # otherwise the type checker may crash as it expects
    # things to be right.
    t.args = tuple(AnyType(TypeOfAny.from_error) for _ in t.type.type_vars)
    t.invalid = True


def expand_type_alias(node: TypeAlias, args: List[Type],
                      fail: MsgCallback, no_args: bool, ctx: Context, *,
                      unexpanded_type: Optional[Type] = None,
                      disallow_any: bool = False) -> Type:
    """Expand a (generic) type alias target following the rules outlined in TypeAlias docstring.

    Here:
        target: original target type (contains unbound type variables)
        alias_tvars: type variable names
        args: types to be substituted in place of type variables
        fail: error reporter callback
        no_args: whether original definition used a bare generic `A = List`
        ctx: context where expansion happens
    """
    exp_len = len(node.alias_tvars)
    act_len = len(args)
    if exp_len > 0 and act_len == 0:
        # Interpret bare Alias same as normal generic, i.e., Alias[Any, Any, ...]
        return set_any_tvars(node, ctx.line, ctx.column,
                             disallow_any=disallow_any, fail=fail,
                             unexpanded_type=unexpanded_type)
    if exp_len == 0 and act_len == 0:
        if no_args:
            assert isinstance(node.target, Instance)  # type: ignore[misc]
            # Note: this is the only case where we use an eager expansion. See more info about
            # no_args aliases like L = List in the docstring for TypeAlias class.
            return Instance(node.target.type, [], line=ctx.line, column=ctx.column)
        return TypeAliasType(node, [], line=ctx.line, column=ctx.column)
    if (exp_len == 0 and act_len > 0
            and isinstance(node.target, Instance)  # type: ignore[misc]
            and no_args):
        tp = Instance(node.target.type, args)
        tp.line = ctx.line
        tp.column = ctx.column
        return tp
    if act_len != exp_len:
        fail('Bad number of arguments for type alias, expected: %s, given: %s'
             % (exp_len, act_len), ctx)
        return set_any_tvars(node, ctx.line, ctx.column, from_error=True)
    typ = TypeAliasType(node, args, ctx.line, ctx.column)
    assert typ.alias is not None
    # HACK: Implement FlexibleAlias[T, typ] by expanding it to typ here.
    if (isinstance(typ.alias.target, Instance)  # type: ignore
            and typ.alias.target.type.fullname == 'mypy_extensions.FlexibleAlias'):
        exp = get_proper_type(typ)
        assert isinstance(exp, Instance)
        return exp.args[-1]
    return typ


def set_any_tvars(node: TypeAlias,
                  newline: int, newcolumn: int, *,
                  from_error: bool = False,
                  disallow_any: bool = False,
                  fail: Optional[MsgCallback] = None,
                  unexpanded_type: Optional[Type] = None) -> Type:
    if from_error or disallow_any:
        type_of_any = TypeOfAny.from_error
    else:
        type_of_any = TypeOfAny.from_omitted_generics
    if disallow_any:
        assert fail is not None
        otype = unexpanded_type or node.target
        type_str = otype.name if isinstance(otype, UnboundType) else format_type_bare(otype)

        fail(message_registry.BARE_GENERIC.format(quote_type_string(type_str)),
             Context(newline, newcolumn), code=codes.TYPE_ARG)
    any_type = AnyType(type_of_any, line=newline, column=newcolumn)
    return TypeAliasType(node, [any_type] * len(node.alias_tvars), newline, newcolumn)


def remove_dups(tvars: Iterable[T]) -> List[T]:
    # Get unique elements in order of appearance
    all_tvars: Set[T] = set()
    new_tvars: List[T] = []
    for t in tvars:
        if t not in all_tvars:
            new_tvars.append(t)
            all_tvars.add(t)
    return new_tvars


def flatten_tvars(ll: Iterable[List[T]]) -> List[T]:
    return remove_dups(chain.from_iterable(ll))


class TypeVarLikeQuery(TypeQuery[TypeVarLikeList]):
    """Find TypeVar and ParamSpec references in an unbound type."""

    def __init__(self,
                 lookup: Callable[[str, Context], Optional[SymbolTableNode]],
                 scope: 'TypeVarLikeScope',
                 *,
                 include_callables: bool = True,
                 include_bound_tvars: bool = False) -> None:
        self.include_callables = include_callables
        self.lookup = lookup
        self.scope = scope
        self.include_bound_tvars = include_bound_tvars
        super().__init__(flatten_tvars)

    def _seems_like_callable(self, type: UnboundType) -> bool:
        if not type.args:
            return False
        if isinstance(type.args[0], (EllipsisType, TypeList)):
            return True
        return False

    def visit_unbound_type(self, t: UnboundType) -> TypeVarLikeList:
        name = t.name
        node = None
        # Special case P.args and P.kwargs for ParamSpecs only.
        if name.endswith('args'):
            if name.endswith('.args') or name.endswith('.kwargs'):
                base = '.'.join(name.split('.')[:-1])
                n = self.lookup(base, t)
                if n is not None and isinstance(n.node, ParamSpecExpr):
                    node = n
                    name = base
        if node is None:
            node = self.lookup(name, t)
        if node and isinstance(node.node, TypeVarLikeExpr) and (
                self.include_bound_tvars or self.scope.get_binding(node) is None):
            assert isinstance(node.node, TypeVarLikeExpr)
            return [(name, node.node)]
        elif not self.include_callables and self._seems_like_callable(t):
            return []
        elif node and node.fullname in ('typing_extensions.Literal', 'typing.Literal'):
            return []
        elif node and node.fullname in ('typing_extensions.Annotated', 'typing.Annotated'):
            # Don't query the second argument to Annotated for TypeVars
            return self.query_types([t.args[0]])
        else:
            return super().visit_unbound_type(t)

    def visit_callable_type(self, t: CallableType) -> TypeVarLikeList:
        if self.include_callables:
            return super().visit_callable_type(t)
        else:
            return []


def check_for_explicit_any(typ: Optional[Type],
                           options: Options,
                           is_typeshed_stub: bool,
                           msg: MessageBuilder,
                           context: Context) -> None:
    if (options.disallow_any_explicit and
            not is_typeshed_stub and
            typ and
            has_explicit_any(typ)):
        msg.explicit_any(context)


def has_explicit_any(t: Type) -> bool:
    """
    Whether this type is or type it contains is an Any coming from explicit type annotation
    """
    return t.accept(HasExplicitAny())


class HasExplicitAny(TypeQuery[bool]):
    def __init__(self) -> None:
        super().__init__(any)

    def visit_any(self, t: AnyType) -> bool:
        return t.type_of_any == TypeOfAny.explicit

    def visit_typeddict_type(self, t: TypedDictType) -> bool:
        # typeddict is checked during TypedDict declaration, so don't typecheck it here.
        return False


def has_any_from_unimported_type(t: Type) -> bool:
    """Return true if this type is Any because an import was not followed.

    If type t is such Any type or has type arguments that contain such Any type
    this function will return true.
    """
    return t.accept(HasAnyFromUnimportedType())


class HasAnyFromUnimportedType(TypeQuery[bool]):
    def __init__(self) -> None:
        super().__init__(any)

    def visit_any(self, t: AnyType) -> bool:
        return t.type_of_any == TypeOfAny.from_unimported_type

    def visit_typeddict_type(self, t: TypedDictType) -> bool:
        # typeddict is checked during TypedDict declaration, so don't typecheck it here
        return False


def collect_any_types(t: Type) -> List[AnyType]:
    """Return all inner `AnyType`s of type t"""
    return t.accept(CollectAnyTypesQuery())


class CollectAnyTypesQuery(TypeQuery[List[AnyType]]):
    def __init__(self) -> None:
        super().__init__(self.combine_lists_strategy)

    def visit_any(self, t: AnyType) -> List[AnyType]:
        return [t]

    @classmethod
    def combine_lists_strategy(cls, it: Iterable[List[AnyType]]) -> List[AnyType]:
        result: List[AnyType] = []
        for l in it:
            result.extend(l)
        return result


def collect_all_inner_types(t: Type) -> List[Type]:
    """
    Return all types that `t` contains
    """
    return t.accept(CollectAllInnerTypesQuery())


class CollectAllInnerTypesQuery(TypeQuery[List[Type]]):
    def __init__(self) -> None:
        super().__init__(self.combine_lists_strategy)

    def query_types(self, types: Iterable[Type]) -> List[Type]:
        return self.strategy([t.accept(self) for t in types]) + list(types)

    @classmethod
    def combine_lists_strategy(cls, it: Iterable[List[Type]]) -> List[Type]:
        return list(itertools.chain.from_iterable(it))


def make_optional_type(t: Type) -> Type:
    """Return the type corresponding to Optional[t].

    Note that we can't use normal union simplification, since this function
    is called during semantic analysis and simplification only works during
    type checking.
    """
    t = get_proper_type(t)
    if isinstance(t, NoneType):
        return t
    elif isinstance(t, UnionType):
        items = [item for item in union_items(t)
                 if not isinstance(item, NoneType)]
        return UnionType(items + [NoneType()], t.line, t.column)
    else:
        return UnionType([t, NoneType()], t.line, t.column)


def fix_instance_types(t: Type, fail: MsgCallback, note: MsgCallback,
                       python_version: Tuple[int, int]) -> None:
    """Recursively fix all instance types (type argument count) in a given type.

    For example 'Union[Dict, List[str, int]]' will be transformed into
    'Union[Dict[Any, Any], List[Any]]' in place.
    """
    t.accept(InstanceFixer(fail, note, python_version))


class InstanceFixer(TypeTraverserVisitor):
    def __init__(
        self, fail: MsgCallback, note: MsgCallback, python_version: Tuple[int, int]
    ) -> None:
        self.fail = fail
        self.note = note
        self.python_version = python_version

    def visit_instance(self, typ: Instance) -> None:
        super().visit_instance(typ)
        if len(typ.args) != len(typ.type.type_vars):
            fix_instance(typ, self.fail, self.note, disallow_any=False,
                         python_version=self.python_version, use_generic_error=True)
