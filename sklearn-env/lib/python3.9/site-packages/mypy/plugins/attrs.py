"""Plugin for supporting the attrs library (http://www.attrs.org)"""

from mypy.backports import OrderedDict

from typing import Optional, Dict, List, cast, Tuple, Iterable
from typing_extensions import Final

import mypy.plugin  # To avoid circular imports.
from mypy.exprtotype import expr_to_unanalyzed_type, TypeTranslationError
from mypy.lookup import lookup_fully_qualified
from mypy.nodes import (
    Context, Argument, Var, ARG_OPT, ARG_POS, TypeInfo, AssignmentStmt,
    TupleExpr, ListExpr, NameExpr, CallExpr, RefExpr, FuncDef,
    is_class_var, TempNode, Decorator, MemberExpr, Expression,
    SymbolTableNode, MDEF, JsonDict, OverloadedFuncDef, ARG_NAMED_OPT, ARG_NAMED,
    TypeVarExpr, PlaceholderNode
)
from mypy.plugin import SemanticAnalyzerPluginInterface
from mypy.plugins.common import (
    _get_argument, _get_bool_argument, _get_decorator_bool_argument, add_method,
    deserialize_and_fixup_type
)
from mypy.types import (
    TupleType, Type, AnyType, TypeOfAny, CallableType, NoneType, TypeVarType,
    Overloaded, UnionType, FunctionLike, get_proper_type,
)
from mypy.typeops import make_simplified_union, map_type_from_supertype
from mypy.typevars import fill_typevars
from mypy.util import unmangle
from mypy.server.trigger import make_wildcard_trigger

KW_ONLY_PYTHON_2_UNSUPPORTED = "kw_only is not supported in Python 2"

# The names of the different functions that create classes or arguments.
attr_class_makers: Final = {
    'attr.s',
    'attr.attrs',
    'attr.attributes',
}
attr_dataclass_makers: Final = {
    'attr.dataclass',
}
attr_frozen_makers: Final = {"attr.frozen"}
attr_define_makers: Final = {"attr.define", "attr.mutable"}
attr_attrib_makers: Final = {
    'attr.ib',
    'attr.attrib',
    'attr.attr',
    'attr.field',
}

SELF_TVAR_NAME: Final = "_AT"


class Converter:
    """Holds information about a `converter=` argument"""

    def __init__(self,
                 name: Optional[str] = None,
                 is_attr_converters_optional: bool = False) -> None:
        self.name = name
        self.is_attr_converters_optional = is_attr_converters_optional


class Attribute:
    """The value of an attr.ib() call."""

    def __init__(self, name: str, info: TypeInfo,
                 has_default: bool, init: bool, kw_only: bool, converter: Converter,
                 context: Context,
                 init_type: Optional[Type]) -> None:
        self.name = name
        self.info = info
        self.has_default = has_default
        self.init = init
        self.kw_only = kw_only
        self.converter = converter
        self.context = context
        self.init_type = init_type

    def argument(self, ctx: 'mypy.plugin.ClassDefContext') -> Argument:
        """Return this attribute as an argument to __init__."""
        assert self.init

        init_type = self.init_type or self.info[self.name].type

        if self.converter.name:
            # When a converter is set the init_type is overridden by the first argument
            # of the converter method.
            converter = lookup_fully_qualified(self.converter.name, ctx.api.modules,
                                               raise_on_missing=False)
            if not converter:
                # The converter may be a local variable. Check there too.
                converter = ctx.api.lookup_qualified(self.converter.name, self.info, True)

            # Get the type of the converter.
            converter_type: Optional[Type] = None
            if converter and isinstance(converter.node, TypeInfo):
                from mypy.checkmember import type_object_type  # To avoid import cycle.
                converter_type = type_object_type(converter.node, ctx.api.named_type)
            elif converter and isinstance(converter.node, OverloadedFuncDef):
                converter_type = converter.node.type
            elif converter and converter.type:
                converter_type = converter.type

            init_type = None
            converter_type = get_proper_type(converter_type)
            if isinstance(converter_type, CallableType) and converter_type.arg_types:
                init_type = ctx.api.anal_type(converter_type.arg_types[0])
            elif isinstance(converter_type, Overloaded):
                types: List[Type] = []
                for item in converter_type.items:
                    # Walk the overloads looking for methods that can accept one argument.
                    num_arg_types = len(item.arg_types)
                    if not num_arg_types:
                        continue
                    if num_arg_types > 1 and any(kind == ARG_POS for kind in item.arg_kinds[1:]):
                        continue
                    types.append(item.arg_types[0])
                # Make a union of all the valid types.
                if types:
                    args = make_simplified_union(types)
                    init_type = ctx.api.anal_type(args)

            if self.converter.is_attr_converters_optional and init_type:
                # If the converter was attr.converter.optional(type) then add None to
                # the allowed init_type.
                init_type = UnionType.make_union([init_type, NoneType()])

            if not init_type:
                ctx.api.fail("Cannot determine __init__ type from converter", self.context)
                init_type = AnyType(TypeOfAny.from_error)
        elif self.converter.name == '':
            # This means we had a converter but it's not of a type we can infer.
            # Error was shown in _get_converter_name
            init_type = AnyType(TypeOfAny.from_error)

        if init_type is None:
            if ctx.api.options.disallow_untyped_defs:
                # This is a compromise.  If you don't have a type here then the
                # __init__ will be untyped. But since the __init__ is added it's
                # pointing at the decorator. So instead we also show the error in the
                # assignment, which is where you would fix the issue.
                node = self.info[self.name].node
                assert node is not None
                ctx.api.msg.need_annotation_for_var(node, self.context)

            # Convert type not set to Any.
            init_type = AnyType(TypeOfAny.unannotated)

        if self.kw_only:
            arg_kind = ARG_NAMED_OPT if self.has_default else ARG_NAMED
        else:
            arg_kind = ARG_OPT if self.has_default else ARG_POS

        # Attrs removes leading underscores when creating the __init__ arguments.
        return Argument(Var(self.name.lstrip("_"), init_type), init_type,
                        None,
                        arg_kind)

    def serialize(self) -> JsonDict:
        """Serialize this object so it can be saved and restored."""
        return {
            'name': self.name,
            'has_default': self.has_default,
            'init': self.init,
            'kw_only': self.kw_only,
            'converter_name': self.converter.name,
            'converter_is_attr_converters_optional': self.converter.is_attr_converters_optional,
            'context_line': self.context.line,
            'context_column': self.context.column,
            'init_type': self.init_type.serialize() if self.init_type else None,
        }

    @classmethod
    def deserialize(cls, info: TypeInfo,
                    data: JsonDict,
                    api: SemanticAnalyzerPluginInterface) -> 'Attribute':
        """Return the Attribute that was serialized."""
        raw_init_type = data['init_type']
        init_type = deserialize_and_fixup_type(raw_init_type, api) if raw_init_type else None

        return Attribute(data['name'],
            info,
            data['has_default'],
            data['init'],
            data['kw_only'],
            Converter(data['converter_name'], data['converter_is_attr_converters_optional']),
            Context(line=data['context_line'], column=data['context_column']),
            init_type)

    def expand_typevar_from_subtype(self, sub_type: TypeInfo) -> None:
        """Expands type vars in the context of a subtype when an attribute is inherited
        from a generic super type."""
        if not isinstance(self.init_type, TypeVarType):
            return

        self.init_type = map_type_from_supertype(self.init_type, sub_type, self.info)


def _determine_eq_order(ctx: 'mypy.plugin.ClassDefContext') -> bool:
    """
    Validate the combination of *cmp*, *eq*, and *order*. Derive the effective
    value of order.
    """
    cmp = _get_decorator_optional_bool_argument(ctx, 'cmp')
    eq = _get_decorator_optional_bool_argument(ctx, 'eq')
    order = _get_decorator_optional_bool_argument(ctx, 'order')

    if cmp is not None and any((eq is not None, order is not None)):
        ctx.api.fail('Don\'t mix "cmp" with "eq" and "order"', ctx.reason)

    # cmp takes precedence due to bw-compatibility.
    if cmp is not None:
        return cmp

    # If left None, equality is on and ordering mirrors equality.
    if eq is None:
        eq = True

    if order is None:
        order = eq

    if eq is False and order is True:
        ctx.api.fail('eq must be True if order is True', ctx.reason)

    return order


def _get_decorator_optional_bool_argument(
    ctx: 'mypy.plugin.ClassDefContext',
    name: str,
    default: Optional[bool] = None,
) -> Optional[bool]:
    """Return the Optional[bool] argument for the decorator.

    This handles both @decorator(...) and @decorator.
    """
    if isinstance(ctx.reason, CallExpr):
        attr_value = _get_argument(ctx.reason, name)
        if attr_value:
            if isinstance(attr_value, NameExpr):
                if attr_value.fullname == 'builtins.True':
                    return True
                if attr_value.fullname == 'builtins.False':
                    return False
                if attr_value.fullname == 'builtins.None':
                    return None
            ctx.api.fail('"{}" argument must be True or False.'.format(name), ctx.reason)
            return default
        return default
    else:
        return default


def attr_class_maker_callback(ctx: 'mypy.plugin.ClassDefContext',
                              auto_attribs_default: Optional[bool] = False,
                              frozen_default: bool = False) -> None:
    """Add necessary dunder methods to classes decorated with attr.s.

    attrs is a package that lets you define classes without writing dull boilerplate code.

    At a quick glance, the decorator searches the class body for assignments of `attr.ib`s (or
    annotated variables if auto_attribs=True), then depending on how the decorator is called,
    it will add an __init__ or all the __cmp__ methods.  For frozen=True it will turn the attrs
    into properties.

    See http://www.attrs.org/en/stable/how-does-it-work.html for information on how attrs works.
    """
    info = ctx.cls.info

    init = _get_decorator_bool_argument(ctx, 'init', True)
    frozen = _get_frozen(ctx, frozen_default)
    order = _determine_eq_order(ctx)
    slots = _get_decorator_bool_argument(ctx, 'slots', False)

    auto_attribs = _get_decorator_optional_bool_argument(ctx, 'auto_attribs', auto_attribs_default)
    kw_only = _get_decorator_bool_argument(ctx, 'kw_only', False)

    if ctx.api.options.python_version[0] < 3:
        if auto_attribs:
            ctx.api.fail("auto_attribs is not supported in Python 2", ctx.reason)
            return
        if not info.defn.base_type_exprs:
            # Note: This will not catch subclassing old-style classes.
            ctx.api.fail("attrs only works with new-style classes", info.defn)
            return
        if kw_only:
            ctx.api.fail(KW_ONLY_PYTHON_2_UNSUPPORTED, ctx.reason)
            return

    attributes = _analyze_class(ctx, auto_attribs, kw_only)

    # Check if attribute types are ready.
    for attr in attributes:
        node = info.get(attr.name)
        if node is None:
            # This name is likely blocked by a star import. We don't need to defer because
            # defer() is already called by mark_incomplete().
            return
        if node.type is None and not ctx.api.final_iteration:
            ctx.api.defer()
            return

    _add_attrs_magic_attribute(ctx, raw_attr_types=[info[attr.name].type for attr in attributes])
    if slots:
        _add_slots(ctx, attributes)

    # Save the attributes so that subclasses can reuse them.
    ctx.cls.info.metadata['attrs'] = {
        'attributes': [attr.serialize() for attr in attributes],
        'frozen': frozen,
    }

    adder = MethodAdder(ctx)
    if init:
        _add_init(ctx, attributes, adder)
    if order:
        _add_order(ctx, adder)
    if frozen:
        _make_frozen(ctx, attributes)


def _get_frozen(ctx: 'mypy.plugin.ClassDefContext', frozen_default: bool) -> bool:
    """Return whether this class is frozen."""
    if _get_decorator_bool_argument(ctx, 'frozen', frozen_default):
        return True
    # Subclasses of frozen classes are frozen so check that.
    for super_info in ctx.cls.info.mro[1:-1]:
        if 'attrs' in super_info.metadata and super_info.metadata['attrs']['frozen']:
            return True
    return False


def _analyze_class(ctx: 'mypy.plugin.ClassDefContext',
                   auto_attribs: Optional[bool],
                   kw_only: bool) -> List[Attribute]:
    """Analyze the class body of an attr maker, its parents, and return the Attributes found.

    auto_attribs=True means we'll generate attributes from type annotations also.
    auto_attribs=None means we'll detect which mode to use.
    kw_only=True means that all attributes created here will be keyword only args in __init__.
    """
    own_attrs: OrderedDict[str, Attribute] = OrderedDict()
    if auto_attribs is None:
        auto_attribs = _detect_auto_attribs(ctx)

    # Walk the body looking for assignments and decorators.
    for stmt in ctx.cls.defs.body:
        if isinstance(stmt, AssignmentStmt):
            for attr in _attributes_from_assignment(ctx, stmt, auto_attribs, kw_only):
                # When attrs are defined twice in the same body we want to use the 2nd definition
                # in the 2nd location. So remove it from the OrderedDict.
                # Unless it's auto_attribs in which case we want the 2nd definition in the
                # 1st location.
                if not auto_attribs and attr.name in own_attrs:
                    del own_attrs[attr.name]
                own_attrs[attr.name] = attr
        elif isinstance(stmt, Decorator):
            _cleanup_decorator(stmt, own_attrs)

    for attribute in own_attrs.values():
        # Even though these look like class level assignments we want them to look like
        # instance level assignments.
        if attribute.name in ctx.cls.info.names:
            node = ctx.cls.info.names[attribute.name].node
            if isinstance(node, PlaceholderNode):
                # This node is not ready yet.
                continue
            assert isinstance(node, Var)
            node.is_initialized_in_class = False

    # Traverse the MRO and collect attributes from the parents.
    taken_attr_names = set(own_attrs)
    super_attrs = []
    for super_info in ctx.cls.info.mro[1:-1]:
        if 'attrs' in super_info.metadata:
            # Each class depends on the set of attributes in its attrs ancestors.
            ctx.api.add_plugin_dependency(make_wildcard_trigger(super_info.fullname))

            for data in super_info.metadata['attrs']['attributes']:
                # Only add an attribute if it hasn't been defined before.  This
                # allows for overwriting attribute definitions by subclassing.
                if data['name'] not in taken_attr_names:
                    a = Attribute.deserialize(super_info, data, ctx.api)
                    a.expand_typevar_from_subtype(ctx.cls.info)
                    super_attrs.append(a)
                    taken_attr_names.add(a.name)
    attributes = super_attrs + list(own_attrs.values())

    # Check the init args for correct default-ness.  Note: This has to be done after all the
    # attributes for all classes have been read, because subclasses can override parents.
    last_default = False

    for i, attribute in enumerate(attributes):
        if not attribute.init:
            continue

        if attribute.kw_only:
            # Keyword-only attributes don't care whether they are default or not.
            continue

        # If the issue comes from merging different classes, report it
        # at the class definition point.
        context = attribute.context if i >= len(super_attrs) else ctx.cls

        if not attribute.has_default and last_default:
            ctx.api.fail(
                "Non-default attributes not allowed after default attributes.",
                context)
        last_default |= attribute.has_default

    return attributes


def _detect_auto_attribs(ctx: 'mypy.plugin.ClassDefContext') -> bool:
    """Return whether auto_attribs should be enabled or disabled.

    It's disabled if there are any unannotated attribs()
    """
    for stmt in ctx.cls.defs.body:
        if isinstance(stmt, AssignmentStmt):
            for lvalue in stmt.lvalues:
                lvalues, rvalues = _parse_assignments(lvalue, stmt)

                if len(lvalues) != len(rvalues):
                    # This means we have some assignment that isn't 1 to 1.
                    # It can't be an attrib.
                    continue

                for lhs, rvalue in zip(lvalues, rvalues):
                    # Check if the right hand side is a call to an attribute maker.
                    if (isinstance(rvalue, CallExpr)
                            and isinstance(rvalue.callee, RefExpr)
                            and rvalue.callee.fullname in attr_attrib_makers
                            and not stmt.new_syntax):
                        # This means we have an attrib without an annotation and so
                        # we can't do auto_attribs=True
                        return False
    return True


def _attributes_from_assignment(ctx: 'mypy.plugin.ClassDefContext',
                                stmt: AssignmentStmt, auto_attribs: bool,
                                kw_only: bool) -> Iterable[Attribute]:
    """Return Attribute objects that are created by this assignment.

    The assignments can look like this:
        x = attr.ib()
        x = y = attr.ib()
        x, y = attr.ib(), attr.ib()
    or if auto_attribs is enabled also like this:
        x: type
        x: type = default_value
    """
    for lvalue in stmt.lvalues:
        lvalues, rvalues = _parse_assignments(lvalue, stmt)

        if len(lvalues) != len(rvalues):
            # This means we have some assignment that isn't 1 to 1.
            # It can't be an attrib.
            continue

        for lhs, rvalue in zip(lvalues, rvalues):
            # Check if the right hand side is a call to an attribute maker.
            if (isinstance(rvalue, CallExpr)
                    and isinstance(rvalue.callee, RefExpr)
                    and rvalue.callee.fullname in attr_attrib_makers):
                attr = _attribute_from_attrib_maker(ctx, auto_attribs, kw_only, lhs, rvalue, stmt)
                if attr:
                    yield attr
            elif auto_attribs and stmt.type and stmt.new_syntax and not is_class_var(lhs):
                yield _attribute_from_auto_attrib(ctx, kw_only, lhs, rvalue, stmt)


def _cleanup_decorator(stmt: Decorator, attr_map: Dict[str, Attribute]) -> None:
    """Handle decorators in class bodies.

    `x.default` will set a default value on x
    `x.validator` and `x.default` will get removed to avoid throwing a type error.
    """
    remove_me = []
    for func_decorator in stmt.decorators:
        if (isinstance(func_decorator, MemberExpr)
                and isinstance(func_decorator.expr, NameExpr)
                and func_decorator.expr.name in attr_map):

            if func_decorator.name == 'default':
                attr_map[func_decorator.expr.name].has_default = True

            if func_decorator.name in ('default', 'validator'):
                # These are decorators on the attrib object that only exist during
                # class creation time.  In order to not trigger a type error later we
                # just remove them.  This might leave us with a Decorator with no
                # decorators (Emperor's new clothes?)
                # TODO: It would be nice to type-check these rather than remove them.
                #       default should be Callable[[], T]
                #       validator should be Callable[[Any, 'Attribute', T], Any]
                #       where T is the type of the attribute.
                remove_me.append(func_decorator)
    for dec in remove_me:
        stmt.decorators.remove(dec)


def _attribute_from_auto_attrib(ctx: 'mypy.plugin.ClassDefContext',
                                kw_only: bool,
                                lhs: NameExpr,
                                rvalue: Expression,
                                stmt: AssignmentStmt) -> Attribute:
    """Return an Attribute for a new type assignment."""
    name = unmangle(lhs.name)
    # `x: int` (without equal sign) assigns rvalue to TempNode(AnyType())
    has_rhs = not isinstance(rvalue, TempNode)
    sym = ctx.cls.info.names.get(name)
    init_type = sym.type if sym else None
    return Attribute(name, ctx.cls.info, has_rhs, True, kw_only, Converter(), stmt, init_type)


def _attribute_from_attrib_maker(ctx: 'mypy.plugin.ClassDefContext',
                                 auto_attribs: bool,
                                 kw_only: bool,
                                 lhs: NameExpr,
                                 rvalue: CallExpr,
                                 stmt: AssignmentStmt) -> Optional[Attribute]:
    """Return an Attribute from the assignment or None if you can't make one."""
    if auto_attribs and not stmt.new_syntax:
        # auto_attribs requires an annotation on *every* attr.ib.
        assert lhs.node is not None
        ctx.api.msg.need_annotation_for_var(lhs.node, stmt)
        return None

    if len(stmt.lvalues) > 1:
        ctx.api.fail("Too many names for one attribute", stmt)
        return None

    # This is the type that belongs in the __init__ method for this attrib.
    init_type = stmt.type

    # Read all the arguments from the call.
    init = _get_bool_argument(ctx, rvalue, 'init', True)
    # Note: If the class decorator says kw_only=True the attribute is ignored.
    # See https://github.com/python-attrs/attrs/issues/481 for explanation.
    kw_only |= _get_bool_argument(ctx, rvalue, 'kw_only', False)
    if kw_only and ctx.api.options.python_version[0] < 3:
        ctx.api.fail(KW_ONLY_PYTHON_2_UNSUPPORTED, stmt)
        return None

    # TODO: Check for attr.NOTHING
    attr_has_default = bool(_get_argument(rvalue, 'default'))
    attr_has_factory = bool(_get_argument(rvalue, 'factory'))

    if attr_has_default and attr_has_factory:
        ctx.api.fail('Can\'t pass both "default" and "factory".', rvalue)
    elif attr_has_factory:
        attr_has_default = True

    # If the type isn't set through annotation but is passed through `type=` use that.
    type_arg = _get_argument(rvalue, 'type')
    if type_arg and not init_type:
        try:
            un_type = expr_to_unanalyzed_type(type_arg, ctx.api.options, ctx.api.is_stub_file)
        except TypeTranslationError:
            ctx.api.fail('Invalid argument to type', type_arg)
        else:
            init_type = ctx.api.anal_type(un_type)
            if init_type and isinstance(lhs.node, Var) and not lhs.node.type:
                # If there is no annotation, add one.
                lhs.node.type = init_type
                lhs.is_inferred_def = False

    # Note: convert is deprecated but works the same as converter.
    converter = _get_argument(rvalue, 'converter')
    convert = _get_argument(rvalue, 'convert')
    if convert and converter:
        ctx.api.fail('Can\'t pass both "convert" and "converter".', rvalue)
    elif convert:
        ctx.api.fail("convert is deprecated, use converter", rvalue)
        converter = convert
    converter_info = _parse_converter(ctx, converter)

    name = unmangle(lhs.name)
    return Attribute(name, ctx.cls.info, attr_has_default, init,
                     kw_only, converter_info, stmt, init_type)


def _parse_converter(ctx: 'mypy.plugin.ClassDefContext',
                     converter: Optional[Expression]) -> Converter:
    """Return the Converter object from an Expression."""
    # TODO: Support complex converters, e.g. lambdas, calls, etc.
    if converter:
        if isinstance(converter, RefExpr) and converter.node:
            if (isinstance(converter.node, FuncDef)
                    and converter.node.type
                    and isinstance(converter.node.type, FunctionLike)):
                return Converter(converter.node.fullname)
            elif (isinstance(converter.node, OverloadedFuncDef)
                    and is_valid_overloaded_converter(converter.node)):
                return Converter(converter.node.fullname)
            elif isinstance(converter.node, TypeInfo):
                return Converter(converter.node.fullname)

        if (isinstance(converter, CallExpr)
                and isinstance(converter.callee, RefExpr)
                and converter.callee.fullname == "attr.converters.optional"
                and converter.args
                and converter.args[0]):
            # Special handling for attr.converters.optional(type)
            # We extract the type and add make the init_args Optional in Attribute.argument
            argument = _parse_converter(ctx, converter.args[0])
            argument.is_attr_converters_optional = True
            return argument

        # Signal that we have an unsupported converter.
        ctx.api.fail(
            "Unsupported converter, only named functions and types are currently supported",
            converter
        )
        return Converter('')
    return Converter(None)


def is_valid_overloaded_converter(defn: OverloadedFuncDef) -> bool:
    return all((not isinstance(item, Decorator) or isinstance(item.func.type, FunctionLike))
               for item in defn.items)


def _parse_assignments(
        lvalue: Expression,
        stmt: AssignmentStmt) -> Tuple[List[NameExpr], List[Expression]]:
    """Convert a possibly complex assignment expression into lists of lvalues and rvalues."""
    lvalues: List[NameExpr] = []
    rvalues: List[Expression] = []
    if isinstance(lvalue, (TupleExpr, ListExpr)):
        if all(isinstance(item, NameExpr) for item in lvalue.items):
            lvalues = cast(List[NameExpr], lvalue.items)
        if isinstance(stmt.rvalue, (TupleExpr, ListExpr)):
            rvalues = stmt.rvalue.items
    elif isinstance(lvalue, NameExpr):
        lvalues = [lvalue]
        rvalues = [stmt.rvalue]
    return lvalues, rvalues


def _add_order(ctx: 'mypy.plugin.ClassDefContext', adder: 'MethodAdder') -> None:
    """Generate all the ordering methods for this class."""
    bool_type = ctx.api.named_type('builtins.bool')
    object_type = ctx.api.named_type('builtins.object')
    # Make the types be:
    #    AT = TypeVar('AT')
    #    def __lt__(self: AT, other: AT) -> bool
    # This way comparisons with subclasses will work correctly.
    tvd = TypeVarType(SELF_TVAR_NAME, ctx.cls.info.fullname + '.' + SELF_TVAR_NAME,
                     -1, [], object_type)
    self_tvar_expr = TypeVarExpr(SELF_TVAR_NAME, ctx.cls.info.fullname + '.' + SELF_TVAR_NAME,
                                 [], object_type)
    ctx.cls.info.names[SELF_TVAR_NAME] = SymbolTableNode(MDEF, self_tvar_expr)

    args = [Argument(Var('other', tvd), tvd, None, ARG_POS)]
    for method in ['__lt__', '__le__', '__gt__', '__ge__']:
        adder.add_method(method, args, bool_type, self_type=tvd, tvd=tvd)


def _make_frozen(ctx: 'mypy.plugin.ClassDefContext', attributes: List[Attribute]) -> None:
    """Turn all the attributes into properties to simulate frozen classes."""
    for attribute in attributes:
        if attribute.name in ctx.cls.info.names:
            # This variable belongs to this class so we can modify it.
            node = ctx.cls.info.names[attribute.name].node
            assert isinstance(node, Var)
            node.is_property = True
        else:
            # This variable belongs to a super class so create new Var so we
            # can modify it.
            var = Var(attribute.name, ctx.cls.info[attribute.name].type)
            var.info = ctx.cls.info
            var._fullname = '%s.%s' % (ctx.cls.info.fullname, var.name)
            ctx.cls.info.names[var.name] = SymbolTableNode(MDEF, var)
            var.is_property = True


def _add_init(ctx: 'mypy.plugin.ClassDefContext', attributes: List[Attribute],
              adder: 'MethodAdder') -> None:
    """Generate an __init__ method for the attributes and add it to the class."""
    # Convert attributes to arguments with kw_only arguments at the  end of
    # the argument list
    pos_args = []
    kw_only_args = []
    for attribute in attributes:
        if not attribute.init:
            continue
        if attribute.kw_only:
            kw_only_args.append(attribute.argument(ctx))
        else:
            pos_args.append(attribute.argument(ctx))
    args = pos_args + kw_only_args
    if all(
        # We use getattr rather than instance checks because the variable.type
        # might be wrapped into a Union or some other type, but even non-Any
        # types reliably track the fact that the argument was not annotated.
        getattr(arg.variable.type, "type_of_any", None) == TypeOfAny.unannotated
        for arg in args
    ):
        # This workaround makes --disallow-incomplete-defs usable with attrs,
        # but is definitely suboptimal as a long-term solution.
        # See https://github.com/python/mypy/issues/5954 for discussion.
        for a in args:
            a.variable.type = AnyType(TypeOfAny.implementation_artifact)
            a.type_annotation = AnyType(TypeOfAny.implementation_artifact)
    adder.add_method('__init__', args, NoneType())


def _add_attrs_magic_attribute(ctx: 'mypy.plugin.ClassDefContext',
                               raw_attr_types: 'List[Optional[Type]]') -> None:
    attr_name = '__attrs_attrs__'
    any_type = AnyType(TypeOfAny.explicit)
    attributes_types: 'List[Type]' = [
        ctx.api.named_type_or_none('attr.Attribute', [attr_type or any_type]) or any_type
        for attr_type in raw_attr_types
    ]
    fallback_type = ctx.api.named_type('builtins.tuple', [
        ctx.api.named_type_or_none('attr.Attribute', [any_type]) or any_type,
    ])
    var = Var(name=attr_name, type=TupleType(attributes_types, fallback=fallback_type))
    var.info = ctx.cls.info
    var._fullname = ctx.cls.info.fullname + '.' + attr_name
    ctx.cls.info.names[attr_name] = SymbolTableNode(
        kind=MDEF,
        node=var,
        plugin_generated=True,
    )


def _add_slots(ctx: 'mypy.plugin.ClassDefContext',
               attributes: List[Attribute]) -> None:
    # Unlike `@dataclasses.dataclass`, `__slots__` is rewritten here.
    ctx.cls.info.slots = {attr.name for attr in attributes}


class MethodAdder:
    """Helper to add methods to a TypeInfo.

    ctx: The ClassDefCtx we are using on which we will add methods.
    """

    # TODO: Combine this with the code build_namedtuple_typeinfo to support both.

    def __init__(self, ctx: 'mypy.plugin.ClassDefContext') -> None:
        self.ctx = ctx
        self.self_type = fill_typevars(ctx.cls.info)

    def add_method(self,
                   method_name: str, args: List[Argument], ret_type: Type,
                   self_type: Optional[Type] = None,
                   tvd: Optional[TypeVarType] = None) -> None:
        """Add a method: def <method_name>(self, <args>) -> <ret_type>): ... to info.

        self_type: The type to use for the self argument or None to use the inferred self type.
        tvd: If the method is generic these should be the type variables.
        """
        self_type = self_type if self_type is not None else self.self_type
        add_method(self.ctx, method_name, args, ret_type, self_type, tvd)
