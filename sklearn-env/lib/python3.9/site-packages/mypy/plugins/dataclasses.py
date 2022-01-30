"""Plugin that provides support for dataclasses."""

from typing import Dict, List, Set, Tuple, Optional
from typing_extensions import Final

from mypy.nodes import (
    ARG_OPT, ARG_NAMED, ARG_NAMED_OPT, ARG_POS, MDEF, Argument, AssignmentStmt, CallExpr,
    Context, Expression, JsonDict, NameExpr, RefExpr,
    SymbolTableNode, TempNode, TypeInfo, Var, TypeVarExpr, PlaceholderNode
)
from mypy.plugin import ClassDefContext, SemanticAnalyzerPluginInterface
from mypy.plugins.common import (
    add_method, _get_decorator_bool_argument, deserialize_and_fixup_type,
)
from mypy.typeops import map_type_from_supertype
from mypy.types import (
    Type, Instance, NoneType, TypeVarType, CallableType, get_proper_type,
    AnyType, TypeOfAny,
)
from mypy.server.trigger import make_wildcard_trigger

# The set of decorators that generate dataclasses.
dataclass_makers: Final = {
    'dataclass',
    'dataclasses.dataclass',
}
# The set of functions that generate dataclass fields.
field_makers: Final = {
    'dataclasses.field',
}


SELF_TVAR_NAME: Final = "_DT"


class DataclassAttribute:
    def __init__(
            self,
            name: str,
            is_in_init: bool,
            is_init_var: bool,
            has_default: bool,
            line: int,
            column: int,
            type: Optional[Type],
            info: TypeInfo,
            kw_only: bool,
    ) -> None:
        self.name = name
        self.is_in_init = is_in_init
        self.is_init_var = is_init_var
        self.has_default = has_default
        self.line = line
        self.column = column
        self.type = type
        self.info = info
        self.kw_only = kw_only

    def to_argument(self) -> Argument:
        arg_kind = ARG_POS
        if self.kw_only and self.has_default:
            arg_kind = ARG_NAMED_OPT
        elif self.kw_only and not self.has_default:
            arg_kind = ARG_NAMED
        elif not self.kw_only and self.has_default:
            arg_kind = ARG_OPT
        return Argument(
            variable=self.to_var(),
            type_annotation=self.type,
            initializer=None,
            kind=arg_kind,
        )

    def to_var(self) -> Var:
        return Var(self.name, self.type)

    def serialize(self) -> JsonDict:
        assert self.type
        return {
            'name': self.name,
            'is_in_init': self.is_in_init,
            'is_init_var': self.is_init_var,
            'has_default': self.has_default,
            'line': self.line,
            'column': self.column,
            'type': self.type.serialize(),
            'kw_only': self.kw_only,
        }

    @classmethod
    def deserialize(
        cls, info: TypeInfo, data: JsonDict, api: SemanticAnalyzerPluginInterface
    ) -> 'DataclassAttribute':
        data = data.copy()
        if data.get('kw_only') is None:
            data['kw_only'] = False
        typ = deserialize_and_fixup_type(data.pop('type'), api)
        return cls(type=typ, info=info, **data)

    def expand_typevar_from_subtype(self, sub_type: TypeInfo) -> None:
        """Expands type vars in the context of a subtype when an attribute is inherited
        from a generic super type."""
        if not isinstance(self.type, TypeVarType):
            return

        self.type = map_type_from_supertype(self.type, sub_type, self.info)


class DataclassTransformer:
    def __init__(self, ctx: ClassDefContext) -> None:
        self._ctx = ctx

    def transform(self) -> None:
        """Apply all the necessary transformations to the underlying
        dataclass so as to ensure it is fully type checked according
        to the rules in PEP 557.
        """
        ctx = self._ctx
        info = self._ctx.cls.info
        attributes = self.collect_attributes()
        if attributes is None:
            # Some definitions are not ready, defer() should be already called.
            return
        for attr in attributes:
            if attr.type is None:
                ctx.api.defer()
                return
        decorator_arguments = {
            'init': _get_decorator_bool_argument(self._ctx, 'init', True),
            'eq': _get_decorator_bool_argument(self._ctx, 'eq', True),
            'order': _get_decorator_bool_argument(self._ctx, 'order', False),
            'frozen': _get_decorator_bool_argument(self._ctx, 'frozen', False),
            'slots': _get_decorator_bool_argument(self._ctx, 'slots', False),
        }
        py_version = self._ctx.api.options.python_version

        # If there are no attributes, it may be that the semantic analyzer has not
        # processed them yet. In order to work around this, we can simply skip generating
        # __init__ if there are no attributes, because if the user truly did not define any,
        # then the object default __init__ with an empty signature will be present anyway.
        if (decorator_arguments['init'] and
                ('__init__' not in info.names or info.names['__init__'].plugin_generated) and
                attributes):
            add_method(
                ctx,
                '__init__',
                args=[attr.to_argument() for attr in attributes if attr.is_in_init
                      and not self._is_kw_only_type(attr.type)],
                return_type=NoneType(),
            )

        if (decorator_arguments['eq'] and info.get('__eq__') is None or
                decorator_arguments['order']):
            # Type variable for self types in generated methods.
            obj_type = ctx.api.named_type('builtins.object')
            self_tvar_expr = TypeVarExpr(SELF_TVAR_NAME, info.fullname + '.' + SELF_TVAR_NAME,
                                         [], obj_type)
            info.names[SELF_TVAR_NAME] = SymbolTableNode(MDEF, self_tvar_expr)

        # Add <, >, <=, >=, but only if the class has an eq method.
        if decorator_arguments['order']:
            if not decorator_arguments['eq']:
                ctx.api.fail('eq must be True if order is True', ctx.cls)

            for method_name in ['__lt__', '__gt__', '__le__', '__ge__']:
                # Like for __eq__ and __ne__, we want "other" to match
                # the self type.
                obj_type = ctx.api.named_type('builtins.object')
                order_tvar_def = TypeVarType(SELF_TVAR_NAME, info.fullname + '.' + SELF_TVAR_NAME,
                                            -1, [], obj_type)
                order_return_type = ctx.api.named_type('builtins.bool')
                order_args = [
                    Argument(Var('other', order_tvar_def), order_tvar_def, None, ARG_POS)
                ]

                existing_method = info.get(method_name)
                if existing_method is not None and not existing_method.plugin_generated:
                    assert existing_method.node
                    ctx.api.fail(
                        'You may not have a custom %s method when order=True' % method_name,
                        existing_method.node,
                    )

                add_method(
                    ctx,
                    method_name,
                    args=order_args,
                    return_type=order_return_type,
                    self_type=order_tvar_def,
                    tvar_def=order_tvar_def,
                )

        if decorator_arguments['frozen']:
            self._freeze(attributes)
        else:
            self._propertize_callables(attributes)

        if decorator_arguments['slots']:
            self.add_slots(info, attributes, correct_version=py_version >= (3, 10))

        self.reset_init_only_vars(info, attributes)

        self._add_dataclass_fields_magic_attribute()

        info.metadata['dataclass'] = {
            'attributes': [attr.serialize() for attr in attributes],
            'frozen': decorator_arguments['frozen'],
        }

    def add_slots(self,
                  info: TypeInfo,
                  attributes: List[DataclassAttribute],
                  *,
                  correct_version: bool) -> None:
        if not correct_version:
            # This means that version is lower than `3.10`,
            # it is just a non-existent argument for `dataclass` function.
            self._ctx.api.fail(
                'Keyword argument "slots" for "dataclass" '
                'is only valid in Python 3.10 and higher',
                self._ctx.reason,
            )
            return

        generated_slots = {attr.name for attr in attributes}
        if ((info.slots is not None and info.slots != generated_slots)
                or info.names.get('__slots__')):
            # This means we have a slots conflict.
            # Class explicitly specifies a different `__slots__` field.
            # And `@dataclass(slots=True)` is used.
            # In runtime this raises a type error.
            self._ctx.api.fail(
                '"{}" both defines "__slots__" and is used with "slots=True"'.format(
                    self._ctx.cls.name,
                ),
                self._ctx.cls,
            )
            return

        info.slots = generated_slots

    def reset_init_only_vars(self, info: TypeInfo, attributes: List[DataclassAttribute]) -> None:
        """Remove init-only vars from the class and reset init var declarations."""
        for attr in attributes:
            if attr.is_init_var:
                if attr.name in info.names:
                    del info.names[attr.name]
                else:
                    # Nodes of superclass InitVars not used in __init__ cannot be reached.
                    assert attr.is_init_var
                for stmt in info.defn.defs.body:
                    if isinstance(stmt, AssignmentStmt) and stmt.unanalyzed_type:
                        lvalue = stmt.lvalues[0]
                        if isinstance(lvalue, NameExpr) and lvalue.name == attr.name:
                            # Reset node so that another semantic analysis pass will
                            # recreate a symbol node for this attribute.
                            lvalue.node = None

    def collect_attributes(self) -> Optional[List[DataclassAttribute]]:
        """Collect all attributes declared in the dataclass and its parents.

        All assignments of the form

          a: SomeType
          b: SomeOtherType = ...

        are collected.
        """
        # First, collect attributes belonging to the current class.
        ctx = self._ctx
        cls = self._ctx.cls
        attrs: List[DataclassAttribute] = []
        known_attrs: Set[str] = set()
        kw_only = _get_decorator_bool_argument(ctx, 'kw_only', False)
        for stmt in cls.defs.body:
            # Any assignment that doesn't use the new type declaration
            # syntax can be ignored out of hand.
            if not (isinstance(stmt, AssignmentStmt) and stmt.new_syntax):
                continue

            # a: int, b: str = 1, 'foo' is not supported syntax so we
            # don't have to worry about it.
            lhs = stmt.lvalues[0]
            if not isinstance(lhs, NameExpr):
                continue

            sym = cls.info.names.get(lhs.name)
            if sym is None:
                # This name is likely blocked by a star import. We don't need to defer because
                # defer() is already called by mark_incomplete().
                continue

            node = sym.node
            if isinstance(node, PlaceholderNode):
                # This node is not ready yet.
                return None
            assert isinstance(node, Var)

            # x: ClassVar[int] is ignored by dataclasses.
            if node.is_classvar:
                continue

            # x: InitVar[int] is turned into x: int and is removed from the class.
            is_init_var = False
            node_type = get_proper_type(node.type)
            if (isinstance(node_type, Instance) and
                    node_type.type.fullname == 'dataclasses.InitVar'):
                is_init_var = True
                node.type = node_type.args[0]

            if self._is_kw_only_type(node_type):
                kw_only = True

            has_field_call, field_args = _collect_field_args(stmt.rvalue, ctx)

            is_in_init_param = field_args.get('init')
            if is_in_init_param is None:
                is_in_init = True
            else:
                is_in_init = bool(ctx.api.parse_bool(is_in_init_param))

            has_default = False
            # Ensure that something like x: int = field() is rejected
            # after an attribute with a default.
            if has_field_call:
                has_default = 'default' in field_args or 'default_factory' in field_args

            # All other assignments are already type checked.
            elif not isinstance(stmt.rvalue, TempNode):
                has_default = True

            if not has_default:
                # Make all non-default attributes implicit because they are de-facto set
                # on self in the generated __init__(), not in the class body.
                sym.implicit = True

            is_kw_only = kw_only
            # Use the kw_only field arg if it is provided. Otherwise use the
            # kw_only value from the decorator parameter.
            field_kw_only_param = field_args.get('kw_only')
            if field_kw_only_param is not None:
                is_kw_only = bool(ctx.api.parse_bool(field_kw_only_param))

            known_attrs.add(lhs.name)
            attrs.append(DataclassAttribute(
                name=lhs.name,
                is_in_init=is_in_init,
                is_init_var=is_init_var,
                has_default=has_default,
                line=stmt.line,
                column=stmt.column,
                type=sym.type,
                info=cls.info,
                kw_only=is_kw_only,
            ))

        # Next, collect attributes belonging to any class in the MRO
        # as long as those attributes weren't already collected.  This
        # makes it possible to overwrite attributes in subclasses.
        # copy() because we potentially modify all_attrs below and if this code requires debugging
        # we'll have unmodified attrs laying around.
        all_attrs = attrs.copy()
        for info in cls.info.mro[1:-1]:
            if 'dataclass' not in info.metadata:
                continue

            super_attrs = []
            # Each class depends on the set of attributes in its dataclass ancestors.
            ctx.api.add_plugin_dependency(make_wildcard_trigger(info.fullname))

            for data in info.metadata["dataclass"]["attributes"]:
                name: str = data["name"]
                if name not in known_attrs:
                    attr = DataclassAttribute.deserialize(info, data, ctx.api)
                    attr.expand_typevar_from_subtype(ctx.cls.info)
                    known_attrs.add(name)
                    super_attrs.append(attr)
                elif all_attrs:
                    # How early in the attribute list an attribute appears is determined by the
                    # reverse MRO, not simply MRO.
                    # See https://docs.python.org/3/library/dataclasses.html#inheritance for
                    # details.
                    for attr in all_attrs:
                        if attr.name == name:
                            all_attrs.remove(attr)
                            super_attrs.append(attr)
                            break
            all_attrs = super_attrs + all_attrs
            all_attrs.sort(key=lambda a: a.kw_only)

        # Ensure that arguments without a default don't follow
        # arguments that have a default.
        found_default = False
        # Ensure that the KW_ONLY sentinel is only provided once
        found_kw_sentinel = False
        for attr in all_attrs:
            # If we find any attribute that is_in_init, not kw_only, and that
            # doesn't have a default after one that does have one,
            # then that's an error.
            if found_default and attr.is_in_init and not attr.has_default and not attr.kw_only:
                # If the issue comes from merging different classes, report it
                # at the class definition point.
                context = (Context(line=attr.line, column=attr.column) if attr in attrs
                           else ctx.cls)
                ctx.api.fail(
                    'Attributes without a default cannot follow attributes with one',
                    context,
                )

            found_default = found_default or (attr.has_default and attr.is_in_init)
            if found_kw_sentinel and self._is_kw_only_type(attr.type):
                context = (Context(line=attr.line, column=attr.column) if attr in attrs
                           else ctx.cls)
                ctx.api.fail(
                    'There may not be more than one field with the KW_ONLY type',
                    context,
                )
            found_kw_sentinel = found_kw_sentinel or self._is_kw_only_type(attr.type)

        return all_attrs

    def _freeze(self, attributes: List[DataclassAttribute]) -> None:
        """Converts all attributes to @property methods in order to
        emulate frozen classes.
        """
        info = self._ctx.cls.info
        for attr in attributes:
            sym_node = info.names.get(attr.name)
            if sym_node is not None:
                var = sym_node.node
                assert isinstance(var, Var)
                var.is_property = True
            else:
                var = attr.to_var()
                var.info = info
                var.is_property = True
                var._fullname = info.fullname + '.' + var.name
                info.names[var.name] = SymbolTableNode(MDEF, var)

    def _propertize_callables(self, attributes: List[DataclassAttribute]) -> None:
        """Converts all attributes with callable types to @property methods.

        This avoids the typechecker getting confused and thinking that
        `my_dataclass_instance.callable_attr(foo)` is going to receive a
        `self` argument (it is not).

        """
        info = self._ctx.cls.info
        for attr in attributes:
            if isinstance(get_proper_type(attr.type), CallableType):
                var = attr.to_var()
                var.info = info
                var.is_property = True
                var.is_settable_property = True
                var._fullname = info.fullname + '.' + var.name
                info.names[var.name] = SymbolTableNode(MDEF, var)

    def _is_kw_only_type(self, node: Optional[Type]) -> bool:
        """Checks if the type of the node is the KW_ONLY sentinel value."""
        if node is None:
            return False
        node_type = get_proper_type(node)
        if not isinstance(node_type, Instance):
            return False
        return node_type.type.fullname == 'dataclasses.KW_ONLY'

    def _add_dataclass_fields_magic_attribute(self) -> None:
        attr_name = '__dataclass_fields__'
        any_type = AnyType(TypeOfAny.explicit)
        field_type = self._ctx.api.named_type_or_none('dataclasses.Field', [any_type]) or any_type
        attr_type = self._ctx.api.named_type('builtins.dict', [
            self._ctx.api.named_type('builtins.str'),
            field_type,
        ])
        var = Var(name=attr_name, type=attr_type)
        var.info = self._ctx.cls.info
        var._fullname = self._ctx.cls.info.fullname + '.' + attr_name
        self._ctx.cls.info.names[attr_name] = SymbolTableNode(
            kind=MDEF,
            node=var,
            plugin_generated=True,
        )


def dataclass_class_maker_callback(ctx: ClassDefContext) -> None:
    """Hooks into the class typechecking process to add support for dataclasses.
    """
    transformer = DataclassTransformer(ctx)
    transformer.transform()


def _collect_field_args(expr: Expression,
                        ctx: ClassDefContext) -> Tuple[bool, Dict[str, Expression]]:
    """Returns a tuple where the first value represents whether or not
    the expression is a call to dataclass.field and the second is a
    dictionary of the keyword arguments that field() was called with.
    """
    if (
            isinstance(expr, CallExpr) and
            isinstance(expr.callee, RefExpr) and
            expr.callee.fullname in field_makers
    ):
        # field() only takes keyword arguments.
        args = {}
        for name, arg, kind in zip(expr.arg_names, expr.args, expr.arg_kinds):
            if not kind.is_named():
                if kind.is_named(star=True):
                    # This means that `field` is used with `**` unpacking,
                    # the best we can do for now is not to fail.
                    # TODO: we can infer what's inside `**` and try to collect it.
                    message = 'Unpacking **kwargs in "field()" is not supported'
                else:
                    message = '"field()" does not accept positional arguments'
                ctx.api.fail(message, expr)
                return True, {}
            assert name is not None
            args[name] = arg
        return True, args
    return False, {}
