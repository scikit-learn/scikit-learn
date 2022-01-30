"""Semantic analysis of call-based Enum definitions.

This is conceptually part of mypy.semanal (semantic analyzer pass 2).
"""

from typing import List, Tuple, Optional, Union, cast
from typing_extensions import Final

from mypy.nodes import (
    Expression, Context, TypeInfo, AssignmentStmt, NameExpr, CallExpr, RefExpr, StrExpr,
    UnicodeExpr, TupleExpr, ListExpr, DictExpr, Var, SymbolTableNode, MDEF, ARG_POS,
    ARG_NAMED, EnumCallExpr, MemberExpr
)
from mypy.semanal_shared import SemanticAnalyzerInterface
from mypy.options import Options

# Note: 'enum.EnumMeta' is deliberately excluded from this list. Classes that directly use
# enum.EnumMeta do not necessarily automatically have the 'name' and 'value' attributes.
ENUM_BASES: Final = frozenset((
    'enum.Enum', 'enum.IntEnum', 'enum.Flag', 'enum.IntFlag',
))
ENUM_SPECIAL_PROPS: Final = frozenset((
    'name', 'value', '_name_', '_value_', '_order_', '__order__',
    # Also attributes from `object`:
    '__module__', '__annotations__', '__doc__', '__slots__', '__dict__',
))


class EnumCallAnalyzer:
    def __init__(self, options: Options, api: SemanticAnalyzerInterface) -> None:
        self.options = options
        self.api = api

    def process_enum_call(self, s: AssignmentStmt, is_func_scope: bool) -> bool:
        """Check if s defines an Enum; if yes, store the definition in symbol table.

        Return True if this looks like an Enum definition (but maybe with errors),
        otherwise return False.
        """
        if len(s.lvalues) != 1 or not isinstance(s.lvalues[0], (NameExpr, MemberExpr)):
            return False
        lvalue = s.lvalues[0]
        name = lvalue.name
        enum_call = self.check_enum_call(s.rvalue, name, is_func_scope)
        if enum_call is None:
            return False
        if isinstance(lvalue, MemberExpr):
            self.fail("Enum type as attribute is not supported", lvalue)
            return False
        # Yes, it's a valid Enum definition. Add it to the symbol table.
        self.api.add_symbol(name, enum_call, s)
        return True

    def check_enum_call(self,
                        node: Expression,
                        var_name: str,
                        is_func_scope: bool) -> Optional[TypeInfo]:
        """Check if a call defines an Enum.

        Example:

          A = enum.Enum('A', 'foo bar')

        is equivalent to:

          class A(enum.Enum):
              foo = 1
              bar = 2
        """
        if not isinstance(node, CallExpr):
            return None
        call = node
        callee = call.callee
        if not isinstance(callee, RefExpr):
            return None
        fullname = callee.fullname
        if fullname not in ENUM_BASES:
            return None
        items, values, ok = self.parse_enum_call_args(call, fullname.split('.')[-1])
        if not ok:
            # Error. Construct dummy return value.
            info = self.build_enum_call_typeinfo(var_name, [], fullname, node.line)
        else:
            name = cast(Union[StrExpr, UnicodeExpr], call.args[0]).value
            if name != var_name or is_func_scope:
                # Give it a unique name derived from the line number.
                name += '@' + str(call.line)
            info = self.build_enum_call_typeinfo(name, items, fullname, call.line)
            # Store generated TypeInfo under both names, see semanal_namedtuple for more details.
            if name != var_name or is_func_scope:
                self.api.add_symbol_skip_local(name, info)
        call.analyzed = EnumCallExpr(info, items, values)
        call.analyzed.set_line(call.line, call.column)
        info.line = node.line
        return info

    def build_enum_call_typeinfo(self, name: str, items: List[str], fullname: str,
                                 line: int) -> TypeInfo:
        base = self.api.named_type_or_none(fullname)
        assert base is not None
        info = self.api.basic_new_typeinfo(name, base, line)
        info.metaclass_type = info.calculate_metaclass_type()
        info.is_enum = True
        for item in items:
            var = Var(item)
            var.info = info
            var.is_property = True
            var._fullname = '{}.{}'.format(info.fullname, item)
            info.names[item] = SymbolTableNode(MDEF, var)
        return info

    def parse_enum_call_args(self, call: CallExpr,
                             class_name: str) -> Tuple[List[str],
                                                       List[Optional[Expression]], bool]:
        """Parse arguments of an Enum call.

        Return a tuple of fields, values, was there an error.
        """
        args = call.args
        if not all([arg_kind in [ARG_POS, ARG_NAMED] for arg_kind in call.arg_kinds]):
            return self.fail_enum_call_arg("Unexpected arguments to %s()" % class_name, call)
        if len(args) < 2:
            return self.fail_enum_call_arg("Too few arguments for %s()" % class_name, call)
        if len(args) > 6:
            return self.fail_enum_call_arg("Too many arguments for %s()" % class_name, call)
        valid_name = [None, 'value', 'names', 'module', 'qualname', 'type', 'start']
        for arg_name in call.arg_names:
            if arg_name not in valid_name:
                self.fail_enum_call_arg('Unexpected keyword argument "{}"'.format(arg_name), call)
        value, names = None, None
        for arg_name, arg in zip(call.arg_names, args):
            if arg_name == 'value':
                value = arg
            if arg_name == 'names':
                names = arg
        if value is None:
            value = args[0]
        if names is None:
            names = args[1]
        if not isinstance(value, (StrExpr, UnicodeExpr)):
            return self.fail_enum_call_arg(
                "%s() expects a string literal as the first argument" % class_name, call)
        items = []
        values: List[Optional[Expression]] = []
        if isinstance(names, (StrExpr, UnicodeExpr)):
            fields = names.value
            for field in fields.replace(',', ' ').split():
                items.append(field)
        elif isinstance(names, (TupleExpr, ListExpr)):
            seq_items = names.items
            if all(isinstance(seq_item, (StrExpr, UnicodeExpr)) for seq_item in seq_items):
                items = [cast(Union[StrExpr, UnicodeExpr], seq_item).value
                         for seq_item in seq_items]
            elif all(isinstance(seq_item, (TupleExpr, ListExpr))
                     and len(seq_item.items) == 2
                     and isinstance(seq_item.items[0], (StrExpr, UnicodeExpr))
                     for seq_item in seq_items):
                for seq_item in seq_items:
                    assert isinstance(seq_item, (TupleExpr, ListExpr))
                    name, value = seq_item.items
                    assert isinstance(name, (StrExpr, UnicodeExpr))
                    items.append(name.value)
                    values.append(value)
            else:
                return self.fail_enum_call_arg(
                    "%s() with tuple or list expects strings or (name, value) pairs" %
                    class_name,
                    call)
        elif isinstance(names, DictExpr):
            for key, value in names.items:
                if not isinstance(key, (StrExpr, UnicodeExpr)):
                    return self.fail_enum_call_arg(
                        "%s() with dict literal requires string literals" % class_name, call)
                items.append(key.value)
                values.append(value)
        else:
            # TODO: Allow dict(x=1, y=2) as a substitute for {'x': 1, 'y': 2}?
            return self.fail_enum_call_arg(
                "%s() expects a string, tuple, list or dict literal as the second argument" %
                class_name,
                call)
        if len(items) == 0:
            return self.fail_enum_call_arg("%s() needs at least one item" % class_name, call)
        if not values:
            values = [None] * len(items)
        assert len(items) == len(values)
        return items, values, True

    def fail_enum_call_arg(self, message: str,
                           context: Context) -> Tuple[List[str],
                                                      List[Optional[Expression]], bool]:
        self.fail(message, context)
        return [], [], False

    # Helpers

    def fail(self, msg: str, ctx: Context) -> None:
        self.api.fail(msg, ctx)
