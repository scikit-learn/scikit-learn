from typing import List, Optional, Union

from mypy.nodes import (
    ARG_POS, MDEF, Argument, Block, CallExpr, ClassDef, Expression, SYMBOL_FUNCBASE_TYPES,
    FuncDef, PassStmt, RefExpr, SymbolTableNode, Var, JsonDict,
)
from mypy.plugin import CheckerPluginInterface, ClassDefContext, SemanticAnalyzerPluginInterface
from mypy.semanal import set_callable_name
from mypy.types import (
    CallableType, Overloaded, Type, TypeVarType, deserialize_type, get_proper_type,
)
from mypy.typevars import fill_typevars
from mypy.util import get_unique_redefinition_name
from mypy.typeops import try_getting_str_literals  # noqa: F401  # Part of public API
from mypy.fixup import TypeFixer


def _get_decorator_bool_argument(
        ctx: ClassDefContext,
        name: str,
        default: bool,
) -> bool:
    """Return the bool argument for the decorator.

    This handles both @decorator(...) and @decorator.
    """
    if isinstance(ctx.reason, CallExpr):
        return _get_bool_argument(ctx, ctx.reason, name, default)
    else:
        return default


def _get_bool_argument(ctx: ClassDefContext, expr: CallExpr,
                       name: str, default: bool) -> bool:
    """Return the boolean value for an argument to a call or the
    default if it's not found.
    """
    attr_value = _get_argument(expr, name)
    if attr_value:
        ret = ctx.api.parse_bool(attr_value)
        if ret is None:
            ctx.api.fail('"{}" argument must be True or False.'.format(name), expr)
            return default
        return ret
    return default


def _get_argument(call: CallExpr, name: str) -> Optional[Expression]:
    """Return the expression for the specific argument."""
    # To do this we use the CallableType of the callee to find the FormalArgument,
    # then walk the actual CallExpr looking for the appropriate argument.
    #
    # Note: I'm not hard-coding the index so that in the future we can support other
    # attrib and class makers.
    if not isinstance(call.callee, RefExpr):
        return None

    callee_type = None
    callee_node = call.callee.node
    if (isinstance(callee_node, (Var, SYMBOL_FUNCBASE_TYPES))
            and callee_node.type):
        callee_node_type = get_proper_type(callee_node.type)
        if isinstance(callee_node_type, Overloaded):
            # We take the last overload.
            callee_type = callee_node_type.items[-1]
        elif isinstance(callee_node_type, CallableType):
            callee_type = callee_node_type

    if not callee_type:
        return None

    argument = callee_type.argument_by_name(name)
    if not argument:
        return None
    assert argument.name

    for i, (attr_name, attr_value) in enumerate(zip(call.arg_names, call.args)):
        if argument.pos is not None and not attr_name and i == argument.pos:
            return attr_value
        if attr_name == argument.name:
            return attr_value
    return None


def add_method(
        ctx: ClassDefContext,
        name: str,
        args: List[Argument],
        return_type: Type,
        self_type: Optional[Type] = None,
        tvar_def: Optional[TypeVarType] = None,
) -> None:
    """
    Adds a new method to a class.
    Deprecated, use add_method_to_class() instead.
    """
    add_method_to_class(ctx.api, ctx.cls,
                        name=name,
                        args=args,
                        return_type=return_type,
                        self_type=self_type,
                        tvar_def=tvar_def)


def add_method_to_class(
        api: Union[SemanticAnalyzerPluginInterface, CheckerPluginInterface],
        cls: ClassDef,
        name: str,
        args: List[Argument],
        return_type: Type,
        self_type: Optional[Type] = None,
        tvar_def: Optional[TypeVarType] = None,
) -> None:
    """Adds a new method to a class definition."""
    info = cls.info

    # First remove any previously generated methods with the same name
    # to avoid clashes and problems in the semantic analyzer.
    if name in info.names:
        sym = info.names[name]
        if sym.plugin_generated and isinstance(sym.node, FuncDef):
            cls.defs.body.remove(sym.node)

    self_type = self_type or fill_typevars(info)
    if isinstance(api, SemanticAnalyzerPluginInterface):
        function_type = api.named_type('builtins.function')
    else:
        function_type = api.named_generic_type('builtins.function', [])

    args = [Argument(Var('self'), self_type, None, ARG_POS)] + args
    arg_types, arg_names, arg_kinds = [], [], []
    for arg in args:
        assert arg.type_annotation, 'All arguments must be fully typed.'
        arg_types.append(arg.type_annotation)
        arg_names.append(arg.variable.name)
        arg_kinds.append(arg.kind)

    signature = CallableType(arg_types, arg_kinds, arg_names, return_type, function_type)
    if tvar_def:
        signature.variables = [tvar_def]

    func = FuncDef(name, args, Block([PassStmt()]))
    func.info = info
    func.type = set_callable_name(signature, func)
    func._fullname = info.fullname + '.' + name
    func.line = info.line

    # NOTE: we would like the plugin generated node to dominate, but we still
    # need to keep any existing definitions so they get semantically analyzed.
    if name in info.names:
        # Get a nice unique name instead.
        r_name = get_unique_redefinition_name(name, info.names)
        info.names[r_name] = info.names[name]

    info.names[name] = SymbolTableNode(MDEF, func, plugin_generated=True)
    info.defn.defs.body.append(func)


def deserialize_and_fixup_type(
    data: Union[str, JsonDict], api: SemanticAnalyzerPluginInterface
) -> Type:
    typ = deserialize_type(data)
    typ.accept(TypeFixer(api.modules, allow_missing=False))
    return typ
