"""
ExtendedSyntaxCheck performs various syntax checks on the pythran AST.
"""

from pythran.passmanager import ModuleAnalysis
from pythran.analyses import StrictAliases, ArgumentEffects
from pythran.syntax import PythranSyntaxError
from pythran.intrinsic import ConstantIntr, FunctionIntr
from pythran import metadata

import gast as ast


def is_global_constant(node):
    if isinstance(node, ConstantIntr):
        return True

    if not isinstance(node, ast.FunctionDef):
        return False

    return metadata.get(node.body[0], metadata.StaticReturn)


def is_global(node):
    return (isinstance(node, (FunctionIntr, ast.FunctionDef)) or
            is_global_constant(node))


class ExtendedSyntaxCheck(ModuleAnalysis[StrictAliases, ArgumentEffects]):
    """
    Perform advanced syntax checking, based on strict aliases analysis:
    - is there a function redefinition?
    - is there a function call that does not match the called expression arity?
    - is there an operation that updates a global variable?
    """

    ResultType = type(None)
    def __init__(self):
        super().__init__()
        self.inassert = False
        self.functions = set()

    def check_global_with_side_effect(self, node, arg):
        if not isinstance(arg, ast.Call):
            return
        try:
            aliases = self.strict_aliases[arg.func]
        except KeyError:
            return

        for alias in aliases:
            if is_global_constant(alias):
                raise PythranSyntaxError(
                    ("Cannot modify '{}': global variables are constant "
                     "in pythran.").format(alias.name),
                    arg.func)

    def visit_FunctionDef(self, node):
        if node.name in self.functions:
            raise PythranSyntaxError("Function {} redefined".format(
                node.name),
                node)
        else:
            self.functions.add(node.name)
        self.generic_visit(node)

    def check_assert_with_side_effect(self, node, arg):
        if self.inassert:
            raise PythranSyntaxError("Cannot call a function with side effect "
                                     "in an assert", node)

    def visit_Assert(self, node):
        self.inassert = True
        self.generic_visit(node)
        self.inassert = False

    def is_immutable_constant(self, node):
        if isinstance(node, ast.Constant):
            return True

        if isinstance(node, ast.Tuple):
            return all(self.is_immutable_constant(elt) for elt in node.elts)

        if isinstance(node, ast.UnaryOp):
            return self.is_immutable_constant(node.operand)

        if isinstance(node, ast.Call):
            target = getattr(node, 'func', node)
            try:
                aliases = self.strict_aliases[target]
            except KeyError:
                return False

            if not aliases:
                return False

            if all(is_global_constant(alias) for alias in aliases):
                return True

        if isinstance(node, ast.Attribute):
            target = getattr(node, 'func', node)
            try:
                aliases = self.strict_aliases[target]
            except KeyError:
                return False

            if not aliases:
                return False

            if all(is_global(alias) for alias in aliases):
                return True

        if isinstance(node, ast.Name):
            try:
                aliases = self.strict_aliases[node]
            except KeyError:
                return False

            if all(isinstance(alias, ast.FunctionDef) for alias in aliases):
                return True

        return False

    def visit_arguments(self, node):
        self.generic_visit(node)
        for arg_default in node.defaults:
            if not self.is_immutable_constant(arg_default):
                raise PythranSyntaxError(
                    "Pythran does not support mutable default values. Use a "
                    "`None' default and set the value at runtime instead.",
                    arg_default)

    def visit_Call(self, node):
        self.generic_visit(node)
        func = node.func
        try:
            aliases = self.strict_aliases[func]
        except KeyError:
            raise PythranSyntaxError(
                "Call to unknown function `{}`, it's a trap!"
                .format(getattr(func, 'id', None) or func),
                node)

        argument_effects = set()

        for alias in aliases:
            # look for effect on arguments to prepare check on globals
            try:
                func_aes = self.argument_effects[alias]
                for i, effect in enumerate(func_aes):
                    if effect:
                        argument_effects.add(i)
            except KeyError:
                pass

            if not isinstance(alias, ast.FunctionDef):
                continue
            ubound = len(alias.args.args)
            lbound = ubound - len(alias.args.defaults)
            call_args_count = len(node.args) + len(node.keywords)
            if lbound <= call_args_count <= ubound:
                continue

            if lbound == ubound:
                msg = 'Invalid call to {}: expected {} arguments, got {}'
                msg = msg.format(alias.name,
                                 len(alias.args.args),
                                 len(node.args)
                                 )
            else:
                msg = ('Invalid {} call: '
                       'expected between {} and {} arguments, got {}')
                msg = msg.format(alias.name,
                                 lbound, ubound,
                                 len(node.args)
                                 )
            raise PythranSyntaxError(msg, node)

        # check for effects on globals
        for i, arg in enumerate(node.args):
            if i not in argument_effects:
                continue
            self.check_global_with_side_effect(node, arg)
            self.check_assert_with_side_effect(node, arg)
