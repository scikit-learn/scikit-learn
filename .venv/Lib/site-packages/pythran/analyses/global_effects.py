""" GlobalEffects computes function effect on global state. """

from pythran.analyses.aliases import Aliases
from pythran.analyses.intrinsics import Intrinsics
from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.passmanager import ModuleAnalysis
from pythran.tables import MODULES
from pythran.graph import DiGraph
import pythran.intrinsic as intrinsic

import gast as ast
from functools import reduce

class FunctionEffect(object):
    def __init__(self, node):
        self.func = node
        if isinstance(node, ast.FunctionDef):
            self.global_effect = False
        elif isinstance(node, intrinsic.Intrinsic):
            self.global_effect = node.global_effects
        elif isinstance(node, ast.alias):
            self.global_effect = False
        elif isinstance(node, str):
            self.global_effect = False
        elif isinstance(node, intrinsic.Class):
            self.global_effect = False
        elif isinstance(node, intrinsic.UnboundValueType):
            self.global_effect = True  # conservative choice
        else:
            raise NotImplementedError

# Compute the intrinsic effects only once
IntrinsicGlobalEffects = {}

def save_global_effects(module):
    """ Recursively save globals effect for all functions. """
    for intr in module.values():
        if isinstance(intr, dict):  # Submodule case
            save_global_effects(intr)
        else:
            fe = FunctionEffect(intr)
            IntrinsicGlobalEffects[intr] = fe
            if isinstance(intr, intrinsic.Class):
                save_global_effects(intr.fields)

for module in MODULES.values():
    save_global_effects(module)

class GlobalEffectsHelper(ModuleAnalysis[Aliases, GlobalDeclarations, Intrinsics]):

    """Add a flag on each function that updates a global variable."""

    ResultType = DiGraph
    def __init__(self):
        super().__init__()
        self.node_to_functioneffect = dict()

    def prepare(self, node):
        """
        Initialise globals effects as this analyse is inter-procedural.

        Initialisation done for Pythonic functions and default value set for
        user defined functions.
        """
        super().prepare(node)

        for i in self.intrinsics:
            fe = IntrinsicGlobalEffects[i]
            self.node_to_functioneffect[i] = fe
            self.result.add_node(fe)

        for n in self.global_declarations.values():
            fe = FunctionEffect(n)
            self.node_to_functioneffect[n] = fe
            self.result.add_node(fe)

        self.node_to_functioneffect[intrinsic.UnboundValue] = \
            FunctionEffect(intrinsic.UnboundValue)

    def visit_FunctionDef(self, node):
        self.current_function = self.node_to_functioneffect[node]
        assert self.current_function in self.result
        self.generic_visit(node)

    def visit_Print(self, _):
        self.current_function.global_effect = True

    def visit_Call(self, node):
        # try to get all aliases of the function, if possible
        # else use [] as a fallback
        func_aliases = self.aliases[node.func]
        # expand argument if any
        func_aliases = reduce(
            # all funcs
            lambda x, y: x + (list(self.node_to_functioneffect.keys())
                              if isinstance(y, ast.Name) else [y]),
            func_aliases,
            list())
        for func_alias in func_aliases:
            # special hook for bound functions
            if isinstance(func_alias, ast.Call):
                fake_call = ast.Call(func_alias.args[0],
                                     func_alias.args[1:], [])
                self.visit(fake_call)
                continue

            # conservative choice
            if func_alias not in self.node_to_functioneffect:
                func_alias = intrinsic.UnboundValue

            func_alias = self.node_to_functioneffect[func_alias]
            self.result.add_edge(self.current_function, func_alias)
        self.generic_visit(node)


class GlobalEffects(ModuleAnalysis[GlobalEffectsHelper]):

    """Add a flag on each function that updates a global variable."""

    ResultType = dict

    def visit_Module(self, node):
        result = self.global_effects_helper
        keep_going = True
        while keep_going:
            keep_going = False
            for function in result:
                if function.global_effect:
                    for pred in result.predecessors(function):
                        if not pred.global_effect:
                            keep_going = pred.global_effect = True
        self.result = {f.func for f in result if f.global_effect}
