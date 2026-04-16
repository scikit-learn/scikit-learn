""" ArgumentEffects computes write effect on arguments. """

from pythran.analyses.aliases import Aliases
from pythran.analyses.intrinsics import Intrinsics
from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.passmanager import ModuleAnalysis
from pythran.tables import MODULES
from pythran.graph import DiGraph
# FIXME: investigate why we need to import it that way
from pythran import intrinsic

import gast as ast
from functools import reduce


class FunctionEffects(object):
    def __init__(self, node):
        self.func = node
        if isinstance(node, ast.FunctionDef):
            self.update_effects = [False] * len(node.args.args)
        elif isinstance(node, intrinsic.Intrinsic):
            self.update_effects = [isinstance(x, intrinsic.UpdateEffect)
                                   for x in node.argument_effects]
        elif isinstance(node, ast.alias):
            self.update_effects = []
        elif isinstance(node, intrinsic.Class):
            self.update_effects = []
        else:
            raise NotImplementedError


# Compute the intrinsic effects only once
IntrinsicArgumentEffects = {}


def save_function_effect(module):
    """ Recursively save function effect for pythonic functions. """
    for intr in module.values():
        if isinstance(intr, dict):  # Submodule case
            save_function_effect(intr)
        else:
            fe = FunctionEffects(intr)
            IntrinsicArgumentEffects[intr] = fe
            if isinstance(intr, intrinsic.Class):
                save_function_effect(intr.fields)

for module in MODULES.values():
    save_function_effect(module)


class ArgumentEffectsHelper(ModuleAnalysis[Aliases, GlobalDeclarations, Intrinsics]):

    """Gathers inter-procedural effects on function arguments."""

    ResultType = DiGraph

    def __init__(self):
        # There's an edge between src and dest if a parameter of dest is
        # modified by src
        super().__init__()
        self.node_to_functioneffect = {}

    def prepare(self, node):
        """
        Initialise arguments effects as this analyse is inter-procedural.

        Initialisation done for Pythonic functions and default value set for
        user defined functions.
        """
        super().prepare(node)
        for i in self.intrinsics:
            fe = IntrinsicArgumentEffects[i]
            self.node_to_functioneffect[i] = fe
            self.result.add_node(fe)

        for n in self.global_declarations.values():
            fe = FunctionEffects(n)
            self.node_to_functioneffect[n] = fe
            self.result.add_node(fe)

    def argument_index(self, node):
        while isinstance(node, ast.Subscript):
            node = node.value
        for node_alias in self.aliases[node]:
            while isinstance(node_alias, ast.Subscript):
                node_alias = node_alias.value
            if node_alias in self.current_arguments:
                return self.current_arguments[node_alias]
            if node_alias in self.current_subscripted_arguments:
                return self.current_subscripted_arguments[node_alias]
        return -1

    def visit_FunctionDef(self, node):
        self.current_function = self.node_to_functioneffect[node]
        self.current_arguments = {arg: i
                                  for i, arg
                                  in enumerate(node.args.args)}
        self.current_subscripted_arguments = dict()
        assert self.current_function in self.result
        self.generic_visit(node)

    def visit_For(self, node):
        ai = self.argument_index(node.iter)
        if ai >= 0:
            self.current_subscripted_arguments[node.target] = ai
        self.generic_visit(node)

    def visit_AugAssign(self, node):
        n = self.argument_index(node.target)
        if n >= 0:
            self.current_function.update_effects[n] = True
        self.generic_visit(node)

    def visit_Assign(self, node):
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for t in targets:
            if isinstance(t, ast.Subscript):
                n = self.argument_index(t)
                if n >= 0:
                    self.current_function.update_effects[n] = True
        self.generic_visit(node)

    visit_AnnAssign = visit_Assign

    def visit_Call(self, node):
        for i, arg in enumerate(node.args):
            n = self.argument_index(arg)
            if n >= 0:
                func_aliases = self.aliases[node.func]

                # pessimistic case: no alias found
                if func_aliases is None:
                    self.current_function.update_effects[n] = True
                    continue

                # expand argument if any
                func_aliases = reduce(
                    lambda x, y: x + (
                        # all functions
                        list(self.node_to_functioneffect.keys())
                        if (isinstance(y, ast.Name) and
                            self.argument_index(y) >= 0)
                        else [y]),
                    func_aliases,
                    list())

                for func_alias in func_aliases:
                    # special hook for binded functions
                    if isinstance(func_alias, ast.Call):
                        bound_name = func_alias.args[0].id
                        func_alias = self.global_declarations[bound_name]
                    if func_alias is intrinsic.UnboundValue:
                        continue

                    if func_alias not in self.node_to_functioneffect:
                        continue

                    if func_alias is MODULES['functools']['partial']:
                        base_func_aliases = self.aliases[node.args[0]]
                        fe = self.node_to_functioneffect[func_alias]
                        if len(base_func_aliases) == 1:
                            base_func_alias = next(iter(base_func_aliases))
                            fe = self.node_to_functioneffect.get(
                                base_func_alias,
                                fe)
                    else:
                        fe = self.node_to_functioneffect[func_alias]

                    if not self.result.has_edge(fe, self.current_function):
                        self.result.add_edge(
                            fe,
                            self.current_function,
                            effective_parameters=[],
                            formal_parameters=[])
                    edge = self.result.edges[fe, self.current_function]
                    edge["effective_parameters"].append(n)
                    edge["formal_parameters"].append(i)
        self.generic_visit(node)

class ArgumentEffects(ModuleAnalysis[ArgumentEffectsHelper]):

    """Gathers inter-procedural effects on function arguments."""

    ResultType = dict

    def visit_Module(self, node):
        result = self.argument_effects_helper
        candidates = set(result)
        while candidates:
            function = candidates.pop()
            for ue in enumerate(function.update_effects):
                update_effect_idx, update_effect = ue
                if not update_effect:
                    continue
                for pred in result.successors(function):
                    edge = result.edges[function, pred]
                    for fp in enumerate(edge["formal_parameters"]):
                        i, formal_parameter_idx = fp
                        # propagate the impurity backward if needed.
                        # Afterward we may need another graph iteration
                        ith_effectiv = edge["effective_parameters"][i]
                        if(formal_parameter_idx == update_effect_idx and
                           not pred.update_effects[ith_effectiv]):
                            pred.update_effects[ith_effectiv] = True
                            candidates.add(pred)
        self.result = {f.func: f.update_effects for f in result}
