"""
PureExpressions detects expressions without side-effects.
"""

from pythran.analyses.aliases import Aliases
from pythran.analyses.argument_effects import ArgumentEffects
from pythran.analyses.global_effects import GlobalEffects
from pythran.analyses.pure_functions import PureFunctions
from pythran.passmanager import ModuleAnalysis
from pythran.intrinsic import Intrinsic

import gast as ast


class PureExpressions(ModuleAnalysis[ArgumentEffects, GlobalEffects, Aliases, PureFunctions]):
    '''Yields the set of pure expressions'''

    ResultType = set

    def visit_FunctionDef(self, node):
        if node in self.pure_functions:
            self.result.add(node)

        # do not visit arguments
        for stmt in node.body:
            self.visit(stmt)

    def generic_visit(self, node):
        is_pure = all([self.visit(x) for x in ast.iter_child_nodes(node)])
        if is_pure:
            self.result.add(node)
        return is_pure

    def visit_Call(self, node):
        # check if all arguments are Pures
        is_pure = all([self.visit(arg) for arg in node.args])

        # check all possible function called
        func_aliases = self.aliases[node.func]
        if func_aliases:
            for func_alias in func_aliases:
                # does the function have a global effect?
                if isinstance(func_alias, Intrinsic):
                    is_pure &= not func_alias.global_effects
                else:
                    is_pure &= func_alias in self.pure_functions

                # does the function have an argument effect ?
                # trivial arguments can be ignored
                if func_alias in self.argument_effects:
                    func_aes = self.argument_effects[func_alias]
                    for arg, ae in zip(node.args, func_aes):
                        if ae:
                            try:
                                ast.literal_eval(arg)
                            except ValueError:
                                is_pure = False
                else:
                    is_pure = False
        else:
            is_pure = False  # conservative choice

        # check for chained call
        is_pure &= self.visit(node.func)
        if is_pure:
            self.result.add(node)
        return is_pure
