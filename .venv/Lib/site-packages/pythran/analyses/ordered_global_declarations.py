""" OrderedGlobalDeclarations orders all global functions. """

from pythran.analyses.aliases import StrictAliases
from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.passmanager import ModuleAnalysis

import gast as ast


class OrderedGlobalDeclarationsHelper(ModuleAnalysis[StrictAliases, GlobalDeclarations]):
    '''Order all global functions according to their callgraph depth'''
    ResultType = dict

    def visit_FunctionDef(self, node):
        self.curr = node
        self.result[node] = set()
        self.generic_visit(node)

    def visit_Name(self, node):
        if node in self.strict_aliases:
            for alias in self.strict_aliases[node]:
                if isinstance(alias, ast.FunctionDef):
                    self.result[self.curr].add(alias)
                elif isinstance(alias, ast.Call):  # this is a bind
                    for alias in self.strict_aliases[alias.args[0]]:
                        if alias in self.global_declarations:
                            self.result[self.curr].add(alias)


class OrderedGlobalDeclarations(ModuleAnalysis[OrderedGlobalDeclarationsHelper]):
    '''Order all global functions according to their callgraph depth'''
    ResultType = list

    def visit_Module(self, node):
        # compute the weight of each function
        # the weight of a function is the number functions it references
        result = self.ordered_global_declarations_helper
        old_count = -1
        new_count = 0
        # iteratively propagate weights
        while new_count != old_count:
            for v in result.values():
                v.update(*[result[f] for f in v])
            old_count = new_count
            new_count = sum(len(value) for value in result.values())
        # return functions, the one with the greatest weight first
        self.result = sorted(result.keys(), reverse=True,
                             key=lambda s: len(result[s]))
        return self.result
