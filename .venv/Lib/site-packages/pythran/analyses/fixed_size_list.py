"""
Whether a list usage makes it a candidate for fixed-size-list

This could be a type information, but it seems easier to implement it that way
"""
from pythran.passmanager import FunctionAnalysis
from pythran.tables import MODULES
from pythran.analyses import Aliases, Ancestors
from pythran.analyses.use_def_chain import DefUseChains
from pythran.analyses import ArgumentEffects

import gast as ast


class FixedSizeList(FunctionAnalysis[Aliases, DefUseChains, Ancestors, ArgumentEffects]):

    ResultType = set

    def is_fixed_size_list_def(self, node):
        if isinstance(node, ast.List):
            return True

        if not isinstance(node, ast.Call):
            return False

        return all(alias == MODULES['builtins']['list']
                   for alias in self.aliases[node.func])

    def is_safe_call(self, node, index):
        func_aliases = list(self.aliases[node])
        for alias in func_aliases:
            if isinstance(alias, ast.Call):
                if not self.is_safe_call(alias.args[0],
                                         index + len(alias.args) - 1):
                    return False

            if alias in self.argument_effects:
                func_aes = self.argument_effects[alias]
                if func_aes[index]:
                    return False
        return True

    def is_safe_use(self, use):
        parent = self.ancestors[use.node][-1]

        OK = ast.Subscript, ast.BinOp
        if isinstance(parent, OK):
            return True

        if isinstance(parent, ast.Call):
            n = parent.args.index(use.node)
            return self.is_safe_call(parent.func, n)

        return False

    def visit_Assign(self, node):
        self.generic_visit(node)
        if not node.value:
            return
        if not self.is_fixed_size_list_def(node.value):
            return

        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            def_ = self.def_use_chains.chains[target]
            if any(not self.is_safe_use(u) for u in def_.users()):
                break
            if not isinstance(target, ast.Name):
                continue
            if sum(1 for d in self.def_use_chains.locals[self.ctx.function]
                   if d.name() == target.id) > 1:
                break
        else:
            self.result.add(node.value)

    visit_AnnAssign = visit_Assign

    def visit_Call(self, node):
        self.generic_visit(node)
        for i, arg in enumerate(node.args):
            if not self.is_fixed_size_list_def(arg):
                continue
            if self.is_safe_call(node.func, i):
                self.result.add(arg)
