"""
PotentialIterator finds if it is possible to use an iterator.
"""

from pythran.analyses.aliases import Aliases
from pythran.analyses.argument_read_once import ArgumentReadOnce
from pythran.passmanager import NodeAnalysis

import gast as ast


class PotentialIterator(NodeAnalysis[Aliases, ArgumentReadOnce]):
    """Find whether an expression can be replaced with an iterator."""
    ResultType = set

    def visit_For(self, node):
        self.result.add(node.iter)
        self.generic_visit(node)

    def visit_Compare(self, node):
        if isinstance(node.ops[0], (ast.In, ast.NotIn)):
            self.result.update(node.comparators)
        self.generic_visit(node)

    def visit_Call(self, node):
        for i, arg in enumerate(node.args):

            def isReadOnce(f, i):
                return (f in self.argument_read_once and
                        self.argument_read_once[f][i] <= 1)
            if all(isReadOnce(alias, i) for alias in self.aliases[node.func]):
                self.result.add(arg)
        self.generic_visit(node)
