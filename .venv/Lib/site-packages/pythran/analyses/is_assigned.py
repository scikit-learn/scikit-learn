""" Gathers variables that have value modification in the given node. """

from pythran.passmanager import NodeAnalysis

import gast as ast


class IsAssigned(NodeAnalysis):

    """
    Gather variable that change in given node.

    It doesn't check constness as it is use for integer so we don't care about
    arguments effects as it is use by value.
    """

    ResultType = list

    def visit_Name(self, node):
        """ Stored variable have new value. """
        if isinstance(node.ctx, ast.Store):
            self.result.append(node)

    def visit_Tuple(self, node):
        if isinstance(node.ctx, ast.Store):

            def rec(n):
                if isinstance(n, ast.Name):
                    self.result.append(n)
                elif isinstance(n, ast.Tuple):
                    for elt in n.elts:
                        rec(elt)
            rec(node)
