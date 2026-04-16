"""
Literals lists nodes that are only literals
"""

from pythran.passmanager import FunctionAnalysis

import gast as ast


class Literals(FunctionAnalysis):
    """
        Store variable that save only Literals (with no construction cost)
    """
    ResultType = set

    def visit_Assign(self, node):
        # list, dict, set and other are not considered as Literals as they have
        # a constructor which may be costly and they can be updated using
        # function call
        if isinstance(node.value, (ast.Constant, ast.Lambda)):
            targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
            targets = [target for target in targets
                       if isinstance(target, ast.Name)]
            self.result.update(targets)
    visit_AnnAssign = visit_Assign
