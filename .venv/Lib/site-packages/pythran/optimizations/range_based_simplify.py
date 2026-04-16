''' Simplify expressions based on range information when possible'''

from pythran.analyses import RangeValues
from pythran.passmanager import Transformation

import gast as ast
from math import isinf
from copy import deepcopy


class RangeBasedSimplify(Transformation[RangeValues]):
    '''
    Simplify expressions based on range analysis

    >>> import gast as ast
    >>> from pythran import passmanager, backend

    >>> node = ast.parse("def any():\\n for x in builtins.range(10): y=x%8")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RangeBasedSimplify, node)
    >>> print(pm.dump(backend.Python, node))
    def any():
        for x in builtins.range(10):
            y = (x if (x < 8) else (x - 8))

    >>> node = ast.parse("def any(): x = 1 or 2; return 3 == x")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RangeBasedSimplify, node)
    >>> print(pm.dump(backend.Python, node))
    def any():
        x = (1 or 2)
        return 0

    >>> node = ast.parse("def a(i): x = 1,1,2; return x[2], x[0 if i else 1]")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RangeBasedSimplify, node)
    >>> print(pm.dump(backend.Python, node))
    def a(i):
        x = (1, 1, 2)
        return (2, 1)
    '''

    def visit_OMPDirective(self, node):
        return node

    def visit_BinOp(self, node):
        node = self.generic_visit(node)
        if not isinstance(node.op, ast.Mod):
            return node

        right_range = self.range_values[node.right]
        left_range = self.range_values[node.left]

        if right_range.low < 0 or isinf(right_range.high):
            return node

        if left_range.low < -right_range.low:
            return node
        if left_range.high > right_range.high * 2:
            return node

        cleft0, cleft1 = deepcopy(node.left), deepcopy(node.left)
        cright = deepcopy(node.right)
        self.update = True
        return ast.IfExp(ast.Compare(node.left, [ast.Lt()], [node.right]),
                         cleft0,
                         ast.BinOp(cleft1, ast.Sub(), cright))

    def visit_range(self, node):
        range_value = self.range_values[node]
        if isinf(range_value.high):
            return self.generic_visit(node)
        elif range_value.low == range_value.high:
            self.update = True
            return ast.Constant(range_value.low, None)
        else:
            return self.generic_visit(node)

    visit_Compare = visit_range

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load):
            return self.visit_range(node)
        return self.generic_visit(node)

    visit_Subscript = visit_Name
