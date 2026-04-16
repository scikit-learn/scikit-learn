""" NormalizeIsNone detects is None patterns. """

from pythran.passmanager import Transformation
from pythran.analyses import Ancestors
from pythran.syntax import PythranSyntaxError
from functools import reduce

import gast as ast


def is_none(expr):
    return isinstance(expr, ast.Constant) and expr.value is None


def is_is_none(expr):
    if not isinstance(expr, ast.Compare):
        return None

    if len(expr.ops) != 1:
        exprs = [expr.left] + expr.comparators
        if any(is_none(expr) for expr in exprs):
            raise PythranSyntaxError("is None in complex condition", expr)
        return None

    if not isinstance(expr.ops[0], (ast.Eq, ast.Is)):
        return None

    if is_none(expr.left):
        return expr.comparators[0]

    if is_none(expr.comparators[0]):
        return expr.left

    return None


def is_is_not_none(expr):
    if not isinstance(expr, ast.Compare):
        return None

    if len(expr.ops) != 1:
        exprs = [expr.left] + expr.comparators
        if any(is_none(expr) for expr in exprs):
            raise PythranSyntaxError("is None in complex condition", expr)
        return None

    if not isinstance(expr.ops[0], (ast.NotEq, ast.IsNot)):
        return None

    if is_none(expr.left):
        return expr.comparators[0]

    if is_none(expr.comparators[0]):
        return expr.left

    return None


class NormalizeIsNone(Transformation[Ancestors]):

    table = {ast.And: ast.BitAnd, ast.Or: ast.BitOr}

    @staticmethod
    def match_is_none(node):
        noned_var = is_is_none(node)
        if noned_var is None:
            noned_var = is_is_not_none(node)
            negated = noned_var is not None
        else:
            negated = False
        return noned_var, negated

    def visit_BoolOp(self, node):
        values = list(node.values)
        self.generic_visit(node)
        if any(x != y for x, y in zip(values, node.values)):
            self.update = True
            expr = reduce(lambda x, y:
                          ast.BinOp(x,
                                    NormalizeIsNone.table[type(node.op)](), y),
                          node.values)
            return expr
        else:
            return node

    def visit_Compare(self, node):
        self.generic_visit(node)
        noned_var, negated = self.match_is_none(node)
        if noned_var is None:
            return node
        call = ast.Call(
            ast.Attribute(
                ast.Attribute(
                    ast.Name('builtins', ast.Load(), None, None),
                    'pythran',
                    ast.Load()
                ),
                'is_none',
                ast.Load()),
            [noned_var], [])

        self.update = True

        if negated:
            return ast.UnaryOp(ast.Not(), call)
        else:
            return call
