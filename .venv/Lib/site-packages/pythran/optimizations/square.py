""" Replaces **2 by a call to numpy.square. """

from pythran.passmanager import Transformation
from pythran.analyses.ast_matcher import ASTMatcher, AST_any
from pythran.conversion import mangle
from pythran.utils import isnum

import gast as ast
import copy


class Square(Transformation):

    """
    Replaces **2 by a call to numpy.square.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('a**2')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(Square, node)
    >>> print(pm.dump(backend.Python, node))
    import numpy as __pythran_import_numpy
    __pythran_import_numpy.square(a)
    >>> node = ast.parse('__pythran_import_numpy.power(a,2)')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(Square, node)
    >>> print(pm.dump(backend.Python, node))
    import numpy as __pythran_import_numpy
    __pythran_import_numpy.square(a)
    """

    POW_PATTERN = ast.BinOp(AST_any(), ast.Pow(), ast.Constant(2, None))
    POWER_PATTERN = ast.Call(
        ast.Attribute(ast.Name(mangle('numpy'), ast.Load(), None, None),
                      'power',
                      ast.Load()),
        [AST_any(), ast.Constant(2, None)],
        [])

    def replace(self, value):
        self.update = self.need_import = True
        module_name = ast.Name(mangle('numpy'), ast.Load(), None, None)
        return ast.Call(ast.Attribute(module_name, 'square', ast.Load()),
                        [value], [])

    def visit_Module(self, node):
        self.need_import = False
        self.generic_visit(node)
        if self.need_import:
            import_alias = ast.alias(name='numpy', asname=mangle('numpy'))
            importIt = ast.Import(names=[import_alias])
            node.body.insert(0, importIt)
        return node

    def expand_pow(self, node, n):
        if n == 0:
            return ast.Constant(1, None)
        elif n == 1:
            return node
        else:
            node_square = self.replace(node)
            node_pow = self.expand_pow(node_square, n >> 1)
            if n & 1:
                return ast.BinOp(node_pow, ast.Mult(), copy.deepcopy(node))
            else:
                return node_pow

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if ASTMatcher(Square.POW_PATTERN).match(node):
            return self.replace(node.left)
        elif isinstance(node.op, ast.Pow) and isnum(node.right):
            n = node.right.value
            if int(n) == n and n > 0:
                return self.expand_pow(node.left, int(n))
            else:
                return node
        else:
            return node

    def visit_Call(self, node):
        self.generic_visit(node)
        if ASTMatcher(Square.POWER_PATTERN).match(node):
            return self.replace(node.args[0])
        else:
            return node
