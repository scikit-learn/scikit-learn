""" Replaces a[:] = b by a call to numpy.copyto. """

from pythran.passmanager import Transformation
from pythran.analyses.ast_matcher import ASTMatcher, AST_any
from pythran.conversion import mangle
from pythran.utils import isnum

import gast as ast
import copy


class CopyTo(Transformation):

    """
    Replaces a[:] = b by a call to numpy.copyto.

    This is a slight extension to numpy.copyto as it assumes it also supports
    string and list as first argument.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('a[:] = b')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(CopyTo, node)
    >>> print(pm.dump(backend.Python, node))
    import numpy as __pythran_import_numpy
    __pythran_import_numpy.copyto(a, b)
    >>> node = ast.parse('a[:] = b[:]')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(CopyTo, node)
    >>> print(pm.dump(backend.Python, node))
    import numpy as __pythran_import_numpy
    __pythran_import_numpy.copyto(a, b)
    """

    def isNone(self, node):
        if node is None:
            return True
        return isinstance(node, ast.Constant) and node.value is None

    def is_full_slice(self, node):
        # FIXME: could accept a call to len for node.upper
        return (
                isinstance(node, ast.Slice) and
                (node.lower == 0 or self.isNone(node.lower)) and
                (self.isNone(node.upper)) and
                (self.isNone(node.step) or node.step == 1)
                )

    def is_fully_sliced(self, node):
        if not isinstance(node, ast.Subscript):
            return False
        if not isinstance(node.value, ast.Name):
            return False
        if self.is_full_slice(node.slice):
            return True
        elif isinstance(node.slice, ast.Tuple):
            return all(self.is_full_slice(elt) for elt in node.slice.elts)
        else:
            return False

    def visit_Module(self, node):
        self.generic_visit(node)
        if self.update:
            import_alias = ast.alias(name='numpy', asname=mangle('numpy'))
            importIt = ast.Import(names=[import_alias])
            node.body.insert(0, importIt)
        return node

    def visit_Assign(self, node):
        if not node.value:
            return node
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        if len(targets) != 1:
            return node
        target, = targets
        if not self.is_fully_sliced(target):
            return node
        if self.is_fully_sliced(node.value):
            value = node.value.value
        else:
            value = node.value

        self.update = True
        return ast.Expr(
                ast.Call(
                    ast.Attribute(ast.Name(mangle('numpy'), ast.Load(), None, None),
                                  'copyto',
                                  ast.Load()),
                    [target.value, value],
                    [])
                )
    visit_AnnAssign = visit_Assign


