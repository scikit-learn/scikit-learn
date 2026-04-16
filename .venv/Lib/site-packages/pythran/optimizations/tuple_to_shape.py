""" TupleToShap transforms some Tuple node into shape nodes when relevant. """

from pythran.analyses import Aliases
from pythran.tables import MODULES
from pythran.passmanager import Transformation
from pythran.utils import pythran_builtin_attr

import gast as ast

patterns = (MODULES['numpy']['full'],
            MODULES['numpy']['ones'],
            MODULES['numpy']['zeros'],
            MODULES['numpy']['empty'],
            )
reshape_patterns = MODULES['numpy']['ndarray']['reshape'],


def istuple(node):
    return isinstance(node, ast.Tuple)


def toshape(node):
    b = pythran_builtin_attr("make_shape")
    return ast.Call(b, node.elts, [])


class TupleToShape(Transformation[Aliases]):

    """
    Replace tuple nodes by shape when relevant

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(n): import numpy; return numpy.ones((n,4))")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(TupleToShape, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n):
        import numpy
        return numpy.ones(builtins.pythran.make_shape(n, 4))
    """
    def visit_Call(self, node):
        func_aliases = self.aliases.get(node.func, None)
        if func_aliases is not None:
            if func_aliases.issubset(patterns):
                if istuple(node.args[0]):
                    self.update = True
                    node.args[0] = toshape(node.args[0])
            elif func_aliases.issubset(reshape_patterns):
                if len(node.args) > 2:
                    self.update = True
                    node.args[1:] = [toshape(ast.List(node.args[1:],
                                                      ast.Load()))]
        return self.generic_visit(node)
