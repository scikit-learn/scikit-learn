""" ExpandBuiltins replaces builtins by their full paths. """

from pythran.analyses import Globals, Locals
from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError
from pythran.tables import MODULES

import builtins
import gast as ast


class ExpandBuiltins(Transformation[Locals, Globals]):

    """
    Expands all builtins into full paths.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(): return list()")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ExpandBuiltins, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        return builtins.list()
    """

    def visit_Name(self, node):
        s = node.id
        if s in ('None', 'True', 'False'):
            self.update = True
            return ast.Constant(getattr(builtins, s), None)
        if(isinstance(node.ctx, ast.Load) and
           s not in self.locals[node] and
           s not in self.globals and
           s in MODULES['builtins']):
            if s == 'getattr':
                raise PythranSyntaxError("You fool! Trying a getattr?", node)
            self.update = True
            return ast.Attribute(
                ast.Name('builtins', ast.Load(), None, None),
                s,
                node.ctx)
        else:
            return node
