""" ListToTuple transforms some List node into more Efficient Tuple nodes. """

from pythran.analyses import Aliases, FixedSizeList
from pythran.tables import MODULES
from pythran.passmanager import Transformation
from pythran.utils import path_to_attr

import gast as ast

patterns = {MODULES['numpy']['full'],
            MODULES['numpy']['ones'],
            MODULES['numpy']['zeros'],
            MODULES['numpy']['empty'],
            MODULES['numpy']['concatenate'],
            MODULES['builtins']['tuple'],
           }


def islist(node):
    return isinstance(node, ast.List)


def totuple(node):
    return ast.Tuple(node.elts, node.ctx)


class ListToTuple(Transformation[Aliases, FixedSizeList]):

    """
    Replace list nodes by tuple nodes when possible

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(n): import numpy; return numpy.ones([n,n])")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ListToTuple, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n):
        import numpy
        return numpy.ones((n, n))
    """

    def visit_AugAssign(self, node):
        if not islist(node.value):
            return self.generic_visit(node)
        node.value = totuple(node.value)
        self.update = True
        return self.generic_visit(node)

    def visit_Call(self, node):
        if node.args and islist(node.args[0]):
            func_aliases = self.aliases.get(node.func, set())
            if func_aliases.issubset(patterns):
                self.update = True
                node.args[0] = totuple(node.args[0])
        return self.generic_visit(node)

    def visit_List(self, node):
        self.generic_visit(node)
        if node in self.fixed_size_list:
            return self.convert(node)
        else:
            return node

    def visit_Assign(self, node):
        """
        Replace list calls by static_list calls when possible

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(n):\\n"
        ...                  "    x = builtins.list(n)\\n"
        ...                  "    x[0] = 0\\n"
        ...                  "    return builtins.tuple(x)")
        >>> pm = passmanager.PassManager("test")
        >>> _, node = pm.apply(ListToTuple, node)
        >>> print(pm.dump(backend.Python, node))
        def foo(n):
            x = builtins.pythran.static_list(n)
            x[0] = 0
            return builtins.tuple(x)

        >>> node = ast.parse("def foo(n):\\n"
        ...                  "    x = builtins.list(n)\\n"
        ...                  "    x[0] = 0\\n"
        ...                  "    return x")
        >>> pm = passmanager.PassManager("test")
        >>> _, node = pm.apply(ListToTuple, node)
        >>> print(pm.dump(backend.Python, node))
        def foo(n):
            x = builtins.list(n)
            x[0] = 0
            return x
        """
        self.generic_visit(node)
        if node.value not in self.fixed_size_list:
            return node

        node.value = self.convert(node.value)
        return node

    visit_AnnAssign = visit_Assign

    def convert(self, node):
        self.update = True

        if isinstance(node, ast.Call):
            if not node.args:
                node = ast.Tuple([])
            else:
                node = node.args[0]
        elif isinstance(node, ast.List):
            node = ast.Tuple(node.elts, ast.Load())

        return ast.Call(path_to_attr(('builtins', 'pythran', 'static_list')),
                        [node], [])
