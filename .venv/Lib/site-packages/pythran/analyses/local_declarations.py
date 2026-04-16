"""
LocalNameDeclarations gathers name of declarations local to a node.
LocalNodeDeclarations gathers node of declarations local to a node.
"""

from pythran.passmanager import NodeAnalysis

import gast as ast


class LocalNodeDeclarations(NodeAnalysis):

    """
    Gathers all local symbols from a function.

    It should not be use from outside a function, but can be used on a function
    (but in that case, parameters are not taken into account)

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... def foo(a):
    ...     b = a + 1''')
    >>> pm = passmanager.PassManager("test")
    >>> [name.id for name in pm.gather(LocalNodeDeclarations, node)]
    ['b']
    >>> node = ast.parse('''
    ... for c in range(n):
    ...     b = a + 1''')
    >>> pm = passmanager.PassManager("test")
    >>> sorted([name.id for name in pm.gather(LocalNodeDeclarations, node)])
    ['b', 'c']
    """

    ResultType = set

    def visit_Name(self, node):
        """ Any node with Store context is a new declaration. """
        if isinstance(node.ctx, ast.Store):
            self.result.add(node)


class LocalNameDeclarations(NodeAnalysis):

    """
    Gathers all local identifiers from a node.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... def foo(a):
    ...     b = a + 1''')
    >>> pm = passmanager.PassManager("test")
    >>> sorted(pm.gather(LocalNameDeclarations, node))
    ['a', 'b', 'foo']
    """

    """ Initialize empty set as the result. """
    ResultType = set

    def visit_Name(self, node):
        """ Any node with Store or Param context is a new identifier. """
        if isinstance(node.ctx, (ast.Store, ast.Param)):
            self.result.add(node.id)

    def visit_FunctionDef(self, node):
        """ Function name is a possible identifier. """
        self.result.add(node.name)
        self.generic_visit(node)
