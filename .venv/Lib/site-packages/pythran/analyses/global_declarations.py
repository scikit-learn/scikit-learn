""" GlobalDeclarations gathers top-level declarations. """

import gast as ast
from pythran.passmanager import ModuleAnalysis

class GlobalDeclarations(ModuleAnalysis):

    """ Gather all kind of identifier defined at global scope.

    >>> import gast as ast
    >>> from pythran import passmanager
    >>> from pythran.analyses import GlobalDeclarations
    >>> node = ast.parse('''
    ... import math
    ... import math as maths
    ... from math import cos
    ... c = 12
    ... def foo(a):
    ...     b = a + 1''')
    >>> pm = passmanager.PassManager("test")
    >>> sorted(pm.gather(GlobalDeclarations, node).keys())
    ['c', 'cos', 'foo', 'math', 'maths']

    """

    ResultType = dict

    def visit_FunctionDef(self, node):
        """ Import module define a new variable name. """
        self.result[node.name] = node

    def visit_Import(self, node):
        for alias in node.names:
            self.result[alias.asname or alias.name] = alias

    def visit_ImportFrom(self, node):
        for alias in node.names:
            self.result[alias.asname or alias.name] = alias

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Store):
            self.result[node.id] = node


class NonlocalDeclarations(ModuleAnalysis):

    """ Transitively gather nonlocal declarations, per function

    >>> import gast as ast
    >>> from pythran import passmanager
    >>> from pythran.analyses import NonlocalDeclarations
    >>> node = ast.parse('''
    ... def foo(a):
    ...   def t():
    ...     def bar():
    ...       nonlocal a
    ...       a = 1
    ...     bar()
    ...   t()''')
    >>> pm = passmanager.PassManager("test")
    >>> [(n.name, l) for n, l in pm.gather(NonlocalDeclarations, node).items()]
    [('foo', set()), ('t', {'a'}), ('bar', {'a'})]

    """

    ResultType = dict

    def __init__(self):
        super().__init__()
        self.context = []
        self.locals = []
        self.nested = [set()]

    def visit_FunctionDef(self, node):
        self.nested[-1].add(node)

        self.nested.append(set())
        self.locals.append(set())

        self.context.append(node)

        self.result[node] = set()
        self.generic_visit(node)

        self.context.pop()

        locals_ = self.locals.pop()
        nested = self.nested.pop()

        for f in nested:
            self.result[node].update(self.result[f] - locals_)

    def visit_Name(self, node):
        if isinstance(node.ctx, (ast.Store, ast.Param)):
            self.locals[-1].add(node.id)

    def visit_Nonlocal(self, node):
        self.result[self.context[-1]].update(node.names)
