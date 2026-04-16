"""
ExtractDocStrings fills a dictionary with doc strings for each function.
"""

from pythran.passmanager import Transformation
from pythran.utils import isstr

import gast as ast


class ExtractDocStrings(Transformation):
    '''
    Extract Python Doc Strings, removing them from the AST and putting them in
    a dictionary for later use.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(): 'my doc is cool' ; pass")
    >>> pm = passmanager.PassManager("test")
    >>> _ = pm.apply(ExtractDocStrings, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        pass
    '''

    def __init__(self):
        super(ExtractDocStrings, self).__init__()
        self.docstrings = dict()

    def run(self, node):
        super(ExtractDocStrings, self).run(node)
        return self.docstrings

    def visit_Expr(self, node):
        'Remove other top-level strings'
        if isstr(node.value):
            return None
        return node

    def visit_documented_node(self, key, node):
        if node.body:
            first_stmt = node.body[0]
            if isinstance(first_stmt, ast.Expr):
                if isstr(first_stmt.value):
                    self.update = True
                    docstring = first_stmt.value.value
                    self.docstrings[key] = docstring
                    node.body.pop(0)
        return self.generic_visit(node)

    def visit_Module(self, node):
        return self.visit_documented_node(None, node)

    def visit_FunctionDef(self, node):
        return self.visit_documented_node(node.name, node)
