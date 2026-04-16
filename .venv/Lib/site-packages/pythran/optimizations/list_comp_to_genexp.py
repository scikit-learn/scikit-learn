""" ListCompToGenexp transforms list comprehension into genexp. """

from pythran.analyses import PotentialIterator
from pythran.passmanager import Transformation

import gast as ast


class ListCompToGenexp(Transformation[PotentialIterator]):
    '''
    Transforms list comprehension into genexp
    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("""                   \\n\
def foo(l):                                    \\n\
    return builtins.sum(l)                     \\n\
def bar(n):                                    \\n\
    return foo([x for x in builtins.range(n)]) \
""")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ListCompToGenexp, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(l):
        return builtins.sum(l)
    def bar(n):
        return foo((x for x in builtins.range(n)))
    '''

    def visit_ListComp(self, node):
        self.generic_visit(node)
        if node in self.potential_iterator:
            self.update = True
            return ast.GeneratorExp(node.elt, node.generators)
        else:
            return node
