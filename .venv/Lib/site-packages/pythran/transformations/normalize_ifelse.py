""" NormalizeIfElse transform early exit in if into if-else. """

from pythran.analyses import Ancestors
from pythran.passmanager import Transformation

import gast as ast


class NormalizeIfElse(Transformation[Ancestors]):
    '''

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("""
    ... def foo(y):
    ...  if y: return 1
    ...  return 2""")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeIfElse, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(y):
        if y:
            return 1
        else:
            return 2

    >>> node = ast.parse("""
    ... def foo(y):
    ...  if y:
    ...    z = y + 1
    ...    if z:
    ...      return 1
    ...    else:
    ...      return 3
    ...  return 2""")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeIfElse, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(y):
        if y:
            z = (y + 1)
            if z:
                return 1
            else:
                return 3
        else:
            return 2
    '''

    def check_lasts(self, node):
        if isinstance(node, (ast.Return, ast.Break, ast.Return)):
            return True
        if isinstance(node, ast.If):
            if not self.check_lasts(node.body[-1]):
                return False
            return node.orelse and self.check_lasts(node.orelse[-1])

    def visit_If(self, node):
        self.generic_visit(node)
        if not self.check_lasts(node.body[-1]):
            return node
        parent = self.ancestors[node][-1]
        for attr in ('body', 'orelse', 'finalbody'):
            try:
                body = getattr(parent, attr)
                index = body.index(node)
                if index == len(body) - 1:
                    return node
                if not node.orelse:
                    node.orelse = []
                node.orelse.extend(body[index + 1:])
                body[index + 1:] = []
                self.update = True
                return node
            except ValueError:
                continue
            except AttributeError:
                continue
        return node
