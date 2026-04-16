"""
NormalizeTypeIs transforms type(x) == y into isinstance(y, type(x)).

FIXME: This is a small deviation from the original semantic of type... is as
isinstance honors subtyping.
"""

from pythran.passmanager import Transformation

import gast as ast
from functools import reduce


def istypecall(node):
    # FIXME: could use alias analysis instead, but this pattern is rather
    # idiomatic.
    if not isinstance(node, ast.Call):
        return False
    return getattr(node.func, 'attr', None) == "type"


class NormalizeTypeIs(Transformation):
    '''

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("""
    ... def foo(y):
    ...  return builtins.int == builtins.type(y)""")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeTypeIs, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(y):
        return builtins.isinstance(y, builtins.int)
    '''

    def visit_Compare(self, node):
        self.generic_visit(node)
        if not all(isinstance(op, (ast.Is, ast.Eq)) for op in node.ops):
            return node

        comparators = node.comparators + [node.left]
        comparators_are_typecall = [istypecall(expr) for expr in comparators]
        if not any(comparators_are_typecall):
            return node

        self.update = True

        typecall_index = comparators_are_typecall.index(True)
        checked_value = comparators[typecall_index].args[0]
        typecalls = [expr for i, expr in enumerate(comparators)
                     if i != typecall_index]
        return reduce(lambda x, y:
                      ast.BinOp(x, ast.BitAnd(),
                                ast.Call(
                                    ast.Attribute(
                                        ast.Name('builtins', ast.Load(),
                                                 None, None),
                                        'isinstance',
                                        ast.Load()),
                                    [checked_value, y],
                                    [])
                                ), typecalls[1:],
                      ast.Call(
                          ast.Attribute(ast.Name('builtins', ast.Load(), None, None),
                                        'isinstance',
                                        ast.Load()),
                          [checked_value, typecalls[0]],
                          [])
                      )
