""" NormalizeException simplifies try blocks. """

from pythran.passmanager import Transformation

import gast as ast


class NormalizeException(Transformation):
    '''
    Transform else statement in try except block in nested try except.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("try:print('t')\\nexcept: print('x')\\n\
else: print('e')")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeException, node)
    >>> print(pm.dump(backend.Python, node))
    try:
        print('t')
        try:
            print('e')
        except:
            pass
    except:
        print('x')
    '''
    # FIXME : The transformation is incorrect. Else statement should propagate
    # exception
    def visit_Try(self, node):
        if node.orelse:
            node.body.append(
                ast.Try(
                    node.orelse,
                    [ast.ExceptHandler(None, None, [ast.Pass()])],
                    [], []
                    )
                )
            node.orelse = []
            self.update = True
        if node.finalbody:
            node.body.extend(node.finalbody)
            node.finalbody.append(ast.Raise(None, None))
            self.update = True
            node = ast.Try(
                node.body,
                [ast.ExceptHandler(None, None, node.finalbody)],
                [], [])
            node.finalbody = []
        return node
