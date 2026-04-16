""" Replaces **2 by a call to numpy.square. """

from pythran.passmanager import Transformation

import gast as ast

def getid(node):
    if isinstance(node, ast.Attribute):
        return getid(node.value), node.attr
    if isinstance(node, ast.Name):
        return node.id
    return node

class SimplifyExcept(Transformation):

    """
    Remove redundant except clauses

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('try: pass\\nexcept (OSError, OSError): pass')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(SimplifyExcept, node)
    >>> print(pm.dump(backend.Python, node))
    try:
        pass
    except OSError:
        pass
    """

    def visit_ExceptHandler(self, node):
        if isinstance(node.type, ast.Tuple):
            all_ids = {getid(elt) for elt in node.type.elts}
            to_remove = []
            for i, elt in enumerate(node.type.elts):
                eltid = getid(elt)
                if eltid in all_ids:
                    all_ids.remove(eltid)
                else:
                    to_remove.append(i)
            for i in reversed(to_remove):
                node.type.elts.pop(i)
            if len(node.type.elts) == 1:
                node.type = node.type.elts[0]
                self.update = True
            self.update |= bool(to_remove)
        return node
