"""
UnshadowParameters prevents the shadow parameter phenomenon
"""

from pythran.analyses import Identifiers
from pythran.passmanager import Transformation

import gast as ast


class UnshadowParameters(Transformation[Identifiers]):
    '''
    Prevents parameter shadowing by creating new variable.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(a): a = 1")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(UnshadowParameters, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(a):
        a_ = a
        a_ = 1
    '''

    def visit_FunctionDef(self, node):
        self.argsid = {arg.id for arg in node.args.args}
        self.renaming = {}
        [self.visit(n) for n in node.body]
        # do it twice to make sure all renaming are done
        [self.visit(n) for n in node.body]
        for k, v in self.renaming.items():
            node.body.insert(
                0,
                ast.Assign(
                    [ast.Name(v, ast.Store(), None, None)],
                    ast.Name(k, ast.Load(), None, None),
                    None)
                )
        self.update |= bool(self.renaming)
        return node

    def update_name(self, node):
        if isinstance(node, ast.Name) and node.id in self.argsid:
            if node.id not in self.renaming:
                new_name = node.id
                while new_name in self.identifiers:
                    new_name = new_name + "_"
                self.renaming[node.id] = new_name

    def visit_Assign(self, node):
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            self.update_name(target)
        try:
            self.generic_visit(node)
        except AttributeError:
            pass
        return node
    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        self.update_name(node.target)
        return self.generic_visit(node)

    def visit_Name(self, node):
        if node.id in self.renaming:
            node.id = self.renaming[node.id]
        return node
