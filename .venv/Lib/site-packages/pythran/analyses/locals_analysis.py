"""
Locals computes the value of locals()
"""

from pythran.passmanager import ModuleAnalysis
import pythran.metadata as md

import gast as ast


class Locals(ModuleAnalysis):
    """
    Statically compute the value of locals() before each statement

    Yields a dictionary binding every node to the set of variable names defined
    *before* this node.

    Following snippet illustrates its behavior:
    >>> import gast as ast
    >>> from pythran import passmanager
    >>> pm = passmanager.PassManager('test')
    >>> code = '''
    ... def b(n):
    ...     m = n + 1
    ...     def b(n):
    ...         return n + 1
    ...     return b(m)'''
    >>> tree = ast.parse(code)
    >>> l = pm.gather(Locals, tree)
    >>> sorted(l[tree.body[0].body[0]])
    ['n']
    >>> sorted(l[tree.body[0].body[1]])
    ['b', 'm', 'n']
    """

    ResultType = dict

    def __init__(self):
        super().__init__()
        self.locals = set()
        self.nesting = 0

    def generic_visit(self, node):
        super(Locals, self).generic_visit(node)
        if node not in self.result:
            self.result[node] = self.result[self.expr_parent]

    def store_and_visit(self, node):
        self.expr_parent = node
        self.result[node] = self.locals.copy()
        self.generic_visit(node)

    def visit_Module(self, node):
        self.expr_parent = node
        self.result[node] = self.locals
        self.generic_visit(node)

    def visit_FunctionDef(self, node):
        # top-level OMP statements attached to that function

        md.visit(self, node)
        # special case for nested functions
        if self.nesting:
            self.locals.add(node.name)
        self.nesting += 1
        self.expr_parent = node
        self.result[node] = self.locals.copy()
        parent_locals = self.locals.copy()
        for default in node.args.defaults:
            self.visit(default)
        for arg in node.args.args:
            if arg.annotation:
                self.visit(arg.annotation)
        if node.returns:
            self.visit(node.returns)
        self.locals.update(arg.id for arg in node.args.args)
        for stmt in node.body:
            self.visit(stmt)
        self.locals = parent_locals
        self.nesting -= 1

    def visit_Assign(self, node):
        self.expr_parent = node
        self.result[node] = self.locals.copy()
        md.visit(self, node)
        if node.value:
            self.visit(node.value)
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        self.locals.update(t.id for t in targets
                           if isinstance(t, ast.Name))
        for target in targets:
            self.visit(target)

    def visit_AnnAssign(self, node):
        self.visit_Assign(node)
        self.visit(node.annotation)

    def visit_For(self, node):
        self.expr_parent = node
        self.result[node] = self.locals.copy()
        md.visit(self, node)
        self.visit(node.iter)
        self.locals.add(node.target.id)
        for stmt in node.body:
            self.visit(stmt)
        for stmt in node.orelse:
            self.visit(stmt)

    def visit_Import(self, node):
        self.result[node] = self.locals.copy()
        self.locals.update(alias.name for alias in node.names)

    def visit_ImportFrom(self, node):
        self.result[node] = self.locals.copy()
        self.locals.update(alias.name for alias in node.names)

    def visit_ExceptHandler(self, node):
        self.expr_parent = node
        self.result[node] = self.locals.copy()
        if node.name:
            self.locals.add(node.name.id)
        node.type and self.visit(node.type)
        for stmt in node.body:
            self.visit(stmt)

    # statements that do not define a new variable
    visit_Return = store_and_visit
    visit_Yield = store_and_visit
    visit_Try = store_and_visit
    visit_AugAssign = store_and_visit
    visit_Print = store_and_visit
    visit_While = store_and_visit
    visit_If = store_and_visit
    visit_Raise = store_and_visit
    visit_Assert = store_and_visit
    visit_Expr = store_and_visit
    visit_Pass = store_and_visit
    visit_Break = store_and_visit
    visit_Continue = store_and_visit
