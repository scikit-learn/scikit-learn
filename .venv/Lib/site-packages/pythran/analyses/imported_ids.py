""" ImportedIds gathers identifiers imported by a node. """

from pythran.analyses.globals_analysis import Globals
from pythran.analyses.locals_analysis import Locals
from pythran.passmanager import NodeAnalysis
import pythran.metadata as md

import gast as ast


class ImportedIds(NodeAnalysis[Globals, Locals]):

    """
    Gather ids referenced by a node and not declared locally.

    >>> import gast as ast
    >>> from pythran import passmanager
    >>> from pythran.analyses import ImportedIds
    >>> node = ast.parse('''
    ... def foo():
    ...   def t():
    ...     nonlocal g
    ...     g = k
    ...   t()''')
    >>> pm = passmanager.PassManager("test")
    >>> sorted(pm.gather(ImportedIds, node))
    ['g', 'k']
    """

    ResultType = set

    def __init__(self):
        super().__init__()
        self.current_locals = set()
        self.current_nonlocals = set()
        self.in_augassign = False

    def visit_Name(self, node):
        if node.id in self.current_nonlocals:
            self.result.add(node.id)
        elif isinstance(node.ctx, ast.Store) and not self.in_augassign:
            self.current_locals.add(node.id)
        elif (node.id not in self.visible_globals and
              node.id not in self.current_locals):
            self.result.add(node.id)

    def visit_Nonlocal(self, node):
        self.current_nonlocals.update(node.names)

    def visit_FunctionDef(self, node):
        self.current_locals.add(node.name)
        current_locals = self.current_locals.copy()
        self.current_locals.update(arg.id for arg in node.args.args)
        for stmt in node.body:
            self.visit(stmt)
        self.current_locals = current_locals

    def visit_AnyComp(self, node):
        current_locals = self.current_locals.copy()
        for generator in node.generators:
            self.visit(generator)
        self.visit(node.elt)
        self.current_locals = current_locals

    visit_ListComp = visit_AnyComp
    visit_SetComp = visit_AnyComp
    visit_DictComp = visit_AnyComp
    visit_GeneratorExp = visit_AnyComp

    def visit_Assign(self, node):
        # order matter as an assignation
        # is evaluated before being assigned
        md.visit(self, node)
        if node.value:
            self.visit(node.value)
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            self.visit(target)

    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        self.in_augassign = True
        self.generic_visit(node)
        self.in_augassign = False

    def visit_Lambda(self, node):
        current_locals = self.current_locals.copy()
        self.current_locals.update(arg.id for arg in node.args.args)
        self.visit(node.body)
        self.current_locals = current_locals

    def visit_Import(self, node):
        self.current_locals.update(alias.name for alias in node.names)

    def visit_StoredTuple(self, node):
        for elt in node.elts:
            if isinstance(elt, ast.Name):
                self.current_locals.add(elt.id)
                continue
            if isinstance(elt, ast.Subscript):
                self.visit(elt)
            if isinstance(elt, ast.Tuple):
                self.visit_StoredTuple(node)

    def visit_Tuple(self, node):
        if isinstance(node.ctx, ast.Load):
            self.generic_visit(node)
        else:
            self.visit_StoredTuple(node)

    visit_List = visit_Tuple

    def visit_ImportFrom(self, node):
        self.current_locals.update(alias.name for alias in node.names)

    def visit_Attribute(self, node):
        pass

    def prepare(self, node):
        super(ImportedIds, self).prepare(node)
        self.visible_globals = set(self.globals) - self.locals[node]
