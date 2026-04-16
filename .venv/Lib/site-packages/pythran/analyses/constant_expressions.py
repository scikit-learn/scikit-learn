""" ConstantExpressions gathers constant expression.  """

from pythran.analyses.aliases import Aliases
from pythran.analyses.globals_analysis import Globals
from pythran.analyses.locals_analysis import Locals
from pythran.analyses.pure_functions import PureFunctions
from pythran.intrinsic import FunctionIntr
from pythran.passmanager import NodeAnalysis

import gast as ast



class ConstantExpressions(NodeAnalysis[Globals, Locals, Aliases, PureFunctions]):

    """Identify constant expressions."""

    ResultType = set

    def add(self, node):
        self.result.add(node)
        return True

    def visit_BoolOp(self, node):
        return all([self.visit(x) for x in node.values]) and self.add(node)

    def visit_BinOp(self, node):
        rec = all([self.visit(x) for x in (node.left, node.right)])
        return rec and self.add(node)

    def visit_UnaryOp(self, node):
        return self.visit(node.operand) and self.add(node)

    def visit_IfExp(self, node):
        rec = all([self.visit(x) for x in (node.test, node.body, node.orelse)])
        return rec and self.add(node)

    def visit_Compare(self, node):
        rec = all([self.visit(x) for x in ([node.left] + node.comparators)])
        return rec and self.add(node)

    def visit_Call(self, node):
        rec = all([self.visit(x) for x in (node.args + [node.func])])
        return rec and self.add(node)

    visit_Constant = add

    def visit_Subscript(self, node):
        rec = all([self.visit(x) for x in (node.value, node.slice)])
        rec = isinstance(node.ctx, ast.Load) and rec
        return rec and self.add(node)

    def visit_Name(self, node):
        if node in self.aliases:
            # params and store are not constants
            if not isinstance(node.ctx, ast.Load):
                return False
            # if we can alias on multiple value, it is not constant
            elif len(self.aliases[node]) > 1:
                return False
            # if it is not a globals, it depends on variable so it is not
            # constant
            elif node.id not in self.globals:
                return False
            # if it is defined in the current function, it is not constant
            elif node.id in self.locals[node]:
                return False

            def is_function(x):
                return isinstance(x, (FunctionIntr,
                                      ast.FunctionDef,
                                      ast.alias))

            pure_fun = all(alias in self.pure_functions
                           for alias in self.aliases[node])
            return pure_fun
        else:
            return False

    def visit_Attribute(self, node):
        return Aliases.access_path(node).isconst() and self.add(node)

    def visit_Dict(self, node):
        rec = all([self.visit(x) for x in (node.keys + node.values)])
        return rec and self.add(node)

    def visit_List(self, node):
        return all([self.visit(x) for x in node.elts]) and self.add(node)

    visit_Tuple = visit_List
    visit_Set = visit_List

    def visit_Slice(self, node):
        return all([x is None or self.visit(x) for x in
                    (node.lower, node.upper, node.step)]) and self.add(node)

    def visit_Index(self, node):
        return self.visit(node.value) and self.add(node)

    def visit_AnnAssign(self, node):
        self.visit(node.target)
        # skip annotation
        if node.value:
            self.visit(node.value)
