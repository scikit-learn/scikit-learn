"""
Dependencies lists the functions and types required by a function
"""

from pythran.passmanager import ModuleAnalysis
from pythran.conversion import demangle

import gast as ast
import math


class Dependencies(ModuleAnalysis):
    OpMap = {
        # binop
        ast.Add: ('operator', 'add'),
        ast.Sub: ('operator', 'sub'),
        ast.Mult: ('operator', 'mul'),
        ast.Div: ('operator', 'div'),
        ast.Mod: ('operator', 'mod'),
        ast.Pow: ('operator', 'pow'),
        ast.LShift: ('operator', 'lshift'),
        ast.RShift: ('operator', 'rshift'),
        ast.BitOr: ('operator', 'or_'),
        ast.BitXor: ('operator', 'xor_'),
        ast.BitAnd: ('operator', 'and_'),
        ast.MatMult: ('operator', 'matmul'),
        ast.FloorDiv: ('operator', 'floordiv'),
        # unaryop
        ast.Invert: ('operator', 'invert'),
        ast.Not: ('operator', 'not_'),
        ast.UAdd: ('operator', 'pos'),
        ast.USub: ('operator', 'neg'),
        # cmpop
        ast.Eq: ('operator', 'eq'),
        ast.NotEq: ('operator', 'ne'),
        ast.Lt: ('operator', 'lt'),
        ast.LtE: ('operator', 'le'),
        ast.Gt: ('operator', 'gt'),
        ast.GtE: ('operator', 'ge'),
        ast.Is: ('operator', 'is_'),
        ast.IsNot: ('operator', 'is_not'),
        ast.In: ('operator', 'contains'),
        ast.NotIn: ('operator', 'contains'),
    }

    IOpMap = {
        ast.Add: ('operator', 'iadd'),
        ast.Sub: ('operator', 'isub'),
        ast.Mult: ('operator', 'imul'),
        ast.Div: ('operator', 'idiv'),
        ast.Mod: ('operator', 'imod'),
        ast.Pow: ('operator', 'ipow'),
        ast.LShift: ('operator', 'ilshift'),
        ast.RShift: ('operator', 'irshift'),
        ast.BitOr: ('operator', 'ior'),
        ast.BitXor: ('operator', 'ixor'),
        ast.BitAnd: ('operator', 'iand'),
        ast.MatMult: ('operator', 'imatmul'),
        ast.FloorDiv: ('operator', 'ifloordiv'),
    }

    ResultType = set

    def visit_List(self, node):
        self.result.add(('builtins', 'list'))
        self.generic_visit(node)

    def visit_Tuple(self, node):
        self.result.add(('builtins', 'tuple'))
        self.generic_visit(node)

    def visit_Set(self, node):
        self.result.add(('builtins', 'set'))
        self.generic_visit(node)

    def visit_Dict(self, node):
        self.result.add(('builtins', 'dict'))
        self.generic_visit(node)

    def visit_Slice(self, node):
        self.result.add(('types', 'slice'))
        self.generic_visit(node)

    def visit_And(self, node):
        self.result.add(('builtins', 'pythran', 'and'))
        self.generic_visit(node)

    def visit_Or(self, node):
        self.result.add(('builtins', 'pythran', 'or'))
        self.generic_visit(node)

    def visit_BinOp(self, node):
        self.visit(node.left)
        self.result.add(Dependencies.OpMap[type(node.op)])
        self.visit(node.right)

    def visit_UnaryOp(self, node):
        self.result.add(Dependencies.OpMap[type(node.op)])
        self.visit(node.operand)

    def visit_Compare(self, node):
        self.visit(node.left)
        for op in node.ops:
            self.result.add(Dependencies.OpMap[type(op)])
        for comparator in node.comparators:
            self.visit(comparator)

    def visit_AugAssign(self, node):
        self.visit(node.target)
        # because of the way type inference turns augassign into assign
        self.result.add(Dependencies.OpMap[type(node.op)])
        self.result.add(Dependencies.IOpMap[type(node.op)])
        self.visit(node.value)

    def visit_Print(self, node):
        self.result.add(('builtins', 'print'))
        self.generic_visit(node)

    def visit_Assert(self, node):
        self.result.add(('builtins', 'assert'))
        self.generic_visit(node)

    def visit_Yield(self, node):
        self.result.add(('utils', 'yield'))
        self.generic_visit(node)

    def visit_Constant(self, node):
        if node.value is None:
            self.result.add(('builtins', 'None'))
        elif isinstance(node.value, bytes):
            self.result.add(('types', 'str'))  # FIXME: using str as backend
        elif isinstance(node.value, str):
            self.result.add(('types', 'str'))
        elif isinstance(node.value, complex):
            self.result.add(('types', 'complex'))
        elif math.isnan(node.value):
            self.result.add(('numpy', 'nan'))
        elif math.isinf(node.value):
            self.result.add(('numpy', 'inf'))

    def visit_Attribute(self, node):
        def rec(n):
            if isinstance(n, ast.Name):
                return demangle(n.id),
            elif isinstance(n, ast.Attribute):
                return rec(n.value) + (n.attr,)
        attr = rec(node)

        attr and self.result.add(attr)
