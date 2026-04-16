""" Optimize a[...] = b[...] + c when we have no conflicting aliasing """

from pythran.analyses import InterproceduralAliases
from pythran.passmanager import Transformation

import gast as ast


class FastGExpr(Transformation):

    def as_gexpr(self, node):
        if not isinstance(node, ast.Subscript):
            return None
        if not isinstance(node.slice, ast.Slice):
            if not isinstance(node.slice, ast.Tuple):
                return None
            if not any(isinstance(elt, ast.Slice) for elt in node.slice.elts):
                return None

        if not isinstance(node.value, ast.Name):
            return None

        return node.value, node.slice

    def may_alias(self, gexpr, value):
        if isinstance(value, ast.Constant):
            return False
        if isinstance(value, (ast.List, ast.Tuple)):
            return any(self.may_alias(gexpr, elt) for elt in value.elts)
        if isinstance(value, ast.UnaryOp):
            return self.may_alias(gexpr, value.operand)
        if isinstance(value, ast.BinOp):
            return any(self.may_alias(gexpr, elt) for elt in (value.left,
                                                              value.right))
        if isinstance(value, ast.Subscript):
            if not isinstance(value.value, ast.Name):
                return True
            # Lazy loading of interprocedural_aliases as it's a costly analysis
            if not self.interprocedural_aliases:
                self.interprocedural_aliases = self.passmanager.gather(InterproceduralAliases, self.current_module)
            return not self.interprocedural_aliases[gexpr[0]].isdisjoint(self.interprocedural_aliases[value.value])

        return True

    def visit_Module(self, node):
        self.current_module = node
        self.interprocedural_aliases = None
        new_node = self.generic_visit(node)
        self.interprocedural_aliases = None
        self.current_module = None
        return new_node

    def visit_Assign(self, node):
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        if len(targets) > 1:
            return node

        if not node.value:
            return node

        target, = targets
        value = node.value
        gexpr = self.as_gexpr(target)
        if not gexpr:
            return node

        if self.may_alias(gexpr, value):
            return node

        self.update = True

        func = ast.Attribute(
            value=ast.Attribute(value=ast.Name('builtins', ast.Load(),
                                               None, None),
                                attr="pythran", ctx=ast.Load()),
            attr="restrict_assign", ctx=ast.Load())
        return ast.Expr(ast.Call(func, args=[target, value], keywords=[]))
    visit_AnnAssign = visit_Assign

