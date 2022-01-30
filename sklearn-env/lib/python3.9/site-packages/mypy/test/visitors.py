"""Visitor classes pulled out from different tests

These are here because we don't currently support having interpreted
classes subtype compiled ones but pytest grabs the python file
even if the test was compiled.

"""

from typing import Set

from mypy.nodes import (
    NameExpr, TypeVarExpr, CallExpr, Expression, MypyFile, AssignmentStmt, IntExpr
)
from mypy.traverser import TraverserVisitor

from mypy.treetransform import TransformVisitor
from mypy.types import Type


# from testtypegen
class SkippedNodeSearcher(TraverserVisitor):
    def __init__(self) -> None:
        self.nodes: Set[Expression] = set()
        self.is_typing = False

    def visit_mypy_file(self, f: MypyFile) -> None:
        self.is_typing = f.fullname == 'typing' or f.fullname == 'builtins'
        super().visit_mypy_file(f)

    def visit_assignment_stmt(self, s: AssignmentStmt) -> None:
        if s.type or ignore_node(s.rvalue):
            for lvalue in s.lvalues:
                if isinstance(lvalue, NameExpr):
                    self.nodes.add(lvalue)
        super().visit_assignment_stmt(s)

    def visit_name_expr(self, n: NameExpr) -> None:
        self.skip_if_typing(n)

    def visit_int_expr(self, n: IntExpr) -> None:
        self.skip_if_typing(n)

    def skip_if_typing(self, n: Expression) -> None:
        if self.is_typing:
            self.nodes.add(n)


def ignore_node(node: Expression) -> bool:
    """Return True if node is to be omitted from test case output."""

    # We want to get rid of object() expressions in the typing module stub
    # and also TypeVar(...) expressions. Since detecting whether a node comes
    # from the typing module is not easy, we just to strip them all away.
    if isinstance(node, TypeVarExpr):
        return True
    if isinstance(node, NameExpr) and node.fullname == 'builtins.object':
        return True
    if isinstance(node, NameExpr) and node.fullname == 'builtins.None':
        return True
    if isinstance(node, CallExpr) and (ignore_node(node.callee) or
                                       node.analyzed):
        return True

    return False


# from testtransform
class TypeAssertTransformVisitor(TransformVisitor):
    def type(self, type: Type) -> Type:
        assert type is not None
        return type
