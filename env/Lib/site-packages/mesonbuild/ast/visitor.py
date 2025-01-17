# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool
from __future__ import annotations

import typing as T
from itertools import zip_longest

if T.TYPE_CHECKING:
    from .. import mparser

class AstVisitor:
    def __init__(self) -> None:
        pass

    def visit_default_func(self, node: mparser.BaseNode) -> None:
        pass

    def visit_BooleanNode(self, node: mparser.BooleanNode) -> None:
        self.visit_default_func(node)

    def visit_IdNode(self, node: mparser.IdNode) -> None:
        self.visit_default_func(node)

    def visit_NumberNode(self, node: mparser.NumberNode) -> None:
        self.visit_default_func(node)

    def visit_StringNode(self, node: mparser.StringNode) -> None:
        self.visit_default_func(node)

    def visit_ContinueNode(self, node: mparser.ContinueNode) -> None:
        self.visit_default_func(node)

    def visit_BreakNode(self, node: mparser.BreakNode) -> None:
        self.visit_default_func(node)

    def visit_SymbolNode(self, node: mparser.SymbolNode) -> None:
        self.visit_default_func(node)

    def visit_WhitespaceNode(self, node: mparser.WhitespaceNode) -> None:
        self.visit_default_func(node)

    def visit_ArrayNode(self, node: mparser.ArrayNode) -> None:
        self.visit_default_func(node)
        node.args.accept(self)

    def visit_DictNode(self, node: mparser.DictNode) -> None:
        self.visit_default_func(node)
        node.args.accept(self)

    def visit_EmptyNode(self, node: mparser.EmptyNode) -> None:
        self.visit_default_func(node)

    def visit_OrNode(self, node: mparser.OrNode) -> None:
        self.visit_default_func(node)
        node.left.accept(self)
        node.right.accept(self)

    def visit_AndNode(self, node: mparser.AndNode) -> None:
        self.visit_default_func(node)
        node.left.accept(self)
        node.right.accept(self)

    def visit_ComparisonNode(self, node: mparser.ComparisonNode) -> None:
        self.visit_default_func(node)
        node.left.accept(self)
        node.right.accept(self)

    def visit_ArithmeticNode(self, node: mparser.ArithmeticNode) -> None:
        self.visit_default_func(node)
        node.left.accept(self)
        node.right.accept(self)

    def visit_NotNode(self, node: mparser.NotNode) -> None:
        self.visit_default_func(node)
        node.value.accept(self)

    def visit_CodeBlockNode(self, node: mparser.CodeBlockNode) -> None:
        self.visit_default_func(node)
        for i in node.lines:
            i.accept(self)

    def visit_IndexNode(self, node: mparser.IndexNode) -> None:
        self.visit_default_func(node)
        node.iobject.accept(self)
        node.index.accept(self)

    def visit_MethodNode(self, node: mparser.MethodNode) -> None:
        self.visit_default_func(node)
        node.source_object.accept(self)
        node.name.accept(self)
        node.args.accept(self)

    def visit_FunctionNode(self, node: mparser.FunctionNode) -> None:
        self.visit_default_func(node)
        node.func_name.accept(self)
        node.args.accept(self)

    def visit_AssignmentNode(self, node: mparser.AssignmentNode) -> None:
        self.visit_default_func(node)
        node.var_name.accept(self)
        node.value.accept(self)

    def visit_PlusAssignmentNode(self, node: mparser.PlusAssignmentNode) -> None:
        self.visit_default_func(node)
        node.var_name.accept(self)
        node.value.accept(self)

    def visit_ForeachClauseNode(self, node: mparser.ForeachClauseNode) -> None:
        self.visit_default_func(node)
        for varname in node.varnames:
            varname.accept(self)
        node.items.accept(self)
        node.block.accept(self)

    def visit_IfClauseNode(self, node: mparser.IfClauseNode) -> None:
        self.visit_default_func(node)
        for i in node.ifs:
            i.accept(self)
        node.elseblock.accept(self)

    def visit_UMinusNode(self, node: mparser.UMinusNode) -> None:
        self.visit_default_func(node)
        node.value.accept(self)

    def visit_IfNode(self, node: mparser.IfNode) -> None:
        self.visit_default_func(node)
        node.condition.accept(self)
        node.block.accept(self)

    def visit_ElseNode(self, node: mparser.ElseNode) -> None:
        self.visit_default_func(node)
        node.block.accept(self)

    def visit_TernaryNode(self, node: mparser.TernaryNode) -> None:
        self.visit_default_func(node)
        node.condition.accept(self)
        node.trueblock.accept(self)
        node.falseblock.accept(self)

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        self.visit_default_func(node)
        for i in node.arguments:
            i.accept(self)
        for key, val in node.kwargs.items():
            key.accept(self)
            val.accept(self)

    def visit_ParenthesizedNode(self, node: mparser.ParenthesizedNode) -> None:
        self.visit_default_func(node)
        node.inner.accept(self)

class FullAstVisitor(AstVisitor):
    """Visit all nodes, including Symbol and Whitespaces"""

    def enter_node(self, node: mparser.BaseNode) -> None:
        pass

    def exit_node(self, node: mparser.BaseNode) -> None:
        if node.whitespaces:
            node.whitespaces.accept(self)

    def visit_default_func(self, node: mparser.BaseNode) -> None:
        self.enter_node(node)
        self.exit_node(node)

    def visit_UnaryOperatorNode(self, node: mparser.UnaryOperatorNode) -> None:
        self.enter_node(node)
        node.operator.accept(self)
        node.value.accept(self)
        self.exit_node(node)

    def visit_BinaryOperatorNode(self, node: mparser.BinaryOperatorNode) -> None:
        self.enter_node(node)
        node.left.accept(self)
        node.operator.accept(self)
        node.right.accept(self)
        self.exit_node(node)

    def visit_ArrayNode(self, node: mparser.ArrayNode) -> None:
        self.enter_node(node)
        node.lbracket.accept(self)
        node.args.accept(self)
        node.rbracket.accept(self)
        self.exit_node(node)

    def visit_DictNode(self, node: mparser.DictNode) -> None:
        self.enter_node(node)
        node.lcurl.accept(self)
        node.args.accept(self)
        node.rcurl.accept(self)
        self.exit_node(node)

    def visit_OrNode(self, node: mparser.OrNode) -> None:
        self.visit_BinaryOperatorNode(node)

    def visit_AndNode(self, node: mparser.AndNode) -> None:
        self.visit_BinaryOperatorNode(node)

    def visit_ComparisonNode(self, node: mparser.ComparisonNode) -> None:
        self.visit_BinaryOperatorNode(node)

    def visit_ArithmeticNode(self, node: mparser.ArithmeticNode) -> None:
        self.visit_BinaryOperatorNode(node)

    def visit_NotNode(self, node: mparser.NotNode) -> None:
        self.visit_UnaryOperatorNode(node)

    def visit_CodeBlockNode(self, node: mparser.CodeBlockNode) -> None:
        self.enter_node(node)
        if node.pre_whitespaces:
            node.pre_whitespaces.accept(self)
        for i in node.lines:
            i.accept(self)
        self.exit_node(node)

    def visit_IndexNode(self, node: mparser.IndexNode) -> None:
        self.enter_node(node)
        node.iobject.accept(self)
        node.lbracket.accept(self)
        node.index.accept(self)
        node.rbracket.accept(self)
        self.exit_node(node)

    def visit_MethodNode(self, node: mparser.MethodNode) -> None:
        self.enter_node(node)
        node.source_object.accept(self)
        node.dot.accept(self)
        node.name.accept(self)
        node.lpar.accept(self)
        node.args.accept(self)
        node.rpar.accept(self)
        self.exit_node(node)

    def visit_FunctionNode(self, node: mparser.FunctionNode) -> None:
        self.enter_node(node)
        node.func_name.accept(self)
        node.lpar.accept(self)
        node.args.accept(self)
        node.rpar.accept(self)
        self.exit_node(node)

    def visit_AssignmentNode(self, node: mparser.AssignmentNode) -> None:
        self.enter_node(node)
        node.var_name.accept(self)
        node.operator.accept(self)
        node.value.accept(self)
        self.exit_node(node)

    def visit_PlusAssignmentNode(self, node: mparser.PlusAssignmentNode) -> None:
        self.visit_AssignmentNode(node)

    def visit_ForeachClauseNode(self, node: mparser.ForeachClauseNode) -> None:
        self.enter_node(node)
        node.foreach_.accept(self)
        for varname, comma in zip_longest(node.varnames, node.commas):
            varname.accept(self)
            if comma is not None:
                comma.accept(self)
        node.colon.accept(self)
        node.items.accept(self)
        node.block.accept(self)
        node.endforeach.accept(self)
        self.exit_node(node)

    def visit_IfClauseNode(self, node: mparser.IfClauseNode) -> None:
        self.enter_node(node)
        for i in node.ifs:
            i.accept(self)
        node.elseblock.accept(self)
        node.endif.accept(self)
        self.exit_node(node)

    def visit_UMinusNode(self, node: mparser.UMinusNode) -> None:
        self.visit_UnaryOperatorNode(node)

    def visit_IfNode(self, node: mparser.IfNode) -> None:
        self.enter_node(node)
        node.if_.accept(self)
        node.condition.accept(self)
        node.block.accept(self)
        self.exit_node(node)

    def visit_ElseNode(self, node: mparser.ElseNode) -> None:
        self.enter_node(node)
        node.else_.accept(self)
        node.block.accept(self)
        self.exit_node(node)

    def visit_TernaryNode(self, node: mparser.TernaryNode) -> None:
        self.enter_node(node)
        node.condition.accept(self)
        node.questionmark.accept(self)
        node.trueblock.accept(self)
        node.colon.accept(self)
        node.falseblock.accept(self)
        self.exit_node(node)

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        self.enter_node(node)
        commas_iter = iter(node.commas)

        for arg in node.arguments:
            arg.accept(self)
            try:
                comma = next(commas_iter)
                comma.accept(self)
            except StopIteration:
                pass

        assert len(node.colons) == len(node.kwargs)
        for (key, val), colon in zip(node.kwargs.items(), node.colons):
            key.accept(self)
            colon.accept(self)
            val.accept(self)
            try:
                comma = next(commas_iter)
                comma.accept(self)
            except StopIteration:
                pass

        self.exit_node(node)

    def visit_ParenthesizedNode(self, node: mparser.ParenthesizedNode) -> None:
        self.enter_node(node)
        node.lpar.accept(self)
        node.inner.accept(self)
        node.rpar.accept(self)
        self.exit_node(node)
