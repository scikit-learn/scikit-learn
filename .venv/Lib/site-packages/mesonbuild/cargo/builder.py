# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2023 Intel Corporation

"""Provides helpers for building AST

This is meant to make building Meson AST from foreign (largely declarative)
build descriptions easier.
"""

from __future__ import annotations
import dataclasses
import typing as T

from .. import mparser

if T.TYPE_CHECKING:
    import builtins


@dataclasses.dataclass
class Builder:

    filename: str

    def _token(self, tid: str, value: mparser.TV_TokenTypes) -> mparser.Token[mparser.TV_TokenTypes]:
        """Create a Token object, but with the line numbers stubbed out.

        :param tid: the token id (such as string, number, etc)
        :param filename: the filename that the token was generated from
        :param value: the value of the token
        :return: A Token object
        """
        return mparser.Token(tid, self.filename, -1, -1, -1, (-1, -1), value)

    def _symbol(self, val: str) -> mparser.SymbolNode:
        return mparser.SymbolNode(self._token('', val))

    def assign(self, value: mparser.BaseNode, varname: str) -> mparser.AssignmentNode:
        return mparser.AssignmentNode(self.identifier(varname), self._symbol('='), value)

    def string(self, value: str) -> mparser.StringNode:
        """Build A StringNode

        :param value: the value of the string
        :return: A StringNode
        """
        return mparser.StringNode(self._token('string', value), escape=False)

    def number(self, value: int) -> mparser.NumberNode:
        """Build A NumberNode

        :param value: the value of the number
        :return: A NumberNode
        """
        return mparser.NumberNode(self._token('number', str(value)))

    def bool(self, value: builtins.bool) -> mparser.BooleanNode:
        """Build A BooleanNode

        :param value: the value of the boolean
        :return: A BooleanNode
        """
        return mparser.BooleanNode(self._token('bool', value))

    def array(self, value: T.List[mparser.BaseNode]) -> mparser.ArrayNode:
        """Build an Array Node

        :param value: A list of nodes to insert into the array
        :return: An ArrayNode built from the arguments
        """
        args = mparser.ArgumentNode(self._token('array', 'unused'))
        args.arguments = value
        return mparser.ArrayNode(self._symbol('['), args, self._symbol(']'))

    def dict(self, value: T.Dict[mparser.BaseNode, mparser.BaseNode]) -> mparser.DictNode:
        """Build an Dictionary Node

        :param value: A dict of nodes to insert into the dictionary
        :return: An DictNode built from the arguments
        """
        args = mparser.ArgumentNode(self._token('dict', 'unused'))
        for key, val in value.items():
            args.set_kwarg_no_check(key, val)
        return mparser.DictNode(self._symbol('{'), args, self._symbol('}'))

    def identifier(self, value: str) -> mparser.IdNode:
        """Build A IdNode

        :param value: the value of the boolean
        :return: A BooleanNode
        """
        return mparser.IdNode(self._token('id', value))

    def method(self, name: str, id_: mparser.BaseNode,
               pos: T.Optional[T.List[mparser.BaseNode]] = None,
               kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
               ) -> mparser.MethodNode:
        """Create a method call.

        :param name: the name of the method
        :param id_: the object to call the method of
        :param pos: a list of positional arguments, defaults to None
        :param kw: a dictionary of keyword arguments, defaults to None
        :return: a method call object
        """
        args = mparser.ArgumentNode(self._token('array', 'unused'))
        if pos is not None:
            args.arguments = pos
        if kw is not None:
            args.kwargs = {self.identifier(k): v for k, v in kw.items()}
        return mparser.MethodNode(id_, self._symbol('.'), self.identifier(name), self._symbol('('), args, self._symbol(')'))

    def function(self, name: str,
                 pos: T.Optional[T.List[mparser.BaseNode]] = None,
                 kw: T.Optional[T.Mapping[str, mparser.BaseNode]] = None,
                 ) -> mparser.FunctionNode:
        """Create a function call.

        :param name: the name of the function
        :param pos: a list of positional arguments, defaults to None
        :param kw: a dictionary of keyword arguments, defaults to None
        :return: a method call object
        """
        args = mparser.ArgumentNode(self._token('array', 'unused'))
        if pos is not None:
            args.arguments = pos
        if kw is not None:
            args.kwargs = {self.identifier(k): v for k, v in kw.items()}
        return mparser.FunctionNode(self.identifier(name), self._symbol('('), args, self._symbol(')'))

    def equal(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
        """Create an equality operation

        :param lhs: The left hand side of the equal
        :param rhs: the right hand side of the equal
        :return: A comparison node
        """
        return mparser.ComparisonNode('==', lhs, self._symbol('=='), rhs)

    def not_equal(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
        """Create an inequality operation

        :param lhs: The left hand side of the "!="
        :param rhs: the right hand side of the "!="
        :return: A comparison node
        """
        return mparser.ComparisonNode('!=', lhs, self._symbol('!='), rhs)

    def in_(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
        """Create an "in" operation

        :param lhs: The left hand side of the "in"
        :param rhs: the right hand side of the "in"
        :return: A comparison node
        """
        return mparser.ComparisonNode('in', lhs, self._symbol('in'), rhs)

    def not_in(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ComparisonNode:
        """Create an "not in" operation

        :param lhs: The left hand side of the "not in"
        :param rhs: the right hand side of the "not in"
        :return: A comparison node
        """
        return mparser.ComparisonNode('not in', lhs, self._symbol('not in'), rhs)

    def or_(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.OrNode:
        """Create and OrNode

        :param lhs: The Left of the Node
        :param rhs: The Right of the Node
        :return: The OrNode
        """
        return mparser.OrNode(lhs, self._symbol('or'), rhs)

    def and_(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.AndNode:
        """Create an AndNode

        :param lhs: The left of the And
        :param rhs: The right of the And
        :return: The AndNode
        """
        return mparser.AndNode(lhs, self._symbol('and'), rhs)

    def not_(self, value: mparser.BaseNode) -> mparser.NotNode:
        """Create a not node

        :param value: The value to negate
        :return: The NotNode
        """
        return mparser.NotNode(self._token('not', ''), self._symbol('not'), value)

    def block(self, lines: T.List[mparser.BaseNode]) -> mparser.CodeBlockNode:
        block = mparser.CodeBlockNode(self._token('node', ''))
        block.lines = lines
        return block

    def plus(self, lhs: mparser.BaseNode, rhs: mparser.BaseNode) -> mparser.ArithmeticNode:
        """Create an addition node

        :param lhs: The left of the addition
        :param rhs: The right of the addition
        :return: The ArithmeticNode
        """
        return mparser.ArithmeticNode('+', lhs, self._symbol('+'), rhs)

    def plusassign(self, value: mparser.BaseNode, varname: str) -> mparser.PlusAssignmentNode:
        """Create a "+=" node

        :param value: The value to add
        :param varname: The variable to assign
        :return: The PlusAssignmentNode
        """
        return mparser.PlusAssignmentNode(self.identifier(varname), self._symbol('+='), value)

    def if_(self, condition: mparser.BaseNode, block: mparser.CodeBlockNode) -> mparser.IfClauseNode:
        """Create a "if" block

        :param condition: The condition
        :param block: Lines inside the condition
        :return: The IfClauseNode
        """
        clause = mparser.IfClauseNode(condition)
        clause.ifs.append(mparser.IfNode(clause, self._symbol('if'), condition, block))
        clause.elseblock = mparser.EmptyNode(-1, -1, self.filename)
        return clause

    def foreach(self, varnames: T.List[str], items: mparser.BaseNode, block: mparser.CodeBlockNode) -> mparser.ForeachClauseNode:
        """Create a "foreach" loop

        :param varnames: Iterator variable names (one for list, two for dict).
        :param items: The list of dict to iterate
        :param block: Lines inside the loop
        :return: The ForeachClauseNode
        """
        varids = [self.identifier(i) for i in varnames]
        commas = [self._symbol(',') for i in range(len(varnames) - 1)]
        return mparser.ForeachClauseNode(self._symbol('foreach'), varids, commas, self._symbol(':'), items, block, self._symbol('endforeach'))
