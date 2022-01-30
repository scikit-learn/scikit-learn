"""
    sphinx.util.tags
    ~~~~~~~~~~~~~~~~

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from typing import Iterator, List

from jinja2 import nodes
from jinja2.environment import Environment
from jinja2.nodes import Node
from jinja2.parser import Parser

env = Environment()


class BooleanParser(Parser):
    """
    Only allow condition exprs and/or/not operations.
    """

    def parse_compare(self) -> Node:
        node: Node
        token = self.stream.current
        if token.type == 'name':
            if token.value in ('true', 'false', 'True', 'False'):
                node = nodes.Const(token.value in ('true', 'True'),
                                   lineno=token.lineno)
            elif token.value in ('none', 'None'):
                node = nodes.Const(None, lineno=token.lineno)
            else:
                node = nodes.Name(token.value, 'load', lineno=token.lineno)
            next(self.stream)
        elif token.type == 'lparen':
            next(self.stream)
            node = self.parse_expression()
            self.stream.expect('rparen')
        else:
            self.fail("unexpected token '%s'" % (token,), token.lineno)
        return node


class Tags:
    def __init__(self, tags: List[str] = None) -> None:
        self.tags = dict.fromkeys(tags or [], True)

    def has(self, tag: str) -> bool:
        return tag in self.tags

    __contains__ = has

    def __iter__(self) -> Iterator[str]:
        return iter(self.tags)

    def add(self, tag: str) -> None:
        self.tags[tag] = True

    def remove(self, tag: str) -> None:
        self.tags.pop(tag, None)

    def eval_condition(self, condition: str) -> bool:
        # exceptions are handled by the caller
        parser = BooleanParser(env, condition, state='variable')
        expr = parser.parse_expression()
        if not parser.stream.eos:
            raise ValueError('chunk after expression')

        def eval_node(node: Node) -> bool:
            if isinstance(node, nodes.CondExpr):
                if eval_node(node.test):
                    return eval_node(node.expr1)
                else:
                    return eval_node(node.expr2)
            elif isinstance(node, nodes.And):
                return eval_node(node.left) and eval_node(node.right)
            elif isinstance(node, nodes.Or):
                return eval_node(node.left) or eval_node(node.right)
            elif isinstance(node, nodes.Not):
                return not eval_node(node.node)
            elif isinstance(node, nodes.Name):
                return self.tags.get(node.name, False)
            else:
                raise ValueError('invalid node, check parsing')

        return eval_node(expr)
