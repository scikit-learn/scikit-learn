from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import jinja2.environment
import jinja2.nodes
import jinja2.parser

from sphinx.deprecation import RemovedInSphinx90Warning

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Literal

_ENV = jinja2.environment.Environment()


class BooleanParser(jinja2.parser.Parser):
    """Only allow conditional expressions and binary operators."""

    def parse_compare(self) -> jinja2.nodes.Expr:
        node: jinja2.nodes.Expr
        token = self.stream.current
        if token.type == 'name':
            if token.value in {'true', 'True'}:
                node = jinja2.nodes.Const(True, lineno=token.lineno)
            elif token.value in {'false', 'False'}:
                node = jinja2.nodes.Const(False, lineno=token.lineno)
            elif token.value in {'none', 'None'}:
                node = jinja2.nodes.Const(None, lineno=token.lineno)
            else:
                node = jinja2.nodes.Name(token.value, 'load', lineno=token.lineno)
            next(self.stream)
        elif token.type == 'lparen':
            next(self.stream)
            node = self.parse_expression()
            self.stream.expect('rparen')
        else:
            self.fail(f"unexpected token '{token}'", token.lineno)
        return node


class Tags:
    def __init__(self, tags: Sequence[str] = ()) -> None:
        self._tags = set(tags or ())
        self._condition_cache: dict[str, bool] = {}

    def __str__(self) -> str:
        return f'{self.__class__.__name__}({", ".join(sorted(self._tags))})'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({tuple(sorted(self._tags))})'

    def __iter__(self) -> Iterator[str]:
        return iter(self._tags)

    def __contains__(self, tag: str) -> bool:
        return tag in self._tags

    def has(self, tag: str) -> bool:
        return tag in self._tags

    def add(self, tag: str) -> None:
        self._tags.add(tag)

    def remove(self, tag: str) -> None:
        self._tags.discard(tag)

    @property
    def tags(self) -> dict[str, Literal[True]]:
        warnings.warn(
            'Tags.tags is deprecated, use methods on Tags.',
            RemovedInSphinx90Warning,
            stacklevel=2,
        )
        return dict.fromkeys(self._tags, True)

    def eval_condition(self, condition: str) -> bool:
        """Evaluate a boolean condition.

        Only conditional expressions and binary operators (and, or, not)
        are permitted, and operate on tag names, where truthy values mean
        the tag is present and vice versa.
        """
        if condition in self._condition_cache:
            return self._condition_cache[condition]

        # exceptions are handled by the caller
        parser = BooleanParser(_ENV, condition, state='variable')
        expr = parser.parse_expression()
        if not parser.stream.eos:
            msg = 'chunk after expression'
            raise ValueError(msg)

        evaluated = self._condition_cache[condition] = self._eval_node(expr)
        return evaluated

    def _eval_node(self, node: jinja2.nodes.Node | None) -> bool:
        if isinstance(node, jinja2.nodes.CondExpr):
            if self._eval_node(node.test):
                return self._eval_node(node.expr1)
            else:
                return self._eval_node(node.expr2)
        elif isinstance(node, jinja2.nodes.And):
            return self._eval_node(node.left) and self._eval_node(node.right)
        elif isinstance(node, jinja2.nodes.Or):
            return self._eval_node(node.left) or self._eval_node(node.right)
        elif isinstance(node, jinja2.nodes.Not):
            return not self._eval_node(node.node)
        elif isinstance(node, jinja2.nodes.Name):
            return node.name in self._tags
        else:
            msg = 'invalid node, check parsing'
            raise ValueError(msg)
