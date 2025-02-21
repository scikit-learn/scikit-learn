# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2023 Intel Corporation

"""Rust CFG parser.

Rust uses its `cfg()` format in cargo.

This may have the following functions:
 - all()
 - any()
 - not()

And additionally is made up of `identifier [ = str]`. Where the str is optional,
so you could have examples like:
```
[target.`cfg(unix)`.dependencies]
[target.'cfg(target_arch = "x86_64")'.dependencies]
[target.'cfg(all(target_arch = "x86_64", target_arch = "x86"))'.dependencies]
```
"""

from __future__ import annotations
import dataclasses
import enum
import functools
import typing as T


from . import builder
from .. import mparser
from ..mesonlib import MesonBugException

if T.TYPE_CHECKING:
    _T = T.TypeVar('_T')
    _LEX_TOKEN = T.Tuple['TokenType', T.Optional[str]]
    _LEX_STREAM = T.Iterable[_LEX_TOKEN]
    _LEX_STREAM_AH = T.Iterator[T.Tuple[_LEX_TOKEN, T.Optional[_LEX_TOKEN]]]


class TokenType(enum.Enum):

    LPAREN = enum.auto()
    RPAREN = enum.auto()
    STRING = enum.auto()
    IDENTIFIER = enum.auto()
    ALL = enum.auto()
    ANY = enum.auto()
    NOT = enum.auto()
    COMMA = enum.auto()
    EQUAL = enum.auto()


def lexer(raw: str) -> _LEX_STREAM:
    """Lex a cfg() expression.

    :param raw: The raw cfg() expression
    :return: An iterable of tokens
    """
    buffer: T.List[str] = []
    is_string: bool = False
    for s in raw:
        if s.isspace() or s in {')', '(', ',', '='} or (s == '"' and buffer):
            val = ''.join(buffer)
            buffer.clear()
            if is_string:
                yield (TokenType.STRING, val)
            elif val == 'any':
                yield (TokenType.ANY, None)
            elif val == 'all':
                yield (TokenType.ALL, None)
            elif val == 'not':
                yield (TokenType.NOT, None)
            elif val:
                yield (TokenType.IDENTIFIER, val)

            if s == '(':
                yield (TokenType.LPAREN, None)
                continue
            elif s == ')':
                yield (TokenType.RPAREN, None)
                continue
            elif s == ',':
                yield (TokenType.COMMA, None)
                continue
            elif s == '=':
                yield (TokenType.EQUAL, None)
                continue
            elif s.isspace():
                continue

        if s == '"':
            is_string = not is_string
        else:
            buffer.append(s)
    if buffer:
        # This should always be an identifier
        yield (TokenType.IDENTIFIER, ''.join(buffer))


def lookahead(iter: T.Iterator[_T]) -> T.Iterator[T.Tuple[_T, T.Optional[_T]]]:
    """Get the current value of the iterable, and the next if possible.

    :param iter: The iterable to look into
    :yield: A tuple of the current value, and, if possible, the next
    :return: nothing
    """
    current: _T
    next_: T.Optional[_T]
    try:
        next_ = next(iter)
    except StopIteration:
        # This is an empty iterator, there's nothing to look ahead to
        return

    while True:
        current = next_
        try:
            next_ = next(iter)
        except StopIteration:
            next_ = None

        yield current, next_

        if next_ is None:
            break


@dataclasses.dataclass
class IR:

    """Base IR node for Cargo CFG."""


@dataclasses.dataclass
class String(IR):

    value: str


@dataclasses.dataclass
class Identifier(IR):

    value: str


@dataclasses.dataclass
class Equal(IR):

    lhs: IR
    rhs: IR


@dataclasses.dataclass
class Any(IR):

    args: T.List[IR]


@dataclasses.dataclass
class All(IR):

    args: T.List[IR]


@dataclasses.dataclass
class Not(IR):

    value: IR


def _parse(ast: _LEX_STREAM_AH) -> IR:
    (token, value), n_stream = next(ast)
    if n_stream is not None:
        ntoken, _ = n_stream
    else:
        ntoken, _ = (None, None)

    stream: T.List[_LEX_TOKEN]
    if token is TokenType.IDENTIFIER:
        if ntoken is TokenType.EQUAL:
            return Equal(Identifier(value), _parse(ast))
    if token is TokenType.STRING:
        return String(value)
    if token is TokenType.EQUAL:
        # In this case the previous caller already has handled the equal
        return _parse(ast)
    if token in {TokenType.ANY, TokenType.ALL}:
        type_ = All if token is TokenType.ALL else Any
        assert ntoken is TokenType.LPAREN
        next(ast)  # advance the iterator to get rid of the LPAREN
        stream = []
        args: T.List[IR] = []
        while token is not TokenType.RPAREN:
            (token, value), _ = next(ast)
            if token is TokenType.COMMA:
                args.append(_parse(lookahead(iter(stream))))
                stream.clear()
            else:
                stream.append((token, value))
        if stream:
            args.append(_parse(lookahead(iter(stream))))
        return type_(args)
    if token is TokenType.NOT:
        next(ast)  # advance the iterator to get rid of the LPAREN
        stream = []
        # Mypy can't figure out that token is overridden inside the while loop
        while token is not TokenType.RPAREN:  # type: ignore
            (token, value), _ = next(ast)
            stream.append((token, value))
        return Not(_parse(lookahead(iter(stream))))

    raise MesonBugException(f'Unhandled Cargo token: {token}')


def parse(ast: _LEX_STREAM) -> IR:
    """Parse the tokenized list into Meson AST.

    :param ast: An iterable of Tokens
    :return: An mparser Node to be used as a conditional
    """
    ast_i: _LEX_STREAM_AH = lookahead(iter(ast))
    return _parse(ast_i)


@functools.singledispatch
def ir_to_meson(ir: T.Any, build: builder.Builder) -> mparser.BaseNode:
    raise NotImplementedError


@ir_to_meson.register
def _(ir: String, build: builder.Builder) -> mparser.BaseNode:
    return build.string(ir.value)


@ir_to_meson.register
def _(ir: Identifier, build: builder.Builder) -> mparser.BaseNode:
    host_machine = build.identifier('host_machine')
    if ir.value == "target_arch":
        return build.method('cpu_family', host_machine)
    elif ir.value in {"target_os", "target_family"}:
        return build.method('system', host_machine)
    elif ir.value == "target_endian":
        return build.method('endian', host_machine)
    raise MesonBugException(f"Unhandled Cargo identifier: {ir.value}")


@ir_to_meson.register
def _(ir: Equal, build: builder.Builder) -> mparser.BaseNode:
    return build.equal(ir_to_meson(ir.lhs, build), ir_to_meson(ir.rhs, build))


@ir_to_meson.register
def _(ir: Not, build: builder.Builder) -> mparser.BaseNode:
    return build.not_(ir_to_meson(ir.value, build))


@ir_to_meson.register
def _(ir: Any, build: builder.Builder) -> mparser.BaseNode:
    args = iter(reversed(ir.args))
    last = next(args)
    cur = build.or_(ir_to_meson(next(args), build), ir_to_meson(last, build))
    for a in args:
        cur = build.or_(ir_to_meson(a, build), cur)
    return cur


@ir_to_meson.register
def _(ir: All, build: builder.Builder) -> mparser.BaseNode:
    args = iter(reversed(ir.args))
    last = next(args)
    cur = build.and_(ir_to_meson(next(args), build), ir_to_meson(last, build))
    for a in args:
        cur = build.and_(ir_to_meson(a, build), cur)
    return cur
