# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2023 Intel Corporation

"""Rust CFG parser.

Rust uses its `cfg()` format in cargo.
https://doc.rust-lang.org/reference/conditional-compilation.html

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
import typing as T


from ..mesonlib import MesonBugException

if T.TYPE_CHECKING:
    _T = T.TypeVar('_T')
    _LEX_TOKEN = T.Tuple['TokenType', T.Optional[str]]
    _LEX_STREAM = T.Iterator[_LEX_TOKEN]
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
    CFG = enum.auto()


def lexer(raw: str) -> _LEX_STREAM:
    """Lex a cfg() expression.

    :param raw: The raw cfg() expression
    :return: An iterable of tokens
    """
    start: int = 0
    is_string: bool = False
    for i, s in enumerate(raw):
        if s.isspace() or s in {')', '(', ',', '=', '"'}:
            val = raw[start:i]
            start = i + 1
            if s == '"' and is_string:
                yield (TokenType.STRING, val)
                is_string = False
                continue
            elif val == 'any':
                yield (TokenType.ANY, None)
            elif val == 'all':
                yield (TokenType.ALL, None)
            elif val == 'not':
                yield (TokenType.NOT, None)
            elif val == 'cfg':
                yield (TokenType.CFG, None)
            elif val:
                yield (TokenType.IDENTIFIER, val)

            if s == '(':
                yield (TokenType.LPAREN, None)
            elif s == ')':
                yield (TokenType.RPAREN, None)
            elif s == ',':
                yield (TokenType.COMMA, None)
            elif s == '=':
                yield (TokenType.EQUAL, None)
            elif s == '"':
                is_string = True
    val = raw[start:]
    if val:
        # This should always be an identifier
        yield (TokenType.IDENTIFIER, val)


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

    lhs: Identifier
    rhs: String


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

    if token is TokenType.IDENTIFIER:
        assert value
        id_ = Identifier(value)
        if ntoken is TokenType.EQUAL:
            next(ast)
            (token, value), _ = next(ast)
            assert token is TokenType.STRING
            assert value is not None
            return Equal(id_, String(value))
        return id_
    elif token in {TokenType.ANY, TokenType.ALL}:
        type_ = All if token is TokenType.ALL else Any
        args: T.List[IR] = []
        (token, value), n_stream = next(ast)
        assert token is TokenType.LPAREN
        if n_stream and n_stream[0] == TokenType.RPAREN:
            return type_(args)
        while True:
            args.append(_parse(ast))
            (token, value), _ = next(ast)
            if token is TokenType.RPAREN:
                break
            assert token is TokenType.COMMA
        return type_(args)
    elif token in {TokenType.NOT, TokenType.CFG}:
        is_not = token is TokenType.NOT
        (token, value), _ = next(ast)
        assert token is TokenType.LPAREN
        arg = _parse(ast)
        (token, value), _ = next(ast)
        assert token is TokenType.RPAREN
        return Not(arg) if is_not else arg
    else:
        raise MesonBugException(f'Unhandled Cargo token:{token} {value}')


def parse(ast: _LEX_STREAM) -> IR:
    """Parse the tokenized list into Meson AST.

    :param ast: An iterable of Tokens
    :return: An mparser Node to be used as a conditional
    """
    ast_i: _LEX_STREAM_AH = lookahead(ast)
    return _parse(ast_i)


def _eval_cfg(ir: IR, cfgs: T.Dict[str, str]) -> bool:
    if isinstance(ir, Identifier):
        return ir.value in cfgs
    elif isinstance(ir, Equal):
        return cfgs.get(ir.lhs.value) == ir.rhs.value
    elif isinstance(ir, Not):
        return not _eval_cfg(ir.value, cfgs)
    elif isinstance(ir, Any):
        return any(_eval_cfg(i, cfgs) for i in ir.args)
    elif isinstance(ir, All):
        return all(_eval_cfg(i, cfgs) for i in ir.args)
    else:
        raise MesonBugException(f'Unhandled Cargo cfg IR: {ir}')


def eval_cfg(raw: str, cfgs: T.Dict[str, str]) -> bool:
    return _eval_cfg(parse(lexer(raw)), cfgs)
