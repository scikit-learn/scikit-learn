from __future__ import annotations

import sys
import warnings
from typing import TYPE_CHECKING, Any, Union, cast

from docutils import nodes

from sphinx import addnodes
from sphinx.domains.c._ids import _id_prefix, _max_id
from sphinx.util.cfamily import (
    ASTAttributeList,
    ASTBaseBase,
    ASTBaseParenExprList,
    UnsupportedMultiCharacterCharLiteral,
    verify_description_mode,
)

if TYPE_CHECKING:
    from typing import TypeAlias

    from docutils.nodes import Element, Node, TextElement

    from sphinx.domains.c._symbol import Symbol
    from sphinx.environment import BuildEnvironment
    from sphinx.util.cfamily import StringifyTransform

DeclarationType: TypeAlias = Union[  # NoQA: UP007
    "ASTStruct", "ASTUnion", "ASTEnum", "ASTEnumerator",
    "ASTType", "ASTTypeWithInit", "ASTMacro",
]


class ASTBase(ASTBaseBase):
    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))


# Names
################################################################################

class ASTIdentifier(ASTBaseBase):
    def __init__(self, name: str) -> None:
        if not isinstance(name, str) or len(name) == 0:
            raise AssertionError
        self.name = sys.intern(name)
        self.is_anonymous = name[0] == '@'

    # ASTBaseBase already implements this method,
    # but specialising it here improves performance
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTIdentifier):
            return NotImplemented
        return self.name == other.name

    def is_anon(self) -> bool:
        return self.is_anonymous

    # and this is where we finally make a difference between __str__ and the display string

    def __str__(self) -> str:
        return self.name

    def get_display_string(self) -> str:
        return "[anonymous]" if self.is_anonymous else self.name

    def describe_signature(self, signode: TextElement, mode: str, env: BuildEnvironment,
                           prefix: str, symbol: Symbol) -> None:
        # note: slightly different signature of describe_signature due to the prefix
        verify_description_mode(mode)
        if self.is_anonymous:
            node = addnodes.desc_sig_name(text="[anonymous]")
        else:
            node = addnodes.desc_sig_name(self.name, self.name)
        if mode == 'markType':
            targetText = prefix + self.name
            pnode = addnodes.pending_xref('', refdomain='c',
                                          reftype='identifier',
                                          reftarget=targetText, modname=None,
                                          classname=None)
            pnode['c:parent_key'] = symbol.get_lookup_key()
            pnode += node
            signode += pnode
        elif mode == 'lastIsName':
            nameNode = addnodes.desc_name()
            nameNode += node
            signode += nameNode
        elif mode == 'noneIsName':
            signode += node
        else:
            raise Exception('Unknown description mode: %s' % mode)

    @property
    def identifier(self) -> str:
        warnings.warn(
            '`ASTIdentifier.identifier` is deprecated, use `ASTIdentifier.name` instead',
            DeprecationWarning, stacklevel=2,
        )
        return self.name


class ASTNestedName(ASTBase):
    def __init__(self, names: list[ASTIdentifier], rooted: bool) -> None:
        assert len(names) > 0
        self.names = names
        self.rooted = rooted

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNestedName):
            return NotImplemented
        return self.names == other.names and self.rooted == other.rooted

    def __hash__(self) -> int:
        return hash((self.names, self.rooted))

    @property
    def name(self) -> ASTNestedName:
        return self

    def get_id(self, version: int) -> str:
        return '.'.join(str(n) for n in self.names)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = '.'.join(transform(n) for n in self.names)
        if self.rooted:
            return '.' + res
        else:
            return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        # just print the name part, with template args, not template params
        if mode == 'noneIsName':
            if self.rooted:
                unreachable = "Can this happen?"
                raise AssertionError(unreachable)  # TODO
                signode += nodes.Text('.')
            for i in range(len(self.names)):
                if i != 0:
                    unreachable = "Can this happen?"
                    raise AssertionError(unreachable)  # TODO
                    signode += nodes.Text('.')
                n = self.names[i]
                n.describe_signature(signode, mode, env, '', symbol)
        elif mode == 'param':
            assert not self.rooted, str(self)
            assert len(self.names) == 1
            self.names[0].describe_signature(signode, 'noneIsName', env, '', symbol)
        elif mode in ('markType', 'lastIsName', 'markName'):
            # Each element should be a pending xref targeting the complete
            # prefix.
            prefix = ''
            first = True
            names = self.names[:-1] if mode == 'lastIsName' else self.names
            # If lastIsName, then wrap all of the prefix in a desc_addname,
            # else append directly to signode.
            # TODO: also for C?
            #  NOTE: Breathe previously relied on the prefix being in the desc_addname node,
            #       so it can remove it in inner declarations.
            dest = signode
            if mode == 'lastIsName':
                dest = addnodes.desc_addname()
            if self.rooted:
                prefix += '.'
                if mode == 'lastIsName' and len(names) == 0:
                    signode += addnodes.desc_sig_punctuation('.', '.')
                else:
                    dest += addnodes.desc_sig_punctuation('.', '.')
            for i in range(len(names)):
                ident = names[i]
                if not first:
                    dest += addnodes.desc_sig_punctuation('.', '.')
                    prefix += '.'
                first = False
                txt_ident = str(ident)
                if txt_ident != '':
                    ident.describe_signature(dest, 'markType', env, prefix, symbol)
                prefix += txt_ident
            if mode == 'lastIsName':
                if len(self.names) > 1:
                    dest += addnodes.desc_sig_punctuation('.', '.')
                    signode += dest
                self.names[-1].describe_signature(signode, mode, env, '', symbol)
        else:
            raise Exception('Unknown description mode: %s' % mode)


################################################################################
# Expressions
################################################################################

class ASTExpression(ASTBase):
    pass


# Primary expressions
################################################################################

class ASTLiteral(ASTExpression):
    pass


class ASTBooleanLiteral(ASTLiteral):
    def __init__(self, value: bool) -> None:
        self.value = value

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTBooleanLiteral):
            return NotImplemented
        return self.value == other.value

    def __hash__(self) -> int:
        return hash(self.value)

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.value:
            return 'true'
        else:
            return 'false'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        txt = str(self)
        signode += addnodes.desc_sig_keyword(txt, txt)


class ASTNumberLiteral(ASTLiteral):
    def __init__(self, data: str) -> None:
        self.data = data

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNumberLiteral):
            return NotImplemented
        return self.data == other.data

    def __hash__(self) -> int:
        return hash(self.data)

    def _stringify(self, transform: StringifyTransform) -> str:
        return self.data

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        txt = str(self)
        signode += addnodes.desc_sig_literal_number(txt, txt)


class ASTCharLiteral(ASTLiteral):
    def __init__(self, prefix: str, data: str) -> None:
        self.prefix = prefix  # may be None when no prefix
        self.data = data
        decoded = data.encode().decode('unicode-escape')
        if len(decoded) == 1:
            self.value = ord(decoded)
        else:
            raise UnsupportedMultiCharacterCharLiteral(decoded)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTCharLiteral):
            return NotImplemented
        return (
            self.prefix == other.prefix
            and self.value == other.value
        )

    def __hash__(self) -> int:
        return hash((self.prefix, self.value))

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.prefix is None:
            return "'" + self.data + "'"
        else:
            return self.prefix + "'" + self.data + "'"

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        txt = str(self)
        signode += addnodes.desc_sig_literal_char(txt, txt)


class ASTStringLiteral(ASTLiteral):
    def __init__(self, data: str) -> None:
        self.data = data

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTStringLiteral):
            return NotImplemented
        return self.data == other.data

    def __hash__(self) -> int:
        return hash(self.data)

    def _stringify(self, transform: StringifyTransform) -> str:
        return self.data

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        txt = str(self)
        signode += addnodes.desc_sig_literal_string(txt, txt)


class ASTIdExpression(ASTExpression):
    def __init__(self, name: ASTNestedName) -> None:
        # note: this class is basically to cast a nested name as an expression
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTIdExpression):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.name)

    def get_id(self, version: int) -> str:
        return self.name.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.name.describe_signature(signode, mode, env, symbol)


class ASTParenExpr(ASTExpression):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTParenExpr):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return '(' + transform(self.expr) + ')'

    def get_id(self, version: int) -> str:
        return self.expr.get_id(version)  # type: ignore[attr-defined]

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


# Postfix expressions
################################################################################

class ASTPostfixOp(ASTBase):
    pass


class ASTPostfixCallExpr(ASTPostfixOp):
    def __init__(self, lst: ASTParenExprList | ASTBracedInitList) -> None:
        self.lst = lst

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPostfixCallExpr):
            return NotImplemented
        return self.lst == other.lst

    def __hash__(self) -> int:
        return hash(self.lst)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.lst)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.lst.describe_signature(signode, mode, env, symbol)


class ASTPostfixArray(ASTPostfixOp):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPostfixArray):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return '[' + transform(self.expr) + ']'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('[', '[')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(']', ']')


class ASTPostfixInc(ASTPostfixOp):
    def _stringify(self, transform: StringifyTransform) -> str:
        return '++'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('++', '++')


class ASTPostfixDec(ASTPostfixOp):
    def _stringify(self, transform: StringifyTransform) -> str:
        return '--'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('--', '--')


class ASTPostfixMemberOfPointer(ASTPostfixOp):
    def __init__(self, name: ASTNestedName) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPostfixMemberOfPointer):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def _stringify(self, transform: StringifyTransform) -> str:
        return '->' + transform(self.name)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('->', '->')
        self.name.describe_signature(signode, 'noneIsName', env, symbol)


class ASTPostfixExpr(ASTExpression):
    def __init__(self, prefix: ASTExpression, postFixes: list[ASTPostfixOp]) -> None:
        self.prefix = prefix
        self.postFixes = postFixes

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPostfixExpr):
            return NotImplemented
        return self.prefix == other.prefix and self.postFixes == other.postFixes

    def __hash__(self) -> int:
        return hash((self.prefix, self.postFixes))

    def _stringify(self, transform: StringifyTransform) -> str:
        return ''.join([transform(self.prefix), *(transform(p) for p in self.postFixes)])

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.prefix.describe_signature(signode, mode, env, symbol)
        for p in self.postFixes:
            p.describe_signature(signode, mode, env, symbol)


# Unary expressions
################################################################################

class ASTUnaryOpExpr(ASTExpression):
    def __init__(self, op: str, expr: ASTExpression) -> None:
        self.op = op
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTUnaryOpExpr):
            return NotImplemented
        return self.op == other.op and self.expr == other.expr

    def __hash__(self) -> int:
        return hash((self.op, self.expr))

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.op[0] in 'cn':
            return self.op + " " + transform(self.expr)
        else:
            return self.op + transform(self.expr)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.op[0] in 'cn':
            signode += addnodes.desc_sig_keyword(self.op, self.op)
            signode += addnodes.desc_sig_space()
        else:
            signode += addnodes.desc_sig_operator(self.op, self.op)
        self.expr.describe_signature(signode, mode, env, symbol)


class ASTSizeofType(ASTExpression):
    def __init__(self, typ: ASTType) -> None:
        self.typ = typ

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTSizeofType):
            return NotImplemented
        return self.typ == other.typ

    def __hash__(self) -> int:
        return hash(self.typ)

    def _stringify(self, transform: StringifyTransform) -> str:
        return "sizeof(" + transform(self.typ) + ")"

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('sizeof', 'sizeof')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTSizeofExpr(ASTExpression):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTSizeofExpr):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return "sizeof " + transform(self.expr)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('sizeof', 'sizeof')
        signode += addnodes.desc_sig_space()
        self.expr.describe_signature(signode, mode, env, symbol)


class ASTAlignofExpr(ASTExpression):
    def __init__(self, typ: ASTType) -> None:
        self.typ = typ

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTAlignofExpr):
            return NotImplemented
        return self.typ == other.typ

    def __hash__(self) -> int:
        return hash(self.typ)

    def _stringify(self, transform: StringifyTransform) -> str:
        return "alignof(" + transform(self.typ) + ")"

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('alignof', 'alignof')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


# Other expressions
################################################################################

class ASTCastExpr(ASTExpression):
    def __init__(self, typ: ASTType, expr: ASTExpression) -> None:
        self.typ = typ
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTCastExpr):
            return NotImplemented
        return (
            self.typ == other.typ
            and self.expr == other.expr
        )

    def __hash__(self) -> int:
        return hash((self.typ, self.expr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['(']
        res.append(transform(self.typ))
        res.append(')')
        res.append(transform(self.expr))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')
        self.expr.describe_signature(signode, mode, env, symbol)


class ASTBinOpExpr(ASTBase):
    def __init__(self, exprs: list[ASTExpression], ops: list[str]) -> None:
        assert len(exprs) > 0
        assert len(exprs) == len(ops) + 1
        self.exprs = exprs
        self.ops = ops

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTBinOpExpr):
            return NotImplemented
        return (
            self.exprs == other.exprs
            and self.ops == other.ops
        )

    def __hash__(self) -> int:
        return hash((self.exprs, self.ops))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.exprs[0]))
        for i in range(1, len(self.exprs)):
            res.append(' ')
            res.append(self.ops[i - 1])
            res.append(' ')
            res.append(transform(self.exprs[i]))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.exprs[0].describe_signature(signode, mode, env, symbol)
        for i in range(1, len(self.exprs)):
            signode += addnodes.desc_sig_space()
            op = self.ops[i - 1]
            if ord(op[0]) >= ord('a') and ord(op[0]) <= ord('z'):
                signode += addnodes.desc_sig_keyword(op, op)
            else:
                signode += addnodes.desc_sig_operator(op, op)
            signode += addnodes.desc_sig_space()
            self.exprs[i].describe_signature(signode, mode, env, symbol)


class ASTAssignmentExpr(ASTExpression):
    def __init__(self, exprs: list[ASTExpression], ops: list[str]) -> None:
        assert len(exprs) > 0
        assert len(exprs) == len(ops) + 1
        self.exprs = exprs
        self.ops = ops

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTAssignmentExpr):
            return NotImplemented
        return (
            self.exprs == other.exprs
            and self.ops == other.ops
        )

    def __hash__(self) -> int:
        return hash((self.exprs, self.ops))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.exprs[0]))
        for i in range(1, len(self.exprs)):
            res.append(' ')
            res.append(self.ops[i - 1])
            res.append(' ')
            res.append(transform(self.exprs[i]))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.exprs[0].describe_signature(signode, mode, env, symbol)
        for i in range(1, len(self.exprs)):
            signode += addnodes.desc_sig_space()
            op = self.ops[i - 1]
            if ord(op[0]) >= ord('a') and ord(op[0]) <= ord('z'):
                signode += addnodes.desc_sig_keyword(op, op)
            else:
                signode += addnodes.desc_sig_operator(op, op)
            signode += addnodes.desc_sig_space()
            self.exprs[i].describe_signature(signode, mode, env, symbol)


class ASTFallbackExpr(ASTExpression):
    def __init__(self, expr: str) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTFallbackExpr):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return self.expr

    def get_id(self, version: int) -> str:
        return str(self.expr)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += nodes.literal(self.expr, self.expr)


################################################################################
# Types
################################################################################

class ASTTrailingTypeSpec(ASTBase):
    pass


class ASTTrailingTypeSpecFundamental(ASTTrailingTypeSpec):
    def __init__(self, names: list[str]) -> None:
        assert len(names) != 0
        self.names = names

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTrailingTypeSpecFundamental):
            return NotImplemented
        return self.names == other.names

    def __hash__(self) -> int:
        return hash(self.names)

    def _stringify(self, transform: StringifyTransform) -> str:
        return ' '.join(self.names)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        first = True
        for n in self.names:
            if not first:
                signode += addnodes.desc_sig_space()
            else:
                first = False
            signode += addnodes.desc_sig_keyword_type(n, n)


class ASTTrailingTypeSpecName(ASTTrailingTypeSpec):
    def __init__(self, prefix: str, nestedName: ASTNestedName) -> None:
        self.prefix = prefix
        self.nestedName = nestedName

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTrailingTypeSpecName):
            return NotImplemented
        return (
            self.prefix == other.prefix
            and self.nestedName == other.nestedName
        )

    def __hash__(self) -> int:
        return hash((self.prefix, self.nestedName))

    @property
    def name(self) -> ASTNestedName:
        return self.nestedName

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.prefix:
            res.append(self.prefix)
            res.append(' ')
        res.append(transform(self.nestedName))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.prefix:
            signode += addnodes.desc_sig_keyword(self.prefix, self.prefix)
            signode += addnodes.desc_sig_space()
        self.nestedName.describe_signature(signode, mode, env, symbol=symbol)


class ASTFunctionParameter(ASTBase):
    def __init__(self, arg: ASTTypeWithInit | None, ellipsis: bool = False) -> None:
        self.arg = arg
        self.ellipsis = ellipsis

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTFunctionParameter):
            return NotImplemented
        return self.arg == other.arg and self.ellipsis == other.ellipsis

    def __hash__(self) -> int:
        return hash((self.arg, self.ellipsis))

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        # the anchor will be our parent
        return symbol.parent.declaration.get_id(version, prefixed=False)

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.ellipsis:
            return '...'
        else:
            return transform(self.arg)

    def describe_signature(self, signode: Any, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.ellipsis:
            signode += addnodes.desc_sig_punctuation('...', '...')
        else:
            self.arg.describe_signature(signode, mode, env, symbol=symbol)


class ASTParameters(ASTBase):
    def __init__(self, args: list[ASTFunctionParameter], attrs: ASTAttributeList) -> None:
        self.args = args
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTParameters):
            return NotImplemented
        return self.args == other.args and self.attrs == other.attrs

    def __hash__(self) -> int:
        return hash((self.args, self.attrs))

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.args

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append('(')
        first = True
        for a in self.args:
            if not first:
                res.append(', ')
            first = False
            res.append(str(a))
        res.append(')')
        if len(self.attrs) != 0:
            res.append(' ')
            res.append(transform(self.attrs))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        multi_line_parameter_list = False
        test_node: Element = signode
        while test_node.parent:
            if not isinstance(test_node, addnodes.desc_signature):
                test_node = test_node.parent
                continue
            multi_line_parameter_list = test_node.get('multi_line_parameter_list', False)
            break

        # only use the desc_parameterlist for the outer list, not for inner lists
        if mode == 'lastIsName':
            paramlist = addnodes.desc_parameterlist()
            paramlist['multi_line_parameter_list'] = multi_line_parameter_list
            for arg in self.args:
                param = addnodes.desc_parameter('', '', noemph=True)
                arg.describe_signature(param, 'param', env, symbol=symbol)
                paramlist += param
            signode += paramlist
        else:
            signode += addnodes.desc_sig_punctuation('(', '(')
            first = True
            for arg in self.args:
                if not first:
                    signode += addnodes.desc_sig_punctuation(',', ',')
                    signode += addnodes.desc_sig_space()
                first = False
                arg.describe_signature(signode, 'markType', env, symbol=symbol)
            signode += addnodes.desc_sig_punctuation(')', ')')

        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
            self.attrs.describe_signature(signode)


class ASTDeclSpecsSimple(ASTBaseBase):
    def __init__(self, storage: str, threadLocal: str, inline: bool,
                 restrict: bool, volatile: bool, const: bool, attrs: ASTAttributeList) -> None:
        self.storage = storage
        self.threadLocal = threadLocal
        self.inline = inline
        self.restrict = restrict
        self.volatile = volatile
        self.const = const
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclSpecsSimple):
            return NotImplemented
        return (
            self.storage == other.storage
            and self.threadLocal == other.threadLocal
            and self.inline == other.inline
            and self.restrict == other.restrict
            and self.volatile == other.volatile
            and self.const == other.const
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((
            self.storage,
            self.threadLocal,
            self.inline,
            self.restrict,
            self.volatile,
            self.const,
            self.attrs,
        ))

    def mergeWith(self, other: ASTDeclSpecsSimple) -> ASTDeclSpecsSimple:
        if not other:
            return self
        return ASTDeclSpecsSimple(self.storage or other.storage,
                                  self.threadLocal or other.threadLocal,
                                  self.inline or other.inline,
                                  self.volatile or other.volatile,
                                  self.const or other.const,
                                  self.restrict or other.restrict,
                                  self.attrs + other.attrs)

    def _stringify(self, transform: StringifyTransform) -> str:
        res: list[str] = []
        if len(self.attrs) != 0:
            res.append(transform(self.attrs))
        if self.storage:
            res.append(self.storage)
        if self.threadLocal:
            res.append(self.threadLocal)
        if self.inline:
            res.append('inline')
        if self.restrict:
            res.append('restrict')
        if self.volatile:
            res.append('volatile')
        if self.const:
            res.append('const')
        return ' '.join(res)

    def describe_signature(self, modifiers: list[Node]) -> None:
        def _add(modifiers: list[Node], text: str) -> None:
            if len(modifiers) != 0:
                modifiers.append(addnodes.desc_sig_space())
            modifiers.append(addnodes.desc_sig_keyword(text, text))

        if len(modifiers) != 0 and len(self.attrs) != 0:
            modifiers.append(addnodes.desc_sig_space())
        tempNode = nodes.TextElement()
        self.attrs.describe_signature(tempNode)
        modifiers.extend(tempNode.children)
        if self.storage:
            _add(modifiers, self.storage)
        if self.threadLocal:
            _add(modifiers, self.threadLocal)
        if self.inline:
            _add(modifiers, 'inline')
        if self.restrict:
            _add(modifiers, 'restrict')
        if self.volatile:
            _add(modifiers, 'volatile')
        if self.const:
            _add(modifiers, 'const')


class ASTDeclSpecs(ASTBase):
    def __init__(self, outer: str,
                 leftSpecs: ASTDeclSpecsSimple,
                 rightSpecs: ASTDeclSpecsSimple,
                 trailing: ASTTrailingTypeSpec) -> None:
        # leftSpecs and rightSpecs are used for output
        # allSpecs are used for id generation TODO: remove?
        self.outer = outer
        self.leftSpecs = leftSpecs
        self.rightSpecs = rightSpecs
        self.allSpecs = self.leftSpecs.mergeWith(self.rightSpecs)
        self.trailingTypeSpec = trailing

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclSpecs):
            return NotImplemented
        return (
            self.outer == other.outer
            and self.leftSpecs == other.leftSpecs
            and self.rightSpecs == other.rightSpecs
            and self.trailingTypeSpec == other.trailingTypeSpec
        )

    def __hash__(self) -> int:
        return hash((
            self.outer,
            self.leftSpecs,
            self.rightSpecs,
            self.trailingTypeSpec,
        ))

    def _stringify(self, transform: StringifyTransform) -> str:
        res: list[str] = []
        l = transform(self.leftSpecs)
        if len(l) > 0:
            res.append(l)
        if self.trailingTypeSpec:
            if len(res) > 0:
                res.append(" ")
            res.append(transform(self.trailingTypeSpec))
            r = str(self.rightSpecs)
            if len(r) > 0:
                if len(res) > 0:
                    res.append(" ")
                res.append(r)
        return "".join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        modifiers: list[Node] = []

        self.leftSpecs.describe_signature(modifiers)

        for m in modifiers:
            signode += m
        if self.trailingTypeSpec:
            if len(modifiers) > 0:
                signode += addnodes.desc_sig_space()
            self.trailingTypeSpec.describe_signature(signode, mode, env,
                                                     symbol=symbol)
            modifiers = []
            self.rightSpecs.describe_signature(modifiers)
            if len(modifiers) > 0:
                signode += addnodes.desc_sig_space()
            for m in modifiers:
                signode += m


# Declarator
################################################################################

class ASTArray(ASTBase):
    def __init__(self, static: bool, const: bool, volatile: bool, restrict: bool,
                 vla: bool, size: ASTExpression) -> None:
        self.static = static
        self.const = const
        self.volatile = volatile
        self.restrict = restrict
        self.vla = vla
        self.size = size
        if vla:
            assert size is None
        if size is not None:
            assert not vla

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTArray):
            return NotImplemented
        return (
            self.static == other.static
            and self.const == other.const
            and self.volatile == other.volatile
            and self.restrict == other.restrict
            and self.vla == other.vla
            and self.size == other.size
        )

    def __hash__(self) -> int:
        return hash((
            self.static,
            self.const,
            self.volatile,
            self.restrict,
            self.vla,
            self.size,
        ))

    def _stringify(self, transform: StringifyTransform) -> str:
        el = []
        if self.static:
            el.append('static')
        if self.restrict:
            el.append('restrict')
        if self.volatile:
            el.append('volatile')
        if self.const:
            el.append('const')
        if self.vla:
            return '[' + ' '.join(el) + '*]'
        elif self.size:
            el.append(transform(self.size))
        return '[' + ' '.join(el) + ']'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('[', '[')
        addSpace = False

        def _add(signode: TextElement, text: str) -> bool:
            if addSpace:
                signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_keyword(text, text)
            return True

        if self.static:
            addSpace = _add(signode, 'static')
        if self.restrict:
            addSpace = _add(signode, 'restrict')
        if self.volatile:
            addSpace = _add(signode, 'volatile')
        if self.const:
            addSpace = _add(signode, 'const')
        if self.vla:
            signode += addnodes.desc_sig_punctuation('*', '*')
        elif self.size:
            if addSpace:
                signode += addnodes.desc_sig_space()
            self.size.describe_signature(signode, 'markType', env, symbol)
        signode += addnodes.desc_sig_punctuation(']', ']')


class ASTDeclarator(ASTBase):
    @property
    def name(self) -> ASTNestedName:
        raise NotImplementedError(repr(self))

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        raise NotImplementedError(repr(self))

    def require_space_after_declSpecs(self) -> bool:
        raise NotImplementedError(repr(self))


class ASTDeclaratorNameParam(ASTDeclarator):
    def __init__(self, declId: ASTNestedName,
                 arrayOps: list[ASTArray], param: ASTParameters) -> None:
        self.declId = declId
        self.arrayOps = arrayOps
        self.param = param

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorNameParam):
            return NotImplemented
        return (
            self.declId == other.declId
            and self.arrayOps == other.arrayOps
            and self.param == other.param
        )

    def __hash__(self) -> int:
        return hash((self.declId, self.arrayOps, self.param))

    @property
    def name(self) -> ASTNestedName:
        return self.declId

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.param.function_params

    # ------------------------------------------------------------------------

    def require_space_after_declSpecs(self) -> bool:
        return self.declId is not None

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.declId:
            res.append(transform(self.declId))
        res.extend(transform(op) for op in self.arrayOps)
        if self.param:
            res.append(transform(self.param))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.declId:
            self.declId.describe_signature(signode, mode, env, symbol)
        for op in self.arrayOps:
            op.describe_signature(signode, mode, env, symbol)
        if self.param:
            self.param.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorNameBitField(ASTDeclarator):
    def __init__(self, declId: ASTNestedName, size: ASTExpression) -> None:
        self.declId = declId
        self.size = size

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorNameBitField):
            return NotImplemented
        return self.declId == other.declId and self.size == other.size

    def __hash__(self) -> int:
        return hash((self.declId, self.size))

    @property
    def name(self) -> ASTNestedName:
        return self.declId

    # ------------------------------------------------------------------------

    def require_space_after_declSpecs(self) -> bool:
        return self.declId is not None

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.declId:
            res.append(transform(self.declId))
        res.append(" : ")
        res.append(transform(self.size))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.declId:
            self.declId.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_space()
        signode += addnodes.desc_sig_punctuation(':', ':')
        signode += addnodes.desc_sig_space()
        self.size.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorPtr(ASTDeclarator):
    def __init__(self, next: ASTDeclarator, restrict: bool, volatile: bool, const: bool,
                 attrs: ASTAttributeList) -> None:
        assert next
        self.next = next
        self.restrict = restrict
        self.volatile = volatile
        self.const = const
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorPtr):
            return NotImplemented
        return (
            self.next == other.next
            and self.restrict == other.restrict
            and self.volatile == other.volatile
            and self.const == other.const
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((self.next, self.restrict, self.volatile, self.const, self.attrs))

    @property
    def name(self) -> ASTNestedName:
        return self.next.name

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.next.function_params

    def require_space_after_declSpecs(self) -> bool:
        return self.const or self.volatile or self.restrict or \
            len(self.attrs) > 0 or \
            self.next.require_space_after_declSpecs()

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['*']
        res.append(transform(self.attrs))
        if len(self.attrs) != 0 and (self.restrict or self.volatile or self.const):
            res.append(' ')
        if self.restrict:
            res.append('restrict')
        if self.volatile:
            if self.restrict:
                res.append(' ')
            res.append('volatile')
        if self.const:
            if self.restrict or self.volatile:
                res.append(' ')
            res.append('const')
        if self.const or self.volatile or self.restrict or len(self.attrs) > 0:
            if self.next.require_space_after_declSpecs():
                res.append(' ')
        res.append(transform(self.next))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('*', '*')
        self.attrs.describe_signature(signode)
        if len(self.attrs) != 0 and (self.restrict or self.volatile or self.const):
            signode += addnodes.desc_sig_space()

        def _add_anno(signode: TextElement, text: str) -> None:
            signode += addnodes.desc_sig_keyword(text, text)

        if self.restrict:
            _add_anno(signode, 'restrict')
        if self.volatile:
            if self.restrict:
                signode += addnodes.desc_sig_space()
            _add_anno(signode, 'volatile')
        if self.const:
            if self.restrict or self.volatile:
                signode += addnodes.desc_sig_space()
            _add_anno(signode, 'const')
        if self.const or self.volatile or self.restrict or len(self.attrs) > 0:
            if self.next.require_space_after_declSpecs():
                signode += addnodes.desc_sig_space()
        self.next.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorParen(ASTDeclarator):
    def __init__(self, inner: ASTDeclarator, next: ASTDeclarator) -> None:
        assert inner
        assert next
        self.inner = inner
        self.next = next
        # TODO: we assume the name and params are in inner

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorParen):
            return NotImplemented
        return self.inner == other.inner and self.next == other.next

    def __hash__(self) -> int:
        return hash((self.inner, self.next))

    @property
    def name(self) -> ASTNestedName:
        return self.inner.name

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.inner.function_params

    def require_space_after_declSpecs(self) -> bool:
        return True

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['(']
        res.append(transform(self.inner))
        res.append(')')
        res.append(transform(self.next))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.inner.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')
        self.next.describe_signature(signode, "noneIsName", env, symbol)


# Initializer
################################################################################

class ASTParenExprList(ASTBaseParenExprList):
    def __init__(self, exprs: list[ASTExpression]) -> None:
        self.exprs = exprs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTParenExprList):
            return NotImplemented
        return self.exprs == other.exprs

    def __hash__(self) -> int:
        return hash(self.exprs)

    def _stringify(self, transform: StringifyTransform) -> str:
        exprs = [transform(e) for e in self.exprs]
        return '(%s)' % ', '.join(exprs)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('(', '(')
        first = True
        for e in self.exprs:
            if not first:
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            else:
                first = False
            e.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTBracedInitList(ASTBase):
    def __init__(self, exprs: list[ASTExpression], trailingComma: bool) -> None:
        self.exprs = exprs
        self.trailingComma = trailingComma

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTBracedInitList):
            return NotImplemented
        return self.exprs == other.exprs and self.trailingComma == other.trailingComma

    def __hash__(self) -> int:
        return hash((self.exprs, self.trailingComma))

    def _stringify(self, transform: StringifyTransform) -> str:
        exprs = ', '.join(transform(e) for e in self.exprs)
        trailingComma = ',' if self.trailingComma else ''
        return f'{{{exprs}{trailingComma}}}'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('{', '{')
        first = True
        for e in self.exprs:
            if not first:
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            else:
                first = False
            e.describe_signature(signode, mode, env, symbol)
        if self.trailingComma:
            signode += addnodes.desc_sig_punctuation(',', ',')
        signode += addnodes.desc_sig_punctuation('}', '}')


class ASTInitializer(ASTBase):
    def __init__(self, value: ASTBracedInitList | ASTExpression,
                 hasAssign: bool = True) -> None:
        self.value = value
        self.hasAssign = hasAssign

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTInitializer):
            return NotImplemented
        return self.value == other.value and self.hasAssign == other.hasAssign

    def __hash__(self) -> int:
        return hash((self.value, self.hasAssign))

    def _stringify(self, transform: StringifyTransform) -> str:
        val = transform(self.value)
        if self.hasAssign:
            return ' = ' + val
        else:
            return val

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.hasAssign:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('=', '=')
            signode += addnodes.desc_sig_space()
        self.value.describe_signature(signode, 'markType', env, symbol)


class ASTType(ASTBase):
    def __init__(self, declSpecs: ASTDeclSpecs, decl: ASTDeclarator) -> None:
        assert declSpecs
        assert decl
        self.declSpecs = declSpecs
        self.decl = decl

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTType):
            return NotImplemented
        return self.declSpecs == other.declSpecs and self.decl == other.decl

    def __hash__(self) -> int:
        return hash((self.declSpecs, self.decl))

    @property
    def name(self) -> ASTNestedName:
        return self.decl.name

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.decl.function_params

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        declSpecs = transform(self.declSpecs)
        res.append(declSpecs)
        if self.decl.require_space_after_declSpecs() and len(declSpecs) > 0:
            res.append(' ')
        res.append(transform(self.decl))
        return ''.join(res)

    def get_type_declaration_prefix(self) -> str:
        if self.declSpecs.trailingTypeSpec:
            return 'typedef'
        else:
            return 'type'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.declSpecs.describe_signature(signode, 'markType', env, symbol)
        if (self.decl.require_space_after_declSpecs() and
                len(str(self.declSpecs)) > 0):
            signode += addnodes.desc_sig_space()
        # for parameters that don't really declare new names we get 'markType',
        # this should not be propagated, but be 'noneIsName'.
        if mode == 'markType':
            mode = 'noneIsName'
        self.decl.describe_signature(signode, mode, env, symbol)


class ASTTypeWithInit(ASTBase):
    def __init__(self, type: ASTType, init: ASTInitializer) -> None:
        self.type = type
        self.init = init

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTypeWithInit):
            return NotImplemented
        return self.type == other.type and self.init == other.init

    def __hash__(self) -> int:
        return hash((self.type, self.init))

    @property
    def name(self) -> ASTNestedName:
        return self.type.name

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return self.type.get_id(version, objectType, symbol)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.type))
        if self.init:
            res.append(transform(self.init))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.type.describe_signature(signode, mode, env, symbol)
        if self.init:
            self.init.describe_signature(signode, mode, env, symbol)


class ASTMacroParameter(ASTBase):
    def __init__(self, arg: ASTNestedName | None, ellipsis: bool = False,
                 variadic: bool = False) -> None:
        self.arg = arg
        self.ellipsis = ellipsis
        self.variadic = variadic

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTMacroParameter):
            return NotImplemented
        return (
            self.arg == other.arg
            and self.ellipsis == other.ellipsis
            and self.variadic == other.variadic
        )

    def __hash__(self) -> int:
        return hash((self.arg, self.ellipsis, self.variadic))

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.ellipsis:
            return '...'
        elif self.variadic:
            return transform(self.arg) + '...'
        else:
            return transform(self.arg)

    def describe_signature(self, signode: Any, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.ellipsis:
            signode += addnodes.desc_sig_punctuation('...', '...')
        elif self.variadic:
            name = str(self)
            signode += addnodes.desc_sig_name(name, name)
        else:
            self.arg.describe_signature(signode, mode, env, symbol=symbol)


class ASTMacro(ASTBase):
    def __init__(self, ident: ASTNestedName, args: list[ASTMacroParameter] | None) -> None:
        self.ident = ident
        self.args = args

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTMacro):
            return NotImplemented
        return self.ident == other.ident and self.args == other.args

    def __hash__(self) -> int:
        return hash((self.ident, self.args))

    @property
    def name(self) -> ASTNestedName:
        return self.ident

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.ident))
        if self.args is not None:
            res.append('(')
            first = True
            for arg in self.args:
                if not first:
                    res.append(', ')
                first = False
                res.append(transform(arg))
            res.append(')')
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.ident.describe_signature(signode, mode, env, symbol)
        if self.args is None:
            return
        paramlist = addnodes.desc_parameterlist()
        for arg in self.args:
            param = addnodes.desc_parameter('', '', noemph=True)
            arg.describe_signature(param, 'param', env, symbol=symbol)
            paramlist += param
        signode += paramlist


class ASTStruct(ASTBase):
    def __init__(self, name: ASTNestedName) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTStruct):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.name)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.name.describe_signature(signode, mode, env, symbol=symbol)


class ASTUnion(ASTBase):
    def __init__(self, name: ASTNestedName) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTUnion):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.name)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.name.describe_signature(signode, mode, env, symbol=symbol)


class ASTEnum(ASTBase):
    def __init__(self, name: ASTNestedName) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTEnum):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.name)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.name.describe_signature(signode, mode, env, symbol=symbol)


class ASTEnumerator(ASTBase):
    def __init__(self, name: ASTNestedName, init: ASTInitializer | None,
                 attrs: ASTAttributeList) -> None:
        self.name = name
        self.init = init
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTEnumerator):
            return NotImplemented
        return (
            self.name == other.name
            and self.init == other.init
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((self.name, self.init, self.attrs))

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.name))
        if len(self.attrs) != 0:
            res.append(' ')
            res.append(transform(self.attrs))
        if self.init:
            res.append(transform(self.init))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.name.describe_signature(signode, mode, env, symbol)
        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
            self.attrs.describe_signature(signode)
        if self.init:
            self.init.describe_signature(signode, 'markType', env, symbol)


class ASTDeclaration(ASTBaseBase):
    def __init__(self, objectType: str, directiveType: str | None,
                 declaration: DeclarationType | ASTFunctionParameter,
                 semicolon: bool = False) -> None:
        self.objectType = objectType
        self.directiveType = directiveType
        self.declaration = declaration
        self.semicolon = semicolon

        self.symbol: Symbol | None = None
        # set by CObject._add_enumerator_to_parent
        self.enumeratorScopedSymbol: Symbol | None = None

        # the cache assumes that by the time get_newest_id is called, no
        # further changes will be made to this object
        self._newest_id_cache: str | None = None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaration):
            return NotImplemented
        return (
            self.objectType == other.objectType
            and self.directiveType == other.directiveType
            and self.declaration == other.declaration
            and self.semicolon == other.semicolon
            and self.symbol == other.symbol
            and self.enumeratorScopedSymbol == other.enumeratorScopedSymbol
        )

    def clone(self) -> ASTDeclaration:
        return ASTDeclaration(self.objectType, self.directiveType,
                              self.declaration.clone(), self.semicolon)

    @property
    def name(self) -> ASTNestedName:
        decl = cast(DeclarationType, self.declaration)
        return decl.name

    @property
    def function_params(self) -> list[ASTFunctionParameter] | None:
        if self.objectType != 'function':
            return None
        decl = cast(ASTType, self.declaration)
        return decl.function_params

    def get_id(self, version: int, prefixed: bool = True) -> str:
        if self.objectType == 'enumerator' and self.enumeratorScopedSymbol:
            return self.enumeratorScopedSymbol.declaration.get_id(version, prefixed)
        id_ = self.declaration.get_id(version, self.objectType, self.symbol)
        if prefixed:
            return _id_prefix[version] + id_
        else:
            return id_

    def get_newest_id(self) -> str:
        if self._newest_id_cache is None:
            self._newest_id_cache = self.get_id(_max_id, True)
        return self._newest_id_cache

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.declaration)
        if self.semicolon:
            res += ';'
        return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, options: dict[str, bool]) -> None:
        verify_description_mode(mode)
        assert self.symbol
        # The caller of the domain added a desc_signature node.
        # Always enable multiline:
        signode['is_multiline'] = True
        # Put each line in a desc_signature_line node.
        mainDeclNode = addnodes.desc_signature_line()
        mainDeclNode.sphinx_line_type = 'declarator'
        mainDeclNode['add_permalink'] = not self.symbol.isRedeclaration
        signode += mainDeclNode

        if self.objectType in {'member', 'function', 'macro'}:
            pass
        elif self.objectType == 'struct':
            mainDeclNode += addnodes.desc_sig_keyword('struct', 'struct')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'union':
            mainDeclNode += addnodes.desc_sig_keyword('union', 'union')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'enum':
            mainDeclNode += addnodes.desc_sig_keyword('enum', 'enum')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'enumerator':
            mainDeclNode += addnodes.desc_sig_keyword('enumerator', 'enumerator')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'type':
            decl = cast(ASTType, self.declaration)
            prefix = decl.get_type_declaration_prefix()
            mainDeclNode += addnodes.desc_sig_keyword(prefix, prefix)
            mainDeclNode += addnodes.desc_sig_space()
        else:
            raise AssertionError
        self.declaration.describe_signature(mainDeclNode, mode, env, self.symbol)
        if self.semicolon:
            mainDeclNode += addnodes.desc_sig_punctuation(';', ';')
