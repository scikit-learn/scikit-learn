from __future__ import annotations

import sys
import warnings
from typing import TYPE_CHECKING, Any, ClassVar, Literal

from docutils import nodes

from sphinx import addnodes
from sphinx.domains.cpp._ids import (
    _id_char_from_prefix,
    _id_explicit_cast,
    _id_fundamental_v1,
    _id_fundamental_v2,
    _id_operator_unary_v2,
    _id_operator_v1,
    _id_operator_v2,
    _id_prefix,
    _id_shorthands_v1,
    _max_id,
)
from sphinx.util.cfamily import (
    ASTAttributeList,
    ASTBaseBase,
    ASTBaseParenExprList,
    NoOldIdError,
    UnsupportedMultiCharacterCharLiteral,
    verify_description_mode,
)

if TYPE_CHECKING:
    from docutils.nodes import Element, TextElement

    from sphinx.addnodes import desc_signature
    from sphinx.domains.cpp._symbol import Symbol
    from sphinx.environment import BuildEnvironment
    from sphinx.util.cfamily import StringifyTransform


class ASTBase(ASTBaseBase):
    pass


# Names
################################################################################

class ASTIdentifier(ASTBase):
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

    def __hash__(self) -> int:
        return hash(self.name)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.name)

    def is_anon(self) -> bool:
        return self.is_anonymous

    def get_id(self, version: int) -> str:
        if self.is_anonymous and version < 3:
            raise NoOldIdError
        if version == 1:
            if self.name == 'size_t':
                return 's'
            else:
                return self.name
        if self.name == "std":
            return 'St'
        elif self.name[0] == "~":
            # a destructor, just use an arbitrary version of dtors
            return 'D0'
        else:
            if self.is_anonymous:
                return 'Ut%d_%s' % (len(self.name) - 1, self.name[1:])
            else:
                return str(len(self.name)) + self.name

    # and this is where we finally make a difference between __str__ and the display string

    def __str__(self) -> str:
        return self.name

    def get_display_string(self) -> str:
        return "[anonymous]" if self.is_anonymous else self.name

    def describe_signature(self, signode: TextElement, mode: str, env: BuildEnvironment,
                           prefix: str, templateArgs: str, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.is_anonymous:
            node = addnodes.desc_sig_name(text="[anonymous]")
        else:
            node = addnodes.desc_sig_name(self.name, self.name)
        if mode == 'markType':
            targetText = prefix + self.name + templateArgs
            pnode = addnodes.pending_xref('', refdomain='cpp',
                                          reftype='identifier',
                                          reftarget=targetText, modname=None,
                                          classname=None)
            pnode['cpp:parent_key'] = symbol.get_lookup_key()
            pnode += node
            signode += pnode
        elif mode == 'lastIsName':
            nameNode = addnodes.desc_name()
            nameNode += node
            signode += nameNode
        elif mode == 'noneIsName':
            signode += node
        elif mode == 'param':
            node['classes'].append('sig-param')
            signode += node
        elif mode == 'udl':
            # the target is 'operator""id' instead of just 'id'
            assert len(prefix) == 0
            assert len(templateArgs) == 0
            assert not self.is_anonymous
            targetText = 'operator""' + self.name
            pnode = addnodes.pending_xref('', refdomain='cpp',
                                          reftype='identifier',
                                          reftarget=targetText, modname=None,
                                          classname=None)
            pnode['cpp:parent_key'] = symbol.get_lookup_key()
            pnode += node
            signode += pnode
        else:
            raise Exception('Unknown description mode: %s' % mode)

    @property
    def identifier(self) -> str:
        warnings.warn(
            '`ASTIdentifier.identifier` is deprecated, use `ASTIdentifier.name` instead',
            DeprecationWarning, stacklevel=2,
        )
        return self.name


class ASTNestedNameElement(ASTBase):
    def __init__(self, identOrOp: ASTIdentifier | ASTOperator,
                 templateArgs: ASTTemplateArgs | None) -> None:
        self.identOrOp = identOrOp
        self.templateArgs = templateArgs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNestedNameElement):
            return NotImplemented
        return self.identOrOp == other.identOrOp and self.templateArgs == other.templateArgs

    def __hash__(self) -> int:
        return hash((self.identOrOp, self.templateArgs))

    def is_operator(self) -> bool:
        return False

    def get_id(self, version: int) -> str:
        res = self.identOrOp.get_id(version)
        if self.templateArgs:
            res += self.templateArgs.get_id(version)
        return res

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.identOrOp)
        if self.templateArgs:
            res += transform(self.templateArgs)
        return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, prefix: str, symbol: Symbol) -> None:
        tArgs = str(self.templateArgs) if self.templateArgs is not None else ''
        self.identOrOp.describe_signature(signode, mode, env, prefix, tArgs, symbol)
        if self.templateArgs is not None:
            self.templateArgs.describe_signature(signode, 'markType', env, symbol)


class ASTNestedName(ASTBase):
    def __init__(self, names: list[ASTNestedNameElement],
                 templates: list[bool], rooted: bool) -> None:
        assert len(names) > 0
        self.names = names
        self.templates = templates
        assert len(self.names) == len(self.templates)
        self.rooted = rooted

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNestedName):
            return NotImplemented
        return (
            self.names == other.names
            and self.templates == other.templates
            and self.rooted == other.rooted
        )

    def __hash__(self) -> int:
        return hash((self.names, self.templates, self.rooted))

    @property
    def name(self) -> ASTNestedName:
        return self

    def num_templates(self) -> int:
        count = 0
        for n in self.names:
            if n.is_operator():
                continue
            if n.templateArgs:
                count += 1
        return count

    def get_id(self, version: int, modifiers: str = '') -> str:
        if version == 1:
            tt = str(self)
            if tt in _id_shorthands_v1:
                return _id_shorthands_v1[tt]
            else:
                return '::'.join(n.get_id(version) for n in self.names)

        res = []
        if len(self.names) > 1 or len(modifiers) > 0:
            res.append('N')
        res.append(modifiers)
        res.extend(n.get_id(version) for n in self.names)
        if len(self.names) > 1 or len(modifiers) > 0:
            res.append('E')
        return ''.join(res)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.rooted:
            res.append('')
        for i in range(len(self.names)):
            n = self.names[i]
            if self.templates[i]:
                res.append("template " + transform(n))
            else:
                res.append(transform(n))
        return '::'.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        # just print the name part, with template args, not template params
        if mode == 'noneIsName':
            if self.rooted:
                unreachable = "Can this happen?"
                raise AssertionError(unreachable)  # TODO
                signode += nodes.Text('::')
            for i in range(len(self.names)):
                if i != 0:
                    unreachable = "Can this happen?"
                    raise AssertionError(unreachable)  # TODO
                    signode += nodes.Text('::blah')
                n = self.names[i]
                if self.templates[i]:
                    unreachable = "Can this happen?"
                    raise AssertionError(unreachable)  # TODO
                    signode += nodes.Text("template")
                    signode += nodes.Text(" ")
                n.describe_signature(signode, mode, env, '', symbol)
        elif mode == 'param':
            assert not self.rooted, str(self)
            assert len(self.names) == 1
            assert not self.templates[0]
            self.names[0].describe_signature(signode, 'param', env, '', symbol)
        elif mode in ('markType', 'lastIsName', 'markName'):
            # Each element should be a pending xref targeting the complete
            # prefix. however, only the identifier part should be a link, such
            # that template args can be a link as well.
            # For 'lastIsName' we should also prepend template parameter lists.
            templateParams: list[Any] = []
            if mode == 'lastIsName':
                assert symbol is not None
                if symbol.declaration.templatePrefix is not None:
                    templateParams = symbol.declaration.templatePrefix.templates
            iTemplateParams = 0
            templateParamsPrefix = ''
            prefix = ''
            first = True
            names = self.names[:-1] if mode == 'lastIsName' else self.names
            # If lastIsName, then wrap all of the prefix in a desc_addname,
            # else append directly to signode.
            # NOTE: Breathe previously relied on the prefix being in the desc_addname node,
            #       so it can remove it in inner declarations.
            dest = signode
            if mode == 'lastIsName':
                dest = addnodes.desc_addname()
            if self.rooted:
                prefix += '::'
                if mode == 'lastIsName' and len(names) == 0:
                    signode += addnodes.desc_sig_punctuation('::', '::')
                else:
                    dest += addnodes.desc_sig_punctuation('::', '::')
            for i in range(len(names)):
                nne = names[i]
                template = self.templates[i]
                if not first:
                    dest += addnodes.desc_sig_punctuation('::', '::')
                    prefix += '::'
                if template:
                    dest += addnodes.desc_sig_keyword('template', 'template')
                    dest += addnodes.desc_sig_space()
                first = False
                txt_nne = str(nne)
                if txt_nne != '':
                    if nne.templateArgs and iTemplateParams < len(templateParams):
                        templateParamsPrefix += str(templateParams[iTemplateParams])
                        iTemplateParams += 1
                    nne.describe_signature(dest, 'markType',
                                           env, templateParamsPrefix + prefix, symbol)
                prefix += txt_nne
            if mode == 'lastIsName':
                if len(self.names) > 1:
                    dest += addnodes.desc_sig_punctuation('::', '::')
                    signode += dest
                if self.templates[-1]:
                    signode += addnodes.desc_sig_keyword('template', 'template')
                    signode += addnodes.desc_sig_space()
                self.names[-1].describe_signature(signode, mode, env, '', symbol)
        else:
            raise Exception('Unknown description mode: %s' % mode)


################################################################################
# Expressions
################################################################################

class ASTExpression(ASTBase):
    def get_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))


# Primary expressions
################################################################################

class ASTLiteral(ASTExpression):
    pass


class ASTPointerLiteral(ASTLiteral):
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ASTPointerLiteral)

    def __hash__(self) -> int:
        return hash('nullptr')

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'nullptr'

    def get_id(self, version: int) -> str:
        return 'LDnE'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('nullptr', 'nullptr')


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

    def get_id(self, version: int) -> str:
        if self.value:
            return 'L1E'
        else:
            return 'L0E'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword(str(self), str(self))


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

    def get_id(self, version: int) -> str:
        # TODO: floats should be mangled by writing the hex of the binary representation
        return "L%sE" % self.data.replace("'", "")

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_literal_number(self.data, self.data)


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

    def get_id(self, version: int) -> str:
        # note: the length is not really correct with escaping
        return "LA%d_KcE" % (len(self.data) - 2)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_literal_string(self.data, self.data)


class ASTCharLiteral(ASTLiteral):
    def __init__(self, prefix: str, data: str) -> None:
        self.prefix = prefix  # may be None when no prefix
        self.data = data
        assert prefix in _id_char_from_prefix
        self.type = _id_char_from_prefix[prefix]
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

    def get_id(self, version: int) -> str:
        # TODO: the ID should be have L E around it
        return self.type + str(self.value)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.prefix is not None:
            signode += addnodes.desc_sig_keyword(self.prefix, self.prefix)
        txt = "'" + self.data + "'"
        signode += addnodes.desc_sig_literal_char(txt, txt)


class ASTUserDefinedLiteral(ASTLiteral):
    def __init__(self, literal: ASTLiteral, ident: ASTIdentifier) -> None:
        self.literal = literal
        self.ident = ident

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTUserDefinedLiteral):
            return NotImplemented
        return self.literal == other.literal and self.ident == other.ident

    def __hash__(self) -> int:
        return hash((self.literal, self.ident))

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.literal) + transform(self.ident)

    def get_id(self, version: int) -> str:
        # mangle as if it was a function call: ident(literal)
        return f'clL_Zli{self.ident.get_id(version)}E{self.literal.get_id(version)}E'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.literal.describe_signature(signode, mode, env, symbol)
        self.ident.describe_signature(signode, "udl", env, "", "", symbol)


################################################################################

class ASTThisLiteral(ASTExpression):
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ASTThisLiteral)

    def __hash__(self) -> int:
        return hash("this")

    def _stringify(self, transform: StringifyTransform) -> str:
        return "this"

    def get_id(self, version: int) -> str:
        return "fpT"

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('this', 'this')


class ASTFoldExpr(ASTExpression):
    def __init__(self, leftExpr: ASTExpression | None,
                 op: str, rightExpr: ASTExpression | None) -> None:
        assert leftExpr is not None or rightExpr is not None
        self.leftExpr = leftExpr
        self.op = op
        self.rightExpr = rightExpr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTFoldExpr):
            return NotImplemented
        return (
            self.leftExpr == other.leftExpr
            and self.op == other.op
            and self.rightExpr == other.rightExpr
        )

    def __hash__(self) -> int:
        return hash((self.leftExpr, self.op, self.rightExpr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['(']
        if self.leftExpr:
            res.append(transform(self.leftExpr))
            res.append(' ')
            res.append(self.op)
            res.append(' ')
        res.append('...')
        if self.rightExpr:
            res.append(' ')
            res.append(self.op)
            res.append(' ')
            res.append(transform(self.rightExpr))
        res.append(')')
        return ''.join(res)

    def get_id(self, version: int) -> str:
        assert version >= 3
        if version == 3:
            return str(self)
        # https://github.com/itanium-cxx-abi/cxx-abi/pull/67
        res = []
        if self.leftExpr is None:  # (... op expr)
            res.append('fl')
        elif self.rightExpr is None:  # (expr op ...)
            res.append('fr')
        else:  # (expr op ... op expr)
            # we don't check where the parameter pack is,
            # we just always call this a binary left fold
            res.append('fL')
        res.append(_id_operator_v2[self.op])
        if self.leftExpr:
            res.append(self.leftExpr.get_id(version))
        if self.rightExpr:
            res.append(self.rightExpr.get_id(version))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('(', '(')
        if self.leftExpr:
            self.leftExpr.describe_signature(signode, mode, env, symbol)
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_operator(self.op, self.op)
            signode += addnodes.desc_sig_space()
        signode += addnodes.desc_sig_punctuation('...', '...')
        if self.rightExpr:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_operator(self.op, self.op)
            signode += addnodes.desc_sig_space()
            self.rightExpr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


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
        return self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


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


# Postfix expressions
################################################################################

class ASTPostfixOp(ASTBase):
    def get_id(self, idPrefix: str, version: int) -> str:
        raise NotImplementedError(repr(self))

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))


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

    def get_id(self, idPrefix: str, version: int) -> str:
        return 'ix' + idPrefix + self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('[', '[')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(']', ']')


class ASTPostfixMember(ASTPostfixOp):
    def __init__(self, name: ASTNestedName) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPostfixMember):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def _stringify(self, transform: StringifyTransform) -> str:
        return '.' + transform(self.name)

    def get_id(self, idPrefix: str, version: int) -> str:
        return 'dt' + idPrefix + self.name.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('.', '.')
        self.name.describe_signature(signode, 'noneIsName', env, symbol)


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

    def get_id(self, idPrefix: str, version: int) -> str:
        return 'pt' + idPrefix + self.name.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('->', '->')
        self.name.describe_signature(signode, 'noneIsName', env, symbol)


class ASTPostfixInc(ASTPostfixOp):
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ASTPostfixInc)

    def __hash__(self) -> int:
        return hash('++')

    def _stringify(self, transform: StringifyTransform) -> str:
        return '++'

    def get_id(self, idPrefix: str, version: int) -> str:
        return 'pp' + idPrefix

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('++', '++')


class ASTPostfixDec(ASTPostfixOp):
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ASTPostfixDec)

    def __hash__(self) -> int:
        return hash('--')

    def _stringify(self, transform: StringifyTransform) -> str:
        return '--'

    def get_id(self, idPrefix: str, version: int) -> str:
        return 'mm' + idPrefix

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_operator('--', '--')


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

    def get_id(self, idPrefix: str, version: int) -> str:
        return ''.join([
            'cl',
            idPrefix,
            *(e.get_id(version) for e in self.lst.exprs),
            'E',
        ])

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.lst.describe_signature(signode, mode, env, symbol)


class ASTPostfixExpr(ASTExpression):
    def __init__(self, prefix: ASTType, postFixes: list[ASTPostfixOp]) -> None:
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

    def get_id(self, version: int) -> str:
        id = self.prefix.get_id(version)
        for p in self.postFixes:
            id = p.get_id(id, version)
        return id

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.prefix.describe_signature(signode, mode, env, symbol)
        for p in self.postFixes:
            p.describe_signature(signode, mode, env, symbol)


class ASTExplicitCast(ASTExpression):
    def __init__(self, cast: str, typ: ASTType, expr: ASTExpression) -> None:
        assert cast in _id_explicit_cast
        self.cast = cast
        self.typ = typ
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTExplicitCast):
            return NotImplemented
        return self.cast == other.cast and self.typ == other.typ and self.expr == other.expr

    def __hash__(self) -> int:
        return hash((self.cast, self.typ, self.expr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = [self.cast]
        res.append('<')
        res.append(transform(self.typ))
        res.append('>(')
        res.append(transform(self.expr))
        res.append(')')
        return ''.join(res)

    def get_id(self, version: int) -> str:
        return (_id_explicit_cast[self.cast] +
                self.typ.get_id(version) +
                self.expr.get_id(version))

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword(self.cast, self.cast)
        signode += addnodes.desc_sig_punctuation('<', '<')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation('>', '>')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTTypeId(ASTExpression):
    def __init__(self, typeOrExpr: ASTType | ASTExpression, isType: bool) -> None:
        self.typeOrExpr = typeOrExpr
        self.isType = isType

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTypeId):
            return NotImplemented
        return self.typeOrExpr == other.typeOrExpr and self.isType == other.isType

    def __hash__(self) -> int:
        return hash((self.typeOrExpr, self.isType))

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'typeid(' + transform(self.typeOrExpr) + ')'

    def get_id(self, version: int) -> str:
        prefix = 'ti' if self.isType else 'te'
        return prefix + self.typeOrExpr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('typeid', 'typeid')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typeOrExpr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


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

    def get_id(self, version: int) -> str:
        return _id_operator_unary_v2[self.op] + self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.op[0] in 'cn':
            signode += addnodes.desc_sig_keyword(self.op, self.op)
            signode += addnodes.desc_sig_space()
        else:
            signode += addnodes.desc_sig_operator(self.op, self.op)
        self.expr.describe_signature(signode, mode, env, symbol)


class ASTSizeofParamPack(ASTExpression):
    def __init__(self, identifier: ASTIdentifier) -> None:
        self.identifier = identifier

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTSizeofParamPack):
            return NotImplemented
        return self.identifier == other.identifier

    def __hash__(self) -> int:
        return hash(self.identifier)

    def _stringify(self, transform: StringifyTransform) -> str:
        return "sizeof...(" + transform(self.identifier) + ")"

    def get_id(self, version: int) -> str:
        return 'sZ' + self.identifier.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('sizeof', 'sizeof')
        signode += addnodes.desc_sig_punctuation('...', '...')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.identifier.describe_signature(signode, 'markType', env,
                                           symbol=symbol, prefix="", templateArgs="")
        signode += addnodes.desc_sig_punctuation(')', ')')


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

    def get_id(self, version: int) -> str:
        return 'st' + self.typ.get_id(version)

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

    def get_id(self, version: int) -> str:
        return 'sz' + self.expr.get_id(version)

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

    def get_id(self, version: int) -> str:
        return 'at' + self.typ.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('alignof', 'alignof')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTNoexceptExpr(ASTExpression):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNoexceptExpr):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'noexcept(' + transform(self.expr) + ')'

    def get_id(self, version: int) -> str:
        return 'nx' + self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('noexcept', 'noexcept')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTNewExpr(ASTExpression):
    def __init__(self, rooted: bool, isNewTypeId: bool, typ: ASTType,
                 initList: ASTParenExprList | ASTBracedInitList) -> None:
        self.rooted = rooted
        self.isNewTypeId = isNewTypeId
        self.typ = typ
        self.initList = initList

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNewExpr):
            return NotImplemented
        return (
            self.rooted == other.rooted
            and self.isNewTypeId == other.isNewTypeId
            and self.typ == other.typ
            and self.initList == other.initList
        )

    def __hash__(self) -> int:
        return hash((self.rooted, self.isNewTypeId, self.typ, self.initList))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.rooted:
            res.append('::')
        res.append('new ')
        # TODO: placement
        if self.isNewTypeId:
            res.append(transform(self.typ))
        else:
            raise AssertionError
        if self.initList is not None:
            res.append(transform(self.initList))
        return ''.join(res)

    def get_id(self, version: int) -> str:
        # the array part will be in the type mangling, so na is not used
        res = ['nw']
        # TODO: placement
        res.append('_')
        res.append(self.typ.get_id(version))
        if self.initList is not None:
            res.append(self.initList.get_id(version))
        else:
            res.append('E')
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.rooted:
            signode += addnodes.desc_sig_punctuation('::', '::')
        signode += addnodes.desc_sig_keyword('new', 'new')
        signode += addnodes.desc_sig_space()
        # TODO: placement
        if self.isNewTypeId:
            self.typ.describe_signature(signode, mode, env, symbol)
        else:
            raise AssertionError
        if self.initList is not None:
            self.initList.describe_signature(signode, mode, env, symbol)


class ASTDeleteExpr(ASTExpression):
    def __init__(self, rooted: bool, array: bool, expr: ASTExpression) -> None:
        self.rooted = rooted
        self.array = array
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeleteExpr):
            return NotImplemented
        return (
            self.rooted == other.rooted
            and self.array == other.array
            and self.expr == other.expr
        )

    def __hash__(self) -> int:
        return hash((self.rooted, self.array, self.expr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.rooted:
            res.append('::')
        res.append('delete ')
        if self.array:
            res.append('[] ')
        res.append(transform(self.expr))
        return ''.join(res)

    def get_id(self, version: int) -> str:
        if self.array:
            id = "da"
        else:
            id = "dl"
        return id + self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.rooted:
            signode += addnodes.desc_sig_punctuation('::', '::')
        signode += addnodes.desc_sig_keyword('delete', 'delete')
        signode += addnodes.desc_sig_space()
        if self.array:
            signode += addnodes.desc_sig_punctuation('[]', '[]')
            signode += addnodes.desc_sig_space()
        self.expr.describe_signature(signode, mode, env, symbol)


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

    def get_id(self, version: int) -> str:
        return 'cv' + self.typ.get_id(version) + self.expr.get_id(version)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.typ.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')
        self.expr.describe_signature(signode, mode, env, symbol)


class ASTBinOpExpr(ASTExpression):
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

    def get_id(self, version: int) -> str:
        assert version >= 2
        res = []
        for i in range(len(self.ops)):
            res.append(_id_operator_v2[self.ops[i]])
            res.append(self.exprs[i].get_id(version))
        res.append(self.exprs[-1].get_id(version))
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


class ASTConditionalExpr(ASTExpression):
    def __init__(self, ifExpr: ASTExpression, thenExpr: ASTExpression,
                 elseExpr: ASTExpression) -> None:
        self.ifExpr = ifExpr
        self.thenExpr = thenExpr
        self.elseExpr = elseExpr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTConditionalExpr):
            return NotImplemented
        return (
            self.ifExpr == other.ifExpr
            and self.thenExpr == other.thenExpr
            and self.elseExpr == other.elseExpr
        )

    def __hash__(self) -> int:
        return hash((self.ifExpr, self.thenExpr, self.elseExpr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.ifExpr))
        res.append(' ? ')
        res.append(transform(self.thenExpr))
        res.append(' : ')
        res.append(transform(self.elseExpr))
        return ''.join(res)

    def get_id(self, version: int) -> str:
        assert version >= 2
        res = []
        res.append(_id_operator_v2['?'])
        res.append(self.ifExpr.get_id(version))
        res.append(self.thenExpr.get_id(version))
        res.append(self.elseExpr.get_id(version))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.ifExpr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_space()
        signode += addnodes.desc_sig_operator('?', '?')
        signode += addnodes.desc_sig_space()
        self.thenExpr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_space()
        signode += addnodes.desc_sig_operator(':', ':')
        signode += addnodes.desc_sig_space()
        self.elseExpr.describe_signature(signode, mode, env, symbol)


class ASTBracedInitList(ASTBase):
    def __init__(self, exprs: list[ASTExpression | ASTBracedInitList],
                 trailingComma: bool) -> None:
        self.exprs = exprs
        self.trailingComma = trailingComma

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTBracedInitList):
            return NotImplemented
        return self.exprs == other.exprs and self.trailingComma == other.trailingComma

    def __hash__(self) -> int:
        return hash((self.exprs, self.trailingComma))

    def get_id(self, version: int) -> str:
        return "il%sE" % ''.join(e.get_id(version) for e in self.exprs)

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


class ASTAssignmentExpr(ASTExpression):
    def __init__(self, leftExpr: ASTExpression, op: str,
                 rightExpr: ASTExpression | ASTBracedInitList) -> None:
        self.leftExpr = leftExpr
        self.op = op
        self.rightExpr = rightExpr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTAssignmentExpr):
            return NotImplemented
        return (
            self.leftExpr == other.leftExpr
            and self.op == other.op
            and self.rightExpr == other.rightExpr
        )

    def __hash__(self) -> int:
        return hash((self.leftExpr, self.op, self.rightExpr))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.leftExpr))
        res.append(' ')
        res.append(self.op)
        res.append(' ')
        res.append(transform(self.rightExpr))
        return ''.join(res)

    def get_id(self, version: int) -> str:
        # we end up generating the ID from left to right, instead of right to left
        res = []
        res.append(_id_operator_v2[self.op])
        res.append(self.leftExpr.get_id(version))
        res.append(self.rightExpr.get_id(version))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.leftExpr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_space()
        if ord(self.op[0]) >= ord('a') and ord(self.op[0]) <= ord('z'):
            signode += addnodes.desc_sig_keyword(self.op, self.op)
        else:
            signode += addnodes.desc_sig_operator(self.op, self.op)
        signode += addnodes.desc_sig_space()
        self.rightExpr.describe_signature(signode, mode, env, symbol)


class ASTCommaExpr(ASTExpression):
    def __init__(self, exprs: list[ASTExpression]) -> None:
        assert len(exprs) > 0
        self.exprs = exprs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTCommaExpr):
            return NotImplemented
        return self.exprs == other.exprs

    def __hash__(self) -> int:
        return hash(self.exprs)

    def _stringify(self, transform: StringifyTransform) -> str:
        return ', '.join(transform(e) for e in self.exprs)

    def get_id(self, version: int) -> str:
        id_ = _id_operator_v2[',']
        res = []
        for i in range(len(self.exprs) - 1):
            res.append(id_)
            res.append(self.exprs[i].get_id(version))
        res.append(self.exprs[-1].get_id(version))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.exprs[0].describe_signature(signode, mode, env, symbol)
        for i in range(1, len(self.exprs)):
            signode += addnodes.desc_sig_punctuation(',', ',')
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

# Things for ASTNestedName
################################################################################

class ASTOperator(ASTBase):
    is_anonymous: ClassVar[Literal[False]] = False

    def __eq__(self, other: object) -> bool:
        raise NotImplementedError(repr(self))

    def __hash__(self) -> int:
        raise NotImplementedError(repr(self))

    def is_anon(self) -> bool:
        return self.is_anonymous

    def is_operator(self) -> bool:
        return True

    def get_id(self, version: int) -> str:
        raise NotImplementedError

    def _describe_identifier(self, signode: TextElement, identnode: TextElement,
                             env: BuildEnvironment, symbol: Symbol) -> None:
        """Render the prefix into signode, and the last part into identnode."""
        raise NotImplementedError

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, prefix: str, templateArgs: str,
                           symbol: Symbol) -> None:
        verify_description_mode(mode)
        if mode == 'lastIsName':
            mainName = addnodes.desc_name()
            self._describe_identifier(mainName, mainName, env, symbol)
            signode += mainName
        elif mode == 'markType':
            targetText = prefix + str(self) + templateArgs
            pnode = addnodes.pending_xref('', refdomain='cpp',
                                          reftype='identifier',
                                          reftarget=targetText, modname=None,
                                          classname=None)
            pnode['cpp:parent_key'] = symbol.get_lookup_key()
            # Render the identifier part, but collapse it into a string
            # and make that the a link to this operator.
            # E.g., if it is 'operator SomeType', then 'SomeType' becomes
            # a link to the operator, not to 'SomeType'.
            container = nodes.literal()
            self._describe_identifier(signode, container, env, symbol)
            txt = container.astext()
            pnode += addnodes.desc_name(txt, txt)
            signode += pnode
        else:
            addName = addnodes.desc_addname()
            self._describe_identifier(addName, addName, env, symbol)
            signode += addName


class ASTOperatorBuildIn(ASTOperator):
    def __init__(self, op: str) -> None:
        self.op = op

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTOperatorBuildIn):
            return NotImplemented
        return self.op == other.op

    def __hash__(self) -> int:
        return hash(self.op)

    def get_id(self, version: int) -> str:
        if version == 1:
            ids = _id_operator_v1
            if self.op not in ids:
                raise NoOldIdError
        else:
            ids = _id_operator_v2
        if self.op not in ids:
            raise Exception('Internal error: Built-in operator "%s" can not '
                            'be mapped to an id.' % self.op)
        return ids[self.op]

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.op in ('new', 'new[]', 'delete', 'delete[]') or self.op[0] in "abcnox":
            return 'operator ' + self.op
        else:
            return 'operator' + self.op

    def _describe_identifier(self, signode: TextElement, identnode: TextElement,
                             env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('operator', 'operator')
        if self.op in ('new', 'new[]', 'delete', 'delete[]') or self.op[0] in "abcnox":
            signode += addnodes.desc_sig_space()
        identnode += addnodes.desc_sig_operator(self.op, self.op)


class ASTOperatorLiteral(ASTOperator):
    def __init__(self, identifier: ASTIdentifier) -> None:
        self.identifier = identifier

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTOperatorLiteral):
            return NotImplemented
        return self.identifier == other.identifier

    def __hash__(self) -> int:
        return hash(self.identifier)

    def get_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError
        return 'li' + self.identifier.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'operator""' + transform(self.identifier)

    def _describe_identifier(self, signode: TextElement, identnode: TextElement,
                             env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('operator', 'operator')
        signode += addnodes.desc_sig_literal_string('""', '""')
        self.identifier.describe_signature(identnode, 'markType', env, '', '', symbol)


class ASTOperatorType(ASTOperator):
    def __init__(self, type: ASTType) -> None:
        self.type = type

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTOperatorType):
            return NotImplemented
        return self.type == other.type

    def __hash__(self) -> int:
        return hash(self.type)

    def get_id(self, version: int) -> str:
        if version == 1:
            return 'castto-%s-operator' % self.type.get_id(version)
        else:
            return 'cv' + self.type.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return f'operator {transform(self.type)}'

    def get_name_no_template(self) -> str:
        return str(self)

    def _describe_identifier(self, signode: TextElement, identnode: TextElement,
                             env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('operator', 'operator')
        signode += addnodes.desc_sig_space()
        self.type.describe_signature(identnode, 'markType', env, symbol)


class ASTTemplateArgConstant(ASTBase):
    def __init__(self, value: ASTExpression) -> None:
        self.value = value

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateArgConstant):
            return NotImplemented
        return self.value == other.value

    def __hash__(self) -> int:
        return hash(self.value)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.value)

    def get_id(self, version: int) -> str:
        if version == 1:
            return str(self).replace(' ', '-')
        if version == 2:
            return 'X' + str(self) + 'E'
        return 'X' + self.value.get_id(version) + 'E'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.value.describe_signature(signode, mode, env, symbol)


class ASTTemplateArgs(ASTBase):
    def __init__(self, args: list[ASTType | ASTTemplateArgConstant],
                 packExpansion: bool) -> None:
        assert args is not None
        self.args = args
        self.packExpansion = packExpansion

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateArgs):
            return NotImplemented
        return self.args == other.args and self.packExpansion == other.packExpansion

    def __hash__(self) -> int:
        return hash((self.args, self.packExpansion))

    def get_id(self, version: int) -> str:
        if version == 1:
            res = []
            res.append(':')
            res.append('.'.join(a.get_id(version) for a in self.args))
            res.append(':')
            return ''.join(res)

        res = []
        res.append('I')
        if len(self.args) > 0:
            for a in self.args[:-1]:
                res.append(a.get_id(version))
            if self.packExpansion:
                res.append('J')
            res.append(self.args[-1].get_id(version))
            if self.packExpansion:
                res.append('E')
        res.append('E')
        return ''.join(res)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ', '.join(transform(a) for a in self.args)
        if self.packExpansion:
            res += '...'
        return '<' + res + '>'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('<', '<')
        first = True
        for a in self.args:
            if not first:
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            first = False
            a.describe_signature(signode, 'markType', env, symbol=symbol)
        if self.packExpansion:
            signode += addnodes.desc_sig_punctuation('...', '...')
        signode += addnodes.desc_sig_punctuation('>', '>')


# Main part of declarations
################################################################################

class ASTTrailingTypeSpec(ASTBase):
    def get_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))


class ASTTrailingTypeSpecFundamental(ASTTrailingTypeSpec):
    def __init__(self, names: list[str], canonNames: list[str]) -> None:
        assert len(names) != 0
        assert len(names) == len(canonNames), (names, canonNames)
        self.names = names
        # the canonical name list is for ID lookup
        self.canonNames = canonNames

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTrailingTypeSpecFundamental):
            return NotImplemented
        return self.names == other.names and self.canonNames == other.canonNames

    def __hash__(self) -> int:
        return hash((self.names, self.canonNames))

    def _stringify(self, transform: StringifyTransform) -> str:
        return ' '.join(self.names)

    def get_id(self, version: int) -> str:
        if version == 1:
            res = []
            for a in self.canonNames:
                if a in _id_fundamental_v1:
                    res.append(_id_fundamental_v1[a])
                else:
                    res.append(a)
            return '-'.join(res)

        txt = ' '.join(self.canonNames)
        if txt not in _id_fundamental_v2:
            raise Exception(
                'Semi-internal error: Fundamental type "%s" can not be mapped '
                'to an ID. Is it a true fundamental type? If not so, the '
                'parser should have rejected it.' % txt)
        return _id_fundamental_v2[txt]

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        first = True
        for n in self.names:
            if not first:
                signode += addnodes.desc_sig_space()
            else:
                first = False
            signode += addnodes.desc_sig_keyword_type(n, n)


class ASTTrailingTypeSpecDecltypeAuto(ASTTrailingTypeSpec):
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ASTTrailingTypeSpecDecltypeAuto)

    def __hash__(self) -> int:
        return hash('decltype(auto)')

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'decltype(auto)'

    def get_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError
        return 'Dc'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('decltype', 'decltype')
        signode += addnodes.desc_sig_punctuation('(', '(')
        signode += addnodes.desc_sig_keyword('auto', 'auto')
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTTrailingTypeSpecDecltype(ASTTrailingTypeSpec):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTrailingTypeSpecDecltype):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'decltype(' + transform(self.expr) + ')'

    def get_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError
        return 'DT' + self.expr.get_id(version) + "E"

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('decltype', 'decltype')
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')


class ASTTrailingTypeSpecName(ASTTrailingTypeSpec):
    def __init__(self, prefix: str, nestedName: ASTNestedName,
                 placeholderType: str | None) -> None:
        self.prefix = prefix
        self.nestedName = nestedName
        self.placeholderType = placeholderType

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTrailingTypeSpecName):
            return NotImplemented
        return (
            self.prefix == other.prefix
            and self.nestedName == other.nestedName
            and self.placeholderType == other.placeholderType
        )

    def __hash__(self) -> int:
        return hash((self.prefix, self.nestedName, self.placeholderType))

    @property
    def name(self) -> ASTNestedName:
        return self.nestedName

    def get_id(self, version: int) -> str:
        return self.nestedName.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.prefix:
            res.append(self.prefix)
            res.append(' ')
        res.append(transform(self.nestedName))
        if self.placeholderType is not None:
            res.append(' ')
            res.append(self.placeholderType)
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.prefix:
            signode += addnodes.desc_sig_keyword(self.prefix, self.prefix)
            signode += addnodes.desc_sig_space()
        self.nestedName.describe_signature(signode, mode, env, symbol=symbol)
        if self.placeholderType is not None:
            signode += addnodes.desc_sig_space()
            if self.placeholderType == 'auto':
                signode += addnodes.desc_sig_keyword('auto', 'auto')
            elif self.placeholderType == 'decltype(auto)':
                signode += addnodes.desc_sig_keyword('decltype', 'decltype')
                signode += addnodes.desc_sig_punctuation('(', '(')
                signode += addnodes.desc_sig_keyword('auto', 'auto')
                signode += addnodes.desc_sig_punctuation(')', ')')
            else:
                raise AssertionError(self.placeholderType)


class ASTFunctionParameter(ASTBase):
    def __init__(self, arg: ASTTypeWithInit | ASTTemplateParamConstrainedTypeWithInit,
                 ellipsis: bool = False) -> None:
        self.arg = arg
        self.ellipsis = ellipsis

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTFunctionParameter):
            return NotImplemented
        return self.arg == other.arg and self.ellipsis == other.ellipsis

    def __hash__(self) -> int:
        return hash((self.arg, self.ellipsis))

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        # this is not part of the normal name mangling in C++
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=False)
        # else, do the usual
        if self.ellipsis:
            return 'z'
        else:
            return self.arg.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.ellipsis:
            return '...'
        else:
            return transform(self.arg)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.ellipsis:
            signode += addnodes.desc_sig_punctuation('...', '...')
        else:
            self.arg.describe_signature(signode, mode, env, symbol=symbol)


class ASTNoexceptSpec(ASTBase):
    def __init__(self, expr: ASTExpression | None) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNoexceptSpec):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.expr:
            return 'noexcept(' + transform(self.expr) + ')'
        return 'noexcept'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('noexcept', 'noexcept')
        if self.expr:
            signode += addnodes.desc_sig_punctuation('(', '(')
            self.expr.describe_signature(signode, 'markType', env, symbol)
            signode += addnodes.desc_sig_punctuation(')', ')')


class ASTParametersQualifiers(ASTBase):
    def __init__(self, args: list[ASTFunctionParameter], volatile: bool, const: bool,
                 refQual: str | None, exceptionSpec: ASTNoexceptSpec,
                 trailingReturn: ASTType,
                 override: bool, final: bool, attrs: ASTAttributeList,
                 initializer: str | None) -> None:
        self.args = args
        self.volatile = volatile
        self.const = const
        self.refQual = refQual
        self.exceptionSpec = exceptionSpec
        self.trailingReturn = trailingReturn
        self.override = override
        self.final = final
        self.attrs = attrs
        self.initializer = initializer

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTParametersQualifiers):
            return NotImplemented
        return (
            self.args == other.args
            and self.volatile == other.volatile
            and self.const == other.const
            and self.refQual == other.refQual
            and self.exceptionSpec == other.exceptionSpec
            and self.trailingReturn == other.trailingReturn
            and self.override == other.override
            and self.final == other.final
            and self.attrs == other.attrs
            and self.initializer == other.initializer
        )

    def __hash__(self) -> int:
        return hash((
            self.args, self.volatile, self.const, self.refQual, self.exceptionSpec,
            self.trailingReturn, self.override, self.final, self.attrs, self.initializer
        ))

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.args

    def get_modifiers_id(self, version: int) -> str:
        res = []
        if self.volatile:
            res.append('V')
        if self.const:
            if version == 1:
                res.append('C')
            else:
                res.append('K')
        if self.refQual == '&&':
            res.append('O')
        elif self.refQual == '&':
            res.append('R')
        return ''.join(res)

    def get_param_id(self, version: int) -> str:
        if version == 1:
            if len(self.args) == 0:
                return ''
            else:
                return '__' + '.'.join(a.get_id(version) for a in self.args)
        if len(self.args) == 0:
            return 'v'
        else:
            return ''.join(a.get_id(version) for a in self.args)

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
        if self.volatile:
            res.append(' volatile')
        if self.const:
            res.append(' const')
        if self.refQual:
            res.append(' ')
            res.append(self.refQual)
        if self.exceptionSpec:
            res.append(' ')
            res.append(transform(self.exceptionSpec))
        if self.trailingReturn:
            res.append(' -> ')
            res.append(transform(self.trailingReturn))
        if self.final:
            res.append(' final')
        if self.override:
            res.append(' override')
        if len(self.attrs) != 0:
            res.append(' ')
            res.append(transform(self.attrs))
        if self.initializer:
            res.append(' = ')
            res.append(self.initializer)
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

        def _add_anno(signode: TextElement, text: str) -> None:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_keyword(text, text)

        if self.volatile:
            _add_anno(signode, 'volatile')
        if self.const:
            _add_anno(signode, 'const')
        if self.refQual:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation(self.refQual, self.refQual)
        if self.exceptionSpec:
            signode += addnodes.desc_sig_space()
            self.exceptionSpec.describe_signature(signode, mode, env, symbol)
        if self.trailingReturn:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_operator('->', '->')
            signode += addnodes.desc_sig_space()
            self.trailingReturn.describe_signature(signode, mode, env, symbol)
        if self.final:
            _add_anno(signode, 'final')
        if self.override:
            _add_anno(signode, 'override')
        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
            self.attrs.describe_signature(signode)
        if self.initializer:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('=', '=')
            signode += addnodes.desc_sig_space()
            assert self.initializer in ('0', 'delete', 'default')
            if self.initializer == '0':
                signode += addnodes.desc_sig_literal_number('0', '0')
            else:
                signode += addnodes.desc_sig_keyword(self.initializer, self.initializer)


class ASTExplicitSpec(ASTBase):
    def __init__(self, expr: ASTExpression | None) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTExplicitSpec):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['explicit']
        if self.expr is not None:
            res.append('(')
            res.append(transform(self.expr))
            res.append(')')
        return ''.join(res)

    def describe_signature(self, signode: TextElement,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('explicit', 'explicit')
        if self.expr is not None:
            signode += addnodes.desc_sig_punctuation('(', '(')
            self.expr.describe_signature(signode, 'markType', env, symbol)
            signode += addnodes.desc_sig_punctuation(')', ')')


class ASTDeclSpecsSimple(ASTBase):
    def __init__(self, storage: str, threadLocal: bool, inline: bool, virtual: bool,
                 explicitSpec: ASTExplicitSpec | None,
                 consteval: bool, constexpr: bool, constinit: bool,
                 volatile: bool, const: bool, friend: bool,
                 attrs: ASTAttributeList) -> None:
        self.storage = storage
        self.threadLocal = threadLocal
        self.inline = inline
        self.virtual = virtual
        self.explicitSpec = explicitSpec
        self.consteval = consteval
        self.constexpr = constexpr
        self.constinit = constinit
        self.volatile = volatile
        self.const = const
        self.friend = friend
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclSpecsSimple):
            return NotImplemented
        return (
            self.storage == other.storage
            and self.threadLocal == other.threadLocal
            and self.inline == other.inline
            and self.virtual == other.virtual
            and self.explicitSpec == other.explicitSpec
            and self.consteval == other.consteval
            and self.constexpr == other.constexpr
            and self.constinit == other.constinit
            and self.volatile == other.volatile
            and self.const == other.const
            and self.friend == other.friend
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((
            self.storage,
            self.threadLocal,
            self.inline,
            self.virtual,
            self.explicitSpec,
            self.consteval,
            self.constexpr,
            self.constinit,
            self.volatile,
            self.const,
            self.friend,
            self.attrs,
        ))

    def mergeWith(self, other: ASTDeclSpecsSimple) -> ASTDeclSpecsSimple:
        if not other:
            return self
        return ASTDeclSpecsSimple(self.storage or other.storage,
                                  self.threadLocal or other.threadLocal,
                                  self.inline or other.inline,
                                  self.virtual or other.virtual,
                                  self.explicitSpec or other.explicitSpec,
                                  self.consteval or other.consteval,
                                  self.constexpr or other.constexpr,
                                  self.constinit or other.constinit,
                                  self.volatile or other.volatile,
                                  self.const or other.const,
                                  self.friend or other.friend,
                                  self.attrs + other.attrs)

    def _stringify(self, transform: StringifyTransform) -> str:
        res: list[str] = []
        if len(self.attrs) != 0:
            res.append(transform(self.attrs))
        if self.storage:
            res.append(self.storage)
        if self.threadLocal:
            res.append('thread_local')
        if self.inline:
            res.append('inline')
        if self.friend:
            res.append('friend')
        if self.virtual:
            res.append('virtual')
        if self.explicitSpec:
            res.append(transform(self.explicitSpec))
        if self.consteval:
            res.append('consteval')
        if self.constexpr:
            res.append('constexpr')
        if self.constinit:
            res.append('constinit')
        if self.volatile:
            res.append('volatile')
        if self.const:
            res.append('const')
        return ' '.join(res)

    def describe_signature(self, signode: TextElement,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.attrs.describe_signature(signode)
        addSpace = len(self.attrs) != 0

        def _add(signode: TextElement, text: str) -> bool:
            if addSpace:
                signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_keyword(text, text)
            return True

        if self.storage:
            addSpace = _add(signode, self.storage)
        if self.threadLocal:
            addSpace = _add(signode, 'thread_local')
        if self.inline:
            addSpace = _add(signode, 'inline')
        if self.friend:
            addSpace = _add(signode, 'friend')
        if self.virtual:
            addSpace = _add(signode, 'virtual')
        if self.explicitSpec:
            if addSpace:
                signode += addnodes.desc_sig_space()
            self.explicitSpec.describe_signature(signode, env, symbol)
            addSpace = True
        if self.consteval:
            addSpace = _add(signode, 'consteval')
        if self.constexpr:
            addSpace = _add(signode, 'constexpr')
        if self.constinit:
            addSpace = _add(signode, 'constinit')
        if self.volatile:
            addSpace = _add(signode, 'volatile')
        if self.const:
            addSpace = _add(signode, 'const')


class ASTDeclSpecs(ASTBase):
    def __init__(self, outer: str,
                 leftSpecs: ASTDeclSpecsSimple, rightSpecs: ASTDeclSpecsSimple,
                 trailing: ASTTrailingTypeSpec) -> None:
        # leftSpecs and rightSpecs are used for output
        # allSpecs are used for id generation
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

    def get_id(self, version: int) -> str:
        if version == 1:
            res = []
            res.append(self.trailingTypeSpec.get_id(version))
            if self.allSpecs.volatile:
                res.append('V')
            if self.allSpecs.const:
                res.append('C')
            return ''.join(res)
        res = []
        if self.allSpecs.volatile:
            res.append('V')
        if self.allSpecs.const:
            res.append('K')
        if self.trailingTypeSpec is not None:
            res.append(self.trailingTypeSpec.get_id(version))
        return ''.join(res)

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
        numChildren = len(signode)
        self.leftSpecs.describe_signature(signode, env, symbol)
        addSpace = len(signode) != numChildren

        if self.trailingTypeSpec:
            if addSpace:
                signode += addnodes.desc_sig_space()
            numChildren = len(signode)
            self.trailingTypeSpec.describe_signature(signode, mode, env,
                                                     symbol=symbol)
            addSpace = len(signode) != numChildren

            if len(str(self.rightSpecs)) > 0:
                if addSpace:
                    signode += addnodes.desc_sig_space()
                self.rightSpecs.describe_signature(signode, env, symbol)


# Declarator
################################################################################

class ASTArray(ASTBase):
    def __init__(self, size: ASTExpression) -> None:
        self.size = size

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTArray):
            return NotImplemented
        return self.size == other.size

    def __hash__(self) -> int:
        return hash(self.size)

    def _stringify(self, transform: StringifyTransform) -> str:
        if self.size:
            return '[' + transform(self.size) + ']'
        else:
            return '[]'

    def get_id(self, version: int) -> str:
        if version == 1:
            return 'A'
        if version == 2:
            if self.size:
                return 'A' + str(self.size) + '_'
            else:
                return 'A_'
        if self.size:
            return 'A' + self.size.get_id(version) + '_'
        else:
            return 'A_'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('[', '[')
        if self.size:
            self.size.describe_signature(signode, 'markType', env, symbol)
        signode += addnodes.desc_sig_punctuation(']', ']')


class ASTDeclarator(ASTBase):
    @property
    def name(self) -> ASTNestedName:
        raise NotImplementedError(repr(self))

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        raise NotImplementedError(repr(self))

    @property
    def isPack(self) -> bool:
        raise NotImplementedError(repr(self))

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        raise NotImplementedError(repr(self))

    @property
    def trailingReturn(self) -> ASTType:
        raise NotImplementedError(repr(self))

    def require_space_after_declSpecs(self) -> bool:
        raise NotImplementedError(repr(self))

    def get_modifiers_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def get_param_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def get_ptr_suffix_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        raise NotImplementedError(repr(self))

    def is_function_type(self) -> bool:
        raise NotImplementedError(repr(self))

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))


class ASTDeclaratorNameParamQual(ASTDeclarator):
    def __init__(self, declId: ASTNestedName,
                 arrayOps: list[ASTArray],
                 paramQual: ASTParametersQualifiers) -> None:
        self.declId = declId
        self.arrayOps = arrayOps
        self.paramQual = paramQual

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorNameParamQual):
            return NotImplemented
        return (
            self.declId == other.declId
            and self.arrayOps == other.arrayOps
            and self.paramQual == other.paramQual
        )

    def __hash__(self) -> int:
        return hash((self.declId, self.arrayOps, self.paramQual))

    @property
    def name(self) -> ASTNestedName:
        return self.declId

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.declId = name

    @property
    def isPack(self) -> bool:
        return False

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.paramQual.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.paramQual.trailingReturn

    # only the modifiers for a function, e.g.,
    def get_modifiers_id(self, version: int) -> str:
        # cv-qualifiers
        if self.paramQual:
            return self.paramQual.get_modifiers_id(version)
        raise Exception("This should only be called on a function: %s" % self)

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        if self.paramQual:
            return self.paramQual.get_param_id(version)
        else:
            return ''

    def get_ptr_suffix_id(self, version: int) -> str:  # only the array specifiers
        return ''.join(a.get_id(version) for a in self.arrayOps)

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        assert version >= 2
        res = []
        # TODO: can we actually have both array ops and paramQual?
        res.append(self.get_ptr_suffix_id(version))
        if self.paramQual:
            res.append(self.get_modifiers_id(version))
            res.append('F')
            res.append(returnTypeId)
            res.append(self.get_param_id(version))
            res.append('E')
        else:
            res.append(returnTypeId)
        return ''.join(res)

    # ------------------------------------------------------------------------

    def require_space_after_declSpecs(self) -> bool:
        return self.declId is not None

    def is_function_type(self) -> bool:
        return self.paramQual is not None

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.declId:
            res.append(transform(self.declId))
        res.extend(transform(op) for op in self.arrayOps)
        if self.paramQual:
            res.append(transform(self.paramQual))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.declId:
            self.declId.describe_signature(signode, mode, env, symbol)
        for op in self.arrayOps:
            op.describe_signature(signode, mode, env, symbol)
        if self.paramQual:
            self.paramQual.describe_signature(signode, mode, env, symbol)


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

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.declId = name

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        return ''

    def get_ptr_suffix_id(self, version: int) -> str:  # only the array specifiers
        return ''

    # ------------------------------------------------------------------------

    def require_space_after_declSpecs(self) -> bool:
        return self.declId is not None

    def is_function_type(self) -> bool:
        return False

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
    def __init__(self, next: ASTDeclarator, volatile: bool, const: bool,
                 attrs: ASTAttributeList) -> None:
        assert next
        self.next = next
        self.volatile = volatile
        self.const = const
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorPtr):
            return NotImplemented
        return (
            self.next == other.next
            and self.volatile == other.volatile
            and self.const == other.const
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((self.next, self.volatile, self.const, self.attrs))

    @property
    def name(self) -> ASTNestedName:
        return self.next.name

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.next.name = name

    @property
    def isPack(self) -> bool:
        return self.next.isPack

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.next.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.next.trailingReturn

    def require_space_after_declSpecs(self) -> bool:
        return self.next.require_space_after_declSpecs()

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['*']
        res.append(transform(self.attrs))
        if len(self.attrs) != 0 and (self.volatile or self.const):
            res.append(' ')
        if self.volatile:
            res.append('volatile')
        if self.const:
            if self.volatile:
                res.append(' ')
            res.append('const')
        if self.const or self.volatile or len(self.attrs) > 0:
            if self.next.require_space_after_declSpecs():
                res.append(' ')
        res.append(transform(self.next))
        return ''.join(res)

    def get_modifiers_id(self, version: int) -> str:
        return self.next.get_modifiers_id(version)

    def get_param_id(self, version: int) -> str:
        return self.next.get_param_id(version)

    def get_ptr_suffix_id(self, version: int) -> str:
        if version == 1:
            res = ['P']
            if self.volatile:
                res.append('V')
            if self.const:
                res.append('C')
            res.append(self.next.get_ptr_suffix_id(version))
            return ''.join(res)

        res = [self.next.get_ptr_suffix_id(version)]
        res.append('P')
        if self.volatile:
            res.append('V')
        if self.const:
            res.append('C')
        return ''.join(res)

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        # ReturnType *next, so we are part of the return type of 'next
        res = ['P']
        if self.volatile:
            res.append('V')
        if self.const:
            res.append('C')
        res.append(returnTypeId)
        return self.next.get_type_id(version, returnTypeId=''.join(res))

    def is_function_type(self) -> bool:
        return self.next.is_function_type()

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('*', '*')
        self.attrs.describe_signature(signode)
        if len(self.attrs) != 0 and (self.volatile or self.const):
            signode += addnodes.desc_sig_space()

        def _add_anno(signode: TextElement, text: str) -> None:
            signode += addnodes.desc_sig_keyword(text, text)
        if self.volatile:
            _add_anno(signode, 'volatile')
        if self.const:
            if self.volatile:
                signode += addnodes.desc_sig_space()
            _add_anno(signode, 'const')
        if self.const or self.volatile or len(self.attrs) > 0:
            if self.next.require_space_after_declSpecs():
                signode += addnodes.desc_sig_space()
        self.next.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorRef(ASTDeclarator):
    def __init__(self, next: ASTDeclarator, attrs: ASTAttributeList) -> None:
        assert next
        self.next = next
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorRef):
            return NotImplemented
        return self.next == other.next and self.attrs == other.attrs

    def __hash__(self) -> int:
        return hash((self.next, self.attrs))

    @property
    def name(self) -> ASTNestedName:
        return self.next.name

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.next.name = name

    @property
    def isPack(self) -> bool:
        return self.next.isPack

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.next.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.next.trailingReturn

    def require_space_after_declSpecs(self) -> bool:
        return self.next.require_space_after_declSpecs()

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['&']
        res.append(transform(self.attrs))
        if len(self.attrs) != 0 and self.next.require_space_after_declSpecs():
            res.append(' ')
        res.append(transform(self.next))
        return ''.join(res)

    def get_modifiers_id(self, version: int) -> str:
        return self.next.get_modifiers_id(version)

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        return self.next.get_param_id(version)

    def get_ptr_suffix_id(self, version: int) -> str:
        if version == 1:
            return 'R' + self.next.get_ptr_suffix_id(version)
        else:
            return self.next.get_ptr_suffix_id(version) + 'R'

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        assert version >= 2
        # ReturnType &next, so we are part of the return type of 'next
        return self.next.get_type_id(version, returnTypeId='R' + returnTypeId)

    def is_function_type(self) -> bool:
        return self.next.is_function_type()

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('&', '&')
        self.attrs.describe_signature(signode)
        if len(self.attrs) > 0 and self.next.require_space_after_declSpecs():
            signode += addnodes.desc_sig_space()
        self.next.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorParamPack(ASTDeclarator):
    def __init__(self, next: ASTDeclarator) -> None:
        assert next
        self.next = next

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorParamPack):
            return NotImplemented
        return self.next == other.next

    def __hash__(self) -> int:
        return hash(self.next)

    @property
    def name(self) -> ASTNestedName:
        return self.next.name

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.next.name = name

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.next.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.next.trailingReturn

    @property
    def isPack(self) -> bool:
        return True

    def require_space_after_declSpecs(self) -> bool:
        return False

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.next)
        if self.next.name:
            res = ' ' + res
        return '...' + res

    def get_modifiers_id(self, version: int) -> str:
        return self.next.get_modifiers_id(version)

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        return self.next.get_param_id(version)

    def get_ptr_suffix_id(self, version: int) -> str:
        if version == 1:
            return 'Dp' + self.next.get_ptr_suffix_id(version)
        else:
            return self.next.get_ptr_suffix_id(version) + 'Dp'

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        assert version >= 2
        # ReturnType... next, so we are part of the return type of 'next
        return self.next.get_type_id(version, returnTypeId='Dp' + returnTypeId)

    def is_function_type(self) -> bool:
        return self.next.is_function_type()

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('...', '...')
        if self.next.name:
            signode += addnodes.desc_sig_space()
        self.next.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorMemPtr(ASTDeclarator):
    def __init__(self, className: ASTNestedName,
                 const: bool, volatile: bool, next: ASTDeclarator) -> None:
        assert className
        assert next
        self.className = className
        self.const = const
        self.volatile = volatile
        self.next = next

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorMemPtr):
            return NotImplemented
        return (
            self.className == other.className
            and self.const == other.const
            and self.volatile == other.volatile
            and self.next == other.next
        )

    def __hash__(self) -> int:
        return hash((self.className, self.const, self.volatile, self.next))

    @property
    def name(self) -> ASTNestedName:
        return self.next.name

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.next.name = name

    @property
    def isPack(self) -> bool:
        return self.next.isPack

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.next.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.next.trailingReturn

    def require_space_after_declSpecs(self) -> bool:
        return True

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.className))
        res.append('::*')
        if self.volatile:
            res.append('volatile')
        if self.const:
            if self.volatile:
                res.append(' ')
            res.append('const')
        if self.next.require_space_after_declSpecs():
            res.append(' ')
        res.append(transform(self.next))
        return ''.join(res)

    def get_modifiers_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError
        return self.next.get_modifiers_id(version)

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        if version == 1:
            raise NoOldIdError
        return self.next.get_param_id(version)

    def get_ptr_suffix_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError
        raise NotImplementedError
        return self.next.get_ptr_suffix_id(version) + 'Dp'

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        assert version >= 2
        # ReturnType name::* next, so we are part of the return type of next
        nextReturnTypeId = ''
        if self.volatile:
            nextReturnTypeId += 'V'
        if self.const:
            nextReturnTypeId += 'K'
        nextReturnTypeId += 'M'
        nextReturnTypeId += self.className.get_id(version)
        nextReturnTypeId += returnTypeId
        return self.next.get_type_id(version, nextReturnTypeId)

    def is_function_type(self) -> bool:
        return self.next.is_function_type()

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.className.describe_signature(signode, 'markType', env, symbol)
        signode += addnodes.desc_sig_punctuation('::', '::')
        signode += addnodes.desc_sig_punctuation('*', '*')

        def _add_anno(signode: TextElement, text: str) -> None:
            signode += addnodes.desc_sig_keyword(text, text)
        if self.volatile:
            _add_anno(signode, 'volatile')
        if self.const:
            if self.volatile:
                signode += addnodes.desc_sig_space()
            _add_anno(signode, 'const')
        if self.next.require_space_after_declSpecs():
            signode += addnodes.desc_sig_space()
        self.next.describe_signature(signode, mode, env, symbol)


class ASTDeclaratorParen(ASTDeclarator):
    def __init__(self, inner: ASTDeclarator, next: ASTDeclarator) -> None:
        assert inner
        assert next
        self.inner = inner
        self.next = next
        # TODO: we assume the name, params, and qualifiers are in inner

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTDeclaratorParen):
            return NotImplemented
        return self.inner == other.inner and self.next == other.next

    def __hash__(self) -> int:
        return hash((self.inner, self.next))

    @property
    def name(self) -> ASTNestedName:
        return self.inner.name

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.inner.name = name

    @property
    def isPack(self) -> bool:
        return self.inner.isPack or self.next.isPack

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.inner.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.inner.trailingReturn

    def require_space_after_declSpecs(self) -> bool:
        return True

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['(']
        res.append(transform(self.inner))
        res.append(')')
        res.append(transform(self.next))
        return ''.join(res)

    def get_modifiers_id(self, version: int) -> str:
        return self.inner.get_modifiers_id(version)

    def get_param_id(self, version: int) -> str:  # only the parameters (if any)
        return self.inner.get_param_id(version)

    def get_ptr_suffix_id(self, version: int) -> str:
        if version == 1:
            raise NoOldIdError  # TODO: was this implemented before?
            return self.next.get_ptr_suffix_id(version) + \
                self.inner.get_ptr_suffix_id(version)
        return self.inner.get_ptr_suffix_id(version) + \
            self.next.get_ptr_suffix_id(version)

    def get_type_id(self, version: int, returnTypeId: str) -> str:
        assert version >= 2
        # ReturnType (inner)next, so 'inner' returns everything outside
        nextId = self.next.get_type_id(version, returnTypeId)
        return self.inner.get_type_id(version, returnTypeId=nextId)

    def is_function_type(self) -> bool:
        return self.inner.is_function_type()

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        signode += addnodes.desc_sig_punctuation('(', '(')
        self.inner.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation(')', ')')
        self.next.describe_signature(signode, "noneIsName", env, symbol)


# Type and initializer stuff
##############################################################################################

class ASTPackExpansionExpr(ASTExpression):
    def __init__(self, expr: ASTExpression | ASTBracedInitList) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTPackExpansionExpr):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.expr) + '...'

    def get_id(self, version: int) -> str:
        id = self.expr.get_id(version)
        return 'sp' + id

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.expr.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation('...', '...')


class ASTParenExprList(ASTBaseParenExprList):
    def __init__(self, exprs: list[ASTExpression | ASTBracedInitList]) -> None:
        self.exprs = exprs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTParenExprList):
            return NotImplemented
        return self.exprs == other.exprs

    def __hash__(self) -> int:
        return hash(self.exprs)

    def get_id(self, version: int) -> str:
        return "pi%sE" % ''.join(e.get_id(version) for e in self.exprs)

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


class ASTInitializer(ASTBase):
    def __init__(self, value: ASTExpression | ASTBracedInitList,
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

    @name.setter
    def name(self, name: ASTNestedName) -> None:
        self.decl.name = name

    @property
    def isPack(self) -> bool:
        return self.decl.isPack

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        return self.decl.function_params

    @property
    def trailingReturn(self) -> ASTType:
        return self.decl.trailingReturn

    def get_id(self, version: int, objectType: str | None = None,
               symbol: Symbol | None = None) -> str:
        if version == 1:
            res = []
            if objectType:  # needs the name
                if objectType == 'function':  # also modifiers
                    res.append(symbol.get_full_nested_name().get_id(version))
                    res.append(self.decl.get_param_id(version))
                    res.append(self.decl.get_modifiers_id(version))
                    if (self.declSpecs.leftSpecs.constexpr or
                            (self.declSpecs.rightSpecs and
                             self.declSpecs.rightSpecs.constexpr)):
                        res.append('CE')
                elif objectType == 'type':  # just the name
                    res.append(symbol.get_full_nested_name().get_id(version))
                else:
                    raise AssertionError(objectType)
            else:  # only type encoding
                if self.decl.is_function_type():
                    raise NoOldIdError
                res.append(self.declSpecs.get_id(version))
                res.append(self.decl.get_ptr_suffix_id(version))
                res.append(self.decl.get_param_id(version))
            return ''.join(res)
        # other versions
        res = []
        if objectType:  # needs the name
            if objectType == 'function':  # also modifiers
                modifiers = self.decl.get_modifiers_id(version)
                res.append(symbol.get_full_nested_name().get_id(version, modifiers))
                if version >= 4:
                    # with templates we need to mangle the return type in as well
                    templ = symbol.declaration.templatePrefix
                    if templ is not None:
                        typeId = self.decl.get_ptr_suffix_id(version)
                        if self.trailingReturn:
                            returnTypeId = self.trailingReturn.get_id(version)
                        else:
                            returnTypeId = self.declSpecs.get_id(version)
                        res.append(typeId)
                        res.append(returnTypeId)
                res.append(self.decl.get_param_id(version))
            elif objectType == 'type':  # just the name
                res.append(symbol.get_full_nested_name().get_id(version))
            else:
                raise AssertionError(objectType)
        else:  # only type encoding
            # the 'returnType' of a non-function type is simply just the last
            # type, i.e., for 'int*' it is 'int'
            returnTypeId = self.declSpecs.get_id(version)
            typeId = self.decl.get_type_id(version, returnTypeId)
            res.append(typeId)
        return ''.join(res)

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


class ASTTemplateParamConstrainedTypeWithInit(ASTBase):
    def __init__(self, type: ASTType, init: ASTType) -> None:
        assert type
        self.type = type
        self.init = init

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateParamConstrainedTypeWithInit):
            return NotImplemented
        return self.type == other.type and self.init == other.init

    def __hash__(self) -> int:
        return hash((self.type, self.init))

    @property
    def name(self) -> ASTNestedName:
        return self.type.name

    @property
    def isPack(self) -> bool:
        return self.type.isPack

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        # this is not part of the normal name mangling in C++
        assert version >= 2
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=False)
        else:
            return self.type.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.type)
        if self.init:
            res += " = "
            res += transform(self.init)
        return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.type.describe_signature(signode, mode, env, symbol)
        if self.init:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('=', '=')
            signode += addnodes.desc_sig_space()
            self.init.describe_signature(signode, mode, env, symbol)


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

    @property
    def isPack(self) -> bool:
        return self.type.isPack

    def get_id(self, version: int, objectType: str | None = None,
               symbol: Symbol | None = None) -> str:
        if objectType != 'member':
            return self.type.get_id(version, objectType)
        if version == 1:
            return (symbol.get_full_nested_name().get_id(version) + '__' +
                    self.type.get_id(version))
        return symbol.get_full_nested_name().get_id(version)

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


class ASTTypeUsing(ASTBase):
    def __init__(self, name: ASTNestedName, type: ASTType | None) -> None:
        self.name = name
        self.type = type

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTypeUsing):
            return NotImplemented
        return self.name == other.name and self.type == other.type

    def __hash__(self) -> int:
        return hash((self.name, self.type))

    def get_id(self, version: int, objectType: str | None = None,
               symbol: Symbol | None = None) -> str:
        if version == 1:
            raise NoOldIdError
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.name))
        if self.type:
            res.append(' = ')
            res.append(transform(self.type))
        return ''.join(res)

    def get_type_declaration_prefix(self) -> str:
        return 'using'

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.name.describe_signature(signode, mode, env, symbol=symbol)
        if self.type:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('=', '=')
            signode += addnodes.desc_sig_space()
            self.type.describe_signature(signode, 'markType', env, symbol=symbol)


# Other declarations
##############################################################################################

class ASTConcept(ASTBase):
    def __init__(self, nestedName: ASTNestedName, initializer: ASTInitializer) -> None:
        self.nestedName = nestedName
        self.initializer = initializer

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTConcept):
            return NotImplemented
        return self.nestedName == other.nestedName and self.initializer == other.initializer

    def __hash__(self) -> int:
        return hash((self.nestedName, self.initializer))

    @property
    def name(self) -> ASTNestedName:
        return self.nestedName

    def get_id(self, version: int, objectType: str | None = None,
               symbol: Symbol | None = None) -> str:
        if version == 1:
            raise NoOldIdError
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.nestedName)
        if self.initializer:
            res += transform(self.initializer)
        return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.nestedName.describe_signature(signode, mode, env, symbol)
        if self.initializer:
            self.initializer.describe_signature(signode, mode, env, symbol)


class ASTBaseClass(ASTBase):
    def __init__(self, name: ASTNestedName, visibility: str,
                 virtual: bool, pack: bool) -> None:
        self.name = name
        self.visibility = visibility
        self.virtual = virtual
        self.pack = pack

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTBaseClass):
            return NotImplemented
        return (
            self.name == other.name
            and self.visibility == other.visibility
            and self.virtual == other.virtual
            and self.pack == other.pack
        )

    def __hash__(self) -> int:
        return hash((self.name, self.visibility, self.virtual, self.pack))

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.visibility is not None:
            res.append(self.visibility)
            res.append(' ')
        if self.virtual:
            res.append('virtual ')
        res.append(transform(self.name))
        if self.pack:
            res.append('...')
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        if self.visibility is not None:
            signode += addnodes.desc_sig_keyword(self.visibility,
                                                 self.visibility)
            signode += addnodes.desc_sig_space()
        if self.virtual:
            signode += addnodes.desc_sig_keyword('virtual', 'virtual')
            signode += addnodes.desc_sig_space()
        self.name.describe_signature(signode, 'markType', env, symbol=symbol)
        if self.pack:
            signode += addnodes.desc_sig_punctuation('...', '...')


class ASTClass(ASTBase):
    def __init__(self, name: ASTNestedName, final: bool, bases: list[ASTBaseClass],
                 attrs: ASTAttributeList) -> None:
        self.name = name
        self.final = final
        self.bases = bases
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTClass):
            return NotImplemented
        return (
            self.name == other.name
            and self.final == other.final
            and self.bases == other.bases
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((self.name, self.final, self.bases, self.attrs))

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.attrs))
        if len(self.attrs) != 0:
            res.append(' ')
        res.append(transform(self.name))
        if self.final:
            res.append(' final')
        if len(self.bases) > 0:
            res.append(' : ')
            first = True
            for b in self.bases:
                if not first:
                    res.append(', ')
                first = False
                res.append(transform(b))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.attrs.describe_signature(signode)
        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
        self.name.describe_signature(signode, mode, env, symbol=symbol)
        if self.final:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_keyword('final', 'final')
        if len(self.bases) > 0:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation(':', ':')
            signode += addnodes.desc_sig_space()
            for b in self.bases:
                b.describe_signature(signode, mode, env, symbol=symbol)
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            signode.pop()
            signode.pop()


class ASTUnion(ASTBase):
    def __init__(self, name: ASTNestedName, attrs: ASTAttributeList) -> None:
        self.name = name
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTUnion):
            return NotImplemented
        return self.name == other.name and self.attrs == other.attrs

    def __hash__(self) -> int:
        return hash((self.name, self.attrs))

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        if version == 1:
            raise NoOldIdError
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.attrs))
        if len(self.attrs) != 0:
            res.append(' ')
        res.append(transform(self.name))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        self.attrs.describe_signature(signode)
        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
        self.name.describe_signature(signode, mode, env, symbol=symbol)


class ASTEnum(ASTBase):
    def __init__(self, name: ASTNestedName, scoped: str, underlyingType: ASTType,
                 attrs: ASTAttributeList) -> None:
        self.name = name
        self.scoped = scoped
        self.underlyingType = underlyingType
        self.attrs = attrs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTEnum):
            return NotImplemented
        return (
            self.name == other.name
            and self.scoped == other.scoped
            and self.underlyingType == other.underlyingType
            and self.attrs == other.attrs
        )

    def __hash__(self) -> int:
        return hash((self.name, self.scoped, self.underlyingType, self.attrs))

    def get_id(self, version: int, objectType: str, symbol: Symbol) -> str:
        if version == 1:
            raise NoOldIdError
        return symbol.get_full_nested_name().get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.scoped:
            res.append(self.scoped)
            res.append(' ')
        res.append(transform(self.attrs))
        if len(self.attrs) != 0:
            res.append(' ')
        res.append(transform(self.name))
        if self.underlyingType:
            res.append(' : ')
            res.append(transform(self.underlyingType))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        verify_description_mode(mode)
        # self.scoped has been done by the CPPEnumObject
        self.attrs.describe_signature(signode)
        if len(self.attrs) != 0:
            signode += addnodes.desc_sig_space()
        self.name.describe_signature(signode, mode, env, symbol=symbol)
        if self.underlyingType:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation(':', ':')
            signode += addnodes.desc_sig_space()
            self.underlyingType.describe_signature(signode, 'noneIsName',
                                                   env, symbol=symbol)


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
        if version == 1:
            raise NoOldIdError
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


################################################################################
# Templates
################################################################################

# Parameters
################################################################################

class ASTTemplateParam(ASTBase):
    def get_identifier(self) -> ASTIdentifier:
        raise NotImplementedError(repr(self))

    def get_id(self, version: int) -> str:
        raise NotImplementedError(repr(self))

    def describe_signature(self, parentNode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        raise NotImplementedError(repr(self))

    @property
    def isPack(self) -> bool:
        raise NotImplementedError(repr(self))

    @property
    def name(self) -> ASTNestedName:
        raise NotImplementedError(repr(self))


class ASTTemplateKeyParamPackIdDefault(ASTTemplateParam):
    def __init__(self, key: str, identifier: ASTIdentifier,
                 parameterPack: bool, default: ASTType) -> None:
        assert key
        if parameterPack:
            assert default is None
        self.key = key
        self.identifier = identifier
        self.parameterPack = parameterPack
        self.default = default

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateKeyParamPackIdDefault):
            return NotImplemented
        return (
            self.key == other.key
            and self.identifier == other.identifier
            and self.parameterPack == other.parameterPack
            and self.default == other.default
        )

    def __hash__(self) -> int:
        return hash((self.key, self.identifier, self.parameterPack, self.default))

    def get_identifier(self) -> ASTIdentifier:
        return self.identifier

    def get_id(self, version: int) -> str:
        assert version >= 2
        # this is not part of the normal name mangling in C++
        res = []
        if self.parameterPack:
            res.append('Dp')
        else:
            res.append('0')  # we need to put something
        return ''.join(res)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = [self.key]
        if self.parameterPack:
            if self.identifier:
                res.append(' ')
            res.append('...')
        if self.identifier:
            if not self.parameterPack:
                res.append(' ')
            res.append(transform(self.identifier))
        if self.default:
            res.append(' = ')
            res.append(transform(self.default))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword(self.key, self.key)
        if self.parameterPack:
            if self.identifier:
                signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('...', '...')
        if self.identifier:
            if not self.parameterPack:
                signode += addnodes.desc_sig_space()
            self.identifier.describe_signature(signode, mode, env, '', '', symbol)
        if self.default:
            signode += addnodes.desc_sig_space()
            signode += addnodes.desc_sig_punctuation('=', '=')
            signode += addnodes.desc_sig_space()
            self.default.describe_signature(signode, 'markType', env, symbol)


class ASTTemplateParamType(ASTTemplateParam):
    def __init__(self, data: ASTTemplateKeyParamPackIdDefault) -> None:
        assert data
        self.data = data

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateParamType):
            return NotImplemented
        return self.data == other.data

    def __hash__(self) -> int:
        return hash(self.data)

    @property
    def name(self) -> ASTNestedName:
        id = self.get_identifier()
        return ASTNestedName([ASTNestedNameElement(id, None)], [False], rooted=False)

    @property
    def isPack(self) -> bool:
        return self.data.parameterPack

    def get_identifier(self) -> ASTIdentifier:
        return self.data.get_identifier()

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        # this is not part of the normal name mangling in C++
        assert version >= 2
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=False)
        else:
            return self.data.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.data)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.data.describe_signature(signode, mode, env, symbol)


class ASTTemplateParamTemplateType(ASTTemplateParam):
    def __init__(self, nestedParams: ASTTemplateParams,
                 data: ASTTemplateKeyParamPackIdDefault) -> None:
        assert nestedParams
        assert data
        self.nestedParams = nestedParams
        self.data = data

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateParamTemplateType):
            return NotImplemented
        return (
            self.nestedParams == other.nestedParams
            and self.data == other.data
        )

    def __hash__(self) -> int:
        return hash((self.nestedParams, self.data))

    @property
    def name(self) -> ASTNestedName:
        id = self.get_identifier()
        return ASTNestedName([ASTNestedNameElement(id, None)], [False], rooted=False)

    @property
    def isPack(self) -> bool:
        return self.data.parameterPack

    def get_identifier(self) -> ASTIdentifier:
        return self.data.get_identifier()

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        assert version >= 2
        # this is not part of the normal name mangling in C++
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=None)
        else:
            return self.nestedParams.get_id(version) + self.data.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        return transform(self.nestedParams) + transform(self.data)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.nestedParams.describe_signature(signode, 'noneIsName', env, symbol)
        signode += addnodes.desc_sig_space()
        self.data.describe_signature(signode, mode, env, symbol)


class ASTTemplateParamNonType(ASTTemplateParam):
    def __init__(self,
                 param: ASTTypeWithInit | ASTTemplateParamConstrainedTypeWithInit,
                 parameterPack: bool = False) -> None:
        assert param
        self.param = param
        self.parameterPack = parameterPack

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateParamNonType):
            return NotImplemented
        return (
            self.param == other.param
            and self.parameterPack == other.parameterPack
        )

    @property
    def name(self) -> ASTNestedName:
        id = self.get_identifier()
        return ASTNestedName([ASTNestedNameElement(id, None)], [False], rooted=False)

    @property
    def isPack(self) -> bool:
        return self.param.isPack or self.parameterPack

    def get_identifier(self) -> ASTIdentifier:
        name = self.param.name
        if name:
            assert len(name.names) == 1
            assert name.names[0].identOrOp
            assert not name.names[0].templateArgs
            res = name.names[0].identOrOp
            assert isinstance(res, ASTIdentifier)
            return res
        else:
            return None

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        assert version >= 2
        # this is not part of the normal name mangling in C++
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=None)
        else:
            res = '_'
            if self.parameterPack:
                res += 'Dp'
            return res + self.param.get_id(version)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = transform(self.param)
        if self.parameterPack:
            res += '...'
        return res

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        self.param.describe_signature(signode, mode, env, symbol)
        if self.parameterPack:
            signode += addnodes.desc_sig_punctuation('...', '...')


class ASTTemplateParams(ASTBase):
    def __init__(self, params: list[ASTTemplateParam],
                 requiresClause: ASTRequiresClause | None) -> None:
        assert params is not None
        self.params = params
        self.requiresClause = requiresClause

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateParams):
            return NotImplemented
        return self.params == other.params and self.requiresClause == other.requiresClause

    def __hash__(self) -> int:
        return hash((self.params, self.requiresClause))

    def get_id(self, version: int, excludeRequires: bool = False) -> str:
        assert version >= 2
        res = []
        res.append("I")
        res.extend(param.get_id(version) for param in self.params)
        res.append("E")
        if not excludeRequires and self.requiresClause:
            res.extend(['IQ', self.requiresClause.expr.get_id(version), 'E'])
        return ''.join(res)

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append("template<")
        res.append(", ".join(transform(a) for a in self.params))
        res.append("> ")
        if self.requiresClause is not None:
            res.append(transform(self.requiresClause))
            res.append(" ")
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('template', 'template')
        signode += addnodes.desc_sig_punctuation('<', '<')
        first = True
        for param in self.params:
            if not first:
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            first = False
            param.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation('>', '>')
        if self.requiresClause is not None:
            signode += addnodes.desc_sig_space()
            self.requiresClause.describe_signature(signode, mode, env, symbol)

    def describe_signature_as_introducer(
            self, parentNode: desc_signature, mode: str, env: BuildEnvironment,
            symbol: Symbol, lineSpec: bool) -> None:
        def makeLine(parentNode: desc_signature) -> addnodes.desc_signature_line:
            signode = addnodes.desc_signature_line()
            parentNode += signode
            signode.sphinx_line_type = 'templateParams'
            return signode
        lineNode = makeLine(parentNode)
        lineNode += addnodes.desc_sig_keyword('template', 'template')
        lineNode += addnodes.desc_sig_punctuation('<', '<')
        first = True
        for param in self.params:
            if not first:
                lineNode += addnodes.desc_sig_punctuation(',', ',')
                lineNode += addnodes.desc_sig_space()
            first = False
            if lineSpec:
                lineNode = makeLine(parentNode)
            param.describe_signature(lineNode, mode, env, symbol)
        if lineSpec and not first:
            lineNode = makeLine(parentNode)
        lineNode += addnodes.desc_sig_punctuation('>', '>')
        if self.requiresClause:
            reqNode = addnodes.desc_signature_line()
            reqNode.sphinx_line_type = 'requiresClause'
            parentNode += reqNode
            self.requiresClause.describe_signature(reqNode, 'markType', env, symbol)


# Template introducers
################################################################################

class ASTTemplateIntroductionParameter(ASTBase):
    def __init__(self, identifier: ASTIdentifier, parameterPack: bool) -> None:
        self.identifier = identifier
        self.parameterPack = parameterPack

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateIntroductionParameter):
            return NotImplemented
        return (
            self.identifier == other.identifier
            and self.parameterPack == other.parameterPack
        )

    def __hash__(self) -> int:
        return hash((self.identifier, self.parameterPack))

    @property
    def name(self) -> ASTNestedName:
        id = self.get_identifier()
        return ASTNestedName([ASTNestedNameElement(id, None)], [False], rooted=False)

    @property
    def isPack(self) -> bool:
        return self.parameterPack

    def get_identifier(self) -> ASTIdentifier:
        return self.identifier

    def get_id(
        self, version: int, objectType: str | None = None, symbol: Symbol | None = None,
    ) -> str:
        assert version >= 2
        # this is not part of the normal name mangling in C++
        if symbol:
            # the anchor will be our parent
            return symbol.parent.declaration.get_id(version, prefixed=None)
        else:
            if self.parameterPack:
                return 'Dp'
            else:
                return '0'  # we need to put something

    def get_id_as_arg(self, version: int) -> str:
        assert version >= 2
        # used for the implicit requires clause
        res = self.identifier.get_id(version)
        if self.parameterPack:
            return 'sp' + res
        else:
            return res

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.parameterPack:
            res.append('...')
        res.append(transform(self.identifier))
        return ''.join(res)

    def describe_signature(self, signode: TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        if self.parameterPack:
            signode += addnodes.desc_sig_punctuation('...', '...')
        self.identifier.describe_signature(signode, mode, env, '', '', symbol)


class ASTTemplateIntroduction(ASTBase):
    def __init__(self, concept: ASTNestedName,
                 params: list[ASTTemplateIntroductionParameter]) -> None:
        assert len(params) > 0
        self.concept = concept
        self.params = params

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateIntroduction):
            return NotImplemented
        return self.concept == other.concept and self.params == other.params

    def __hash__(self) -> int:
        return hash((self.concept, self.params))

    def get_id(self, version: int) -> str:
        assert version >= 2
        return ''.join([
            # first do the same as a normal template parameter list
            "I",
            *(param.get_id(version) for param in self.params),
            "E",
            # let's use X expr E, which is otherwise for constant template args
            "X",
            self.concept.get_id(version),
            "I",
            *(param.get_id_as_arg(version) for param in self.params),
            "E",
            "E",
        ])

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        res.append(transform(self.concept))
        res.append('{')
        res.append(', '.join(transform(param) for param in self.params))
        res.append('} ')
        return ''.join(res)

    def describe_signature_as_introducer(
            self, parentNode: desc_signature, mode: str,
            env: BuildEnvironment, symbol: Symbol, lineSpec: bool) -> None:
        # Note: 'lineSpec' has no effect on template introductions.
        signode = addnodes.desc_signature_line()
        parentNode += signode
        signode.sphinx_line_type = 'templateIntroduction'
        self.concept.describe_signature(signode, 'markType', env, symbol)
        signode += addnodes.desc_sig_punctuation('{', '{')
        first = True
        for param in self.params:
            if not first:
                signode += addnodes.desc_sig_punctuation(',', ',')
                signode += addnodes.desc_sig_space()
            first = False
            param.describe_signature(signode, mode, env, symbol)
        signode += addnodes.desc_sig_punctuation('}', '}')


################################################################################

class ASTTemplateDeclarationPrefix(ASTBase):
    def __init__(self,
                 templates: list[ASTTemplateParams | ASTTemplateIntroduction] | None) -> None:
        # templates is None means it's an explicit instantiation of a variable
        self.templates = templates

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTTemplateDeclarationPrefix):
            return NotImplemented
        return self.templates == other.templates

    def __hash__(self) -> int:
        return hash(self.templates)

    def get_requires_clause_in_last(self) -> ASTRequiresClause | None:
        if self.templates is None:
            return None
        lastList = self.templates[-1]
        if not isinstance(lastList, ASTTemplateParams):
            return None
        return lastList.requiresClause  # which may be None

    def get_id_except_requires_clause_in_last(self, version: int) -> str:
        assert version >= 2
        # This is not part of the Itanium ABI mangling system.
        res = []
        lastIndex = len(self.templates) - 1
        for i, t in enumerate(self.templates):
            if isinstance(t, ASTTemplateParams):
                res.append(t.get_id(version, excludeRequires=(i == lastIndex)))
            else:
                res.append(t.get_id(version))
        return ''.join(res)

    def _stringify(self, transform: StringifyTransform) -> str:
        return ''.join(map(transform, self.templates))

    def describe_signature(self, signode: desc_signature, mode: str,
                           env: BuildEnvironment, symbol: Symbol, lineSpec: bool) -> None:
        verify_description_mode(mode)
        for t in self.templates:
            t.describe_signature_as_introducer(signode, 'lastIsName', env, symbol, lineSpec)


class ASTRequiresClause(ASTBase):
    def __init__(self, expr: ASTExpression) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTRequiresClause):
            return NotImplemented
        return self.expr == other.expr

    def __hash__(self) -> int:
        return hash(self.expr)

    def _stringify(self, transform: StringifyTransform) -> str:
        return 'requires ' + transform(self.expr)

    def describe_signature(self, signode: nodes.TextElement, mode: str,
                           env: BuildEnvironment, symbol: Symbol) -> None:
        signode += addnodes.desc_sig_keyword('requires', 'requires')
        signode += addnodes.desc_sig_space()
        self.expr.describe_signature(signode, mode, env, symbol)


################################################################################
################################################################################

class ASTDeclaration(ASTBase):
    def __init__(self, objectType: str, directiveType: str | None = None,
                 visibility: str | None = None,
                 templatePrefix: ASTTemplateDeclarationPrefix | None = None,
                 declaration: Any = None,
                 trailingRequiresClause: ASTRequiresClause | None = None,
                 semicolon: bool = False) -> None:
        self.objectType = objectType
        self.directiveType = directiveType
        self.visibility = visibility
        self.templatePrefix = templatePrefix
        self.declaration = declaration
        self.trailingRequiresClause = trailingRequiresClause
        self.semicolon = semicolon

        self.symbol: Symbol | None = None
        # set by CPPObject._add_enumerator_to_parent
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
            and self.visibility == other.visibility
            and self.templatePrefix == other.templatePrefix
            and self.declaration == other.declaration
            and self.trailingRequiresClause == other.trailingRequiresClause
            and self.semicolon == other.semicolon
            and self.symbol == other.symbol
            and self.enumeratorScopedSymbol == other.enumeratorScopedSymbol
        )

    def clone(self) -> ASTDeclaration:
        templatePrefixClone = self.templatePrefix.clone() if self.templatePrefix else None
        trailingRequiresClasueClone = self.trailingRequiresClause.clone() \
            if self.trailingRequiresClause else None
        return ASTDeclaration(self.objectType, self.directiveType, self.visibility,
                              templatePrefixClone,
                              self.declaration.clone(), trailingRequiresClasueClone,
                              self.semicolon)

    @property
    def name(self) -> ASTNestedName:
        return self.declaration.name

    @property
    def function_params(self) -> list[ASTFunctionParameter]:
        if self.objectType != 'function':
            return None
        return self.declaration.function_params

    def get_id(self, version: int, prefixed: bool = True) -> str:
        if version == 1:
            if self.templatePrefix or self.trailingRequiresClause:
                raise NoOldIdError
            if self.objectType == 'enumerator' and self.enumeratorScopedSymbol:
                return self.enumeratorScopedSymbol.declaration.get_id(version)
            return self.declaration.get_id(version, self.objectType, self.symbol)
        # version >= 2
        if self.objectType == 'enumerator' and self.enumeratorScopedSymbol:
            return self.enumeratorScopedSymbol.declaration.get_id(version, prefixed)
        if prefixed:
            res = [_id_prefix[version]]
        else:
            res = []
        # (See also https://github.com/sphinx-doc/sphinx/pull/10286#issuecomment-1168102147)
        # The first implementation of requires clauses only supported a single clause after the
        # template prefix, and no trailing clause. It put the ID after the template parameter
        # list, i.e.,
        #    "I" + template_parameter_list_id + "E" + "IQ" + requires_clause_id + "E"
        # but the second implementation associates the requires clause with each list, i.e.,
        #    "I" + template_parameter_list_id + "IQ" + requires_clause_id + "E" + "E"
        # To avoid making a new ID version, we make an exception for the last requires clause
        # in the template prefix, and still put it in the end.
        # As we now support trailing requires clauses we add that as if it was a conjunction.
        if self.templatePrefix is not None:
            res.append(self.templatePrefix.get_id_except_requires_clause_in_last(version))
            requiresClauseInLast = self.templatePrefix.get_requires_clause_in_last()
        else:
            requiresClauseInLast = None

        if requiresClauseInLast or self.trailingRequiresClause:
            if version < 4:
                raise NoOldIdError
            res.append('IQ')
            if requiresClauseInLast and self.trailingRequiresClause:
                # make a conjunction of them
                res.append('aa')
            if requiresClauseInLast:
                res.append(requiresClauseInLast.expr.get_id(version))
            if self.trailingRequiresClause:
                res.append(self.trailingRequiresClause.expr.get_id(version))
            res.append('E')
        res.append(self.declaration.get_id(version, self.objectType, self.symbol))
        return ''.join(res)

    def get_newest_id(self) -> str:
        if self._newest_id_cache is None:
            self._newest_id_cache = self.get_id(_max_id, True)
        return self._newest_id_cache

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.visibility and self.visibility != "public":
            res.append(self.visibility)
            res.append(' ')
        if self.templatePrefix:
            res.append(transform(self.templatePrefix))
        res.append(transform(self.declaration))
        if self.trailingRequiresClause:
            res.append(' ')
            res.append(transform(self.trailingRequiresClause))
        if self.semicolon:
            res.append(';')
        return ''.join(res)

    def describe_signature(self, signode: desc_signature, mode: str,
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

        if self.templatePrefix:
            self.templatePrefix.describe_signature(signode, mode, env,
                                                   symbol=self.symbol,
                                                   lineSpec=options.get('tparam-line-spec'))
        signode += mainDeclNode
        if self.visibility and self.visibility != "public":
            mainDeclNode += addnodes.desc_sig_keyword(self.visibility, self.visibility)
            mainDeclNode += addnodes.desc_sig_space()
        if self.objectType == 'type':
            prefix = self.declaration.get_type_declaration_prefix()
            mainDeclNode += addnodes.desc_sig_keyword(prefix, prefix)
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'concept':
            mainDeclNode += addnodes.desc_sig_keyword('concept', 'concept')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType in {'member', 'function'}:
            pass
        elif self.objectType == 'class':
            assert self.directiveType in ('class', 'struct')
            mainDeclNode += addnodes.desc_sig_keyword(self.directiveType, self.directiveType)
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'union':
            mainDeclNode += addnodes.desc_sig_keyword('union', 'union')
            mainDeclNode += addnodes.desc_sig_space()
        elif self.objectType == 'enum':
            mainDeclNode += addnodes.desc_sig_keyword('enum', 'enum')
            mainDeclNode += addnodes.desc_sig_space()
            if self.directiveType == 'enum-class':
                mainDeclNode += addnodes.desc_sig_keyword('class', 'class')
                mainDeclNode += addnodes.desc_sig_space()
            elif self.directiveType == 'enum-struct':
                mainDeclNode += addnodes.desc_sig_keyword('struct', 'struct')
                mainDeclNode += addnodes.desc_sig_space()
            else:
                assert self.directiveType == 'enum', self.directiveType
        elif self.objectType == 'enumerator':
            mainDeclNode += addnodes.desc_sig_keyword('enumerator', 'enumerator')
            mainDeclNode += addnodes.desc_sig_space()
        else:
            raise AssertionError(self.objectType)
        self.declaration.describe_signature(mainDeclNode, mode, env, self.symbol)
        lastDeclNode = mainDeclNode
        if self.trailingRequiresClause:
            trailingReqNode = addnodes.desc_signature_line()
            trailingReqNode.sphinx_line_type = 'trailingRequiresClause'
            signode.append(trailingReqNode)
            lastDeclNode = trailingReqNode
            self.trailingRequiresClause.describe_signature(
                trailingReqNode, 'markType', env, self.symbol)
        if self.semicolon:
            lastDeclNode += addnodes.desc_sig_punctuation(';', ';')


class ASTNamespace(ASTBase):
    def __init__(self, nestedName: ASTNestedName,
                 templatePrefix: ASTTemplateDeclarationPrefix) -> None:
        self.nestedName = nestedName
        self.templatePrefix = templatePrefix

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ASTNamespace):
            return NotImplemented
        return (
            self.nestedName == other.nestedName
            and self.templatePrefix == other.templatePrefix
        )

    def _stringify(self, transform: StringifyTransform) -> str:
        res = []
        if self.templatePrefix:
            res.append(transform(self.templatePrefix))
        res.append(transform(self.nestedName))
        return ''.join(res)
