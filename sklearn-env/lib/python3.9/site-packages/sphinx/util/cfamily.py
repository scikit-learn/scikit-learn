"""
    sphinx.util.cfamily
    ~~~~~~~~~~~~~~~~~~~

    Utility functions common to the C and C++ domains.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
from copy import deepcopy
from typing import Any, Callable, List, Match, Optional, Pattern, Tuple, Union

from docutils import nodes
from docutils.nodes import TextElement

from sphinx.config import Config
from sphinx.util import logging

logger = logging.getLogger(__name__)

StringifyTransform = Callable[[Any], str]


_whitespace_re = re.compile(r'(?u)\s+')
anon_identifier_re = re.compile(r'(@[a-zA-Z0-9_])[a-zA-Z0-9_]*\b')
identifier_re = re.compile(r'''(?x)
    (   # This 'extends' _anon_identifier_re with the ordinary identifiers,
        # make sure they are in sync.
        (~?\b[a-zA-Z_])  # ordinary identifiers
    |   (@[a-zA-Z0-9_])  # our extension for names of anonymous entities
    )
    [a-zA-Z0-9_]*\b
''')
integer_literal_re = re.compile(r'[1-9][0-9]*(\'[0-9]+)*')
octal_literal_re = re.compile(r'0[0-7]*(\'[0-7]+)*')
hex_literal_re = re.compile(r'0[xX][0-9a-fA-F]+(\'[0-9a-fA-F]+)*')
binary_literal_re = re.compile(r'0[bB][01]+(\'[01]+)*')
integers_literal_suffix_re = re.compile(r'''(?x)
    # unsigned and/or (long) long, in any order, but at least one of them
    (
        ([uU]    ([lL]  |  (ll)  |  (LL))?)
        |
        (([lL]  |  (ll)  |  (LL))    [uU]?)
    )\b
    # the ending word boundary is important for distinguishing
    # between suffixes and UDLs in C++
''')
float_literal_re = re.compile(r'''(?x)
    [+-]?(
    # decimal
      ([0-9]+(\'[0-9]+)*[eE][+-]?[0-9]+(\'[0-9]+)*)
    | (([0-9]+(\'[0-9]+)*)?\.[0-9]+(\'[0-9]+)*([eE][+-]?[0-9]+(\'[0-9]+)*)?)
    | ([0-9]+(\'[0-9]+)*\.([eE][+-]?[0-9]+(\'[0-9]+)*)?)
    # hex
    | (0[xX][0-9a-fA-F]+(\'[0-9a-fA-F]+)*[pP][+-]?[0-9a-fA-F]+(\'[0-9a-fA-F]+)*)
    | (0[xX]([0-9a-fA-F]+(\'[0-9a-fA-F]+)*)?\.
        [0-9a-fA-F]+(\'[0-9a-fA-F]+)*([pP][+-]?[0-9a-fA-F]+(\'[0-9a-fA-F]+)*)?)
    | (0[xX][0-9a-fA-F]+(\'[0-9a-fA-F]+)*\.([pP][+-]?[0-9a-fA-F]+(\'[0-9a-fA-F]+)*)?)
    )
''')
float_literal_suffix_re = re.compile(r'[fFlL]\b')
# the ending word boundary is important for distinguishing between suffixes and UDLs in C++
char_literal_re = re.compile(r'''(?x)
    ((?:u8)|u|U|L)?
    '(
      (?:[^\\'])
    | (\\(
        (?:['"?\\abfnrtv])
      | (?:[0-7]{1,3})
      | (?:x[0-9a-fA-F]{2})
      | (?:u[0-9a-fA-F]{4})
      | (?:U[0-9a-fA-F]{8})
      ))
    )'
''')


def verify_description_mode(mode: str) -> None:
    if mode not in ('lastIsName', 'noneIsName', 'markType', 'markName', 'param', 'udl'):
        raise Exception("Description mode '%s' is invalid." % mode)


class NoOldIdError(Exception):
    # Used to avoid implementing unneeded id generation for old id schemes.
    pass


class ASTBaseBase:
    def __eq__(self, other: Any) -> bool:
        if type(self) is not type(other):
            return False
        try:
            for key, value in self.__dict__.items():
                if value != getattr(other, key):
                    return False
        except AttributeError:
            return False
        return True

    __hash__: Callable[[], int] = None

    def clone(self) -> Any:
        return deepcopy(self)

    def _stringify(self, transform: StringifyTransform) -> str:
        raise NotImplementedError(repr(self))

    def __str__(self) -> str:
        return self._stringify(lambda ast: str(ast))

    def get_display_string(self) -> str:
        return self._stringify(lambda ast: ast.get_display_string())

    def __repr__(self) -> str:
        return '<%s>' % self.__class__.__name__


################################################################################
# Attributes
################################################################################

class ASTAttribute(ASTBaseBase):
    def describe_signature(self, signode: TextElement) -> None:
        raise NotImplementedError(repr(self))


class ASTCPPAttribute(ASTAttribute):
    def __init__(self, arg: str) -> None:
        self.arg = arg

    def _stringify(self, transform: StringifyTransform) -> str:
        return "[[" + self.arg + "]]"

    def describe_signature(self, signode: TextElement) -> None:
        txt = str(self)
        signode.append(nodes.Text(txt, txt))


class ASTGnuAttribute(ASTBaseBase):
    def __init__(self, name: str, args: Optional["ASTBaseParenExprList"]) -> None:
        self.name = name
        self.args = args

    def _stringify(self, transform: StringifyTransform) -> str:
        res = [self.name]
        if self.args:
            res.append(transform(self.args))
        return ''.join(res)


class ASTGnuAttributeList(ASTAttribute):
    def __init__(self, attrs: List[ASTGnuAttribute]) -> None:
        self.attrs = attrs

    def _stringify(self, transform: StringifyTransform) -> str:
        res = ['__attribute__((']
        first = True
        for attr in self.attrs:
            if not first:
                res.append(', ')
            first = False
            res.append(transform(attr))
        res.append('))')
        return ''.join(res)

    def describe_signature(self, signode: TextElement) -> None:
        txt = str(self)
        signode.append(nodes.Text(txt, txt))


class ASTIdAttribute(ASTAttribute):
    """For simple attributes defined by the user."""

    def __init__(self, id: str) -> None:
        self.id = id

    def _stringify(self, transform: StringifyTransform) -> str:
        return self.id

    def describe_signature(self, signode: TextElement) -> None:
        signode.append(nodes.Text(self.id, self.id))


class ASTParenAttribute(ASTAttribute):
    """For paren attributes defined by the user."""

    def __init__(self, id: str, arg: str) -> None:
        self.id = id
        self.arg = arg

    def _stringify(self, transform: StringifyTransform) -> str:
        return self.id + '(' + self.arg + ')'

    def describe_signature(self, signode: TextElement) -> None:
        txt = str(self)
        signode.append(nodes.Text(txt, txt))


################################################################################

class ASTBaseParenExprList(ASTBaseBase):
    pass


################################################################################

class UnsupportedMultiCharacterCharLiteral(Exception):
    pass


class DefinitionError(Exception):
    pass


class BaseParser:
    def __init__(self, definition: str, *,
                 location: Union[nodes.Node, Tuple[str, int], str],
                 config: "Config") -> None:
        self.definition = definition.strip()
        self.location = location  # for warnings
        self.config = config

        self.pos = 0
        self.end = len(self.definition)
        self.last_match: Match = None
        self._previous_state: Tuple[int, Match] = (0, None)
        self.otherErrors: List[DefinitionError] = []

        # in our tests the following is set to False to capture bad parsing
        self.allowFallbackExpressionParsing = True

    def _make_multi_error(self, errors: List[Any], header: str) -> DefinitionError:
        if len(errors) == 1:
            if len(header) > 0:
                return DefinitionError(header + '\n' + str(errors[0][0]))
            else:
                return DefinitionError(str(errors[0][0]))
        result = [header, '\n']
        for e in errors:
            if len(e[1]) > 0:
                indent = '  '
                result.append(e[1])
                result.append(':\n')
                for line in str(e[0]).split('\n'):
                    if len(line) == 0:
                        continue
                    result.append(indent)
                    result.append(line)
                    result.append('\n')
            else:
                result.append(str(e[0]))
        return DefinitionError(''.join(result))

    @property
    def language(self) -> str:
        raise NotImplementedError

    def status(self, msg: str) -> None:
        # for debugging
        indicator = '-' * self.pos + '^'
        print("%s\n%s\n%s" % (msg, self.definition, indicator))

    def fail(self, msg: str) -> None:
        errors = []
        indicator = '-' * self.pos + '^'
        exMain = DefinitionError(
            'Invalid %s declaration: %s [error at %d]\n  %s\n  %s' %
            (self.language, msg, self.pos, self.definition, indicator))
        errors.append((exMain, "Main error"))
        for err in self.otherErrors:
            errors.append((err, "Potential other error"))
        self.otherErrors = []
        raise self._make_multi_error(errors, '')

    def warn(self, msg: str) -> None:
        logger.warning(msg, location=self.location)

    def match(self, regex: Pattern) -> bool:
        match = regex.match(self.definition, self.pos)
        if match is not None:
            self._previous_state = (self.pos, self.last_match)
            self.pos = match.end()
            self.last_match = match
            return True
        return False

    def skip_string(self, string: str) -> bool:
        strlen = len(string)
        if self.definition[self.pos:self.pos + strlen] == string:
            self.pos += strlen
            return True
        return False

    def skip_word(self, word: str) -> bool:
        return self.match(re.compile(r'\b%s\b' % re.escape(word)))

    def skip_ws(self) -> bool:
        return self.match(_whitespace_re)

    def skip_word_and_ws(self, word: str) -> bool:
        if self.skip_word(word):
            self.skip_ws()
            return True
        return False

    def skip_string_and_ws(self, string: str) -> bool:
        if self.skip_string(string):
            self.skip_ws()
            return True
        return False

    @property
    def eof(self) -> bool:
        return self.pos >= self.end

    @property
    def current_char(self) -> str:
        try:
            return self.definition[self.pos]
        except IndexError:
            return 'EOF'

    @property
    def matched_text(self) -> str:
        if self.last_match is not None:
            return self.last_match.group()
        else:
            return None

    def read_rest(self) -> str:
        rv = self.definition[self.pos:]
        self.pos = self.end
        return rv

    def assert_end(self, *, allowSemicolon: bool = False) -> None:
        self.skip_ws()
        if allowSemicolon:
            if not self.eof and self.definition[self.pos:] != ';':
                self.fail('Expected end of definition or ;.')
        else:
            if not self.eof:
                self.fail('Expected end of definition.')

    ################################################################################

    @property
    def id_attributes(self):
        raise NotImplementedError

    @property
    def paren_attributes(self):
        raise NotImplementedError

    def _parse_balanced_token_seq(self, end: List[str]) -> str:
        # TODO: add handling of string literals and similar
        brackets = {'(': ')', '[': ']', '{': '}'}
        startPos = self.pos
        symbols: List[str] = []
        while not self.eof:
            if len(symbols) == 0 and self.current_char in end:
                break
            if self.current_char in brackets.keys():
                symbols.append(brackets[self.current_char])
            elif len(symbols) > 0 and self.current_char == symbols[-1]:
                symbols.pop()
            elif self.current_char in ")]}":
                self.fail("Unexpected '%s' in balanced-token-seq." % self.current_char)
            self.pos += 1
        if self.eof:
            self.fail("Could not find end of balanced-token-seq starting at %d."
                      % startPos)
        return self.definition[startPos:self.pos]

    def _parse_attribute(self) -> Optional[ASTAttribute]:
        self.skip_ws()
        # try C++11 style
        startPos = self.pos
        if self.skip_string_and_ws('['):
            if not self.skip_string('['):
                self.pos = startPos
            else:
                # TODO: actually implement the correct grammar
                arg = self._parse_balanced_token_seq(end=[']'])
                if not self.skip_string_and_ws(']'):
                    self.fail("Expected ']' in end of attribute.")
                if not self.skip_string_and_ws(']'):
                    self.fail("Expected ']' in end of attribute after [[...]")
                return ASTCPPAttribute(arg)

        # try GNU style
        if self.skip_word_and_ws('__attribute__'):
            if not self.skip_string_and_ws('('):
                self.fail("Expected '(' after '__attribute__'.")
            if not self.skip_string_and_ws('('):
                self.fail("Expected '(' after '__attribute__('.")
            attrs = []
            while 1:
                if self.match(identifier_re):
                    name = self.matched_text
                    exprs = self._parse_paren_expression_list()
                    attrs.append(ASTGnuAttribute(name, exprs))
                if self.skip_string_and_ws(','):
                    continue
                elif self.skip_string_and_ws(')'):
                    break
                else:
                    self.fail("Expected identifier, ')', or ',' in __attribute__.")
            if not self.skip_string_and_ws(')'):
                self.fail("Expected ')' after '__attribute__((...)'")
            return ASTGnuAttributeList(attrs)

        # try the simple id attributes defined by the user
        for id in self.id_attributes:
            if self.skip_word_and_ws(id):
                return ASTIdAttribute(id)

        # try the paren attributes defined by the user
        for id in self.paren_attributes:
            if not self.skip_string_and_ws(id):
                continue
            if not self.skip_string('('):
                self.fail("Expected '(' after user-defined paren-attribute.")
            arg = self._parse_balanced_token_seq(end=[')'])
            if not self.skip_string(')'):
                self.fail("Expected ')' to end user-defined paren-attribute.")
            return ASTParenAttribute(id, arg)

        return None

    def _parse_paren_expression_list(self) -> ASTBaseParenExprList:
        raise NotImplementedError
