"""
`GrammarLexer` is compatible with Pygments lexers and can be used to highlight
the input using a regular grammar with token annotations.
"""
from __future__ import unicode_literals
from prompt_toolkit.document import Document
from prompt_toolkit.layout.lexers import Lexer
from prompt_toolkit.layout.utils import split_lines
from prompt_toolkit.token import Token

from .compiler import _CompiledGrammar
from six.moves import range

__all__ = (
    'GrammarLexer',
)


class GrammarLexer(Lexer):
    """
    Lexer which can be used for highlighting of tokens according to variables in the grammar.

    (It does not actual lexing of the string, but it exposes an API, compatible
    with the Pygments lexer class.)

    :param compiled_grammar: Grammar as returned by the `compile()` function.
    :param lexers: Dictionary mapping variable names of the regular grammar to
                   the lexers that should be used for this part. (This can
                   call other lexers recursively.) If you wish a part of the
                   grammar to just get one token, use a
                   `prompt_toolkit.layout.lexers.SimpleLexer`.
    """
    def __init__(self, compiled_grammar, default_token=None, lexers=None):
        assert isinstance(compiled_grammar, _CompiledGrammar)
        assert default_token is None or isinstance(default_token, tuple)
        assert lexers is None or all(isinstance(v, Lexer) for k, v in lexers.items())
        assert lexers is None or isinstance(lexers, dict)

        self.compiled_grammar = compiled_grammar
        self.default_token = default_token or Token
        self.lexers = lexers or {}

    def _get_tokens(self, cli, text):
        m = self.compiled_grammar.match_prefix(text)

        if m:
            characters = [[self.default_token, c] for c in text]

            for v in m.variables():
                # If we have a `Lexer` instance for this part of the input.
                # Tokenize recursively and apply tokens.
                lexer = self.lexers.get(v.varname)

                if lexer:
                    document = Document(text[v.start:v.stop])
                    lexer_tokens_for_line = lexer.lex_document(cli, document)
                    lexer_tokens = []
                    for i in range(len(document.lines)):
                        lexer_tokens.extend(lexer_tokens_for_line(i))
                        lexer_tokens.append((Token, '\n'))
                    if lexer_tokens:
                        lexer_tokens.pop()

                    i = v.start
                    for t, s in lexer_tokens:
                        for c in s:
                            if characters[i][0] == self.default_token:
                                characters[i][0] = t
                            i += 1

            # Highlight trailing input.
            trailing_input = m.trailing_input()
            if trailing_input:
                for i in range(trailing_input.start, trailing_input.stop):
                    characters[i][0] = Token.TrailingInput

            return characters
        else:
            return [(Token, text)]

    def lex_document(self, cli, document):
        lines = list(split_lines(self._get_tokens(cli, document.text)))

        def get_line(lineno):
            try:
                return lines[lineno]
            except IndexError:
                return []

        return get_line
