"""
    pygments.lexers.rita
    ~~~~~~~~~~~~~~~~~~~~

    Lexers for RITA language

    :copyright: Copyright 2006-2021 by the Pygments team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re

from pygments.lexer import RegexLexer, include, bygroups, using, this, \
    inherit, words
from pygments.token import Comment, Operator, Keyword, Name, Literal, Punctuation, Text, Whitespace

__all__ = ['RitaLexer']


class RitaLexer(RegexLexer):
    """
    Lexer for `RITA <https://github.com/zaibacu/rita-dsl>`_

     .. versionadded:: 2.11
    """
    name = 'Rita'
    filenames = ['*.rita']
    aliases = ['rita']
    mimetypes = ['text/rita']

    tokens = {
        'root': [
            (r'\n', Whitespace),
            (r'\s+', Whitespace),
            (r'#(.*?)\n', Comment.Single),
            (r'@(.*?)\n', Operator),  # Yes, whole line as an operator
            (r'"(\w|\d|\s|(\\")|[\'_\-./,\?\!])+?"', Literal),
            (r'\'(\w|\d|\s|(\\\')|["_\-./,\?\!])+?\'', Literal),
            (r'([A-Z_]+)', Keyword),
            (r'([a-z0-9_]+)', Name),
            (r'((->)|[!?+*|=])', Operator),
            (r'[\(\),\{\}]', Punctuation)
        ]
    }
