"""
The Token class, interchangeable with ``pygments.token``.

A `Token` has some semantics for a piece of text that is given a style through
a :class:`~prompt_toolkit.styles.Style` class. A pygments lexer for instance,
returns a list of (Token, text) tuples. Each fragment of text has a token
assigned, which when combined with a style sheet, will determine the fine
style.
"""

# If we don't need any lexers or style classes from Pygments, we don't want
# Pygments to be installed for only the following 10 lines of code. So, there
# is some duplication, but this should stay compatible with Pygments.

__all__ = (
    'Token',
    'ZeroWidthEscape',
)


class _TokenType(tuple):
    def __getattr__(self, val):
        if not val or not val[0].isupper():
            return tuple.__getattribute__(self, val)

        new = _TokenType(self + (val,))
        setattr(self, val, new)
        return new

    def __repr__(self):
        return 'Token' + (self and '.' or '') + '.'.join(self)


# Prefer the Token class from Pygments. If Pygments is not installed, use our
# minimalistic Token class.
try:
    from pygments.token import Token
except ImportError:
    Token = _TokenType()


# Built-in tokens:

#: `ZeroWidthEscape` can be used for raw VT escape sequences that don't
#: cause the cursor position to move. (E.g. FinalTerm's escape sequences
#: for shell integration.)
ZeroWidthEscape = Token.ZeroWidthEscape
