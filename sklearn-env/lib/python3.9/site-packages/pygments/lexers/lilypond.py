"""
    pygments.lexers.lilypond
    ~~~~~~~~~~~~~~~~~~~~~~~~

    Lexer for LilyPond.

    :copyright: Copyright 2006-2021 by the Pygments team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re

from pygments.lexer import default, inherit, words
from pygments.lexers.lisp import SchemeLexer
from pygments.lexers._lilypond_builtins import (
    keywords, pitch_language_names, clefs, scales, repeat_types, units,
    chord_modifiers, pitches, music_functions, dynamics, articulations,
    music_commands, markup_commands, grobs, translators, contexts,
    context_properties, grob_properties, scheme_functions, paper_variables,
    header_variables
)
from pygments.token import Token

__all__ = ["LilyPondLexer"]

NAME_END_RE = r"(?=\d|[^\w\-]|[\-_][\W\d])"

def builtin_words(names, backslash, suffix=NAME_END_RE):
    prefix = r"[\-_^]?"
    if backslash == "mandatory":
        prefix += r"\\"
    elif backslash == "optional":
        prefix += r"\\?"
    else:
        assert backslash == "disallowed"
    return words(names, prefix, suffix)

class LilyPondLexer(SchemeLexer):
    """
    Lexer for input to `LilyPond <https://lilypond.org>`_, a text-based music typesetter.

    .. important::

       This lexer is meant to be used in conjunction with the ``lilypond`` style.

    .. versionadded:: 2.11
    """
    name = 'LilyPond'
    aliases = ['lilypond']
    filenames = ['*.ly']
    mimetypes = []

    flags = re.DOTALL | re.MULTILINE

    # Because parsing LilyPond input is very tricky (and in fact
    # impossible without executing LilyPond when there is Scheme
    # code in the file), this lexer does not try to recognize
    # lexical modes. Instead, it catches the most frequent pieces
    # of syntax, and, above all, knows about many kinds of builtins.

    # In order to parse embedded Scheme, this lexer subclasses the SchemeLexer.
    # It redefines the 'root' state entirely, and adds a rule for #{ #}
    # to the 'value' state. The latter is used to parse a Scheme expression
    # after #.

    def get_tokens_unprocessed(self, text):
        """Highlight Scheme variables as LilyPond builtins when applicable."""
        for index, token, value in super().get_tokens_unprocessed(text):
            if token is Token.Name.Function or token is Token.Name.Variable:
                if value in scheme_functions:
                    token = Token.Name.Builtin.SchemeFunction
            elif token is Token.Name.Builtin:
                token = Token.Name.Builtin.SchemeBuiltin
            yield index, token, value

    tokens = {
        "root": [
            # Whitespace.
            (r"\s+", Token.Whitespace),

            # Multi-line comment. These are non-nestable.
            (r"%\{.*?%\}", Token.Comment.Multiline),

            # Simple comment.
            (r"%.*?$", Token.Comment.Single),

            # End of embedded LilyPond in Scheme.
            (r"#\}", Token.Punctuation, "#pop"),

            # Embedded Scheme, starting with # ("delayed"),
            # or $ (immediate). #@ and and $@ are the lesser known
            # "list splicing operators".
            (r"[#$]@?", Token.Punctuation, "value"),

            # Any kind of punctuation:
            # - sequential music: { },
            # - parallel music: << >>,
            # - voice separator: << \\ >>,
            # - chord: < >,
            # - bar check: |,
            # - dot in nested properties: \revert NoteHead.color,
            # - equals sign in assignemnts and lists for various commands:
            #   \override Stem.color = red,
            # - comma as alternative syntax for lists: \time 3,3,2 4/4,
            # - colon in tremolos: c:32,
            # - double hyphen in lyrics: li -- ly -- pond,
            (r"\\\\|--|[{}<>=.,:|]", Token.Punctuation),

            # Pitch, with optional octavation marks, octave check,
            # and forced or cautionary accidental.
            (words(pitches, suffix=r"=?[',]*!?\??" + NAME_END_RE), Token.Pitch),

            # String, optionally with direction specifier.
            (r'[\-_^]?"', Token.String, "string"),

            # Numbers.
            (r"-?\d+\.\d+", Token.Number.Float), # 5. and .5 are not allowed
            (r"-?\d+/\d+", Token.Number.Fraction),
            # Integer, or duration with optional augmentation dots. We have no
            # way to distinguish these, so we highlight them all as numbers.
            (r"-?(\d+|\\longa|\\breve)\.*", Token.Number),
            # Separates duration and duration multiplier highlighted as fraction.
            (r"\*", Token.Number),

            # Ties, slurs, manual beams.
            (r"[~()[\]]", Token.Name.Builtin.Articulation),

            # Predefined articulation shortcuts. A direction specifier is
            # required here.
            (r"[\-_^][>^_!.\-+]", Token.Name.Builtin.Articulation),

            # Fingering numbers, string numbers.
            (r"[\-_^]?\\?\d+", Token.Name.Builtin.Articulation),

            # Builtins.
            (builtin_words(keywords, "mandatory"), Token.Keyword),
            (builtin_words(pitch_language_names, "disallowed"), Token.Name.PitchLanguage),
            (builtin_words(clefs, "disallowed"), Token.Name.Builtin.Clef),
            (builtin_words(scales, "mandatory"), Token.Name.Builtin.Scale),
            (builtin_words(repeat_types, "disallowed"), Token.Name.Builtin.RepeatType),
            (builtin_words(units, "mandatory"), Token.Number),
            (builtin_words(chord_modifiers, "disallowed"), Token.ChordModifier),
            (builtin_words(music_functions, "mandatory"), Token.Name.Builtin.MusicFunction),
            (builtin_words(dynamics, "mandatory"), Token.Name.Builtin.Dynamic),
            # Those like slurs that don't take a backslash are covered above.
            (builtin_words(articulations, "mandatory"), Token.Name.Builtin.Articulation),
            (builtin_words(music_commands, "mandatory"), Token.Name.Builtin.MusicCommand),
            (builtin_words(markup_commands, "mandatory"), Token.Name.Builtin.MarkupCommand),
            (builtin_words(grobs, "disallowed"), Token.Name.Builtin.Grob),
            (builtin_words(translators, "disallowed"), Token.Name.Builtin.Translator),
            # Optional backslash because of \layout { \context { \Score ... } }.
            (builtin_words(contexts, "optional"), Token.Name.Builtin.Context),
            (builtin_words(context_properties, "disallowed"), Token.Name.Builtin.ContextProperty),
            (builtin_words(grob_properties, "disallowed"),
             Token.Name.Builtin.GrobProperty,
             "maybe-subproperties"),
            # Optional backslashes here because output definitions are wrappers
            # around modules.  Concretely, you can do, e.g.,
            # \paper { oddHeaderMarkup = \evenHeaderMarkup }
            (builtin_words(paper_variables, "optional"), Token.Name.Builtin.PaperVariable),
            (builtin_words(header_variables, "optional"), Token.Name.Builtin.HeaderVariable),

            # Other backslashed-escaped names (like dereferencing a
            # music variable), possibly with a direction specifier.
            (r"[\-_^]?\\.+?" + NAME_END_RE, Token.Name.BackslashReference),

            # Definition of a variable. Support assignments to alist keys
            # (myAlist.my-key.my-nested-key = \markup \spam \eggs).
            (r"([^\W\d]|-)+(?=([^\W\d]|[\-.])*\s*=)", Token.Name.Lvalue),

            # Virtually everything can appear in markup mode, so we highlight
            # as text.
            (r".", Token.Text),
        ],
        "string": [
            (r'"', Token.String, "#pop"),
            (r'\\.', Token.String.Escape),
            (r'[^\\"]+', Token.String),
        ],
        "value": [
            # Scan a LilyPond value, then pop back since we had a
            # complete expression.
            (r"#\{", Token.Punctuation, ("#pop", "root")),
            inherit,
        ],
        # Grob subproperties are undeclared and it would be tedious
        # to maintain them by hand. Instead, this state allows recognizing
        # everything that looks like a-known-property.foo.bar-baz as
        # one single property name.
        "maybe-subproperties": [
            (r"\.", Token.Punctuation),
            (r"\s+", Token.Whitespace),
            (r"([^\W\d])+" + NAME_END_RE, Token.Name.Builtin.GrobProperty),
            default("#pop"),
        ]
    }
