"""
The default styling.
"""
from __future__ import unicode_literals

from prompt_toolkit.token import Token

__all__ = (
    'DEFAULT_STYLE_EXTENSIONS',
    'default_style_extensions',
)


#: Styling of prompt-toolkit specific tokens, that are not know by the default
#: Pygments style.
DEFAULT_STYLE_EXTENSIONS = {
    # Highlighting of search matches in document.
    Token.SearchMatch:                            'noinherit reverse',
    Token.SearchMatch.Current:                    'noinherit #ffffff bg:#448844 underline',

    # Highlighting of select text in document.
    Token.SelectedText:                           'reverse',

    Token.CursorColumn:                           'bg:#dddddd',
    Token.CursorLine:                             'underline',
    Token.ColorColumn:                            'bg:#ccaacc',

    # Highlighting of matching brackets.
    Token.MatchingBracket:                        '',
    Token.MatchingBracket.Other:                  '#000000 bg:#aacccc',
    Token.MatchingBracket.Cursor:                 '#ff8888 bg:#880000',

    Token.MultipleCursors.Cursor:                 '#000000 bg:#ccccaa',

    # Line numbers.
    Token.LineNumber:                             '#888888',
    Token.LineNumber.Current:                     'bold',
    Token.Tilde:                                  '#8888ff',

    # Default prompt.
    Token.Prompt:                                 '',
    Token.Prompt.Arg:                             'noinherit',
    Token.Prompt.Search:                          'noinherit',
    Token.Prompt.Search.Text:                     '',

    # Search toolbar.
    Token.Toolbar.Search:                         'bold',
    Token.Toolbar.Search.Text:                    'nobold',

    # System toolbar
    Token.Toolbar.System:                         'bold',
    Token.Toolbar.System.Text:                    'nobold',

    # "arg" toolbar.
    Token.Toolbar.Arg:                            'bold',
    Token.Toolbar.Arg.Text:                       'nobold',

    # Validation toolbar.
    Token.Toolbar.Validation:                     'bg:#550000 #ffffff',
    Token.WindowTooSmall:                         'bg:#550000 #ffffff',

    # Completions toolbar.
    Token.Toolbar.Completions:                    'bg:#bbbbbb #000000',
    Token.Toolbar.Completions.Arrow:              'bg:#bbbbbb #000000 bold',
    Token.Toolbar.Completions.Completion:         'bg:#bbbbbb #000000',
    Token.Toolbar.Completions.Completion.Current: 'bg:#444444 #ffffff',

    # Completions menu.
    Token.Menu.Completions:                       'bg:#bbbbbb #000000',
    Token.Menu.Completions.Completion:            '',
    Token.Menu.Completions.Completion.Current:    'bg:#888888 #ffffff',
    Token.Menu.Completions.Meta:                  'bg:#999999 #000000',
    Token.Menu.Completions.Meta.Current:          'bg:#aaaaaa #000000',
    Token.Menu.Completions.MultiColumnMeta:       'bg:#aaaaaa #000000',

    # Scrollbars.
    Token.Scrollbar:                              'bg:#888888',
    Token.Scrollbar.Button:                       'bg:#444444',
    Token.Scrollbar.Arrow:                        'bg:#222222 #888888 bold',

    # Auto suggestion text.
    Token.AutoSuggestion:                         '#666666',

    # Trailing whitespace and tabs.
    Token.TrailingWhiteSpace:                     '#999999',
    Token.Tab:                                    '#999999',

    # When Control-C has been pressed. Grayed.
    Token.Aborted:                                '#888888',

    # Entering a Vi digraph.
    Token.Digraph:                                '#4444ff',
}

default_style_extensions = DEFAULT_STYLE_EXTENSIONS  # Old name.
