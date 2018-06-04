"""
Processors are little transformation blocks that transform the token list from
a buffer before the BufferControl will render it to the screen.

They can insert tokens before or after, or highlight fragments by replacing the
token types.
"""
from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass
from six.moves import range

from prompt_toolkit.cache import SimpleCache
from prompt_toolkit.document import Document
from prompt_toolkit.enums import SEARCH_BUFFER
from prompt_toolkit.filters import to_cli_filter, ViInsertMultipleMode
from prompt_toolkit.layout.utils import token_list_to_text
from prompt_toolkit.reactive import Integer
from prompt_toolkit.token import Token

from .utils import token_list_len, explode_tokens

import re

__all__ = (
    'Processor',
    'Transformation',

    'HighlightSearchProcessor',
    'HighlightSelectionProcessor',
    'PasswordProcessor',
    'HighlightMatchingBracketProcessor',
    'DisplayMultipleCursors',
    'BeforeInput',
    'AfterInput',
    'AppendAutoSuggestion',
    'ConditionalProcessor',
    'ShowLeadingWhiteSpaceProcessor',
    'ShowTrailingWhiteSpaceProcessor',
    'TabsProcessor',
)


class Processor(with_metaclass(ABCMeta, object)):
    """
    Manipulate the tokens for a given line in a
    :class:`~prompt_toolkit.layout.controls.BufferControl`.
    """
    @abstractmethod
    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        """
        Apply transformation.  Returns a :class:`.Transformation` instance.

        :param cli: :class:`.CommandLineInterface` instance.
        :param lineno: The number of the line to which we apply the processor.
        :param source_to_display: A function that returns the position in the
            `tokens` for any position in the source string. (This takes
            previous processors into account.)
        :param tokens: List of tokens that we can transform. (Received from the
            previous processor.)
        """
        return Transformation(tokens)

    def has_focus(self, cli):
        """
        Processors can override the focus.
        (Used for the reverse-i-search prefix in DefaultPrompt.)
        """
        return False


class Transformation(object):
    """
    Transformation result, as returned by :meth:`.Processor.apply_transformation`.

    Important: Always make sure that the length of `document.text` is equal to
               the length of all the text in `tokens`!

    :param tokens: The transformed tokens. To be displayed, or to pass to the
        next processor.
    :param source_to_display: Cursor position transformation from original string to
        transformed string.
    :param display_to_source: Cursor position transformed from source string to
        original string.
    """
    def __init__(self, tokens, source_to_display=None, display_to_source=None):
        self.tokens = tokens
        self.source_to_display = source_to_display or (lambda i: i)
        self.display_to_source = display_to_source or (lambda i: i)


class HighlightSearchProcessor(Processor):
    """
    Processor that highlights search matches in the document.
    Note that this doesn't support multiline search matches yet.

    :param preview_search: A Filter; when active it indicates that we take
        the search text in real time while the user is typing, instead of the
        last active search state.
    """
    def __init__(self, preview_search=False, search_buffer_name=SEARCH_BUFFER,
                 get_search_state=None):
        self.preview_search = to_cli_filter(preview_search)
        self.search_buffer_name = search_buffer_name
        self.get_search_state = get_search_state or (lambda cli: cli.search_state)

    def _get_search_text(self, cli):
        """
        The text we are searching for.
        """
        # When the search buffer has focus, take that text.
        if self.preview_search(cli) and cli.buffers[self.search_buffer_name].text:
            return cli.buffers[self.search_buffer_name].text
        # Otherwise, take the text of the last active search.
        else:
            return self.get_search_state(cli).text

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        search_text = self._get_search_text(cli)
        searchmatch_current_token = (':', ) + Token.SearchMatch.Current
        searchmatch_token = (':', ) + Token.SearchMatch

        if search_text and not cli.is_returning:
            # For each search match, replace the Token.
            line_text = token_list_to_text(tokens)
            tokens = explode_tokens(tokens)

            flags = re.IGNORECASE if cli.is_ignoring_case else 0

            # Get cursor column.
            if document.cursor_position_row == lineno:
                cursor_column = source_to_display(document.cursor_position_col)
            else:
                cursor_column = None

            for match in re.finditer(re.escape(search_text), line_text, flags=flags):
                if cursor_column is not None:
                    on_cursor = match.start() <= cursor_column < match.end()
                else:
                    on_cursor = False

                for i in range(match.start(), match.end()):
                    old_token, text = tokens[i]
                    if on_cursor:
                        tokens[i] = (old_token + searchmatch_current_token, tokens[i][1])
                    else:
                        tokens[i] = (old_token + searchmatch_token, tokens[i][1])

        return Transformation(tokens)


class HighlightSelectionProcessor(Processor):
    """
    Processor that highlights the selection in the document.
    """
    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        selected_token = (':', ) + Token.SelectedText

        # In case of selection, highlight all matches.
        selection_at_line = document.selection_range_at_line(lineno)

        if selection_at_line:
            from_, to = selection_at_line
            from_ = source_to_display(from_)
            to = source_to_display(to)

            tokens = explode_tokens(tokens)

            if from_ == 0 and to == 0 and len(tokens) == 0:
                # When this is an empty line, insert a space in order to
                # visualiase the selection.
                return Transformation([(Token.SelectedText, ' ')])
            else:
                for i in range(from_, to + 1):
                    if i < len(tokens):
                        old_token, old_text = tokens[i]
                        tokens[i] = (old_token + selected_token, old_text)

        return Transformation(tokens)


class PasswordProcessor(Processor):
    """
    Processor that turns masks the input. (For passwords.)

    :param char: (string) Character to be used. "*" by default.
    """
    def __init__(self, char='*'):
        self.char = char

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        tokens = [(token, self.char * len(text)) for token, text in tokens]
        return Transformation(tokens)


class HighlightMatchingBracketProcessor(Processor):
    """
    When the cursor is on or right after a bracket, it highlights the matching
    bracket.

    :param max_cursor_distance: Only highlight matching brackets when the
        cursor is within this distance. (From inside a `Processor`, we can't
        know which lines will be visible on the screen. But we also don't want
        to scan the whole document for matching brackets on each key press, so
        we limit to this value.)
    """
    _closing_braces = '])}>'

    def __init__(self, chars='[](){}<>', max_cursor_distance=1000):
        self.chars = chars
        self.max_cursor_distance = max_cursor_distance

        self._positions_cache = SimpleCache(maxsize=8)

    def _get_positions_to_highlight(self, document):
        """
        Return a list of (row, col) tuples that need to be highlighted.
        """
        # Try for the character under the cursor.
        if document.current_char and document.current_char in self.chars:
            pos = document.find_matching_bracket_position(
                    start_pos=document.cursor_position - self.max_cursor_distance,
                    end_pos=document.cursor_position + self.max_cursor_distance)

        # Try for the character before the cursor.
        elif (document.char_before_cursor and document.char_before_cursor in
              self._closing_braces and document.char_before_cursor in self.chars):
            document = Document(document.text, document.cursor_position - 1)

            pos = document.find_matching_bracket_position(
                    start_pos=document.cursor_position - self.max_cursor_distance,
                    end_pos=document.cursor_position + self.max_cursor_distance)
        else:
            pos = None

        # Return a list of (row, col) tuples that need to be highlighted.
        if pos:
            pos += document.cursor_position  # pos is relative.
            row, col = document.translate_index_to_position(pos)
            return [(row, col), (document.cursor_position_row, document.cursor_position_col)]
        else:
            return []

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Get the highlight positions.
        key = (cli.render_counter, document.text, document.cursor_position)
        positions = self._positions_cache.get(
            key, lambda: self._get_positions_to_highlight(document))

        # Apply if positions were found at this line.
        if positions:
            for row, col in positions:
                if row == lineno:
                    col = source_to_display(col)
                    tokens = explode_tokens(tokens)
                    token, text = tokens[col]

                    if col == document.cursor_position_col:
                        token += (':', ) + Token.MatchingBracket.Cursor
                    else:
                        token += (':', ) + Token.MatchingBracket.Other

                    tokens[col] = (token, text)

        return Transformation(tokens)


class DisplayMultipleCursors(Processor):
    """
    When we're in Vi block insert mode, display all the cursors.
    """
    _insert_multiple =  ViInsertMultipleMode()

    def __init__(self, buffer_name):
        self.buffer_name = buffer_name

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        buff = cli.buffers[self.buffer_name]

        if self._insert_multiple(cli):
            positions = buff.multiple_cursor_positions
            tokens = explode_tokens(tokens)

            # If any cursor appears on the current line, highlight that.
            start_pos = document.translate_row_col_to_index(lineno, 0)
            end_pos = start_pos + len(document.lines[lineno])

            token_suffix = (':', ) + Token.MultipleCursors.Cursor

            for p in positions:
                if start_pos <= p < end_pos:
                    column = source_to_display(p - start_pos)

                    # Replace token.
                    token, text = tokens[column]
                    token += token_suffix
                    tokens[column] = (token, text)
                elif p == end_pos:
                    tokens.append((token_suffix, ' '))

            return Transformation(tokens)
        else:
            return Transformation(tokens)


class BeforeInput(Processor):
    """
    Insert tokens before the input.

    :param get_tokens: Callable that takes a
        :class:`~prompt_toolkit.interface.CommandLineInterface` and returns the
        list of tokens to be inserted.
    """
    def __init__(self, get_tokens):
        assert callable(get_tokens)
        self.get_tokens = get_tokens

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        if lineno == 0:
            tokens_before = self.get_tokens(cli)
            tokens = tokens_before + tokens

            shift_position = token_list_len(tokens_before)
            source_to_display = lambda i: i + shift_position
            display_to_source = lambda i: i - shift_position
        else:
            source_to_display = None
            display_to_source = None

        return Transformation(tokens, source_to_display=source_to_display,
                              display_to_source=display_to_source)

    @classmethod
    def static(cls, text, token=Token):
        """
        Create a :class:`.BeforeInput` instance that always inserts the same
        text.
        """
        def get_static_tokens(cli):
            return [(token, text)]
        return cls(get_static_tokens)

    def __repr__(self):
        return '%s(get_tokens=%r)' % (
            self.__class__.__name__, self.get_tokens)


class AfterInput(Processor):
    """
    Insert tokens after the input.

    :param get_tokens: Callable that takes a
        :class:`~prompt_toolkit.interface.CommandLineInterface` and returns the
        list of tokens to be appended.
    """
    def __init__(self, get_tokens):
        assert callable(get_tokens)
        self.get_tokens = get_tokens

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Insert tokens after the last line.
        if lineno == document.line_count - 1:
            return Transformation(tokens=tokens + self.get_tokens(cli))
        else:
            return Transformation(tokens=tokens)

    @classmethod
    def static(cls, text, token=Token):
        """
        Create a :class:`.AfterInput` instance that always inserts the same
        text.
        """
        def get_static_tokens(cli):
            return [(token, text)]
        return cls(get_static_tokens)

    def __repr__(self):
        return '%s(get_tokens=%r)' % (
            self.__class__.__name__, self.get_tokens)


class AppendAutoSuggestion(Processor):
    """
    Append the auto suggestion to the input.
    (The user can then press the right arrow the insert the suggestion.)

    :param buffer_name: The name of the buffer from where we should take the
        auto suggestion. If not given, we take the current buffer.
    """
    def __init__(self, buffer_name=None, token=Token.AutoSuggestion):
        self.buffer_name = buffer_name
        self.token = token

    def _get_buffer(self, cli):
        if self.buffer_name:
            return cli.buffers[self.buffer_name]
        else:
            return cli.current_buffer

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Insert tokens after the last line.
        if lineno == document.line_count - 1:
            buffer = self._get_buffer(cli)

            if buffer.suggestion and buffer.document.is_cursor_at_the_end:
                suggestion = buffer.suggestion.text
            else:
                suggestion = ''

            return Transformation(tokens=tokens + [(self.token, suggestion)])
        else:
            return Transformation(tokens=tokens)


class ShowLeadingWhiteSpaceProcessor(Processor):
    """
    Make leading whitespace visible.

    :param get_char: Callable that takes a :class:`CommandLineInterface`
        instance and returns one character.
    :param token: Token to be used.
    """
    def __init__(self, get_char=None, token=Token.LeadingWhiteSpace):
        assert get_char is None or callable(get_char)

        if get_char is None:
            def get_char(cli):
                if '\xb7'.encode(cli.output.encoding(), 'replace') == b'?':
                    return '.'
                else:
                    return '\xb7'

        self.token = token
        self.get_char = get_char

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Walk through all te tokens.
        if tokens and token_list_to_text(tokens).startswith(' '):
            t = (self.token, self.get_char(cli))
            tokens = explode_tokens(tokens)

            for i in range(len(tokens)):
                if tokens[i][1] == ' ':
                    tokens[i] = t
                else:
                    break

        return Transformation(tokens)


class ShowTrailingWhiteSpaceProcessor(Processor):
    """
    Make trailing whitespace visible.

    :param get_char: Callable that takes a :class:`CommandLineInterface`
        instance and returns one character.
    :param token: Token to be used.
    """
    def __init__(self, get_char=None, token=Token.TrailingWhiteSpace):
        assert get_char is None or callable(get_char)

        if get_char is None:
            def get_char(cli):
                if '\xb7'.encode(cli.output.encoding(), 'replace') == b'?':
                    return '.'
                else:
                    return '\xb7'

        self.token = token
        self.get_char = get_char


    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        if tokens and tokens[-1][1].endswith(' '):
            t = (self.token, self.get_char(cli))
            tokens = explode_tokens(tokens)

            # Walk backwards through all te tokens and replace whitespace.
            for i in range(len(tokens) - 1, -1, -1):
                char = tokens[i][1]
                if char == ' ':
                    tokens[i] = t
                else:
                    break

        return Transformation(tokens)


class TabsProcessor(Processor):
    """
    Render tabs as spaces (instead of ^I) or make them visible (for instance,
    by replacing them with dots.)

    :param tabstop: (Integer) Horizontal space taken by a tab.
    :param get_char1: Callable that takes a `CommandLineInterface` and return a
        character (text of length one). This one is used for the first space
        taken by the tab.
    :param get_char2: Like `get_char1`, but for the rest of the space.
    """
    def __init__(self, tabstop=4, get_char1=None, get_char2=None, token=Token.Tab):
        assert isinstance(tabstop, Integer)
        assert get_char1 is None or callable(get_char1)
        assert get_char2 is None or callable(get_char2)

        self.get_char1 = get_char1 or get_char2 or (lambda cli: '|')
        self.get_char2 = get_char2 or get_char1 or (lambda cli: '\u2508')
        self.tabstop = tabstop
        self.token = token

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        tabstop = int(self.tabstop)
        token = self.token

        # Create separator for tabs.
        separator1 = self.get_char1(cli)
        separator2 = self.get_char2(cli)

        # Transform tokens.
        tokens = explode_tokens(tokens)

        position_mappings = {}
        result_tokens = []
        pos = 0

        for i, token_and_text in enumerate(tokens):
            position_mappings[i] = pos

            if token_and_text[1] == '\t':
                # Calculate how many characters we have to insert.
                count = tabstop - (pos % tabstop)
                if count == 0:
                    count = tabstop

                # Insert tab.
                result_tokens.append((token, separator1))
                result_tokens.append((token, separator2 * (count - 1)))
                pos += count
            else:
                result_tokens.append(token_and_text)
                pos += 1

        position_mappings[len(tokens)] = pos

        def source_to_display(from_position):
            " Maps original cursor position to the new one. "
            return position_mappings[from_position]

        def display_to_source(display_pos):
            " Maps display cursor position to the original one. "
            position_mappings_reversed = dict((v, k) for k, v in position_mappings.items())

            while display_pos >= 0:
                try:
                    return position_mappings_reversed[display_pos]
                except KeyError:
                    display_pos -= 1
            return 0

        return Transformation(
            result_tokens,
            source_to_display=source_to_display,
            display_to_source=display_to_source)


class ConditionalProcessor(Processor):
    """
    Processor that applies another processor, according to a certain condition.
    Example::

        # Create a function that returns whether or not the processor should
        # currently be applied.
        def highlight_enabled(cli):
            return true_or_false

        # Wrapt it in a `ConditionalProcessor` for usage in a `BufferControl`.
        BufferControl(input_processors=[
            ConditionalProcessor(HighlightSearchProcessor(),
                                 Condition(highlight_enabled))])

    :param processor: :class:`.Processor` instance.
    :param filter: :class:`~prompt_toolkit.filters.CLIFilter` instance.
    """
    def __init__(self, processor, filter):
        assert isinstance(processor, Processor)

        self.processor = processor
        self.filter = to_cli_filter(filter)

    def apply_transformation(self, cli, document, lineno, source_to_display, tokens):
        # Run processor when enabled.
        if self.filter(cli):
            return self.processor.apply_transformation(
                cli, document, lineno, source_to_display, tokens)
        else:
            return Transformation(tokens)

    def has_focus(self, cli):
        if self.filter(cli):
            return self.processor.has_focus(cli)
        else:
            return False

    def __repr__(self):
        return '%s(processor=%r, filter=%r)' % (
            self.__class__.__name__, self.processor, self.filter)
