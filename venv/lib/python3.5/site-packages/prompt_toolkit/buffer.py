"""
Data structures for the Buffer.
It holds the text, cursor position, history, etc...
"""
from __future__ import unicode_literals

from .auto_suggest import AutoSuggest
from .clipboard import ClipboardData
from .completion import Completer, Completion, CompleteEvent
from .document import Document
from .enums import IncrementalSearchDirection
from .filters import to_simple_filter
from .history import History, InMemoryHistory
from .search_state import SearchState
from .selection import SelectionType, SelectionState, PasteMode
from .utils import Event
from .cache import FastDictCache
from .validation import ValidationError

from six.moves import range

import os
import re
import shlex
import six
import subprocess
import tempfile

__all__ = (
    'EditReadOnlyBuffer',
    'AcceptAction',
    'Buffer',
    'indent',
    'unindent',
    'reshape_text',
)


class EditReadOnlyBuffer(Exception):
    " Attempt editing of read-only :class:`.Buffer`. "


class AcceptAction(object):
    """
    What to do when the input is accepted by the user.
    (When Enter was pressed in the command line.)

    :param handler: (optional) A callable which takes a
        :class:`~prompt_toolkit.interface.CommandLineInterface` and
        :class:`~prompt_toolkit.document.Document`. It is called when the user
        accepts input.
    """
    def __init__(self, handler=None):
        assert handler is None or callable(handler)
        self.handler = handler

    @classmethod
    def run_in_terminal(cls, handler, render_cli_done=False):
        """
        Create an :class:`.AcceptAction` that runs the given handler in the
        terminal.

        :param render_cli_done: When True, render the interface in the 'Done'
                state first, then execute the function. If False, erase the
                interface instead.
        """
        def _handler(cli, buffer):
            cli.run_in_terminal(lambda: handler(cli, buffer), render_cli_done=render_cli_done)
        return AcceptAction(handler=_handler)

    @property
    def is_returnable(self):
        """
        True when there is something handling accept.
        """
        return bool(self.handler)

    def validate_and_handle(self, cli, buffer):
        """
        Validate buffer and handle the accept action.
        """
        if buffer.validate():
            if self.handler:
                self.handler(cli, buffer)

            buffer.append_to_history()


def _return_document_handler(cli, buffer):
    # Set return value.
    cli.set_return_value(buffer.document)

    # Make sure that if we run this UI again, that we reset this buffer, next
    # time.
    def reset_this_buffer():
        buffer.reset()
    cli.pre_run_callables.append(reset_this_buffer)


AcceptAction.RETURN_DOCUMENT = AcceptAction(_return_document_handler)
AcceptAction.IGNORE = AcceptAction(handler=None)


class ValidationState(object):
    " The validation state of a buffer. This is set after the validation. "
    VALID = 'VALID'
    INVALID = 'INVALID'
    UNKNOWN = 'UNKNOWN'


class CompletionState(object):
    """
    Immutable class that contains a completion state.
    """
    def __init__(self, original_document, current_completions=None, complete_index=None):
        #: Document as it was when the completion started.
        self.original_document = original_document

        #: List of all the current Completion instances which are possible at
        #: this point.
        self.current_completions = current_completions or []

        #: Position in the `current_completions` array.
        #: This can be `None` to indicate "no completion", the original text.
        self.complete_index = complete_index  # Position in the `_completions` array.

    def __repr__(self):
        return '%s(%r, <%r> completions, index=%r)' % (
            self.__class__.__name__,
            self.original_document, len(self.current_completions), self.complete_index)

    def go_to_index(self, index):
        """
        Create a new :class:`.CompletionState` object with the new index.
        """
        return CompletionState(self.original_document, self.current_completions, complete_index=index)

    def new_text_and_position(self):
        """
        Return (new_text, new_cursor_position) for this completion.
        """
        if self.complete_index is None:
            return self.original_document.text, self.original_document.cursor_position
        else:
            original_text_before_cursor = self.original_document.text_before_cursor
            original_text_after_cursor = self.original_document.text_after_cursor

            c = self.current_completions[self.complete_index]
            if c.start_position == 0:
                before = original_text_before_cursor
            else:
                before = original_text_before_cursor[:c.start_position]

            new_text = before + c.text + original_text_after_cursor
            new_cursor_position = len(before) + len(c.text)
            return new_text, new_cursor_position

    @property
    def current_completion(self):
        """
        Return the current completion, or return `None` when no completion is
        selected.
        """
        if self.complete_index is not None:
            return self.current_completions[self.complete_index]


_QUOTED_WORDS_RE = re.compile(r"""(\s+|".*?"|'.*?')""")


class YankNthArgState(object):
    """
    For yank-last-arg/yank-nth-arg: Keep track of where we are in the history.
    """
    def __init__(self, history_position=0, n=-1, previous_inserted_word=''):
        self.history_position = history_position
        self.previous_inserted_word = previous_inserted_word
        self.n = n

    def __repr__(self):
        return '%s(history_position=%r, n=%r, previous_inserted_word=%r)' % (
            self.__class__.__name__, self.history_position, self.n,
            self.previous_inserted_word)


class Buffer(object):
    """
    The core data structure that holds the text and cursor position of the
    current input line and implements all text manupulations on top of it. It
    also implements the history, undo stack and the completion state.

    :param completer: :class:`~prompt_toolkit.completion.Completer` instance.
    :param history: :class:`~prompt_toolkit.history.History` instance.
    :param tempfile_suffix: Suffix to be appended to the tempfile for the 'open
                           in editor' function.

    Events:

    :param on_text_changed: When the buffer text changes. (Callable on None.)
    :param on_text_insert: When new text is inserted. (Callable on None.)
    :param on_cursor_position_changed: When the cursor moves. (Callable on None.)

    Filters:

    :param is_multiline: :class:`~prompt_toolkit.filters.SimpleFilter` to
        indicate whether we should consider this buffer a multiline input. If
        so, key bindings can decide to insert newlines when pressing [Enter].
        (Instead of accepting the input.)
    :param complete_while_typing: :class:`~prompt_toolkit.filters.SimpleFilter`
        instance. Decide whether or not to do asynchronous autocompleting while
        typing.
    :param enable_history_search: :class:`~prompt_toolkit.filters.SimpleFilter`
        to indicate when up-arrow partial string matching is enabled. It is
        adviced to not enable this at the same time as `complete_while_typing`,
        because when there is an autocompletion found, the up arrows usually
        browse through the completions, rather than through the history.
    :param read_only: :class:`~prompt_toolkit.filters.SimpleFilter`. When True,
        changes will not be allowed.
    """
    def __init__(self, completer=None, auto_suggest=None, history=None,
                 validator=None, tempfile_suffix='',
                 is_multiline=False, complete_while_typing=False,
                 enable_history_search=False, initial_document=None,
                 accept_action=AcceptAction.IGNORE, read_only=False,
                 on_text_changed=None, on_text_insert=None, on_cursor_position_changed=None):

        # Accept both filters and booleans as input.
        enable_history_search = to_simple_filter(enable_history_search)
        is_multiline = to_simple_filter(is_multiline)
        complete_while_typing = to_simple_filter(complete_while_typing)
        read_only = to_simple_filter(read_only)

        # Validate input.
        assert completer is None or isinstance(completer, Completer)
        assert auto_suggest is None or isinstance(auto_suggest, AutoSuggest)
        assert history is None or isinstance(history, History)
        assert on_text_changed is None or callable(on_text_changed)
        assert on_text_insert is None or callable(on_text_insert)
        assert on_cursor_position_changed is None or callable(on_cursor_position_changed)

        self.completer = completer
        self.auto_suggest = auto_suggest
        self.validator = validator
        self.tempfile_suffix = tempfile_suffix
        self.accept_action = accept_action

        # Filters. (Usually, used by the key bindings to drive the buffer.)
        self.is_multiline = is_multiline
        self.complete_while_typing = complete_while_typing
        self.enable_history_search = enable_history_search
        self.read_only = read_only

        # Text width. (For wrapping, used by the Vi 'gq' operator.)
        self.text_width = 0

        #: The command buffer history.
        # Note that we shouldn't use a lazy 'or' here. bool(history) could be
        # False when empty.
        self.history = InMemoryHistory() if history is None else history

        self.__cursor_position = 0

        # Events
        self.on_text_changed = Event(self, on_text_changed)
        self.on_text_insert = Event(self, on_text_insert)
        self.on_cursor_position_changed = Event(self, on_cursor_position_changed)

        # Document cache. (Avoid creating new Document instances.)
        self._document_cache = FastDictCache(Document, size=10)

        self.reset(initial_document=initial_document)

    def reset(self, initial_document=None, append_to_history=False):
        """
        :param append_to_history: Append current input to history first.
        """
        assert initial_document is None or isinstance(initial_document, Document)

        if append_to_history:
            self.append_to_history()

        initial_document = initial_document or Document()

        self.__cursor_position = initial_document.cursor_position

        # `ValidationError` instance. (Will be set when the input is wrong.)
        self.validation_error = None
        self.validation_state = ValidationState.UNKNOWN

        # State of the selection.
        self.selection_state = None

        # Multiple cursor mode. (When we press 'I' or 'A' in visual-block mode,
        # we can insert text on multiple lines at once. This is implemented by
        # using multiple cursors.)
        self.multiple_cursor_positions = []

        # When doing consecutive up/down movements, prefer to stay at this column.
        self.preferred_column = None

        # State of complete browser
        self.complete_state = None  # For interactive completion through Ctrl-N/Ctrl-P.

        # State of Emacs yank-nth-arg completion.
        self.yank_nth_arg_state = None  # for yank-nth-arg.

        # Remember the document that we had *right before* the last paste
        # operation. This is used for rotating through the kill ring.
        self.document_before_paste = None

        # Current suggestion.
        self.suggestion = None

        # The history search text. (Used for filtering the history when we
        # browse through it.)
        self.history_search_text = None

        # Undo/redo stacks
        self._undo_stack = []  # Stack of (text, cursor_position)
        self._redo_stack = []

        #: The working lines. Similar to history, except that this can be
        #: modified. The user can press arrow_up and edit previous entries.
        #: Ctrl-C should reset this, and copy the whole history back in here.
        #: Enter should process the current command and append to the real
        #: history.
        self._working_lines = self.history.strings[:]
        self._working_lines.append(initial_document.text)
        self.__working_index = len(self._working_lines) - 1

    # <getters/setters>

    def _set_text(self, value):
        """ set text at current working_index. Return whether it changed. """
        working_index = self.working_index
        working_lines = self._working_lines

        original_value = working_lines[working_index]
        working_lines[working_index] = value

        # Return True when this text has been changed.
        if len(value) != len(original_value):
            # For Python 2, it seems that when two strings have a different
            # length and one is a prefix of the other, Python still scans
            # character by character to see whether the strings are different.
            # (Some benchmarking showed significant differences for big
            # documents. >100,000 of lines.)
            return True
        elif value != original_value:
            return True
        return False

    def _set_cursor_position(self, value):
        """ Set cursor position. Return whether it changed. """
        original_position = self.__cursor_position
        self.__cursor_position = max(0, value)

        return value != original_position

    @property
    def text(self):
        return self._working_lines[self.working_index]

    @text.setter
    def text(self, value):
        """
        Setting text. (When doing this, make sure that the cursor_position is
        valid for this text. text/cursor_position should be consistent at any time,
        otherwise set a Document instead.)
        """
        assert isinstance(value, six.text_type), 'Got %r' % value
        assert self.cursor_position <= len(value)

        # Don't allow editing of read-only buffers.
        if self.read_only():
            raise EditReadOnlyBuffer()

        changed = self._set_text(value)

        if changed:
            self._text_changed()

            # Reset history search text.
            self.history_search_text = None

    @property
    def cursor_position(self):
        return self.__cursor_position

    @cursor_position.setter
    def cursor_position(self, value):
        """
        Setting cursor position.
        """
        assert isinstance(value, int)
        assert value <= len(self.text)

        changed = self._set_cursor_position(value)

        if changed:
            self._cursor_position_changed()

    @property
    def working_index(self):
        return self.__working_index

    @working_index.setter
    def working_index(self, value):
        if self.__working_index != value:
            self.__working_index = value
            self._text_changed()

    def _text_changed(self):
        # Remove any validation errors and complete state.
        self.validation_error = None
        self.validation_state = ValidationState.UNKNOWN
        self.complete_state = None
        self.yank_nth_arg_state = None
        self.document_before_paste = None
        self.selection_state = None
        self.suggestion = None
        self.preferred_column = None

        # fire 'on_text_changed' event.
        self.on_text_changed.fire()

    def _cursor_position_changed(self):
        # Remove any validation errors and complete state.
        self.validation_error = None
        self.validation_state = ValidationState.UNKNOWN
        self.complete_state = None
        self.yank_nth_arg_state = None
        self.document_before_paste = None

        # Unset preferred_column. (Will be set after the cursor movement, if
        # required.)
        self.preferred_column = None

        # Note that the cursor position can change if we have a selection the
        # new position of the cursor determines the end of the selection.

        # fire 'on_cursor_position_changed' event.
        self.on_cursor_position_changed.fire()

    @property
    def document(self):
        """
        Return :class:`~prompt_toolkit.document.Document` instance from the
        current text, cursor position and selection state.
        """
        return self._document_cache[
            self.text, self.cursor_position, self.selection_state]

    @document.setter
    def document(self, value):
        """
        Set :class:`~prompt_toolkit.document.Document` instance.

        This will set both the text and cursor position at the same time, but
        atomically. (Change events will be triggered only after both have been set.)
        """
        self.set_document(value)

    def set_document(self, value, bypass_readonly=False):
        """
        Set :class:`~prompt_toolkit.document.Document` instance. Like the
        ``document`` property, but accept an ``bypass_readonly`` argument.

        :param bypass_readonly: When True, don't raise an
                                :class:`.EditReadOnlyBuffer` exception, even
                                when the buffer is read-only.
        """
        assert isinstance(value, Document)

        # Don't allow editing of read-only buffers.
        if not bypass_readonly and self.read_only():
            raise EditReadOnlyBuffer()

        # Set text and cursor position first.
        text_changed = self._set_text(value.text)
        cursor_position_changed = self._set_cursor_position(value.cursor_position)

        # Now handle change events. (We do this when text/cursor position is
        # both set and consistent.)
        if text_changed:
            self._text_changed()

        if cursor_position_changed:
            self._cursor_position_changed()

    # End of <getters/setters>

    def save_to_undo_stack(self, clear_redo_stack=True):
        """
        Safe current state (input text and cursor position), so that we can
        restore it by calling undo.
        """
        # Safe if the text is different from the text at the top of the stack
        # is different. If the text is the same, just update the cursor position.
        if self._undo_stack and self._undo_stack[-1][0] == self.text:
            self._undo_stack[-1] = (self._undo_stack[-1][0], self.cursor_position)
        else:
            self._undo_stack.append((self.text, self.cursor_position))

        # Saving anything to the undo stack, clears the redo stack.
        if clear_redo_stack:
            self._redo_stack = []

    def transform_lines(self, line_index_iterator, transform_callback):
        """
        Transforms the text on a range of lines.
        When the iterator yield an index not in the range of lines that the
        document contains, it skips them silently.

        To uppercase some lines::

            new_text = transform_lines(range(5,10), lambda text: text.upper())

        :param line_index_iterator: Iterator of line numbers (int)
        :param transform_callback: callable that takes the original text of a
                                   line, and return the new text for this line.

        :returns: The new text.
        """
        # Split lines
        lines = self.text.split('\n')

        # Apply transformation
        for index in line_index_iterator:
            try:
                lines[index] = transform_callback(lines[index])
            except IndexError:
                pass

        return '\n'.join(lines)

    def transform_current_line(self, transform_callback):
        """
        Apply the given transformation function to the current line.

        :param transform_callback: callable that takes a string and return a new string.
        """
        document = self.document
        a = document.cursor_position + document.get_start_of_line_position()
        b = document.cursor_position + document.get_end_of_line_position()
        self.text = (
            document.text[:a] +
            transform_callback(document.text[a:b]) +
            document.text[b:])

    def transform_region(self, from_, to, transform_callback):
        """
        Transform a part of the input string.

        :param from_: (int) start position.
        :param to: (int) end position.
        :param transform_callback: Callable which accepts a string and returns
            the transformed string.
        """
        assert from_ < to

        self.text = ''.join([
            self.text[:from_] +
            transform_callback(self.text[from_:to]) +
            self.text[to:]
        ])

    def cursor_left(self, count=1):
        self.cursor_position += self.document.get_cursor_left_position(count=count)

    def cursor_right(self, count=1):
        self.cursor_position += self.document.get_cursor_right_position(count=count)

    def cursor_up(self, count=1):
        """ (for multiline edit). Move cursor to the previous line.  """
        original_column = self.preferred_column or self.document.cursor_position_col
        self.cursor_position += self.document.get_cursor_up_position(
            count=count, preferred_column=original_column)

        # Remember the original column for the next up/down movement.
        self.preferred_column = original_column

    def cursor_down(self, count=1):
        """ (for multiline edit). Move cursor to the next line.  """
        original_column = self.preferred_column or self.document.cursor_position_col
        self.cursor_position += self.document.get_cursor_down_position(
            count=count, preferred_column=original_column)

        # Remember the original column for the next up/down movement.
        self.preferred_column = original_column

    def auto_up(self, count=1, go_to_start_of_line_if_history_changes=False):
        """
        If we're not on the first line (of a multiline input) go a line up,
        otherwise go back in history. (If nothing is selected.)
        """
        if self.complete_state:
            self.complete_previous(count=count)
        elif self.document.cursor_position_row > 0:
            self.cursor_up(count=count)
        elif not self.selection_state:
            self.history_backward(count=count)

            # Go to the start of the line?
            if go_to_start_of_line_if_history_changes:
                self.cursor_position += self.document.get_start_of_line_position()

    def auto_down(self, count=1, go_to_start_of_line_if_history_changes=False):
        """
        If we're not on the last line (of a multiline input) go a line down,
        otherwise go forward in history. (If nothing is selected.)
        """
        if self.complete_state:
            self.complete_next(count=count)
        elif self.document.cursor_position_row < self.document.line_count - 1:
            self.cursor_down(count=count)
        elif not self.selection_state:
            self.history_forward(count=count)

            # Go to the start of the line?
            if go_to_start_of_line_if_history_changes:
                self.cursor_position += self.document.get_start_of_line_position()

    def delete_before_cursor(self, count=1):
        """
        Delete specified number of characters before cursor and return the
        deleted text.
        """
        assert count >= 0
        deleted = ''

        if self.cursor_position > 0:
            deleted = self.text[self.cursor_position - count:self.cursor_position]

            new_text = self.text[:self.cursor_position - count] + self.text[self.cursor_position:]
            new_cursor_position = self.cursor_position - len(deleted)

            # Set new Document atomically.
            self.document = Document(new_text, new_cursor_position)

        return deleted

    def delete(self, count=1):
        """
        Delete specified number of characters and Return the deleted text.
        """
        if self.cursor_position < len(self.text):
            deleted = self.document.text_after_cursor[:count]
            self.text = self.text[:self.cursor_position] + \
                self.text[self.cursor_position + len(deleted):]
            return deleted
        else:
            return ''

    def join_next_line(self, separator=' '):
        """
        Join the next line to the current one by deleting the line ending after
        the current line.
        """
        if not self.document.on_last_line:
            self.cursor_position += self.document.get_end_of_line_position()
            self.delete()

            # Remove spaces.
            self.text = (self.document.text_before_cursor + separator +
                         self.document.text_after_cursor.lstrip(' '))

    def join_selected_lines(self, separator=' '):
        """
        Join the selected lines.
        """
        assert self.selection_state

        # Get lines.
        from_, to = sorted([self.cursor_position, self.selection_state.original_cursor_position])

        before = self.text[:from_]
        lines = self.text[from_:to].splitlines()
        after = self.text[to:]

        # Replace leading spaces with just one space.
        lines = [l.lstrip(' ') + separator for l in lines]

        # Set new document.
        self.document = Document(text=before + ''.join(lines) + after,
                                 cursor_position=len(before + ''.join(lines[:-1])) - 1)

    def swap_characters_before_cursor(self):
        """
        Swap the last two characters before the cursor.
        """
        pos = self.cursor_position

        if pos >= 2:
            a = self.text[pos - 2]
            b = self.text[pos - 1]

            self.text = self.text[:pos-2] + b + a + self.text[pos:]

    def go_to_history(self, index):
        """
        Go to this item in the history.
        """
        if index < len(self._working_lines):
            self.working_index = index
            self.cursor_position = len(self.text)

    def complete_next(self, count=1, disable_wrap_around=False):
        """
        Browse to the next completions.
        (Does nothing if there are no completion.)
        """
        if self.complete_state:
            completions_count = len(self.complete_state.current_completions)

            if self.complete_state.complete_index is None:
                index = 0
            elif self.complete_state.complete_index == completions_count - 1:
                index = None

                if disable_wrap_around:
                    return
            else:
                index = min(completions_count-1, self.complete_state.complete_index + count)
            self.go_to_completion(index)

    def complete_previous(self, count=1, disable_wrap_around=False):
        """
        Browse to the previous completions.
        (Does nothing if there are no completion.)
        """
        if self.complete_state:
            if self.complete_state.complete_index == 0:
                index = None

                if disable_wrap_around:
                    return
            elif self.complete_state.complete_index is None:
                index = len(self.complete_state.current_completions) - 1
            else:
                index = max(0, self.complete_state.complete_index - count)

            self.go_to_completion(index)

    def cancel_completion(self):
        """
        Cancel completion, go back to the original text.
        """
        if self.complete_state:
            self.go_to_completion(None)
            self.complete_state = None

    def set_completions(self, completions, go_to_first=True, go_to_last=False):
        """
        Start completions. (Generate list of completions and initialize.)
        """
        assert not (go_to_first and go_to_last)

        # Generate list of all completions.
        if completions is None:
            if self.completer:
                completions = list(self.completer.get_completions(
                    self.document,
                    CompleteEvent(completion_requested=True)
                ))
            else:
                completions = []

        # Set `complete_state`.
        if completions:
            self.complete_state = CompletionState(
                original_document=self.document,
                current_completions=completions)
            if go_to_first:
                self.go_to_completion(0)
            elif go_to_last:
                self.go_to_completion(len(completions) - 1)
            else:
                self.go_to_completion(None)

        else:
            self.complete_state = None

    def start_history_lines_completion(self):
        """
        Start a completion based on all the other lines in the document and the
        history.
        """
        found_completions = set()
        completions = []

        # For every line of the whole history, find matches with the current line.
        current_line = self.document.current_line_before_cursor.lstrip()

        for i, string in enumerate(self._working_lines):
            for j, l in enumerate(string.split('\n')):
                l = l.strip()
                if l and l.startswith(current_line):
                    # When a new line has been found.
                    if l not in found_completions:
                        found_completions.add(l)

                        # Create completion.
                        if i == self.working_index:
                            display_meta = "Current, line %s" % (j+1)
                        else:
                            display_meta = "History %s, line %s" % (i+1, j+1)

                        completions.append(Completion(
                            l,
                            start_position=-len(current_line),
                            display_meta=display_meta))

        self.set_completions(completions=completions[::-1])

    def go_to_completion(self, index):
        """
        Select a completion from the list of current completions.
        """
        assert index is None or isinstance(index, int)
        assert self.complete_state

        # Set new completion
        state = self.complete_state.go_to_index(index)

        # Set text/cursor position
        new_text, new_cursor_position = state.new_text_and_position()
        self.document = Document(new_text, new_cursor_position)

        # (changing text/cursor position will unset complete_state.)
        self.complete_state = state

    def apply_completion(self, completion):
        """
        Insert a given completion.
        """
        assert isinstance(completion, Completion)

        # If there was already a completion active, cancel that one.
        if self.complete_state:
            self.go_to_completion(None)
        self.complete_state = None

        # Insert text from the given completion.
        self.delete_before_cursor(-completion.start_position)
        self.insert_text(completion.text)

    def _set_history_search(self):
        """ Set `history_search_text`. """
        if self.enable_history_search():
            if self.history_search_text is None:
                self.history_search_text = self.document.text_before_cursor
        else:
            self.history_search_text = None

    def _history_matches(self, i):
        """
        True when the current entry matches the history search.
        (when we don't have history search, it's also True.)
        """
        return (self.history_search_text is None or
                self._working_lines[i].startswith(self.history_search_text))

    def history_forward(self, count=1):
        """
        Move forwards through the history.

        :param count: Amount of items to move forward.
        """
        self._set_history_search()

        # Go forward in history.
        found_something = False

        for i in range(self.working_index + 1, len(self._working_lines)):
            if self._history_matches(i):
                self.working_index = i
                count -= 1
                found_something = True
            if count == 0:
                break

        # If we found an entry, move cursor to the end of the first line.
        if found_something:
            self.cursor_position = 0
            self.cursor_position += self.document.get_end_of_line_position()

    def history_backward(self, count=1):
        """
        Move backwards through history.
        """
        self._set_history_search()

        # Go back in history.
        found_something = False

        for i in range(self.working_index - 1, -1, -1):
            if self._history_matches(i):
                self.working_index = i
                count -= 1
                found_something = True
            if count == 0:
                break

        # If we move to another entry, move cursor to the end of the line.
        if found_something:
            self.cursor_position = len(self.text)

    def yank_nth_arg(self, n=None, _yank_last_arg=False):
        """
        Pick nth word from previous history entry (depending on current
        `yank_nth_arg_state`) and insert it at current position. Rotate through
        history if called repeatedly. If no `n` has been given, take the first
        argument. (The second word.)

        :param n: (None or int), The index of the word from the previous line
            to take.
        """
        assert n is None or isinstance(n, int)

        if not len(self.history):
            return

        # Make sure we have a `YankNthArgState`.
        if self.yank_nth_arg_state is None:
            state = YankNthArgState(n=-1 if _yank_last_arg else 1)
        else:
            state = self.yank_nth_arg_state

        if n is not None:
            state.n = n

        # Get new history position.
        new_pos = state.history_position - 1
        if -new_pos > len(self.history):
            new_pos = -1

        # Take argument from line.
        line = self.history[new_pos]

        words = [w.strip() for w in _QUOTED_WORDS_RE.split(line)]
        words = [w for w in words if w]
        try:
            word = words[state.n]
        except IndexError:
            word = ''

        # Insert new argument.
        if state.previous_inserted_word:
            self.delete_before_cursor(len(state.previous_inserted_word))
        self.insert_text(word)

        # Save state again for next completion. (Note that the 'insert'
        # operation from above clears `self.yank_nth_arg_state`.)
        state.previous_inserted_word = word
        state.history_position = new_pos
        self.yank_nth_arg_state = state

    def yank_last_arg(self, n=None):
        """
        Like `yank_nth_arg`, but if no argument has been given, yank the last
        word by default.
        """
        self.yank_nth_arg(n=n, _yank_last_arg=True)

    def start_selection(self, selection_type=SelectionType.CHARACTERS):
        """
        Take the current cursor position as the start of this selection.
        """
        self.selection_state = SelectionState(self.cursor_position, selection_type)

    def copy_selection(self, _cut=False):
        """
        Copy selected text and return :class:`.ClipboardData` instance.
        """
        new_document, clipboard_data = self.document.cut_selection()
        if _cut:
            self.document = new_document

        self.selection_state = None
        return clipboard_data

    def cut_selection(self):
        """
        Delete selected text and return :class:`.ClipboardData` instance.
        """
        return self.copy_selection(_cut=True)

    def paste_clipboard_data(self, data, paste_mode=PasteMode.EMACS, count=1):
        """
        Insert the data from the clipboard.
        """
        assert isinstance(data, ClipboardData)
        assert paste_mode in (PasteMode.VI_BEFORE, PasteMode.VI_AFTER, PasteMode.EMACS)

        original_document = self.document
        self.document = self.document.paste_clipboard_data(data, paste_mode=paste_mode, count=count)

        # Remember original document. This assignment should come at the end,
        # because assigning to 'document' will erase it.
        self.document_before_paste = original_document

    def newline(self, copy_margin=True):
        """
        Insert a line ending at the current position.
        """
        if copy_margin:
            self.insert_text('\n' + self.document.leading_whitespace_in_current_line)
        else:
            self.insert_text('\n')

    def insert_line_above(self, copy_margin=True):
        """
        Insert a new line above the current one.
        """
        if copy_margin:
            insert = self.document.leading_whitespace_in_current_line + '\n'
        else:
            insert = '\n'

        self.cursor_position += self.document.get_start_of_line_position()
        self.insert_text(insert)
        self.cursor_position -= 1

    def insert_line_below(self, copy_margin=True):
        """
        Insert a new line below the current one.
        """
        if copy_margin:
            insert = '\n' + self.document.leading_whitespace_in_current_line
        else:
            insert = '\n'

        self.cursor_position += self.document.get_end_of_line_position()
        self.insert_text(insert)

    def insert_text(self, data, overwrite=False, move_cursor=True, fire_event=True):
        """
        Insert characters at cursor position.

        :param fire_event: Fire `on_text_insert` event. This is mainly used to
            trigger autocompletion while typing.
        """
        # Original text & cursor position.
        otext = self.text
        ocpos = self.cursor_position

        # In insert/text mode.
        if overwrite:
            # Don't overwrite the newline itself. Just before the line ending,
            # it should act like insert mode.
            overwritten_text = otext[ocpos:ocpos + len(data)]
            if '\n' in overwritten_text:
                overwritten_text = overwritten_text[:overwritten_text.find('\n')]

            self.text = otext[:ocpos] + data + otext[ocpos + len(overwritten_text):]
        else:
            self.text = otext[:ocpos] + data + otext[ocpos:]

        if move_cursor:
            self.cursor_position += len(data)

        # Fire 'on_text_insert' event.
        if fire_event:
            self.on_text_insert.fire()

    def undo(self):
        # Pop from the undo-stack until we find a text that if different from
        # the current text. (The current logic of `save_to_undo_stack` will
        # cause that the top of the undo stack is usually the same as the
        # current text, so in that case we have to pop twice.)
        while self._undo_stack:
            text, pos = self._undo_stack.pop()

            if text != self.text:
                # Push current text to redo stack.
                self._redo_stack.append((self.text, self.cursor_position))

                # Set new text/cursor_position.
                self.document = Document(text, cursor_position=pos)
                break

    def redo(self):
        if self._redo_stack:
            # Copy current state on undo stack.
            self.save_to_undo_stack(clear_redo_stack=False)

            # Pop state from redo stack.
            text, pos = self._redo_stack.pop()
            self.document = Document(text, cursor_position=pos)

    def validate(self):
        """
        Returns `True` if valid.
        """
        # Don't call the validator again, if it was already called for the
        # current input.
        if self.validation_state != ValidationState.UNKNOWN:
            return self.validation_state == ValidationState.VALID

        # Validate first. If not valid, set validation exception.
        if self.validator:
            try:
                self.validator.validate(self.document)
            except ValidationError as e:
                # Set cursor position (don't allow invalid values.)
                cursor_position = e.cursor_position
                self.cursor_position = min(max(0, cursor_position), len(self.text))

                self.validation_state = ValidationState.INVALID
                self.validation_error = e
                return False

        self.validation_state = ValidationState.VALID
        self.validation_error = None
        return True

    def append_to_history(self):
        """
        Append the current input to the history.
        (Only if valid input.)
        """
        # Validate first. If not valid, set validation exception.
        if not self.validate():
            return

        # Save at the tail of the history. (But don't if the last entry the
        # history is already the same.)
        if self.text and (not len(self.history) or self.history[-1] != self.text):
            self.history.append(self.text)

    def _search(self, search_state, include_current_position=False, count=1):
        """
        Execute search. Return (working_index, cursor_position) tuple when this
        search is applied. Returns `None` when this text cannot be found.
        """
        assert isinstance(search_state, SearchState)
        assert isinstance(count, int) and count > 0

        text = search_state.text
        direction = search_state.direction
        ignore_case = search_state.ignore_case()

        def search_once(working_index, document):
            """
            Do search one time.
            Return (working_index, document) or `None`
            """
            if direction == IncrementalSearchDirection.FORWARD:
                # Try find at the current input.
                new_index = document.find(
                   text, include_current_position=include_current_position,
                   ignore_case=ignore_case)

                if new_index is not None:
                    return (working_index,
                            Document(document.text, document.cursor_position + new_index))
                else:
                    # No match, go forward in the history. (Include len+1 to wrap around.)
                    # (Here we should always include all cursor positions, because
                    # it's a different line.)
                    for i in range(working_index + 1, len(self._working_lines) + 1):
                        i %= len(self._working_lines)

                        document = Document(self._working_lines[i], 0)
                        new_index = document.find(text, include_current_position=True,
                                                  ignore_case=ignore_case)
                        if new_index is not None:
                            return (i, Document(document.text, new_index))
            else:
                # Try find at the current input.
                new_index = document.find_backwards(
                    text, ignore_case=ignore_case)

                if new_index is not None:
                    return (working_index,
                            Document(document.text, document.cursor_position + new_index))
                else:
                    # No match, go back in the history. (Include -1 to wrap around.)
                    for i in range(working_index - 1, -2, -1):
                        i %= len(self._working_lines)

                        document = Document(self._working_lines[i], len(self._working_lines[i]))
                        new_index = document.find_backwards(
                            text, ignore_case=ignore_case)
                        if new_index is not None:
                            return (i, Document(document.text, len(document.text) + new_index))

        # Do 'count' search iterations.
        working_index = self.working_index
        document = self.document
        for _ in range(count):
            result = search_once(working_index, document)
            if result is None:
                return  # Nothing found.
            else:
                working_index, document = result

        return (working_index, document.cursor_position)

    def document_for_search(self, search_state):
        """
        Return a :class:`~prompt_toolkit.document.Document` instance that has
        the text/cursor position for this search, if we would apply it. This
        will be used in the
        :class:`~prompt_toolkit.layout.controls.BufferControl` to display
        feedback while searching.
        """
        search_result = self._search(search_state, include_current_position=True)

        if search_result is None:
            return self.document
        else:
            working_index, cursor_position = search_result

            # Keep selection, when `working_index` was not changed.
            if working_index == self.working_index:
                selection = self.selection_state
            else:
                selection = None

            return Document(self._working_lines[working_index],
                            cursor_position, selection=selection)

    def get_search_position(self, search_state, include_current_position=True, count=1):
        """
        Get the cursor position for this search.
        (This operation won't change the `working_index`. It's won't go through
        the history. Vi text objects can't span multiple items.)
        """
        search_result = self._search(
            search_state, include_current_position=include_current_position, count=count)

        if search_result is None:
            return self.cursor_position
        else:
            working_index, cursor_position = search_result
            return cursor_position

    def apply_search(self, search_state, include_current_position=True, count=1):
        """
        Apply search. If something is found, set `working_index` and
        `cursor_position`.
        """
        search_result = self._search(
            search_state, include_current_position=include_current_position, count=count)

        if search_result is not None:
            working_index, cursor_position = search_result
            self.working_index = working_index
            self.cursor_position = cursor_position

    def exit_selection(self):
        self.selection_state = None

    def open_in_editor(self, cli):
        """
        Open code in editor.

        :param cli: :class:`~prompt_toolkit.interface.CommandLineInterface`
            instance.
        """
        if self.read_only():
            raise EditReadOnlyBuffer()

        # Write to temporary file
        descriptor, filename = tempfile.mkstemp(self.tempfile_suffix)
        os.write(descriptor, self.text.encode('utf-8'))
        os.close(descriptor)

        # Open in editor
        # (We need to use `cli.run_in_terminal`, because not all editors go to
        # the alternate screen buffer, and some could influence the cursor
        # position.)
        succes = cli.run_in_terminal(lambda: self._open_file_in_editor(filename))

        # Read content again.
        if succes:
            with open(filename, 'rb') as f:
                text = f.read().decode('utf-8')

                # Drop trailing newline. (Editors are supposed to add it at the
                # end, but we don't need it.)
                if text.endswith('\n'):
                    text = text[:-1]

                self.document = Document(
                    text=text,
                    cursor_position=len(text))

        # Clean up temp file.
        os.remove(filename)

    def _open_file_in_editor(self, filename):
        """
        Call editor executable.

        Return True when we received a zero return code.
        """
        # If the 'VISUAL' or 'EDITOR' environment variable has been set, use that.
        # Otherwise, fall back to the first available editor that we can find.
        visual = os.environ.get('VISUAL')
        editor = os.environ.get('EDITOR')

        editors = [
            visual,
            editor,

            # Order of preference.
            '/usr/bin/editor',
            '/usr/bin/nano',
            '/usr/bin/pico',
            '/usr/bin/vi',
            '/usr/bin/emacs',
        ]

        for e in editors:
            if e:
                try:
                    # Use 'shlex.split()', because $VISUAL can contain spaces
                    # and quotes.
                    returncode = subprocess.call(shlex.split(e) + [filename])
                    return returncode == 0

                except OSError:
                    # Executable does not exist, try the next one.
                    pass

        return False


def indent(buffer, from_row, to_row, count=1):
    """
    Indent text of a :class:`.Buffer` object.
    """
    current_row = buffer.document.cursor_position_row
    line_range = range(from_row, to_row)

    # Apply transformation.
    new_text = buffer.transform_lines(line_range, lambda l: '    ' * count + l)
    buffer.document = Document(
        new_text,
        Document(new_text).translate_row_col_to_index(current_row, 0))

    # Go to the start of the line.
    buffer.cursor_position += buffer.document.get_start_of_line_position(after_whitespace=True)


def unindent(buffer, from_row, to_row, count=1):
    """
    Unindent text of a :class:`.Buffer` object.
    """
    current_row = buffer.document.cursor_position_row
    line_range = range(from_row, to_row)

    def transform(text):
        remove = '    ' * count
        if text.startswith(remove):
            return text[len(remove):]
        else:
            return text.lstrip()

    # Apply transformation.
    new_text = buffer.transform_lines(line_range, transform)
    buffer.document = Document(
        new_text,
        Document(new_text).translate_row_col_to_index(current_row, 0))

    # Go to the start of the line.
    buffer.cursor_position += buffer.document.get_start_of_line_position(after_whitespace=True)


def reshape_text(buffer, from_row, to_row):
    """
    Reformat text, taking the width into account.
    `to_row` is included.
    (Vi 'gq' operator.)
    """
    lines = buffer.text.splitlines(True)
    lines_before = lines[:from_row]
    lines_after = lines[to_row + 1:]
    lines_to_reformat = lines[from_row:to_row + 1]

    if lines_to_reformat:
        # Take indentation from the first line.
        length = re.search(r'^\s*', lines_to_reformat[0]).end()
        indent = lines_to_reformat[0][:length].replace('\n', '')

        # Now, take all the 'words' from the lines to be reshaped.
        words = ''.join(lines_to_reformat).split()

        # And reshape.
        width = (buffer.text_width or 80) - len(indent)
        reshaped_text = [indent]
        current_width = 0
        for w in words:
            if current_width:
                if len(w) + current_width + 1 > width:
                    reshaped_text.append('\n')
                    reshaped_text.append(indent)
                    current_width = 0
                else:
                    reshaped_text.append(' ')
                    current_width += 1

            reshaped_text.append(w)
            current_width += len(w)

        if reshaped_text[-1] != '\n':
            reshaped_text.append('\n')

        # Apply result.
        buffer.document = Document(
            text=''.join(lines_before + reshaped_text + lines_after),
            cursor_position=len(''.join(lines_before + reshaped_text)))
