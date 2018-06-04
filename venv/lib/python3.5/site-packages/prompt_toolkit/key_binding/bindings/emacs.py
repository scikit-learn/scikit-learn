# pylint: disable=function-redefined
from __future__ import unicode_literals
from prompt_toolkit.buffer import SelectionType, indent, unindent
from prompt_toolkit.keys import Keys
from prompt_toolkit.enums import IncrementalSearchDirection, SEARCH_BUFFER, SYSTEM_BUFFER
from prompt_toolkit.filters import Condition, EmacsMode, HasSelection, EmacsInsertMode, HasFocus, HasArg
from prompt_toolkit.completion import CompleteEvent

from .scroll import scroll_page_up, scroll_page_down
from .named_commands import get_by_name
from ..registry import Registry, ConditionalRegistry

__all__ = (
    'load_emacs_bindings',
    'load_emacs_search_bindings',
    'load_emacs_system_bindings',
    'load_extra_emacs_page_navigation_bindings',
)


def load_emacs_bindings():
    """
    Some e-macs extensions.
    """
    # Overview of Readline emacs commands:
    # http://www.catonmat.net/download/readline-emacs-editing-mode-cheat-sheet.pdf
    registry = ConditionalRegistry(Registry(), EmacsMode())
    handle = registry.add_binding

    insert_mode = EmacsInsertMode()
    has_selection = HasSelection()

    @handle(Keys.Escape)
    def _(event):
        """
        By default, ignore escape key.

        (If we don't put this here, and Esc is followed by a key which sequence
        is not handled, we'll insert an Escape character in the input stream.
        Something we don't want and happens to easily in emacs mode.
        Further, people can always use ControlQ to do a quoted insert.)
        """
        pass

    handle(Keys.ControlA)(get_by_name('beginning-of-line'))
    handle(Keys.ControlB)(get_by_name('backward-char'))
    handle(Keys.ControlDelete, filter=insert_mode)(get_by_name('kill-word'))
    handle(Keys.ControlE)(get_by_name('end-of-line'))
    handle(Keys.ControlF)(get_by_name('forward-char'))
    handle(Keys.ControlLeft)(get_by_name('backward-word'))
    handle(Keys.ControlRight)(get_by_name('forward-word'))
    handle(Keys.ControlX, 'r', 'y', filter=insert_mode)(get_by_name('yank'))
    handle(Keys.ControlY, filter=insert_mode)(get_by_name('yank'))
    handle(Keys.Escape, 'b')(get_by_name('backward-word'))
    handle(Keys.Escape, 'c', filter=insert_mode)(get_by_name('capitalize-word'))
    handle(Keys.Escape, 'd', filter=insert_mode)(get_by_name('kill-word'))
    handle(Keys.Escape, 'f')(get_by_name('forward-word'))
    handle(Keys.Escape, 'l', filter=insert_mode)(get_by_name('downcase-word'))
    handle(Keys.Escape, 'u', filter=insert_mode)(get_by_name('uppercase-word'))
    handle(Keys.Escape, 'y', filter=insert_mode)(get_by_name('yank-pop'))
    handle(Keys.Escape, Keys.ControlH, filter=insert_mode)(get_by_name('backward-kill-word'))
    handle(Keys.Escape, Keys.Backspace, filter=insert_mode)(get_by_name('backward-kill-word'))
    handle(Keys.Escape, '\\', filter=insert_mode)(get_by_name('delete-horizontal-space'))

    handle(Keys.ControlUnderscore, save_before=(lambda e: False), filter=insert_mode)(
        get_by_name('undo'))

    handle(Keys.ControlX, Keys.ControlU, save_before=(lambda e: False), filter=insert_mode)(
        get_by_name('undo'))


    handle(Keys.Escape, '<', filter= ~has_selection)(get_by_name('beginning-of-history'))
    handle(Keys.Escape, '>', filter= ~has_selection)(get_by_name('end-of-history'))

    handle(Keys.Escape, '.', filter=insert_mode)(get_by_name('yank-last-arg'))
    handle(Keys.Escape, '_', filter=insert_mode)(get_by_name('yank-last-arg'))
    handle(Keys.Escape, Keys.ControlY, filter=insert_mode)(get_by_name('yank-nth-arg'))
    handle(Keys.Escape, '#', filter=insert_mode)(get_by_name('insert-comment'))
    handle(Keys.ControlO)(get_by_name('operate-and-get-next'))

    # ControlQ does a quoted insert. Not that for vt100 terminals, you have to
    # disable flow control by running ``stty -ixon``, otherwise Ctrl-Q and
    # Ctrl-S are captured by the terminal.
    handle(Keys.ControlQ, filter= ~has_selection)(get_by_name('quoted-insert'))

    handle(Keys.ControlX, '(')(get_by_name('start-kbd-macro'))
    handle(Keys.ControlX, ')')(get_by_name('end-kbd-macro'))
    handle(Keys.ControlX, 'e')(get_by_name('call-last-kbd-macro'))

    @handle(Keys.ControlN)
    def _(event):
        " Next line. "
        event.current_buffer.auto_down()

    @handle(Keys.ControlP)
    def _(event):
        " Previous line. "
        event.current_buffer.auto_up(count=event.arg)

    def handle_digit(c):
        """
        Handle input of arguments.
        The first number needs to be preceeded by escape.
        """
        @handle(c, filter=HasArg())
        @handle(Keys.Escape, c)
        def _(event):
            event.append_to_arg_count(c)

    for c in '0123456789':
        handle_digit(c)

    @handle(Keys.Escape, '-', filter=~HasArg())
    def _(event):
        """
        """
        if event._arg is None:
            event.append_to_arg_count('-')

    @handle('-', filter=Condition(lambda cli: cli.input_processor.arg == '-'))
    def _(event):
        """
        When '-' is typed again, after exactly '-' has been given as an
        argument, ignore this.
        """
        event.cli.input_processor.arg = '-'

    is_returnable = Condition(
        lambda cli: cli.current_buffer.accept_action.is_returnable)

    # Meta + Newline: always accept input.
    handle(Keys.Escape, Keys.ControlJ, filter=insert_mode & is_returnable)(
        get_by_name('accept-line'))

    def character_search(buff, char, count):
        if count < 0:
            match = buff.document.find_backwards(char, in_current_line=True, count=-count)
        else:
            match = buff.document.find(char, in_current_line=True, count=count)

        if match is not None:
            buff.cursor_position += match

    @handle(Keys.ControlSquareClose, Keys.Any)
    def _(event):
        " When Ctl-] + a character is pressed. go to that character. "
        # Also named 'character-search'
        character_search(event.current_buffer, event.data, event.arg)

    @handle(Keys.Escape, Keys.ControlSquareClose, Keys.Any)
    def _(event):
        " Like Ctl-], but backwards. "
        # Also named 'character-search-backward'
        character_search(event.current_buffer, event.data, -event.arg)

    @handle(Keys.Escape, 'a')
    def _(event):
        " Previous sentence. "
        # TODO:

    @handle(Keys.Escape, 'e')
    def _(event):
        " Move to end of sentence. "
        # TODO:

    @handle(Keys.Escape, 't', filter=insert_mode)
    def _(event):
        """
        Swap the last two words before the cursor.
        """
        # TODO

    @handle(Keys.Escape, '*', filter=insert_mode)
    def _(event):
        """
        `meta-*`: Insert all possible completions of the preceding text.
        """
        buff = event.current_buffer

        # List all completions.
        complete_event = CompleteEvent(text_inserted=False, completion_requested=True)
        completions = list(buff.completer.get_completions(buff.document, complete_event))

        # Insert them.
        text_to_insert = ' '.join(c.text for c in completions)
        buff.insert_text(text_to_insert)

    @handle(Keys.ControlX, Keys.ControlX)
    def _(event):
        """
        Move cursor back and forth between the start and end of the current
        line.
        """
        buffer = event.current_buffer

        if buffer.document.is_cursor_at_the_end_of_line:
            buffer.cursor_position += buffer.document.get_start_of_line_position(after_whitespace=False)
        else:
            buffer.cursor_position += buffer.document.get_end_of_line_position()

    @handle(Keys.ControlSpace)
    def _(event):
        """
        Start of the selection (if the current buffer is not empty).
        """
        # Take the current cursor position as the start of this selection.
        buff = event.current_buffer
        if buff.text:
            buff.start_selection(selection_type=SelectionType.CHARACTERS)

    @handle(Keys.ControlG, filter= ~has_selection)
    def _(event):
        """
        Control + G: Cancel completion menu and validation state.
        """
        event.current_buffer.complete_state = None
        event.current_buffer.validation_error = None

    @handle(Keys.ControlG, filter=has_selection)
    def _(event):
        """
        Cancel selection.
        """
        event.current_buffer.exit_selection()

    @handle(Keys.ControlW, filter=has_selection)
    @handle(Keys.ControlX, 'r', 'k', filter=has_selection)
    def _(event):
        """
        Cut selected text.
        """
        data = event.current_buffer.cut_selection()
        event.cli.clipboard.set_data(data)

    @handle(Keys.Escape, 'w', filter=has_selection)
    def _(event):
        """
        Copy selected text.
        """
        data = event.current_buffer.copy_selection()
        event.cli.clipboard.set_data(data)

    @handle(Keys.Escape, Keys.Left)
    def _(event):
        """
        Cursor to start of previous word.
        """
        buffer = event.current_buffer
        buffer.cursor_position += buffer.document.find_previous_word_beginning(count=event.arg) or 0

    @handle(Keys.Escape, Keys.Right)
    def _(event):
        """
        Cursor to start of next word.
        """
        buffer = event.current_buffer
        buffer.cursor_position += buffer.document.find_next_word_beginning(count=event.arg) or \
            buffer.document.get_end_of_document_position()

    @handle(Keys.Escape, '/', filter=insert_mode)
    def _(event):
        """
        M-/: Complete.
        """
        b = event.current_buffer
        if b.complete_state:
            b.complete_next()
        else:
            event.cli.start_completion(select_first=True)

    @handle(Keys.ControlC, '>', filter=has_selection)
    def _(event):
        """
        Indent selected text.
        """
        buffer = event.current_buffer

        buffer.cursor_position += buffer.document.get_start_of_line_position(after_whitespace=True)

        from_, to = buffer.document.selection_range()
        from_, _ = buffer.document.translate_index_to_position(from_)
        to, _ = buffer.document.translate_index_to_position(to)

        indent(buffer, from_, to + 1, count=event.arg)

    @handle(Keys.ControlC, '<', filter=has_selection)
    def _(event):
        """
        Unindent selected text.
        """
        buffer = event.current_buffer

        from_, to = buffer.document.selection_range()
        from_, _ = buffer.document.translate_index_to_position(from_)
        to, _ = buffer.document.translate_index_to_position(to)

        unindent(buffer, from_, to + 1, count=event.arg)

    return registry


def load_emacs_open_in_editor_bindings():
    """
    Pressing C-X C-E will open the buffer in an external editor.
    """
    registry = Registry()

    registry.add_binding(Keys.ControlX, Keys.ControlE,
                         filter=EmacsMode() & ~HasSelection())(
         get_by_name('edit-and-execute-command'))

    return registry


def load_emacs_system_bindings():
    registry = ConditionalRegistry(Registry(), EmacsMode())
    handle = registry.add_binding

    has_focus = HasFocus(SYSTEM_BUFFER)

    @handle(Keys.Escape, '!', filter= ~has_focus)
    def _(event):
        """
        M-'!' opens the system prompt.
        """
        event.cli.push_focus(SYSTEM_BUFFER)

    @handle(Keys.Escape, filter=has_focus)
    @handle(Keys.ControlG, filter=has_focus)
    @handle(Keys.ControlC, filter=has_focus)
    def _(event):
        """
        Cancel system prompt.
        """
        event.cli.buffers[SYSTEM_BUFFER].reset()
        event.cli.pop_focus()

    @handle(Keys.ControlJ, filter=has_focus)
    def _(event):
        """
        Run system command.
        """
        system_line = event.cli.buffers[SYSTEM_BUFFER]
        event.cli.run_system_command(system_line.text)
        system_line.reset(append_to_history=True)

        # Focus previous buffer again.
        event.cli.pop_focus()

    return registry


def load_emacs_search_bindings(get_search_state=None):
    registry = ConditionalRegistry(Registry(), EmacsMode())
    handle = registry.add_binding

    has_focus = HasFocus(SEARCH_BUFFER)

    assert get_search_state is None or callable(get_search_state)

    if not get_search_state:
        def get_search_state(cli): return cli.search_state

    @handle(Keys.ControlG, filter=has_focus)
    @handle(Keys.ControlC, filter=has_focus)
    # NOTE: the reason for not also binding Escape to this one, is that we want
    #       Alt+Enter to accept input directly in incremental search mode.
    def _(event):
        """
        Abort an incremental search and restore the original line.
        """
        search_buffer = event.cli.buffers[SEARCH_BUFFER]

        search_buffer.reset()
        event.cli.pop_focus()

    @handle(Keys.ControlJ, filter=has_focus)
    @handle(Keys.Escape, filter=has_focus, eager=True)
    def _(event):
        """
        When enter pressed in isearch, quit isearch mode. (Multiline
        isearch would be too complicated.)
        """
        input_buffer = event.cli.buffers.previous(event.cli)
        search_buffer = event.cli.buffers[SEARCH_BUFFER]

        # Update search state.
        if search_buffer.text:
            get_search_state(event.cli).text = search_buffer.text

        # Apply search.
        input_buffer.apply_search(get_search_state(event.cli), include_current_position=True)

        # Add query to history of search line.
        search_buffer.append_to_history()
        search_buffer.reset()

        # Focus previous document again.
        event.cli.pop_focus()

    @handle(Keys.ControlR, filter= ~has_focus)
    def _(event):
        get_search_state(event.cli).direction = IncrementalSearchDirection.BACKWARD
        event.cli.push_focus(SEARCH_BUFFER)

    @handle(Keys.ControlS, filter= ~has_focus)
    def _(event):
        get_search_state(event.cli).direction = IncrementalSearchDirection.FORWARD
        event.cli.push_focus(SEARCH_BUFFER)

    def incremental_search(cli, direction, count=1):
        " Apply search, but keep search buffer focussed. "
        # Update search_state.
        search_state = get_search_state(cli)
        direction_changed = search_state.direction != direction

        search_state.text = cli.buffers[SEARCH_BUFFER].text
        search_state.direction = direction

        # Apply search to current buffer.
        if not direction_changed:
            input_buffer = cli.buffers.previous(cli)
            input_buffer.apply_search(search_state,
                                      include_current_position=False, count=count)

    @handle(Keys.ControlR, filter=has_focus)
    @handle(Keys.Up, filter=has_focus)
    def _(event):
        incremental_search(event.cli, IncrementalSearchDirection.BACKWARD, count=event.arg)

    @handle(Keys.ControlS, filter=has_focus)
    @handle(Keys.Down, filter=has_focus)
    def _(event):
        incremental_search(event.cli, IncrementalSearchDirection.FORWARD, count=event.arg)

    return registry


def load_extra_emacs_page_navigation_bindings():
    """
    Key bindings, for scrolling up and down through pages.
    This are separate bindings, because GNU readline doesn't have them.
    """
    registry = ConditionalRegistry(Registry(), EmacsMode())
    handle = registry.add_binding

    handle(Keys.ControlV)(scroll_page_down)
    handle(Keys.PageDown)(scroll_page_down)
    handle(Keys.Escape, 'v')(scroll_page_up)
    handle(Keys.PageUp)(scroll_page_up)

    return registry
