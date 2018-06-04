"""
Module to define and register Terminal IPython shortcuts with
:mod:`prompt_toolkit`
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import warnings
import signal
import sys
from typing import Callable


from prompt_toolkit.enums import DEFAULT_BUFFER, SEARCH_BUFFER
from prompt_toolkit.filters import (HasFocus, HasSelection, Condition,
    ViInsertMode, EmacsInsertMode, HasCompletions)
from prompt_toolkit.filters.cli import ViMode, ViNavigationMode
from prompt_toolkit.keys import Keys
from prompt_toolkit.key_binding.bindings.completion import display_completions_like_readline

from IPython.utils.decorators import undoc

@undoc
@Condition
def cursor_in_leading_ws(cli):
    before = cli.application.buffer.document.current_line_before_cursor
    return (not before) or before.isspace()

def register_ipython_shortcuts(registry, shell):
    """Set up the prompt_toolkit keyboard shortcuts for IPython"""
    insert_mode = ViInsertMode() | EmacsInsertMode()

    if getattr(shell, 'handle_return', None):
        return_handler = shell.handle_return(shell)
    else:
        return_handler = newline_or_execute_outer(shell)

    # Ctrl+J == Enter, seemingly
    registry.add_binding(Keys.ControlJ,
                         filter=(HasFocus(DEFAULT_BUFFER)
                                 & ~HasSelection()
                                 & insert_mode
                        ))(return_handler)

    registry.add_binding(Keys.ControlBackslash)(force_exit)

    registry.add_binding(Keys.ControlP,
                         filter=(ViInsertMode() & HasFocus(DEFAULT_BUFFER)
                        ))(previous_history_or_previous_completion)

    registry.add_binding(Keys.ControlN,
                         filter=(ViInsertMode() & HasFocus(DEFAULT_BUFFER)
                        ))(next_history_or_next_completion)

    registry.add_binding(Keys.ControlG,
                         filter=(HasFocus(DEFAULT_BUFFER) & HasCompletions()
                        ))(dismiss_completion)

    registry.add_binding(Keys.ControlC, filter=HasFocus(DEFAULT_BUFFER)
                        )(reset_buffer)

    registry.add_binding(Keys.ControlC, filter=HasFocus(SEARCH_BUFFER)
                        )(reset_search_buffer)

    supports_suspend = Condition(lambda cli: hasattr(signal, 'SIGTSTP'))
    registry.add_binding(Keys.ControlZ, filter=supports_suspend
                        )(suspend_to_bg)

    # Ctrl+I == Tab
    registry.add_binding(Keys.ControlI,
                         filter=(HasFocus(DEFAULT_BUFFER)
                                 & ~HasSelection()
                                 & insert_mode
                                 & cursor_in_leading_ws
                        ))(indent_buffer)

    registry.add_binding(Keys.ControlO,
                         filter=(HasFocus(DEFAULT_BUFFER)
                                & EmacsInsertMode()))(newline_autoindent_outer(shell.input_splitter))

    registry.add_binding(Keys.F2,
                         filter=HasFocus(DEFAULT_BUFFER)
                        )(open_input_in_editor)

    if shell.display_completions == 'readlinelike':
        registry.add_binding(Keys.ControlI,
                             filter=(HasFocus(DEFAULT_BUFFER)
                                     & ~HasSelection()
                                     & insert_mode
                                     & ~cursor_in_leading_ws
                            ))(display_completions_like_readline)

    if sys.platform == 'win32':
        registry.add_binding(Keys.ControlV,
                             filter=(
                             HasFocus(
                             DEFAULT_BUFFER) & ~ViMode()
                            ))(win_paste)


def newline_or_execute_outer(shell):
    def newline_or_execute(event):
        """When the user presses return, insert a newline or execute the code."""
        b = event.current_buffer
        d = b.document

        if b.complete_state:
            cc = b.complete_state.current_completion
            if cc:
                b.apply_completion(cc)
            else:
                b.cancel_completion()
            return

        # If there's only one line, treat it as if the cursor is at the end.
        # See https://github.com/ipython/ipython/issues/10425
        if d.line_count == 1:
            check_text = d.text
        else:
            check_text = d.text[:d.cursor_position]
        status, indent = shell.input_splitter.check_complete(check_text + '\n')

        if not (d.on_last_line or
                d.cursor_position_row >= d.line_count - d.empty_line_count_at_the_end()
                ):
            b.insert_text('\n' + (' ' * (indent or 0)))
            return

        if (status != 'incomplete') and b.accept_action.is_returnable:
            b.accept_action.validate_and_handle(event.cli, b)
        else:
            b.insert_text('\n' + (' ' * (indent or 0)))
    return newline_or_execute


def previous_history_or_previous_completion(event):
    """
    Control-P in vi edit mode on readline is history next, unlike default prompt toolkit.

    If completer is open this still select previous completion.
    """
    event.current_buffer.auto_up()


def next_history_or_next_completion(event):
    """
    Control-N in vi edit mode on readline is history previous, unlike default prompt toolkit.

    If completer is open this still select next completion.
    """
    event.current_buffer.auto_down()


def dismiss_completion(event):
    b = event.current_buffer
    if b.complete_state:
        b.cancel_completion()


def reset_buffer(event):
    b = event.current_buffer
    if b.complete_state:
        b.cancel_completion()
    else:
        b.reset()


def reset_search_buffer(event):
    if event.current_buffer.document.text:
        event.current_buffer.reset()
    else:
        event.cli.push_focus(DEFAULT_BUFFER)

def suspend_to_bg(event):
    event.cli.suspend_to_background()

def force_exit(event):
    """
    Force exit (with a non-zero return value)
    """
    sys.exit("Quit")

def indent_buffer(event):
    event.current_buffer.insert_text(' ' * 4)

@undoc
def newline_with_copy_margin(event):
    """
    DEPRECATED since IPython 6.0

    See :any:`newline_autoindent_outer` for a replacement.

    Preserve margin and cursor position when using
    Control-O to insert a newline in EMACS mode
    """
    warnings.warn("`newline_with_copy_margin(event)` is deprecated since IPython 6.0. "
      "see `newline_autoindent_outer(shell)(event)` for a replacement.",
                  DeprecationWarning, stacklevel=2)

    b = event.current_buffer
    cursor_start_pos = b.document.cursor_position_col
    b.newline(copy_margin=True)
    b.cursor_up(count=1)
    cursor_end_pos = b.document.cursor_position_col
    if cursor_start_pos != cursor_end_pos:
        pos_diff = cursor_start_pos - cursor_end_pos
        b.cursor_right(count=pos_diff)

def newline_autoindent_outer(inputsplitter) -> Callable[..., None]:
    """
    Return a function suitable for inserting a indented newline after the cursor.

    Fancier version of deprecated ``newline_with_copy_margin`` which should
    compute the correct indentation of the inserted line. That is to say, indent
    by 4 extra space after a function definition, class definition, context
    manager... And dedent by 4 space after ``pass``, ``return``, ``raise ...``.
    """

    def newline_autoindent(event):
        """insert a newline after the cursor indented appropriately."""
        b = event.current_buffer
        d = b.document

        if b.complete_state:
            b.cancel_completion()
        text = d.text[:d.cursor_position] + '\n'
        _, indent = inputsplitter.check_complete(text)
        b.insert_text('\n' + (' ' * (indent or 0)), move_cursor=False)

    return newline_autoindent


def open_input_in_editor(event):
    event.cli.current_buffer.tempfile_suffix = ".py"
    event.cli.current_buffer.open_in_editor(event.cli)


if sys.platform == 'win32':
    from IPython.core.error import TryNext
    from IPython.lib.clipboard import (ClipboardEmpty,
                                       win32_clipboard_get,
                                       tkinter_clipboard_get)

    @undoc
    def win_paste(event):
        try:
            text = win32_clipboard_get()
        except TryNext:
            try:
                text = tkinter_clipboard_get()
            except (TryNext, ClipboardEmpty):
                return
        except ClipboardEmpty:
            return
        event.current_buffer.insert_text(text.replace('\t', ' ' * 4))
