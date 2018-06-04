# pylint: disable=function-redefined
from __future__ import unicode_literals

from prompt_toolkit.enums import DEFAULT_BUFFER
from prompt_toolkit.filters import HasSelection, Condition, EmacsInsertMode, ViInsertMode
from prompt_toolkit.keys import Keys
from prompt_toolkit.layout.screen import Point
from prompt_toolkit.mouse_events import MouseEventType, MouseEvent
from prompt_toolkit.renderer import HeightIsUnknownError
from prompt_toolkit.utils import suspend_to_background_supported, is_windows

from .named_commands import get_by_name
from ..registry import Registry


__all__ = (
    'load_basic_bindings',
    'load_abort_and_exit_bindings',
    'load_basic_system_bindings',
    'load_auto_suggestion_bindings',
)

def if_no_repeat(event):
    """ Callable that returns True when the previous event was delivered to
    another handler. """
    return not event.is_repeat


def load_basic_bindings():
    registry = Registry()
    insert_mode = ViInsertMode() | EmacsInsertMode()
    handle = registry.add_binding
    has_selection = HasSelection()

    @handle(Keys.ControlA)
    @handle(Keys.ControlB)
    @handle(Keys.ControlC)
    @handle(Keys.ControlD)
    @handle(Keys.ControlE)
    @handle(Keys.ControlF)
    @handle(Keys.ControlG)
    @handle(Keys.ControlH)
    @handle(Keys.ControlI)
    @handle(Keys.ControlJ)
    @handle(Keys.ControlK)
    @handle(Keys.ControlL)
    @handle(Keys.ControlM)
    @handle(Keys.ControlN)
    @handle(Keys.ControlO)
    @handle(Keys.ControlP)
    @handle(Keys.ControlQ)
    @handle(Keys.ControlR)
    @handle(Keys.ControlS)
    @handle(Keys.ControlT)
    @handle(Keys.ControlU)
    @handle(Keys.ControlV)
    @handle(Keys.ControlW)
    @handle(Keys.ControlX)
    @handle(Keys.ControlY)
    @handle(Keys.ControlZ)
    @handle(Keys.F1)
    @handle(Keys.F2)
    @handle(Keys.F3)
    @handle(Keys.F4)
    @handle(Keys.F5)
    @handle(Keys.F6)
    @handle(Keys.F7)
    @handle(Keys.F8)
    @handle(Keys.F9)
    @handle(Keys.F10)
    @handle(Keys.F11)
    @handle(Keys.F12)
    @handle(Keys.F13)
    @handle(Keys.F14)
    @handle(Keys.F15)
    @handle(Keys.F16)
    @handle(Keys.F17)
    @handle(Keys.F18)
    @handle(Keys.F19)
    @handle(Keys.F20)
    @handle(Keys.ControlSpace)
    @handle(Keys.ControlBackslash)
    @handle(Keys.ControlSquareClose)
    @handle(Keys.ControlCircumflex)
    @handle(Keys.ControlUnderscore)
    @handle(Keys.Backspace)
    @handle(Keys.Up)
    @handle(Keys.Down)
    @handle(Keys.Right)
    @handle(Keys.Left)
    @handle(Keys.ShiftUp)
    @handle(Keys.ShiftDown)
    @handle(Keys.ShiftRight)
    @handle(Keys.ShiftLeft)
    @handle(Keys.Home)
    @handle(Keys.End)
    @handle(Keys.Delete)
    @handle(Keys.ShiftDelete)
    @handle(Keys.ControlDelete)
    @handle(Keys.PageUp)
    @handle(Keys.PageDown)
    @handle(Keys.BackTab)
    @handle(Keys.Tab)
    @handle(Keys.ControlLeft)
    @handle(Keys.ControlRight)
    @handle(Keys.ControlUp)
    @handle(Keys.ControlDown)
    @handle(Keys.Insert)
    @handle(Keys.Ignore)
    def _(event):
        """
        First, for any of these keys, Don't do anything by default. Also don't
        catch them in the 'Any' handler which will insert them as data.

        If people want to insert these characters as a literal, they can always
        do by doing a quoted insert. (ControlQ in emacs mode, ControlV in Vi
        mode.)
        """
        pass

    # Readline-style bindings.
    handle(Keys.Home)(get_by_name('beginning-of-line'))
    handle(Keys.End)(get_by_name('end-of-line'))
    handle(Keys.Left)(get_by_name('backward-char'))
    handle(Keys.Right)(get_by_name('forward-char'))
    handle(Keys.ControlUp)(get_by_name('previous-history'))
    handle(Keys.ControlDown)(get_by_name('next-history'))
    handle(Keys.ControlL)(get_by_name('clear-screen'))

    handle(Keys.ControlK, filter=insert_mode)(get_by_name('kill-line'))
    handle(Keys.ControlU, filter=insert_mode)(get_by_name('unix-line-discard'))
    handle(Keys.ControlH, filter=insert_mode, save_before=if_no_repeat)(
        get_by_name('backward-delete-char'))
    handle(Keys.Backspace, filter=insert_mode, save_before=if_no_repeat)(
        get_by_name('backward-delete-char'))
    handle(Keys.Delete, filter=insert_mode, save_before=if_no_repeat)(
        get_by_name('delete-char'))
    handle(Keys.ShiftDelete, filter=insert_mode, save_before=if_no_repeat)(
        get_by_name('delete-char'))
    handle(Keys.Any, filter=insert_mode, save_before=if_no_repeat)(
        get_by_name('self-insert'))
    handle(Keys.ControlT, filter=insert_mode)(get_by_name('transpose-chars'))
    handle(Keys.ControlW, filter=insert_mode)(get_by_name('unix-word-rubout'))
    handle(Keys.ControlI, filter=insert_mode)(get_by_name('menu-complete'))
    handle(Keys.BackTab, filter=insert_mode)(get_by_name('menu-complete-backward'))

    handle(Keys.PageUp, filter= ~has_selection)(get_by_name('previous-history'))
    handle(Keys.PageDown, filter= ~has_selection)(get_by_name('next-history'))

    # CTRL keys.

    text_before_cursor = Condition(lambda cli: cli.current_buffer.text)
    handle(Keys.ControlD, filter=text_before_cursor & insert_mode)(get_by_name('delete-char'))

    is_multiline = Condition(lambda cli: cli.current_buffer.is_multiline())
    is_returnable = Condition(lambda cli: cli.current_buffer.accept_action.is_returnable)

    @handle(Keys.ControlJ, filter=is_multiline & insert_mode)
    def _(event):
        " Newline (in case of multiline input. "
        event.current_buffer.newline(copy_margin=not event.cli.in_paste_mode)

    @handle(Keys.ControlJ, filter=~is_multiline & is_returnable)
    def _(event):
        " Enter, accept input. "
        buff = event.current_buffer
        buff.accept_action.validate_and_handle(event.cli, buff)

    # Delete the word before the cursor.

    @handle(Keys.Up)
    def _(event):
        event.current_buffer.auto_up(count=event.arg)

    @handle(Keys.Down)
    def _(event):
        event.current_buffer.auto_down(count=event.arg)

    @handle(Keys.Delete, filter=has_selection)
    def _(event):
        data = event.current_buffer.cut_selection()
        event.cli.clipboard.set_data(data)

    # Global bindings.

    @handle(Keys.ControlZ)
    def _(event):
        """
        By default, control-Z should literally insert Ctrl-Z.
        (Ansi Ctrl-Z, code 26 in MSDOS means End-Of-File.
        In a Python REPL for instance, it's possible to type
        Control-Z followed by enter to quit.)

        When the system bindings are loaded and suspend-to-background is
        supported, that will override this binding.
        """
        event.current_buffer.insert_text(event.data)

    @handle(Keys.CPRResponse, save_before=lambda e: False)
    def _(event):
        """
        Handle incoming Cursor-Position-Request response.
        """
        # The incoming data looks like u'\x1b[35;1R'
        # Parse row/col information.
        row, col = map(int, event.data[2:-1].split(';'))

        # Report absolute cursor position to the renderer.
        event.cli.renderer.report_absolute_cursor_row(row)

    @handle(Keys.BracketedPaste)
    def _(event):
        " Pasting from clipboard. "
        data = event.data

        # Be sure to use \n as line ending.
        # Some terminals (Like iTerm2) seem to paste \r\n line endings in a
        # bracketed paste. See: https://github.com/ipython/ipython/issues/9737
        data = data.replace('\r\n', '\n')
        data = data.replace('\r', '\n')

        event.current_buffer.insert_text(data)

    @handle(Keys.Any, filter=Condition(lambda cli: cli.quoted_insert), eager=True)
    def _(event):
        """
        Handle quoted insert.
        """
        event.current_buffer.insert_text(event.data, overwrite=False)
        event.cli.quoted_insert = False

    return registry


def load_mouse_bindings():
    """
    Key bindings, required for mouse support.
    (Mouse events enter through the key binding system.)
    """
    registry = Registry()

    @registry.add_binding(Keys.Vt100MouseEvent)
    def _(event):
        """
        Handling of incoming mouse event.
        """
        # Typical:   "Esc[MaB*"
        # Urxvt:     "Esc[96;14;13M"
        # Xterm SGR: "Esc[<64;85;12M"

        # Parse incoming packet.
        if event.data[2] == 'M':
            # Typical.
            mouse_event, x, y = map(ord, event.data[3:])
            mouse_event = {
                32: MouseEventType.MOUSE_DOWN,
                35: MouseEventType.MOUSE_UP,
                96: MouseEventType.SCROLL_UP,
                97: MouseEventType.SCROLL_DOWN,
            }.get(mouse_event)

            # Handle situations where `PosixStdinReader` used surrogateescapes.
            if x >= 0xdc00: x-= 0xdc00
            if y >= 0xdc00: y-= 0xdc00

            x -= 32
            y -= 32
        else:
            # Urxvt and Xterm SGR.
            # When the '<' is not present, we are not using the Xterm SGR mode,
            # but Urxvt instead.
            data = event.data[2:]
            if data[:1] == '<':
                sgr = True
                data = data[1:]
            else:
                sgr = False

            # Extract coordinates.
            mouse_event, x, y = map(int, data[:-1].split(';'))
            m = data[-1]

            # Parse event type.
            if sgr:
                mouse_event = {
                    (0, 'M'): MouseEventType.MOUSE_DOWN,
                    (0, 'm'): MouseEventType.MOUSE_UP,
                    (64, 'M'): MouseEventType.SCROLL_UP,
                    (65, 'M'): MouseEventType.SCROLL_DOWN,
                }.get((mouse_event, m))
            else:
                mouse_event = {
                    32: MouseEventType.MOUSE_DOWN,
                    35: MouseEventType.MOUSE_UP,
                    96: MouseEventType.SCROLL_UP,
                    97: MouseEventType.SCROLL_DOWN,
                    }.get(mouse_event)

        x -= 1
        y -= 1

        # Only handle mouse events when we know the window height.
        if event.cli.renderer.height_is_known and mouse_event is not None:
            # Take region above the layout into account. The reported
            # coordinates are absolute to the visible part of the terminal.
            try:
                y -= event.cli.renderer.rows_above_layout
            except HeightIsUnknownError:
                return

            # Call the mouse handler from the renderer.
            handler = event.cli.renderer.mouse_handlers.mouse_handlers[x,y]
            handler(event.cli, MouseEvent(position=Point(x=x, y=y),
                                          event_type=mouse_event))

    @registry.add_binding(Keys.WindowsMouseEvent)
    def _(event):
        """
        Handling of mouse events for Windows.
        """
        assert is_windows()  # This key binding should only exist for Windows.

        # Parse data.
        event_type, x, y = event.data.split(';')
        x = int(x)
        y = int(y)

        # Make coordinates absolute to the visible part of the terminal.
        screen_buffer_info = event.cli.renderer.output.get_win32_screen_buffer_info()
        rows_above_cursor = screen_buffer_info.dwCursorPosition.Y - event.cli.renderer._cursor_pos.y
        y -= rows_above_cursor

        # Call the mouse event handler.
        handler = event.cli.renderer.mouse_handlers.mouse_handlers[x,y]
        handler(event.cli, MouseEvent(position=Point(x=x, y=y),
                                      event_type=event_type))

    return registry


def load_abort_and_exit_bindings():
    """
    Basic bindings for abort (Ctrl-C) and exit (Ctrl-D).
    """
    registry = Registry()
    handle = registry.add_binding

    @handle(Keys.ControlC)
    def _(event):
        " Abort when Control-C has been pressed. "
        event.cli.abort()

    @Condition
    def ctrl_d_condition(cli):
        """ Ctrl-D binding is only active when the default buffer is selected
        and empty. """
        return (cli.current_buffer_name == DEFAULT_BUFFER and
                not cli.current_buffer.text)

    handle(Keys.ControlD, filter=ctrl_d_condition)(get_by_name('end-of-file'))

    return registry


def load_basic_system_bindings():
    """
    Basic system bindings (For both Emacs and Vi mode.)
    """
    registry = Registry()

    suspend_supported = Condition(
        lambda cli: suspend_to_background_supported())

    @registry.add_binding(Keys.ControlZ, filter=suspend_supported)
    def _(event):
        """
        Suspend process to background.
        """
        event.cli.suspend_to_background()

    return registry


def load_auto_suggestion_bindings():
    """
    Key bindings for accepting auto suggestion text.
    """
    registry = Registry()
    handle = registry.add_binding

    suggestion_available = Condition(
        lambda cli:
            cli.current_buffer.suggestion is not None and
            cli.current_buffer.document.is_cursor_at_the_end)

    @handle(Keys.ControlF, filter=suggestion_available)
    @handle(Keys.ControlE, filter=suggestion_available)
    @handle(Keys.Right, filter=suggestion_available)
    def _(event):
        " Accept suggestion. "
        b = event.current_buffer
        suggestion = b.suggestion

        if suggestion:
            b.insert_text(suggestion.text)

    return registry
