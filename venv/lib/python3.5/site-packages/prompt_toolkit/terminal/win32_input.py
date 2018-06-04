from __future__ import unicode_literals
from ctypes import windll, pointer
from ctypes.wintypes import DWORD
from six.moves import range

from prompt_toolkit.key_binding.input_processor import KeyPress
from prompt_toolkit.keys import Keys
from prompt_toolkit.mouse_events import MouseEventType
from prompt_toolkit.win32_types import EventTypes, KEY_EVENT_RECORD, MOUSE_EVENT_RECORD, INPUT_RECORD, STD_INPUT_HANDLE

import msvcrt
import os
import sys
import six

__all__ = (
    'ConsoleInputReader',
    'raw_mode',
    'cooked_mode'
)


class ConsoleInputReader(object):
    """
    :param recognize_paste: When True, try to discover paste actions and turn
        the event into a BracketedPaste.
    """
    # Keys with character data.
    mappings = {
        b'\x1b': Keys.Escape,

        b'\x00': Keys.ControlSpace,  # Control-Space (Also for Ctrl-@)
        b'\x01': Keys.ControlA,  # Control-A (home)
        b'\x02': Keys.ControlB,  # Control-B (emacs cursor left)
        b'\x03': Keys.ControlC,  # Control-C (interrupt)
        b'\x04': Keys.ControlD,  # Control-D (exit)
        b'\x05': Keys.ControlE,  # Contrel-E (end)
        b'\x06': Keys.ControlF,  # Control-F (cursor forward)
        b'\x07': Keys.ControlG,  # Control-G
        b'\x08': Keys.ControlH,  # Control-H (8) (Identical to '\b')
        b'\x09': Keys.ControlI,  # Control-I (9) (Identical to '\t')
        b'\x0a': Keys.ControlJ,  # Control-J (10) (Identical to '\n')
        b'\x0b': Keys.ControlK,  # Control-K (delete until end of line; vertical tab)
        b'\x0c': Keys.ControlL,  # Control-L (clear; form feed)
        b'\x0d': Keys.ControlJ,  # Control-J NOTE: Windows sends \r instead of
                                 #   \n when pressing enter. We turn it into \n
                                 #   to be compatible with other platforms.
        b'\x0e': Keys.ControlN,  # Control-N (14) (history forward)
        b'\x0f': Keys.ControlO,  # Control-O (15)
        b'\x10': Keys.ControlP,  # Control-P (16) (history back)
        b'\x11': Keys.ControlQ,  # Control-Q
        b'\x12': Keys.ControlR,  # Control-R (18) (reverse search)
        b'\x13': Keys.ControlS,  # Control-S (19) (forward search)
        b'\x14': Keys.ControlT,  # Control-T
        b'\x15': Keys.ControlU,  # Control-U
        b'\x16': Keys.ControlV,  # Control-V
        b'\x17': Keys.ControlW,  # Control-W
        b'\x18': Keys.ControlX,  # Control-X
        b'\x19': Keys.ControlY,  # Control-Y (25)
        b'\x1a': Keys.ControlZ,  # Control-Z

        b'\x1c': Keys.ControlBackslash,  # Both Control-\ and Ctrl-|
        b'\x1d': Keys.ControlSquareClose,  # Control-]
        b'\x1e': Keys.ControlCircumflex,  # Control-^
        b'\x1f': Keys.ControlUnderscore,  # Control-underscore (Also for Ctrl-hypen.)
        b'\x7f': Keys.Backspace,  # (127) Backspace
    }

    # Keys that don't carry character data.
    keycodes = {
        # Home/End
        33: Keys.PageUp,
        34: Keys.PageDown,
        35: Keys.End,
        36: Keys.Home,

        # Arrows
        37: Keys.Left,
        38: Keys.Up,
        39: Keys.Right,
        40: Keys.Down,

        45: Keys.Insert,
        46: Keys.Delete,

        # F-keys.
        112: Keys.F1,
        113: Keys.F2,
        114: Keys.F3,
        115: Keys.F4,
        116: Keys.F5,
        117: Keys.F6,
        118: Keys.F7,
        119: Keys.F8,
        120: Keys.F9,
        121: Keys.F10,
        122: Keys.F11,
        123: Keys.F12,
    }

    LEFT_ALT_PRESSED = 0x0002
    RIGHT_ALT_PRESSED = 0x0001
    SHIFT_PRESSED = 0x0010
    LEFT_CTRL_PRESSED = 0x0008
    RIGHT_CTRL_PRESSED = 0x0004

    def __init__(self, recognize_paste=True):
        self._fdcon = None
        self.recognize_paste = recognize_paste

        # When stdin is a tty, use that handle, otherwise, create a handle from
        # CONIN$.
        if sys.stdin.isatty():
            self.handle = windll.kernel32.GetStdHandle(STD_INPUT_HANDLE)
        else:
            self._fdcon = os.open('CONIN$', os.O_RDWR | os.O_BINARY)
            self.handle = msvcrt.get_osfhandle(self._fdcon)

    def close(self):
        " Close fdcon. "
        if self._fdcon is not None:
            os.close(self._fdcon)

    def read(self):
        """
        Return a list of `KeyPress` instances. It won't return anything when
        there was nothing to read.  (This function doesn't block.)

        http://msdn.microsoft.com/en-us/library/windows/desktop/ms684961(v=vs.85).aspx
        """
        max_count = 2048  # Max events to read at the same time.

        read = DWORD(0)
        arrtype = INPUT_RECORD * max_count
        input_records = arrtype()

        # Get next batch of input event.
        windll.kernel32.ReadConsoleInputW(
            self.handle, pointer(input_records), max_count, pointer(read))

        # First, get all the keys from the input buffer, in order to determine
        # whether we should consider this a paste event or not.
        all_keys = list(self._get_keys(read, input_records))

        if self.recognize_paste and self._is_paste(all_keys):
            gen = iter(all_keys)
            for k in gen:
                # Pasting: if the current key consists of text or \n, turn it
                # into a BracketedPaste.
                data = []
                while k and (isinstance(k.key, six.text_type) or
                             k.key == Keys.ControlJ):
                    data.append(k.data)
                    try:
                        k = next(gen)
                    except StopIteration:
                        k = None

                if data:
                    yield KeyPress(Keys.BracketedPaste, ''.join(data))
                if k is not None:
                    yield k
        else:
            for k in all_keys:
                yield k

    def _get_keys(self, read, input_records):
        """
        Generator that yields `KeyPress` objects from the input records.
        """
        for i in range(read.value):
            ir = input_records[i]

            # Get the right EventType from the EVENT_RECORD.
            # (For some reason the Windows console application 'cmder'
            # [http://gooseberrycreative.com/cmder/] can return '0' for
            # ir.EventType. -- Just ignore that.)
            if ir.EventType in EventTypes:
                ev = getattr(ir.Event, EventTypes[ir.EventType])

                # Process if this is a key event. (We also have mouse, menu and
                # focus events.)
                if type(ev) == KEY_EVENT_RECORD and ev.KeyDown:
                    for key_press in self._event_to_key_presses(ev):
                        yield key_press

                elif type(ev) == MOUSE_EVENT_RECORD:
                    for key_press in self._handle_mouse(ev):
                        yield key_press

    @staticmethod
    def _is_paste(keys):
        """
        Return `True` when we should consider this list of keys as a paste
        event. Pasted text on windows will be turned into a
        `Keys.BracketedPaste` event. (It's not 100% correct, but it is probably
        the best possible way to detect pasting of text and handle that
        correctly.)
        """
        # Consider paste when it contains at least one newline and at least one
        # other character.
        text_count = 0
        newline_count = 0

        for k in keys:
            if isinstance(k.key, six.text_type):
                text_count += 1
            if k.key == Keys.ControlJ:
                newline_count += 1

        return newline_count >= 1 and text_count > 1

    def _event_to_key_presses(self, ev):
        """
        For this `KEY_EVENT_RECORD`, return a list of `KeyPress` instances.
        """
        assert type(ev) == KEY_EVENT_RECORD and ev.KeyDown

        result = None

        u_char = ev.uChar.UnicodeChar
        ascii_char = u_char.encode('utf-8')

        # NOTE: We don't use `ev.uChar.AsciiChar`. That appears to be latin-1
        #       encoded. See also:
        # https://github.com/ipython/ipython/issues/10004
        # https://github.com/jonathanslenders/python-prompt-toolkit/issues/389

        if u_char == '\x00':
            if ev.VirtualKeyCode in self.keycodes:
                result = KeyPress(self.keycodes[ev.VirtualKeyCode], '')
        else:
            if ascii_char in self.mappings:
                if self.mappings[ascii_char] == Keys.ControlJ:
                    u_char = '\n'  # Windows sends \n, turn into \r for unix compatibility.
                result = KeyPress(self.mappings[ascii_char], u_char)
            else:
                result = KeyPress(u_char, u_char)

        # Correctly handle Control-Arrow keys.
        if (ev.ControlKeyState & self.LEFT_CTRL_PRESSED or
                ev.ControlKeyState & self.RIGHT_CTRL_PRESSED) and result:
            if result.key == Keys.Left:
                result.key = Keys.ControlLeft

            if result.key == Keys.Right:
                result.key = Keys.ControlRight

            if result.key == Keys.Up:
                result.key = Keys.ControlUp

            if result.key == Keys.Down:
                result.key = Keys.ControlDown

        # Turn 'Tab' into 'BackTab' when shift was pressed.
        if ev.ControlKeyState & self.SHIFT_PRESSED and result:
            if result.key == Keys.Tab:
                result.key = Keys.BackTab

        # Turn 'Space' into 'ControlSpace' when control was pressed.
        if (ev.ControlKeyState & self.LEFT_CTRL_PRESSED or
                ev.ControlKeyState & self.RIGHT_CTRL_PRESSED) and result and result.data == ' ':
            result = KeyPress(Keys.ControlSpace, ' ')

        # Turn Control-Enter into META-Enter. (On a vt100 terminal, we cannot
        # detect this combination. But it's really practical on Windows.)
        if (ev.ControlKeyState & self.LEFT_CTRL_PRESSED or
                ev.ControlKeyState & self.RIGHT_CTRL_PRESSED) and result and \
                result.key == Keys.ControlJ:
            return [KeyPress(Keys.Escape, ''), result]

        # Return result. If alt was pressed, prefix the result with an
        # 'Escape' key, just like unix VT100 terminals do.

        # NOTE: Only replace the left alt with escape. The right alt key often
        #       acts as altgr and is used in many non US keyboard layouts for
        #       typing some special characters, like a backslash. We don't want
        #       all backslashes to be prefixed with escape. (Esc-\ has a
        #       meaning in E-macs, for instance.)
        if result:
            meta_pressed = ev.ControlKeyState & self.LEFT_ALT_PRESSED

            if meta_pressed:
                return [KeyPress(Keys.Escape, ''), result]
            else:
                return [result]

        else:
            return []

    def _handle_mouse(self, ev):
        """
        Handle mouse events. Return a list of KeyPress instances.
        """
        FROM_LEFT_1ST_BUTTON_PRESSED = 0x1

        result = []

        # Check event type.
        if ev.ButtonState == FROM_LEFT_1ST_BUTTON_PRESSED:
            # On a key press, generate both the mouse down and up event.
            for event_type in [MouseEventType.MOUSE_DOWN, MouseEventType.MOUSE_UP]:
                data = ';'.join([
                   event_type,
                   str(ev.MousePosition.X),
                   str(ev.MousePosition.Y)
                ])
                result.append(KeyPress(Keys.WindowsMouseEvent, data))

        return result


class raw_mode(object):
    """
    ::

        with raw_mode(stdin):
            ''' the windows terminal is now in 'raw' mode. '''

    The ``fileno`` attribute is ignored. This is to be compatble with the
    `raw_input` method of `.vt100_input`.
    """
    def __init__(self, fileno=None):
        self.handle = windll.kernel32.GetStdHandle(STD_INPUT_HANDLE)

    def __enter__(self):
        # Remember original mode.
        original_mode = DWORD()
        windll.kernel32.GetConsoleMode(self.handle, pointer(original_mode))
        self.original_mode = original_mode

        self._patch()

    def _patch(self):
        # Set raw
        ENABLE_ECHO_INPUT = 0x0004
        ENABLE_LINE_INPUT = 0x0002
        ENABLE_PROCESSED_INPUT = 0x0001

        windll.kernel32.SetConsoleMode(
            self.handle, self.original_mode.value &
            ~(ENABLE_ECHO_INPUT | ENABLE_LINE_INPUT | ENABLE_PROCESSED_INPUT))

    def __exit__(self, *a, **kw):
        # Restore original mode
        windll.kernel32.SetConsoleMode(self.handle, self.original_mode)


class cooked_mode(raw_mode):
    """
    ::

        with cooked_mode(stdin):
            ''' the pseudo-terminal stdin is now used in raw mode '''
    """
    def _patch(self):
        # Set cooked.
        ENABLE_ECHO_INPUT = 0x0004
        ENABLE_LINE_INPUT = 0x0002
        ENABLE_PROCESSED_INPUT = 0x0001

        windll.kernel32.SetConsoleMode(
            self.handle, self.original_mode.value |
            (ENABLE_ECHO_INPUT | ENABLE_LINE_INPUT | ENABLE_PROCESSED_INPUT))
