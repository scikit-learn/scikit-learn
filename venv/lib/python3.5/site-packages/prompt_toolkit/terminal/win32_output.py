from __future__ import unicode_literals

from ctypes import windll, byref, ArgumentError, c_char, c_long, c_ulong, c_uint, pointer
from ctypes.wintypes import DWORD

from prompt_toolkit.renderer import Output
from prompt_toolkit.styles import ANSI_COLOR_NAMES
from prompt_toolkit.win32_types import CONSOLE_SCREEN_BUFFER_INFO, STD_OUTPUT_HANDLE, STD_INPUT_HANDLE, COORD, SMALL_RECT

import os
import six

__all__ = (
    'Win32Output',
)


def _coord_byval(coord):
    """
    Turns a COORD object into a c_long.
    This will cause it to be passed by value instead of by reference. (That is what I think at least.)

    When runing ``ptipython`` is run (only with IPython), we often got the following error::

         Error in 'SetConsoleCursorPosition'.
         ArgumentError("argument 2: <class 'TypeError'>: wrong type",)
     argument 2: <class 'TypeError'>: wrong type

    It was solved by turning ``COORD`` parameters into a ``c_long`` like this.

    More info: http://msdn.microsoft.com/en-us/library/windows/desktop/ms686025(v=vs.85).aspx
    """
    return c_long(coord.Y * 0x10000 | coord.X & 0xFFFF)


#: If True: write the output of the renderer also to the following file. This
#: is very useful for debugging. (e.g.: to see that we don't write more bytes
#: than required.)
_DEBUG_RENDER_OUTPUT = False
_DEBUG_RENDER_OUTPUT_FILENAME = r'prompt-toolkit-windows-output.log'


class NoConsoleScreenBufferError(Exception):
    """
    Raised when the application is not running inside a Windows Console, but
    the user tries to instantiate Win32Output.
    """
    def __init__(self):
        # Are we running in 'xterm' on Windows, like git-bash for instance?
        xterm = 'xterm' in os.environ.get('TERM', '')

        if xterm:
            message = ('Found %s, while expecting a Windows console. '
                       'Maybe try to run this program using "winpty" '
                       'or run it in cmd.exe instead. Or otherwise, '
                       'in case of Cygwin, use the Python executable '
                       'that is compiled for Cygwin.' % os.environ['TERM'])
        else:
            message = 'No Windows console found. Are you running cmd.exe?'
        super(NoConsoleScreenBufferError, self).__init__(message)


class Win32Output(Output):
    """
    I/O abstraction for rendering to Windows consoles.
    (cmd.exe and similar.)
    """
    def __init__(self, stdout, use_complete_width=False):
        self.use_complete_width = use_complete_width

        self._buffer = []
        self.stdout = stdout
        self.hconsole = windll.kernel32.GetStdHandle(STD_OUTPUT_HANDLE)

        self._in_alternate_screen = False

        self.color_lookup_table = ColorLookupTable()

        # Remember the default console colors.
        info = self.get_win32_screen_buffer_info()
        self.default_attrs = info.wAttributes if info else 15

        if _DEBUG_RENDER_OUTPUT:
            self.LOG = open(_DEBUG_RENDER_OUTPUT_FILENAME, 'ab')

    def fileno(self):
        " Return file descriptor. "
        return self.stdout.fileno()

    def encoding(self):
        " Return encoding used for stdout. "
        return self.stdout.encoding

    def write(self, data):
        self._buffer.append(data)

    def write_raw(self, data):
        " For win32, there is no difference between write and write_raw. "
        self.write(data)

    def get_size(self):
        from prompt_toolkit.layout.screen import Size
        info = self.get_win32_screen_buffer_info()

        # We take the width of the *visible* region as the size. Not the width
        # of the complete screen buffer. (Unless use_complete_width has been
        # set.)
        if self.use_complete_width:
            width = info.dwSize.X
        else:
            width = info.srWindow.Right - info.srWindow.Left

        height = info.srWindow.Bottom - info.srWindow.Top + 1

        # We avoid the right margin, windows will wrap otherwise.
        maxwidth = info.dwSize.X - 1
        width = min(maxwidth, width)

        # Create `Size` object.
        return Size(rows=height, columns=width)

    def _winapi(self, func, *a, **kw):
        """
        Flush and call win API function.
        """
        self.flush()

        if _DEBUG_RENDER_OUTPUT:
            self.LOG.write(('%r' % func.__name__).encode('utf-8') + b'\n')
            self.LOG.write(b'     ' + ', '.join(['%r' % i for i in a]).encode('utf-8') + b'\n')
            self.LOG.write(b'     ' + ', '.join(['%r' % type(i) for i in a]).encode('utf-8') + b'\n')
            self.LOG.flush()

        try:
            return func(*a, **kw)
        except ArgumentError as e:
            if _DEBUG_RENDER_OUTPUT:
                self.LOG.write(('    Error in %r %r %s\n' % (func.__name__, e, e)).encode('utf-8'))

    def get_win32_screen_buffer_info(self):
        """
        Return Screen buffer info.
        """
        # NOTE: We don't call the `GetConsoleScreenBufferInfo` API through
        #     `self._winapi`. Doing so causes Python to crash on certain 64bit
        #     Python versions. (Reproduced with 64bit Python 2.7.6, on Windows
        #     10). It is not clear why. Possibly, it has to do with passing
        #     these objects as an argument, or through *args.

        # The Python documentation contains the following - possibly related - warning:
        #     ctypes does not support passing unions or structures with
        #     bit-fields to functions by value. While this may work on 32-bit
        #     x86, it's not guaranteed by the library to work in the general
        #     case. Unions and structures with bit-fields should always be
        #     passed to functions by pointer.

        # Also see:
        #    - https://github.com/ipython/ipython/issues/10070
        #    - https://github.com/jonathanslenders/python-prompt-toolkit/issues/406
        #    - https://github.com/jonathanslenders/python-prompt-toolkit/issues/86

        self.flush()
        sbinfo = CONSOLE_SCREEN_BUFFER_INFO()
        success = windll.kernel32.GetConsoleScreenBufferInfo(self.hconsole, byref(sbinfo))

        # success = self._winapi(windll.kernel32.GetConsoleScreenBufferInfo,
        #                        self.hconsole, byref(sbinfo))

        if success:
            return sbinfo
        else:
            raise NoConsoleScreenBufferError

    def set_title(self, title):
        """
        Set terminal title.
        """
        assert isinstance(title, six.text_type)
        self._winapi(windll.kernel32.SetConsoleTitleW, title)

    def clear_title(self):
        self._winapi(windll.kernel32.SetConsoleTitleW, '')

    def erase_screen(self):
        start = COORD(0, 0)
        sbinfo = self.get_win32_screen_buffer_info()
        length = sbinfo.dwSize.X * sbinfo.dwSize.Y

        self.cursor_goto(row=0, column=0)
        self._erase(start, length)

    def erase_down(self):
        sbinfo = self.get_win32_screen_buffer_info()
        size = sbinfo.dwSize

        start = sbinfo.dwCursorPosition
        length = ((size.X - size.X) + size.X * (size.Y - sbinfo.dwCursorPosition.Y))

        self._erase(start, length)

    def erase_end_of_line(self):
        """
        """
        sbinfo = self.get_win32_screen_buffer_info()
        start = sbinfo.dwCursorPosition
        length = sbinfo.dwSize.X - sbinfo.dwCursorPosition.X

        self._erase(start, length)

    def _erase(self, start, length):
        chars_written = c_ulong()

        self._winapi(windll.kernel32.FillConsoleOutputCharacterA,
                     self.hconsole, c_char(b' '), DWORD(length), _coord_byval(start),
                     byref(chars_written))

        # Reset attributes.
        sbinfo = self.get_win32_screen_buffer_info()
        self._winapi(windll.kernel32.FillConsoleOutputAttribute,
                     self.hconsole, sbinfo.wAttributes, length, _coord_byval(start),
                     byref(chars_written))

    def reset_attributes(self):
        " Reset the console foreground/background color. "
        self._winapi(windll.kernel32.SetConsoleTextAttribute, self.hconsole,
                     self.default_attrs)

    def set_attributes(self, attrs):
        fgcolor, bgcolor, bold, underline, italic, blink, reverse = attrs

        # Start from the default attributes.
        attrs = self.default_attrs

        # Override the last four bits: foreground color.
        if fgcolor is not None:
            attrs = attrs & ~0xf
            attrs |= self.color_lookup_table.lookup_fg_color(fgcolor)

        # Override the next four bits: background color.
        if bgcolor is not None:
            attrs = attrs & ~0xf0
            attrs |= self.color_lookup_table.lookup_bg_color(bgcolor)

        # Reverse: swap these four bits groups.
        if reverse:
            attrs = (attrs & ~0xff) | ((attrs & 0xf) << 4) | ((attrs & 0xf0) >> 4)

        self._winapi(windll.kernel32.SetConsoleTextAttribute, self.hconsole, attrs)

    def disable_autowrap(self):
        # Not supported by Windows.
        pass

    def enable_autowrap(self):
        # Not supported by Windows.
        pass

    def cursor_goto(self, row=0, column=0):
        pos = COORD(x=column, y=row)
        self._winapi(windll.kernel32.SetConsoleCursorPosition, self.hconsole, _coord_byval(pos))

    def cursor_up(self, amount):
        sr = self.get_win32_screen_buffer_info().dwCursorPosition
        pos = COORD(sr.X, sr.Y - amount)
        self._winapi(windll.kernel32.SetConsoleCursorPosition, self.hconsole, _coord_byval(pos))

    def cursor_down(self, amount):
        self.cursor_up(-amount)

    def cursor_forward(self, amount):
        sr = self.get_win32_screen_buffer_info().dwCursorPosition
#        assert sr.X + amount >= 0, 'Negative cursor position: x=%r amount=%r' % (sr.X, amount)

        pos = COORD(max(0, sr.X + amount), sr.Y)
        self._winapi(windll.kernel32.SetConsoleCursorPosition, self.hconsole, _coord_byval(pos))

    def cursor_backward(self, amount):
        self.cursor_forward(-amount)

    def flush(self):
        """
        Write to output stream and flush.
        """
        if not self._buffer:
            # Only flush stdout buffer. (It could be that Python still has
            # something in its buffer. -- We want to be sure to print that in
            # the correct color.)
            self.stdout.flush()
            return

        data = ''.join(self._buffer)

        if _DEBUG_RENDER_OUTPUT:
            self.LOG.write(('%r' % data).encode('utf-8') + b'\n')
            self.LOG.flush()

        # Print characters one by one. This appears to be the best soluton
        # in oder to avoid traces of vertical lines when the completion
        # menu disappears.
        for b in data:
            written = DWORD()

            retval = windll.kernel32.WriteConsoleW(self.hconsole, b, 1, byref(written), None)
            assert retval != 0

        self._buffer = []

    def get_rows_below_cursor_position(self):
        info = self.get_win32_screen_buffer_info()
        return info.srWindow.Bottom - info.dwCursorPosition.Y + 1

    def scroll_buffer_to_prompt(self):
        """
        To be called before drawing the prompt. This should scroll the console
        to left, with the cursor at the bottom (if possible).
        """
        # Get current window size
        info = self.get_win32_screen_buffer_info()
        sr = info.srWindow
        cursor_pos = info.dwCursorPosition

        result = SMALL_RECT()

        # Scroll to the left.
        result.Left = 0
        result.Right = sr.Right - sr.Left

        # Scroll vertical
        win_height = sr.Bottom - sr.Top
        if 0 < sr.Bottom - cursor_pos.Y < win_height - 1:
            # no vertical scroll if cursor already on the screen
            result.Bottom = sr.Bottom
        else:
            result.Bottom = max(win_height, cursor_pos.Y)
        result.Top = result.Bottom - win_height

        # Scroll API
        self._winapi(windll.kernel32.SetConsoleWindowInfo, self.hconsole, True, byref(result))

    def enter_alternate_screen(self):
        """
        Go to alternate screen buffer.
        """
        if not self._in_alternate_screen:
            GENERIC_READ = 0x80000000
            GENERIC_WRITE = 0x40000000

            # Create a new console buffer and activate that one.
            handle = self._winapi(windll.kernel32.CreateConsoleScreenBuffer, GENERIC_READ|GENERIC_WRITE,
                                  DWORD(0), None, DWORD(1), None)

            self._winapi(windll.kernel32.SetConsoleActiveScreenBuffer, handle)
            self.hconsole = handle
            self._in_alternate_screen = True

    def quit_alternate_screen(self):
        """
        Make stdout again the active buffer.
        """
        if self._in_alternate_screen:
            stdout = self._winapi(windll.kernel32.GetStdHandle, STD_OUTPUT_HANDLE)
            self._winapi(windll.kernel32.SetConsoleActiveScreenBuffer, stdout)
            self._winapi(windll.kernel32.CloseHandle, self.hconsole)
            self.hconsole = stdout
            self._in_alternate_screen = False

    def enable_mouse_support(self):
        ENABLE_MOUSE_INPUT = 0x10
        handle = windll.kernel32.GetStdHandle(STD_INPUT_HANDLE)

        original_mode = DWORD()
        self._winapi(windll.kernel32.GetConsoleMode, handle, pointer(original_mode))
        self._winapi(windll.kernel32.SetConsoleMode, handle, original_mode.value | ENABLE_MOUSE_INPUT)

    def disable_mouse_support(self):
        ENABLE_MOUSE_INPUT = 0x10
        handle = windll.kernel32.GetStdHandle(STD_INPUT_HANDLE)

        original_mode = DWORD()
        self._winapi(windll.kernel32.GetConsoleMode, handle, pointer(original_mode))
        self._winapi(windll.kernel32.SetConsoleMode, handle, original_mode.value & ~ ENABLE_MOUSE_INPUT)

    def hide_cursor(self):
        pass

    def show_cursor(self):
        pass

    @classmethod
    def win32_refresh_window(cls):
        """
        Call win32 API to refresh the whole Window.

        This is sometimes necessary when the application paints background
        for completion menus. When the menu disappears, it leaves traces due
        to a bug in the Windows Console. Sending a repaint request solves it.
        """
        # Get console handle
        handle = windll.kernel32.GetConsoleWindow()

        RDW_INVALIDATE = 0x0001
        windll.user32.RedrawWindow(handle, None, None, c_uint(RDW_INVALIDATE))


class FOREGROUND_COLOR:
    BLACK     = 0x0000
    BLUE      = 0x0001
    GREEN     = 0x0002
    CYAN      = 0x0003
    RED       = 0x0004
    MAGENTA   = 0x0005
    YELLOW    = 0x0006
    GRAY      = 0x0007
    INTENSITY = 0x0008  # Foreground color is intensified.


class BACKROUND_COLOR:
    BLACK     = 0x0000
    BLUE      = 0x0010
    GREEN     = 0x0020
    CYAN      = 0x0030
    RED       = 0x0040
    MAGENTA   = 0x0050
    YELLOW    = 0x0060
    GRAY      = 0x0070
    INTENSITY = 0x0080  # Background color is intensified.


def _create_ansi_color_dict(color_cls):
    " Create a table that maps the 16 named ansi colors to their Windows code. "
    return {
        'ansidefault':   color_cls.BLACK,
        'ansiblack':     color_cls.BLACK,
        'ansidarkgray':  color_cls.BLACK | color_cls.INTENSITY,
        'ansilightgray': color_cls.GRAY,
        'ansiwhite':     color_cls.GRAY | color_cls.INTENSITY,
        
        # Low intensity.
        'ansidarkred':     color_cls.RED,
        'ansidarkgreen':   color_cls.GREEN,
        'ansibrown':       color_cls.YELLOW,
        'ansidarkblue':    color_cls.BLUE,
        'ansipurple':      color_cls.MAGENTA,
        'ansiteal':        color_cls.CYAN,

        # High intensity.
        'ansired':        color_cls.RED | color_cls.INTENSITY,
        'ansigreen':      color_cls.GREEN | color_cls.INTENSITY,
        'ansiyellow':     color_cls.YELLOW | color_cls.INTENSITY,
        'ansiblue':       color_cls.BLUE | color_cls.INTENSITY,
        'ansifuchsia':    color_cls.MAGENTA | color_cls.INTENSITY,
        'ansiturquoise':  color_cls.CYAN | color_cls.INTENSITY,
    }

FG_ANSI_COLORS = _create_ansi_color_dict(FOREGROUND_COLOR)
BG_ANSI_COLORS = _create_ansi_color_dict(BACKROUND_COLOR)

assert set(FG_ANSI_COLORS) == set(ANSI_COLOR_NAMES)
assert set(BG_ANSI_COLORS) == set(ANSI_COLOR_NAMES)


class ColorLookupTable(object):
    """
    Inspired by pygments/formatters/terminal256.py
    """
    def __init__(self):
        self._win32_colors = self._build_color_table()
        self.best_match = {}  # Cache

    @staticmethod
    def _build_color_table():
        """
        Build an RGB-to-256 color conversion table
        """
        FG = FOREGROUND_COLOR
        BG = BACKROUND_COLOR

        return [
            (0x00, 0x00, 0x00, FG.BLACK, BG.BLACK),
            (0x00, 0x00, 0xaa, FG.BLUE, BG.BLUE),
            (0x00, 0xaa, 0x00, FG.GREEN, BG.GREEN),
            (0x00, 0xaa, 0xaa, FG.CYAN, BG.CYAN),
            (0xaa, 0x00, 0x00, FG.RED, BG.RED),
            (0xaa, 0x00, 0xaa, FG.MAGENTA, BG.MAGENTA),
            (0xaa, 0xaa, 0x00, FG.YELLOW, BG.YELLOW),
            (0x88, 0x88, 0x88, FG.GRAY, BG.GRAY),

            (0x44, 0x44, 0xff, FG.BLUE | FG.INTENSITY, BG.BLUE | BG.INTENSITY),
            (0x44, 0xff, 0x44, FG.GREEN | FG.INTENSITY, BG.GREEN | BG.INTENSITY),
            (0x44, 0xff, 0xff, FG.CYAN | FG.INTENSITY, BG.CYAN | BG.INTENSITY),
            (0xff, 0x44, 0x44, FG.RED | FG.INTENSITY, BG.RED | BG.INTENSITY),
            (0xff, 0x44, 0xff, FG.MAGENTA | FG.INTENSITY, BG.MAGENTA | BG.INTENSITY),
            (0xff, 0xff, 0x44, FG.YELLOW | FG.INTENSITY, BG.YELLOW | BG.INTENSITY),

            (0x44, 0x44, 0x44, FG.BLACK | FG.INTENSITY, BG.BLACK | BG.INTENSITY),
            (0xff, 0xff, 0xff, FG.GRAY | FG.INTENSITY, BG.GRAY | BG.INTENSITY),
        ]

    def _closest_color(self, r, g, b):
        distance = 257 * 257 * 3  # "infinity" (>distance from #000000 to #ffffff)
        fg_match = 0
        bg_match = 0

        for r_, g_, b_, fg_, bg_ in self._win32_colors:
            rd = r - r_
            gd = g - g_
            bd = b - b_

            d = rd * rd + gd * gd + bd * bd

            if d < distance:
                fg_match = fg_
                bg_match = bg_
                distance = d
        return fg_match, bg_match

    def _color_indexes(self, color):
        indexes = self.best_match.get(color, None)
        if indexes is None:
            try:
                rgb = int(str(color), 16)
            except ValueError:
                rgb = 0

            r = (rgb >> 16) & 0xff
            g = (rgb >> 8) & 0xff
            b = rgb & 0xff
            indexes = self._closest_color(r, g, b)
            self.best_match[color] = indexes
        return indexes

    def lookup_fg_color(self, fg_color):
        """
        Return the color for use in the
        `windll.kernel32.SetConsoleTextAttribute` API call.

        :param fg_color: Foreground as text. E.g. 'ffffff' or 'red'
        """
        # Foreground.
        if fg_color in FG_ANSI_COLORS:
            return FG_ANSI_COLORS[fg_color]
        else:
            return self._color_indexes(fg_color)[0]

    def lookup_bg_color(self, bg_color):
        """
        Return the color for use in the
        `windll.kernel32.SetConsoleTextAttribute` API call.

        :param bg_color: Background as text. E.g. 'ffffff' or 'red'
        """
        # Background.
        if bg_color in BG_ANSI_COLORS:
            return BG_ANSI_COLORS[bg_color]
        else:
            return self._color_indexes(bg_color)[1]
