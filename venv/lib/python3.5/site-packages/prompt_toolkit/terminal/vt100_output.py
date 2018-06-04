"""
Output for vt100 terminals.

A lot of thanks, regarding outputting of colors, goes to the Pygments project:
(We don't rely on Pygments anymore, because many things are very custom, and
everything has been highly optimized.)
http://pygments.org/
"""
from __future__ import unicode_literals

from prompt_toolkit.filters import to_simple_filter, Condition
from prompt_toolkit.layout.screen import Size
from prompt_toolkit.renderer import Output
from prompt_toolkit.styles import ANSI_COLOR_NAMES

from six.moves import range
import array
import errno
import os
import six

__all__ = (
    'Vt100_Output',
)


FG_ANSI_COLORS = {
    'ansidefault': 39,

    # Low intensity.
    'ansiblack':       30,
    'ansidarkred':     31,
    'ansidarkgreen':   32,
    'ansibrown':       33,
    'ansidarkblue':    34,
    'ansipurple':      35,
    'ansiteal':        36,
    'ansilightgray':   37,

    # High intensity.
    'ansidarkgray':    90,
    'ansired':         91,
    'ansigreen':       92,
    'ansiyellow':      93,
    'ansiblue':        94,
    'ansifuchsia':     95,
    'ansiturquoise':   96,
    'ansiwhite':       97,
}

BG_ANSI_COLORS = {
    'ansidefault':     49,

    # Low intensity.
    'ansiblack':       40,
    'ansidarkred':     41,
    'ansidarkgreen':   42,
    'ansibrown':       43,
    'ansidarkblue':    44,
    'ansipurple':      45,
    'ansiteal':        46,
    'ansilightgray':   47,

    # High intensity.
    'ansidarkgray':    100,
    'ansired':         101,
    'ansigreen':       102,
    'ansiyellow':      103,
    'ansiblue':        104,
    'ansifuchsia':     105,
    'ansiturquoise':   106,
    'ansiwhite':       107,
}


ANSI_COLORS_TO_RGB = {
    'ansidefault':   (0x00, 0x00, 0x00),  # Don't use, 'default' doesn't really have a value.
    'ansiblack':     (0x00, 0x00, 0x00),
    'ansidarkgray':  (0x7f, 0x7f, 0x7f), 
    'ansiwhite':     (0xff, 0xff, 0xff),
    'ansilightgray': (0xe5, 0xe5, 0xe5),

    # Low intensity.
    'ansidarkred':     (0xcd, 0x00, 0x00),
    'ansidarkgreen':   (0x00, 0xcd, 0x00),
    'ansibrown':       (0xcd, 0xcd, 0x00),
    'ansidarkblue':    (0x00, 0x00, 0xcd),
    'ansipurple':      (0xcd, 0x00, 0xcd),
    'ansiteal':        (0x00, 0xcd, 0xcd),

    # High intensity.
    'ansired':         (0xff, 0x00, 0x00),
    'ansigreen':       (0x00, 0xff, 0x00),
    'ansiyellow':      (0xff, 0xff, 0x00),
    'ansiblue':        (0x00, 0x00, 0xff),
    'ansifuchsia':     (0xff, 0x00, 0xff),
    'ansiturquoise':   (0x00, 0xff, 0xff),
}


assert set(FG_ANSI_COLORS) == set(ANSI_COLOR_NAMES)
assert set(BG_ANSI_COLORS) == set(ANSI_COLOR_NAMES)
assert set(ANSI_COLORS_TO_RGB) == set(ANSI_COLOR_NAMES)


def _get_closest_ansi_color(r, g, b, exclude=()):
    """
    Find closest ANSI color. Return it by name.

    :param r: Red (Between 0 and 255.)
    :param g: Green (Between 0 and 255.)
    :param b: Blue (Between 0 and 255.)
    :param exclude: A tuple of color names to exclude. (E.g. ``('ansired', )``.)
    """
    assert isinstance(exclude, tuple)

    # When we have a bit of saturation, avoid the gray-like colors, otherwise,
    # too often the distance to the gray color is less.
    saturation = abs(r - g) + abs(g - b) + abs(b - r)  # Between 0..510

    if saturation > 30:
        exclude += ('ansilightgray', 'ansidarkgray', 'ansiwhite', 'ansiblack')

    # Take the closest color.
    # (Thanks to Pygments for this part.)
    distance = 257*257*3  # "infinity" (>distance from #000000 to #ffffff)
    match = 'ansidefault'

    for name, (r2, g2, b2) in ANSI_COLORS_TO_RGB.items():
        if name != 'ansidefault' and name not in exclude:
            d = (r - r2) ** 2 + (g - g2) ** 2 + (b - b2) ** 2

            if d < distance:
                match = name
                distance = d

    return match


class _16ColorCache(dict):
    """
    Cache which maps (r, g, b) tuples to 16 ansi colors.

    :param bg: Cache for background colors, instead of foreground.
    """
    def __init__(self, bg=False):
        assert isinstance(bg, bool)
        self.bg = bg

    def get_code(self, value, exclude=()):
        """
        Return a (ansi_code, ansi_name) tuple. (E.g. ``(44, 'ansiblue')``.) for
        a given (r,g,b) value.
        """
        key = (value, exclude)
        if key not in self:
            self[key] = self._get(value, exclude)
        return self[key]

    def _get(self, value, exclude=()):
        r, g, b = value
        match = _get_closest_ansi_color(r, g, b, exclude=exclude)

        # Turn color name into code.
        if self.bg:
            code = BG_ANSI_COLORS[match]
        else:
            code = FG_ANSI_COLORS[match]

        self[value] = code
        return code, match


class _256ColorCache(dict):
    """
    Cach which maps (r, g, b) tuples to 256 colors.
    """
    def __init__(self):
        # Build color table.
        colors = []

        # colors 0..15: 16 basic colors
        colors.append((0x00, 0x00, 0x00))  # 0
        colors.append((0xcd, 0x00, 0x00))  # 1
        colors.append((0x00, 0xcd, 0x00))  # 2
        colors.append((0xcd, 0xcd, 0x00))  # 3
        colors.append((0x00, 0x00, 0xee))  # 4
        colors.append((0xcd, 0x00, 0xcd))  # 5
        colors.append((0x00, 0xcd, 0xcd))  # 6
        colors.append((0xe5, 0xe5, 0xe5))  # 7
        colors.append((0x7f, 0x7f, 0x7f))  # 8
        colors.append((0xff, 0x00, 0x00))  # 9
        colors.append((0x00, 0xff, 0x00))  # 10
        colors.append((0xff, 0xff, 0x00))  # 11
        colors.append((0x5c, 0x5c, 0xff))  # 12
        colors.append((0xff, 0x00, 0xff))  # 13
        colors.append((0x00, 0xff, 0xff))  # 14
        colors.append((0xff, 0xff, 0xff))  # 15

        # colors 16..232: the 6x6x6 color cube
        valuerange = (0x00, 0x5f, 0x87, 0xaf, 0xd7, 0xff)

        for i in range(217):
            r = valuerange[(i // 36) % 6]
            g = valuerange[(i // 6) % 6]
            b = valuerange[i % 6]
            colors.append((r, g, b))

        # colors 233..253: grayscale
        for i in range(1, 22):
            v = 8 + i * 10
            colors.append((v, v, v))

        self.colors = colors

    def __missing__(self, value):
        r, g, b = value

        # Find closest color.
        # (Thanks to Pygments for this!)
        distance = 257*257*3  # "infinity" (>distance from #000000 to #ffffff)
        match = 0

        for i, (r2, g2, b2) in enumerate(self.colors):
            d = (r - r2) ** 2 + (g - g2) ** 2 + (b - b2) ** 2

            if d < distance:
                match = i
                distance = d

        # Turn color name into code.
        self[value] = match
        return match


_16_fg_colors = _16ColorCache(bg=False)
_16_bg_colors = _16ColorCache(bg=True)
_256_colors = _256ColorCache()


class _EscapeCodeCache(dict):
    """
    Cache for VT100 escape codes. It maps
    (fgcolor, bgcolor, bold, underline, reverse) tuples to VT100 escape sequences.

    :param true_color: When True, use 24bit colors instead of 256 colors.
    """
    def __init__(self, true_color=False, ansi_colors_only=False):
        assert isinstance(true_color, bool)
        self.true_color = true_color
        self.ansi_colors_only = to_simple_filter(ansi_colors_only)

    def __missing__(self, attrs):
        fgcolor, bgcolor, bold, underline, italic, blink, reverse = attrs
        parts = []

        parts.extend(self._colors_to_code(fgcolor, bgcolor))

        if bold:
            parts.append('1')
        if italic:
            parts.append('3')
        if blink:
            parts.append('5')
        if underline:
            parts.append('4')
        if reverse:
            parts.append('7')

        if parts:
            result = '\x1b[0;' + ';'.join(parts) + 'm'
        else:
            result = '\x1b[0m'

        self[attrs] = result
        return result

    def _color_name_to_rgb(self, color):
        " Turn 'ffffff', into (0xff, 0xff, 0xff). "
        try:
            rgb = int(color, 16)
        except ValueError:
            raise
        else:
            r = (rgb >> 16) & 0xff
            g = (rgb >> 8) & 0xff
            b = rgb & 0xff
            return r, g, b

    def _colors_to_code(self, fg_color, bg_color):
        " Return a tuple with the vt100 values  that represent this color. "
        # When requesting ANSI colors only, and both fg/bg color were converted
        # to ANSI, ensure that the foreground and background color are not the
        # same. (Unless they were explicitely defined to be the same color.)
        fg_ansi = [()]

        def get(color, bg):
            table = BG_ANSI_COLORS if bg else FG_ANSI_COLORS

            if color is None:
                return ()

            # 16 ANSI colors. (Given by name.)
            elif color in table:
                return (table[color], )

            # RGB colors. (Defined as 'ffffff'.)
            else:
                try:
                    rgb = self._color_name_to_rgb(color)
                except ValueError:
                    return ()

                # When only 16 colors are supported, use that.
                if self.ansi_colors_only():
                    if bg:  # Background.
                        if fg_color != bg_color:
                            exclude = (fg_ansi[0], )
                        else:
                            exclude = ()
                        code, name = _16_bg_colors.get_code(rgb, exclude=exclude)
                        return (code, )
                    else:  # Foreground.
                        code, name = _16_fg_colors.get_code(rgb)
                        fg_ansi[0] = name
                        return (code, )

                # True colors. (Only when this feature is enabled.)
                elif self.true_color:
                    r, g, b = rgb
                    return (48 if bg else 38, 2, r, g, b)

                # 256 RGB colors.
                else:
                    return (48 if bg else 38, 5, _256_colors[rgb])

        result = []
        result.extend(get(fg_color, False))
        result.extend(get(bg_color, True))

        return map(six.text_type, result)


def _get_size(fileno):
    # Thanks to fabric (fabfile.org), and
    # http://sqizit.bartletts.id.au/2011/02/14/pseudo-terminals-in-python/
    """
    Get the size of this pseudo terminal.

    :param fileno: stdout.fileno()
    :returns: A (rows, cols) tuple.
    """
    # Inline imports, because these modules are not available on Windows.
    # (This file is used by ConEmuOutput, which is used on Windows.)
    import fcntl
    import termios

    # Buffer for the C call
    buf = array.array(b'h' if six.PY2 else u'h', [0, 0, 0, 0])

    # Do TIOCGWINSZ (Get)
    # Note: We should not pass 'True' as a fourth parameter to 'ioctl'. (True
    #       is the default.) This causes segmentation faults on some systems.
    #       See: https://github.com/jonathanslenders/python-prompt-toolkit/pull/364
    fcntl.ioctl(fileno, termios.TIOCGWINSZ, buf)

    # Return rows, cols
    return buf[0], buf[1]


class Vt100_Output(Output):
    """
    :param get_size: A callable which returns the `Size` of the output terminal.
    :param stdout: Any object with has a `write` and `flush` method + an 'encoding' property.
    :param true_color: Use 24bit color instead of 256 colors. (Can be a :class:`SimpleFilter`.)
        When `ansi_colors_only` is set, only 16 colors are used.
    :param ansi_colors_only: Restrict to 16 ANSI colors only.
    :param term: The terminal environment variable. (xterm, xterm-256color, linux, ...)
    :param write_binary: Encode the output before writing it. If `True` (the
        default), the `stdout` object is supposed to expose an `encoding` attribute.
    """
    def __init__(self, stdout, get_size, true_color=False,
                 ansi_colors_only=None, term=None, write_binary=True):
        assert callable(get_size)
        assert term is None or isinstance(term, six.text_type)
        assert all(hasattr(stdout, a) for a in ('write', 'flush'))

        if write_binary:
            assert hasattr(stdout, 'encoding')

        self._buffer = []
        self.stdout = stdout
        self.write_binary = write_binary
        self.get_size = get_size
        self.true_color = to_simple_filter(true_color)
        self.term = term or 'xterm'

        # ANSI colors only?
        if ansi_colors_only is None:
            # When not given, use the following default.
            ANSI_COLORS_ONLY = bool(os.environ.get(
                'PROMPT_TOOLKIT_ANSI_COLORS_ONLY', False))

            @Condition
            def ansi_colors_only():
                return ANSI_COLORS_ONLY or term in ('linux', 'eterm-color')
        else:
            ansi_colors_only = to_simple_filter(ansi_colors_only)

        self.ansi_colors_only = ansi_colors_only

        # Cache for escape codes.
        self._escape_code_cache = _EscapeCodeCache(ansi_colors_only=ansi_colors_only)
        self._escape_code_cache_true_color = _EscapeCodeCache(
            true_color=True, ansi_colors_only=ansi_colors_only)

    @classmethod
    def from_pty(cls, stdout, true_color=False, ansi_colors_only=None, term=None):
        """
        Create an Output class from a pseudo terminal.
        (This will take the dimensions by reading the pseudo
        terminal attributes.)
        """
        assert stdout.isatty()
        def get_size():
            rows, columns = _get_size(stdout.fileno())
            # If terminal (incorrectly) reports its size as 0, pick a reasonable default.
            # See https://github.com/ipython/ipython/issues/10071
            return Size(rows=(rows or 24), columns=(columns or 80))

        return cls(stdout, get_size, true_color=true_color,
                   ansi_colors_only=ansi_colors_only, term=term)

    def fileno(self):
        " Return file descriptor. "
        return self.stdout.fileno()

    def encoding(self):
        " Return encoding used for stdout. "
        return self.stdout.encoding

    def write_raw(self, data):
        """
        Write raw data to output.
        """
        self._buffer.append(data)

    def write(self, data):
        """
        Write text to output.
        (Removes vt100 escape codes. -- used for safely writing text.)
        """
        self._buffer.append(data.replace('\x1b', '?'))

    def set_title(self, title):
        """
        Set terminal title.
        """
        if self.term not in ('linux', 'eterm-color'):  # Not supported by the Linux console.
            self.write_raw('\x1b]2;%s\x07' % title.replace('\x1b', '').replace('\x07', ''))

    def clear_title(self):
        self.set_title('')

    def erase_screen(self):
        """
        Erases the screen with the background colour and moves the cursor to
        home.
        """
        self.write_raw('\x1b[2J')

    def enter_alternate_screen(self):
        self.write_raw('\x1b[?1049h\x1b[H')

    def quit_alternate_screen(self):
        self.write_raw('\x1b[?1049l')

    def enable_mouse_support(self):
        self.write_raw('\x1b[?1000h')

        # Enable urxvt Mouse mode. (For terminals that understand this.)
        self.write_raw('\x1b[?1015h')

        # Also enable Xterm SGR mouse mode. (For terminals that understand this.)
        self.write_raw('\x1b[?1006h')

        # Note: E.g. lxterminal understands 1000h, but not the urxvt or sgr
        #       extensions.

    def disable_mouse_support(self):
        self.write_raw('\x1b[?1000l')
        self.write_raw('\x1b[?1015l')
        self.write_raw('\x1b[?1006l')

    def erase_end_of_line(self):
        """
        Erases from the current cursor position to the end of the current line.
        """
        self.write_raw('\x1b[K')

    def erase_down(self):
        """
        Erases the screen from the current line down to the bottom of the
        screen.
        """
        self.write_raw('\x1b[J')

    def reset_attributes(self):
        self.write_raw('\x1b[0m')

    def set_attributes(self, attrs):
        """
        Create new style and output.

        :param attrs: `Attrs` instance.
        """
        if self.true_color() and not self.ansi_colors_only():
            self.write_raw(self._escape_code_cache_true_color[attrs])
        else:
            self.write_raw(self._escape_code_cache[attrs])

    def disable_autowrap(self):
        self.write_raw('\x1b[?7l')

    def enable_autowrap(self):
        self.write_raw('\x1b[?7h')

    def enable_bracketed_paste(self):
        self.write_raw('\x1b[?2004h')

    def disable_bracketed_paste(self):
        self.write_raw('\x1b[?2004l')

    def cursor_goto(self, row=0, column=0):
        """ Move cursor position. """
        self.write_raw('\x1b[%i;%iH' % (row, column))

    def cursor_up(self, amount):
        if amount == 0:
            pass
        elif amount == 1:
            self.write_raw('\x1b[A')
        else:
            self.write_raw('\x1b[%iA' % amount)

    def cursor_down(self, amount):
        if amount == 0:
            pass
        elif amount == 1:
            # Note: Not the same as '\n', '\n' can cause the window content to
            #       scroll.
            self.write_raw('\x1b[B')
        else:
            self.write_raw('\x1b[%iB' % amount)

    def cursor_forward(self, amount):
        if amount == 0:
            pass
        elif amount == 1:
            self.write_raw('\x1b[C')
        else:
            self.write_raw('\x1b[%iC' % amount)

    def cursor_backward(self, amount):
        if amount == 0:
            pass
        elif amount == 1:
            self.write_raw('\b')  # '\x1b[D'
        else:
            self.write_raw('\x1b[%iD' % amount)

    def hide_cursor(self):
        self.write_raw('\x1b[?25l')

    def show_cursor(self):
        self.write_raw('\x1b[?12l\x1b[?25h')  # Stop blinking cursor and show.

    def flush(self):
        """
        Write to output stream and flush.
        """
        if not self._buffer:
            return

        data = ''.join(self._buffer)

        try:
            # (We try to encode ourself, because that way we can replace
            # characters that don't exist in the character set, avoiding
            # UnicodeEncodeError crashes. E.g. u'\xb7' does not appear in 'ascii'.)
            # My Arch Linux installation of july 2015 reported 'ANSI_X3.4-1968'
            # for sys.stdout.encoding in xterm.
            if self.write_binary:
                if hasattr(self.stdout, 'buffer'):
                    out = self.stdout.buffer  # Py3.
                else:
                    out = self.stdout
                out.write(data.encode(self.stdout.encoding or 'utf-8', 'replace'))
            else:
                self.stdout.write(data)

            self.stdout.flush()
        except IOError as e:
            if e.args and e.args[0] == errno.EINTR:
                # Interrupted system call. Can happpen in case of a window
                # resize signal. (Just ignore. The resize handler will render
                # again anyway.)
                pass
            elif e.args and e.args[0] == 0:
                # This can happen when there is a lot of output and the user
                # sends a KeyboardInterrupt by pressing Control-C. E.g. in
                # a Python REPL when we execute "while True: print('test')".
                # (The `ptpython` REPL uses this `Output` class instead of
                # `stdout` directly -- in order to be network transparent.)
                # So, just ignore.
                pass
            else:
                raise

        self._buffer = []

    def ask_for_cpr(self):
        """
        Asks for a cursor position report (CPR).
        """
        self.write_raw('\x1b[6n')
        self.flush()

    def bell(self):
        " Sound bell. "
        self.write_raw('\a')
        self.flush()
