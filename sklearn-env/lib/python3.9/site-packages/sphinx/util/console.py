"""
    sphinx.util.console
    ~~~~~~~~~~~~~~~~~~~

    Format colored console output.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import os
import re
import sys
from typing import Dict, Pattern

try:
    # check if colorama is installed to support color on Windows
    import colorama
except ImportError:
    colorama = None


_ansi_re: Pattern = re.compile('\x1b\\[(\\d\\d;){0,2}\\d\\dm')
codes: Dict[str, str] = {}


def terminal_safe(s: str) -> str:
    """Safely encode a string for printing to the terminal."""
    return s.encode('ascii', 'backslashreplace').decode('ascii')


def get_terminal_width() -> int:
    """Borrowed from the py lib."""
    try:
        import fcntl
        import struct
        import termios
        call = fcntl.ioctl(0, termios.TIOCGWINSZ, struct.pack('hhhh', 0, 0, 0, 0))
        height, width = struct.unpack('hhhh', call)[:2]
        terminal_width = width
    except Exception:
        # FALLBACK
        terminal_width = int(os.environ.get('COLUMNS', "80")) - 1
    return terminal_width


_tw: int = get_terminal_width()


def term_width_line(text: str) -> str:
    if not codes:
        # if no coloring, don't output fancy backspaces
        return text + '\n'
    else:
        # codes are not displayed, this must be taken into account
        return text.ljust(_tw + len(text) - len(_ansi_re.sub('', text))) + '\r'


def color_terminal() -> bool:
    if sys.platform == 'win32' and colorama is not None:
        colorama.init()
        return True
    if not hasattr(sys.stdout, 'isatty'):
        return False
    if not sys.stdout.isatty():
        return False
    if 'COLORTERM' in os.environ:
        return True
    term = os.environ.get('TERM', 'dumb').lower()
    if term in ('xterm', 'linux') or 'color' in term:
        return True
    return False


def nocolor() -> None:
    if sys.platform == 'win32' and colorama is not None:
        colorama.deinit()
    codes.clear()


def coloron() -> None:
    codes.update(_orig_codes)


def colorize(name: str, text: str, input_mode: bool = False) -> str:
    def escseq(name: str) -> str:
        # Wrap escape sequence with ``\1`` and ``\2`` to let readline know
        # it is non-printable characters
        # ref: https://tiswww.case.edu/php/chet/readline/readline.html
        #
        # Note: This hack does not work well in Windows (see #5059)
        escape = codes.get(name, '')
        if input_mode and escape and sys.platform != 'win32':
            return '\1' + escape + '\2'
        else:
            return escape

    return escseq(name) + text + escseq('reset')


def strip_colors(s: str) -> str:
    return re.compile('\x1b.*?m').sub('', s)


def create_color_func(name: str) -> None:
    def inner(text: str) -> str:
        return colorize(name, text)
    globals()[name] = inner


_attrs = {
    'reset':     '39;49;00m',
    'bold':      '01m',
    'faint':     '02m',
    'standout':  '03m',
    'underline': '04m',
    'blink':     '05m',
}

for _name, _value in _attrs.items():
    codes[_name] = '\x1b[' + _value

_colors = [
    ('black',     'darkgray'),
    ('darkred',   'red'),
    ('darkgreen', 'green'),
    ('brown',     'yellow'),
    ('darkblue',  'blue'),
    ('purple',    'fuchsia'),
    ('turquoise', 'teal'),
    ('lightgray', 'white'),
]

for i, (dark, light) in enumerate(_colors, 30):
    codes[dark] = '\x1b[%im' % i
    codes[light] = '\x1b[%im' % (i + 60)

_orig_codes = codes.copy()

for _name in codes:
    create_color_func(_name)
