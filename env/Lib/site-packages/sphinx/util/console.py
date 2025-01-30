"""Format colored console output."""

from __future__ import annotations

import os
import re
import shutil
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Final

    # fmt: off
    def reset(text: str) -> str: ...  # NoQA: E704
    def bold(text: str) -> str: ...  # NoQA: E704
    def faint(text: str) -> str: ...  # NoQA: E704
    def standout(text: str) -> str: ...  # NoQA: E704
    def underline(text: str) -> str: ...  # NoQA: E704
    def blink(text: str) -> str: ...  # NoQA: E704

    def black(text: str) -> str: ...  # NoQA: E704
    def white(text: str) -> str: ...  # NoQA: E704
    def red(text: str) -> str: ...  # NoQA: E704
    def green(text: str) -> str: ...  # NoQA: E704
    def yellow(text: str) -> str: ...  # NoQA: E704
    def blue(text: str) -> str: ...  # NoQA: E704
    def fuchsia(text: str) -> str: ...  # NoQA: E704
    def teal(text: str) -> str: ...  # NoQA: E704

    def darkgray(text: str) -> str: ...  # NoQA: E704
    def lightgray(text: str) -> str: ...  # NoQA: E704
    def darkred(text: str) -> str: ...  # NoQA: E704
    def darkgreen(text: str) -> str: ...  # NoQA: E704
    def brown(text: str) -> str: ...  # NoQA: E704
    def darkblue(text: str) -> str: ...  # NoQA: E704
    def purple(text: str) -> str: ...  # NoQA: E704
    def turquoise(text: str) -> str: ...  # NoQA: E704
    # fmt: on

try:
    # check if colorama is installed to support color on Windows
    import colorama

    COLORAMA_AVAILABLE = True
except ImportError:
    COLORAMA_AVAILABLE = False

_CSI: Final[str] = re.escape('\x1b[')  # 'ESC [': Control Sequence Introducer

# Pattern matching ANSI control sequences containing colors.
_ansi_color_re: Final[re.Pattern[str]] = re.compile(r'\x1b\[(?:\d+;){0,2}\d*m')

_ansi_re: Final[re.Pattern[str]] = re.compile(
    _CSI
    + r"""
    (?:
      (?:\d+;){0,2}\d*m     # ANSI color code    ('m' is equivalent to '0m')
    |
      [012]?K               # ANSI Erase in Line ('K' is equivalent to '0K')
    )""",
    re.VERBOSE | re.ASCII,
)
"""Pattern matching ANSI CSI colors (SGR) and erase line (EL) sequences.

See :func:`strip_escape_sequences` for details.
"""

codes: dict[str, str] = {}


def terminal_safe(s: str) -> str:
    """Safely encode a string for printing to the terminal."""
    return s.encode('ascii', 'backslashreplace').decode('ascii')


def get_terminal_width() -> int:
    """Return the width of the terminal in columns."""
    return shutil.get_terminal_size().columns - 1


_tw: int = get_terminal_width()


def term_width_line(text: str) -> str:
    if not codes:
        # if no coloring, don't output fancy backspaces
        return text + '\n'
    else:
        # codes are not displayed, this must be taken into account
        return text.ljust(_tw + len(text) - len(strip_escape_sequences(text))) + '\r'


def color_terminal() -> bool:
    if 'NO_COLOR' in os.environ:
        return False
    if sys.platform == 'win32' and COLORAMA_AVAILABLE:
        colorama.just_fix_windows_console()
        return True
    if 'FORCE_COLOR' in os.environ:
        return True
    if not hasattr(sys.stdout, 'isatty'):
        return False
    if not sys.stdout.isatty():
        return False
    if 'COLORTERM' in os.environ:
        return True
    term = os.environ.get('TERM', 'dumb').lower()
    return term in ('xterm', 'linux') or 'color' in term


def nocolor() -> None:
    if sys.platform == 'win32' and COLORAMA_AVAILABLE:
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
    """Remove the ANSI color codes in a string *s*.

    .. caution::

       This function is not meant to be used in production and should only
       be used for testing Sphinx's output messages.

    .. seealso:: :func:`strip_escape_sequences`
    """
    return _ansi_color_re.sub('', s)


def strip_escape_sequences(text: str, /) -> str:
    r"""Remove the ANSI CSI colors and "erase in line" sequences.

    Other `escape sequences `__ (e.g., VT100-specific functions) are not
    supported and only control sequences *natively* known to Sphinx (i.e.,
    colors declared in this module and "erase entire line" (``'\x1b[2K'``))
    are eliminated by this function.

    .. caution::

       This function is not meant to be used in production and should only
       be used for testing Sphinx's output messages that were not tempered
       with by third-party extensions.

    .. versionadded:: 7.3

       This function is added as an *experimental* feature.

    __ https://en.wikipedia.org/wiki/ANSI_escape_code
    """
    return _ansi_re.sub('', text)


def create_color_func(name: str) -> None:
    def inner(text: str) -> str:
        return colorize(name, text)

    globals()[name] = inner


_attrs = {
    'reset': '39;49;00m',
    'bold': '01m',
    'faint': '02m',
    'standout': '03m',
    'underline': '04m',
    'blink': '05m',
}

for __name, __value in _attrs.items():
    codes[__name] = '\x1b[' + __value

_colors = [
    ('black', 'darkgray'),
    ('darkred', 'red'),
    ('darkgreen', 'green'),
    ('brown', 'yellow'),
    ('darkblue', 'blue'),
    ('purple', 'fuchsia'),
    ('turquoise', 'teal'),
    ('lightgray', 'white'),
]

for __i, (__dark, __light) in enumerate(_colors, 30):
    codes[__dark] = '\x1b[%im' % __i
    codes[__light] = '\x1b[%im' % (__i + 60)

_orig_codes = codes.copy()

for _name in codes:
    create_color_func(_name)
