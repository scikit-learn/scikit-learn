from __future__ import unicode_literals

from prompt_toolkit.cache import FastDictCache
from prompt_toolkit.token import Token
from prompt_toolkit.utils import get_cwidth

from collections import defaultdict, namedtuple

__all__ = (
    'Point',
    'Size',
    'Screen',
    'Char',
)


Point = namedtuple('Point', 'y x')
Size = namedtuple('Size', 'rows columns')


class Char(object):
    """
    Represent a single character in a :class:`.Screen`.

    This should be considered immutable.
    """
    __slots__ = ('char', 'token', 'width')

    # If we end up having one of these special control sequences in the input string,
    # we should display them as follows:
    # Usually this happens after a "quoted insert".
    display_mappings = {
        '\x00': '^@',  # Control space
        '\x01': '^A',
        '\x02': '^B',
        '\x03': '^C',
        '\x04': '^D',
        '\x05': '^E',
        '\x06': '^F',
        '\x07': '^G',
        '\x08': '^H',
        '\x09': '^I',
        '\x0a': '^J',
        '\x0b': '^K',
        '\x0c': '^L',
        '\x0d': '^M',
        '\x0e': '^N',
        '\x0f': '^O',
        '\x10': '^P',
        '\x11': '^Q',
        '\x12': '^R',
        '\x13': '^S',
        '\x14': '^T',
        '\x15': '^U',
        '\x16': '^V',
        '\x17': '^W',
        '\x18': '^X',
        '\x19': '^Y',
        '\x1a': '^Z',
        '\x1b': '^[',  # Escape
        '\x1c': '^\\',
        '\x1d': '^]',
        '\x1f': '^_',
        '\x7f': '^?',  # Backspace
    }

    def __init__(self, char=' ', token=Token):
        # If this character has to be displayed otherwise, take that one.
        char = self.display_mappings.get(char, char)

        self.char = char
        self.token = token

        # Calculate width. (We always need this, so better to store it directly
        # as a member for performance.)
        self.width = get_cwidth(char)

    def __eq__(self, other):
        return self.char == other.char and self.token == other.token

    def __ne__(self, other):
        # Not equal: We don't do `not char.__eq__` here, because of the
        # performance of calling yet another function.
        return self.char != other.char or self.token != other.token

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.char, self.token)


_CHAR_CACHE = FastDictCache(Char, size=1000 * 1000)
Transparent = Token.Transparent


class Screen(object):
    """
    Two dimentional buffer of :class:`.Char` instances.
    """
    def __init__(self, default_char=None, initial_width=0, initial_height=0):
        if default_char is None:
            default_char = _CHAR_CACHE[' ', Transparent]

        self.data_buffer = defaultdict(lambda: defaultdict(lambda: default_char))

        #: Escape sequences to be injected.
        self.zero_width_escapes = defaultdict(lambda: defaultdict(lambda: ''))

        #: Position of the cursor.
        self.cursor_position = Point(y=0, x=0)

        #: Visibility of the cursor.
        self.show_cursor = True

        #: (Optional) Where to position the menu. E.g. at the start of a completion.
        #: (We can't use the cursor position, because we don't want the
        #: completion menu to change its position when we browse through all the
        #: completions.)
        self.menu_position = None

        #: Currently used width/height of the screen. This will increase when
        #: data is written to the screen.
        self.width = initial_width or 0
        self.height = initial_height or 0

    def replace_all_tokens(self, token):
        """
        For all the characters in the screen. Set the token to the given `token`.
        """
        b = self.data_buffer

        for y, row in b.items():
            for x, char in row.items():
                b[y][x] = _CHAR_CACHE[char.char, token]


class WritePosition(object):
    def __init__(self, xpos, ypos, width, height, extended_height=None):
        assert height >= 0
        assert extended_height is None or extended_height >= 0
        assert width >= 0
        # xpos and ypos can be negative. (A float can be partially visible.)

        self.xpos = xpos
        self.ypos = ypos
        self.width = width
        self.height = height
        self.extended_height = extended_height or height

    def __repr__(self):
        return '%s(%r, %r, %r, %r, %r)' % (
            self.__class__.__name__,
            self.xpos, self.ypos, self.width, self.height, self.extended_height)
