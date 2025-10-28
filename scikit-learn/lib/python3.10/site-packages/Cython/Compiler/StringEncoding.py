#
#   Cython -- encoding related tools
#


import re
import sys

join_bytes = b''.join


class UnicodeLiteralBuilder:
    """Assemble a unicode string.
    """
    def __init__(self):
        self.chars = []

    def append(self, characters):
        assert isinstance(characters, str), f"Expected str, got {type(characters)}"
        self.chars.append(characters)

    if sys.maxunicode == 65535:
        def append_charval(self, char_number):
            if char_number > 65535:
                # wide Unicode character on narrow platform => replace
                # by surrogate pair
                char_number -= 0x10000
                self.chars.append( chr((char_number // 1024) + 0xD800) )
                self.chars.append( chr((char_number  % 1024) + 0xDC00) )
            else:
                self.chars.append( chr(char_number) )
    else:
        def append_charval(self, char_number):
            self.chars.append( chr(char_number) )

    def append_uescape(self, char_number, escape_string):
        self.append_charval(char_number)

    def getstring(self):
        return EncodedString(''.join(self.chars))

    def getstrings(self):
        return (None, self.getstring())


class BytesLiteralBuilder:
    """Assemble a byte string or char value.
    """
    def __init__(self, target_encoding):
        self.chars = []
        self.target_encoding = target_encoding

    def append(self, characters):
        if isinstance(characters, str):
            characters = characters.encode(self.target_encoding)
        assert isinstance(characters, bytes), str(type(characters))
        self.chars.append(characters)

    def append_charval(self, char_number):
        self.chars.append( chr(char_number).encode('ISO-8859-1') )

    def append_uescape(self, char_number, escape_string):
        self.append(escape_string)

    def getstring(self):
        # this *must* return a byte string!
        return bytes_literal(join_bytes(self.chars), self.target_encoding)

    def getchar(self):
        # this *must* return a byte string!
        return self.getstring()

    def getstrings(self):
        return (self.getstring(), None)


class StrLiteralBuilder:
    """Assemble both a bytes and a unicode representation of a string.
    """
    def __init__(self, target_encoding):
        self._bytes   = BytesLiteralBuilder(target_encoding)
        self._unicode = UnicodeLiteralBuilder()

    def append(self, characters):
        self._bytes.append(characters)
        self._unicode.append(characters)

    def append_charval(self, char_number):
        self._bytes.append_charval(char_number)
        self._unicode.append_charval(char_number)

    def append_uescape(self, char_number, escape_string):
        self._bytes.append(escape_string)
        self._unicode.append_charval(char_number)

    def getstrings(self):
        return (self._bytes.getstring(), self._unicode.getstring())


class EncodedString(str):
    # unicode string subclass to keep track of the original encoding.
    # 'encoding' is None for unicode strings and the source encoding
    # otherwise
    encoding = None

    def __deepcopy__(self, memo):
        return self

    def byteencode(self):
        assert self.encoding is not None
        return self.encode(self.encoding)

    def utf8encode(self):
        assert self.encoding is None
        return self.encode("UTF-8")

    @property
    def is_unicode(self):
        return self.encoding is None

    def contains_surrogates(self):
        return string_contains_surrogates(self)

    def as_utf8_string(self):
        return bytes_literal(self.utf8encode(), 'utf8')

    def as_c_string_literal(self):
        # first encodes the string then produces a c string literal
        if self.encoding is None:
            s = self.as_utf8_string()
        else:
            s = bytes_literal(self.byteencode(), self.encoding)
        return s.as_c_string_literal()


def string_contains_surrogates(ustring):
    """
    Check if the unicode string contains surrogate code points
    on a CPython platform with wide (UCS-4) or narrow (UTF-16)
    Unicode, i.e. characters that would be spelled as two
    separate code units on a narrow platform.
    """
    for c in map(ord, ustring):
        if c > 65535:  # can only happen on wide platforms
            return True
        if 0xD800 <= c <= 0xDFFF:
            return True
    return False


def string_contains_lone_surrogates(ustring):
    """
    Check if the unicode string contains lone surrogate code points
    on a CPython platform with wide (UCS-4) or narrow (UTF-16)
    Unicode, i.e. characters that would be spelled as two
    separate code units on a narrow platform, but that do not form a pair.
    """
    last_was_start = False
    unicode_uses_surrogate_encoding = sys.maxunicode == 65535
    for c in map(ord, ustring):
        # surrogates tend to be rare
        if c < 0xD800 or c > 0xDFFF:
            if last_was_start:
                return True
        elif not unicode_uses_surrogate_encoding:
            # on 32bit Unicode platforms, there is never a pair
            return True
        elif c <= 0xDBFF:
            if last_was_start:
                return True  # lone start
            last_was_start = True
        else:
            if not last_was_start:
                return True  # lone end
            last_was_start = False
    return last_was_start


class BytesLiteral(bytes):
    # bytes subclass that is compatible with EncodedString
    encoding = None

    def __deepcopy__(self, memo):
        return self

    def byteencode(self):
        return bytes(self)

    def utf8encode(self):
        assert False, "this is not a unicode string: %r" % self

    def __str__(self):
        """Fake-decode the byte string to unicode to support %
        formatting of unicode strings.
        """
        return self.decode('ISO-8859-1')

    is_unicode = False

    def as_c_string_literal(self):
        value = split_string_literal(escape_byte_string(self))
        return '"%s"' % value


def bytes_literal(s, encoding):
    assert isinstance(s, bytes)
    s = BytesLiteral(s)
    s.encoding = encoding
    return s


def encoded_string(s, encoding):
    assert isinstance(s, (str, bytes))
    s = EncodedString(s)
    if encoding is not None:
        s.encoding = encoding
    return s

def encoded_string_or_bytes_literal(s, encoding):
    if isinstance(s, bytes):
        return bytes_literal(s, encoding)
    else:
        return encoded_string(s, encoding)


char_from_escape_sequence = {
    r'\a' : '\a',
    r'\b' : '\b',
    r'\f' : '\f',
    r'\n' : '\n',
    r'\r' : '\r',
    r'\t' : '\t',
    r'\v' : '\v',
    }.get

_c_special = ('\\', '??', '"') + tuple(map(chr, range(32)))


def _to_escape_sequence(s):
    if s in '\n\r\t':
        return repr(s)[1:-1]
    elif s == '"':
        return r'\"'
    elif s == '\\':
        return r'\\'
    else:
        # within a character sequence, oct passes much better than hex
        return ''.join([f'\\{ord(c):03o}' for c in s])


def _build_specials_replacer():
    subexps = []
    replacements = {}
    for special in _c_special:
        regexp = ''.join(['[%s]' % c.replace('\\', '\\\\') for c in special])
        subexps.append(regexp)
        replacements[special.encode('ASCII')] = _to_escape_sequence(special).encode('ASCII')
    sub = re.compile(('(%s)' % '|'.join(subexps)).encode('ASCII')).sub
    def replace_specials(m):
        return replacements[m.group(1)]
    def replace(s):
        return sub(replace_specials, s)
    return replace

_replace_specials = _build_specials_replacer()


def escape_char(c):
    c = c.decode('ISO-8859-1')
    if c in '\n\r\t\\':
        return repr(c)[1:-1]
    elif c == "'":
        return "\\'"
    n = ord(c)
    if n < 32 or n >= 127:
        # hex works well for characters
        return "\\x%02X" % n
    else:
        # strictly £, @ and ` (which fall in this list) are only allowed
        # in C23. But practically they're well-supported earlier.
        return c

def escape_byte_string(s):
    """Escape a byte string so that it can be written into C code.
    Note that this returns a Unicode string instead which, when
    encoded as ASCII, will result in the correct byte sequence
    being written.
    """
    s = _replace_specials(s)
    try:
        return s.decode("ASCII")  #  trial decoding: plain ASCII => done
    except UnicodeDecodeError:
        pass
    s_new = bytearray()
    append, extend = s_new.append, s_new.extend
    for b in s:
        if b >= 127:
            extend(b'\\%03o' % b)
        else:
            append(b)
    return s_new.decode('ASCII')

def split_string_literal(s, limit=2000):
    # MSVC can't handle long string literals.
    if len(s) < limit:
        return s
    else:
        start = 0
        chunks = []
        while start < len(s):
            end = start + limit
            if len(s) > end-4 and '\\' in s[end-4:end]:
                end -= 4 - s[end-4:end].find('\\')  # just before the backslash
                while s[end-1] == '\\':
                    end -= 1
                    if end == start:
                        # must have been a long line of backslashes
                        end = start + limit - (limit % 2) - 4
                        break
            chunks.append(s[start:end])
            start = end
        return '""'.join(chunks)

def encode_pyunicode_string(s):
    """Create Py_UNICODE[] representation of a given unicode string.
    """
    s = list(map(ord, s)) + [0]

    if sys.maxunicode >= 0x10000:  # Wide build or Py3.3
        utf16, utf32 = [], s
        for code_point in s:
            if code_point >= 0x10000:  # outside of BMP
                high, low = divmod(code_point - 0x10000, 1024)
                utf16.append(high + 0xD800)
                utf16.append(low + 0xDC00)
            else:
                utf16.append(code_point)
    else:
        utf16, utf32 = s, []
        for code_unit in s:
            if 0xDC00 <= code_unit <= 0xDFFF and utf32 and 0xD800 <= utf32[-1] <= 0xDBFF:
                high, low = utf32[-1], code_unit
                utf32[-1] = ((high & 0x3FF) << 10) + (low & 0x3FF) + 0x10000
            else:
                utf32.append(code_unit)

    if utf16 == utf32:
        utf16 = []
    return ",".join(map(str, utf16)), ",".join(map(str, utf32))
