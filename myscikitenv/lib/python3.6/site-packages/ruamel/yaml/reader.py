# coding: utf-8

from __future__ import absolute_import

# This module contains abstractions for the input stream. You don't have to
# looks further, there are no pretty code.
#
# We define two classes here.
#
#   Mark(source, line, column)
# It's just a record and its only use is producing nice error messages.
# Parser does not use it for any other purposes.
#
#   Reader(source, data)
# Reader determines the encoding of `data` and converts it to unicode.
# Reader provides the following methods and attributes:
#   reader.peek(length=1) - return the next `length` characters
#   reader.forward(length=1) - move the current position to `length`
#      characters.
#   reader.index - the number of the current character.
#   reader.line, stream.column - the line and the column of the current
#      character.

import codecs

from ruamel.yaml.error import YAMLError, FileMark, StringMark, YAMLStreamError
from ruamel.yaml.compat import text_type, binary_type, PY3, UNICODE_SIZE
from ruamel.yaml.util import RegExp

if False:  # MYPY
    from typing import Any, Dict, Optional, List, Union, Text, Tuple, Optional  # NOQA
#    from ruamel.yaml.compat import StreamTextType  # NOQA

__all__ = ['Reader', 'ReaderError']


class ReaderError(YAMLError):
    def __init__(self, name, position, character, encoding, reason):
        # type: (Any, Any, Any, Any, Any) -> None
        self.name = name
        self.character = character
        self.position = position
        self.encoding = encoding
        self.reason = reason

    def __str__(self):
        # type: () -> str
        if isinstance(self.character, binary_type):
            return "'%s' codec can't decode byte #x%02x: %s\n" '  in "%s", position %d' % (
                self.encoding,
                ord(self.character),
                self.reason,
                self.name,
                self.position,
            )
        else:
            return 'unacceptable character #x%04x: %s\n' '  in "%s", position %d' % (
                self.character,
                self.reason,
                self.name,
                self.position,
            )


class Reader(object):
    # Reader:
    # - determines the data encoding and converts it to a unicode string,
    # - checks if characters are in allowed range,
    # - adds '\0' to the end.

    # Reader accepts
    #  - a `str` object (PY2) / a `bytes` object (PY3),
    #  - a `unicode` object (PY2) / a `str` object (PY3),
    #  - a file-like object with its `read` method returning `str`,
    #  - a file-like object with its `read` method returning `unicode`.

    # Yeah, it's ugly and slow.

    def __init__(self, stream, loader=None):
        # type: (Any, Any) -> None
        self.loader = loader
        if self.loader is not None and getattr(self.loader, '_reader', None) is None:
            self.loader._reader = self
        self.reset_reader()
        self.stream = stream  # type: Any  # as .read is called

    def reset_reader(self):
        # type: () -> None
        self.name = None  # type: Any
        self.stream_pointer = 0
        self.eof = True
        self.buffer = ""
        self.pointer = 0
        self.raw_buffer = None  # type: Any
        self.raw_decode = None
        self.encoding = None  # type: Optional[Text]
        self.index = 0
        self.line = 0
        self.column = 0

    @property
    def stream(self):
        # type: () -> Any
        try:
            return self._stream
        except AttributeError:
            raise YAMLStreamError('input stream needs to specified')

    @stream.setter
    def stream(self, val):
        # type: (Any) -> None
        if val is None:
            return
        self._stream = None
        if isinstance(val, text_type):
            self.name = '<unicode string>'
            self.check_printable(val)
            self.buffer = val + u'\0'  # type: ignore
        elif isinstance(val, binary_type):
            self.name = '<byte string>'
            self.raw_buffer = val
            self.determine_encoding()
        else:
            if not hasattr(val, 'read'):
                raise YAMLStreamError('stream argument needs to have a read() method')
            self._stream = val
            self.name = getattr(self.stream, 'name', '<file>')
            self.eof = False
            self.raw_buffer = None
            self.determine_encoding()

    def peek(self, index=0):
        # type: (int) -> Text
        try:
            return self.buffer[self.pointer + index]
        except IndexError:
            self.update(index + 1)
            return self.buffer[self.pointer + index]

    def prefix(self, length=1):
        # type: (int) -> Any
        if self.pointer + length >= len(self.buffer):
            self.update(length)
        return self.buffer[self.pointer : self.pointer + length]

    def forward_1_1(self, length=1):
        # type: (int) -> None
        if self.pointer + length + 1 >= len(self.buffer):
            self.update(length + 1)
        while length != 0:
            ch = self.buffer[self.pointer]
            self.pointer += 1
            self.index += 1
            if ch in u'\n\x85\u2028\u2029' or (
                ch == u'\r' and self.buffer[self.pointer] != u'\n'
            ):
                self.line += 1
                self.column = 0
            elif ch != u'\uFEFF':
                self.column += 1
            length -= 1

    def forward(self, length=1):
        # type: (int) -> None
        if self.pointer + length + 1 >= len(self.buffer):
            self.update(length + 1)
        while length != 0:
            ch = self.buffer[self.pointer]
            self.pointer += 1
            self.index += 1
            if ch == u'\n' or (ch == u'\r' and self.buffer[self.pointer] != u'\n'):
                self.line += 1
                self.column = 0
            elif ch != u'\uFEFF':
                self.column += 1
            length -= 1

    def get_mark(self):
        # type: () -> Any
        if self.stream is None:
            return StringMark(
                self.name, self.index, self.line, self.column, self.buffer, self.pointer
            )
        else:
            return FileMark(self.name, self.index, self.line, self.column)

    def determine_encoding(self):
        # type: () -> None
        while not self.eof and (self.raw_buffer is None or len(self.raw_buffer) < 2):
            self.update_raw()
        if isinstance(self.raw_buffer, binary_type):
            if self.raw_buffer.startswith(codecs.BOM_UTF16_LE):
                self.raw_decode = codecs.utf_16_le_decode  # type: ignore
                self.encoding = 'utf-16-le'
            elif self.raw_buffer.startswith(codecs.BOM_UTF16_BE):
                self.raw_decode = codecs.utf_16_be_decode  # type: ignore
                self.encoding = 'utf-16-be'
            else:
                self.raw_decode = codecs.utf_8_decode  # type: ignore
                self.encoding = 'utf-8'
        self.update(1)

    if UNICODE_SIZE == 2:
        NON_PRINTABLE = RegExp(
            u'[^\x09\x0A\x0D\x20-\x7E\x85' u'\xA0-\uD7FF' u'\uE000-\uFFFD' u']'
        )
    else:
        NON_PRINTABLE = RegExp(
            u'[^\x09\x0A\x0D\x20-\x7E\x85'
            u'\xA0-\uD7FF'
            u'\uE000-\uFFFD'
            u'\U00010000-\U0010FFFF'
            u']'
        )

    _printable_ascii = ('\x09\x0A\x0D' + "".join(map(chr, range(0x20, 0x7F)))).encode('ascii')

    @classmethod
    def _get_non_printable_ascii(cls, data):  # type: ignore
        # type: (Text, bytes) -> Optional[Tuple[int, Text]]
        ascii_bytes = data.encode('ascii')
        non_printables = ascii_bytes.translate(None, cls._printable_ascii)  # type: ignore
        if not non_printables:
            return None
        non_printable = non_printables[:1]
        return ascii_bytes.index(non_printable), non_printable.decode('ascii')

    @classmethod
    def _get_non_printable_regex(cls, data):
        # type: (Text) -> Optional[Tuple[int, Text]]
        match = cls.NON_PRINTABLE.search(data)
        if not bool(match):
            return None
        return match.start(), match.group()

    @classmethod
    def _get_non_printable(cls, data):
        # type: (Text) -> Optional[Tuple[int, Text]]
        try:
            return cls._get_non_printable_ascii(data)  # type: ignore
        except UnicodeEncodeError:
            return cls._get_non_printable_regex(data)

    def check_printable(self, data):
        # type: (Any) -> None
        non_printable_match = self._get_non_printable(data)
        if non_printable_match is not None:
            start, character = non_printable_match
            position = self.index + (len(self.buffer) - self.pointer) + start
            raise ReaderError(
                self.name,
                position,
                ord(character),
                'unicode',
                'special characters are not allowed',
            )

    def update(self, length):
        # type: (int) -> None
        if self.raw_buffer is None:
            return
        self.buffer = self.buffer[self.pointer :]
        self.pointer = 0
        while len(self.buffer) < length:
            if not self.eof:
                self.update_raw()
            if self.raw_decode is not None:
                try:
                    data, converted = self.raw_decode(self.raw_buffer, 'strict', self.eof)
                except UnicodeDecodeError as exc:
                    if PY3:
                        character = self.raw_buffer[exc.start]
                    else:
                        character = exc.object[exc.start]
                    if self.stream is not None:
                        position = self.stream_pointer - len(self.raw_buffer) + exc.start
                    elif self.stream is not None:
                        position = self.stream_pointer - len(self.raw_buffer) + exc.start
                    else:
                        position = exc.start
                    raise ReaderError(self.name, position, character, exc.encoding, exc.reason)
            else:
                data = self.raw_buffer
                converted = len(data)
            self.check_printable(data)
            self.buffer += data
            self.raw_buffer = self.raw_buffer[converted:]
            if self.eof:
                self.buffer += '\0'
                self.raw_buffer = None
                break

    def update_raw(self, size=None):
        # type: (Optional[int]) -> None
        if size is None:
            size = 4096 if PY3 else 1024
        data = self.stream.read(size)
        if self.raw_buffer is None:
            self.raw_buffer = data
        else:
            self.raw_buffer += data
        self.stream_pointer += len(data)
        if not data:
            self.eof = True


# try:
#     import psyco
#     psyco.bind(Reader)
# except ImportError:
#     pass
