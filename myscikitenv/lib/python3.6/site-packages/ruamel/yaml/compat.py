# coding: utf-8

from __future__ import print_function

# partially from package six by Benjamin Peterson

import sys
import os
import types
import traceback
from abc import abstractmethod


# fmt: off
if False:  # MYPY
    from typing import Any, Dict, Optional, List, Union, BinaryIO, IO, Text, Tuple, Optional  # NOQA
# fmt: on

_DEFAULT_YAML_VERSION = (1, 2)

try:
    from ruamel.ordereddict import ordereddict
except:  # NOQA
    try:
        from collections import OrderedDict
    except ImportError:
        from ordereddict import OrderedDict  # type: ignore
    # to get the right name import ... as ordereddict doesn't do that

    class ordereddict(OrderedDict):  # type: ignore
        if not hasattr(OrderedDict, 'insert'):

            def insert(self, pos, key, value):
                # type: (int, Any, Any) -> None
                if pos >= len(self):
                    self[key] = value
                    return
                od = ordereddict()
                od.update(self)
                for k in od:
                    del self[k]
                for index, old_key in enumerate(od):
                    if pos == index:
                        self[key] = value
                    self[old_key] = od[old_key]


PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3


if PY3:

    def utf8(s):
        # type: (str) -> str
        return s

    def to_str(s):
        # type: (str) -> str
        return s

    def to_unicode(s):
        # type: (str) -> str
        return s


else:
    if False:
        unicode = str

    def utf8(s):
        # type: (unicode) -> str
        return s.encode('utf-8')

    def to_str(s):
        # type: (str) -> str
        return str(s)

    def to_unicode(s):
        # type: (str) -> unicode
        return unicode(s)  # NOQA


if PY3:
    string_types = str
    integer_types = int
    class_types = type
    text_type = str
    binary_type = bytes

    MAXSIZE = sys.maxsize
    unichr = chr
    import io

    StringIO = io.StringIO
    BytesIO = io.BytesIO
    # have unlimited precision
    no_limit_int = int
    from collections.abc import Hashable, MutableSequence, MutableMapping, Mapping  # NOQA

else:
    string_types = basestring  # NOQA
    integer_types = (int, long)  # NOQA
    class_types = (type, types.ClassType)
    text_type = unicode  # NOQA
    binary_type = str

    # to allow importing
    unichr = unichr  # type: ignore
    from StringIO import StringIO as _StringIO

    StringIO = _StringIO
    import cStringIO

    BytesIO = cStringIO.StringIO
    # have unlimited precision
    no_limit_int = long  # NOQA not available on Python 3
    from collections import Hashable, MutableSequence, MutableMapping, Mapping  # NOQA

if False:  # MYPY
    # StreamType = Union[BinaryIO, IO[str], IO[unicode],  StringIO]
    # StreamType = Union[BinaryIO, IO[str], StringIO]  # type: ignore
    StreamType = Any

    StreamTextType = Union[Text, StreamType]
    VersionType = Union[List[int], str, Tuple[int, int]]

if PY3:
    builtins_module = 'builtins'
else:
    builtins_module = '__builtin__'

UNICODE_SIZE = 4 if sys.maxunicode > 65535 else 2


def with_metaclass(meta, *bases):
    # type: (Any, Any) -> Any
    """Create a base class with a metaclass."""
    return meta('NewBase', bases, {})


DBG_TOKEN = 1
DBG_EVENT = 2
DBG_NODE = 4


_debug = None  # type: Optional[int]
if 'RUAMELDEBUG' in os.environ:
    _debugx = os.environ.get('RUAMELDEBUG')
    if _debugx is None:
        _debug = 0
    else:
        _debug = int(_debugx)


if bool(_debug):

    class ObjectCounter(object):
        def __init__(self):
            # type: () -> None
            self.map = {}  # type: Dict[Any, Any]

        def __call__(self, k):
            # type: (Any) -> None
            self.map[k] = self.map.get(k, 0) + 1

        def dump(self):
            # type: () -> None
            for k in sorted(self.map):
                sys.stdout.write('{} -> {}'.format(k, self.map[k]))

    object_counter = ObjectCounter()


# used from yaml util when testing
def dbg(val=None):
    # type: (Any) -> Any
    global _debug
    if _debug is None:
        # set to true or false
        _debugx = os.environ.get('YAMLDEBUG')
        if _debugx is None:
            _debug = 0
        else:
            _debug = int(_debugx)
    if val is None:
        return _debug
    return _debug & val


class Nprint(object):
    def __init__(self, file_name=None):
        # type: (Any) -> None
        self._max_print = None  # type: Any
        self._count = None  # type: Any
        self._file_name = file_name

    def __call__(self, *args, **kw):
        # type: (Any, Any) -> None
        if not bool(_debug):
            return
        out = sys.stdout if self._file_name is None else open(self._file_name, 'a')
        dbgprint = print  # to fool checking for print statements by dv utility
        kw1 = kw.copy()
        kw1['file'] = out
        dbgprint(*args, **kw1)
        out.flush()
        if self._max_print is not None:
            if self._count is None:
                self._count = self._max_print
            self._count -= 1
            if self._count == 0:
                dbgprint('forced exit\n')
                traceback.print_stack()
                out.flush()
                sys.exit(0)
        if self._file_name:
            out.close()

    def set_max_print(self, i):
        # type: (int) -> None
        self._max_print = i
        self._count = None


nprint = Nprint()
nprintf = Nprint('/var/tmp/ruamel.yaml.log')

# char checkers following production rules


def check_namespace_char(ch):
    # type: (Any) -> bool
    if u'\x21' <= ch <= u'\x7E':  # ! to ~
        return True
    if u'\xA0' <= ch <= u'\uD7FF':
        return True
    if (u'\uE000' <= ch <= u'\uFFFD') and ch != u'\uFEFF':  # excl. byte order mark
        return True
    if u'\U00010000' <= ch <= u'\U0010FFFF':
        return True
    return False


def check_anchorname_char(ch):
    # type: (Any) -> bool
    if ch in u',[]{}':
        return False
    return check_namespace_char(ch)


def version_tnf(t1, t2=None):
    # type: (Any, Any) -> Any
    from ruamel.yaml import version_info  # NOQA

    if version_info < t1:
        return True
    if t2 is not None and version_info < t2:
        return None
    return False


class MutableSliceableSequence(MutableSequence):  # type: ignore
    __slots__ = ()

    def __getitem__(self, index):
        # type: (Any) -> Any
        if not isinstance(index, slice):
            return self.__getsingleitem__(index)
        return type(self)([self[i] for i in range(*index.indices(len(self)))])  # type: ignore

    def __setitem__(self, index, value):
        # type: (Any, Any) -> None
        if not isinstance(index, slice):
            return self.__setsingleitem__(index, value)
        assert iter(value)
        # nprint(index.start, index.stop, index.step, index.indices(len(self)))
        if index.step is None:
            del self[index.start : index.stop]
            for elem in reversed(value):
                self.insert(0 if index.start is None else index.start, elem)
        else:
            range_parms = index.indices(len(self))
            nr_assigned_items = (range_parms[1] - range_parms[0] - 1) // range_parms[2] + 1
            # need to test before changing, in case TypeError is caught
            if nr_assigned_items < len(value):
                raise TypeError(
                    'too many elements in value {} < {}'.format(nr_assigned_items, len(value))
                )
            elif nr_assigned_items > len(value):
                raise TypeError(
                    'not enough elements in value {} > {}'.format(
                        nr_assigned_items, len(value)
                    )
                )
            for idx, i in enumerate(range(*range_parms)):
                self[i] = value[idx]

    def __delitem__(self, index):
        # type: (Any) -> None
        if not isinstance(index, slice):
            return self.__delsingleitem__(index)
        # nprint(index.start, index.stop, index.step, index.indices(len(self)))
        for i in reversed(range(*index.indices(len(self)))):
            del self[i]

    @abstractmethod
    def __getsingleitem__(self, index):
        # type: (Any) -> Any
        raise IndexError

    @abstractmethod
    def __setsingleitem__(self, index, value):
        # type: (Any, Any) -> None
        raise IndexError

    @abstractmethod
    def __delsingleitem__(self, index):
        # type: (Any) -> None
        raise IndexError
