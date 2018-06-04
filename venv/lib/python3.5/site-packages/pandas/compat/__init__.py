"""
compat
======

Cross-compatible functions for Python 2 and 3.

Key items to import for 2/3 compatible code:
* iterators: range(), map(), zip(), filter(), reduce()
* lists: lrange(), lmap(), lzip(), lfilter()
* unicode: u() [no unicode builtin in Python 3]
* longs: long (int in Python 3)
* callable
* iterable method compatibility: iteritems, iterkeys, itervalues
  * Uses the original method if available, otherwise uses items, keys, values.
* types:
    * text_type: unicode in Python 2, str in Python 3
    * binary_type: str in Python 2, bytes in Python 3
    * string_types: basestring in Python 2, str in Python 3
* bind_method: binds functions to classes
* add_metaclass(metaclass) - class decorator that recreates class with with the
  given metaclass instead (and avoids intermediary class creation)

Other items:
* platform checker
"""
# pylint disable=W0611
# flake8: noqa

import re
import functools
import itertools
from distutils.version import LooseVersion
from itertools import product
import sys
import platform
import types
from unicodedata import east_asian_width
import struct
import inspect
from collections import namedtuple

PY2 = sys.version_info[0] == 2
PY3 = (sys.version_info[0] >= 3)
PY35 = (sys.version_info >= (3, 5))
PY36 = (sys.version_info >= (3, 6))
PYPY = (platform.python_implementation() == 'PyPy')

try:
    import __builtin__ as builtins
    # not writeable when instantiated with string, doesn't handle unicode well
    from cStringIO import StringIO as cStringIO
    # always writeable
    from StringIO import StringIO
    BytesIO = StringIO
    import cPickle
    import httplib
except ImportError:
    import builtins
    from io import StringIO, BytesIO
    cStringIO = StringIO
    import pickle as cPickle
    import http.client as httplib

from pandas.compat.chainmap import DeepChainMap


if PY3:
    def isidentifier(s):
        return s.isidentifier()

    def str_to_bytes(s, encoding=None):
        return s.encode(encoding or 'ascii')

    def bytes_to_str(b, encoding=None):
        return b.decode(encoding or 'utf-8')

    # The signature version below is directly copied from Django,
    # https://github.com/django/django/pull/4846
    def signature(f):
        sig = inspect.signature(f)
        args = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
        ]
        varargs = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.VAR_POSITIONAL
        ]
        varargs = varargs[0] if varargs else None
        keywords = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.VAR_KEYWORD
        ]
        keywords = keywords[0] if keywords else None
        defaults = [
            p.default for p in sig.parameters.values()
            if p.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
            and p.default is not p.empty
        ] or None
        argspec = namedtuple('Signature', ['args', 'defaults',
                                           'varargs', 'keywords'])
        return argspec(args, defaults, varargs, keywords)

    def get_range_parameters(data):
        """Gets the start, stop, and step parameters from a range object"""
        return data.start, data.stop, data.step

    # have to explicitly put builtins into the namespace
    range = range
    map = map
    zip = zip
    filter = filter
    intern = sys.intern
    reduce = functools.reduce
    long = int
    unichr = chr

    # This was introduced in Python 3.3, but we don't support
    # Python 3.x < 3.5, so checking PY3 is safe.
    FileNotFoundError = FileNotFoundError

    # list-producing versions of the major Python iterating functions
    def lrange(*args, **kwargs):
        return list(range(*args, **kwargs))

    def lzip(*args, **kwargs):
        return list(zip(*args, **kwargs))

    def lmap(*args, **kwargs):
        return list(map(*args, **kwargs))

    def lfilter(*args, **kwargs):
        return list(filter(*args, **kwargs))

    from importlib import reload
    reload = reload

else:
    # Python 2
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")

    FileNotFoundError = IOError

    def isidentifier(s, dotted=False):
        return bool(_name_re.match(s))

    def str_to_bytes(s, encoding='ascii'):
        return s

    def bytes_to_str(b, encoding='ascii'):
        return b

    def signature(f):
        return inspect.getargspec(f)

    def get_range_parameters(data):
        """Gets the start, stop, and step parameters from a range object"""
        # seems we only have indexing ops to infer
        # rather than direct accessors
        if len(data) > 1:
            step = data[1] - data[0]
            stop = data[-1] + step
            start = data[0]
        elif len(data):
            start = data[0]
            stop = data[0] + 1
            step = 1
        else:
            start = stop = 0
            step = 1

        return start, stop, step

    # import iterator versions of these functions
    range = xrange
    intern = intern
    zip = itertools.izip
    filter = itertools.ifilter
    map = itertools.imap
    reduce = reduce
    long = long
    unichr = unichr

    # Python 2-builtin ranges produce lists
    lrange = builtins.range
    lzip = builtins.zip
    lmap = builtins.map
    lfilter = builtins.filter

    reload = builtins.reload

if PY2:
    def iteritems(obj, **kw):
        return obj.iteritems(**kw)

    def iterkeys(obj, **kw):
        return obj.iterkeys(**kw)

    def itervalues(obj, **kw):
        return obj.itervalues(**kw)

    next = lambda it: it.next()
else:
    def iteritems(obj, **kw):
        return iter(obj.items(**kw))

    def iterkeys(obj, **kw):
        return iter(obj.keys(**kw))

    def itervalues(obj, **kw):
        return iter(obj.values(**kw))

    next = next


def bind_method(cls, name, func):
    """Bind a method to class, python 2 and python 3 compatible.

    Parameters
    ----------

    cls : type
        class to receive bound method
    name : basestring
        name of method on class instance
    func : function
        function to be bound as method


    Returns
    -------
    None
    """
    # only python 2 has bound/unbound method issue
    if not PY3:
        setattr(cls, name, types.MethodType(func, None, cls))
    else:
        setattr(cls, name, func)
# ----------------------------------------------------------------------------
# functions largely based / taken from the six module

# Much of the code in this module comes from Benjamin Peterson's six library.
# The license for this library can be found in LICENSES/SIX and the code can be
# found at https://bitbucket.org/gutworth/six

# Definition of East Asian Width
# http://unicode.org/reports/tr11/
# Ambiguous width can be changed by option
_EAW_MAP = {'Na': 1, 'N': 1, 'W': 2, 'F': 2, 'H': 1}

if PY3:
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str
    binary_type = bytes

    def u(s):
        return s

    def u_safe(s):
        return s

    def to_str(s):
        """
        Convert bytes and non-string into Python 3 str
        """
        if isinstance(s, binary_type):
            s = bytes_to_str(s)
        elif not isinstance(s, string_types):
            s = str(s)
        return s

    def strlen(data, encoding=None):
        # encoding is for compat with PY2
        return len(data)

    def east_asian_len(data, encoding=None, ambiguous_width=1):
        """
        Calculate display width considering unicode East Asian Width
        """
        if isinstance(data, text_type):
            return sum(_EAW_MAP.get(east_asian_width(c), ambiguous_width) for c in data)
        else:
            return len(data)

    def import_lzma():
        """ import lzma from the std library """
        import lzma
        return lzma

    def set_function_name(f, name, cls):
        """ Bind the name/qualname attributes of the function """
        f.__name__ = name
        f.__qualname__ = '{klass}.{name}'.format(
            klass=cls.__name__,
            name=name)
        f.__module__ = cls.__module__
        return f

    ResourceWarning = ResourceWarning

else:
    string_types = basestring,
    integer_types = (int, long)
    class_types = (type, types.ClassType)
    text_type = unicode
    binary_type = str

    def u(s):
        return unicode(s, "unicode_escape")

    def u_safe(s):
        try:
            return unicode(s, "unicode_escape")
        except:
            return s

    def to_str(s):
        """
        Convert unicode and non-string into Python 2 str
        """
        if not isinstance(s, string_types):
            s = str(s)
        return s

    def strlen(data, encoding=None):
        try:
            data = data.decode(encoding)
        except UnicodeError:
            pass
        return len(data)

    def east_asian_len(data, encoding=None, ambiguous_width=1):
        """
        Calculate display width considering unicode East Asian Width
        """
        if isinstance(data, text_type):
            try:
                data = data.decode(encoding)
            except UnicodeError:
                pass
            return sum(_EAW_MAP.get(east_asian_width(c), ambiguous_width) for c in data)
        else:
            return len(data)

    def import_lzma():
        """ import the backported lzma library
        or raise ImportError if not available """
        from backports import lzma
        return lzma

    def set_function_name(f, name, cls):
        """ Bind the name attributes of the function """
        f.__name__ = name
        return f

    class ResourceWarning(Warning):
        pass

string_and_binary_types = string_types + (binary_type,)


try:
    # callable reintroduced in later versions of Python
    callable = callable
except NameError:
    def callable(obj):
        return any("__call__" in klass.__dict__ for klass in type(obj).__mro__)


if PY2:
    # In PY2 functools.wraps doesn't provide metadata pytest needs to generate
    # decorated tests using parametrization. See pytest GH issue #2782
    def wraps(wrapped, assigned=functools.WRAPPER_ASSIGNMENTS,
              updated=functools.WRAPPER_UPDATES):
        def wrapper(f):
            f = functools.wraps(wrapped, assigned, updated)(f)
            f.__wrapped__ = wrapped
            return f
        return wrapper
else:
    wraps = functools.wraps


def add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        for slots_var in orig_vars.get('__slots__', ()):
            orig_vars.pop(slots_var)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper

from collections import OrderedDict, Counter

if PY3:
    def raise_with_traceback(exc, traceback=Ellipsis):
        if traceback == Ellipsis:
            _, _, traceback = sys.exc_info()
        raise exc.with_traceback(traceback)
else:
    # this version of raise is a syntax error in Python 3
    exec("""
def raise_with_traceback(exc, traceback=Ellipsis):
    if traceback == Ellipsis:
        _, _, traceback = sys.exc_info()
    raise exc, None, traceback
""")

raise_with_traceback.__doc__ = """Raise exception with existing traceback.
If traceback is not passed, uses sys.exc_info() to get traceback."""


# dateutil minimum version
import dateutil

if LooseVersion(dateutil.__version__) < LooseVersion('2.5'):
    raise ImportError('dateutil 2.5.0 is the minimum required version')
from dateutil import parser as _date_parser
parse_date = _date_parser.parse


# In Python 3.7, the private re._pattern_type is removed.
# Python 3.5+ have typing.re.Pattern
if PY35:
    import typing
    re_type = typing.re.Pattern
else:
    re_type = type(re.compile(''))

# https://github.com/pandas-dev/pandas/pull/9123
def is_platform_little_endian():
    """ am I little endian """
    return sys.byteorder == 'little'


def is_platform_windows():
    return sys.platform == 'win32' or sys.platform == 'cygwin'


def is_platform_linux():
    return sys.platform == 'linux2'


def is_platform_mac():
    return sys.platform == 'darwin'


def is_platform_32bit():
    return struct.calcsize("P") * 8 < 64
