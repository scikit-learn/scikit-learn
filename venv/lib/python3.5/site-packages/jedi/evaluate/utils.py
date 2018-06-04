""" A universal module with functions / classes without dependencies. """
import sys
import contextlib
import functools
import re
import os

from jedi._compatibility import reraise


_sep = os.path.sep
if os.path.altsep is not None:
    _sep += os.path.altsep
_path_re = re.compile('(?:\.[^{0}]+|[{0}]__init__\.py)$'.format(re.escape(_sep)))
del _sep


def to_list(func):
    def wrapper(*args, **kwargs):
        return list(func(*args, **kwargs))
    return wrapper


def unite(iterable):
    """Turns a two dimensional array into a one dimensional."""
    return set(typ for types in iterable for typ in types)


class UncaughtAttributeError(Exception):
    """
    Important, because `__getattr__` and `hasattr` catch AttributeErrors
    implicitly. This is really evil (mainly because of `__getattr__`).
    `hasattr` in Python 2 is even more evil, because it catches ALL exceptions.
    Therefore this class originally had to be derived from `BaseException`
    instead of `Exception`.  But because I removed relevant `hasattr` from
    the code base, we can now switch back to `Exception`.

    :param base: return values of sys.exc_info().
    """


def safe_property(func):
    return property(reraise_uncaught(func))


def reraise_uncaught(func):
    """
    Re-throw uncaught `AttributeError`.

    Usage:  Put ``@rethrow_uncaught`` in front of the function
    which does **not** suppose to raise `AttributeError`.

    AttributeError is easily get caught by `hasattr` and another
    ``except AttributeError`` clause.  This becomes problem when you use
    a lot of "dynamic" attributes (e.g., using ``@property``) because you
    can't distinguish if the property does not exist for real or some code
    inside of the "dynamic" attribute through that error.  In a well
    written code, such error should not exist but getting there is very
    difficult.  This decorator is to help us getting there by changing
    `AttributeError` to `UncaughtAttributeError` to avoid unexpected catch.
    This helps us noticing bugs earlier and facilitates debugging.

    .. note:: Treating StopIteration here is easy.
              Add that feature when needed.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwds):
        try:
            return func(*args, **kwds)
        except AttributeError:
            exc_info = sys.exc_info()
            reraise(UncaughtAttributeError(exc_info[1]), exc_info[2])
    return wrapper


class PushBackIterator(object):
    def __init__(self, iterator):
        self.pushes = []
        self.iterator = iterator
        self.current = None

    def push_back(self, value):
        self.pushes.append(value)

    def __iter__(self):
        return self

    def next(self):
        """ Python 2 Compatibility """
        return self.__next__()

    def __next__(self):
        if self.pushes:
            self.current = self.pushes.pop()
        else:
            self.current = next(self.iterator)
        return self.current


@contextlib.contextmanager
def ignored(*exceptions):
    """
    Context manager that ignores all of the specified exceptions. This will
    be in the standard library starting with Python 3.4.
    """
    try:
        yield
    except exceptions:
        pass


def indent_block(text, indention='    '):
    """This function indents a text block with a default of four spaces."""
    temp = ''
    while text and text[-1] == '\n':
        temp += text[-1]
        text = text[:-1]
    lines = text.split('\n')
    return '\n'.join(map(lambda s: indention + s, lines)) + temp


def dotted_from_fs_path(fs_path, sys_path):
    """
    Changes `/usr/lib/python3.4/email/utils.py` to `email.utils`.  I.e.
    compares the path with sys.path and then returns the dotted_path. If the
    path is not in the sys.path, just returns None.
    """
    if os.path.basename(fs_path).startswith('__init__.'):
        # We are calculating the path. __init__ files are not interesting.
        fs_path = os.path.dirname(fs_path)

    # prefer
    #   - UNIX
    #     /path/to/pythonX.Y/lib-dynload
    #     /path/to/pythonX.Y/site-packages
    #   - Windows
    #     C:\path\to\DLLs
    #     C:\path\to\Lib\site-packages
    # over
    #   - UNIX
    #     /path/to/pythonX.Y
    #   - Windows
    #     C:\path\to\Lib
    path = ''
    for s in sys_path:
        if (fs_path.startswith(s) and len(path) < len(s)):
            path = s

    # - Window
    # X:\path\to\lib-dynload/datetime.pyd => datetime
    module_path = fs_path[len(path):].lstrip(os.path.sep).lstrip('/')
    # - Window
    # Replace like X:\path\to\something/foo/bar.py
    return _path_re.sub('', module_path).replace(os.path.sep, '.').replace('/', '.')
