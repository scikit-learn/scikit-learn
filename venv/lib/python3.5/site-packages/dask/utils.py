from __future__ import absolute_import, division, print_function

import functools
import inspect
import os
import shutil
import sys
import tempfile
import re
from errno import ENOENT
from collections import Iterator
from contextlib import contextmanager
from importlib import import_module
from numbers import Integral
from threading import Lock
import multiprocessing as mp
import uuid
from weakref import WeakValueDictionary

from .compatibility import get_named_args, getargspec, PY3, unicode, bind_method
from .core import get_deps
from .context import _globals
from .optimization import key_split    # noqa: F401


system_encoding = sys.getdefaultencoding()
if system_encoding == 'ascii':
    system_encoding = 'utf-8'


def deepmap(func, *seqs):
    """ Apply function inside nested lists

    >>> inc = lambda x: x + 1
    >>> deepmap(inc, [[1, 2], [3, 4]])
    [[2, 3], [4, 5]]

    >>> add = lambda x, y: x + y
    >>> deepmap(add, [[1, 2], [3, 4]], [[10, 20], [30, 40]])
    [[11, 22], [33, 44]]
    """
    if isinstance(seqs[0], (list, Iterator)):
        return [deepmap(func, *items) for items in zip(*seqs)]
    else:
        return func(*seqs)


def homogeneous_deepmap(func, seq):
    if not seq:
        return seq
    n = 0
    tmp = seq
    while isinstance(tmp, list):
        n += 1
        tmp = tmp[0]

    return ndeepmap(n, func, seq)


def ndeepmap(n, func, seq):
    """ Call a function on every element within a nested container

    >>> def inc(x):
    ...     return x + 1
    >>> L = [[1, 2], [3, 4, 5]]
    >>> ndeepmap(2, inc, L)
    [[2, 3], [4, 5, 6]]
    """
    if n == 1:
        return [func(item) for item in seq]
    elif n > 1:
        return [ndeepmap(n - 1, func, item) for item in seq]
    elif isinstance(seq, list):
        return func(seq[0])
    else:
        return func(seq)


@contextmanager
def ignoring(*exceptions):
    try:
        yield
    except exceptions:
        pass


def import_required(mod_name, error_msg):
    """Attempt to import a required dependency.

    Raises a RuntimeError if the requested module is not available.
    """
    try:
        return import_module(mod_name)
    except ImportError:
        raise RuntimeError(error_msg)


@contextmanager
def tmpfile(extension='', dir=None):
    extension = '.' + extension.lstrip('.')
    handle, filename = tempfile.mkstemp(extension, dir=dir)
    os.close(handle)
    os.remove(filename)

    try:
        yield filename
    finally:
        if os.path.exists(filename):
            if os.path.isdir(filename):
                shutil.rmtree(filename)
            else:
                with ignoring(OSError):
                    os.remove(filename)


@contextmanager
def tmpdir(dir=None):
    dirname = tempfile.mkdtemp(dir=dir)

    try:
        yield dirname
    finally:
        if os.path.exists(dirname):
            if os.path.isdir(dirname):
                with ignoring(OSError):
                    shutil.rmtree(dirname)
            else:
                with ignoring(OSError):
                    os.remove(dirname)


@contextmanager
def filetext(text, extension='', open=open, mode='w'):
    with tmpfile(extension=extension) as filename:
        f = open(filename, mode=mode)
        try:
            f.write(text)
        finally:
            try:
                f.close()
            except AttributeError:
                pass

        yield filename


@contextmanager
def changed_cwd(new_cwd):
    old_cwd = os.getcwd()
    os.chdir(new_cwd)
    try:
        yield
    finally:
        os.chdir(old_cwd)


@contextmanager
def tmp_cwd(dir=None):
    with tmpdir(dir) as dirname:
        with changed_cwd(dirname):
            yield dirname


@contextmanager
def noop_context():
    yield


class IndexCallable(object):
    """ Provide getitem syntax for functions

    >>> def inc(x):
    ...     return x + 1

    >>> I = IndexCallable(inc)
    >>> I[3]
    4
    """
    __slots__ = 'fn',

    def __init__(self, fn):
        self.fn = fn

    def __getitem__(self, key):
        return self.fn(key)


@contextmanager
def filetexts(d, open=open, mode='t', use_tmpdir=True):
    """ Dumps a number of textfiles to disk

    d - dict
        a mapping from filename to text like {'a.csv': '1,1\n2,2'}

    Since this is meant for use in tests, this context manager will
    automatically switch to a temporary current directory, to avoid
    race conditions when running tests in parallel.
    """
    with (tmp_cwd() if use_tmpdir else noop_context()):
        for filename, text in d.items():
            f = open(filename, 'w' + mode)
            try:
                f.write(text)
            finally:
                try:
                    f.close()
                except AttributeError:
                    pass

        yield list(d)

        for filename in d:
            if os.path.exists(filename):
                with ignoring(OSError):
                    os.remove(filename)


def concrete(seq):
    """ Make nested iterators concrete lists

    >>> data = [[1, 2], [3, 4]]
    >>> seq = iter(map(iter, data))
    >>> concrete(seq)
    [[1, 2], [3, 4]]
    """
    if isinstance(seq, Iterator):
        seq = list(seq)
    if isinstance(seq, (tuple, list)):
        seq = list(map(concrete, seq))
    return seq


def pseudorandom(n, p, random_state=None):
    """ Pseudorandom array of integer indexes

    >>> pseudorandom(5, [0.5, 0.5], random_state=123)
    array([1, 0, 0, 1, 1], dtype=int8)

    >>> pseudorandom(10, [0.5, 0.2, 0.2, 0.1], random_state=5)
    array([0, 2, 0, 3, 0, 1, 2, 1, 0, 0], dtype=int8)
    """
    import numpy as np
    p = list(p)
    cp = np.cumsum([0] + p)
    assert np.allclose(1, cp[-1])
    assert len(p) < 256

    if not isinstance(random_state, np.random.RandomState):
        random_state = np.random.RandomState(random_state)

    x = random_state.random_sample(n)
    out = np.empty(n, dtype='i1')

    for i, (low, high) in enumerate(zip(cp[:-1], cp[1:])):
        out[(x >= low) & (x < high)] = i
    return out


def random_state_data(n, random_state=None):
    """Return a list of arrays that can initialize
    ``np.random.RandomState``.

    Parameters
    ----------
    n : int
        Number of arrays to return.
    random_state : int or np.random.RandomState, optional
        If an int, is used to seed a new ``RandomState``.
    """
    import numpy as np

    if not isinstance(random_state, np.random.RandomState):
        random_state = np.random.RandomState(random_state)

    random_data = random_state.bytes(624 * n * 4)  # `n * 624` 32-bit integers
    l = list(np.frombuffer(random_data, dtype=np.uint32).reshape((n, -1)))
    assert len(l) == n
    return l


def is_integer(i):
    """
    >>> is_integer(6)
    True
    >>> is_integer(42.0)
    True
    >>> is_integer('abc')
    False
    """
    return isinstance(i, Integral) or (isinstance(i, float) and i.is_integer())


ONE_ARITY_BUILTINS = set([abs, all, any, bool, bytearray, bytes, callable, chr,
                          classmethod, complex, dict, dir, enumerate, eval,
                          float, format, frozenset, hash, hex, id, int, iter,
                          len, list, max, min, next, oct, open, ord, range,
                          repr, reversed, round, set, slice, sorted,
                          staticmethod, str, sum, tuple,
                          type, vars, zip, memoryview])
if PY3:
    ONE_ARITY_BUILTINS.add(ascii)  # noqa: F821
MULTI_ARITY_BUILTINS = set([compile, delattr, divmod, filter, getattr, hasattr,
                            isinstance, issubclass, map, pow, setattr])


def takes_multiple_arguments(func, varargs=True):
    """ Does this function take multiple arguments?

    >>> def f(x, y): pass
    >>> takes_multiple_arguments(f)
    True

    >>> def f(x): pass
    >>> takes_multiple_arguments(f)
    False

    >>> def f(x, y=None): pass
    >>> takes_multiple_arguments(f)
    False

    >>> def f(*args): pass
    >>> takes_multiple_arguments(f)
    True

    >>> class Thing(object):
    ...     def __init__(self, a): pass
    >>> takes_multiple_arguments(Thing)
    False

    """
    if func in ONE_ARITY_BUILTINS:
        return False
    elif func in MULTI_ARITY_BUILTINS:
        return True

    try:
        spec = getargspec(func)
    except Exception:
        return False

    try:
        is_constructor = spec.args[0] == 'self' and isinstance(func, type)
    except Exception:
        is_constructor = False

    if varargs and spec.varargs:
        return True

    ndefaults = 0 if spec.defaults is None else len(spec.defaults)
    return len(spec.args) - ndefaults - is_constructor > 1


class Dispatch(object):
    """Simple single dispatch."""
    def __init__(self, name=None):
        self._lookup = {}
        self._lazy = {}
        if name:
            self.__name__ = name

    def register(self, type, func=None):
        """Register dispatch of `func` on arguments of type `type`"""
        def wrapper(func):
            if isinstance(type, tuple):
                for t in type:
                    self.register(t, func)
            else:
                self._lookup[type] = func
            return func

        return wrapper(func) if func is not None else wrapper

    def register_lazy(self, toplevel, func=None):
        """
        Register a registration function which will be called if the
        *toplevel* module (e.g. 'pandas') is ever loaded.
        """
        def wrapper(func):
            self._lazy[toplevel] = func
            return func

        return wrapper(func) if func is not None else wrapper

    def dispatch(self, cls):
        """Return the function implementation for the given ``cls``"""
        # Fast path with direct lookup on cls
        lk = self._lookup
        try:
            impl = lk[cls]
        except KeyError:
            pass
        else:
            return impl
        # Is a lazy registration function present?
        toplevel, _, _ = cls.__module__.partition('.')
        try:
            register = self._lazy.pop(toplevel)
        except KeyError:
            pass
        else:
            register()
            return self.dispatch(cls) # recurse
        # Walk the MRO and cache the lookup result
        for cls2 in inspect.getmro(cls)[1:]:
            if cls2 in lk:
                lk[cls] = lk[cls2]
                return lk[cls2]
        raise TypeError("No dispatch for {0}".format(cls))

    def __call__(self, arg):
        """
        Call the corresponding method based on type of argument.
        """
        meth = self.dispatch(type(arg))
        return meth(arg)


def ensure_not_exists(filename):
    """
    Ensure that a file does not exist.
    """
    try:
        os.unlink(filename)
    except OSError as e:
        if e.errno != ENOENT:
            raise


def _skip_doctest(line):
    # NumPy docstring contains cursor and comment only example
    stripped = line.strip()
    if stripped == '>>>' or stripped.startswith('>>> #'):
        return stripped
    elif '>>>' in stripped and '+SKIP' not in stripped:
        return line + '  # doctest: +SKIP'
    else:
        return line


def skip_doctest(doc):
    if doc is None:
        return ''
    return '\n'.join([_skip_doctest(line) for line in doc.split('\n')])


def extra_titles(doc):
    lines = doc.split('\n')
    titles = {i: lines[i].strip() for i in range(len(lines) - 1)
              if lines[i + 1] and all(c == '-' for c in lines[i + 1].strip())}

    seen = set()
    for i, title in sorted(titles.items()):
        if title in seen:
            new_title = 'Extra ' + title
            lines[i] = lines[i].replace(title, new_title)
            lines[i + 1] = lines[i + 1].replace('-' * len(title),
                                                '-' * len(new_title))
        else:
            seen.add(title)

    return '\n'.join(lines)


def derived_from(original_klass, version=None, ua_args=[]):
    """Decorator to attach original class's docstring to the wrapped method.

    Parameters
    ----------
    original_klass: type
        Original class which the method is derived from
    version : str
        Original package version which supports the wrapped method
    ua_args : list
        List of keywords which Dask doesn't support. Keywords existing in
        original but not in Dask will automatically be added.
    """
    def wrapper(method):
        method_name = method.__name__

        try:
            # do not use wraps here, as it hides keyword arguments displayed
            # in the doc
            original_method = getattr(original_klass, method_name)
            doc = original_method.__doc__
            if doc is None:
                doc = ''

            try:
                method_args = get_named_args(method)
                original_args = get_named_args(original_method)
                not_supported = [m for m in original_args if m not in method_args]
            except ValueError:
                not_supported = []

            if len(ua_args) > 0:
                not_supported.extend(ua_args)

            if len(not_supported) > 0:
                note = ("\n        Notes\n        -----\n"
                        "        Dask doesn't support the following argument(s).\n\n")
                args = ''.join(['        * {0}\n'.format(a) for a in not_supported])
                doc = doc + note + args
            doc = skip_doctest(doc)
            doc = extra_titles(doc)

            method.__doc__ = doc
            return method

        except AttributeError:
            module_name = original_klass.__module__.split('.')[0]

            @functools.wraps(method)
            def wrapped(*args, **kwargs):
                msg = "Base package doesn't support '{0}'.".format(method_name)
                if version is not None:
                    msg2 = " Use {0} {1} or later to use this method."
                    msg += msg2.format(module_name, version)
                raise NotImplementedError(msg)
            return wrapped
    return wrapper


def funcname(func):
    """Get the name of a function."""
    # functools.partial
    if isinstance(func, functools.partial):
        return funcname(func.func)
    # methodcaller
    if isinstance(func, methodcaller):
        return func.method

    module_name = getattr(func, '__module__', None) or ''
    type_name = getattr(type(func), '__name__', None) or ''

    # toolz.curry
    if 'toolz' in module_name and 'curry' == type_name:
        return func.func_name
    # multipledispatch objects
    if 'multipledispatch' in module_name and 'Dispatcher' == type_name:
        return func.name

    # All other callables
    try:
        name = func.__name__
        if name == '<lambda>':
            return 'lambda'
        return name
    except AttributeError:
        return str(func)


def ensure_bytes(s):
    """ Turn string or bytes to bytes

    >>> ensure_bytes(u'123')
    '123'
    >>> ensure_bytes('123')
    '123'
    >>> ensure_bytes(b'123')
    '123'
    """
    if isinstance(s, bytes):
        return s
    if hasattr(s, 'encode'):
        return s.encode()
    msg = "Object %s is neither a bytes object nor has an encode method"
    raise TypeError(msg % s)


def ensure_unicode(s):
    """ Turn string or bytes to bytes

    >>> ensure_unicode(u'123')
    u'123'
    >>> ensure_unicode('123')
    u'123'
    >>> ensure_unicode(b'123')
    u'123'
    """
    if isinstance(s, unicode):
        return s
    if hasattr(s, 'decode'):
        return s.decode()
    msg = "Object %s is neither a bytes object nor has an encode method"
    raise TypeError(msg % s)


def digit(n, k, base):
    """

    >>> digit(1234, 0, 10)
    4
    >>> digit(1234, 1, 10)
    3
    >>> digit(1234, 2, 10)
    2
    >>> digit(1234, 3, 10)
    1
    """
    return n // base**k % base


def insert(tup, loc, val):
    """

    >>> insert(('a', 'b', 'c'), 0, 'x')
    ('x', 'b', 'c')
    """
    L = list(tup)
    L[loc] = val
    return tuple(L)


def dependency_depth(dsk):
    import toolz

    deps, _ = get_deps(dsk)

    @toolz.memoize
    def max_depth_by_deps(key):
        if not deps[key]:
            return 1

        d = 1 + max(max_depth_by_deps(dep_key) for dep_key in deps[key])
        return d

    return max(max_depth_by_deps(dep_key) for dep_key in deps.keys())


def memory_repr(num):
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0


def asciitable(columns, rows):
    """Formats an ascii table for given columns and rows.

    Parameters
    ----------
    columns : list
        The column names
    rows : list of tuples
        The rows in the table. Each tuple must be the same length as
        ``columns``.
    """
    rows = [tuple(str(i) for i in r) for r in rows]
    columns = tuple(str(i) for i in columns)
    widths = tuple(max(max(map(len, x)), len(c))
                   for x, c in zip(zip(*rows), columns))
    row_template = ('|' + (' %%-%ds |' * len(columns))) % widths
    header = row_template % tuple(columns)
    bar = '+%s+' % '+'.join('-' * (w + 2) for w in widths)
    data = '\n'.join(row_template % r for r in rows)
    return '\n'.join([bar, header, bar, data, bar])


def put_lines(buf, lines):
    if any(not isinstance(x, unicode) for x in lines):
        lines = [unicode(x) for x in lines]
    buf.write('\n'.join(lines))


_method_cache = {}


class methodcaller(object):
    """
    Return a callable object that calls the given method on its operand.

    Unlike the builtin `operator.methodcaller`, instances of this class are
    serializable
    """

    __slots__ = ('method',)
    func = property(lambda self: self.method)  # For `funcname` to work

    def __new__(cls, method):
        if method in _method_cache:
            return _method_cache[method]
        self = object.__new__(cls)
        self.method = method
        _method_cache[method] = self
        return self

    def __call__(self, obj, *args, **kwargs):
        return getattr(obj, self.method)(*args, **kwargs)

    def __reduce__(self):
        return (methodcaller, (self.method,))

    def __str__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.method)

    __repr__ = __str__


class itemgetter(object):
    """
    Return a callable object that gets an item from the operand

    Unlike the builtin `operator.itemgetter`, instances of this class are
    serializable
    """
    __slots__ = ('index',)

    def __init__(self, index):
        self.index = index

    def __call__(self, x):
        return x[self.index]

    def __reduce__(self):
        return (itemgetter, (self.index,))

    def __eq__(self, other):
        return type(self) is type(other) and self.index == other.index


class MethodCache(object):
    """Attribute access on this object returns a methodcaller for that
    attribute.

    Examples
    --------
    >>> a = [1, 3, 3]
    >>> M.count(a, 3) == a.count(3)
    True
    """
    __getattr__ = staticmethod(methodcaller)
    __dir__ = lambda self: list(_method_cache)


M = MethodCache()


class SerializableLock(object):
    _locks = WeakValueDictionary()
    """ A Serializable per-process Lock

    This wraps a normal ``threading.Lock`` object and satisfies the same
    interface.  However, this lock can also be serialized and sent to different
    processes.  It will not block concurrent operations between processes (for
    this you should look at ``multiprocessing.Lock`` or ``locket.lock_file``
    but will consistently deserialize into the same lock.

    So if we make a lock in one process::

        lock = SerializableLock()

    And then send it over to another process multiple times::

        bytes = pickle.dumps(lock)
        a = pickle.loads(bytes)
        b = pickle.loads(bytes)

    Then the deserialized objects will operate as though they were the same
    lock, and collide as appropriate.

    This is useful for consistently protecting resources on a per-process
    level.

    The creation of locks is itself not threadsafe.
    """
    def __init__(self, token=None):
        self.token = token or str(uuid.uuid4())
        if self.token in SerializableLock._locks:
            self.lock = SerializableLock._locks[self.token]
        else:
            self.lock = Lock()
            SerializableLock._locks[self.token] = self.lock

    def acquire(self, *args):
        return self.lock.acquire(*args)

    def release(self, *args):
        return self.lock.release(*args)

    def __enter__(self):
        self.lock.__enter__()

    def __exit__(self, *args):
        self.lock.__exit__(*args)

    @property
    def locked(self):
        return self.locked

    def __getstate__(self):
        return self.token

    def __setstate__(self, token):
        self.__init__(token)

    def __str__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.token)

    __repr__ = __str__


def effective_get(get=None, collection=None):
    """Get the effective get method used in a given situation"""
    return (get or _globals.get('get') or
            getattr(collection, '__dask_scheduler__', None))


def get_scheduler_lock(get=None, collection=None):
    """Get an instance of the appropriate lock for a certain situation based on
       scheduler used."""
    from . import multiprocessing
    actual_get = effective_get(get, collection)

    if actual_get == multiprocessing.get:
        return mp.Manager().Lock()
    return SerializableLock()


def ensure_dict(d):
    if type(d) is dict:
        return d
    elif hasattr(d, 'dicts'):
        result = {}
        for dd in d.dicts.values():
            result.update(dd)
        return result
    return dict(d)


class OperatorMethodMixin(object):
    """A mixin for dynamically implementing operators"""

    @classmethod
    def _bind_operator(cls, op):
        """ bind operator to this class """
        name = op.__name__

        if name.endswith('_'):
            # for and_ and or_
            name = name[:-1]
        elif name == 'inv':
            name = 'invert'

        meth = '__{0}__'.format(name)

        if name in ('abs', 'invert', 'neg', 'pos'):
            bind_method(cls, meth, cls._get_unary_operator(op))
        else:
            bind_method(cls, meth, cls._get_binary_operator(op))

            if name in ('eq', 'gt', 'ge', 'lt', 'le', 'ne', 'getitem'):
                return

            rmeth = '__r{0}__'.format(name)
            bind_method(cls, rmeth, cls._get_binary_operator(op, inv=True))

    @classmethod
    def _get_unary_operator(cls, op):
        """ Must return a method used by unary operator """
        raise NotImplementedError

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        """ Must return a method used by binary operator """
        raise NotImplementedError


def partial_by_order(*args, **kwargs):
    """

    >>> from operator import add
    >>> partial_by_order(5, function=add, other=[(1, 10)])
    15
    """
    function = kwargs.pop('function')
    other = kwargs.pop('other')
    args2 = list(args)
    for i, arg in other:
        args2.insert(i, arg)
    return function(*args2, **kwargs)


def is_arraylike(x):
    """ Is this object a numpy array or something similar?

    Examples
    --------
    >>> import numpy as np
    >>> x = np.ones(5)
    >>> is_arraylike(x)
    True
    >>> is_arraylike(5)
    False
    >>> is_arraylike('cat')
    False
    """
    return (hasattr(x, '__array__') and
            hasattr(x, 'shape') and x.shape and
            hasattr(x, 'dtype'))


def natural_sort_key(s):
    """
    Sorting `key` function for performing a natural sort on a collection of
    strings

    See https://en.wikipedia.org/wiki/Natural_sort_order

    Parameters
    ----------
    s : str
        A string that is an element of the collection being sorted

    Returns
    -------
    tuple[str or int]
        Tuple of the parts of the input string where each part is either a
        string or an integer

    Examples
    --------
    >>> a = ['f0', 'f1', 'f2', 'f8', 'f9', 'f10', 'f11', 'f19', 'f20', 'f21']
    >>> sorted(a)
    ['f0', 'f1', 'f10', 'f11', 'f19', 'f2', 'f20', 'f21', 'f8', 'f9']
    >>> sorted(a, key=natural_sort_key)
    ['f0', 'f1', 'f2', 'f8', 'f9', 'f10', 'f11', 'f19', 'f20', 'f21']
    """
    return [int(part) if part.isdigit() else part
            for part in re.split('(\d+)', s)]
