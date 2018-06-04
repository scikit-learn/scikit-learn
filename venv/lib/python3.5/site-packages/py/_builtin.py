import sys

try:
    reversed = reversed
except NameError:
    def reversed(sequence):
        """reversed(sequence) -> reverse iterator over values of the sequence

        Return a reverse iterator
        """
        if hasattr(sequence, '__reversed__'):
            return sequence.__reversed__()
        if not hasattr(sequence, '__getitem__'):
            raise TypeError("argument to reversed() must be a sequence")
        return reversed_iterator(sequence)

    class reversed_iterator(object):

        def __init__(self, seq):
            self.seq = seq
            self.remaining = len(seq)

        def __iter__(self):
            return self

        def next(self):
            i = self.remaining
            if i > 0:
                i -= 1
                item = self.seq[i]
                self.remaining = i
                return item
            raise StopIteration

        def __length_hint__(self):
            return self.remaining

try:
    any = any
except NameError:
    def any(iterable):
        for x in iterable:
            if x:
                return True
        return False

try:
    all = all
except NameError:
    def all(iterable):
        for x in iterable:
            if not x:
                return False
        return True

try:
    sorted = sorted
except NameError:
    builtin_cmp = cmp # need to use cmp as keyword arg

    def sorted(iterable, cmp=None, key=None, reverse=0):
        use_cmp = None
        if key is not None:
            if cmp is None:
                def use_cmp(x, y):
                    return builtin_cmp(x[0], y[0])
            else:
                def use_cmp(x, y):
                    return cmp(x[0], y[0])
            l = [(key(element), element) for element in iterable]
        else:
            if cmp is not None:
                use_cmp = cmp
            l = list(iterable)
        if use_cmp is not None:
            l.sort(use_cmp)
        else:
            l.sort()
        if reverse:
            l.reverse()
        if key is not None:
            return [element for (_, element) in l]
        return l

try:
    set, frozenset = set, frozenset
except NameError:
    from sets import set, frozenset

# pass through
enumerate = enumerate

try:
    BaseException = BaseException
except NameError:
    BaseException = Exception

try:
    GeneratorExit = GeneratorExit
except NameError:
    class GeneratorExit(Exception):
        """ This exception is never raised, it is there to make it possible to
        write code compatible with CPython 2.5 even in lower CPython
        versions."""
        pass
    GeneratorExit.__module__ = 'exceptions'

_sysex = (KeyboardInterrupt, SystemExit, MemoryError, GeneratorExit)

try:
    callable = callable
except NameError:
    def callable(obj):
        return hasattr(obj, "__call__")

if sys.version_info >= (3, 0):
    exec ("print_ = print ; exec_=exec")
    import builtins

    # some backward compatibility helpers
    _basestring = str
    def _totext(obj, encoding=None, errors=None):
        if isinstance(obj, bytes):
            if errors is None:
                obj = obj.decode(encoding)
            else:
                obj = obj.decode(encoding, errors)
        elif not isinstance(obj, str):
            obj = str(obj)
        return obj

    def _isbytes(x):
        return isinstance(x, bytes)
    def _istext(x):
        return isinstance(x, str)

    text = str
    bytes = bytes


    def _getimself(function):
        return getattr(function, '__self__', None)

    def _getfuncdict(function):
        return getattr(function, "__dict__", None)

    def _getcode(function):
        return getattr(function, "__code__", None)

    def execfile(fn, globs=None, locs=None):
        if globs is None:
            back = sys._getframe(1)
            globs = back.f_globals
            locs = back.f_locals
            del back
        elif locs is None:
            locs = globs
        fp = open(fn, "r")
        try:
            source = fp.read()
        finally:
            fp.close()
        co = compile(source, fn, "exec", dont_inherit=True)
        exec_(co, globs, locs)

else:
    import __builtin__ as builtins
    _totext = unicode
    _basestring = basestring
    text = unicode
    bytes = str
    execfile = execfile
    callable = callable
    def _isbytes(x):
        return isinstance(x, str)
    def _istext(x):
        return isinstance(x, unicode)

    def _getimself(function):
        return getattr(function, 'im_self', None)

    def _getfuncdict(function):
        return getattr(function, "__dict__", None)

    def _getcode(function):
        try:
            return getattr(function, "__code__")
        except AttributeError:
            return getattr(function, "func_code", None)

    def print_(*args, **kwargs):
        """ minimal backport of py3k print statement. """
        sep = ' '
        if 'sep' in kwargs:
            sep = kwargs.pop('sep')
        end = '\n'
        if 'end' in kwargs:
            end = kwargs.pop('end')
        file = 'file' in kwargs and kwargs.pop('file') or sys.stdout
        if kwargs:
            args = ", ".join([str(x) for x in kwargs])
            raise TypeError("invalid keyword arguments: %s" % args)
        at_start = True
        for x in args:
            if not at_start:
                file.write(sep)
            file.write(str(x))
            at_start = False
        file.write(end)

    def exec_(obj, globals=None, locals=None):
        """ minimal backport of py3k exec statement. """
        __tracebackhide__ = True
        if globals is None:
            frame = sys._getframe(1)
            globals = frame.f_globals
            if locals is None:
                locals = frame.f_locals
        elif locals is None:
            locals = globals
        exec2(obj, globals, locals)

if sys.version_info >= (3, 0):
    def _reraise(cls, val, tb):
        __tracebackhide__ = True
        assert hasattr(val, '__traceback__')
        raise cls.with_traceback(val, tb)
else:
    exec ("""
def _reraise(cls, val, tb):
    __tracebackhide__ = True
    raise cls, val, tb
def exec2(obj, globals, locals):
    __tracebackhide__ = True
    exec obj in globals, locals
""")

def _tryimport(*names):
    """ return the first successfully imported module. """
    assert names
    for name in names:
        try:
            __import__(name)
        except ImportError:
            excinfo = sys.exc_info()
        else:
            return sys.modules[name]
    _reraise(*excinfo)
