import sys
import warnings


def saveit(func):
    """A decorator that caches the return value of a function"""

    name = '_' + func.__name__

    def _wrapper(self, *args, **kwds):
        if not hasattr(self, name):
            setattr(self, name, func(self, *args, **kwds))
        return getattr(self, name)
    return _wrapper

cacheit = saveit


def prevent_recursion(default):
    """A decorator that returns the return value of `default` in recursions"""
    def decorator(func):
        name = '_calling_%s_' % func.__name__

        def newfunc(self, *args, **kwds):
            if getattr(self, name, False):
                return default()
            setattr(self, name, True)
            try:
                return func(self, *args, **kwds)
            finally:
                setattr(self, name, False)
        return newfunc
    return decorator


def ignore_exception(exception_class):
    """A decorator that ignores `exception_class` exceptions"""
    def _decorator(func):
        def newfunc(*args, **kwds):
            try:
                return func(*args, **kwds)
            except exception_class:
                pass
        return newfunc
    return _decorator


def deprecated(message=None):
    """A decorator for deprecated functions"""
    def _decorator(func, message=message):
        if message is None:
            message = '%s is deprecated' % func.__name__

        def newfunc(*args, **kwds):
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return func(*args, **kwds)
        return newfunc
    return _decorator


def cached(size):
    """A caching decorator based on parameter objects"""
    def decorator(func):
        cached_func = _Cached(func, size)
        return lambda *a, **kw: cached_func(*a, **kw)
    return decorator


class _Cached(object):

    def __init__(self, func, count):
        self.func = func
        self.cache = []
        self.count = count

    def __call__(self, *args, **kwds):
        key = (args, kwds)
        for cached_key, cached_result in self.cache:
            if cached_key == key:
                return cached_result
        result = self.func(*args, **kwds)
        self.cache.append((key, result))
        if len(self.cache) > self.count:
            del self.cache[0]
        return result


def resolve(str_or_obj):
    """Returns object from string"""
    from rope.base.utils.pycompat import string_types
    if not isinstance(str_or_obj, string_types):
        return str_or_obj
    if '.' not in str_or_obj:
        str_or_obj += '.'
    mod_name, obj_name = str_or_obj.rsplit('.', 1)
    __import__(mod_name)
    mod = sys.modules[mod_name]
    return getattr(mod, obj_name) if obj_name else mod
