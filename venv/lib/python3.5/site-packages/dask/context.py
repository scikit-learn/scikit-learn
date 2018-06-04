"""
Control global computation context
"""
from __future__ import absolute_import, division, print_function

import threading
from functools import partial
from collections import defaultdict

_globals = defaultdict(lambda: None)
_globals['callbacks'] = set()


thread_state = threading.local()


class set_options(object):
    """ Set global state within controlled context

    This lets you specify various global settings in a tightly controlled
    ``with`` block.

    Valid keyword arguments currently include the following::

        get - the scheduler to use
        pool - a thread or process pool
        cache - Cache to use for intermediate results
        func_loads/func_dumps - loads/dumps functions for serialization of data
            likely to contain functions.  Defaults to
            cloudpickle.loads/cloudpickle.dumps
        optimizations - List of additional optimizations to run

    Examples
    --------
    >>> with set_options(get=dask.get):  # doctest: +SKIP
    ...     x = np.array(x)  # uses dask.get internally
    """
    def __init__(self, **kwargs):
        self.old = _globals.copy()
        _globals.update(kwargs)

    def __enter__(self):
        return

    def __exit__(self, type, value, traceback):
        _globals.clear()
        _globals.update(self.old)


def globalmethod(default=None, key=None, falsey=None):
    """ Allow function to be taken over by globals

    This modifies a method so that occurrences of it may be taken over by
    functions registered in the global options. Can be used as a decorator or a
    function.

    Parameters
    ----------
    default : callable
        The default callable to use.
    key : str
        Key under which we register this function in the global parameters
    falsey : callable, None, optional
        A function to use if the option is falsey. If not provided, the default
        is used instead.

    Examples
    --------
    >>> import dask
    >>> class Foo(object):
    ...     @globalmethod(key='bar', falsey=lambda: 3)
    ...     def bar():
    ...         return 1
    >>> f = Foo()
    >>> f.bar()
    1
    >>> with dask.set_options(bar=lambda: 2):
    ...     print(f.bar())
    2
    >>> with dask.set_options(bar=False):
    ...     print(f.bar())
    3
    """
    if default is None:
        return partial(globalmethod, key=key, falsey=falsey)
    return GlobalMethod(default=default, key=key, falsey=falsey)


class GlobalMethod(object):
    def __init__(self, default, key, falsey=None):
        self._default = default
        self._key = key
        self._falsey = falsey

    def __get__(self, instance, owner=None):
        if self._key in _globals:
            if _globals[self._key]:
                return _globals[self._key]
            elif self._falsey is not None:
                return self._falsey
        return self._default
