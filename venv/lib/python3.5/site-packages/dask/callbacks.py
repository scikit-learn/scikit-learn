from __future__ import absolute_import, division, print_function

from contextlib import contextmanager

from .context import _globals

__all__ = ['Callback', 'add_callbacks']


class Callback(object):
    """ Base class for using the callback mechanism

    Create a callback with functions of the following signatures:

    >>> def start(dsk):
    ...     pass
    >>> def start_state(dsk, state):
    ...     pass
    >>> def pretask(key, dsk, state):
    ...     pass
    >>> def posttask(key, result, dsk, state, worker_id):
    ...     pass
    >>> def finish(dsk, state, failed):
    ...     pass

    You may then construct a callback object with any number of them

    >>> cb = Callback(pretask=pretask, finish=finish)  # doctest: +SKIP

    And use it either as a context manager over a compute/get call

    >>> with cb:  # doctest: +SKIP
    ...     x.compute()  # doctest: +SKIP

    Or globally with the ``register`` method

    >>> cb.register()  # doctest: +SKIP
    >>> cb.unregister()  # doctest: +SKIP

    Alternatively subclass the ``Callback`` class with your own methods.

    >>> class PrintKeys(Callback):
    ...     def _pretask(self, key, dask, state):
    ...         print("Computing: {0}!".format(repr(key)))

    >>> with PrintKeys():  # doctest: +SKIP
    ...     x.compute()  # doctest: +SKIP
    """

    def __init__(self, start=None, start_state=None, pretask=None, posttask=None, finish=None):
        if start:
            self._start = start
        if start_state:
            self._start_state = start_state
        if pretask:
            self._pretask = pretask
        if posttask:
            self._posttask = posttask
        if finish:
            self._finish = finish

    @property
    def _callback(self):
        fields = ['_start', '_start_state', '_pretask', '_posttask', '_finish']
        return tuple(getattr(self, i, None) for i in fields)

    def __enter__(self):
        self._cm = add_callbacks(self)
        self._cm.__enter__()
        return self

    def __exit__(self, *args):
        self._cm.__exit__(*args)

    def register(self):
        _globals['callbacks'].add(self._callback)

    def unregister(self):
        _globals['callbacks'].remove(self._callback)


def unpack_callbacks(cbs):
    """Take an iterable of callbacks, return a list of each callback."""
    if cbs:
        return [[i for i in f if i] for f in zip(*cbs)]
    else:
        return [(), (), (), (), ()]


@contextmanager
def local_callbacks(callbacks=None):
    """Allows callbacks to work with nested schedulers.

    Callbacks will only be used by the first started scheduler they encounter.
    This means that only the outermost scheduler will use global callbacks."""
    global_callbacks = callbacks is None
    if global_callbacks:
        callbacks, _globals['callbacks'] = _globals['callbacks'], set()
    try:
        yield callbacks or ()
    finally:
        if global_callbacks:
            _globals['callbacks'] = callbacks


def normalize_callback(cb):
    """Normalizes a callback to a tuple"""
    if isinstance(cb, Callback):
        return cb._callback
    elif isinstance(cb, tuple):
        return cb
    else:
        raise TypeError("Callbacks must be either `Callback` or `tuple`")


class add_callbacks(object):
    """Context manager for callbacks.

    Takes several callbacks and applies them only in the enclosed context.
    Callbacks can either be represented as a ``Callback`` object, or as a tuple
    of length 4.

    Examples
    --------
    >>> def pretask(key, dsk, state):
    ...     print("Now running {0}").format(key)
    >>> callbacks = (None, pretask, None, None)
    >>> with add_callbacks(callbacks):    # doctest: +SKIP
    ...     res.compute()
    """
    def __init__(self, *callbacks):
        self.callbacks = [normalize_callback(c) for c in callbacks]
        _globals['callbacks'].update(self.callbacks)

    def __enter__(self):
        return

    def __exit__(self, type, value, traceback):
        for c in self.callbacks:
            _globals['callbacks'].discard(c)
