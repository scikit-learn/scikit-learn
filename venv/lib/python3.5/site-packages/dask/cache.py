from __future__ import absolute_import, division, print_function

from .callbacks import Callback
from timeit import default_timer
from numbers import Number
import sys

overhead = sys.getsizeof(1.23) * 4 + sys.getsizeof(()) * 4


class Cache(Callback):
    """ Use cache for computation

    Examples
    --------

    >>> cache = Cache(1e9)          # doctest: +SKIP

    The cache can be used locally as a context manager around ``compute`` or
    ``get`` calls:

    >>> with cache:                 # doctest: +SKIP
    ...     result = x.compute()

    You can also register a cache globally, so that it works for all
    computations:

    >>> cache.register()            # doctest: +SKIP
    >>> cache.unregister()          # doctest: +SKIP
    """

    def __init__(self, cache, *args, **kwargs):
        try:
            import cachey
        except ImportError as ex:
            raise ImportError('Cache requires cachey, "{ex}" problem '
                              'importing'.format(ex=str(ex)))
        self._nbytes = cachey.nbytes
        if isinstance(cache, Number):
            cache = cachey.Cache(cache, *args, **kwargs)
        else:
            assert not args and not kwargs
        self.cache = cache
        self.starttimes = dict()

    def _start(self, dsk):
        self.durations = dict()
        overlap = set(dsk) & set(self.cache.data)
        for key in overlap:
            dsk[key] = self.cache.data[key]

    def _pretask(self, key, dsk, state):
        self.starttimes[key] = default_timer()

    def _posttask(self, key, value, dsk, state, id):
        duration = default_timer() - self.starttimes[key]
        deps = state['dependencies'][key]
        if deps:
            duration += max(self.durations.get(k, 0) for k in deps)
        self.durations[key] = duration
        nb = self._nbytes(value) + overhead + sys.getsizeof(key) * 4
        self.cache.put(key, value, cost=duration / nb / 1e9, nbytes=nb)

    def _finish(self, dsk, state, errored):
        self.starttimes.clear()
        self.durations.clear()
