"""
A threaded shared-memory scheduler

See scheduler.py
"""
from __future__ import absolute_import, division, print_function

import sys
from collections import defaultdict
from multiprocessing.pool import ThreadPool
import threading
from threading import current_thread, Lock

from .local import get_async
from .context import _globals
from .utils_test import inc, add  # noqa: F401


def _thread_get_id():
    return current_thread().ident


main_thread = current_thread()
default_pool = None
pools = defaultdict(dict)
pools_lock = Lock()


def pack_exception(e, dumps):
    return e, sys.exc_info()[2]


def get(dsk, result, cache=None, num_workers=None, **kwargs):
    """ Threaded cached implementation of dask.get

    Parameters
    ----------

    dsk: dict
        A dask dictionary specifying a workflow
    result: key or list of keys
        Keys corresponding to desired data
    num_workers: integer of thread count
        The number of threads to use in the ThreadPool that will actually execute tasks
    cache: dict-like (optional)
        Temporary storage of results

    Examples
    --------

    >>> dsk = {'x': 1, 'y': 2, 'z': (inc, 'x'), 'w': (add, 'z', 'y')}
    >>> get(dsk, 'w')
    4
    >>> get(dsk, ['w', 'y'])
    (4, 2)
    """
    global default_pool
    pool = _globals['pool']
    thread = current_thread()

    with pools_lock:
        if pool is None:
            if num_workers is None and thread is main_thread:
                if default_pool is None:
                    default_pool = ThreadPool()
                pool = default_pool
            elif thread in pools and num_workers in pools[thread]:
                pool = pools[thread][num_workers]
            else:
                pool = ThreadPool(num_workers)
                pools[thread][num_workers] = pool

    results = get_async(pool.apply_async, len(pool._pool), dsk, result,
                        cache=cache, get_id=_thread_get_id,
                        pack_exception=pack_exception, **kwargs)

    # Cleanup pools associated to dead threads
    with pools_lock:
        active_threads = set(threading.enumerate())
        if thread is not main_thread:
            for t in list(pools):
                if t not in active_threads:
                    for p in pools.pop(t).values():
                        p.close()

    return results
