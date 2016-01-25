"""
Backends for embarrassingly parallel code.
"""

import gc
import os
import warnings
import threading

from ._multiprocessing_helpers import mp
if mp is not None:
    from .pool import MemmapingPool
    from multiprocessing.pool import ThreadPool

# Environment variables to protect against bad situations when nesting
JOBLIB_SPAWNED_PROCESS = "__JOBLIB_SPAWNED_PARALLEL__"

# Under Linux or OS X the default start method of multiprocessing
# can cause third party libraries to crash. Under Python 3.4+ it is possible
# to set an environment variable to switch the default start method from
# 'fork' to 'forkserver' or 'spawn' to avoid this issue albeit at the cost
# of causing semantic changes and some additional pool instanciation overhead.
if hasattr(mp, 'get_context'):
    method = os.environ.get('JOBLIB_START_METHOD', '').strip() or None
    DEFAULT_MP_CONTEXT = mp.get_context(method=method)
else:
    DEFAULT_MP_CONTEXT = None


###############################################################################
class ThreadingBackend(object):

    def effective_n_jobs(self, n_jobs):
        """ Determine the number of jobs which are going to run in parallel """
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        elif mp is None or n_jobs is None:
            # multiprocessing is not available or disabled, fallback
            # to sequential mode
            return 1
        elif n_jobs < 0:
            n_jobs = max(mp.cpu_count() + 1 + n_jobs, 1)
        return n_jobs

    def initialize(self, n_jobs, poolargs):
        self._pool = ThreadPool(n_jobs)
        return n_jobs

    def terminate(self):
        if self._pool is not None:
            self._pool.close()
            self._pool.terminate()  # terminate does a join()
            self._pool = None

    def apply_async(self, *args, **kwargs):
        return self._pool.apply_async(*args, **kwargs)

    def get_exceptions(self):
        return []


###############################################################################
class MultiProcessingBackend(ThreadingBackend):

    def effective_n_jobs(self, n_jobs):
        """ Determine the number of jobs which are going to run in parallel. Will
        also check if we're attempting to create a nested parallel loop. """
        if mp.current_process().daemon:
            # Daemonic processes cannot have children
            self._pool = None
            warnings.warn(
                'Multiprocessing-backed parallel loops cannot be nested,'
                ' setting n_jobs=1',
                stacklevel=3)
            return 1

        elif threading.current_thread().name != 'MainThread':
            # Prevent posix fork inside in non-main posix threads
            self._pool = None
            warnings.warn(
                'Multiprocessing backed parallel loops cannot be nested'
                ' below threads, setting n_jobs=1',
                stacklevel=3)
            return 1

        return super(MultiProcessingBackend, self).effective_n_jobs(n_jobs)

    def initialize(self, n_jobs, poolargs):
        already_forked = int(os.environ.get(JOBLIB_SPAWNED_PROCESS, 0))
        if already_forked:
            raise ImportError('[joblib] Attempting to do parallel computing '
                    'without protecting your import on a system that does '
                    'not support forking. To use parallel-computing in a '
                    'script, you must protect your main loop using "if '
                    "__name__ == '__main__'"
                    '". Please see the joblib documentation on Parallel '
                    'for more information'
                )
        # Set an environment variable to avoid infinite loops
        os.environ[JOBLIB_SPAWNED_PROCESS] = '1'

        # Make sure to free as much memory as possible before forking
        gc.collect()
        self._pool = MemmapingPool(n_jobs, **poolargs)

        return n_jobs

    def terminate(self):
        super(MultiProcessingBackend, self).terminate()
        if JOBLIB_SPAWNED_PROCESS in os.environ:
            del os.environ[JOBLIB_SPAWNED_PROCESS]

    def get_exceptions(self):
        # We are using multiprocessing, we also want to capture
        # KeyboardInterrupts
        from .parallel import WorkerInterrupt
        return [KeyboardInterrupt, WorkerInterrupt]
