"""
Backends for embarrassingly parallel code.
"""

import gc
import os
import sys
import warnings
import threading
import functools
import contextlib
from abc import ABCMeta, abstractmethod

from .format_stack import format_exc
from .my_exceptions import WorkerInterrupt, TransportableException
from ._multiprocessing_helpers import mp
from ._compat import with_metaclass, PY27
if mp is not None:
    from .disk import delete_folder
    from .pool import MemmappingPool
    from multiprocessing.pool import ThreadPool
    from .executor import get_memmapping_executor

    # Compat between concurrent.futures and multiprocessing TimeoutError
    from multiprocessing import TimeoutError
    from .externals.loky._base import TimeoutError as LokyTimeoutError
    from .externals.loky import process_executor, cpu_count


class ParallelBackendBase(with_metaclass(ABCMeta)):
    """Helper abc which defines all methods a ParallelBackend must implement"""

    supports_timeout = False
    nesting_level = 0

    def __init__(self, nesting_level=0):
        self.nesting_level = nesting_level

    SUPPORTED_CLIB_VARS = [
        'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'
    ]

    @abstractmethod
    def effective_n_jobs(self, n_jobs):
        """Determine the number of jobs that can actually run in parallel

        n_jobs is the number of workers requested by the callers. Passing
        n_jobs=-1 means requesting all available workers for instance matching
        the number of CPU cores on the worker host(s).

        This method should return a guesstimate of the number of workers that
        can actually perform work concurrently. The primary use case is to make
        it possible for the caller to know in how many chunks to slice the
        work.

        In general working on larger data chunks is more efficient (less
        scheduling overhead and better use of CPU cache prefetching heuristics)
        as long as all the workers have enough work to do.
        """

    @abstractmethod
    def apply_async(self, func, callback=None):
        """Schedule a func to be run"""

    def configure(self, n_jobs=1, parallel=None, prefer=None, require=None,
                  **backend_args):
        """Reconfigure the backend and return the number of workers.

        This makes it possible to reuse an existing backend instance for
        successive independent calls to Parallel with different parameters.
        """
        self.parallel = parallel
        return self.effective_n_jobs(n_jobs)

    def start_call(self):
        """Call-back method called at the beginning of a Parallel call"""

    def stop_call(self):
        """Call-back method called at the end of a Parallel call"""

    def terminate(self):
        """Shutdown the workers and free the shared memory."""

    def compute_batch_size(self):
        """Determine the optimal batch size"""
        return 1

    def batch_completed(self, batch_size, duration):
        """Callback indicate how long it took to run a batch"""

    def get_exceptions(self):
        """List of exception types to be captured."""
        return []

    def abort_everything(self, ensure_ready=True):
        """Abort any running tasks

        This is called when an exception has been raised when executing a tasks
        and all the remaining tasks will be ignored and can therefore be
        aborted to spare computation resources.

        If ensure_ready is True, the backend should be left in an operating
        state as future tasks might be re-submitted via that same backend
        instance.

        If ensure_ready is False, the implementer of this method can decide
        to leave the backend in a closed / terminated state as no new task
        are expected to be submitted to this backend.

        Setting ensure_ready to False is an optimization that can be leveraged
        when aborting tasks via killing processes from a local process pool
        managed by the backend it-self: if we expect no new tasks, there is no
        point in re-creating new workers.
        """
        # Does nothing by default: to be overridden in subclasses when
        # canceling tasks is possible.
        pass

    def get_nested_backend(self):
        """Backend instance to be used by nested Parallel calls.

        By default a thread-based backend is used for the first level of
        nesting. Beyond, switch to sequential backend to avoid spawning too
        many threads on the host.
        """
        nesting_level = getattr(self, 'nesting_level', 0) + 1
        if nesting_level > 1:
            return SequentialBackend(nesting_level=nesting_level)
        else:
            return ThreadingBackend(nesting_level=nesting_level)

    @contextlib.contextmanager
    def retrieval_context(self):
        """Context manager to manage an execution context.

        Calls to Parallel.retrieve will be made inside this context.

        By default, this does nothing. It may be useful for subclasses to
        handle nested parallelism. In particular, it may be required to avoid
        deadlocks if a backend manages a fixed number of workers, when those
        workers may be asked to do nested Parallel calls. Without
        'retrieval_context' this could lead to deadlock, as all the workers
        managed by the backend may be "busy" waiting for the nested parallel
        calls to finish, but the backend has no free workers to execute those
        tasks.
        """
        yield

    @classmethod
    def limit_clib_threads(cls, n_threads=1):
        """Initializer to limit the number of threads used by some C-libraries.

        This function set the number of threads to `n_threads` for OpenMP, MKL,
        Accelerated and OpenBLAS libraries, that can be used with scientific
        computing tools like numpy.
        """
        for var in cls.SUPPORTED_CLIB_VARS:
            var_value = os.environ.get(var, None)
            if var_value is None:
                os.environ[var] = str(n_threads)


class SequentialBackend(ParallelBackendBase):
    """A ParallelBackend which will execute all batches sequentially.

    Does not use/create any threading objects, and hence has minimal
    overhead. Used when n_jobs == 1.
    """

    uses_threads = True
    supports_sharedmem = True

    def effective_n_jobs(self, n_jobs):
        """Determine the number of jobs which are going to run in parallel"""
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        return 1

    def apply_async(self, func, callback=None):
        """Schedule a func to be run"""
        result = ImmediateResult(func)
        if callback:
            callback(result)
        return result

    def get_nested_backend(self):
        nested_level = getattr(self, 'nesting_level', 0) + 1
        return SequentialBackend(nesting_level=nested_level)


class PoolManagerMixin(object):
    """A helper class for managing pool of workers."""

    _pool = None

    def effective_n_jobs(self, n_jobs):
        """Determine the number of jobs which are going to run in parallel"""
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        elif mp is None or n_jobs is None:
            # multiprocessing is not available or disabled, fallback
            # to sequential mode
            return 1
        elif n_jobs < 0:
            n_jobs = max(cpu_count() + 1 + n_jobs, 1)
        return n_jobs

    def terminate(self):
        """Shutdown the process or thread pool"""
        if self._pool is not None:
            self._pool.close()
            self._pool.terminate()  # terminate does a join()
            self._pool = None

    def _get_pool(self):
        """Used by apply_async to make it possible to implement lazy init"""
        return self._pool

    def apply_async(self, func, callback=None):
        """Schedule a func to be run"""
        return self._get_pool().apply_async(
            SafeFunction(func), callback=callback)

    def abort_everything(self, ensure_ready=True):
        """Shutdown the pool and restart a new one with the same parameters"""
        self.terminate()
        if ensure_ready:
            self.configure(n_jobs=self.parallel.n_jobs, parallel=self.parallel,
                           **self.parallel._backend_args)


class AutoBatchingMixin(object):
    """A helper class for automagically batching jobs."""

    # In seconds, should be big enough to hide multiprocessing dispatching
    # overhead.
    # This settings was found by running benchmarks/bench_auto_batching.py
    # with various parameters on various platforms.
    MIN_IDEAL_BATCH_DURATION = .2

    # Should not be too high to avoid stragglers: long jobs running alone
    # on a single worker while other workers have no work to process any more.
    MAX_IDEAL_BATCH_DURATION = 2

    # Batching counters default values
    _DEFAULT_EFFECTIVE_BATCH_SIZE = 1
    _DEFAULT_SMOOTHED_BATCH_DURATION = 0.0

    def __init__(self):
        self._effective_batch_size = self._DEFAULT_EFFECTIVE_BATCH_SIZE
        self._smoothed_batch_duration = self._DEFAULT_SMOOTHED_BATCH_DURATION

    def compute_batch_size(self):
        """Determine the optimal batch size"""
        old_batch_size = self._effective_batch_size
        batch_duration = self._smoothed_batch_duration
        if (batch_duration > 0 and
                batch_duration < self.MIN_IDEAL_BATCH_DURATION):
            # The current batch size is too small: the duration of the
            # processing of a batch of task is not large enough to hide
            # the scheduling overhead.
            ideal_batch_size = int(old_batch_size *
                                   self.MIN_IDEAL_BATCH_DURATION /
                                   batch_duration)
            # Multiply by two to limit oscilations between min and max.
            batch_size = max(2 * ideal_batch_size, 1)
            self._effective_batch_size = batch_size
            if self.parallel.verbose >= 10:
                self.parallel._print(
                    "Batch computation too fast (%.4fs.) "
                    "Setting batch_size=%d.", (batch_duration, batch_size))
        elif (batch_duration > self.MAX_IDEAL_BATCH_DURATION and
              old_batch_size >= 2):
            # The current batch size is too big. If we schedule overly long
            # running batches some CPUs might wait with nothing left to do
            # while a couple of CPUs a left processing a few long running
            # batches. Better reduce the batch size a bit to limit the
            # likelihood of scheduling such stragglers.
            batch_size = old_batch_size // 2
            self._effective_batch_size = batch_size
            if self.parallel.verbose >= 10:
                self.parallel._print(
                    "Batch computation too slow (%.4fs.) "
                    "Setting batch_size=%d.", (batch_duration, batch_size))
        else:
            # No batch size adjustment
            batch_size = old_batch_size

        if batch_size != old_batch_size:
            # Reset estimation of the smoothed mean batch duration: this
            # estimate is updated in the multiprocessing apply_async
            # CallBack as long as the batch_size is constant. Therefore
            # we need to reset the estimate whenever we re-tune the batch
            # size.
            self._smoothed_batch_duration = \
                self._DEFAULT_SMOOTHED_BATCH_DURATION

        return batch_size

    def batch_completed(self, batch_size, duration):
        """Callback indicate how long it took to run a batch"""
        if batch_size == self._effective_batch_size:
            # Update the smoothed streaming estimate of the duration of a batch
            # from dispatch to completion
            old_duration = self._smoothed_batch_duration
            if old_duration == self._DEFAULT_SMOOTHED_BATCH_DURATION:
                # First record of duration for this batch size after the last
                # reset.
                new_duration = duration
            else:
                # Update the exponentially weighted average of the duration of
                # batch for the current effective size.
                new_duration = 0.8 * old_duration + 0.2 * duration
            self._smoothed_batch_duration = new_duration

    def reset_batch_stats(self):
        """Reset batch statistics to default values.

        This avoids interferences with future jobs.
        """
        self._effective_batch_size = self._DEFAULT_EFFECTIVE_BATCH_SIZE
        self._smoothed_batch_duration = self._DEFAULT_SMOOTHED_BATCH_DURATION


class ThreadingBackend(PoolManagerMixin, ParallelBackendBase):
    """A ParallelBackend which will use a thread pool to execute batches in.

    This is a low-overhead backend but it suffers from the Python Global
    Interpreter Lock if the called function relies a lot on Python objects.
    Mostly useful when the execution bottleneck is a compiled extension that
    explicitly releases the GIL (for instance a Cython loop wrapped in a "with
    nogil" block or an expensive call to a library such as NumPy).

    The actual thread pool is lazily initialized: the actual thread pool
    construction is delayed to the first call to apply_async.

    ThreadingBackend is used as the default backend for nested calls.
    """

    supports_timeout = True
    uses_threads = True
    supports_sharedmem = True

    def configure(self, n_jobs=1, parallel=None, **backend_args):
        """Build a process or thread pool and return the number of workers"""
        n_jobs = self.effective_n_jobs(n_jobs)
        if n_jobs == 1:
            # Avoid unnecessary overhead and use sequential backend instead.
            raise FallbackToBackend(
                SequentialBackend(nesting_level=self.nesting_level))
        self.parallel = parallel
        self._n_jobs = n_jobs
        return n_jobs

    def _get_pool(self):
        """Lazily initialize the thread pool

        The actual pool of worker threads is only initialized at the first
        call to apply_async.
        """
        if self._pool is None:
            self._pool = ThreadPool(self._n_jobs)
        return self._pool


class MultiprocessingBackend(PoolManagerMixin, AutoBatchingMixin,
                             ParallelBackendBase):
    """A ParallelBackend which will use a multiprocessing.Pool.

    Will introduce some communication and memory overhead when exchanging
    input and output data with the with the worker Python processes.
    However, does not suffer from the Python Global Interpreter Lock.
    """

    # Environment variables to protect against bad situations when nesting
    JOBLIB_SPAWNED_PROCESS = "__JOBLIB_SPAWNED_PARALLEL__"

    supports_timeout = True

    def effective_n_jobs(self, n_jobs):
        """Determine the number of jobs which are going to run in parallel.

        This also checks if we are attempting to create a nested parallel
        loop.
        """
        if mp is None:
            return 1

        if mp.current_process().daemon:
            # Daemonic processes cannot have children
            if n_jobs != 1:
                warnings.warn(
                    'Multiprocessing-backed parallel loops cannot be nested,'
                    ' setting n_jobs=1',
                    stacklevel=3)
            return 1

        if process_executor._CURRENT_DEPTH > 0:
            # Mixing loky and multiprocessing in nested loop is not supported
            if n_jobs != 1:
                warnings.warn(
                    'Multiprocessing-backed parallel loops cannot be nested,'
                    ' below loky, setting n_jobs=1',
                    stacklevel=3)
            return 1

        if not isinstance(threading.current_thread(), threading._MainThread):
            # Prevent posix fork inside in non-main posix threads
            if n_jobs != 1:
                warnings.warn(
                    'Multiprocessing-backed parallel loops cannot be nested'
                    ' below threads, setting n_jobs=1',
                    stacklevel=3)
            return 1

        return super(MultiprocessingBackend, self).effective_n_jobs(n_jobs)

    def configure(self, n_jobs=1, parallel=None, prefer=None, require=None,
                  **memmappingpool_args):
        """Build a process or thread pool and return the number of workers"""
        n_jobs = self.effective_n_jobs(n_jobs)
        if n_jobs == 1:
            raise FallbackToBackend(
                SequentialBackend(nesting_level=self.nesting_level))

        already_forked = int(os.environ.get(self.JOBLIB_SPAWNED_PROCESS, 0))
        if already_forked:
            raise ImportError(
                '[joblib] Attempting to do parallel computing '
                'without protecting your import on a system that does '
                'not support forking. To use parallel-computing in a '
                'script, you must protect your main loop using "if '
                "__name__ == '__main__'"
                '". Please see the joblib documentation on Parallel '
                'for more information')
        # Set an environment variable to avoid infinite loops
        os.environ[self.JOBLIB_SPAWNED_PROCESS] = '1'

        # Make sure to free as much memory as possible before forking
        gc.collect()
        self._pool = MemmappingPool(
            n_jobs, initializer=self.limit_clib_threads, **memmappingpool_args)
        self.parallel = parallel
        return n_jobs

    def terminate(self):
        """Shutdown the process or thread pool"""
        super(MultiprocessingBackend, self).terminate()
        if self.JOBLIB_SPAWNED_PROCESS in os.environ:
            del os.environ[self.JOBLIB_SPAWNED_PROCESS]

        self.reset_batch_stats()


class LokyBackend(AutoBatchingMixin, ParallelBackendBase):
    """Managing pool of workers with loky instead of multiprocessing."""

    supports_timeout = True

    def configure(self, n_jobs=1, parallel=None, prefer=None, require=None,
                  idle_worker_timeout=300, **memmappingexecutor_args):
        """Build a process executor and return the number of workers"""
        n_jobs = self.effective_n_jobs(n_jobs)
        if n_jobs == 1:
            raise FallbackToBackend(
                SequentialBackend(nesting_level=self.nesting_level))

        self._workers = get_memmapping_executor(
            n_jobs, timeout=idle_worker_timeout,
            initializer=self.limit_clib_threads,
            **memmappingexecutor_args)
        self.parallel = parallel
        return n_jobs

    def effective_n_jobs(self, n_jobs):
        """Determine the number of jobs which are going to run in parallel"""
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        elif mp is None or n_jobs is None:
            # multiprocessing is not available or disabled, fallback
            # to sequential mode
            return 1
        elif mp.current_process().daemon:
            # Daemonic processes cannot have children
            if n_jobs != 1:
                warnings.warn(
                    'Loky-backed parallel loops cannot be called in a'
                    ' multiprocessing, setting n_jobs=1',
                    stacklevel=3)
            return 1
        elif not isinstance(threading.current_thread(), threading._MainThread):
            # Prevent posix fork inside in non-main posix threads
            if n_jobs != 1:
                warnings.warn(
                    'Loky-backed parallel loops cannot be nested below '
                    'threads, setting n_jobs=1',
                    stacklevel=3)
            return 1
        elif n_jobs < 0:
            n_jobs = max(cpu_count() + 1 + n_jobs, 1)
        return n_jobs

    def apply_async(self, func, callback=None):
        """Schedule a func to be run"""
        future = self._workers.submit(SafeFunction(func))
        future.get = functools.partial(self.wrap_future_result, future)
        if callback is not None:
            future.add_done_callback(callback)
        return future

    @staticmethod
    def wrap_future_result(future, timeout=None):
        """Wrapper for Future.result to implement the same behaviour as
        AsyncResults.get from multiprocessing."""
        try:
            return future.result(timeout=timeout)
        except LokyTimeoutError:
            raise TimeoutError()

    def terminate(self):
        if self._workers is not None:
            # Terminate does not shutdown the workers as we want to reuse them
            # in latter calls but we free as much memory as we can by deleting
            # the shared memory
            delete_folder(self._workers._temp_folder)
            self._workers = None

        self.reset_batch_stats()

    def abort_everything(self, ensure_ready=True):
        """Shutdown the workers and restart a new one with the same parameters
        """
        self._workers.shutdown(kill_workers=True)
        delete_folder(self._workers._temp_folder)
        self._workers = None
        if ensure_ready:
            self.configure(n_jobs=self.parallel.n_jobs, parallel=self.parallel)


class ImmediateResult(object):
    def __init__(self, batch):
        # Don't delay the application, to avoid keeping the input
        # arguments in memory
        self.results = batch()

    def get(self):
        return self.results


class SafeFunction(object):
    """Wrapper that handles the serialization of exception tracebacks.

    If an exception is triggered when calling the inner function, a copy of
    the full traceback is captured to make it possible to serialize
    it so that it can be rendered in a different Python process.
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        try:
            return self.func(*args, **kwargs)
        except KeyboardInterrupt:
            # We capture the KeyboardInterrupt and reraise it as
            # something different, as multiprocessing does not
            # interrupt processing for a KeyboardInterrupt
            raise WorkerInterrupt()
        except BaseException:
            if PY27:
                # Capture the traceback of the worker to make it part of
                # the final exception message.
                e_type, e_value, e_tb = sys.exc_info()
                text = format_exc(e_type, e_value, e_tb, context=10,
                                  tb_offset=1)
                raise TransportableException(text, e_type)
            else:
                # Rely on Python 3 built-in Remote Traceback reporting
                raise


class FallbackToBackend(Exception):
    """Raised when configuration should fallback to another backend"""

    def __init__(self, backend):
        self.backend = backend
