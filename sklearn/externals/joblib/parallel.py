"""
Helpers for embarrassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

from __future__ import division

import os
import sys
import gc
import warnings
from math import sqrt
import functools
import time
import threading
import itertools
from numbers import Integral
try:
    import cPickle as pickle
except:
    import pickle

from ._multiprocessing_helpers import mp
if mp is not None:
    from .pool import MemmapingPool
    from multiprocessing.pool import ThreadPool

from .format_stack import format_exc, format_outer_frames
from .logger import Logger, short_format_time
from .my_exceptions import TransportableException, _mk_exception
from .disk import memstr_to_kbytes
from ._compat import _basestring


VALID_BACKENDS = ['multiprocessing', 'threading']

# Environment variables to protect against bad situations when nesting
JOBLIB_SPAWNED_PROCESS = "__JOBLIB_SPAWNED_PARALLEL__"

# In seconds, should be big enough to hide multiprocessing dispatching
# overhead.
# This settings was found by running benchmarks/bench_auto_batching.py
# with various parameters on various platforms.
MIN_IDEAL_BATCH_DURATION = .2

# Should not be too high to avoid stragglers: long jobs running alone
# on a single worker while other workers have no work to process any more.
MAX_IDEAL_BATCH_DURATION = 2

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


class BatchedCalls(object):
    """Wrap a sequence of (func, args, kwargs) tuples as a single callable"""

    def __init__(self, iterator_slice):
        self.items = list(iterator_slice)
        self._size = len(self.items)

    def __call__(self):
        return [func(*args, **kwargs) for func, args, kwargs in self.items]

    def __len__(self):
        return self._size


###############################################################################
# CPU count that works also when multiprocessing has been disabled via
# the JOBLIB_MULTIPROCESSING environment variable
def cpu_count():
    """ Return the number of CPUs.
    """
    if mp is None:
        return 1
    return mp.cpu_count()


###############################################################################
# For verbosity

def _verbosity_filter(index, verbose):
    """ Returns False for indices increasingly apart, the distance
        depending on the value of verbose.

        We use a lag increasing as the square of index
    """
    if not verbose:
        return True
    elif verbose > 10:
        return False
    if index == 0:
        return False
    verbose = .5 * (11 - verbose) ** 2
    scale = sqrt(index / verbose)
    next_scale = sqrt((index + 1) / verbose)
    return (int(next_scale) == int(scale))


###############################################################################
class WorkerInterrupt(Exception):
    """ An exception that is not KeyboardInterrupt to allow subprocesses
        to be interrupted.
    """
    pass


###############################################################################
class SafeFunction(object):
    """ Wraps a function to make it exception with full traceback in
        their representation.
        Useful for parallel computing with multiprocessing, for which
        exceptions cannot be captured.
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
        except:
            e_type, e_value, e_tb = sys.exc_info()
            text = format_exc(e_type, e_value, e_tb, context=10,
                              tb_offset=1)
            raise TransportableException(text, e_type)


###############################################################################
def delayed(function, check_pickle=True):
    """Decorator used to capture the arguments of a function.

    Pass `check_pickle=False` when:

    - performing a possibly repeated check is too costly and has been done
      already once outside of the call to delayed.

    - when used in conjunction `Parallel(backend='threading')`.

    """
    # Try to pickle the input function, to catch the problems early when
    # using with multiprocessing:
    if check_pickle:
        pickle.dumps(function)

    def delayed_function(*args, **kwargs):
        return function, args, kwargs
    try:
        delayed_function = functools.wraps(function)(delayed_function)
    except AttributeError:
        " functools.wraps fails on some callable objects "
    return delayed_function


###############################################################################
class ImmediateComputeBatch(object):
    """Sequential computation of a batch of tasks.

    This replicates the async computation API but actually does not delay
    the computations when joblib.Parallel runs in sequential mode.

    """
    def __init__(self, batch):
        # Don't delay the application, to avoid keeping the input
        # arguments in memory
        self.results = batch()

    def get(self):
        return self.results


###############################################################################
class BatchCompletionCallBack(object):
    """Callback used by joblib.Parallel's multiprocessing backend.

    This callable is executed by the parent process whenever a worker process
    has returned the results of a batch of tasks.

    It is used for progress reporting, to update estimate of the batch
    processing duration and to schedule the next batch of tasks to be
    processed.

    """
    def __init__(self, dispatch_timestamp, batch_size, parallel):
        self.dispatch_timestamp = dispatch_timestamp
        self.batch_size = batch_size
        self.parallel = parallel

    def __call__(self, out):
        self.parallel.n_completed_tasks += self.batch_size
        this_batch_duration = time.time() - self.dispatch_timestamp

        if (self.parallel.batch_size == 'auto'
                and self.batch_size == self.parallel._effective_batch_size):
            # Update the smoothed streaming estimate of the duration of a batch
            # from dispatch to completion
            old_duration = self.parallel._smoothed_batch_duration
            if old_duration == 0:
                # First record of duration for this batch size after the last
                # reset.
                new_duration = this_batch_duration
            else:
                # Update the exponentially weighted average of the duration of
                # batch for the current effective size.
                new_duration = 0.8 * old_duration + 0.2 * this_batch_duration
            self.parallel._smoothed_batch_duration = new_duration

        self.parallel.print_progress()
        if self.parallel._original_iterator is not None:
            self.parallel.dispatch_next()


###############################################################################
class Parallel(Logger):
    ''' Helper class for readable parallel mapping.

        Parameters
        -----------
        n_jobs: int, default: 1
            The maximum number of concurrently running jobs, such as the number
            of Python worker processes when backend="multiprocessing"
            or the size of the thread-pool when backend="threading".
            If -1 all CPUs are used. If 1 is given, no parallel computing code
            is used at all, which is useful for debugging. For n_jobs below -1,
            (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all
            CPUs but one are used.
        backend: str or None, default: 'multiprocessing'
            Specify the parallelization backend implementation.
            Supported backends are:
              - "multiprocessing" used by default, can induce some
                communication and memory overhead when exchanging input and
                output data with the with the worker Python processes.
              - "threading" is a very low-overhead backend but it suffers
                from the Python Global Interpreter Lock if the called function
                relies a lot on Python objects. "threading" is mostly useful
                when the execution bottleneck is a compiled extension that
                explicitly releases the GIL (for instance a Cython loop wrapped
                in a "with nogil" block or an expensive call to a library such
                as NumPy).
        verbose: int, optional
            The verbosity level: if non zero, progress messages are
            printed. Above 50, the output is sent to stdout.
            The frequency of the messages increases with the verbosity level.
            If it more than 10, all iterations are reported.
        pre_dispatch: {'all', integer, or expression, as in '3*n_jobs'}
            The number of batches (of tasks) to be pre-dispatched.
            Default is '2*n_jobs'. When batch_size="auto" this is reasonable
            default and the multiprocessing workers shoud never starve.
        batch_size: int or 'auto', default: 'auto'
            The number of atomic tasks to dispatch at once to each
            worker. When individual evaluations are very fast, multiprocessing
            can be slower than sequential computation because of the overhead.
            Batching fast computations together can mitigate this.
            The ``'auto'`` strategy keeps track of the time it takes for a batch
            to complete, and dynamically adjusts the batch size to keep the time
            on the order of half a second, using a heuristic. The initial batch
            size is 1.
            ``batch_size="auto"`` with ``backend="threading"`` will dispatch
            batches of a single task at a time as the threading backend has
            very little overhead and using larger batch size has not proved to
            bring any gain in that case.
        temp_folder: str, optional
            Folder to be used by the pool for memmaping large arrays
            for sharing memory with worker processes. If None, this will try in
            order:
            - a folder pointed by the JOBLIB_TEMP_FOLDER environment variable,
            - /dev/shm if the folder exists and is writable: this is a RAMdisk
              filesystem available by default on modern Linux distributions,
            - the default system temporary folder that can be overridden
              with TMP, TMPDIR or TEMP environment variables, typically /tmp
              under Unix operating systems.
            Only active when backend="multiprocessing".
        max_nbytes int, str, or None, optional, 1M by default
            Threshold on the size of arrays passed to the workers that
            triggers automated memory mapping in temp_folder. Can be an int
            in Bytes, or a human-readable string, e.g., '1M' for 1 megabyte.
            Use None to disable memmaping of large arrays.
            Only active when backend="multiprocessing".

        Notes
        -----

        This object uses the multiprocessing module to compute in
        parallel the application of a function to many different
        arguments. The main functionality it brings in addition to
        using the raw multiprocessing API are (see examples for details):

            * More readable code, in particular since it avoids
              constructing list of arguments.

            * Easier debugging:
                - informative tracebacks even when the error happens on
                  the client side
                - using 'n_jobs=1' enables to turn off parallel computing
                  for debugging without changing the codepath
                - early capture of pickling errors

            * An optional progress meter.

            * Interruption of multiprocesses jobs with 'Ctrl-C'

            * Flexible pickling control for the communication to and from
              the worker processes.

            * Ability to use shared memory efficiently with worker
              processes for large numpy-based datastructures.

        Examples
        --------

        A simple example:

        >>> from math import sqrt
        >>> from sklearn.externals.joblib import Parallel, delayed
        >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

        Reshaping the output when the function has several return
        values:

        >>> from math import modf
        >>> from sklearn.externals.joblib import Parallel, delayed
        >>> r = Parallel(n_jobs=1)(delayed(modf)(i/2.) for i in range(10))
        >>> res, i = zip(*r)
        >>> res
        (0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5)
        >>> i
        (0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0)

        The progress meter: the higher the value of `verbose`, the more
        messages::

            >>> from time import sleep
            >>> from sklearn.externals.joblib import Parallel, delayed
            >>> r = Parallel(n_jobs=2, verbose=5)(delayed(sleep)(.1) for _ in range(10)) #doctest: +SKIP
            [Parallel(n_jobs=2)]: Done   1 out of  10 | elapsed:    0.1s remaining:    0.9s
            [Parallel(n_jobs=2)]: Done   3 out of  10 | elapsed:    0.2s remaining:    0.5s
            [Parallel(n_jobs=2)]: Done   6 out of  10 | elapsed:    0.3s remaining:    0.2s
            [Parallel(n_jobs=2)]: Done   9 out of  10 | elapsed:    0.5s remaining:    0.1s
            [Parallel(n_jobs=2)]: Done  10 out of  10 | elapsed:    0.5s finished

        Traceback example, note how the line of the error is indicated
        as well as the values of the parameter passed to the function that
        triggered the exception, even though the traceback happens in the
        child process::

         >>> from heapq import nlargest
         >>> from sklearn.externals.joblib import Parallel, delayed
         >>> Parallel(n_jobs=2)(delayed(nlargest)(2, n) for n in (range(4), 'abcde', 3)) #doctest: +SKIP
         #...
         ---------------------------------------------------------------------------
         Sub-process traceback:
         ---------------------------------------------------------------------------
         TypeError                                          Mon Nov 12 11:37:46 2012
         PID: 12934                                    Python 2.7.3: /usr/bin/python
         ...........................................................................
         /usr/lib/python2.7/heapq.pyc in nlargest(n=2, iterable=3, key=None)
             419         if n >= size:
             420             return sorted(iterable, key=key, reverse=True)[:n]
             421
             422     # When key is none, use simpler decoration
             423     if key is None:
         --> 424         it = izip(iterable, count(0,-1))                    # decorate
             425         result = _nlargest(n, it)
             426         return map(itemgetter(0), result)                   # undecorate
             427
             428     # General case, slowest method

         TypeError: izip argument #1 must support iteration
         ___________________________________________________________________________


        Using pre_dispatch in a producer/consumer situation, where the
        data is generated on the fly. Note how the producer is first
        called a 3 times before the parallel loop is initiated, and then
        called to generate new data on the fly. In this case the total
        number of iterations cannot be reported in the progress messages::

         >>> from math import sqrt
         >>> from sklearn.externals.joblib import Parallel, delayed

         >>> def producer():
         ...     for i in range(6):
         ...         print('Produced %s' % i)
         ...         yield i

         >>> out = Parallel(n_jobs=2, verbose=100, pre_dispatch='1.5*n_jobs')(
         ...                         delayed(sqrt)(i) for i in producer()) #doctest: +SKIP
         Produced 0
         Produced 1
         Produced 2
         [Parallel(n_jobs=2)]: Done   1 jobs       | elapsed:    0.0s
         Produced 3
         [Parallel(n_jobs=2)]: Done   2 jobs       | elapsed:    0.0s
         Produced 4
         [Parallel(n_jobs=2)]: Done   3 jobs       | elapsed:    0.0s
         Produced 5
         [Parallel(n_jobs=2)]: Done   4 jobs       | elapsed:    0.0s
         [Parallel(n_jobs=2)]: Done   5 out of   6 | elapsed:    0.0s remaining:    0.0s
         [Parallel(n_jobs=2)]: Done   6 out of   6 | elapsed:    0.0s finished
    '''
    def __init__(self, n_jobs=1, backend='multiprocessing', verbose=0,
                 pre_dispatch='2 * n_jobs', batch_size='auto',
                 temp_folder=None, max_nbytes='1M', mmap_mode='r'):
        self.verbose = verbose
        self._mp_context = DEFAULT_MP_CONTEXT
        if backend is None:
            # `backend=None` was supported in 0.8.2 with this effect
            backend = "multiprocessing"
        elif hasattr(backend, 'Pool') and hasattr(backend, 'Lock'):
            # Make it possible to pass a custom multiprocessing context as
            # backend to change the start method to forkserver or spawn or
            # preload modules on the forkserver helper process.
            self._mp_context = backend
            backend = "multiprocessing"
        if backend not in VALID_BACKENDS:
            raise ValueError("Invalid backend: %s, expected one of %r"
                             % (backend, VALID_BACKENDS))
        self.backend = backend
        self.n_jobs = n_jobs
        if (batch_size == 'auto'
                or isinstance(batch_size, Integral) and batch_size > 0):
            self.batch_size = batch_size
        else:
            raise ValueError(
                "batch_size must be 'auto' or a positive integer, got: %r"
                % batch_size)

        self.pre_dispatch = pre_dispatch
        self._temp_folder = temp_folder
        if isinstance(max_nbytes, _basestring):
            self._max_nbytes = 1024 * memstr_to_kbytes(max_nbytes)
        else:
            self._max_nbytes = max_nbytes
        self._mmap_mode = mmap_mode
        # Not starting the pool in the __init__ is a design decision, to be
        # able to close it ASAP, and not burden the user with closing it
        # unless they choose to use the context manager API with a with block.
        self._pool = None
        self._output = None
        self._jobs = list()
        self._managed_pool = False

        # This lock is used coordinate the main thread of this process with
        # the async callback thread of our the pool.
        self._lock = threading.Lock()

    def __enter__(self):
        self._managed_pool = True
        self._initialize_pool()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._terminate_pool()
        self._managed_pool = False

    def _effective_n_jobs(self):
        n_jobs = self.n_jobs
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        elif mp is None or n_jobs is None:
            # multiprocessing is not available or disabled, fallback
            # to sequential mode
            return 1
        elif n_jobs < 0:
            n_jobs = max(mp.cpu_count() + 1 + n_jobs, 1)
        return n_jobs

    def _initialize_pool(self):
        """Build a process or thread pool and return the number of workers"""
        n_jobs = self._effective_n_jobs()
        # The list of exceptions that we will capture
        self.exceptions = [TransportableException]

        if n_jobs == 1:
            # Sequential mode: do not use a pool instance to avoid any
            # useless dispatching overhead
            self._pool = None
        elif self.backend == 'threading':
            self._pool = ThreadPool(n_jobs)
        elif self.backend == 'multiprocessing':
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
            else:
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
                poolargs = dict(
                    max_nbytes=self._max_nbytes,
                    mmap_mode=self._mmap_mode,
                    temp_folder=self._temp_folder,
                    verbose=max(0, self.verbose - 50),
                )
                if self._mp_context is not None:
                    # Use Python 3.4+ multiprocessing context isolation
                    poolargs['context'] = self._mp_context
                self._pool = MemmapingPool(n_jobs, **poolargs)

                # We are using multiprocessing, we also want to capture
                # KeyboardInterrupts
                self.exceptions.extend([KeyboardInterrupt, WorkerInterrupt])
        else:
            raise ValueError("Unsupported backend: %s" % self.backend)
        return n_jobs

    def _terminate_pool(self):
        if self._pool is not None:
            self._pool.close()
            self._pool.terminate()  # terminate does a join()
            self._pool = None
            if self.backend == 'multiprocessing':
                os.environ.pop(JOBLIB_SPAWNED_PROCESS, 0)

    def _dispatch(self, batch):
        """Queue the batch for computing, with or without multiprocessing

        WARNING: this method is not thread-safe: it should be only called
        indirectly via dispatch_one_batch.

        """
        # If job.get() catches an exception, it closes the queue:
        if self._aborting:
            return

        if self._pool is None:
            job = ImmediateComputeBatch(batch)
            self._jobs.append(job)
            self.n_dispatched_batches += 1
            self.n_dispatched_tasks += len(batch)
            self.n_completed_tasks += len(batch)
            if not _verbosity_filter(self.n_dispatched_batches, self.verbose):
                self._print('Done %3i tasks       | elapsed: %s',
                        (self.n_completed_tasks,
                            short_format_time(time.time() - self._start_time)
                        ))
        else:
            dispatch_timestamp = time.time()
            cb = BatchCompletionCallBack(dispatch_timestamp, len(batch), self)
            job = self._pool.apply_async(SafeFunction(batch), callback=cb)
            self._jobs.append(job)
            self.n_dispatched_tasks += len(batch)
            self.n_dispatched_batches += 1

    def dispatch_next(self):
        """Dispatch more data for parallel processing

        This method is meant to be called concurrently by the multiprocessing
        callback. We rely on the thread-safety of dispatch_one_batch to protect
        against concurrent consumption of the unprotected iterator.

        """
        if not self.dispatch_one_batch(self._original_iterator):
            self._iterating = False
            self._original_iterator = None

    def dispatch_one_batch(self, iterator):
        """Prefetch the tasks for the next batch and dispatch them.

        The effective size of the batch is computed here.
        If there are no more jobs to dispatch, return False, else return True.

        The iterator consumption and dispatching is protected by the same
        lock so calling this function should be thread safe.

        """
        if self.batch_size == 'auto' and self.backend == 'threading':
            # Batching is never beneficial with the threading backend
            batch_size = 1
        elif self.batch_size == 'auto':
            old_batch_size = self._effective_batch_size
            batch_duration = self._smoothed_batch_duration
            if (batch_duration > 0 and
                    batch_duration < MIN_IDEAL_BATCH_DURATION):
                # The current batch size is too small: the duration of the
                # processing of a batch of task is not large enough to hide
                # the scheduling overhead.
                ideal_batch_size = int(
                    old_batch_size * MIN_IDEAL_BATCH_DURATION / batch_duration)
                # Multiply by two to limit oscilations between min and max.
                batch_size = max(2 * ideal_batch_size, 1)
                self._effective_batch_size = batch_size
                if self.verbose >= 10:
                    self._print("Batch computation too fast (%.4fs.) "
                                "Setting batch_size=%d.", (
                                    batch_duration, batch_size))
            elif (batch_duration > MAX_IDEAL_BATCH_DURATION and
                  old_batch_size >= 2):
                # The current batch size is too big. If we schedule overly long
                # running batches some CPUs might wait with nothing left to do
                # while a couple of CPUs a left processing a few long running
                # batches. Better reduce the batch size a bit to limit the
                # likelihood of scheduling such stragglers.
                self._effective_batch_size = batch_size = old_batch_size // 2
                if self.verbose >= 10:
                    self._print("Batch computation too slow (%.2fs.) "
                                "Setting batch_size=%d.", (
                                    batch_duration, batch_size))
            else:
                # No batch size adjustment
                batch_size = old_batch_size

            if batch_size != old_batch_size:
                # Reset estimation of the smoothed mean batch duration: this
                # estimate is updated in the multiprocessing apply_async
                # CallBack as long as the batch_size is constant. Therefore
                # we need to reset the estimate whenever we re-tune the batch
                # size.
                self._smoothed_batch_duration = 0
        else:
            # Fixed batch size strategy
            batch_size = self.batch_size
        with self._lock:
            tasks = BatchedCalls(itertools.islice(iterator, batch_size))
            if not tasks:
                # No more tasks available in the iterator: tell caller to stop.
                return False
            else:
                self._dispatch(tasks)
                return True

    def _print(self, msg, msg_args):
        """Display the message on stout or stderr depending on verbosity"""
        # XXX: Not using the logger framework: need to
        # learn to use logger better.
        if not self.verbose:
            return
        if self.verbose < 50:
            writer = sys.stderr.write
        else:
            writer = sys.stdout.write
        msg = msg % msg_args
        writer('[%s]: %s\n' % (self, msg))

    def print_progress(self):
        """Display the process of the parallel execution only a fraction
           of time, controlled by self.verbose.
        """
        if not self.verbose:
            return
        elapsed_time = time.time() - self._start_time

        # This is heuristic code to print only 'verbose' times a messages
        # The challenge is that we may not know the queue length
        if self._original_iterator:
            if _verbosity_filter(self.n_dispatched_batches, self.verbose):
                return
            self._print('Done %3i tasks      | elapsed: %s',
                        (self.n_completed_tasks,
                         short_format_time(elapsed_time),
                        ))
        else:
            index = self.n_dispatched_batches
            # We are finished dispatching
            total_tasks = self.n_dispatched_tasks
            # We always display the first loop
            if not index == 0:
                # Display depending on the number of remaining items
                # A message as soon as we finish dispatching, cursor is 0
                cursor = (total_tasks - index + 1
                          - self._pre_dispatch_amount)
                frequency = (total_tasks // self.verbose) + 1
                is_last_item = (index + 1 == total_tasks)
                if (is_last_item or cursor % frequency):
                    return
            remaining_time = (elapsed_time / (index + 1) *
                              (self.n_dispatched_tasks - index - 1.))
            self._print('Done %3i out of %3i | elapsed: %s remaining: %s',
                        (index + 1,
                         total_tasks,
                         short_format_time(elapsed_time),
                         short_format_time(remaining_time),
                        ))

    def retrieve(self):
        self._output = list()
        while self._iterating or len(self._jobs) > 0:
            if len(self._jobs) == 0:
                # Wait for an async callback to dispatch new jobs
                time.sleep(0.01)
                continue
            # We need to be careful: the job list can be filling up as
            # we empty it and Python list are not thread-safe by default hence
            # the use of the lock
            with self._lock:
                job = self._jobs.pop(0)
            try:
                self._output.extend(job.get())
            except tuple(self.exceptions) as exception:
                # Stop dispatching any new job in the async callback thread
                self._aborting = True

                if isinstance(exception, TransportableException):
                    # Capture exception to add information on the local
                    # stack in addition to the distant stack
                    this_report = format_outer_frames(context=10,
                                                      stack_start=1)
                    report = """Multiprocessing exception:
%s
---------------------------------------------------------------------------
Sub-process traceback:
---------------------------------------------------------------------------
%s""" % (this_report, exception.message)
                    # Convert this to a JoblibException
                    exception_type = _mk_exception(exception.etype)[0]
                    exception = exception_type(report)

                # Kill remaining running processes without waiting for
                # the results as we will raise the exception we got back
                # to the caller instead of returning any result.
                self._terminate_pool()
                if self._managed_pool:
                    # In case we had to terminate a managed pool, let
                    # us start a new one to ensure that subsequent calls
                    # to __call__ on the same Parallel instance will get
                    # a working pool as they expect.
                    self._initialize_pool()
                raise exception

    def __call__(self, iterable):
        if self._jobs:
            raise ValueError('This Parallel instance is already running')
        # A flag used to abort the dispatching of jobs in case an
        # exception is found
        self._aborting = False
        if not self._managed_pool:
            n_jobs = self._initialize_pool()
        else:
            n_jobs = self._effective_n_jobs()

        if self.batch_size == 'auto':
            self._effective_batch_size = 1

        iterator = iter(iterable)
        pre_dispatch = self.pre_dispatch

        if pre_dispatch == 'all' or n_jobs == 1:
            # prevent further dispatch via multiprocessing callback thread
            self._original_iterator = None
            self._pre_dispatch_amount = 0
        else:
            self._original_iterator = iterator
            if hasattr(pre_dispatch, 'endswith'):
                pre_dispatch = eval(pre_dispatch)
            self._pre_dispatch_amount = pre_dispatch = int(pre_dispatch)

            # The main thread will consume the first pre_dispatch items and
            # the remaining items will later be lazily dispatched by async
            # callbacks upon task completions.
            iterator = itertools.islice(iterator, pre_dispatch)

        self._start_time = time.time()
        self.n_dispatched_batches = 0
        self.n_dispatched_tasks = 0
        self.n_completed_tasks = 0
        self._smoothed_batch_duration = 0.0
        try:
            # Only set self._iterating to True if at least a batch
            # was dispatched. In particular this covers the edge
            # case of Parallel used with an exhausted iterator.
            while self.dispatch_one_batch(iterator):
                self._iterating = True
            else:
                self._iterating = False

            if pre_dispatch == "all" or n_jobs == 1:
                # The iterable was consumed all at once by the above for loop.
                # No need to wait for async callbacks to trigger to
                # consumption.
                self._iterating = False
            self.retrieve()
            # Make sure that we get a last message telling us we are done
            elapsed_time = time.time() - self._start_time
            self._print('Done %3i out of %3i | elapsed: %s finished',
                        (len(self._output), len(self._output),
                         short_format_time(elapsed_time)))
        finally:
            if not self._managed_pool:
                self._terminate_pool()
            self._jobs = list()
        output = self._output
        self._output = None
        return output

    def __repr__(self):
        return '%s(n_jobs=%s)' % (self.__class__.__name__, self.n_jobs)
