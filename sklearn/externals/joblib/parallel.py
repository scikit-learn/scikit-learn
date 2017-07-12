"""
Helpers for embarrassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

from __future__ import division

import os
import sys
from math import sqrt
import functools
import time
import threading
import itertools
from numbers import Integral
from contextlib import contextmanager
import warnings
try:
    import cPickle as pickle
except ImportError:
    import pickle

from ._multiprocessing_helpers import mp

from .format_stack import format_outer_frames
from .logger import Logger, short_format_time
from .my_exceptions import TransportableException, _mk_exception
from .disk import memstr_to_bytes
from ._parallel_backends import (FallbackToBackend, MultiprocessingBackend,
                                 ThreadingBackend, SequentialBackend)
from ._compat import _basestring

# Make sure that those two classes are part of the public joblib.parallel API
# so that 3rd party backend implementers can import them from here.
from ._parallel_backends import AutoBatchingMixin  # noqa
from ._parallel_backends import ParallelBackendBase  # noqa

BACKENDS = {
    'multiprocessing': MultiprocessingBackend,
    'threading': ThreadingBackend,
    'sequential': SequentialBackend,
}

# name of the backend used by default by Parallel outside of any context
# managed by ``parallel_backend``.
DEFAULT_BACKEND = 'multiprocessing'
DEFAULT_N_JOBS = 1

# Thread local value that can be overridden by the ``parallel_backend`` context
# manager
_backend = threading.local()


def get_active_backend():
    """Return the active default backend"""
    active_backend_and_jobs = getattr(_backend, 'backend_and_jobs', None)
    if active_backend_and_jobs is not None:
        return active_backend_and_jobs
    # We are outside of the scope of any parallel_backend context manager,
    # create the default backend instance now
    active_backend = BACKENDS[DEFAULT_BACKEND]()
    return active_backend, DEFAULT_N_JOBS


@contextmanager
def parallel_backend(backend, n_jobs=-1, **backend_params):
    """Change the default backend used by Parallel inside a with block.

    If ``backend`` is a string it must match a previously registered
    implementation using the ``register_parallel_backend`` function.

    Alternatively backend can be passed directly as an instance.

    By default all available workers will be used (``n_jobs=-1``) unless the
    caller passes an explicit value for the ``n_jobs`` parameter.

    This is an alternative to passing a ``backend='backend_name'`` argument to
    the ``Parallel`` class constructor. It is particularly useful when calling
    into library code that uses joblib internally but does not expose the
    backend argument in its own API.

    >>> from operator import neg
    >>> with parallel_backend('threading'):
    ...     print(Parallel()(delayed(neg)(i + 1) for i in range(5)))
    ...
    [-1, -2, -3, -4, -5]

    Warning: this function is experimental and subject to change in a future
    version of joblib.

    .. versionadded:: 0.10

    """
    if isinstance(backend, _basestring):
        backend = BACKENDS[backend](**backend_params)
    old_backend_and_jobs = getattr(_backend, 'backend_and_jobs', None)
    try:
        _backend.backend_and_jobs = (backend, n_jobs)
        # return the backend instance to make it easier to write tests
        yield backend, n_jobs
    finally:
        if old_backend_and_jobs is None:
            if getattr(_backend, 'backend_and_jobs', None) is not None:
                del _backend.backend_and_jobs
        else:
            _backend.backend_and_jobs = old_backend_and_jobs


# Under Linux or OS X the default start method of multiprocessing
# can cause third party libraries to crash. Under Python 3.4+ it is possible
# to set an environment variable to switch the default start method from
# 'fork' to 'forkserver' or 'spawn' to avoid this issue albeit at the cost
# of causing semantic changes and some additional pool instantiation overhead.
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
    """Return the number of CPUs."""
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

        self.parallel._backend.batch_completed(self.batch_size,
                                               this_batch_duration)
        self.parallel.print_progress()
        if self.parallel._original_iterator is not None:
            self.parallel.dispatch_next()


###############################################################################
def register_parallel_backend(name, factory, make_default=False):
    """Register a new Parallel backend factory.

    The new backend can then be selected by passing its name as the backend
    argument to the Parallel class. Moreover, the default backend can be
    overwritten globally by setting make_default=True.

    The factory can be any callable that takes no argument and return an
    instance of ``ParallelBackendBase``.

    Warning: this function is experimental and subject to change in a future
    version of joblib.

    .. versionadded:: 0.10

    """
    BACKENDS[name] = factory
    if make_default:
        global DEFAULT_BACKEND
        DEFAULT_BACKEND = name


def effective_n_jobs(n_jobs=-1):
    """Determine the number of jobs that can actually run in parallel

    n_jobs is the is the number of workers requested by the callers.
    Passing n_jobs=-1 means requesting all available workers for instance
    matching the number of CPU cores on the worker host(s).

    This method should return a guesstimate of the number of workers that can
    actually perform work concurrently with the currently enabled default
    backend. The primary use case is to make it possible for the caller to know
    in how many chunks to slice the work.

    In general working on larger data chunks is more efficient (less
    scheduling overhead and better use of CPU cache prefetching heuristics)
    as long as all the workers have enough work to do.

    Warning: this function is experimental and subject to change in a future
    version of joblib.

    .. versionadded:: 0.10

    """
    backend, _ = get_active_backend()
    return backend.effective_n_jobs(n_jobs=n_jobs)


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
        backend: str, ParallelBackendBase instance or None, \
                default: 'multiprocessing'
            Specify the parallelization backend implementation.
            Supported backends are:

            - "multiprocessing" used by default, can induce some
              communication and memory overhead when exchanging input and
              output data with the worker Python processes.
            - "threading" is a very low-overhead backend but it suffers
              from the Python Global Interpreter Lock if the called function
              relies a lot on Python objects. "threading" is mostly useful
              when the execution bottleneck is a compiled extension that
              explicitly releases the GIL (for instance a Cython loop wrapped
              in a "with nogil" block or an expensive call to a library such
              as NumPy).
            - finally, you can register backends by calling
              register_parallel_backend. This will allow you to implement
              a backend of your liking.
        verbose: int, optional
            The verbosity level: if non zero, progress messages are
            printed. Above 50, the output is sent to stdout.
            The frequency of the messages increases with the verbosity level.
            If it more than 10, all iterations are reported.
        timeout: float, optional
            Timeout limit for each task to complete.  If any task takes longer
            a TimeOutError will be raised. Only applied when n_jobs != 1
        pre_dispatch: {'all', integer, or expression, as in '3*n_jobs'}
            The number of batches (of tasks) to be pre-dispatched.
            Default is '2*n_jobs'. When batch_size="auto" this is reasonable
            default and the multiprocessing workers should never starve.
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

            - a folder pointed by the JOBLIB_TEMP_FOLDER environment
              variable,
            - /dev/shm if the folder exists and is writable: this is a
              RAMdisk filesystem available by default on modern Linux
              distributions,
            - the default system temporary folder that can be
              overridden with TMP, TMPDIR or TEMP environment
              variables, typically /tmp under Unix operating systems.

            Only active when backend="multiprocessing".
        max_nbytes int, str, or None, optional, 1M by default
            Threshold on the size of arrays passed to the workers that
            triggers automated memory mapping in temp_folder. Can be an int
            in Bytes, or a human-readable string, e.g., '1M' for 1 megabyte.
            Use None to disable memmaping of large arrays.
            Only active when backend="multiprocessing".
        mmap_mode: {None, 'r+', 'r', 'w+', 'c'}
            Memmapping mode for numpy arrays passed to workers.
            See 'max_nbytes' parameter documentation for more details.

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
        messages:

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
        child process:

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
        called 3 times before the parallel loop is initiated, and then
        called to generate new data on the fly. In this case the total
        number of iterations cannot be reported in the progress messages:

        >>> from math import sqrt
        >>> from sklearn.externals.joblib import Parallel, delayed
        >>> def producer():
        ...     for i in range(6):
        ...         print('Produced %s' % i)
        ...         yield i
        >>> out = Parallel(n_jobs=2, verbose=100, pre_dispatch='1.5*n_jobs')(
        ...                delayed(sqrt)(i) for i in producer()) #doctest: +SKIP
        Produced 0
        Produced 1
        Produced 2
        [Parallel(n_jobs=2)]: Done 1 jobs     | elapsed:  0.0s
        Produced 3
        [Parallel(n_jobs=2)]: Done 2 jobs     | elapsed:  0.0s
        Produced 4
        [Parallel(n_jobs=2)]: Done 3 jobs     | elapsed:  0.0s
        Produced 5
        [Parallel(n_jobs=2)]: Done 4 jobs     | elapsed:  0.0s
        [Parallel(n_jobs=2)]: Done 5 out of 6 | elapsed:  0.0s remaining: 0.0s
        [Parallel(n_jobs=2)]: Done 6 out of 6 | elapsed:  0.0s finished

    '''
    def __init__(self, n_jobs=1, backend=None, verbose=0, timeout=None,
                 pre_dispatch='2 * n_jobs', batch_size='auto',
                 temp_folder=None, max_nbytes='1M', mmap_mode='r'):
        active_backend, default_n_jobs = get_active_backend()
        if backend is None and n_jobs == 1:
            # If we are under a parallel_backend context manager, look up
            # the default number of jobs and use that instead:
            n_jobs = default_n_jobs
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.timeout = timeout
        self.pre_dispatch = pre_dispatch

        if isinstance(max_nbytes, _basestring):
            max_nbytes = memstr_to_bytes(max_nbytes)

        self._backend_args = dict(
            max_nbytes=max_nbytes,
            mmap_mode=mmap_mode,
            temp_folder=temp_folder,
            verbose=max(0, self.verbose - 50),
        )
        if DEFAULT_MP_CONTEXT is not None:
            self._backend_args['context'] = DEFAULT_MP_CONTEXT

        if backend is None:
            backend = active_backend
        elif isinstance(backend, ParallelBackendBase):
            # Use provided backend as is
            pass
        elif hasattr(backend, 'Pool') and hasattr(backend, 'Lock'):
            # Make it possible to pass a custom multiprocessing context as
            # backend to change the start method to forkserver or spawn or
            # preload modules on the forkserver helper process.
            self._backend_args['context'] = backend
            backend = MultiprocessingBackend()
        else:
            try:
                backend_factory = BACKENDS[backend]
            except KeyError:
                raise ValueError("Invalid backend: %s, expected one of %r"
                                 % (backend, sorted(BACKENDS.keys())))
            backend = backend_factory()

        if (batch_size == 'auto' or isinstance(batch_size, Integral) and
                batch_size > 0):
            self.batch_size = batch_size
        else:
            raise ValueError(
                "batch_size must be 'auto' or a positive integer, got: %r"
                % batch_size)

        self._backend = backend
        self._output = None
        self._jobs = list()
        self._managed_backend = False

        # This lock is used coordinate the main thread of this process with
        # the async callback thread of our the pool.
        self._lock = threading.Lock()

    def __enter__(self):
        self._managed_backend = True
        self._initialize_backend()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._terminate_backend()
        self._managed_backend = False

    def _initialize_backend(self):
        """Build a process or thread pool and return the number of workers"""
        try:
            n_jobs = self._backend.configure(n_jobs=self.n_jobs, parallel=self,
                                             **self._backend_args)
            if self.timeout is not None and not self._backend.supports_timeout:
                warnings.warn(
                    'The backend class {!r} does not support timeout. '
                    "You have set 'timeout={}' in Parallel but "
                    "the 'timeout' parameter will not be used.".format(
                        self._backend.__class__.__name__,
                        self.timeout))

        except FallbackToBackend as e:
            # Recursively initialize the backend in case of requested fallback.
            self._backend = e.backend
            n_jobs = self._initialize_backend()

        return n_jobs

    def _effective_n_jobs(self):
        if self._backend:
            return self._backend.effective_n_jobs(self.n_jobs)
        return 1

    def _terminate_backend(self):
        if self._backend is not None:
            self._backend.terminate()

    def _dispatch(self, batch):
        """Queue the batch for computing, with or without multiprocessing

        WARNING: this method is not thread-safe: it should be only called
        indirectly via dispatch_one_batch.

        """
        # If job.get() catches an exception, it closes the queue:
        if self._aborting:
            return

        self.n_dispatched_tasks += len(batch)
        self.n_dispatched_batches += 1

        dispatch_timestamp = time.time()
        cb = BatchCompletionCallBack(dispatch_timestamp, len(batch), self)
        job = self._backend.apply_async(batch, callback=cb)
        self._jobs.append(job)

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
        if self.batch_size == 'auto':
            batch_size = self._backend.compute_batch_size()
        else:
            # Fixed batch size strategy
            batch_size = self.batch_size

        with self._lock:
            tasks = BatchedCalls(itertools.islice(iterator, batch_size))
            if len(tasks) == 0:
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

        # Original job iterator becomes None once it has been fully
        # consumed : at this point we know the total number of jobs and we are
        # able to display an estimation of the remaining time based on already
        # completed jobs. Otherwise, we simply display the number of completed
        # tasks.
        if self._original_iterator is not None:
            if _verbosity_filter(self.n_dispatched_batches, self.verbose):
                return
            self._print('Done %3i tasks      | elapsed: %s',
                        (self.n_completed_tasks,
                         short_format_time(elapsed_time), ))
        else:
            index = self.n_completed_tasks
            # We are finished dispatching
            total_tasks = self.n_dispatched_tasks
            # We always display the first loop
            if not index == 0:
                # Display depending on the number of remaining items
                # A message as soon as we finish dispatching, cursor is 0
                cursor = (total_tasks - index + 1 -
                          self._pre_dispatch_amount)
                frequency = (total_tasks // self.verbose) + 1
                is_last_item = (index + 1 == total_tasks)
                if (is_last_item or cursor % frequency):
                    return
            remaining_time = (elapsed_time / index) * \
                             (self.n_dispatched_tasks - index * 1.0)
            # only display status if remaining time is greater or equal to 0
            self._print('Done %3i out of %3i | elapsed: %s remaining: %s',
                        (index,
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
                if getattr(self._backend, 'supports_timeout', False):
                    self._output.extend(job.get(timeout=self.timeout))
                else:
                    self._output.extend(job.get())

            except BaseException as exception:
                # Note: we catch any BaseException instead of just Exception
                # instances to also include KeyboardInterrupt.

                # Stop dispatching any new job in the async callback thread
                self._aborting = True

                # If the backend allows it, cancel or kill remaining running
                # tasks without waiting for the results as we will raise
                # the exception we got back to the caller instead of returning
                # any result.
                backend = self._backend
                if (backend is not None and
                        hasattr(backend, 'abort_everything')):
                    # If the backend is managed externally we need to make sure
                    # to leave it in a working state to allow for future jobs
                    # scheduling.
                    ensure_ready = self._managed_backend
                    backend.abort_everything(ensure_ready=ensure_ready)

                if not isinstance(exception, TransportableException):
                    raise
                else:
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

                    raise exception

    def __call__(self, iterable):
        if self._jobs:
            raise ValueError('This Parallel instance is already running')
        # A flag used to abort the dispatching of jobs in case an
        # exception is found
        self._aborting = False
        if not self._managed_backend:
            n_jobs = self._initialize_backend()
        else:
            n_jobs = self._effective_n_jobs()

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
            if not self._managed_backend:
                self._terminate_backend()
            self._jobs = list()
        output = self._output
        self._output = None
        return output

    def __repr__(self):
        return '%s(n_jobs=%s)' % (self.__class__.__name__, self.n_jobs)
