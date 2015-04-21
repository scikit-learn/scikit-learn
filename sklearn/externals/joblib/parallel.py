"""
Helpers for embarrassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

import os
import sys
import gc
import warnings
from collections import Sized
from math import sqrt
import functools
import time
import threading
import itertools

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


###############################################################################
# CPU that works also when multiprocessing is not installed (python2.5)
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
class ImmediateApply(object):
    """ A non-delayed apply function.
    """
    def __init__(self, func, args, kwargs):
        # Don't delay the application, to avoid keeping the input
        # arguments in memory
        self.results = func(*args, **kwargs)

    def get(self):
        return self.results


###############################################################################
class CallBack(object):
    """ Callback used by parallel: it is used for progress reporting, and
        to add data to be processed
    """
    def __init__(self, index, parallel):
        self.parallel = parallel
        self.index = index

    def __call__(self, out):
        self.parallel.print_progress(self.index)
        if self.parallel._original_iterable:
            self.parallel.dispatch_next()


class LockedIterator(object):
    """Wrapper to protect a thread-unsafe iterable against concurrent access.

    A Python generator is not thread-safe by default and will raise
    ValueError("generator already executing") if two threads consume it
    concurrently.

    In joblib this could typically happen when the passed iterator is a
    generator expression and pre_dispatch != 'all'. In that case a callback is
    passed to the multiprocessing apply_async call and helper threads will
    trigger the consumption of the source iterable in the dispatch_next
    method.

    """
    def __init__(self, it):
        self._lock = threading.Lock()
        self._it = iter(it)

    def __iter__(self):
        return self

    def next(self):
        with self._lock:
            return next(self._it)

    # For Python 3 compat
    __next__ = next



###############################################################################
class Parallel(Logger):
    ''' Helper class for readable parallel mapping.

        Parameters
        -----------
        n_jobs : int
            The number of jobs to use for the computation. If -1 all CPUs
            are used. If 1 is given, no parallel computing code is used
            at all, which is useful for debugging. For n_jobs below -1,
            (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all
            CPUs but one are used.
        backend : str or None
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
        verbose : int, optional
            The verbosity level: if non zero, progress messages are
            printed. Above 50, the output is sent to stdout.
            The frequency of the messages increases with the verbosity level.
            If it more than 10, all iterations are reported.
        pre_dispatch : {'all', integer, or expression, as in '3*n_jobs'}
            The amount of jobs to be pre-dispatched. Default is 'all',
            but it may be memory consuming, for instance if each job
            involves a lot of a data.
        temp_folder : str, optional
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
        max_nbytes : int, str, or None, optional, 100e6 (100MB) by default
            Threshold on the size of arrays passed to the workers that
            triggers automated memory mapping in temp_folder. Can be an int
            in Bytes, or a human-readable string, e.g., '1M' for 1 megabyte.
            Use None to disable memmaping of large arrays.
            Only active when backend="multiprocessing".
        mmap_mode : 'r', 'r+' or 'c'
            Mode for the created memmap datastructure. See the documentation of
            numpy.memmap for more details. Note: 'w+' is coerced to 'r+'
            automatically to avoid zeroing the data on unpickling.

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
    def __init__(self, n_jobs=1, backend=None, verbose=0, pre_dispatch='all',
                 temp_folder=None, max_nbytes=100e6, mmap_mode='r'):
        self.verbose = verbose
        self._mp_context = None
        if backend is None:
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
        self.pre_dispatch = pre_dispatch
        self._pool = None
        self._temp_folder = temp_folder
        if isinstance(max_nbytes, _basestring):
            self._max_nbytes = 1024 * memstr_to_kbytes(max_nbytes)
        else:
            self._max_nbytes = max_nbytes
        self._mmap_mode = mmap_mode
        # Not starting the pool in the __init__ is a design decision, to be
        # able to close it ASAP, and not burden the user with closing it.
        self._output = None
        self._jobs = list()
        # A flag used to abort the dispatching of jobs in case an
        # exception is found
        self._aborting = False

    def dispatch(self, func, args, kwargs):
        """ Queue the function for computing, with or without multiprocessing
        """
        if self._pool is None:
            job = ImmediateApply(func, args, kwargs)
            index = len(self._jobs)
            if not _verbosity_filter(index, self.verbose):
                self._print('Done %3i jobs       | elapsed: %s',
                        (index + 1,
                            short_format_time(time.time() - self._start_time)
                        ))
            self._jobs.append(job)
            self.n_dispatched += 1
        else:
            # If job.get() catches an exception, it closes the queue:
            if self._aborting:
                return
            try:
                self._lock.acquire()
                job = self._pool.apply_async(SafeFunction(func), args,
                            kwargs, callback=CallBack(self.n_dispatched, self))
                self._jobs.append(job)
                self.n_dispatched += 1
            except AssertionError:
                print('[Parallel] Pool seems closed')
            finally:
                self._lock.release()

    def dispatch_next(self):
        """ Dispatch more data for parallel processing
        """
        self._dispatch_amount += 1
        while self._dispatch_amount:
            try:
                # XXX: possible race condition shuffling the order of
                # dispatches in the next two lines.
                func, args, kwargs = next(self._original_iterable)
                self.dispatch(func, args, kwargs)
                self._dispatch_amount -= 1
            except ValueError:
                """ Race condition in accessing a generator, we skip,
                    the dispatch will be done later.
                """
            except StopIteration:
                self._iterating = False
                self._original_iterable = None
                return

    def _print(self, msg, msg_args):
        """ Display the message on stout or stderr depending on verbosity
        """
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

    def print_progress(self, index):
        """Display the process of the parallel execution only a fraction
           of time, controlled by self.verbose.
        """
        if not self.verbose:
            return
        elapsed_time = time.time() - self._start_time

        # This is heuristic code to print only 'verbose' times a messages
        # The challenge is that we may not know the queue length
        if self._original_iterable:
            if _verbosity_filter(index, self.verbose):
                return
            self._print('Done %3i jobs       | elapsed: %s',
                        (index + 1,
                         short_format_time(elapsed_time),
                        ))
        else:
            # We are finished dispatching
            queue_length = self.n_dispatched
            # We always display the first loop
            if not index == 0:
                # Display depending on the number of remaining items
                # A message as soon as we finish dispatching, cursor is 0
                cursor = (queue_length - index + 1
                          - self._pre_dispatch_amount)
                frequency = (queue_length // self.verbose) + 1
                is_last_item = (index + 1 == queue_length)
                if (is_last_item or cursor % frequency):
                    return
            remaining_time = (elapsed_time / (index + 1) *
                        (self.n_dispatched - index - 1.))
            self._print('Done %3i out of %3i | elapsed: %s remaining: %s',
                        (index + 1,
                         queue_length,
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
            # We need to be careful: the job queue can be filling up as
            # we empty it
            if hasattr(self, '_lock'):
                self._lock.acquire()
            job = self._jobs.pop(0)
            if hasattr(self, '_lock'):
                self._lock.release()
            try:
                self._output.append(job.get())
            except tuple(self.exceptions) as exception:
                try:
                    self._aborting = True
                    self._lock.acquire()
                    if isinstance(exception,
                            (KeyboardInterrupt, WorkerInterrupt)):
                        # We have captured a user interruption, clean up
                        # everything
                        if hasattr(self, '_pool'):
                            self._pool.close()
                            self._pool.terminate()
                            # We can now allow subprocesses again
                            os.environ.pop('__JOBLIB_SPAWNED_PARALLEL__', 0)
                        raise exception
                    elif isinstance(exception, TransportableException):
                        # Capture exception to add information on the local
                        # stack in addition to the distant stack
                        this_report = format_outer_frames(context=10,
                                                          stack_start=1)
                        report = """Multiprocessing exception:
    %s
    ---------------------------------------------------------------------------
    Sub-process traceback:
    ---------------------------------------------------------------------------
    %s""" % (
                                this_report,
                                exception.message,
                            )
                        # Convert this to a JoblibException
                        exception_type = _mk_exception(exception.etype)[0]
                        raise exception_type(report)
                    raise exception
                finally:
                    self._lock.release()

    def __call__(self, iterable):
        if self._jobs:
            raise ValueError('This Parallel instance is already running')
        n_jobs = self.n_jobs
        if n_jobs == 0:
            raise ValueError('n_jobs == 0 in Parallel has no meaning')
        if n_jobs < 0 and mp is not None:
            n_jobs = max(mp.cpu_count() + 1 + n_jobs, 1)

        # The list of exceptions that we will capture
        self.exceptions = [TransportableException]
        self._lock = threading.Lock()

        # Whether or not to set an environment flag to track
        # multiple process spawning
        set_environ_flag = False
        if (n_jobs is None or mp is None or n_jobs == 1):
            n_jobs = 1
            self._pool = None
        elif self.backend == 'threading':
            self._pool = ThreadPool(n_jobs)
        elif self.backend == 'multiprocessing':
            if mp.current_process().daemon:
                # Daemonic processes cannot have children
                n_jobs = 1
                self._pool = None
                warnings.warn(
                    'Multiprocessing-backed parallel loops cannot be nested,'
                    ' setting n_jobs=1',
                    stacklevel=2)
            elif threading.current_thread().name != 'MainThread':
                # Prevent posix fork inside in non-main posix threads
                n_jobs = 1
                self._pool = None
                warnings.warn(
                    'Multiprocessing backed parallel loops cannot be nested'
                    ' below threads, setting n_jobs=1',
                    stacklevel=2)
            else:
                already_forked = int(os.environ.get('__JOBLIB_SPAWNED_PARALLEL__', 0))
                if already_forked:
                    raise ImportError('[joblib] Attempting to do parallel computing '
                            'without protecting your import on a system that does '
                            'not support forking. To use parallel-computing in a '
                            'script, you must protect your main loop using "if '
                            "__name__ == '__main__'"
                            '". Please see the joblib documentation on Parallel '
                            'for more information'
                        )

                # Make sure to free as much memory as possible before forking
                gc.collect()

                # Set an environment variable to avoid infinite loops
                set_environ_flag = True
                poolargs = dict(
                    max_nbytes=self._max_nbytes,
                    mmap_mode=self._mmap_mode,
                    temp_folder=self._temp_folder,
                    verbose=max(0, self.verbose - 50),
                    context_id=0,  # the pool is used only for one call
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

        pre_dispatch = self.pre_dispatch
        if isinstance(iterable, Sized):
            # We are given a sized (an object with len). No need to be lazy.
            pre_dispatch = 'all'

        if pre_dispatch == 'all' or n_jobs == 1:
            self._original_iterable = None
            self._pre_dispatch_amount = 0
        else:
            # The dispatch mechanism relies on multiprocessing helper threads
            # to dispatch tasks from the original iterable concurrently upon
            # job completions. As Python generators are not thread-safe we
            # need to wrap it with a lock
            iterable = LockedIterator(iterable)
            self._original_iterable = iterable
            self._dispatch_amount = 0
            if hasattr(pre_dispatch, 'endswith'):
                pre_dispatch = eval(pre_dispatch)
            self._pre_dispatch_amount = pre_dispatch = int(pre_dispatch)

            # The main thread will consume the first pre_dispatch items and
            # the remaining items will later be lazily dispatched by async
            # callbacks upon task completions
            iterable = itertools.islice(iterable, pre_dispatch)

        self._start_time = time.time()
        self.n_dispatched = 0
        try:
            if set_environ_flag:
                # Set an environment variable to avoid infinite loops
                os.environ[JOBLIB_SPAWNED_PROCESS] = '1'
            self._iterating = True
            for function, args, kwargs in iterable:
                self.dispatch(function, args, kwargs)

            if pre_dispatch == "all" or n_jobs == 1:
                # The iterable was consumed all at once by the above for loop.
                # No need to wait for async callbacks to trigger to
                # consumption.
                self._iterating = False
            self.retrieve()
            # Make sure that we get a last message telling us we are done
            elapsed_time = time.time() - self._start_time
            self._print('Done %3i out of %3i | elapsed: %s finished',
                        (len(self._output),
                         len(self._output),
                            short_format_time(elapsed_time)
                        ))

        finally:
            if n_jobs > 1:
                self._pool.close()
                self._pool.terminate()  # terminate does a join()
                if self.backend == 'multiprocessing':
                    os.environ.pop(JOBLIB_SPAWNED_PROCESS, 0)
            self._jobs = list()
        output = self._output
        self._output = None
        return output

    def __repr__(self):
        return '%s(n_jobs=%s)' % (self.__class__.__name__, self.n_jobs)
