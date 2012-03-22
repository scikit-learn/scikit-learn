"""
Helpers for embarrassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

import os
import sys
import warnings
from math import sqrt
import functools
import time
import threading
import itertools
try:
    import cPickle as pickle
except:
    import pickle

# Obtain possible configuration from the environment, assuming 1 (on)
# by default, upon 0 set to None. Should instructively fail if some non
# 0/1 value is set.
multiprocessing = int(os.environ.get('JOBLIB_MULTIPROCESSING', 1)) or None
if multiprocessing:
    try:
        import multiprocessing
    except ImportError:
        multiprocessing = None

from .format_stack import format_exc, format_outer_frames
from .logger import Logger, short_format_time
from .my_exceptions import TransportableException, _mk_exception


###############################################################################
# CPU that works also when multiprocessing is not installed (python2.5)
def cpu_count():
    """ Return the number of CPUs.
    """
    if multiprocessing is None:
        return 1
    return multiprocessing.cpu_count()


###############################################################################
# For verbosity

def _verbosity_filter(index, verbose):
    """ Returns False for indices increasingly appart, the distance
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
def delayed(function):
    """ Decorator used to capture the arguments of a function.
    """
    # Try to pickle the input function, to catch the problems early when
    # using with multiprocessing
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
        if self.parallel._iterable:
            self.parallel.dispatch_next()


###############################################################################
class Parallel(Logger):
    ''' Helper class for readable parallel mapping.

        Parameters
        -----------
        n_jobs: int
            The number of jobs to use for the computation. If -1 all CPUs
            are used. If 1 is given, no parallel computing code is used
            at all, which is useful for debuging. For n_jobs below -1,
            (n_cpus + 1 - n_jobs) are used. Thus for n_jobs = -2, all
            CPUs but one are used.
        verbose: int, optional
            The verbosity level: if non zero, progress messages are
            printed. Above 50, the output is sent to stdout.
            The frequency of the messages increases with the verbosity level.
            If it more than 10, all iterations are reported.
        pre_dispatch: {'all', integer, or expression, as in '3*n_jobs'}
            The amount of jobs to be pre-dispatched. Default is 'all',
            but it may be memory consuming, for instance if each job
            involves a lot of a data.

        Notes
        -----

        This object uses the multiprocessing module to compute in
        parallel the application of a function to many different
        arguments. The main functionality it brings in addition to
        using the raw multiprocessing API are (see examples for details):

            * More readable code, in particular since it avoids
              constructing list of arguments.

            * Easier debuging:
                - informative tracebacks even when the error happens on
                  the client side
                - using 'n_jobs=1' enables to turn off parallel computing
                  for debuging without changing the codepath
                - early capture of pickling errors

            * An optional progress meter.

            * Interruption of multiprocesses jobs with 'Ctrl-C'

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

         >>> from string import atoi
         >>> from sklearn.externals.joblib import Parallel, delayed
         >>> Parallel(n_jobs=2)(delayed(atoi)(n) for n in ('1', '300', 30)) #doctest: +SKIP
         #...
         ---------------------------------------------------------------------------
         Sub-process traceback:
         ---------------------------------------------------------------------------
         TypeError                                          Fri Jul  2 20:32:05 2010
         PID: 4151                                     Python 2.6.5: /usr/bin/python
         ...........................................................................
         /usr/lib/python2.6/string.pyc in atoi(s=30, base=10)
             398     is chosen from the leading characters of s, 0 for octal, 0x or
             399     0X for hexadecimal.  If base is 16, a preceding 0x or 0X is
             400     accepted.
             401
             402     """
         --> 403     return _int(s, base)
             404
             405
             406 # Convert string to long integer
             407 def atol(s, base=10):

         TypeError: int() can't convert non-string with explicit base
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
         ...         print 'Produced %s' % i
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
    def __init__(self, n_jobs=None, verbose=0, pre_dispatch='all'):
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self._pool = None
        # Not starting the pool in the __init__ is a design decision, to be
        # able to close it ASAP, and not burden the user with closing it.
        self._output = None
        self._jobs = list()

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
            self._lock.acquire()
            # If job.get() catches an exception, it closes the queue:
            try:
                job = self._pool.apply_async(SafeFunction(func), args,
                            kwargs, callback=CallBack(self.n_dispatched, self))
                self._jobs.append(job)
                self.n_dispatched += 1
            except AssertionError:
                print '[Parallel] Pool seems closed'
            finally:
                self._lock.release()

    def dispatch_next(self):
        """ Dispatch more data for parallel processing
        """
        self._dispatch_amount += 1
        while self._dispatch_amount:
            try:
                # XXX: possible race condition shuffling the order of
                # dispatchs in the next two lines.
                func, args, kwargs = self._iterable.next()
                self.dispatch(func, args, kwargs)
                self._dispatch_amount -= 1
            except ValueError:
                """ Race condition in accessing a generator, we skip,
                    the dispatch will be done later.
                """
            except StopIteration:
                self._iterable = None
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
           of time, controled by self.verbose.
        """
        if not self.verbose:
            return
        elapsed_time = time.time() - self._start_time

        # This is heuristic code to print only 'verbose' times a messages
        # The challenge is that we may not know the queue length
        if self._iterable:
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
        while self._jobs:
            # We need to be careful: the job queue can be filling up as
            # we empty it
            if hasattr(self, '_lock'):
                self._lock.acquire()
            job = self._jobs.pop(0)
            if hasattr(self, '_lock'):
                self._lock.release()
            try:
                self._output.append(job.get())
            except tuple(self.exceptions), exception:
                if isinstance(exception,
                        (KeyboardInterrupt, WorkerInterrupt)):
                    # We have captured a user interruption, clean up
                    # everything
                    if hasattr(self, '_pool'):
                        self._pool.close()
                        self._pool.terminate()
                    raise exception
                elif isinstance(exception, TransportableException):
                    # Capture exception to add information on the local stack
                    # in addition to the distant stack
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

    def __call__(self, iterable):
        if self._jobs:
            raise ValueError('This Parallel instance is already running')
        n_jobs = self.n_jobs
        if n_jobs < 0 and multiprocessing is not None:
            n_jobs = max(multiprocessing.cpu_count() + 1 + n_jobs, 1)

        # The list of exceptions that we will capture
        self.exceptions = [TransportableException]
        if n_jobs is None or multiprocessing is None or n_jobs == 1:
            n_jobs = 1
            self._pool = None
        else:
            if multiprocessing.current_process()._daemonic:
                # Daemonic processes cannot have children
                n_jobs = 1
                self._pool = None
                warnings.warn(
                    'Parallel loops cannot be nested, setting n_jobs=1',
                    stacklevel=2)
            else:
                self._pool = multiprocessing.Pool(n_jobs)
                self._lock = threading.Lock()
                # We are using multiprocessing, we also want to capture
                # KeyboardInterrupts
                self.exceptions.extend([KeyboardInterrupt, WorkerInterrupt])

        if self.pre_dispatch == 'all' or n_jobs == 1:
            self._iterable = None
            self._pre_dispatch_amount = 0
        else:
            self._iterable = iterable
            self._dispatch_amount = 0
            pre_dispatch = self.pre_dispatch
            if hasattr(pre_dispatch, 'endswith'):
                pre_dispatch = eval(pre_dispatch)
            self._pre_dispatch_amount = pre_dispatch = int(pre_dispatch)
            iterable = itertools.islice(iterable, pre_dispatch)

        self._start_time = time.time()
        self.n_dispatched = 0
        try:
            for function, args, kwargs in iterable:
                self.dispatch(function, args, kwargs)

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
                self._pool.join()
            self._jobs = list()
        output = self._output
        self._output = None
        return output

    def __repr__(self):
        return '%s(n_jobs=%s)' % (self.__class__.__name__, self.n_jobs)
