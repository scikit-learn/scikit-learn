"""
Helpers for embarassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

import sys
import functools
import time
import threading
import itertools
try:
    import cPickle as pickle
except:
    import pickle

try:
    import multiprocessing
except ImportError:
    multiprocessing = None

from .format_stack import format_exc, format_outer_frames
from .logger import Logger, short_format_time
from .my_exceptions import TransportableException, _mk_exception


################################################################################
class WorkerInterrupt(Exception):
    """ An exception that is not KeyboardInterrupt to allow subprocesses
        to be interrupted.
    """
    pass


################################################################################
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


################################################################################
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


################################################################################
class ImmediateApply(object):
    """ A non-delayed apply function.
    """
    def __init__ (self, func, args, kwargs):
        # Don't delay the application, to avoid keeping the input
        # arguments in memory
        self.results = func(*args, **kwargs)

    def get (self):
        return self.results


################################################################################
class CallBack(object):
    """ Callback used by parallel: it is used for progress reporting, and 
        to add data to be processed
    """
    def __init__(self, index, parallel):
        self.parallel = parallel
        self.index = index

    def __call__(self, out):
        if self.parallel.verbose:
            self.print_progress()
        if self.parallel._iterable:
            self.parallel.dispatch_next()

    def print_progress(self):
        # XXX: Not using the logger framework: need to
        # learn to use logger better.
        n_jobs = len(self.parallel._pool._pool)
        if self.parallel.n_dispatched > 2*n_jobs:
            # Report less often
            if not self.index % n_jobs == 0:
                return
        elapsed_time = time.time() - self.parallel._start_time
        remaining_time = (elapsed_time/(self.index + 1)*
                    (self.parallel.n_dispatched - self.index - 1.))
        if self.parallel._iterable:
            # The object is still building it's job list
            total = "%3i+" % self.parallel.n_dispatched
        else:
            total = "%3i " % self.parallel.n_dispatched

        if self.parallel.verbose < 50:
            writer = sys.stderr.write
        else:
            writer = sys.stdout.write
        writer('[%s]: Done %3i out of %s |elapsed: %s remaining: %s\n'
                % (self.parallel,
                    self.index+1, 
                    total, 
                    short_format_time(elapsed_time),
                    short_format_time(remaining_time),
                    ))


################################################################################
class Parallel(Logger):
    ''' Helper class for readable parallel mapping.

        Parameters
        -----------
        n_jobs: int
            The number of jobs to use for the computation. If -1 all CPUs
            are used. If 1 is given, no parallel computing code is used
            at all, which is useful for debuging.
        verbose: int, optional
            The verbosity level. If 1 is given, the elapsed time as well
            as the estimated remaining time are displayed. Above 100, the
            output is sent to stdout.
        pre_dispatch: {'all', integer, or expression, as in '3*n_jobs'}
            The amount of jobs to be pre-dispatched. Default is 'all',
            but it may be memory consuming, for instance if each job 
            involves a lot of a data.
        
        Notes
        -----

        This object uses the multiprocessing module to compute in
        parallel the application of a function to many different
        arguments. The main functionnality it brings in addition to 
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
        >>> from scikits.learn.externals.joblib import Parallel, delayed
        >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

        Reshaping the output when the function has several return
        values:
        
        >>> from math import modf
        >>> from scikits.learn.externals.joblib import Parallel, delayed
        >>> r = Parallel(n_jobs=1)(delayed(modf)(i/2.) for i in range(10))
        >>> res, i = zip(*r)
        >>> res
        (0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5)
        >>> i
        (0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0)
       
        The progress meter::

            >>> from time import sleep
            >>> from scikits.learn.externals.joblib import Parallel, delayed
            >>> r = Parallel(n_jobs=2, verbose=1)(delayed(sleep)(.1) for _ in range(10)) #doctest: +SKIP
            [Parallel(n_jobs=2)]: Done   1 out of  10 |elapsed:    0.1s remaining:    0.9s
            [Parallel(n_jobs=2)]: Done   3 out of  10 |elapsed:    0.2s remaining:    0.5s
            [Parallel(n_jobs=2)]: Done   5 out of  10 |elapsed:    0.3s remaining:    0.3s
            [Parallel(n_jobs=2)]: Done   7 out of  10 |elapsed:    0.4s remaining:    0.2s
            [Parallel(n_jobs=2)]: Done   9 out of  10 |elapsed:    0.5s remaining:    0.1s

        Traceback example, note how the ligne of the error is indicated 
        as well as the values of the parameter passed to the function that
        triggered the exception, eventhough the traceback happens in the 
        child process::

         >>> from string import atoi
         >>> from scikits.learn.externals.joblib import Parallel, delayed
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
        number of iterations reported is underestimated::

         >>> from math import sqrt
         >>> from scikits.learn.externals.joblib import Parallel, delayed

         >>> def producer():
         ...     for i in range(6):
         ...         print 'Produced %s' % i
         ...         yield i
         
         >>> out = Parallel(n_jobs=2, verbose=100, pre_dispatch='1.5*n_jobs')(
         ...                         delayed(sqrt)(i) for i in producer()) #doctest: +ELLIPSIS
         Produced 0
         Produced 1
         Produced 2
         [Parallel(n_jobs=2)]: Done   1 out of   3+ |elapsed:   ...s remaining:   ...s
         Produced 3
         [Parallel(n_jobs=2)]: Done ... out of   4+ |elapsed:   ...s remaining:   ...s
         ...

    '''
    def __init__(self, n_jobs=None, verbose=0, pre_dispatch='all'):
        self.verbose = verbose
        self.n_jobs  = n_jobs
        self.pre_dispatch   = pre_dispatch
        self._pool   = None
        # Not starting the pool in the __init__ is a design decision, to
        # be able to close it ASAP, and not burden the user with closing
        # it.
        self._output = None
        self._jobs   = list()


    def dispatch(self, func, args, kwargs):
        """ Queue the function for computing, with or without multiprocessing 
        """
        if self._pool is None:
            job = ImmediateApply(func, args, kwargs)
            if self.verbose:
                print '[%s]: Done job %3i | elapsed: %s' % (
                        self, len(self._jobs)+1,
                        short_format_time(time.time() - self._start_time)
                    )
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
                    # Capture exception to add information on 
                    # the local stack in addition to the distant
                    # stack
                    this_report = format_outer_frames(
                                            context=10,
                                            stack_start=1,
                                            )
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
        if n_jobs == -1 and multiprocessing is not None:
            n_jobs = multiprocessing.cpu_count()

        # The list of exceptions that we will capture
        self.exceptions = [TransportableException]
        if n_jobs is None or multiprocessing is None or n_jobs == 1:
            n_jobs = 1
            self._pool = None
        else:
            self._pool = multiprocessing.Pool(n_jobs)
            self._lock = threading.Lock()
            # We are using multiprocessing, we also want to capture
            # KeyboardInterrupts
            self.exceptions.extend([KeyboardInterrupt, WorkerInterrupt])

        if self.pre_dispatch == 'all' or n_jobs == 1:
            self._iterable = None
        else:
            self._iterable = iterable
            self._dispatch_amount = 0
            pre_dispatch = self.pre_dispatch
            if hasattr(pre_dispatch, 'endswith'):
                pre_dispatch = eval(pre_dispatch)
            pre_dispatch = int(pre_dispatch)
            iterable = itertools.islice(iterable, pre_dispatch)

        self._start_time = time.time()
        self.n_dispatched = 0
        try:
            for function, args, kwargs in iterable:
                self.dispatch(function, args, kwargs)
            
            self.retrieve()
        finally:
            if n_jobs > 1:
                self._pool.close()
                self._pool.join()
            self._jobs = list()
        output = self._output
        self._output = None
        return output


    def __repr__(self):
        return '%s(n_jobs=%s)' % (
                    self.__class__.__name__,
                    self.n_jobs,
                )



