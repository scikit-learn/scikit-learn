"""
Helpers for embarassingly parallel code.
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

import sys
import functools
import time
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
from .my_exceptions import JoblibException, _mk_exception

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
        except:
            e_type, e_value, e_tb = sys.exc_info()
            text = format_exc(e_type, e_value, e_tb, context=10,
                             tb_offset=1)
            exception = _mk_exception(e_type)[0]
            raise exception(text)

def print_progress(msg, index, total, start_time, n_jobs=1):
    # XXX: Not using the logger framework: need to
    # learn to use logger better.
    if total > 2*n_jobs:
        # Report less often
        if not index % n_jobs == 0:
            return
    elapsed_time = time.time() - start_time
    remaining_time = (elapsed_time/(index + 1)*
                (total - index - 1.))
    sys.stderr.write('[%s]: Done %3i out of %3i |elapsed: %s remaining: %s\n'
            % (msg,
                index+1,
                total,
                short_format_time(elapsed_time),
                short_format_time(remaining_time),
                ))


################################################################################
def delayed(function):
    """ Decorator used to capture the arguments of a function.
    """
    # Try to pickle the input function, to catch the problems early when
    # using with multiprocessing
    pickle.dumps(function)

    @functools.wraps(function)
    def delayed_function(*args, **kwargs):
        return function, args, kwargs
    return delayed_function


class LazyApply (object):
    """
    Lazy version of the apply builtin function.
    """
    def __init__ (self, func, args, kwargs):
        self.func   = func
        self.args   = args
        self.kwargs = kwargs

    def get (self):
        return self.func(*self.args, **self.kwargs)


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
            as the estimated remaining time are displayed.

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

        Examples
        --------

        A simple example:

        >>> from math import sqrt
        >>> from joblib import Parallel, delayed
        >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

        Reshaping the output when the function has several return
        values:

        >>> from math import modf
        >>> from joblib import Parallel, delayed
        >>> r = Parallel(n_jobs=1)(delayed(modf)(i/2.) for i in range(10))
        >>> res, i = zip(*r)
        >>> res
        (0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5)
        >>> i
        (0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0)

        The progress meter::

            >>> from time import sleep
            >>> from joblib import Parallel, delayed
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
         >>> from joblib import Parallel, delayed
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

    '''
    def __init__(self, n_jobs=None, verbose=0):
        self.verbose = verbose
        self.n_jobs = n_jobs
        # Not starting the pool in the __init__ is a design decision, to
        # be able to close it ASAP, and not burden the user with closing
        # it.


    def __call__(self, iterable):
        n_jobs = self.n_jobs
        if n_jobs == -1:
            if multiprocessing is None:
                 n_jobs = 1
            else:
                n_jobs = multiprocessing.cpu_count()

        if n_jobs is None or multiprocessing is None or n_jobs == 1:
            n_jobs = 1
            apply = LazyApply
        else:
            pool = multiprocessing.Pool(n_jobs)
            apply = pool.apply_async

        output = list()
        start_time = time.time()
        try:
            for index, (function, args, kwargs) in enumerate(iterable):
                function = SafeFunction(function)
                output.append(apply(function, args, kwargs))
                if self.verbose and n_jobs == 1:
                    print '[%s]: Done job %3i | elapsed: %s' % (
                            self, index,
                            short_format_time(time.time() - start_time)
                        )

            start_time = time.time()
            jobs = output
            output = list()
            for index, job in enumerate(jobs):
                try:
                    output.append(job.get())
                    if self.verbose:
                        print_progress(self, index, len(jobs), start_time,
                                       n_jobs=n_jobs)
                except JoblibException, exception:
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
                    # No need to convert this to a JoblibException, the
                    # SafeFunction already did it
                    raise exception.__class__(report)
        finally:
            if n_jobs > 1:
                pool.close()
                pool.join()
        return output


    def __repr__(self):
        return '%s(n_jobs=%s)' % (
                    self.__class__.__name__,
                    self.n_jobs,
                )



