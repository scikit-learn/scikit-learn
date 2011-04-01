"""
Test the parallel module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2010-2011 Gael Varoquaux
# License: BSD Style, 3 clauses.

import time
try:
    import cPickle as pickle
    PickleError = TypeError
except:
    import pickle
    PickleError = pickle.PicklingError

from ..parallel import Parallel, delayed, SafeFunction, WorkerInterrupt
from ..my_exceptions import JoblibException

import nose

################################################################################

def division(x, y):
    return x/y

def square(x):
    return x**2

def exception_raiser(x):
    if x == 7:
        raise ValueError
    return x

def interrupt_raiser(x):
    time.sleep(.05)
    raise KeyboardInterrupt

def f(x, y=0, z=0):
    """ A module-level function so that it can be spawn with
    multiprocessing.
    """
    return x**2 + y + z

################################################################################
# Test parallel
def test_simple_parallel():
    X = range(10)
    for n_jobs in (1, 2, -1):
        yield (nose.tools.assert_equal, [square(x) for x in X],
                        Parallel(n_jobs=-1)(delayed(square)(x) for x in X))


def test_parallel_kwargs():
    """ Check the keyword argument processing of pmap.
    """
    lst = range(10)
    for n_jobs in (1, 4):
        yield (nose.tools.assert_equal, 
               [f(x, y=1) for x in lst], 
               Parallel(n_jobs=n_jobs)(delayed(f)(x, y=1) for x in lst)
              )

        
def test_parallel_pickling():
    """ Check that pmap captures the errors when it is passed an object
        that cannot be pickled.
    """
    def g(x):
        return x**2
    nose.tools.assert_raises(PickleError,
                             Parallel(), 
                             (delayed(g)(x) for x in range(10))
                            )


def test_error_capture():
    """ Check that error are captured, and that correct exceptions
        are raised.
    """
    try:
        import multiprocessing
    except ImportError:
        multiprocessing = None
    if multiprocessing is not None:
        # A JoblibException will be raised only if there is indeed
        # multiprocessing
        nose.tools.assert_raises(JoblibException,
                                Parallel(n_jobs=2),
                    [delayed(division)(x, y) for x, y in zip((0, 1), (1, 0))],
                        )
        nose.tools.assert_raises(WorkerInterrupt,
                                    Parallel(n_jobs=2),
                        [delayed(interrupt_raiser)(x) for x in (1, 0)],
                            )
    else:
        nose.tools.assert_raises(KeyboardInterrupt,
                                    Parallel(n_jobs=2),
                        [delayed(interrupt_raiser)(x) for x in (1, 0)],
                            )
    nose.tools.assert_raises(ZeroDivisionError,
                                Parallel(n_jobs=2),
                    [delayed(division)(x, y) for x, y in zip((0, 1), (1, 0))],
                        )
    try:
        Parallel(n_jobs=1)(
                    delayed(division)(x, y) for x, y in zip((0, 1), (1, 0)))
    except Exception, e:
        pass
    nose.tools.assert_false(isinstance(e, JoblibException))


class Counter(object):
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2

    def __call__(self, i):
        self.list1.append(i)
        nose.tools.assert_equal(len(self.list1), len(self.list2))


def test_dispatch_one_job():
    """ Test that with only one job, Parallel does act as a iterator.
    """
    produced = list()
    consumed = list()
    def producer():
        for i in range(6):
            produced.append(i)
            yield i
    consumer = Counter(list1=consumed, list2=consumed)

    Parallel(n_jobs=1)(delayed(consumer)(x) for x in producer())


def test_exception_dispatch():
    "Make sure that exception raised during dispatch are indeed captured"
    nose.tools.assert_raises(
            ValueError,
            Parallel(n_jobs=6, pre_dispatch=16, verbose=0),
                    (delayed(exception_raiser)(i) for i in range(30)),
            )


################################################################################
# Test helpers
def test_joblib_exception():
    # Smoke-test the custom exception
    e = JoblibException('foobar')
    # Test the repr
    repr(e)
    # Test the pickle
    pickle.dumps(e)


def test_safe_function():
    safe_division = SafeFunction(division)
    nose.tools.assert_raises(JoblibException, safe_division, 1, 0)

