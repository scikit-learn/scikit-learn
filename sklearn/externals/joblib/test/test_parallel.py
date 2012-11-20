"""
Test the parallel module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2010-2011 Gael Varoquaux
# License: BSD Style, 3 clauses.

import time
import sys
import io
import os
try:
    import cPickle as pickle
    PickleError = TypeError
except:
    import pickle
    PickleError = pickle.PicklingError


if sys.version_info[0] == 3:
    PickleError = pickle.PicklingError

try:
    # Python 2/Python 3 compat
    unicode('str')
except NameError:
    unicode = lambda s: s


from ..parallel import Parallel, delayed, SafeFunction, WorkerInterrupt, \
        multiprocessing, cpu_count
from ..my_exceptions import JoblibException

import nose


###############################################################################

def division(x, y):
    return x / y


def square(x):
    return x ** 2


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
    return x ** 2 + y + z


###############################################################################
def test_cpu_count():
    assert cpu_count() > 0


###############################################################################
# Test parallel
def test_simple_parallel():
    X = range(5)
    for n_jobs in (1, 2, -1, -2):
        yield (nose.tools.assert_equal, [square(x) for x in X],
                Parallel(n_jobs=-1)(
                        delayed(square)(x) for x in X))
    try:
        # To smoke-test verbosity, we capture stdout
        orig_stdout = sys.stdout
        orig_stderr = sys.stdout
        if sys.version_info[0] == 3:
            sys.stderr = io.StringIO()
            sys.stderr = io.StringIO()
        else:
            sys.stdout = io.BytesIO()
            sys.stderr = io.BytesIO()
        for verbose in (2, 11, 100):
                Parallel(n_jobs=-1, verbose=verbose)(
                        delayed(square)(x) for x in X)
                Parallel(n_jobs=1, verbose=verbose)(
                        delayed(square)(x) for x in X)
                Parallel(n_jobs=2, verbose=verbose, pre_dispatch=2)(
                        delayed(square)(x) for x in X)
    except Exception as e:
        my_stdout = sys.stdout
        my_stderr = sys.stderr
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        print(unicode(my_stdout.getvalue()))
        print(unicode(my_stderr.getvalue()))
        raise e
    finally:
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr


def nested_loop():
    Parallel(n_jobs=2)(delayed(square)(.01) for _ in range(2))


def test_nested_loop():
    Parallel(n_jobs=2)(delayed(nested_loop)() for _ in range(2))


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
        return x ** 2
    nose.tools.assert_raises(PickleError,
                             Parallel(),
                             (delayed(g)(x) for x in range(10))
                            )


def test_error_capture():
    # Check that error are captured, and that correct exceptions
    # are raised.
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
        ex = JoblibException
        Parallel(n_jobs=1)(
                    delayed(division)(x, y) for x, y in zip((0, 1), (1, 0)))
    except Exception:
        # Cannot use 'except as' to maintain Python 2.5 compatibility
        ex = sys.exc_info()[1]
    nose.tools.assert_false(isinstance(ex, JoblibException))


class Counter(object):
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2

    def __call__(self, i):
        self.list1.append(i)
        nose.tools.assert_equal(len(self.list1), len(self.list2))


def consumer(queue, item):
    queue.append('Consumed %s' % item)


def test_dispatch_one_job():
    """ Test that with only one job, Parallel does act as a iterator.
    """
    queue = list()

    def producer():
        for i in range(6):
            queue.append('Produced %i' % i)
            yield i

    Parallel(n_jobs=1)(delayed(consumer)(queue, x) for x in producer())
    nose.tools.assert_equal(queue,
                              ['Produced 0', 'Consumed 0',
                               'Produced 1', 'Consumed 1',
                               'Produced 2', 'Consumed 2',
                               'Produced 3', 'Consumed 3',
                               'Produced 4', 'Consumed 4',
                               'Produced 5', 'Consumed 5']
                               )
    nose.tools.assert_equal(len(queue), 12)


def test_dispatch_multiprocessing():
    """ Check that using pre_dispatch Parallel does indeed dispatch items
        lazily.
    """
    if multiprocessing is None:
        return
    manager = multiprocessing.Manager()
    queue = manager.list()

    def producer():
        for i in range(6):
            queue.append('Produced %i' % i)
            yield i

    Parallel(n_jobs=2, pre_dispatch=3)(delayed(consumer)(queue, i)
                                       for i in producer())
    nose.tools.assert_equal(list(queue)[:4],
            ['Produced 0', 'Produced 1', 'Produced 2',
             'Consumed 0', ])
    nose.tools.assert_equal(len(queue), 12)


def test_exception_dispatch():
    "Make sure that exception raised during dispatch are indeed captured"
    nose.tools.assert_raises(
            ValueError,
            Parallel(n_jobs=6, pre_dispatch=16, verbose=0),
                    (delayed(exception_raiser)(i) for i in range(30)),
            )


def _reload_joblib():
    # Retrieve the path of the parallel module in a robust way
    joblib_path = Parallel.__module__.split(os.sep)
    joblib_path = joblib_path[:1]
    joblib_path.append('parallel.py')
    joblib_path = '/'.join(joblib_path)
    module = __import__(joblib_path)
    # Reload the module. This should trigger a fail
    reload(module)


def test_multiple_spawning():
    # Test that attempting to launch a new Python after spawned
    # subprocesses will raise an error, to avoid infinite loops on
    # systems that do not support fork
    if not int(os.environ.get('JOBLIB_MULTIPROCESSING', 1)):
        raise nose.SkipTest()
    nose.tools.assert_raises(ImportError, Parallel(n_jobs=2),
                    [delayed(_reload_joblib)() for i in range(10)])


###############################################################################
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
