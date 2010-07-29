"""
Test the parallel module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2010 Gael Varoquaux
# License: BSD Style, 3 clauses.

try:
    import cPickle as pickle
    PickleError = TypeError
except:
    import pickle
    PickleError = pickle.PicklingError

from ..parallel import Parallel, delayed, JoblibException, SafeFunction

import nose

################################################################################

def division(x, y):
    return x/y

def square(x):
    return x**2

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
    """ Check that error are captured
    """
    nose.tools.assert_raises(JoblibException,
                                Parallel(n_jobs=2),
                    [delayed(division)(x, y) for x, y in zip((0, 1), (1, 0))],
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

