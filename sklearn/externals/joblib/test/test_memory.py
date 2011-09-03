"""
Test the memory module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

from __future__ import with_statement

import shutil
import os
from tempfile import mkdtemp
import pickle
import warnings

import nose

from ..memory import Memory, MemorizedFunc
from .common import with_numpy, np

################################################################################
# Module-level variables for the tests
def f(x, y=1):
    """ A module-level function for testing purposes.
    """
    return x**2 + y


################################################################################
# Test fixtures
env = dict()

def setup_module():
    """ Test setup.
    """
    cachedir = mkdtemp()
    #cachedir = 'foobar'
    env['dir'] = cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    # Don't make the cachedir, Memory should be able to do that on the
    # fly
    print 80*'_'
    print 'test_memory setup'
    print 80*'_'
    
def _rmtree_onerror(func, path, excinfo):
    print '!'*79
    print 'os function failed:', repr(func)
    print 'file to be removed:', path
    print 'exception was:', excinfo[1]
    print '!'*79

def teardown_module():
    """ Test teardown.
    """
    shutil.rmtree(env['dir'], False, _rmtree_onerror)
    print 80*'_'
    print 'test_memory teardown'
    print 80*'_'


################################################################################
# Helper function for the tests
def check_identity_lazy(func, accumulator):
    """ Given a function and an accumulator (a list that grows every
        time the function is called, check that the function can be
        decorated by memory to be a lazy identity.
    """
    # Call each function with several arguments, and check that it is
    # evaluated only once per argument.
    memory = Memory(cachedir=env['dir'], verbose=0)
    memory.clear(warn=False)
    func = memory.cache(func)
    for i in range(3):
        for _ in range(2):
            yield nose.tools.assert_equal, func(i), i
            yield nose.tools.assert_equal, len(accumulator), i + 1


################################################################################
# Tests
def test_memory_integration():
    """ Simple test of memory lazy evaluation.
    """
    accumulator = list()
    # Rmk: this function has the same name than a module-level function,
    # thus it serves as a test to see that both are identified
    # as different.
    def f(l):
        accumulator.append(1)
        return l

    for test in check_identity_lazy(f, accumulator):
        yield test

    # Now test clearing
    memory = Memory(cachedir=env['dir'], verbose=0)
    # First clear the cache directory, to check that our code can
    # handle that
    # NOTE: this line would raise an exception, as the database file is still
    # open; we ignore the error since we want to test what happens if the
    # directory disappears
    shutil.rmtree(env['dir'], ignore_errors=True)
    g = memory.cache(f)
    g(1)
    g.clear(warn=False)
    current_accumulator = len(accumulator)
    out = g(1)
    yield nose.tools.assert_equal, len(accumulator), \
                current_accumulator + 1
    # Also, check that Memory.eval works similarly
    yield nose.tools.assert_equal, memory.eval(f, 1), out
    yield nose.tools.assert_equal, len(accumulator), \
                current_accumulator + 1


def test_no_memory():
    """ Test memory with cachedir=None: no memoize
    """
    accumulator = list()
    def ff(l):
        accumulator.append(1)
        return l
    mem = Memory(cachedir=None, verbose=0)
    gg = mem.cache(ff)
    for _ in range(4):
        current_accumulator = len(accumulator)
        gg(1)
        yield nose.tools.assert_equal, len(accumulator), \
                    current_accumulator + 1


def test_memory_kwarg():
    " Test memory with a function with keyword arguments."
    accumulator = list()
    def g(l=None, m=1):
        accumulator.append(1)
        return l

    for test in check_identity_lazy(g, accumulator):
        yield test

    memory = Memory(cachedir=env['dir'], verbose=0)
    g = memory.cache(g)
    # Smoke test with an explicit keyword argument:
    nose.tools.assert_equal(g(l=30, m=2), 30)


def test_memory_lambda():
    " Test memory with a function with a lambda."
    accumulator = list()
    def helper(x):
        """ A helper function to define l as a lambda.
        """
        accumulator.append(1)
        return x

    l = lambda x: helper(x)

    for test in check_identity_lazy(l, accumulator):
        yield test


def test_memory_name_collision():
    " Check that name collisions with functions will raise warnings"
    memory = Memory(cachedir=env['dir'], verbose=0)

    @memory.cache
    def name_collision(x):
        """ A first function called name_collision 
        """
        return x

    a = name_collision

    @memory.cache
    def name_collision(x):
        """ A second function called name_collision 
        """
        return x

    b = name_collision

    if not hasattr(warnings, 'catch_warnings'):
        # catch_warnings is new in Python 2.6
        return

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        a(1)
        b(1)

        yield nose.tools.assert_equal, len(w), 1
        yield nose.tools.assert_true, "collision" in str(w[-1].message)


def test_memory_warning_lambda_collisions():
    " Check that multiple use of lambda will raise collisions"
    memory = Memory(cachedir=env['dir'], verbose=0)
    a = lambda x: x
    a = memory.cache(a)
    b = lambda x: x+1
    b = memory.cache(b)

    if not hasattr(warnings, 'catch_warnings'):
        # catch_warnings is new in Python 2.6
        return

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        a(1)
        b(1)

        yield nose.tools.assert_equal, len(w), 2
        yield nose.tools.assert_true, "collision" in str(w[-1].message)
        yield nose.tools.assert_true, "collision" in str(w[-2].message)


def test_memory_warning_collision_detection():
    """ Check that collisions impossible to detect will raise appropriate 
        warnings.
    """
    memory = Memory(cachedir=env['dir'], verbose=0)
    a = eval('lambda x: x')
    a = memory.cache(a)
    b = eval('lambda x: x+1')
    b = memory.cache(b)

    if not hasattr(warnings, 'catch_warnings'):
        # catch_warnings is new in Python 2.6
        return

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        a(1)
        b(1)

        yield nose.tools.assert_equal, len(w), 1
        yield nose.tools.assert_true, \
                "cannot detect" in str(w[-1].message).lower()


def test_memory_partial():
    " Test memory with functools.partial."
    accumulator = list()
    def func(x, y):
        """ A helper function to define l as a lambda.
        """
        accumulator.append(1)
        return y

    import functools
    function = functools.partial(func, 1)

    for test in check_identity_lazy(function, accumulator):
        yield test


def test_memory_eval():
    " Smoke test memory with a function with a function defined in an eval."
    memory = Memory(cachedir=env['dir'], verbose=0)

    m = eval('lambda x: x')
    mm = memory.cache(m)

    yield nose.tools.assert_equal, 1, mm(1)


def count_and_append(x=[]):
    """ A function with a side effect in its arguments. 

        Return the lenght of its argument and append one element.
    """
    len_x = len(x)
    x.append(None)
    return len_x

def test_argument_change():
    """ Check that if a function has a side effect in its arguments, it
        should use the hash of changing arguments.
    """
    mem = Memory(cachedir=env['dir'], verbose=0)
    func = mem.cache(count_and_append)
    # call the function for the first time, is should cache it with
    # argument x=[]
    assert func() == 0
    # the second time the argument is x=[None], which is not cached
    # yet, so the functions should be called a second time
    assert func() == 1


@with_numpy
def test_memory_numpy():
    " Test memory with a function with numpy arrays."
    # Check with memmapping and without.
    for mmap_mode in (None, 'r'):
        accumulator = list()
        def n(l=None):
            accumulator.append(1)
            return l

        memory = Memory(cachedir=env['dir'], mmap_mode=mmap_mode,
                            verbose=0)
        memory.clear(warn=False)
        cached_n = memory.cache(n)
        for i in range(3):
            a = np.random.random((10, 10))
            for _ in range(3):
                yield nose.tools.assert_true, np.all(cached_n(a) == a)
                yield nose.tools.assert_equal, len(accumulator), i + 1


def test_memory_exception():
    """ Smoketest the exception handling of Memory. 
    """
    memory = Memory(cachedir=env['dir'], verbose=0)
    class MyException(Exception):
        pass

    @memory.cache
    def h(exc=0):
        if exc:
            raise MyException

    # Call once, to initialise the cache
    h()

    for _ in range(3):
        # Call 3 times, to be sure that the Exception is always raised
        yield nose.tools.assert_raises, MyException, h, 1


def test_memory_ignore():
    " Test the ignore feature of memory "
    memory = Memory(cachedir=env['dir'], verbose=0)
    accumulator = list()

    @memory.cache(ignore=['y'])
    def z(x, y=1):
        accumulator.append(1)

    yield nose.tools.assert_equal, z.ignore, ['y']

    z(0, y=1)
    yield nose.tools.assert_equal, len(accumulator), 1
    z(0, y=1)
    yield nose.tools.assert_equal, len(accumulator), 1
    z(0, y=2)
    yield nose.tools.assert_equal, len(accumulator), 1


def test_func_dir():
    """ Test the creation of the memory cache directory for the function.
    """
    memory = Memory(cachedir=env['dir'], verbose=0)
    path = __name__.split('.')
    path.append('f')
    path = os.path.join(env['dir'], 'joblib', *path)

    g = memory.cache(f)
    # Test that the function directory is created on demand
    yield nose.tools.assert_equal, g._get_func_dir(), path
    yield nose.tools.assert_true, os.path.exists(path)

    # Test that the code is stored.
    yield nose.tools.assert_false, \
        g._check_previous_func_code()
    yield nose.tools.assert_true, \
            os.path.exists(os.path.join(path, 'func_code.py'))
    yield nose.tools.assert_true, \
        g._check_previous_func_code()

    # Test the robustness to failure of loading previous results.
    dir, _ = g.get_output_dir(1)
    a = g(1)
    yield nose.tools.assert_true, os.path.exists(dir)
    os.remove(os.path.join(dir, 'output.pkl'))
    yield nose.tools.assert_equal, a, g(1)



def test_persistence():
    """ Test the memorized functions can be pickled and restored.
    """
    memory = Memory(cachedir=env['dir'], verbose=0)
    g = memory.cache(f)
    output = g(1)

    h = pickle.loads(pickle.dumps(g))

    output_dir, _ = g.get_output_dir(1)
    yield nose.tools.assert_equal, output, h.load_output(output_dir)


def test_format_signature():
    """ Test the signature formatting.
    """
    func = MemorizedFunc(f, cachedir=env['dir'])
    path, sgn = func.format_signature(f, range(10))
    yield nose.tools.assert_equal, \
                sgn, \
                'f([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])'
    path, sgn = func.format_signature(f, range(10), y=range(10))
    yield nose.tools.assert_equal, \
                sgn, \
        'f([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])'


@with_numpy
def test_format_signature_numpy():
    """ Test the format signature formatting with numpy.
    """

