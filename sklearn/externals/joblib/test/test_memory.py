"""
Test the memory module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import shutil
import os
import os.path
from tempfile import mkdtemp
import pickle
import warnings
import io
import sys
import time

import nose

from ..memory import Memory, MemorizedFunc, NotMemorizedFunc, MemorizedResult
from ..memory import NotMemorizedResult, _FUNCTION_HASHES
from .common import with_numpy, np


###############################################################################
# Module-level variables for the tests
def f(x, y=1):
    """ A module-level function for testing purposes.
    """
    return x ** 2 + y


###############################################################################
# Test fixtures
env = dict()


def setup_module():
    """ Test setup.
    """
    cachedir = mkdtemp()
    env['dir'] = cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    # Don't make the cachedir, Memory should be able to do that on the fly
    print(80 * '_')
    print('test_memory setup (%s)' % env['dir'])
    print(80 * '_')


def _rmtree_onerror(func, path, excinfo):
    print('!' * 79)
    print('os function failed: %r' % func)
    print('file to be removed: %s' % path)
    print('exception was: %r' % excinfo[1])
    print('!' * 79)


def teardown_module():
    """ Test teardown.
    """
    shutil.rmtree(env['dir'], False, _rmtree_onerror)
    print(80 * '_')
    print('test_memory teardown (%s)' % env['dir'])
    print(80 * '_')


###############################################################################
# Helper function for the tests
def check_identity_lazy(func, accumulator):
    """ Given a function and an accumulator (a list that grows every
        time the function is called), check that the function can be
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


###############################################################################
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
    for compress in (False, True):
     for mmap_mode in ('r', None):
        # We turn verbosity on to smoke test the verbosity code, however,
        # we capture it, as it is ugly
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

            memory = Memory(cachedir=env['dir'], verbose=10,
                            mmap_mode=mmap_mode, compress=compress)
            # First clear the cache directory, to check that our code can
            # handle that
            # NOTE: this line would raise an exception, as the database file is
            # still open; we ignore the error since we want to test what
            # happens if the directory disappears
            shutil.rmtree(env['dir'], ignore_errors=True)
            g = memory.cache(f)
            g(1)
            g.clear(warn=False)
            current_accumulator = len(accumulator)
            out = g(1)
        finally:
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr

        yield nose.tools.assert_equal, len(accumulator), \
                    current_accumulator + 1
        # Also, check that Memory.eval works similarly
        yield nose.tools.assert_equal, memory.eval(f, 1), out
        yield nose.tools.assert_equal, len(accumulator), \
                    current_accumulator + 1

    # Now do a smoke test with a function defined in __main__, as the name
    # mangling rules are more complex
    f.__module__ = '__main__'
    memory = Memory(cachedir=env['dir'], verbose=0)
    memory.cache(f)(1)


def test_no_memory():
    """ Test memory with cachedir=None: no memoize """
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
    # Check that multiple use of lambda will raise collisions
    memory = Memory(cachedir=env['dir'], verbose=0)
    # For isolation with other tests
    memory.clear()
    a = lambda x: x
    a = memory.cache(a)
    b = lambda x: x + 1
    b = memory.cache(b)

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        nose.tools.assert_equal(0, a(0))
        nose.tools.assert_equal(2, b(1))
        nose.tools.assert_equal(1, a(1))

    # In recent Python versions, we can retrieve the code of lambdas,
    # thus nothing is raised
    nose.tools.assert_equal(len(w), 4)


def test_memory_warning_collision_detection():
    # Check that collisions impossible to detect will raise appropriate
    # warnings.
    memory = Memory(cachedir=env['dir'], verbose=0)
    # For isolation with other tests
    memory.clear()
    a1 = eval('lambda x: x')
    a1 = memory.cache(a1)
    b1 = eval('lambda x: x+1')
    b1 = memory.cache(b1)

    if not hasattr(warnings, 'catch_warnings'):
        # catch_warnings is new in Python 2.6
        return

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        a1(1)
        b1(1)
        a1(0)

        yield nose.tools.assert_equal, len(w), 2
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

        rnd = np.random.RandomState(0)
        for i in range(3):
            a = rnd.random_sample((10, 10))
            for _ in range(3):
                yield nose.tools.assert_true, np.all(cached_n(a) == a)
                yield nose.tools.assert_equal, len(accumulator), i + 1


@with_numpy
def test_memory_numpy_check_mmap_mode():
    """Check that mmap_mode is respected even at the first call"""

    memory = Memory(cachedir=env['dir'], mmap_mode='r', verbose=0)
    memory.clear(warn=False)

    @memory.cache()
    def twice(a):
        return a * 2

    a = np.ones(3)

    b = twice(a)
    c = twice(a)

    nose.tools.assert_true(isinstance(c, np.memmap))
    nose.tools.assert_equal(c.mode, 'r')

    nose.tools.assert_true(isinstance(b, np.memmap))
    nose.tools.assert_equal(b.mode, 'r')


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


def test_partial_decoration():
    "Check cache may be called with kwargs before decorating"
    memory = Memory(cachedir=env['dir'], verbose=0)

    test_values = [
        (['x'], 100, 'r'),
        ([], 10, None),
    ]
    for ignore, verbose, mmap_mode in test_values:
        @memory.cache(ignore=ignore, verbose=verbose, mmap_mode=mmap_mode)
        def z(x):
            pass

        yield nose.tools.assert_equal, z.ignore, ignore
        yield nose.tools.assert_equal, z._verbose, verbose
        yield nose.tools.assert_equal, z.mmap_mode, mmap_mode


def test_func_dir():
    # Test the creation of the memory cache directory for the function.
    memory = Memory(cachedir=env['dir'], verbose=0)
    memory.clear()
    path = __name__.split('.')
    path.append('f')
    path = os.path.join(env['dir'], 'joblib', *path)

    g = memory.cache(f)
    # Test that the function directory is created on demand
    yield nose.tools.assert_equal, g._get_func_dir(), path
    yield nose.tools.assert_true, os.path.exists(path)

    # Test that the code is stored.
    # For the following test to be robust to previous execution, we clear
    # the in-memory store
    _FUNCTION_HASHES.clear()
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
    # Test the memorized functions can be pickled and restored.
    memory = Memory(cachedir=env['dir'], verbose=0)
    g = memory.cache(f)
    output = g(1)

    h = pickle.loads(pickle.dumps(g))

    output_dir, _ = g.get_output_dir(1)
    yield nose.tools.assert_equal, output, h.load_output(output_dir)
    memory2 = pickle.loads(pickle.dumps(memory))
    yield nose.tools.assert_equal, memory.cachedir, memory2.cachedir

    # Smoke test that pickling a memory with cachedir=None works
    memory = Memory(cachedir=None, verbose=0)
    pickle.loads(pickle.dumps(memory))
    g = memory.cache(f)
    gp = pickle.loads(pickle.dumps(g))
    gp(1)


def test_call_and_shelve():
    """Test MemorizedFunc outputting a reference to cache.
    """

    for func, Result in zip((MemorizedFunc(f, env['dir']),
                             NotMemorizedFunc(f),
                             Memory(cachedir=env['dir']).cache(f),
                             Memory(cachedir=None).cache(f),
                             ),
                            (MemorizedResult, NotMemorizedResult,
                             MemorizedResult, NotMemorizedResult)):
        nose.tools.assert_equal(func(2), 5)
        result = func.call_and_shelve(2)
        nose.tools.assert_true(isinstance(result, Result))
        nose.tools.assert_equal(result.get(), 5)

        result.clear()
        nose.tools.assert_raises(KeyError, result.get)
        result.clear()  # Do nothing if there is no cache.


def test_memorized_pickling():
    for func in (MemorizedFunc(f, env['dir']), NotMemorizedFunc(f)):
        filename = os.path.join(env['dir'], 'pickling_test.dat')
        result = func.call_and_shelve(2)
        with open(filename, 'wb') as fp:
            pickle.dump(result, fp)
        with open(filename, 'rb') as fp:
            result2 = pickle.load(fp)
        nose.tools.assert_equal(result2.get(), result.get())
        os.remove(filename)


def test_memorized_repr():
    func = MemorizedFunc(f, env['dir'])
    result = func.call_and_shelve(2)

    func2 = MemorizedFunc(f, env['dir'])
    result2 = func2.call_and_shelve(2)
    nose.tools.assert_equal(result.get(), result2.get())
    nose.tools.assert_equal(repr(func), repr(func2))

    # Smoke test on deprecated methods
    func.format_signature(2)
    func.format_call(2)

    # Smoke test with NotMemorizedFunc
    func = NotMemorizedFunc(f)
    repr(func)
    repr(func.call_and_shelve(2))

    # Smoke test for message output (increase code coverage)
    func = MemorizedFunc(f, env['dir'], verbose=11, timestamp=time.time())
    result = func.call_and_shelve(11)
    result.get()

    func = MemorizedFunc(f, env['dir'], verbose=11)
    result = func.call_and_shelve(11)
    result.get()

    func = MemorizedFunc(f, env['dir'], verbose=5, timestamp=time.time())
    result = func.call_and_shelve(11)
    result.get()

    func = MemorizedFunc(f, env['dir'], verbose=5)
    result = func.call_and_shelve(11)
    result.get()


def test_memory_file_modification():
    # Test that modifying a Python file after loading it does not lead to
    # Recomputation
    dir_name = os.path.join(env['dir'], 'tmp_import')
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    filename = os.path.join(dir_name, 'tmp_joblib_.py')
    content = 'def f(x):\n    print(x)\n    return x\n'
    with open(filename, 'w') as module_file:
        module_file.write(content)

    # Load the module:
    sys.path.append(dir_name)
    import tmp_joblib_ as tmp

    mem = Memory(cachedir=env['dir'], verbose=0)
    f = mem.cache(tmp.f)
    # Capture sys.stdout to count how many time f is called
    orig_stdout = sys.stdout
    if sys.version_info[0] == 3:
        my_stdout = io.StringIO()
    else:
        my_stdout = io.BytesIO()

    try:
        sys.stdout = my_stdout

        # First call f a few times
        f(1)
        f(2)
        f(1)

        # Now modify the module where f is stored without modifying f
        with open(filename, 'w') as module_file:
            module_file.write('\n\n' + content)

        # And call f a couple more times
        f(1)
        f(1)

        # Flush the .pyc files
        shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        # Now modify the module where f is stored, modifying f
        content = 'def f(x):\n    print("x=%s" % x)\n    return x\n'
        with open(filename, 'w') as module_file:
            module_file.write(content)

        # And call f more times prior to reloading: the cache should not be
        # invalidated at this point as the active function definition has not
        # changed in memory yet.
        f(1)
        f(1)

        # Now reload
        my_stdout.write('Reloading\n')
        sys.modules.pop('tmp_joblib_')
        import tmp_joblib_ as tmp
        f = mem.cache(tmp.f)

        # And call f more times
        f(1)
        f(1)

    finally:
        sys.stdout = orig_stdout
    nose.tools.assert_equal(my_stdout.getvalue(), '1\n2\nReloading\nx=1\n')


def _function_to_cache(a, b):
    # Just a place holder function to be mutated by tests
    pass


def _sum(a, b):
    return a + b


def _product(a, b):
    return a * b


def test_memory_in_memory_function_code_change():
    _function_to_cache.__code__ = _sum.__code__

    mem = Memory(cachedir=env['dir'], verbose=0)
    f = mem.cache(_function_to_cache)

    nose.tools.assert_equal(f(1, 2), 3)
    nose.tools.assert_equal(f(1, 2), 3)

    with warnings.catch_warnings(record=True):
        # ignore name collision warnings
        warnings.simplefilter("always")

        # Check that inline function modification triggers a cache invalidation

        _function_to_cache.__code__ = _product.__code__
        nose.tools.assert_equal(f(1, 2), 2)
        nose.tools.assert_equal(f(1, 2), 2)
