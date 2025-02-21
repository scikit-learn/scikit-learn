import threading

import numpy as np
import pytest

from numpy.testing import IS_WASM
from numpy.testing._private.utils import run_threaded

if IS_WASM:
    pytest.skip(allow_module_level=True, reason="no threading support in wasm")


def test_parallel_randomstate_creation():
    # if the coercion cache is enabled and not thread-safe, creating
    # RandomState instances simultaneously leads to a data race
    def func(seed):
        np.random.RandomState(seed)

    run_threaded(func, 500, pass_count=True)


def test_parallel_ufunc_execution():
    # if the loop data cache or dispatch cache are not thread-safe
    # computing ufuncs simultaneously in multiple threads leads
    # to a data race that causes crashes or spurious exceptions
    def func():
        arr = np.random.random((25,))
        np.isnan(arr)

    run_threaded(func, 500)

    # see gh-26690
    NUM_THREADS = 50

    a = np.ones(1000)

    def f(b):
        b.wait()
        return a.sum()

    run_threaded(f, NUM_THREADS, max_workers=NUM_THREADS, pass_barrier=True)


def test_temp_elision_thread_safety():
    amid = np.ones(50000)
    bmid = np.ones(50000)
    alarge = np.ones(1000000)
    blarge = np.ones(1000000)

    def func(count):
        if count % 4 == 0:
            (amid * 2) + bmid
        elif count % 4 == 1:
            (amid + bmid) - 2
        elif count % 4 == 2:
            (alarge * 2) + blarge
        else:
            (alarge + blarge) - 2

    run_threaded(func, 100, pass_count=True)


def test_eigvalsh_thread_safety():
    # if lapack isn't thread safe this will randomly segfault or error
    # see gh-24512
    rng = np.random.RandomState(873699172)
    matrices = (
        rng.random((5, 10, 10, 3, 3)),
        rng.random((5, 10, 10, 3, 3)),
    )

    run_threaded(lambda i: np.linalg.eigvalsh(matrices[i]), 2,
                 pass_count=True)


def test_printoptions_thread_safety():
    # until NumPy 2.1 the printoptions state was stored in globals
    # this verifies that they are now stored in a context variable
    b = threading.Barrier(2)

    def legacy_113():
        np.set_printoptions(legacy='1.13', precision=12)
        b.wait()
        po = np.get_printoptions()
        assert po['legacy'] == '1.13'
        assert po['precision'] == 12
        orig_linewidth = po['linewidth']
        with np.printoptions(linewidth=34, legacy='1.21'):
            po = np.get_printoptions()
            assert po['legacy'] == '1.21'
            assert po['precision'] == 12
            assert po['linewidth'] == 34
        po = np.get_printoptions()
        assert po['linewidth'] == orig_linewidth
        assert po['legacy'] == '1.13'
        assert po['precision'] == 12

    def legacy_125():
        np.set_printoptions(legacy='1.25', precision=7)
        b.wait()
        po = np.get_printoptions()
        assert po['legacy'] == '1.25'
        assert po['precision'] == 7
        orig_linewidth = po['linewidth']
        with np.printoptions(linewidth=6, legacy='1.13'):
            po = np.get_printoptions()
            assert po['legacy'] == '1.13'
            assert po['precision'] == 7
            assert po['linewidth'] == 6
        po = np.get_printoptions()
        assert po['linewidth'] == orig_linewidth
        assert po['legacy'] == '1.25'
        assert po['precision'] == 7

    task1 = threading.Thread(target=legacy_113)
    task2 = threading.Thread(target=legacy_125)

    task1.start()
    task2.start()


def test_parallel_reduction():
    # gh-28041
    NUM_THREADS = 50

    x = np.arange(1000)

    def closure(b):
        b.wait()
        np.sum(x)

    run_threaded(closure, NUM_THREADS, max_workers=NUM_THREADS,
                 pass_barrier=True)


def test_parallel_flat_iterator():
    # gh-28042
    x = np.arange(20).reshape(5, 4).T

    def closure(b):
        b.wait()
        for _ in range(100):
            list(x.flat)

    run_threaded(closure, outer_iterations=100, pass_barrier=True)

    # gh-28143
    def prepare_args():
        return [np.arange(10)]

    def closure(x, b):
        b.wait()
        for _ in range(100):
            y = np.arange(10)
            y.flat[x] = x

    run_threaded(closure, pass_barrier=True, prepare_args=prepare_args)
