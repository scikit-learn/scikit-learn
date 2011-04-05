"""
Test the hashing module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import nose
import time
import hashlib
import tempfile
import os
import gc

from ..hashing import hash
from .common import np, with_numpy

################################################################################
# Helper functions for the tests 
def time_func(func, *args):
    """ Time function func on *args.
    """
    times = list()
    for _ in range(3):
        t1 = time.time()
        func(*args)
        times.append(time.time() - t1)
    return min(times)


def relative_time(func1, func2, *args):
    """ Return the relative time between func1 and func2 applied on
        *args.
    """
    time_func1 = time_func(func1, *args)
    time_func2 = time_func(func2, *args)
    relative_diff = 0.5*( abs(time_func1 - time_func2)
                          /  (time_func1 + time_func2) )
    return relative_diff


################################################################################
# Tests

def test_trival_hash():
    """ Smoke test hash on various types.
    """
    obj_list = [1, 1., 1+1j,
                'a', 
                (1, ), [1, ], {1:1},
                None,
               ]
    for obj1 in obj_list:
        for obj2 in obj_list:
            yield nose.tools.assert_equal, hash(obj1) == hash(obj2), \
                obj1 is obj2



@with_numpy
def test_hash_numpy():
    """ Test hashing with numpy arrays.
    """
    arr1 = np.random.random((10, 10))
    arr2 = arr1.copy()
    arr3 = arr2.copy()
    arr3[0] += 1
    obj_list = (arr1, arr2, arr3)
    for obj1 in obj_list:
        for obj2 in obj_list:
            yield nose.tools.assert_equal, hash(obj1) == hash(obj2), \
                np.all(obj1 == obj2)

    d1 = {1:arr1, 2:arr1}
    d2 = {1:arr2, 2:arr2}
    yield nose.tools.assert_equal, hash(d1), hash(d2)

    d3 = {1:arr2, 2:arr3}
    yield nose.tools.assert_not_equal, hash(d1), hash(d3)

    yield nose.tools.assert_not_equal, hash(arr1), hash(arr1.T)


@with_numpy
def test_hash_memmap():
    """ Check that memmap and arrays hash identically if coerce_mmap is
        True.
    """
    filename = tempfile.mktemp()
    try:
        m = np.memmap(filename, shape=(10, 10), mode='w+')
        a = np.asarray(m)
        for coerce_mmap in (False, True):
            yield (nose.tools.assert_equal,
                            hash(a, coerce_mmap=coerce_mmap) 
                                == hash(m, coerce_mmap=coerce_mmap),
                            coerce_mmap)
    finally:
        if 'm' in locals():
            del m
            # Force a garbage-collection cycle, to be certain that the
            # object is delete, and we don't run in a problem under 
            # Windows with a file handle still open.
            gc.collect()
            try:
                os.unlink(filename)
            except OSError, e:
                # Under windows, some files don't get erased.
                if not os.name == 'nt':
                    raise e

@with_numpy
def test_hash_numpy_performance():
    """ Check the performance of hashing numpy arrays:

        In [22]: a = np.random.random(1000000)

        In [23]: %timeit hashlib.md5(a).hexdigest()
        100 loops, best of 3: 20.7 ms per loop

        In [24]: %timeit hashlib.md5(pickle.dumps(a, protocol=2)).hexdigest()
        1 loops, best of 3: 73.1 ms per loop

        In [25]: %timeit hashlib.md5(cPickle.dumps(a, protocol=2)).hexdigest()
        10 loops, best of 3: 53.9 ms per loop

        In [26]: %timeit hash(a)
        100 loops, best of 3: 20.8 ms per loop
    """
    a = np.random.random(1000000)
    md5_hash = lambda x: hashlib.md5(np.getbuffer(x)).hexdigest()

    relative_diff = relative_time(md5_hash, hash, a)
    yield nose.tools.assert_true, relative_diff < 0.1

    # Check that hashing an tuple of 3 arrays takes approximately
    # 3 times as much as hashing one array
    time_hashlib = 3*time_func(md5_hash, a)
    time_hash = time_func(hash, (a, a, a))
    relative_diff = 0.5*( abs(time_hash - time_hashlib)
                          /  (time_hash + time_hashlib) )

    yield nose.tools.assert_true, relative_diff < 0.2

