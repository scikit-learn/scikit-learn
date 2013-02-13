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
import io
import collections

from ..hashing import hash
from ..func_inspect import filter_args
from ..memory import Memory
from .common import np, with_numpy

from .test_memory import env as test_memory_env
from .test_memory import setup_module as test_memory_setup_func
from .test_memory import teardown_module as test_memory_teardown_func

try:
    # Python 2/Python 3 compat
    unicode('str')
except NameError:
    unicode = lambda s: s


###############################################################################
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
    relative_diff = 0.5 * (abs(time_func1 - time_func2)
                           / (time_func1 + time_func2))
    return relative_diff


class Klass(object):

    def f(self, x):
        return x


class KlassWithCachedMethod(object):

    def __init__(self):
        mem = Memory(cachedir=test_memory_env['dir'])
        self.f = mem.cache(self.f)

    def f(self, x):
        return x


###############################################################################
# Tests

def test_trival_hash():
    """ Smoke test hash on various types.
    """
    obj_list = [1, 2, 1., 2., 1 + 1j, 2. + 1j,
                'a', 'b',
                (1, ), (1, 1, ), [1, ], [1, 1, ],
                {1: 1}, {1: 2}, {2: 1},
                None,
                gc.collect,
                [1, ].append,
               ]
    for obj1 in obj_list:
        for obj2 in obj_list:
            # Check that 2 objects have the same hash only if they are
            # the same.
            yield nose.tools.assert_equal, hash(obj1) == hash(obj2), \
                obj1 is obj2


def test_hash_methods():
    # Check that hashing instance methods works
    a = io.StringIO(unicode('a'))
    nose.tools.assert_equal(hash(a.flush), hash(a.flush))
    a1 = collections.deque(range(10))
    a2 = collections.deque(range(9))
    nose.tools.assert_not_equal(hash(a1.extend), hash(a2.extend))


@with_numpy
def test_hash_numpy():
    """ Test hashing with numpy arrays.
    """
    rnd = np.random.RandomState(0)
    arr1 = rnd.random_sample((10, 10))
    arr2 = arr1.copy()
    arr3 = arr2.copy()
    arr3[0] += 1
    obj_list = (arr1, arr2, arr3)
    for obj1 in obj_list:
        for obj2 in obj_list:
            yield nose.tools.assert_equal, hash(obj1) == hash(obj2), \
                np.all(obj1 == obj2)

    d1 = {1: arr1, 2: arr1}
    d2 = {1: arr2, 2: arr2}
    yield nose.tools.assert_equal, hash(d1), hash(d2)

    d3 = {1: arr2, 2: arr3}
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
            except OSError as e:
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
    rnd = np.random.RandomState(0)
    a = rnd.random_sample(1000000)
    if hasattr(np, 'getbuffer'):
        # Under python 3, there is no getbuffer
        getbuffer = np.getbuffer
    else:
        getbuffer = memoryview
    md5_hash = lambda x: hashlib.md5(getbuffer(x)).hexdigest()

    relative_diff = relative_time(md5_hash, hash, a)
    yield nose.tools.assert_true, relative_diff < 0.1

    # Check that hashing an tuple of 3 arrays takes approximately
    # 3 times as much as hashing one array
    time_hashlib = 3 * time_func(md5_hash, a)
    time_hash = time_func(hash, (a, a, a))
    relative_diff = 0.5 * (abs(time_hash - time_hashlib)
                           / (time_hash + time_hashlib))

    yield nose.tools.assert_true, relative_diff < 0.2


def test_bound_methods_hash():
    """ Make sure that calling the same method on two different instances
    of the same class does resolve to the same hashes.
    """
    a = Klass()
    b = Klass()
    nose.tools.assert_equal(hash(filter_args(a.f, [], (1, ))),
                            hash(filter_args(b.f, [], (1, ))))


@nose.tools.with_setup(test_memory_setup_func, test_memory_teardown_func)
def test_bound_cached_methods_hash():
    """ Make sure that calling the same _cached_ method on two different
    instances of the same class does resolve to the same hashes.
    """
    a = KlassWithCachedMethod()
    b = KlassWithCachedMethod()
    nose.tools.assert_equal(hash(filter_args(a.f.func, [], (1, ))),
                            hash(filter_args(b.f.func, [], (1, ))))


@with_numpy
def test_hash_object_dtype():
    """ Make sure that ndarrays with dtype `object' hash correctly."""

    a = np.array([np.arange(i) for i in range(6)], dtype=object)
    b = np.array([np.arange(i) for i in range(6)], dtype=object)

    nose.tools.assert_equal(hash(a),
                            hash(b))


@with_numpy
def test_numpy_scalar():
    # Numpy scalars are built from compiled functions, and lead to
    # strange pickling paths explored, that can give hash collisions
    a = np.float64(2.0)
    b = np.float64(3.0)
    nose.tools.assert_not_equal(hash(a), hash(b))


def test_dict_hash():
    # Check that dictionaries hash consistently, eventhough the ordering
    # of the keys is not garanteed
    k = KlassWithCachedMethod()

    d = {'#s12069__c_maps.nii.gz': [33],
         '#s12158__c_maps.nii.gz': [33],
         '#s12258__c_maps.nii.gz': [33],
         '#s12277__c_maps.nii.gz': [33],
         '#s12300__c_maps.nii.gz': [33],
         '#s12401__c_maps.nii.gz': [33],
         '#s12430__c_maps.nii.gz': [33],
         '#s13817__c_maps.nii.gz': [33],
         '#s13903__c_maps.nii.gz': [33],
         '#s13916__c_maps.nii.gz': [33],
         '#s13981__c_maps.nii.gz': [33],
         '#s13982__c_maps.nii.gz': [33],
         '#s13983__c_maps.nii.gz': [33]}

    a = k.f(d)
    b = k.f(a)

    nose.tools.assert_equal(hash(a),
                            hash(b))


def test_set_hash():
    # Check that sets hash consistently, eventhough their ordering
    # is not garanteed
    k = KlassWithCachedMethod()

    s = set(['#s12069__c_maps.nii.gz',
             '#s12158__c_maps.nii.gz',
             '#s12258__c_maps.nii.gz',
             '#s12277__c_maps.nii.gz',
             '#s12300__c_maps.nii.gz',
             '#s12401__c_maps.nii.gz',
             '#s12430__c_maps.nii.gz',
             '#s13817__c_maps.nii.gz',
             '#s13903__c_maps.nii.gz',
             '#s13916__c_maps.nii.gz',
             '#s13981__c_maps.nii.gz',
             '#s13982__c_maps.nii.gz',
             '#s13983__c_maps.nii.gz'])

    a = k.f(s)
    b = k.f(a)

    nose.tools.assert_equal(hash(a), hash(b))
