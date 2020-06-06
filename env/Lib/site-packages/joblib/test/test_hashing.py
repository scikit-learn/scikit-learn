"""
Test the hashing module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import time
import hashlib
import sys
import os
import gc
import io
import collections
import itertools
import pickle
import random
from decimal import Decimal
import pytest

from joblib.hashing import hash
from joblib.func_inspect import filter_args
from joblib.memory import Memory
from joblib.testing import raises, skipif, fixture, parametrize
from joblib.test.common import np, with_numpy


def unicode(s):
    return s


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

    def __init__(self, cachedir):
        mem = Memory(cachedir=cachedir)
        self.f = mem.cache(self.f)

    def f(self, x):
        return x


###############################################################################
# Tests

input_list = [1, 2, 1., 2., 1 + 1j, 2. + 1j,
              'a', 'b',
              (1,), (1, 1,), [1, ], [1, 1, ],
              {1: 1}, {1: 2}, {2: 1},
              None,
              gc.collect,
              [1, ].append,
              # Next 2 sets have unorderable elements in python 3.
              set(('a', 1)),
              set(('a', 1, ('a', 1))),
              # Next 2 dicts have unorderable type of keys in python 3.
              {'a': 1, 1: 2},
              {'a': 1, 1: 2, 'd': {'a': 1}}]


@parametrize('obj1', input_list)
@parametrize('obj2', input_list)
def test_trivial_hash(obj1, obj2):
    """Smoke test hash on various types."""
    # Check that 2 objects have the same hash only if they are the same.
    are_hashes_equal = hash(obj1) == hash(obj2)
    are_objs_identical = obj1 is obj2
    assert are_hashes_equal == are_objs_identical


def test_hash_methods():
    # Check that hashing instance methods works
    a = io.StringIO(unicode('a'))
    assert hash(a.flush) == hash(a.flush)
    a1 = collections.deque(range(10))
    a2 = collections.deque(range(9))
    assert hash(a1.extend) != hash(a2.extend)


@fixture(scope='function')
@with_numpy
def three_np_arrays():
    rnd = np.random.RandomState(0)
    arr1 = rnd.random_sample((10, 10))
    arr2 = arr1.copy()
    arr3 = arr2.copy()
    arr3[0] += 1
    return arr1, arr2, arr3


def test_hash_numpy_arrays(three_np_arrays):
    arr1, arr2, arr3 = three_np_arrays

    for obj1, obj2 in itertools.product(three_np_arrays, repeat=2):
        are_hashes_equal = hash(obj1) == hash(obj2)
        are_arrays_equal = np.all(obj1 == obj2)
        assert are_hashes_equal == are_arrays_equal

    assert hash(arr1) != hash(arr1.T)


def test_hash_numpy_dict_of_arrays(three_np_arrays):
    arr1, arr2, arr3 = three_np_arrays

    d1 = {1: arr1, 2: arr2}
    d2 = {1: arr2, 2: arr1}
    d3 = {1: arr2, 2: arr3}

    assert hash(d1) == hash(d2)
    assert hash(d1) != hash(d3)


@with_numpy
@parametrize('dtype', ['datetime64[s]', 'timedelta64[D]'])
def test_numpy_datetime_array(dtype):
    # memoryview is not supported for some dtypes e.g. datetime64
    # see https://github.com/joblib/joblib/issues/188 for more details
    a_hash = hash(np.arange(10))
    array = np.arange(0, 10, dtype=dtype)
    assert hash(array) != a_hash


@with_numpy
def test_hash_numpy_noncontiguous():
    a = np.asarray(np.arange(6000).reshape((1000, 2, 3)),
                   order='F')[:, :1, :]
    b = np.ascontiguousarray(a)
    assert hash(a) != hash(b)

    c = np.asfortranarray(a)
    assert hash(a) != hash(c)


@with_numpy
@parametrize('coerce_mmap', [True, False])
def test_hash_memmap(tmpdir, coerce_mmap):
    """Check that memmap and arrays hash identically if coerce_mmap is True."""
    filename = tmpdir.join('memmap_temp').strpath
    try:
        m = np.memmap(filename, shape=(10, 10), mode='w+')
        a = np.asarray(m)
        are_hashes_equal = (hash(a, coerce_mmap=coerce_mmap) ==
                            hash(m, coerce_mmap=coerce_mmap))
        assert are_hashes_equal == coerce_mmap
    finally:
        if 'm' in locals():
            del m
            # Force a garbage-collection cycle, to be certain that the
            # object is delete, and we don't run in a problem under
            # Windows with a file handle still open.
            gc.collect()


@with_numpy
@skipif(sys.platform == 'win32', reason='This test is not stable under windows'
                                        ' for some reason')
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

    def md5_hash(x):
        return hashlib.md5(memoryview(x)).hexdigest()

    relative_diff = relative_time(md5_hash, hash, a)
    assert relative_diff < 0.3

    # Check that hashing an tuple of 3 arrays takes approximately
    # 3 times as much as hashing one array
    time_hashlib = 3 * time_func(md5_hash, a)
    time_hash = time_func(hash, (a, a, a))
    relative_diff = 0.5 * (abs(time_hash - time_hashlib)
                           / (time_hash + time_hashlib))
    assert relative_diff < 0.3


def test_bound_methods_hash():
    """ Make sure that calling the same method on two different instances
    of the same class does resolve to the same hashes.
    """
    a = Klass()
    b = Klass()
    assert (hash(filter_args(a.f, [], (1, ))) ==
            hash(filter_args(b.f, [], (1, ))))


def test_bound_cached_methods_hash(tmpdir):
    """ Make sure that calling the same _cached_ method on two different
    instances of the same class does resolve to the same hashes.
    """
    a = KlassWithCachedMethod(tmpdir.strpath)
    b = KlassWithCachedMethod(tmpdir.strpath)
    assert (hash(filter_args(a.f.func, [], (1, ))) ==
            hash(filter_args(b.f.func, [], (1, ))))


@with_numpy
def test_hash_object_dtype():
    """ Make sure that ndarrays with dtype `object' hash correctly."""

    a = np.array([np.arange(i) for i in range(6)], dtype=object)
    b = np.array([np.arange(i) for i in range(6)], dtype=object)

    assert hash(a) == hash(b)


@with_numpy
def test_numpy_scalar():
    # Numpy scalars are built from compiled functions, and lead to
    # strange pickling paths explored, that can give hash collisions
    a = np.float64(2.0)
    b = np.float64(3.0)
    assert hash(a) != hash(b)


def test_dict_hash(tmpdir):
    # Check that dictionaries hash consistently, eventhough the ordering
    # of the keys is not garanteed
    k = KlassWithCachedMethod(tmpdir.strpath)

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

    assert hash(a) == hash(b)


def test_set_hash(tmpdir):
    # Check that sets hash consistently, even though their ordering
    # is not guaranteed
    k = KlassWithCachedMethod(tmpdir.strpath)

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

    assert hash(a) == hash(b)


def test_set_decimal_hash():
    # Check that sets containing decimals hash consistently, even though
    # ordering is not guaranteed
    assert (hash(set([Decimal(0), Decimal('NaN')])) ==
            hash(set([Decimal('NaN'), Decimal(0)])))


def test_string():
    # Test that we obtain the same hash for object owning several strings,
    # whatever the past of these strings (which are immutable in Python)
    string = 'foo'
    a = {string: 'bar'}
    b = {string: 'bar'}
    c = pickle.loads(pickle.dumps(b))
    assert hash([a, b]) == hash([a, c])


@with_numpy
def test_dtype():
    # Test that we obtain the same hash for object owning several dtype,
    # whatever the past of these dtypes. Catter for cache invalidation with
    # complex dtype
    a = np.dtype([('f1', np.uint), ('f2', np.int32)])
    b = a
    c = pickle.loads(pickle.dumps(a))
    assert hash([a, c]) == hash([a, b])


@parametrize('to_hash,expected',
             [('This is a string to hash',
               '71b3f47df22cb19431d85d92d0b230b2'),
              (u"C'est l\xe9t\xe9",
               '2d8d189e9b2b0b2e384d93c868c0e576'),
              ((123456, 54321, -98765),
               'e205227dd82250871fa25aa0ec690aa3'),
              ([random.Random(42).random() for _ in range(5)],
               'a11ffad81f9682a7d901e6edc3d16c84'),
              ({'abcde': 123, 'sadfas': [-9999, 2, 3]},
                  'aeda150553d4bb5c69f0e69d51b0e2ef')])
def test_hashes_stay_the_same(to_hash, expected):
    # We want to make sure that hashes don't change with joblib
    # version. For end users, that would mean that they have to
    # regenerate their cache from scratch, which potentially means
    # lengthy recomputations.
    # Expected results have been generated with joblib 0.9.2
    assert hash(to_hash) == expected


@with_numpy
def test_hashes_are_different_between_c_and_fortran_contiguous_arrays():
    # We want to be sure that the c-contiguous and f-contiguous versions of the
    # same array produce 2 different hashes.
    rng = np.random.RandomState(0)
    arr_c = rng.random_sample((10, 10))
    arr_f = np.asfortranarray(arr_c)
    assert hash(arr_c) != hash(arr_f)


@with_numpy
def test_0d_array():
    hash(np.array(0))


@with_numpy
def test_0d_and_1d_array_hashing_is_different():
    assert hash(np.array(0)) != hash(np.array([0]))


@with_numpy
def test_hashes_stay_the_same_with_numpy_objects():
    # We want to make sure that hashes don't change with joblib
    # version. For end users, that would mean that they have to
    # regenerate their cache from scratch, which potentially means
    # lengthy recomputations.
    rng = np.random.RandomState(42)
    # Being explicit about dtypes in order to avoid
    # architecture-related differences. Also using 'f4' rather than
    # 'f8' for float arrays because 'f8' arrays generated by
    # rng.random.randn don't seem to be bit-identical on 32bit and
    # 64bit machines.
    to_hash_list = [
        rng.randint(-1000, high=1000, size=50).astype('<i8'),
        tuple(rng.randn(3).astype('<f4') for _ in range(5)),
        [rng.randn(3).astype('<f4') for _ in range(5)],
        {
            -3333: rng.randn(3, 5).astype('<f4'),
            0: [
                rng.randint(10, size=20).astype('<i8'),
                rng.randn(10).astype('<f4')
            ]
        },
        # Non regression cases for https://github.com/joblib/joblib/issues/308.
        # Generated with joblib 0.9.4.
        np.arange(100, dtype='<i8').reshape((10, 10)),
        # Fortran contiguous array
        np.asfortranarray(np.arange(100, dtype='<i8').reshape((10, 10))),
        # Non contiguous array
        np.arange(100, dtype='<i8').reshape((10, 10))[:, :2],
    ]

    # These expected results have been generated with joblib 0.9.0
    expected_hashes = [
        '10a6afc379ca2708acfbaef0ab676eab',
        '988a7114f337f381393025911ebc823b',
        'c6809f4b97e35f2fa0ee8d653cbd025c',
        'b3ad17348e32728a7eb9cda1e7ede438',
        '927b3e6b0b6a037e8e035bda134e0b05',
        '108f6ee98e7db19ea2006ffd208f4bf1',
        'bd48ccaaff28e16e6badee81041b7180'
    ]

    for to_hash, expected in zip(to_hash_list, expected_hashes):
        assert hash(to_hash) == expected


def test_hashing_pickling_error():
    def non_picklable():
        return 42

    with raises(pickle.PicklingError) as excinfo:
        hash(non_picklable)
    excinfo.match('PicklingError while hashing')


def test_wrong_hash_name():
    msg = "Valid options for 'hash_name' are"
    with raises(ValueError, match=msg):
        data = {'foo': 'bar'}
        hash(data, hash_name='invalid')
