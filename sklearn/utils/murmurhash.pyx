"""Cython wrapper for MurmurHash3 non-cryptographic hash function

MurmurHash is an extensively tested and very fast hash function that has
good distribution properties suitable for machine learning use cases
such as feature hashing and random projections.

The original C++ code by Austin Appleby is released the public domain
and can be found here:

  https://code.google.com/p/smhasher/

"""
# Author: Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD 3 clause
#
# cython: language_level=3


cimport cython
cimport numpy as np
import numpy as np

cdef extern from "src/MurmurHash3.h":
    void MurmurHash3_x86_32(void *key, int len, np.uint32_t seed, void *out)
    void MurmurHash3_x86_128(void *key, int len, np.uint32_t seed, void *out)
    void MurmurHash3_x64_128 (void *key, int len, np.uint32_t seed, void *out)

np.import_array()


cpdef np.uint32_t murmurhash3_int_u32(int key, unsigned int seed):
    """Compute the 32bit murmurhash3 of a int key at seed."""
    cdef np.uint32_t out
    MurmurHash3_x86_32(&key, sizeof(int), seed, &out)
    return out


cpdef np.int32_t murmurhash3_int_s32(int key, unsigned int seed):
    """Compute the 32bit murmurhash3 of a int key at seed."""
    cdef np.int32_t out
    MurmurHash3_x86_32(&key, sizeof(int), seed, &out)
    return out


cpdef np.uint32_t murmurhash3_bytes_u32(bytes key, unsigned int seed):
    """Compute the 32bit murmurhash3 of a bytes key at seed."""
    cdef np.uint32_t out
    MurmurHash3_x86_32(<char*> key, len(key), seed, &out)
    return out


cpdef np.int32_t murmurhash3_bytes_s32(bytes key, unsigned int seed):
    """Compute the 32bit murmurhash3 of a bytes key at seed."""
    cdef np.int32_t out
    MurmurHash3_x86_32(<char*> key, len(key), seed, &out)
    return out


@cython.boundscheck(False)
cpdef np.ndarray[np.uint32_t, ndim=1] murmurhash3_bytes_array_u32(
    np.ndarray[np.int32_t] key, unsigned int seed):
    """Compute 32bit murmurhash3 hashes of a key int array at seed."""
    # TODO make it possible to pass preallocated output array
    cdef np.ndarray[np.uint32_t, ndim=1] out = np.zeros(key.size, np.uint32)
    cdef Py_ssize_t i
    for i in range(key.shape[0]):
        out[i] = murmurhash3_int_u32(key[i], seed)
    return out


@cython.boundscheck(False)
cpdef np.ndarray[np.int32_t, ndim=1] murmurhash3_bytes_array_s32(
    np.ndarray[np.int32_t] key, unsigned int seed):
    """Compute 32bit murmurhash3 hashes of a key int array at seed."""
    # TODO make it possible to pass preallocated output array
    cdef np.ndarray[np.int32_t, ndim=1] out = np.zeros(key.size, np.int32)
    cdef Py_ssize_t i
    for i in range(key.shape[0]):
        out[i] = murmurhash3_int_s32(key[i], seed)
    return out


def murmurhash3_32(key, seed=0, positive=False):
    """Compute the 32bit murmurhash3 of key at seed.

    The underlying implementation is MurmurHash3_x86_32 generating low
    latency 32bits hash suitable for implementing lookup tables, Bloom
    filters, count min sketch or feature hashing.

    Parameters
    ----------
    key : int32, bytes, unicode or ndarray with dtype int32
        the physical object to hash

    seed : int, optional default is 0
        integer seed for the hashing algorithm.

    positive : boolean, optional default is False
        True: the results is casted to an unsigned int
          from 0 to 2 ** 32 - 1
        False: the results is casted to a signed int
          from -(2 ** 31) to 2 ** 31 - 1

    """
    if isinstance(key, bytes):
        if positive:
            return murmurhash3_bytes_u32(key, seed)
        else:
            return murmurhash3_bytes_s32(key, seed)
    elif isinstance(key, unicode):
        if positive:
            return murmurhash3_bytes_u32(key.encode('utf-8'), seed)
        else:
            return murmurhash3_bytes_s32(key.encode('utf-8'), seed)
    elif isinstance(key, int) or isinstance(key, np.int32):
        if positive:
            return murmurhash3_int_u32(<np.int32_t>key, seed)
        else:
            return murmurhash3_int_s32(<np.int32_t>key, seed)
    elif isinstance(key, np.ndarray):
        if key.dtype != np.int32:
            raise TypeError(
                "key.dtype should be int32, got %s" % key.dtype)
        if positive:
            return murmurhash3_bytes_array_u32(key.ravel(),
                                               seed).reshape(key.shape)
        else:
            return murmurhash3_bytes_array_s32(key.ravel(),
                                               seed).reshape(key.shape)
    else:
        raise TypeError(
            "key %r with type %s is not supported. "
            "Explicit conversion to bytes is required" % (key, type(key)))
