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
# License: BSD Style.


cdef extern from "MurmurHash3.h":
    void MurmurHash3_x86_32(void* key, int len, unsigned int seed,
                            void* out)


cpdef unsigned int murmurhash3_32_int_uint(int key, unsigned int seed):
    cdef unsigned int out
    MurmurHash3_x86_32(&key, sizeof(int), seed, &out)
    return out


cpdef int murmurhash3_32_int_int(int key, unsigned int seed):
    cdef int out
    MurmurHash3_x86_32(&key, sizeof(int), seed, &out)
    return out


cpdef unsigned int murmurhash3_32_bytes_uint(bytes key, unsigned int seed):
    cdef unsigned int out
    cdef char* key_c = key
    MurmurHash3_x86_32(key_c, len(key), seed, &out)
    return out


cpdef int murmurhash3_32_bytes_int(bytes key, unsigned int seed):
    cdef int out
    cdef char* key_c = key
    MurmurHash3_x86_32(key_c, len(key), seed, &out)
    return out


def murmurhash3_32(key, seed=0, positive=False):
    """Compute the 32bit murmurhash3_32 of key at seed.

    The underlying implementation is MurmurHash3_x86_32 generating low
    latency 32bits hash suitable for implementing lookup tables, Bloom
    filters, count min sketch or feature hashing.

    Parameters
    ----------
    key: int, bytes or unicode
        the physical object to hash

    seed: int, optional default is 0
        integer seed for the hashing algorithm.

    positive: boolean, optional default is False
        True: the results is casted to an unsigned int
          from 0 to 2 ** 32 - 1
        False: the results is casted to a signed int
          from -(2 ** 31) to 2 ** 31 - 1

    """
    if isinstance(key, bytes):
        if positive:
            return murmurhash3_32_bytes_uint(key, seed)
        else:
            return murmurhash3_32_bytes_int(key, seed)
    elif isinstance(key, unicode):
        if positive:
            return murmurhash3_32_bytes_uint(key.encode('utf-8'), seed)
        else:
            return murmurhash3_32_bytes_int(key.encode('utf-8'), seed)
    elif isinstance(key, int):
        if positive:
            return murmurhash3_32_int_uint(key, seed)
        else:
            return murmurhash3_32_int_int(key, seed)
    else:
        raise ValueError(
            "key %r with type %s is not supported. "
            "Explicit conversion to bytes is required" % (key, type(key)))
