# Author: Lars Buitinck
# License: BSD 3 clause
#
# cython: boundscheck=False, cdivision=True

import sys
import array
from cpython cimport array
cimport cython
from libc.stdlib cimport abs
cimport numpy as np
import numpy as np

from ..utils.murmurhash cimport murmurhash3_bytes_s32
from ..utils.fixes import sp_version

np.import_array()


def transform(raw_X, Py_ssize_t n_features, dtype, unsigned int seed=0):
    """Guts of FeatureLightweightRandomIndexing.transform.

    Returns
    -------
    n_samples : integer
    indices, indptr, values : lists
        For constructing a scipy.sparse.csr_matrix.

    """
    assert n_features > 0

    cdef np.int32_t h
    cdef double value

    cdef array.array indices
    cdef array.array indptr
    indices = array.array("i")
    indices_array_dtype = "q"
    indices_np_dtype = np.longlong


    indptr = array.array(indices_array_dtype, [0])

    # Since Python array does not understand Numpy dtypes, we grow the indices
    # and values arrays ourselves. Use a Py_ssize_t capacity for safety.
    cdef Py_ssize_t capacity = 8192     # arbitrary
    cdef np.int64_t size = 0
    cdef np.ndarray values = np.empty(capacity, dtype=dtype)

    for x in raw_X:
        for f, v in x:
            if isinstance(v, (str, unicode)):
                f = "%s%s%s" % (f, '=', v)
                value = 1
            else:
                value = v

            if value == 0:
                continue

            if isinstance(f, unicode):
                f = (<unicode>f).encode("utf-8")
            # Need explicit type check because Murmurhash does not propagate
            # all exceptions. Add "except *" there?
            elif not isinstance(f, bytes):
                raise TypeError("feature names must be strings")

            h0 = murmurhash3_bytes_s32(<bytes>f, seed)

            array.resize_smart(indices, len(indices) + 1)
            indices[len(indices) - 1] = abs(h0) % n_features
            # improve inner product preservation in the hashed space
            value0 = value * ((h0 >= 0) * 2 - 1)
            values[size] = value0
            size += 1

            if size == capacity:
                capacity *= 2
                # can't use resize member because there might be multiple
                # references to the arrays due to Cython's error checking
                values = np.resize(values, capacity)

            h0 = abs(h0) % n_features

            ih1 = h0
            offset = 1
            while ih1==h0:
                h1 = murmurhash3_bytes_s32(<bytes>f, seed+offset)
                ih1 = abs(h1) % n_features
                offset += 1

            array.resize_smart(indices, len(indices) + 1)
            indices[len(indices) - 1] = ih1
            # improve inner product preservation in the hashed space
            value1 = value * ((h1 >= 0) * 2 - 1)
            values[size] = value1
            size += 1

            if size == capacity:
                capacity *= 2
                # can't use resize member because there might be multiple
                # references to the arrays due to Cython's error checking
                values = np.resize(values, capacity)


        array.resize_smart(indptr, len(indptr) + 1)
        indptr[len(indptr) - 1] = size

    indices_a = np.frombuffer(indices, dtype=np.int32)
    indptr_a = np.frombuffer(indptr, dtype=indices_np_dtype)

    if indptr[-1] > np.iinfo(np.int32).max:  # = 2**31 - 1
        if sp_version < (0, 14):
            raise ValueError(('sparse CSR array has {} non-zero '
                              'elements and requires 64 bit indexing, '
                              ' which is unsupported with scipy {}. '
                              'Please upgrade to scipy >=0.14')
                             .format(indptr[-1], '.'.join(sp_version)))
        # both indices and indptr have the same dtype in CSR arrays
        indices_a = indices_a.astype(np.int64, copy=False)
    else:
        indptr_a = indptr_a.astype(np.int32, copy=False)

    return (indices_a, indptr_a, values[:size])
