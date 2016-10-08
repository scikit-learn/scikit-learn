# Author: Lars Buitinck
# License: BSD 3 clause

import array
from cpython cimport array
cimport cython
from libc.stdlib cimport abs
cimport numpy as np
import numpy as np

from ..externals.six import string_types

from sklearn.utils.murmurhash cimport murmurhash3_bytes_s32

np.import_array()


@cython.boundscheck(False)
@cython.cdivision(True)
def transform(raw_X, Py_ssize_t n_features, dtype):
    """Guts of FeatureHasher.transform.

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
    indptr = array.array("i", [0])

    # Since Python array does not understand Numpy dtypes, we grow the indices
    # and values arrays ourselves. Use a Py_ssize_t capacity for safety.
    cdef Py_ssize_t capacity = 8192     # arbitrary
    cdef np.int32_t size = 0
    cdef np.ndarray values = np.empty(capacity, dtype=dtype)

    for x in raw_X:
        for f, v in x:
            if isinstance(v, string_types):
                f = "%s%s%s" % (f, '=', v)
                value = 1
            else:
                value = v

            if value == 0:
                continue

            if isinstance(f, unicode):
                f = f.encode("utf-8")
            # Need explicit type check because Murmurhash does not propagate
            # all exceptions. Add "except *" there?
            elif not isinstance(f, bytes):
                raise TypeError("feature names must be strings")

            h = murmurhash3_bytes_s32(f, 0)

            array.resize_smart(indices, len(indices) + 1)
            indices[len(indices) - 1] = abs(h) % n_features
            value *= (h >= 0) * 2 - 1
            values[size] = value
            size += 1

            if size == capacity:
                capacity *= 2
                # can't use resize member because there might be multiple
                # references to the arrays due to Cython's error checking
                values = np.resize(values, capacity)

        array.resize_smart(indptr, len(indptr) + 1)
        indptr[len(indptr) - 1] = size

    if len(indices):
        indices_a = np.frombuffer(indices, dtype=np.int32)
    else:       # workaround for NumPy < 1.7.0
        indices_a = np.empty(0, dtype=np.int32)
    return (indices_a, np.frombuffer(indptr, dtype=np.int32), values[:size])
