# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: BSD 3 clause

import array
from cpython cimport array
cimport cython
from libc.stdlib cimport abs
cimport numpy as np
import numpy as np

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

    cdef np.int32_t feature_idx
    cdef dict feature_counter

    cdef array.array indices
    cdef array.array indptr
    cdef array.array values
    indices = array.array("i")
    indptr = array.array("i", [0])
    values = array.array("d")

    for x in raw_X:
        feature_counter = {}
        for f, v in x:
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
            feature_idx = abs(h) % n_features 
            value *= (h >= 0) * 2 - 1
            if feature_idx not in feature_counter:
                feature_counter[feature_idx] = value
            else:
                feature_counter[feature_idx] += value

        indices.extend(feature_counter.keys())
        values.extend(feature_counter.values())
        del(feature_counter)

        indptr.append(len(indices))

    if len(indices):
        indices_a = np.frombuffer(indices, dtype=np.int32)
        values_a = np.frombuffer(values, dtype=np.float64)
        values_a = np.array(values, dtype=dtype)
    else:       # workaround for NumPy < 1.7.0
        indices_a = np.empty(0, dtype=np.int32)
        values_a = np.empty(0, dtype=dtype)
    return (indices_a, np.frombuffer(indptr, dtype=np.int32), values_a)
