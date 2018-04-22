# Author: Lars Buitinck
# License: BSD 3 clause

import sys
import array
import numpy as np

from sklearn.utils.murmurhash import murmurhash3_bytes_s32
from sklearn.utils.fixes import sp_version
from sklearn.externals import six


def transform(raw_X, n_features, dtype, alternate_sign=1):
    """Guts of FeatureHasher.transform.

    Returns
    -------
    n_samples : integer
    indices, indptr, values : lists
        For constructing a scipy.sparse.csr_matrix.

    """
    assert n_features > 0

    indices = []
    if sys.version_info >= (3, 3):
        indices_array_dtype = "q"
        indices_np_dtype = np.longlong
    else:
        # On Windows with PY2.7 long int would still correspond to 32 bit.
        indices_array_dtype = "l"
        indices_np_dtype = np.int_

    indptr = [0]

    # Since Python array does not understand Numpy dtypes, we grow the indices
    # and values arrays ourselves. Use a Py_ssize_t capacity for safety.
    capacity = 8192     # arbitrary
    size = 0
    values = np.empty(capacity, dtype=dtype)

    for x in raw_X:
        for f, v in x:
            if isinstance(v, six.string_types):
                f = "%s%s%s" % (f, '=', v)
                value = 1
            else:
                value = v

            if value == 0:
                continue

            if isinstance(f, six.text_type):
                f = f.encode("utf-8")
            # Need explicit type check because Murmurhash does not propagate
            # all exceptions. Add "except *" there?
            elif not isinstance(f, six.binary_type):
                raise TypeError("feature names must be strings")

            h = murmurhash3_bytes_s32(f, 0)

            indices.append(abs(h) % n_features)
            # improve inner product preservation in the hashed space
            if alternate_sign:
                value *= (h >= 0) * 2 - 1
            values[size] = value
            size += 1

            if size == capacity:
                capacity *= 2
                # can't use resize member because there might be multiple
                # references to the arrays due to Cython's error checking
                values = np.resize(values, capacity)

        indptr.append(size)

    indices = array.array("i", indices)
    indptr = array.array(indices_array_dtype, indptr)
    indices_a = np.frombuffer(indices, dtype=np.int32)
    indptr_a = np.frombuffer(indptr, dtype=indices_np_dtype)

    if indptr[-1] > 2147483648:  # = 2**31
        if sp_version < (0, 14):
            raise ValueError(('sparse CSR array has {} non-zero '
                              'elements and requires 64 bit indexing, '
                              ' which is unsupported with scipy {}. '
                              'Please upgrade to scipy >=0.14')
                             .format(indptr[-1], '.'.join(sp_version)))
        # both indices and indptr have the same dtype in CSR arrays
        indices_a = indices_a.astype(np.int64)
    else:
        indptr_a = indptr_a.astype(np.int32)

    return (indices_a, indptr_a, values[:size])
