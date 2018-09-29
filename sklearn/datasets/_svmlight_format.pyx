# Optimized inner loop of load_svmlight_file.
#
# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import array
from cpython cimport array
cimport cython
from libc.string cimport strchr

cimport numpy as np
import numpy as np
import scipy.sparse as sp

from ..externals.six import b

np.import_array()


cdef bytes COMMA = u','.encode('ascii')
cdef bytes COLON = u':'.encode('ascii')


@cython.boundscheck(False)
@cython.wraparound(False)
def _load_svmlight_file(f, dtype, bint multilabel, bint zero_based,
                        bint query_id, long long offset, long long length):
    cdef array.array data, indices, indptr
    cdef bytes line
    cdef char *hash_ptr
    cdef char *line_cstr
    cdef int idx, prev_idx
    cdef Py_ssize_t i
    cdef bytes qid_prefix = b('qid')
    cdef Py_ssize_t n_features
    cdef long long offset_max = offset + length if length > 0 else -1

    # Special-case float32 but use float64 for everything else;
    # the Python code will do further conversions.
    if dtype == np.float32:
        data = array.array("f")
    else:
        dtype = np.float64
        data = array.array("d")
    indices = array.array("i")
    indptr = array.array("i", [0])
    query = np.arange(0, dtype=np.int64)

    if multilabel:
        labels = []
    else:
        labels = array.array("d")

    if offset > 0:
        f.seek(offset)
        # drop the current line that might be truncated and is to be
        # fetched by another call
        f.readline()

    for line in f:
        # skip comments
        line_cstr = line
        hash_ptr = strchr(line_cstr, '#')
        if hash_ptr != NULL:
            line = line[:hash_ptr - line_cstr]

        line_parts = line.split()
        if len(line_parts) == 0:
            continue

        target, features = line_parts[0], line_parts[1:]
        if multilabel:
            if COLON in target:
                target, features = [], line_parts[0:]
            else:
                target = [float(y) for y in target.split(COMMA)]
            target.sort()
            labels.append(tuple(target))
        else:
            array.resize_smart(labels, len(labels) + 1)
            labels[len(labels) - 1] = float(target)

        prev_idx = -1
        n_features = len(features)
        if n_features and features[0].startswith(qid_prefix):
            _, value = features[0].split(COLON, 1)
            if query_id:
                query.resize(len(query) + 1)
                query[len(query) - 1] = np.int64(value)
            features.pop(0)
            n_features -= 1

        for i in xrange(0, n_features):
            idx_s, value = features[i].split(COLON, 1)
            idx = int(idx_s)
            if idx < 0 or not zero_based and idx == 0:
                raise ValueError(
                    "Invalid index %d in SVMlight/LibSVM data file." % idx)
            if idx <= prev_idx:
                raise ValueError("Feature indices in SVMlight/LibSVM data "
                                 "file should be sorted and unique.")

            array.resize_smart(indices, len(indices) + 1)
            indices[len(indices) - 1] = idx

            array.resize_smart(data, len(data) + 1)
            data[len(data) - 1] = float(value)

            prev_idx = idx

        array.resize_smart(indptr, len(indptr) + 1)
        indptr[len(indptr) - 1] = len(data)

        if offset_max != -1 and f.tell() > offset_max:
            # Stop here and let another call deal with the following.
            break

    return (dtype, data, indices, indptr, labels, query)
