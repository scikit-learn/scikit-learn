# Optimized inner loop of load_svmlight_file.
#
# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: Simple BSD.

from libc.string cimport strchr
cimport numpy as np
import numpy as np
import scipy.sparse as sp

from ..utils.arraybuilder import ArrayBuilder


# csr_matrix.indices and .indptr's dtypes are undocumented. We derive them
# empirically.
_temp_csr = sp.csr_matrix(0)
_INDICES_DTYPE = _temp_csr.indices.dtype
_INDPTR_DTYPE = _temp_csr.indptr.dtype
del _temp_csr


cdef bytes COMMA = u','.encode('ascii')
cdef bytes COLON = u':'.encode('ascii')


def _load_svmlight_file(f, dtype, bint multilabel, bint zero_based):
    cdef bytes line
    cdef char *hash_ptr, *line_cstr
    cdef np.int32_t idx, prev_idx
    cdef Py_ssize_t i

    data = ArrayBuilder(dtype=dtype)
    indptr = ArrayBuilder(dtype=_INDPTR_DTYPE)
    indices = ArrayBuilder(dtype=_INDICES_DTYPE)
    if multilabel:
        labels = []
    else:
        labels = ArrayBuilder(dtype=np.double)

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
            target = [float(y) for y in target.split(COMMA)]
            target.sort()
            labels.append(tuple(target))
        else:
            labels.append(float(target))
        indptr.append(len(data))

        prev_idx = -1
        for i in xrange(1, len(line_parts)):
            idx_s, value = line_parts[i].split(COLON, 1)
            # XXX if we replace int with np.int32 in the line below, this
            # function becomes twice as slow.
            idx = int(idx_s)
            if idx < 0 or not zero_based and idx == 0:
                raise ValueError(
                        "Invalid index %d in SVMlight/LibSVM data file." % idx)
            if idx <= prev_idx:
                raise ValueError("Feature ndices in SVMlight/LibSVM data "
                                 "file should be sorted and unique.")
            indices.append(idx)
            data.append(dtype(value))
            prev_idx = idx

    indptr.append(len(data))

    indptr = indptr.get()
    data = data.get()
    indices = indices.get()

    if not multilabel:
        labels = labels.get()

    return data, indices, indptr, labels
