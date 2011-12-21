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


def _load_svmlight_file(f, n_features, dtype):
    cdef bytes line
    cdef char *hash_ptr, *line_cstr
    cdef Py_ssize_t hash_idx

    data = ArrayBuilder(dtype=dtype)
    indptr = ArrayBuilder(dtype=_INDPTR_DTYPE)
    indices = ArrayBuilder(dtype=_INDICES_DTYPE)
    labels = ArrayBuilder(dtype=np.double)

    for line in f:
        # skip comments
        line_cstr = line
        hash_ptr = strchr(line_cstr, '#')
        if hash_ptr == NULL:
            hash_idx = -1           # index of '\n' in line
        else:
            hash_idx = hash_ptr - <char *>line
        line = line[:hash_idx]

        line_parts = line.split()
        if len(line_parts) == 0:
            continue

        y, features = line_parts[0], line_parts[1:]
        labels.append(float(y))
        indptr.append(len(data))

        for i in xrange(1, len(line_parts)):
            idx, value = line_parts[i].split(":", 1)
            # Real programmers count from zero.
            idx = int(idx)
            if idx <= 0:
                raise ValueError(
                        "invalid index %d in SVMlight/LibSVM data file" % idx)
            indices.append(idx - 1)
            data.append(dtype(value))

    indptr.append(len(data))
    indptr = indptr.get()

    if n_features is not None:
        shape = (indptr.shape[0] - 1, n_features)
    else:
        shape = None    # inferred

    X = sp.csr_matrix((np.array(data.get()),
                       np.array(indices.get(), dtype=np.int),
                       indptr), shape)

    return X, labels.get()
