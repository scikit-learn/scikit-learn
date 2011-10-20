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

def _load_svmlight_file(f, n_features, dtype):
    cdef bytes line
    cdef char *hash_ptr, *line_cstr
    cdef Py_ssize_t hash_idx

    data = []
    indptr = []
    indices = []
    labels = []

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
            indices.append(int(idx))
            data.append(dtype(value))

    indptr.append(len(data))
    indptr = np.array(indptr, dtype=np.int)

    if n_features is not None:
        shape = (indptr.shape[0] - 1, n_features)
    else:
        shape = None    # inferred

    X = sp.csr_matrix((np.array(data),
                       np.array(indices, dtype=np.int),
                       indptr), shape)

    return X, np.array(labels, dtype=np.double)
