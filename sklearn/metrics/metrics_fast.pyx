# Author: Mathieu Blondel <mathieu@mblondel.org>
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
import cython

np.import_array()


cdef inline int _sign(double val):
    if val > 0:
        return 1
    elif val < 0:
        return -1
    else:
        return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _pairwise_ranking_accuracy(np.ndarray[double, ndim=1] y_true,
                               np.ndarray[double, ndim=1] y_score):
    cdef int n_samples = y_true.shape[0]

    cdef int n_correct = 0
    cdef int n_total = 0
    cdef int i, j

    # FIXME: the following algorithm can be used for reducing the complexity
    # from O(n^2) to O(n log(n)):
    # David Christensen.
    # Fast algorithms for the calculation of Kendall’s tau.
    # Computational Statistics, 20:51–62, 2005.

    for i in xrange(n_samples):
        for j in xrange(i + 1, n_samples):
            if y_true[i] == y_true[j]:
                continue

            n_total += 1
            if _sign(y_true[i] - y_true[j]) == _sign(y_score[i] - y_score[j]):
                n_correct += 1

    return n_correct, n_total
