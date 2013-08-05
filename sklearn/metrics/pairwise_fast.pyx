# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
import cython

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _chi2_kernel_fast(cython.floating[:,:] X,
                      cython.floating[:,:] Y,
                      cython.floating[:, :] result):
    cdef int i, j, k
    cdef int n_samples_X = X.shape[0]
    cdef int n_samples_Y = Y.shape[0]
    cdef int n_features = X.shape[1]
    cdef double res, nom, denom
    for i in xrange(n_samples_X):
        for j in xrange(n_samples_Y):
            res = 0
            for k in xrange(n_features):
                denom = (X[i, k] - Y[j, k])
                nom = (X[i, k] + Y[j, k])
                if nom != 0:
                    res  += denom * denom / nom
            result[i, j] = -res
