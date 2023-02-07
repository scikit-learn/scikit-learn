# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#         Lars Buitinck
#         Paolo Toccaceli
#
# License: BSD 3 clause

cimport numpy as cnp
from cython cimport floating

cnp.import_array()


def _chi2_kernel_fast(floating[:, :] X,
                      floating[:, :] Y,
                      floating[:, :] result):
    cdef cnp.npy_intp i, j, k
    cdef cnp.npy_intp n_samples_X = X.shape[0]
    cdef cnp.npy_intp n_samples_Y = Y.shape[0]
    cdef cnp.npy_intp n_features = X.shape[1]
    cdef double res, nom, denom

    with nogil:
        for i in range(n_samples_X):
            for j in range(n_samples_Y):
                res = 0
                for k in range(n_features):
                    denom = (X[i, k] - Y[j, k])
                    nom = (X[i, k] + Y[j, k])
                    if nom != 0:
                        res  += denom * denom / nom
                result[i, j] = -res
