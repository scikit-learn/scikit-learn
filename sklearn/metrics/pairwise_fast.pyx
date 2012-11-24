# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# License: Simplified BSD

import numpy as np
cimport numpy as np
import cython

from ..utils.validation import array2d

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def chi2_kernel(X, Y=None):
    """Computes the chi-squared kernel between observations in X and Y.

    The chi-squared kernel is computed between each pair of rows in X and Y.  X
    and Y have to be non-negative. This kernel is most commonly applied to
    histograms.

    The chi-squared kernel is given by::

        k(x, y) = \sum_i (x[i] - y[i]) ** 2 / (x[i] + y[i])

    It can be interpreted as a weighted difference per entry.

    Parameters
    ----------
    X : array-like of shape (n_samples_X, n_features)

    Y : array of shape (n_samples_Y, n_features)

    Returns
    -------
    kernel_matrix : array of shape (n_samples_X, n_samples_Y)

    See also
    --------
    exponential_chi2_kernel : An exponentiated version of this kernel.

    kernel_approximation.AdditiveChi2Sampler : A Fourier approximation to this
        kernel.
    """
    ### once we support sparse matrices, we can use check_pairwise

    if Y is None:
        # optimize this case!
        Y = X
    if X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices: "
                         "X.shape[1] == %d while Y.shape[1] == %d" % (
                             X.shape[1], Y.shape[1]))

    X = array2d(X, dtype=np.float)
    Y = array2d(Y, dtype=np.float)

    if (X < 0).any():
        raise ValueError("X contains negative values.")
    if (Y < 0).any():
        raise ValueError("Y contains negative values.")

    result = np.zeros((X.shape[0], Y.shape[0]))

    cdef double[:, :] result_memview = result
    cdef double[:, :] X_memview = X
    cdef double[:, :] Y_memview = Y
    cdef int i, j, k
    cdef int n_samples_X = X.shape[0]
    cdef int n_samples_Y = Y.shape[0]
    cdef int n_features = X.shape[1]
    cdef double res, nom, denom
    for i in xrange(n_samples_X):
        for j in xrange(n_samples_Y):
            res = 0
            for k in xrange(n_features):
                denom = (X_memview[i, k] - Y_memview[j, k])
                nom = (X_memview[i, k] + Y_memview[j, k] + 1e-10)
                res  += denom * denom / nom
            result_memview[i, j] = res
    return result
