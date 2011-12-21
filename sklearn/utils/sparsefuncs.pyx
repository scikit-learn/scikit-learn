# Author: Mathieu Blondel
#
# License: BSD Style.

cimport numpy as np
import numpy as np
cimport cython

ctypedef np.float64_t DOUBLE

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def mean_variance_axis0(X):
    """Compute mean and variance along axis 0

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    (means, variances):
        means: np.mean(X, axis=0)
        variances: np.var(X, axis=0)
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    # the column indices for row i are stored in:
    #    indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in:
    #    data[indptr[i]:indptr[i+1]]
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int ind
    cdef double diff

    # means[j] contains the mean of feature j
    means = np.asarray(X.mean(axis=0))[0]

    # counts[j] contains the number of samples where feature j is non-zero
    counts = np.zeros_like(means)

    # variances[j] contains the variance of feature j
    variances = np.zeros_like(means)

    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            ind = X_indices[j]
            diff = X_data[j] - means[ind]
            variances[ind] += diff * diff
            counts[ind] += 1

    nz = n_samples - counts
    variances += nz * means ** 2

    variances /= n_samples

    return means, variances
