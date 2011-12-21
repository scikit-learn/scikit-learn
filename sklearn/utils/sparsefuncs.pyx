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
    cdef unsigned int ptr
    cdef unsigned int ind
    cdef double diff

    # store the means in a 1d-array
    means = np.asarray(X.mean(axis=0))[0]

    variances = np.zeros_like(means)

    for i in xrange(n_samples):
        ptr = X_indptr[i]

        # we need to iterate over all features
        # but CSR matrices do not provide O(1)
        # access to features
        for j in xrange(n_features):
            ind = X_indices[ptr]

            if j != ind:
                diff = means[j]
                variances[j] += diff * diff
            else:
                diff = X_data[ptr] - means[ind]
                variances[ind] += diff * diff
                ptr += 1

    variances /= n_samples

    return means, variances
