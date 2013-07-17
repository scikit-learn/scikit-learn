# Authors: Mathieu Blondel
#          Olivier Grisel
#          Lars Buitinck
#
# Licence: BSD 3 clause

from libc.math cimport fabs, sqrt
cimport numpy as np
import numpy as np
import scipy.sparse as sp
cimport cython

np.import_array()


ctypedef np.float64_t DOUBLE

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csr_min_max_axis0(X):
    """Compute minimum and maximum along axis 0 on a CSR matrix

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    min: float array with shape (n_features,)
        Feature-wise minima

    max: float array with shape (n_features,)
        Feature-wise maxima
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int ind
    cdef double v

    # mins[j] contains the minimum of feature j
    cdef np.ndarray[DOUBLE, ndim=1] mins = np.empty((n_features))

    # maxs[j] contains the maximum of feature j
    cdef np.ndarray[DOUBLE, ndim=1] maxs = np.empty((n_features))

    # nvals[j] contains the number of values for feature j
    cdef np.ndarray[np.int64_t, ndim=1] nvals = np.zeros((n_features), dtype=np.int64)

    mins[:] = np.inf
    maxs[:] = -np.inf

    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            ind = X_indices[j]
            v   = X_data[j]
            mins[ind] = min(mins[ind], v)
            maxs[ind] = max(maxs[ind], v)
            nvals[ind] += 1

    # cython cannot do boolean arrays (yet)
    # cdef np.ndarray[bool, ndim=1] not_dense = nvals != n_samples
    not_dense = nvals != n_samples
    mins[not_dense * (mins > 0.0)] = 0.0
    maxs[not_dense * (maxs < 0.0)] = 0.0
    return mins, maxs


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csr_mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR matrix

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int ind
    cdef double diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] means = np.asarray(X.mean(axis=0))[0]

    # variances[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    # counts[j] contains the number of samples where feature j is non-zero
    cdef np.ndarray[np.int64_t, ndim=1] counts = np.zeros_like(means, dtype=np.int64)

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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_row_normalize_l1(X):
    """Inplace row normalize using the l1 norm"""
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
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += fabs(X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_row_normalize_l2(X):
    """Inplace row normalize using the l2 norm"""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        sum_ = sqrt(sum_)

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_column_scale(X, np.ndarray[DOUBLE, ndim=1] scale):
    """Inplace column scaling of a CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i, j
    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] *= scale[X_indices[j]]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csc_min_max_axis0(X):
    """Compute minimum and maximum along axis 0 on a CSC matrix

    Parameters
    ----------
    X: CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    min: float array with shape (n_features,)
        Feature-wise minima

    max: float array with shape (n_features,)
        Feature-wise maxima
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double v

    # mins[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] mins = np.empty((n_features))

    # maxs[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] maxs = np.empty((n_features))

    # nvals[j] contains the number of values for feature j
    cdef np.ndarray[np.int64_t, ndim=1] nvals = np.zeros((n_features), dtype=np.int64)

    for i in xrange(n_features):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            v = X_data[j]
            mins[i] = min(mins[i], v)
            maxs[i] = max(maxs[i], v)
            nvals[i] += 1

    # cython cannot do boolean arrays (yet)
    # cdef np.ndarray[bool, ndim=1] not_dense = nvals != n_samples
    not_dense = nvals != n_samples
    mins[not_dense * (mins > 0.0)] = 0.0
    maxs[not_dense * (maxs < 0.0)] = 0.0
    return mins, maxs


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csc_mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSC matrix

    Parameters
    ----------
    X: CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] means = np.asarray(X.mean(axis=0))[0]

    # variances[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    # counts[j] contains the number of samples where feature j is non-zero
    cdef np.ndarray[np.int64_t, ndim=1] counts = np.zeros_like(means, dtype=np.int64)

    for i in xrange(n_features):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            diff = X_data[j] - means[i]
            variances[i] += diff * diff
            counts[i] += 1

    nz = n_samples - counts
    variances += nz * means ** 2
    variances /= n_samples

    return means, variances


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csc_column_scale(X, np.ndarray[DOUBLE, ndim=1] scale):
    """Inplace column scaling of a CSC matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSC matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i, j
    for i in xrange(n_features):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] *= scale[i]


def min_max_axis0(X):
    """Compute minimum and maximum along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    mins: float array with shape (n_features,)
        Feature-wise minima

    maxs: float array with shape (n_features,)
        Feature-wise maxima

    """
    if isinstance(X, sp.csr_matrix):
        return csr_min_max_axis0(X)
    elif isinstance(X, sp.csc_matrix):
        return csc_min_max_axis0(X)
    else:
        raise TypeError(
                "Unsupported type; expected a CSR or CSC sparse matrix.")

def mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    if isinstance(X, sp.csr_matrix):
        return csr_mean_variance_axis0(X)
    elif isinstance(X, sp.csc_matrix):
        return csc_mean_variance_axis0(X)
    else:
        raise TypeError(
                "Unsupported type; expected a CSR or CSC sparse matrix.")


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void add_row_csr(np.ndarray[np.float64_t, ndim=1] data,
                      np.ndarray[int, ndim=1] indices,
                      np.ndarray[int, ndim=1] indptr,
                      int i, np.ndarray[np.float64_t, ndim=1, mode="c"] out):
    """Add row i of CSR matrix (data, indices, indptr) to array out.

    Equivalent to out += X[i].toarray(). Returns None.
    """
    cdef int ind, j

    for ind in range(indptr[i], indptr[i + 1]):
        j = indices[ind]
        out[j] += data[ind]


@cython.boundscheck(False)
@cython.wraparound(False)
def assign_rows_csr(X,
                    np.ndarray[np.npy_intp, ndim=1] X_rows,
                    np.ndarray[np.npy_intp, ndim=1] out_rows,
                    np.ndarray[np.float64_t, ndim=2, mode="c"] out):
    """Densify selected rows of a CSR matrix into a preallocated array.

    Like out[out_rows] = X[X_rows].toarray() but without copying. Only supported for
    dtype=np.float64.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix, shape=(n_samples, n_features)
    X_rows : array, dtype=np.intp, shape=n_rows
    out_rows : array, dtype=np.intp, shape=n_rows
    out : array, shape=(arbitrary, n_features)
    """
    cdef:
        # npy_intp (np.intc in Python) is what np.where returns,
        # but int is what scipy.sparse uses.
        int i, ind, j
        np.npy_intp rX
        np.ndarray[DOUBLE, ndim=1] data = X.data
        np.ndarray[int, ndim=1] indices = X.indices, indptr = X.indptr

    if X_rows.shape[0] != out_rows.shape[0]:
        raise ValueError("cannot assign %d rows to %d"
                         % (X_rows.shape[0], out_rows.shape[0]))

    out[:] = 0.
    for i in range(X_rows.shape[0]):
        # XXX we could reuse add_row_csr here, but the array slice
        # is not optimized away.
        rX = X_rows[i]
        for ind in range(indptr[rX], indptr[rX + 1]):
            j = indices[ind]
            out[out_rows[i], j] = data[ind]
