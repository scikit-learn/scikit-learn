# Authors: Mathieu Blondel
#          Olivier Grisel
#          Peter Prettenhofer
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
def csr_row_norms(X):
    """L2 norm of each row in CSR matrix X."""
    cdef:
        unsigned int n_samples = X.shape[0]
        unsigned int n_features = X.shape[1]
        np.ndarray[DOUBLE, ndim=1, mode="c"] norms
        np.ndarray[DOUBLE, ndim=1, mode="c"] data
        np.ndarray[int, ndim=1, mode="c"] indices = X.indices
        np.ndarray[int, ndim=1, mode="c"] indptr = X.indptr

        np.npy_intp i, j
        double sum_

    norms = np.zeros(n_samples, dtype=np.float64)
    data = np.asarray(X.data, dtype=np.float64)     # might copy!

    for i in range(n_samples):
        sum_ = 0.0
        for j in range(indptr[i], indptr[i + 1]):
            sum_ += data[j] * data[j]
        norms[i] = sum_

    return norms


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

    cdef unsigned int i
    cdef unsigned int non_zero = X_indices.shape[0]
    cdef unsigned int col_ind
    cdef double diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] means = np.zeros(n_features)

    # variances[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    # counts[j] contains the number of samples where feature j is non-zero
    cdef np.ndarray[int, ndim=1] counts = np.zeros(n_features,
                                                   dtype=np.int32)

    for i in xrange(non_zero):
        col_ind = X_indices[i]
        means[col_ind] += X_data[i]
    
    means /= n_samples

    for i in xrange(non_zero):
        col_ind = X_indices[i]
        diff = X_data[i] - means[col_ind]
        variances[col_ind] += diff * diff
        counts[col_ind] += 1

    for i in xrange(n_features):
        variances[i] += (n_samples - counts[i]) * means[i] ** 2
        variances[i] /= n_samples

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
    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices

    cdef unsigned int i
    cdef unsigned non_zero = len(X_indices)

    for i in xrange(non_zero):
        X_data[i] *= scale[X_indices[i]]


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
    cdef unsigned int counts
    cdef unsigned int startptr
    cdef unsigned int endptr
    cdef double diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] means = np.zeros(n_features)

    # variances[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    for i in xrange(n_features):

        startptr = X_indptr[i]
        endptr = X_indptr[i + 1]
        counts = endptr - startptr

        for j in xrange(startptr, endptr):
            means[i] += X_data[j]
        means[i] /= n_samples

        for j in xrange(startptr, endptr):
            diff = X_data[j] - means[i]
            variances[i] += diff * diff

        variances[i] += (n_samples - counts) * means[i] * means[i]
        variances[i] /= n_samples

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


@cython.boundscheck(False)
@cython.wraparound(False)
def swap_row_csc(X, int m, int n):
    """
    Swaps two rows of a CSC matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first_sample
    n : int, index of second_sample
    """
    cdef np.ndarray[int, ndim=1] indices = X.indices
    cdef unsigned int i
    cdef unsigned nonzero = indices.shape[0]

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    for i in range(nonzero):
        if indices[i] == m:
            indices[i] = n
        elif indices[i] == n:
            indices[i] = m


@cython.boundscheck(False)
@cython.wraparound(False)
def swap_row_csr(X, int m, int n):
    """
    Swaps two rows of a CSC matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first_sample
    n : int, index of second_sample
    """
    cdef np.ndarray[int, ndim=1] indices = X.indices, indptr = X.indptr
    cdef np.ndarray[DOUBLE, ndim=1] data = X.data
    cdef unsigned int n_nz_rm, n_nz_rn, i, j1, j2, j3
    cdef unsigned int mptr, nptr

    # Some copies to retain the original information
    cdef np.ndarray[int, ndim=1] ptrcopy = np.copy(indptr)
    cdef np.ndarray[int, ndim=1] indcopy = np.copy(indices)
    cdef np.ndarray[DOUBLE, ndim=1] datacopy = np.copy(data)

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    if m > n:
        m, n = n, m

    # Modify indptr in the region in between the indices of two rows to
    # be swapped, if the number of non-zero values in both the rows are
    # not equal.
    n_nz_rm = indptr[m + 1] - indptr[m]
    n_nz_rn = indptr[n + 1] - indptr[n]

    if n_nz_rm != n_nz_rn:
        indptr[m + 1] = indptr[m] + n_nz_rn
        for i in range(m + 2, n + 1):
            indptr[i] = indptr[i - 1] + (ptrcopy[i] - ptrcopy[i - 1])

    # Modify indices and data in between the indices of two rows to be swapped
    m_ptr = ptrcopy[m]
    n_ptr = ptrcopy[n]
    if n_nz_rm == n_nz_rn:
        for i in range(n_nz_rm):
            indices[m_ptr + i], indices[n_ptr + i] = \
                indices[n_ptr + i], indices[m_ptr + i]
            data[m_ptr + i], data[n_ptr + i] = \
                data[n_ptr + i], data[m_ptr + i]

    else:
        j1 = 0
        j2 = ptrcopy[m + 1]
        j3 = 0
        for i in range(m_ptr, ptrcopy[n + 1]):
            if i < m_ptr + n_nz_rn:
                indices[i] = indcopy[n_ptr + j1]
                data[i] = datacopy[n_ptr + j1]
                j1 += 1
            elif i >= m_ptr + n_nz_rn and i < ptrcopy[n + 1] - n_nz_rm:
                indices[i] = indcopy[j2]
                data[i] = datacopy[j2]
                j2 += 1
            else:
                indices[i] = indcopy[m_ptr + j3]
                data[i] = datacopy[m_ptr + j3]
                j3 += 1
