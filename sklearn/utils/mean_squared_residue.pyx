# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

import numpy as np

cimport numpy as np
cimport cython

np.import_array()

ctypedef np.float64_t DOUBLE
ctypedef np.int64_t LONG


def mean_squared_residue(long[:] rows,
                         long[:] cols,
                         double[:] row_mean,
                         double[:] col_mean,
                         double bicluster_mean,
                         double[:, :] X):
    """Computes the mean squared reside (MSR) of a bicluster.

    The MSR, which was defined by Cheng and Church (2000), is a
    measure of the additive similarity of the row and column vectors
    of a matrix.

    In plain NumPy code, the array of residues ``R`` for an array
    ``X`` is calculated as:
        R = X - X.mean(axis=1, keepdims=True) - X.mean(axis=0) + X.mean()``

    It measures each element's coherence with the overall mean, row
    mean, and column mean. To get the overall mean squared residue:
        msr = (R * R).mean()

    Similarly, the row and column MSRs are:
        row_msr = (R * R).mean(axis=1)
        col_msr = (R * R).mean(axis=0)

    Arguments
    -----------
    rows : 1-D array of longs
        Rows in the bicluster.

    cols : 1-D array of longs
        Columns in the bicluster.

    row_mean : 1-D array of doubles
        ``row_mean[i]`` is the mean of bicluster row ``rows[i]``.

    col_mean : 1-D array of doubles
        ``col_mean[j]`` is the mean of bicluster column
        ``columns[j]``.

    bicluster_mean : double
        Mean of the entire bicluster.

    X : 2-D array of doubles
        The data in which the bicluster is defined.

    Returns
    -------
    msr : double
        The mean squared residue of the bicluster.

    row_msr : double
        ``row_msr[i]`` is the MSR of row ``rows[i]``.

    col_msr : double
        ``col_msr[j]`` is the MSR of column ``cols[j]``.

    References
    ----------

    - Cheng, Y., & Church, G. M. (2000). `Biclustering of
      expression data
      <ftp://samba.ad.sdsc.edu/pub/sdsc/biology/ISMB00/157.pdf>`__.

    """
    cdef long n_rows = rows.shape[0]
    cdef long n_cols = cols.shape[0]

    row_msr = np.zeros(n_rows, dtype=np.float, order="c")
    col_msr = np.zeros(n_cols, dtype=np.float, order="c")
    cdef double[:] row_msr_view = row_msr
    cdef double[:] col_msr_view = col_msr

    cdef double msr = 0.0
    cdef double val = 0.0

    cdef long i
    cdef long j

    with nogil:
        for i in range(n_rows):
            for j in range(n_cols):
                val = (X[rows[i], cols[j]] - row_mean[i] - col_mean[j] + arr_mean)
                val = val * val
                row_msr_view[i] += val
                col_msr_view[j] += val
                msr += val
        for i in range(n_rows):
            row_msr_view[i] /= n_cols
        for j in range(n_cols):
            col_msr_view[j] /= n_rows
        msr /= (n_rows * n_cols)

    return msr, row_msr, col_msr
