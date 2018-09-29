# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Andrew nystrom <awnystrom@gmail.com>

# When testing, test the following cases:
# A single dense row
# An empty matrix
# Multiple rows, but one is all zeros.
# The first row is all zeros. Tests we're making indptr correctly.

from scipy.sparse import csr_matrix
from numpy cimport ndarray
cimport numpy as np

ctypedef np.int32_t INDEX_T

ctypedef fused DATA_T:
    np.float32_t
    np.float64_t
    np.int32_t
    np.int64_t


def csr_expansion(ndarray[DATA_T, ndim=1] data,
                  ndarray[INDEX_T, ndim=1] indices,
                  ndarray[INDEX_T, ndim=1] indptr,
                  INDEX_T D, INDEX_T interaction_only, INDEX_T degree):
    """
    Perform a second-degree polynomial or interaction expansion on a scipy
    compressed sparse row (CSR) matrix. The method used only takes products of
    non-zero features. For a matrix with density d, this results in a speedup
    on the order of d^k where k is the degree of the expansion, assuming all
    rows are of similar density. This method was laid out in the work
    "Leveraging Sparsity to Speed Up Polynomial Feature Expansions of CSR
    Matrices Using K-Simplex Numbers" by Andrew Nystrom and John Hughes.

    Parameters
    ----------
    data : the "data" attribute of the input CSR matrix.
    indices : the "indices" attribute of the input CSR matrix.
    indptr : the "indptr" attribute of the input CSR matrix.
    D : The dimensionality of the input CSR matrix.
    interaction_only : 0 for a polynomial expansion, 1 for an interaction
                       expansion.
    degree : The degree of the expansion. This must be either 2 or 3.
    """

    assert degree in (2, 3)

    cdef INDEX_T total_nnz = 0, row_i, R

    # Count how many nonzero elements the expanded matrix will contain.
    for row_i in range(indptr.shape[0]-1):
        # R is the number of nonzero elements in this row.
        R = indptr[row_i + 1] - indptr[row_i]
        total_nnz += (R**3 + 3*R**2 + 2*R) / 6 - interaction_only*R**2 \
                     if degree == 3 else                               \
                     (R**2 + R) / 2 - interaction_only*R

    # Make the arrays that will form the CSR matrix of the expansion.
    cdef ndarray[DATA_T, ndim=1] expanded_data = ndarray(shape=total_nnz,
                                                         dtype=data.dtype,
                                                         order='C')
    cdef ndarray[INDEX_T, ndim=1] expanded_indices = ndarray(shape=total_nnz,
                                                             dtype='int32',
                                                             order='C')
    cdef INDEX_T num_rows = indptr.shape[0]
    cdef ndarray[INDEX_T, ndim=1] expanded_indptr = ndarray(shape=num_rows,
                                                              dtype='int32',
                                                              order='C')

    cdef INDEX_T expanded_index = 0, row_starts, row_ends, i, j, k, \
                 col_ptr_1, col_ptr_2, col_ptr_3, num_cols_in_row,  \
                 expanded_column, expanded_dimensionality

    expanded_dimensionality = (D**3 + 3*D**2 + 2*D)/6            \
                               - interaction_only*D**2           \
                              if degree == 3 else                \
                              (D**2 + D) / 2 - interaction_only*D

    expanded_indptr[0] = indptr[0]
    for row_i in range(indptr.shape[0]-1):
        row_starts = indptr[row_i]
        row_ends = indptr[row_i + 1]
        num_cols_in_row = 0
        for col_ptr_1 in range(row_starts, row_ends):
            i = indices[col_ptr_1]
            for col_ptr_2 in range(col_ptr_1 + interaction_only, row_ends):
                j = indices[col_ptr_2]
                if degree == 2:
                    expanded_column = D*i - (i**2 + 3*i) / 2 - 1 + j \
                                      if interaction_only else       \
                                      D*i - (i**2 + i) / 2 + j
                    expanded_indices[expanded_index] = expanded_column
                    expanded_data[expanded_index] = data[col_ptr_1] * \
                                                    data[col_ptr_2]
                    expanded_index += 1
                    num_cols_in_row += 1
                else:
                    # degree == 3
                    for col_ptr_3 in range(col_ptr_2 + interaction_only,
                                           row_ends):
                        k = indices[col_ptr_3]
                        expanded_column = (3*D**2*i - 3*D*i**2 + i**3 + 11*i \
                                           - 3*j**2 - 9*j) / 6 + i**2        \
                                           - 2*D*i + D*j - D + k             \
                                           if interaction_only else          \
                                           (3*D**2*i - 3*D*i**2 + i**3 - i   \
                                            - 3*j**2 - 3*j) / 6 + D*j + k
                        expanded_indices[expanded_index] = expanded_column
                        expanded_data[expanded_index] = data[col_ptr_1] * \
                                                        data[col_ptr_2] * \
                                                        data[col_ptr_3]
                        expanded_index += 1
                        num_cols_in_row += 1

        expanded_indptr[row_i+1] = expanded_indptr[row_i] + num_cols_in_row

    X_expanded = csr_matrix([])
    X_expanded.data = expanded_data
    X_expanded.indices = expanded_indices
    X_expanded.indptr = expanded_indptr
    X_expanded._shape = (indptr.shape[0]-1, expanded_dimensionality)
    return X_expanded
