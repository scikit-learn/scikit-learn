# Author: Andrew nystrom <awnystrom@gmail.com>

from scipy.sparse import csr_matrix
cimport numpy as cnp
cnp.import_array()

ctypedef cnp.int8_t FLAG_T
ctypedef cnp.int32_t INDEX_A

ctypedef fused INDEX_B:
    cnp.int32_t
    cnp.int64_t
ctypedef fused DATA_T:
    cnp.float32_t
    cnp.float64_t
    cnp.int32_t
    cnp.int64_t

cdef inline INDEX_B _deg2_column(INDEX_B d, INDEX_B i, INDEX_B j,
                                 FLAG_T interaction_only) nogil:
    """Compute the index of the column for a degree 2 expansion

    d is the dimensionality of the input data, i and j are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return d * i - (i**2 + 3 * i) / 2 - 1 + j
    else:
        return d * i - (i**2 + i) / 2 + j


cdef inline INDEX_B _deg3_column(INDEX_B d, INDEX_B i, INDEX_B j, INDEX_B k,
                                 FLAG_T interaction_only) nogil:
    """Compute the index of the column for a degree 3 expansion

    d is the dimensionality of the input data, i, j and k are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return ((3 * d**2 * i - 3 * d * i**2 + i**3
                 + 11 * i - 3 * j**2 - 9 * j) / 6
                + i**2 - 2 * d * i + d * j - d + k)
    else:
        return ((3 * d**2 * i - 3 * d * i**2 + i ** 3 - i
                 - 3 * j**2 - 3 * j) / 6
                + d * j + k)


def _csr_polynomial_expansion(cnp.ndarray[DATA_T, ndim=1] data,
                              cnp.ndarray[INDEX_A, ndim=1] indices,
                              cnp.ndarray[INDEX_A, ndim=1] indptr,
                              INDEX_A d,
                              cnp.ndarray[DATA_T, ndim=1] result_data,
                              cnp.ndarray[INDEX_B, ndim=1] result_indices,
                              cnp.ndarray[INDEX_B, ndim=1] result_indptr,
                              FLAG_T interaction_only,
                              FLAG_T degree):
    """
    Perform a second-degree polynomial or interaction expansion on a scipy
    compressed sparse row (CSR) matrix. The method used only takes products of
    non-zero features. For a matrix with density d, this results in a speedup
    on the order of d^k where k is the degree of the expansion, assuming all
    rows are of similar density.

    Parameters
    ----------
    data : nd-array
        The "data" attribute of the input CSR matrix.

    indices : nd-array
        The "indices" attribute of the input CSR matrix.

    indptr : nd-array
        The "indptr" attribute of the input CSR matrix.

    d : int
        The dimensionality of the input CSR matrix.

    result_data : nd-array
         The output CSR matrix's "data" attribute

    result_indices : nd-array
         The output CSR matrix's "indices" attribute

    result_indptr : nd-array
        The output CSR matrix's "indptr" attribute

    interaction_only : int
        0 for a polynomial expansion, 1 for an interaction expansion.

    degree : int
        The degree of the expansion. This must be either 2 or 3.

    References
    ----------
    "Leveraging Sparsity to Speed Up Polynomial Feature Expansions of CSR
    Matrices Using K-Simplex Numbers" by Andrew Nystrom and John Hughes.
    """

    # Make the arrays that will form the CSR matrix of the expansion.
    cdef INDEX_A row_i, row_starts, row_ends, i, j, k, i_ptr, j_ptr, k_ptr

    cdef INDEX_B expanded_index=0, num_cols_in_row, col

    with nogil:
        result_indptr[0] = indptr[0]
        for row_i in range(indptr.shape[0]-1):
            row_starts = indptr[row_i]
            row_ends = indptr[row_i + 1]
            num_cols_in_row = 0
            for i_ptr in range(row_starts, row_ends):
                i = indices[i_ptr]
                for j_ptr in range(i_ptr + interaction_only, row_ends):
                    j = indices[j_ptr]
                    if degree == 2:
                        col = _deg2_column[INDEX_B](d, i, j, interaction_only)
                        result_indices[expanded_index] = col
                        result_data[expanded_index] = (
                            data[i_ptr] * data[j_ptr])
                        expanded_index += 1
                        num_cols_in_row += 1
                    else:
                        # degree == 3
                        for k_ptr in range(j_ptr + interaction_only,
                                            row_ends):
                            k = indices[k_ptr]
                            col = _deg3_column[INDEX_B](d, i, j, k, interaction_only)
                            result_indices[expanded_index] = col
                            result_data[expanded_index] = (
                                data[i_ptr] * data[j_ptr] * data[k_ptr])
                            expanded_index += 1
                            num_cols_in_row += 1

            result_indptr[row_i+1] = result_indptr[row_i] + num_cols_in_row
