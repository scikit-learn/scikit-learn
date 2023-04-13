# Author: Andrew nystrom <awnystrom@gmail.com>

from scipy.sparse import csr_matrix
cimport numpy as cnp
import numpy as np

cnp.import_array()

# TODO: use `cnp.{int,float}{32,64}` when cython#5230 is resolved:
# https://github.com/cython/cython/issues/5230
ctypedef fused DATA_T:
    float
    double
    int
    long


cdef inline cnp.int32_t _deg2_column(
    cnp.int32_t d,
    cnp.int32_t i,
    cnp.int32_t j,
    cnp.int32_t interaction_only,
) noexcept nogil:
    """Compute the index of the column for a degree 2 expansion

    d is the dimensionality of the input data, i and j are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return d * i - (i**2 + 3 * i) / 2 - 1 + j
    else:
        return d * i - (i**2 + i) / 2 + j


cdef inline cnp.int32_t _deg3_column(
    cnp.int32_t d,
    cnp.int32_t i,
    cnp.int32_t j,
    cnp.int32_t k,
    cnp.int32_t interaction_only
) noexcept nogil:
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


def _csr_polynomial_expansion(
    const DATA_T[:] data,
    const cnp.int32_t[:] indices,
    const cnp.int32_t[:] indptr,
    cnp.int32_t d,
    cnp.int32_t interaction_only,
    cnp.int32_t degree
):
    """
    Perform a second-degree polynomial or interaction expansion on a scipy
    compressed sparse row (CSR) matrix. The method used only takes products of
    non-zero features. For a matrix with density d, this results in a speedup
    on the order of d^k where k is the degree of the expansion, assuming all
    rows are of similar density.

    Parameters
    ----------
    data : memory view on nd-array
        The "data" attribute of the input CSR matrix.

    indices : memory view on nd-array
        The "indices" attribute of the input CSR matrix.

    indptr : memory view on nd-array
        The "indptr" attribute of the input CSR matrix.

    d : int
        The dimensionality of the input CSR matrix.

    interaction_only : int
        0 for a polynomial expansion, 1 for an interaction expansion.

    degree : int
        The degree of the expansion. This must be either 2 or 3.

    References
    ----------
    "Leveraging Sparsity to Speed Up Polynomial Feature Expansions of CSR
    Matrices Using K-Simplex Numbers" by Andrew Nystrom and John Hughes.
    """

    assert degree in (2, 3)

    if degree == 2:
        expanded_dimensionality = int((d**2 + d) / 2 - interaction_only*d)
    else:
        expanded_dimensionality = int((d**3 + 3*d**2 + 2*d) / 6
                                      - interaction_only*d**2)
    if expanded_dimensionality == 0:
        return None
    assert expanded_dimensionality > 0

    cdef cnp.int32_t total_nnz = 0, row_i, nnz

    # Count how many nonzero elements the expanded matrix will contain.
    for row_i in range(indptr.shape[0]-1):
        # nnz is the number of nonzero elements in this row.
        nnz = indptr[row_i + 1] - indptr[row_i]
        if degree == 2:
            total_nnz += (nnz ** 2 + nnz) / 2 - interaction_only * nnz
        else:
            total_nnz += ((nnz ** 3 + 3 * nnz ** 2 + 2 * nnz) / 6
                          - interaction_only * nnz ** 2)

    # Make the arrays that will form the CSR matrix of the expansion.
    cdef:
        DATA_T[:] expanded_data = np.empty(
            shape=total_nnz, dtype=data.base.dtype
        )
        cnp.int32_t[:] expanded_indices = np.empty(
            shape=total_nnz, dtype=np.int32
        )
        cnp.int32_t num_rows = indptr.shape[0] - 1
        cnp.int32_t[:] expanded_indptr = np.empty(
            shape=num_rows + 1, dtype=np.int32
        )

        cnp.int32_t expanded_index = 0, row_starts, row_ends, i, j, k, \
                i_ptr, j_ptr, k_ptr, num_cols_in_row,  \
                expanded_column

    with nogil:
        expanded_indptr[0] = indptr[0]
        for row_i in range(indptr.shape[0]-1):
            row_starts = indptr[row_i]
            row_ends = indptr[row_i + 1]
            num_cols_in_row = 0
            for i_ptr in range(row_starts, row_ends):
                i = indices[i_ptr]
                for j_ptr in range(i_ptr + interaction_only, row_ends):
                    j = indices[j_ptr]
                    if degree == 2:
                        col = _deg2_column(d, i, j, interaction_only)
                        expanded_indices[expanded_index] = col
                        expanded_data[expanded_index] = (
                            data[i_ptr] * data[j_ptr])
                        expanded_index += 1
                        num_cols_in_row += 1
                    else:
                        # degree == 3
                        for k_ptr in range(j_ptr + interaction_only,
                                            row_ends):
                            k = indices[k_ptr]
                            col = _deg3_column(d, i, j, k, interaction_only)
                            expanded_indices[expanded_index] = col
                            expanded_data[expanded_index] = (
                                data[i_ptr] * data[j_ptr] * data[k_ptr])
                            expanded_index += 1
                            num_cols_in_row += 1

            expanded_indptr[row_i+1] = expanded_indptr[row_i] + num_cols_in_row

    return csr_matrix((expanded_data, expanded_indices, expanded_indptr),
                      shape=(num_rows, expanded_dimensionality))
