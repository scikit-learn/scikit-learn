# Author: Andrew nystrom <awnystrom@gmail.com>

from scipy.sparse import csr_matrix
from libc.math cimport sqrt, pow

cimport numpy as cnp
cnp.import_array()

cdef extern from "limits.h":
    cdef int INT_MAX
    cdef long LONG_MAX

ctypedef cnp.int8_t FLAG_t

# This is the maximum value before the intermediate computation
# i**2 + 3*i blows out.
cdef cnp.int64_t MAX_SAFE_INDEX_DEG2 = <cnp.int64_t> sqrt(LONG_MAX) - 4

# This is the maximum value before the intermediate computation
# 3 * d**2 * i + i**3 blows out, since d is the upper bound of i.
cdef cnp.int64_t MAX_SAFE_INDEX_DEG3 = <cnp.int64_t> pow(LONG_MAX, 1/3)/4

# INDEX_{A,B}_t are defined to generate a proper Cartesian product
# of types through Cython fused-type expansion.
ctypedef fused INDEX_A_t:
    cnp.int32_t
    cnp.int64_t
ctypedef fused INDEX_B_t:
    cnp.int32_t
    cnp.int64_t

ctypedef fused DATA_t:
    cnp.float32_t
    cnp.float64_t
    cnp.int32_t
    cnp.int64_t

cdef inline cnp.int64_t _deg2_column(
    cnp.int64_t d,
    cnp.int64_t i,
    cnp.int64_t j,
    FLAG_t interaction_only
) nogil:
    """Compute the index of the column for a degree 2 expansion

    d is the dimensionality of the input data, i and j are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return d * i - (i**2 + 3 * i) / 2 - 1 + j
    else:
        return d * i - (i**2 + i) / 2 + j

def py_deg2_column(d, i, j, interaction_only):
    """Compute the index of the column for a degree 2 expansion

    d is the dimensionality of the input data, i and j are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return d * i - (i**2 + 3 * i) // 2 - 1 + j
    else:
        return d * i - (i**2 + i) // 2 + j


cdef inline cnp.int64_t _deg3_column(
    cnp.int64_t d,
    cnp.int64_t i,
    cnp.int64_t j,
    cnp.int64_t k,
    FLAG_t interaction_only
    ) nogil:
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
def py_deg3_column(d, i, j, k, interaction_only):
    """Compute the index of the column for a degree 3 expansion

    d is the dimensionality of the input data, i, j and k are the indices
    for the columns involved in the expansion.
    """
    if interaction_only:
        return ((3 * d**2 * i - 3 * d * i**2 + i**3
                 + 11 * i - 3 * j**2 - 9 * j) // 6
                + i**2 - 2 * d * i + d * j - d + k)
    else:
        return ((3 * d**2 * i - 3 * d * i**2 + i ** 3 - i
                 - 3 * j**2 - 3 * j) // 6
                + d * j + k)

def _csr_polynomial_expansion(
    cnp.ndarray[DATA_t, ndim=1] data,           # TODO: Make const in Cython 3
    cnp.ndarray[INDEX_A_t, ndim=1] indices,     # TODO: Make const in Cython 3
    cnp.ndarray[INDEX_A_t, ndim=1] indptr,      # TODO: Make const in Cython 3
    INDEX_A_t d,
    cnp.ndarray[DATA_t, ndim=1] result_data,
    cnp.ndarray[INDEX_B_t, ndim=1] result_indices,
    cnp.ndarray[INDEX_B_t, ndim=1] result_indptr,
    FLAG_t interaction_only,
    FLAG_t degree
):
    """
    Perform a second-degree polynomial or interaction expansion on a scipy
    compressed sparse row (CSR) matrix. The method used only takes products of
    non-zero features. For a matrix with density :math:`d`, this results in a
    speedup on the order of :math:`(1/d)^k` where :math:`k` is the degree of
    the expansion, assuming all rows are of similar density.

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
    cdef INDEX_A_t row_i, row_starts, row_ends, i, j, k, i_ptr, j_ptr, k_ptr

    cdef INDEX_B_t expanded_index=0, num_cols_in_row, col

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
                        if max(i, j) > MAX_SAFE_INDEX_DEG2:
                            with gil:
                                col = <INDEX_B_t> py_deg2_column(d, i, j, interaction_only)
                        else:
                            col = <INDEX_B_t> _deg2_column(d, i, j, interaction_only)
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
                            if max(i, j, k) > MAX_SAFE_INDEX_DEG3:
                                with gil:
                                    col = <INDEX_B_t> py_deg3_column(d, i, j, k, interaction_only)
                            else:
                                col = <INDEX_B_t> _deg3_column(d, i, j, k, interaction_only)
                            result_indices[expanded_index] = col
                            result_data[expanded_index] = (
                                data[i_ptr] * data[j_ptr] * data[k_ptr])
                            expanded_index += 1
                            num_cols_in_row += 1

            result_indptr[row_i+1] = result_indptr[row_i] + num_cols_in_row
