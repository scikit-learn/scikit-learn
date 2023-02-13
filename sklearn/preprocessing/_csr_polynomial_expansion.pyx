# Authors: Andrew nystrom <awnystrom@gmail.com>
#          Meekail Zain <zainmeekail@gmail.com>
from libc.limits cimport INT_MAX, LONG_MAX
from libc.math cimport sqrt, pow
cimport numpy as cnp
cnp.import_array()

from scipy.sparse import csr_matrix
import numpy as np

ctypedef cnp.int8_t FLAG_t


# INDEX_{A,B}_t are defined to generate a proper Cartesian product
# of types through Cython fused-type expansion.
ctypedef fused INDEX_A_t:
    cnp.int32_t
    cnp.int64_t
ctypedef fused INDEX_B_t:
    cnp.int32_t
    cnp.int64_t

ctypedef fused DATA_t:
    cnp.int32_t
    cnp.int64_t
    cnp.float32_t
    cnp.float64_t

cdef inline cnp.int64_t _deg2_column(
    cnp.int64_t n_features,
    cnp.int64_t i,
    cnp.int64_t j,
    FLAG_t interaction_only
) nogil:
    """Compute the index of the column for a degree 2 expansion

    n_features is the dimensionality of the input data, i and j are the indices
    for the columns involved in the expansion.
    """

    # This is a conservative upper bound for the maximum value before the
    # intermediate computation i**2 + 3*i overflows.
    # Corresponds to to i**2 + 3*i = maxint32
    # i = sqrt(maxint32) - 4
    cdef cnp.int64_t MAX_SAFE_INDEX_DEG2 = 46340

    if max(i, j) > MAX_SAFE_INDEX_DEG2 or i > LONG_MAX // n_features :
        # In this case, the Cython implementation
        # would result in an integer overflow.
        # Here, we take advantage of `PyLong` for arbitrary precision.
        with gil:
            return <cnp.int64_t> py_deg2_column(n_features, i, j, interaction_only)
    else:
        if interaction_only:
            return n_features * i - (i**2 + 3 * i) / 2 - 1 + j
        else:
            return n_features * i - (i**2 + i) / 2 + j

def py_deg2_column(n_features, i, j, interaction_only):
    if interaction_only:
        return n_features * i - (i**2 + 3 * i) // 2 - 1 + j
    else:
        return n_features * i - (i**2 + i) // 2 + j


cdef inline cnp.int64_t _deg3_column(
    cnp.int64_t n_features,
    cnp.int64_t i,
    cnp.int64_t j,
    cnp.int64_t k,
    FLAG_t interaction_only
    ) nogil:
    """Compute the index of the column for a degree 3 expansion

    n_features is the dimensionality of the input data, i, j and k are the indices
    for the columns involved in the expansion.
    """
    # This is a conservative upper bound for the maximum value before the
    # intermediate computation 3 * d**2 * d + d**3 overflows. We leverage the fact
    # that i<=d.
    # Corresponds to 3 * n_features**2 * d + d**3 = maxint32
    # when n_features == d.
    cdef cnp.int64_t MAX_SAFE_INDEX_DEG3 = <cnp.int64_t> pow(LONG_MAX, 1/3)/4

    if max(i, j, k) > MAX_SAFE_INDEX_DEG3 or i // 3 + 1 > LONG_MAX // n_features // n_features:
        # In this case, the Cython implementation
        # would result in an integer overflow.
        # Here, we take advantage of `PyLong` for arbitrary precision.
        with gil:
            return <cnp.int64_t> py_deg3_column(n_features, i, j, k, interaction_only)
    if interaction_only:
        return (
            (3 * n_features**2 * i - 3 * n_features * i**2 + i**3
            + 11 * i - 3 * j**2 - 9 * j) / 6
            + i**2 - 2 * n_features * i + n_features * j - n_features + k
        )
    else:
        return (
            (3 * n_features**2 * i - 3 * n_features * i**2 + i ** 3 - i
            - 3 * j**2 - 3 * j) / 6 + n_features * j + k
        )

def py_deg3_column(n_features, i, j, k, interaction_only):
    if interaction_only:
        return (
            (3 * n_features**2 * i - 3 * n_features * i**2 + i**3
            + 11 * i - 3 * j**2 - 9 * j) // 6
            + i**2 - 2 * n_features * i + n_features * j - n_features + k
        )
    else:
        return (
            (3 * n_features**2 * i - 3 * n_features * i**2 + i ** 3 - i
            - 3 * j**2 - 3 * j) // 6 + n_features * j + k
        )


def py_calc_expanded_nnz_deg2(n, interaction_only):
    return (n**2 + n) // 2 - interaction_only * n

def py_calc_expanded_nnz_deg3(n, interaction_only):
    return (n**3 + 3 * n**2 + 2 * n) // 6 - interaction_only * n**2

cpdef cnp.int64_t _calc_expanded_nnz(
    cnp.int64_t n,
    FLAG_t interaction_only,
    cnp.int64_t degree
):
    """
    Calculates the number of non-zero interaction terms generated by the
    non-zero elements of a single row.
    """
    # This is the maximum value before the intermediate computation
    # d**2 + d overflows
    # Solution to d**2 + d = maxint32
    cdef cnp.int64_t MAX_SAFE_INDEX_CALC_DEG2 = <cnp.int64_t> sqrt(LONG_MAX)

    # This is the maximum value before the intermediate computation
    # d**3 + 3 * d**2 + 2*d overflows
    # Solution to d**3 + 3 * d**2 + 2*d = maxint32
    cdef cnp.int64_t MAX_SAFE_INDEX_CALC_DEG3 = 2097151

    if degree == 2:
        if n <= MAX_SAFE_INDEX_CALC_DEG2:
            return (n**2 + n) / 2 - interaction_only * n
        return <cnp.int64_t> py_calc_expanded_nnz_deg2(n, interaction_only)
    else:
        if n <= MAX_SAFE_INDEX_CALC_DEG3:
            return (n**3 + 3 * n**2 + 2 * n) / 6 - interaction_only * n**2
        return <cnp.int64_t> py_calc_expanded_nnz_deg3(n, interaction_only)

cpdef cnp.int64_t _calc_total_nnz(
    INDEX_A_t[:] indptr,
    FLAG_t interaction_only,
    cnp.int64_t degree,
):
    """
    Calculates the number of non-zero interaction terms generated by the
    non-zero elements across all rows for a single degree.
    """
    cdef cnp.int64_t total_nnz=0
    cdef cnp.intp_t row_idx
    for row_idx in range(len(indptr) - 1):
        total_nnz += _calc_expanded_nnz(
            indptr[row_idx + 1] - indptr[row_idx],
            interaction_only,
            degree
        )
    return total_nnz


cpdef void _csr_polynomial_expansion(
    const DATA_t[:] data,                 # IN READ-ONLY
    const INDEX_A_t[:] indices,           # IN READ-ONLY
    const INDEX_A_t[:] indptr,            # IN READ-ONLY
    INDEX_A_t n_features,
    DATA_t[:] result_data,          # OUT
    INDEX_B_t[:] result_indices,    # OUT
    INDEX_B_t[:] result_indptr,     # OUT
    FLAG_t interaction_only,
    FLAG_t degree
) nogil:
    """
    Perform a second or third degree polynomial or interaction expansion on a
    compressed sparse row (CSR) matrix. The method used only takes products of
    non-zero features. For a matrix with density :math:`d`, this results in a
    speedup on the order of :math:`(1/d)^k` where :math:`k` is the degree of
    the expansion, assuming all rows are of similar density.

    Parameters
    ----------
    data : memory view on nd-array
        The "data" attribute of the input CSR matrix.

    indices : memory view on nd-array
        The "indices" attribute of the input CSR matrix.

    indptr : memory view on nd-array
        The "indptr" attribute of the input CSR matrix.

    n_features : int
        The dimensionality of the input CSR matrix.

    result_data : nd-array
        The output CSR matrix's "data" attribute.
        It is modified by this routine.

    result_indices : nd-array
        The output CSR matrix's "indices" attribute.
        It is modified by this routine.

    result_indptr : nd-array
        The output CSR matrix's "indptr" attribute.
        It is modified by this routine.

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
                        col = <INDEX_B_t> _deg2_column(
                            n_features,
                            i, j,
                            interaction_only
                        )
                        result_indices[expanded_index] = col
                        result_data[expanded_index] = (
                            data[i_ptr] * data[j_ptr]
                        )
                        expanded_index += 1
                        num_cols_in_row += 1
                    else:
                        # degree == 3
                        for k_ptr in range(j_ptr + interaction_only, row_ends):
                            k = indices[k_ptr]
                            col = <INDEX_B_t> _deg3_column(
                                n_features,
                                i, j, k,
                                interaction_only
                            )
                            result_indices[expanded_index] = col
                            result_data[expanded_index] = (
                                data[i_ptr] * data[j_ptr] * data[k_ptr]
                            )
                            expanded_index += 1
                            num_cols_in_row += 1

            result_indptr[row_i+1] = result_indptr[row_i] + num_cols_in_row
    return
