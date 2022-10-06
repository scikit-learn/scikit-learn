# Based on `scipy.sparse.sparsetools.csr.h` under the BSD license.

cimport cython
from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector

cimport numpy as cnp

# Fused type for data
ctypedef fused DATA_t:
    cython.integral
    cython.floating

# Fused types for indices and indptr.
# Need uniquely named fused types to generate full cross-product for all four
# input/output array combinations

# Input CSR indices and indptr type.
# In the Python layer, SciPy's implementation sets indices and indptr
# type to smallest one necessary to maintain indices values
# representability, determined by the maximum value across
# all input blocks after concatenation.
ctypedef fused IND_t:
    cnp.int32_t
    cnp.int64_t

# TODO: `const`-qualify fused-typed memoryviews when Cython>=3.0 is released.
cpdef _csr_hstack(
    const Py_ssize_t n_blocks,      # Number of matrices to stack
    const Py_ssize_t n_rows,        # Number of rows (same across all matrices)
    const cnp.int64_t[:] n_cols,    # Number of columns (one per matrix)
    IND_t[:] indptr_cat,            # Input concatenated array of indptrs
    IND_t[:] indices_cat,           # Input concatenated array of indices
    DATA_t[:] data_cat,             # Input concatenated array of data
    IND_t[:] indptr,                # Output array to write indptr into
    IND_t[:] indices,               # Output array to write indices into
    DATA_t[:] data,                 # Output array to write data into
    ):

    cdef:
        Py_ssize_t offset, row_start, row_end, row_sum
        Py_ssize_t idx, jdx, kdx, cat_index

        vector[Py_ssize_t] col_offset = vector[Py_ssize_t](n_blocks)
        vector[Py_ssize_t] indptr_bound = vector[Py_ssize_t](n_blocks)
        vector[Py_ssize_t] indices_bound = vector[Py_ssize_t](n_blocks)

    # The bounds will store the locations/indices of the flat concatenated
    # arrays that correspond to each block.
    with nogil:

        # We set the initial index values here and update iteratively
        col_offset[0] = 0
        indptr_bound[0] = 0
        indices_bound[0] = 0

        # We populate the bounds iteratively
        for idx in range(1, n_blocks):
            col_offset[idx] = col_offset[idx - 1] + n_cols[idx - 1]
            indptr_bound[idx] = indptr_bound[idx - 1] + n_rows + 1
            indices_bound[idx] = indices_bound[idx - 1] + indptr_cat[indptr_bound[idx - 1] + n_rows]

        # Now we populate the full output matrix
        indptr[0] = 0
        row_sum = 0

        # We iterate across rows for convenience
        for idx in range(n_rows):
            # For each row, in each matrix we find the valid indices to iterate
            for jdx in range(n_blocks):
                row_start = indptr_cat[indptr_bound[jdx] + idx]
                row_end = indptr_cat[indptr_bound[jdx] + idx + 1]

                # We account for the offset due to horizontal concatenation
                offset = col_offset[jdx]

                # We iterate over the valid indices, updating the indices and
                # data
                for kdx in range(row_end - row_start):
                    cat_index = indices_bound[jdx] + row_start + kdx
                    indices[row_sum + kdx] = indices_cat[cat_index] + offset
                    data[row_sum + kdx] = data_cat[cat_index]
                row_sum += row_end - row_start

            # Cumulative row index (indptr)
            indptr[idx + 1] = row_sum
    return
