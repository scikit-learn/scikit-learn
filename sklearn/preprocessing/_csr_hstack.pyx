# distutils: language=c++
# Based on `scipy.sparse.sparsetools.csr.h` under the BSD license.

cimport cython
from libc.stdlib cimport malloc, free

cimport numpy as cnp
cnp.import_array()

# Fused type for data
ctypedef fused DATA_t:
    cython.numeric

# Fused types for indices and indptr.
# Need uniquely named fused types to generate full cross-product for all four
# input/output array combinations

# Input IND type. In the Python layer this will be set to smallest type
# necessary to maintain integrity, determined by the maximum value across
# all input blocks.
ctypedef fused IND_IN_t:
    cnp.int32_t
    cnp.int64_t

# Output IND type. In the Python layer this will be set to smallest type
# necessary to maintain integrity, determined by the maximum post-concatenation
# index/indptr value.
ctypedef fused IND_OUT_t:
    cnp.int32_t
    cnp.int64_t

# TODO Add const keyword to input arrays in Cython 3.X
cpdef void _csr_hstack(
    const Py_ssize_t n_blocks,      # Number of matrices to stack
    const Py_ssize_t n_rows,        # Number of rows (same across all matrices)
    const cnp.int32_t[:] n_cols,    # Number of columns (one per matrix)
    IND_IN_t[:] indptr_cat,         # Input concatenated array of indptrs
    IND_IN_t[:] indices_cat,        # Input concatenated array of indices
    DATA_t[:] data_cat,             # Input concatenated array of data
    IND_OUT_t[:] indptr,            # Output array to write indptr into
    IND_OUT_t[:] indices,           # Output array to write indices into
    DATA_t[:] data,                 # Output array to write data into
    ) nogil:

    cdef:
        Py_ssize_t offset, row_start, row_end, row_sum, idx, jdx, kdx
        IND_IN_t* col_offset
        IND_IN_t** indptr_bound
        IND_IN_t** indices_bound
        DATA_t** data_bound

    # The bounds will store the locations/indices of the flat concatenated
    # arrays that correspond to each block. Need to dynamically allocate their
    # memory.
    with nogil:
        col_offset = <IND_IN_t*> malloc(n_blocks * sizeof(IND_IN_t))
        indptr_bound = <IND_IN_t**> malloc(n_blocks * sizeof(IND_IN_t*))
        indices_bound = <IND_IN_t**> malloc(n_blocks * sizeof(IND_IN_t*))
        data_bound = <DATA_t**> malloc(n_blocks * sizeof(DATA_t*))

        # We set the initial pointers/values here and update iteratively
        col_offset[0] = 0
        indptr_bound[0] = &indptr_cat[0]
        indices_bound[0] = &indices_cat[0]
        data_bound[0] = &data_cat[0]

        # We populate the bounds iteratively
        for idx in range(1, n_blocks):
            col_offset[idx] = col_offset[idx - 1] + n_cols[idx - 1]
            indptr_bound[idx] = indptr_bound[idx - 1] + n_rows + 1
            indices_bound[idx] = indices_bound[idx - 1] + indptr_bound[idx - 1][n_rows]
            data_bound[idx] = data_bound[idx - 1] + indptr_bound[idx - 1][n_rows]

        # Now we populate the full output matrix
        indptr[0] = 0
        row_sum = 0

        # We iterate across rows for convenience
        for idx in range(n_rows):
            # For each row, in each matrix we find the valid indices to iterate
            for jdx in range(n_blocks):
                row_start = indptr_bound[jdx][idx]
                row_end = indptr_bound[jdx][idx+1]

                # We account for the offset due to horizontal concatenation
                offset = col_offset[jdx]

                # We iterate over the valid indices, updating the indices and
                # data
                for kdx in range(row_start, row_end):
                    indices[row_sum + kdx] = indices_bound[jdx][kdx] + offset
                    data[row_sum + kdx] = data_bound[jdx][kdx]
                row_sum += row_end - row_start

            # Cumulative row index (indptr)
            indptr[idx + 1] = row_sum

        # Free the memory allocated earlier
        free(col_offset)
        free(indptr_bound)
        free(indices_bound)
        free(data_bound)
        return
