from ._tree cimport Node

import numpy as np
cimport numpy as np

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef np.npy_uint32 UINT32_t          # Unsigned 32 bit integer

# safe_realloc(&p, n) resizes the allocation of p to n * sizeof(*p) bytes or
# raises a MemoryError. It never calls free, since that's __dealloc__'s job.
#   cdef DTYPE_t *p = NULL
#   safe_realloc(&p, n)
# is equivalent to p = malloc(n * sizeof(*p)) with error checking.
ctypedef fused realloc_ptr:
# Add pointer types here as needed.
    (DTYPE_t*)
    (SIZE_t*)
    (unsigned char*)
    (DOUBLE_t*)
    (DOUBLE_t**)
    (Node*)
    (Node**)


cdef realloc_ptr safe_realloc(realloc_ptr* p, size_t nelems) nogil except *


cdef np.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size)
