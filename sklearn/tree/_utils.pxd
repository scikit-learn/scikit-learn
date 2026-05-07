# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.

cimport numpy as cnp
from sklearn.tree._common cimport node_struct, Y_DTYPE_C
from sklearn.neighbors._quad_tree cimport Cell
from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t


cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    # It corresponds to the maximum representable value for
    # 32-bit signed integers (i.e. 2^31 - 1).
    RAND_R_MAX = 2147483647


# safe_realloc(&p, n) resizes the allocation of p to n * sizeof(*p) bytes or
# raises a MemoryError. It never calls free, since that's __dealloc__'s job.
#   cdef float32_t *p = NULL
#   safe_realloc(&p, n)
# is equivalent to p = malloc(n * sizeof(*p)) with error checking.
ctypedef fused realloc_ptr:
    # Add pointer types here as needed.
    (float32_t*)
    (intp_t*)
    (uint8_t*)
    (float64_t*)
    (float64_t**)
    (node_struct*)
    (Cell*)
    (node_struct**)

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil


cdef cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size)


cdef intp_t rand_int(intp_t low, intp_t high,
                     uint32_t* random_state) noexcept nogil


cdef float64_t rand_uniform(float64_t low, float64_t high,
                            uint32_t* random_state) noexcept nogil


cdef float64_t log(float64_t x) noexcept nogil


cdef class WeightedFenwickTree:
    cdef intp_t size         # number of leaves (ranks)
    cdef Y_DTYPE_C* tree_w   # BIT for weights
    cdef Y_DTYPE_C* tree_wy  # BIT for weighted targets
    cdef intp_t max_pow2     # highest power of two <= n
    cdef Y_DTYPE_C total_w   # running total weight
    cdef Y_DTYPE_C total_wy  # running total weighted target

    cdef void reset(self, intp_t size) noexcept nogil
    cdef void add(self, intp_t idx, Y_DTYPE_C y, Y_DTYPE_C w) noexcept nogil
    cdef intp_t search(
        self,
        Y_DTYPE_C t,
        Y_DTYPE_C* cw_out,
        Y_DTYPE_C* cwy_out,
        intp_t* prev_idx_out,
    ) noexcept nogil
