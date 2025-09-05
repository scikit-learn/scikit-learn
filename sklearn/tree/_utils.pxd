# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.

cimport numpy as cnp
from ._tree cimport Node
from ..neighbors._quad_tree cimport Cell
from ..utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t


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
    (Node*)
    (Cell*)
    (Node**)

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil


cdef cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size)


cdef intp_t rand_int(intp_t low, intp_t high,
                     uint32_t* random_state) noexcept nogil


cdef float64_t rand_uniform(float64_t low, float64_t high,
                            uint32_t* random_state) noexcept nogil

cdef class WeightedHeap:
    cdef intp_t capacity
    cdef intp_t size
    cdef float64_t* heap
    cdef float64_t* weights
    cdef float64_t total_weight
    cdef float64_t weighted_sum
    cdef bint min_heap

    cdef void reset(self) noexcept nogil
    cdef bint is_empty(self) noexcept nogil
    cdef void push(self, float64_t value, float64_t weight) noexcept nogil
    cdef void pop(self, float64_t* value, float64_t* weight) noexcept nogil
    cdef float64_t top_weight(self) noexcept nogil
    cdef float64_t top(self) noexcept nogil
    cdef void _swap(self, intp_t, intp_t) noexcept nogil
    cdef void _perc_up(self, intp_t) noexcept nogil
    cdef void _perc_down(self, intp_t) noexcept nogil

cdef void precompute_absolute_errors(
    const float64_t[:, ::1] ys,
    const float64_t[:] sample_weight,
    const intp_t[:] sample_indices,
    WeightedHeap above,
    WeightedHeap below,
    intp_t k,
    intp_t start,
    intp_t end,
    float64_t[::1] abs_errors,
    float64_t[::1] medians
) noexcept nogil

cdef float64_t log(float64_t x) noexcept nogil
