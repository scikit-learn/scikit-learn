# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.stdlib cimport free
from libc.stdlib cimport realloc
from libc.math cimport log as ln
from libc.math cimport isnan
from libc.math cimport fabs

import numpy as np
cimport numpy as cnp
cnp.import_array()

from ..utils._random cimport our_rand_r

# =============================================================================
# Helper functions
# =============================================================================

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil:
    # sizeof(realloc_ptr[0]) would be more like idiomatic C, but causes Cython
    # 0.20.1 to crash.
    cdef size_t nbytes = nelems * sizeof(p[0][0])
    if nbytes / sizeof(p[0][0]) != nelems:
        # Overflow in the multiplication
        raise MemoryError(f"could not allocate ({nelems} * {sizeof(p[0][0])}) bytes")

    cdef realloc_ptr tmp = <realloc_ptr>realloc(p[0], nbytes)
    if tmp == NULL:
        raise MemoryError(f"could not allocate {nbytes} bytes")

    p[0] = tmp
    return 0


def _realloc_test():
    # Helper for tests. Tries to allocate <size_t>(-1) / 2 * sizeof(size_t)
    # bytes, which will always overflow.
    cdef intp_t* p = NULL
    safe_realloc(&p, <size_t>(-1) / 2)
    if p != NULL:
        free(p)
        assert False


cdef inline cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size):
    """Return copied data as 1D numpy array of intp's."""
    cdef cnp.npy_intp shape[1]
    shape[0] = <cnp.npy_intp> size
    return cnp.PyArray_SimpleNewFromData(1, shape, cnp.NPY_INTP, data).copy()


cdef inline intp_t rand_int(intp_t low, intp_t high,
                            uint32_t* random_state) noexcept nogil:
    """Generate a random integer in [low; end)."""
    return low + our_rand_r(random_state) % (high - low)


cdef inline float64_t rand_uniform(float64_t low, float64_t high,
                                   uint32_t* random_state) noexcept nogil:
    """Generate a random float64_t in [low; high)."""
    return ((high - low) * <float64_t> our_rand_r(random_state) /
            <float64_t> RAND_R_MAX) + low


cdef inline float64_t log(float64_t x) noexcept nogil:
    return ln(x) / ln(2.0)

# =============================================================================
# WeightedHeap data structure
# =============================================================================

cdef class WeightedHeap:
    """Binary heap with per-item weights, supporting min-heap and max-heap modes.

    Values are stored sign-adjusted internally so that the ordering logic
    is always "min-heap" on the stored buffer:
      - if min_heap: store v
      - else (max-heap): store -v

    Attributes
    ----------
    capacity : intp_t
        Allocated capacity for the heap arrays.

    size_ : intp_t
        Current number of elements in the heap.

    heap_ : float64_t*
        Array of (possibly sign-adjusted) values that determines ordering.

    weights_ : float64_t*
        Parallel array of weights.

    total_weight : float64_t
        Sum of all weights currently in the heap.

    weighted_sum : float64_t
        Sum over items of (original_value * weight), i.e. without sign-adjustment.

    min_heap : bint
        If True, behaves as a min-heap; if False, behaves as a max-heap.
    """

    def __cinit__(self, intp_t capacity, bint min_heap=True):
        if capacity <= 0:
            capacity = 1
        self.capacity = capacity
        self.size_ = 0
        self.min_heap = min_heap
        self.total_weight = 0.0
        self.weighted_sum = 0.0
        self.heap_ = NULL
        self.weights_ = NULL
        # safe_realloc can raise MemoryError -> __cinit__ may propagate
        safe_realloc(&self.heap_, capacity)
        safe_realloc(&self.weights_, capacity)

    def __dealloc__(self):
        if self.heap_ != NULL:
            free(self.heap_)
        if self.weights_ != NULL:
            free(self.weights_)

    cdef void reset(self) noexcept nogil:
        """Reset to construction state (keeps capacity)."""
        self.size_ = 0
        self.total_weight = 0.0
        self.weighted_sum = 0.0

    cdef bint is_empty(self) noexcept nogil:
        return self.size_ == 0

    cdef intp_t size(self) noexcept nogil:
        return self.size_

    cdef int push(self, float64_t value, float64_t weight) except -1 nogil:
        """Insert a (value, weight). Returns 0 or raises MemoryError on alloc fail."""
        cdef intp_t n = self.size_
        cdef float64_t stored = value if self.min_heap else -value

        if n >= self.capacity:
            # should never happen as capacity is set to the max possible size
            return -1

        self.heap_[n] = stored
        self.weights_[n] = weight
        self.size_ = n + 1

        self.total_weight += weight
        self.weighted_sum += value * weight

        self._perc_up(n)
        return 0

    cdef int pop(self, float64_t* value, float64_t* weight) noexcept nogil:
        """Pop top element into pointers. Returns 0 on success, -1 if empty."""
        cdef intp_t n = self.size_
        if n == 0:
            return -1

        self._peek_raw(value, weight)

        # Update aggregates with *original* value (undo sign for max-heap)
        cdef float64_t orig_v = value[0]
        cdef float64_t w = weight[0]
        self.total_weight -= w
        self.weighted_sum -= orig_v * w

        # Move last to root and sift down
        n -= 1
        self.size_ = n
        if n > 0:
            self.heap_[0] = self.heap_[n]
            self.weights_[0] = self.weights_[n]
            self._perc_down(0)
        return 0

    cdef float64_t get_total_weight(self) noexcept nogil:
        return self.total_weight

    cdef float64_t get_weighted_sum(self) noexcept nogil:
        return self.weighted_sum

    cdef float64_t top_weight(self) noexcept nogil:
        if self.size_ == 0:
            return 0.0
        return self.weights_[0]

    cdef float64_t top(self) noexcept nogil:
        if self.size_ == 0:
            return 0.0
        cdef float64_t s = self.heap_[0]
        return s if self.min_heap else -s

    # ----------------------------
    # Internal helpers (nogil)
    # ----------------------------

    cdef inline void _peek_raw(self, float64_t* value, float64_t* weight) noexcept nogil:
        """Internal: read top with proper sign restoration."""
        cdef float64_t stored = self.heap_[0]
        value[0] = stored if self.min_heap else -stored
        weight[0] = self.weights_[0]

    cdef inline void _swap(self, intp_t i, intp_t j) noexcept nogil:
        cdef float64_t tv = self.heap_[i]
        cdef float64_t tw = self.weights_[i]
        self.heap_[i] = self.heap_[j]
        self.weights_[i] = self.weights_[j]
        self.heap_[j] = tv
        self.weights_[j] = tw

    cdef inline void _perc_up(self, intp_t i) noexcept nogil:
        cdef intp_t p
        while i > 0:
            p = (i - 1) >> 1
            if self.heap_[i] < self.heap_[p]:
                self._swap(i, p)
                i = p
            else:
                break

    cdef inline void _perc_down(self, intp_t i) noexcept nogil:
        cdef intp_t n = self.size_
        cdef intp_t left, right, mc
        while True:
            left = (i << 1) + 1
            right = left + 1
            if left >= n:
                return
            mc = left
            if right < n and self.heap_[right] < self.heap_[left]:
                mc = right
            if self.heap_[i] > self.heap_[mc]:
                self._swap(i, mc)
                i = mc
            else:
                return


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
) noexcept nogil:
    """Fill `abs_errors` with the optimal AEs for (y[:i], w[:i])
    i in [1, n].

    Parameters
    ----------
    y : 1D float64_t[::1]
        Values.
    w : 1D float64_t[::1]
        Sample weights.
    abs_errors : 1D float64_t[::1]
        Output buffer, must have shape (n,).
    """
    cdef intp_t j, p, i, step, n
    if start < end:
        j = 0
        step = 1
        n = end - start
    else:
        n = start - end
        step = -1
        j = n - 1

    above.reset()
    below.reset()
    cdef float64_t y
    cdef float64_t w = 1.0
    cdef float64_t val = 0.0
    cdef float64_t wt = 0.0
    cdef float64_t median = 0.0
    cdef float64_t half_weight

    p = start
    for _ in range(n):
        i = sample_indices[p]
        if sample_weight is not None:
            w = sample_weight[i]
        y = ys[i, k]

        # Insert into the appropriate heap
        if below.is_empty():
            above.push(y, w)
        else:
            if y > below.top():
                above.push(y, w)
            else:
                below.push(y, w)

        half_weight = (above.get_total_weight() + below.get_total_weight()) / 2.0

        # Rebalance heaps
        while above.get_total_weight() < half_weight and not below.is_empty():
            if below.pop(&val, &wt) == 0:
                above.push(val, wt)
        while (
            not above.is_empty()
            and (above.get_total_weight() - above.top_weight()) >= half_weight
        ):
            if above.pop(&val, &wt) == 0:
                below.push(val, wt)

        # Current median
        if above.get_total_weight() > half_weight + 1e-5 * fabs(half_weight):
            median = above.top()
        else:  # above and below weight are almost exactly equals
            median = (above.top() + below.top()) / 2.
        medians[j] = median
        abs_errors[j] = (
            (below.get_total_weight() - above.get_total_weight()) * median
            - below.get_weighted_sum()
            + above.get_weighted_sum()
        )
        p += step
        j += step


def _py_precompute_absolute_errors(
    const float64_t[:, ::1] ys,
    const float64_t[:] sample_weight,
    const intp_t[:] sample_indices,
    bint suffix=False
):
    """Used for testing precompute_absolute_errors"""
    cdef:
        intp_t n = sample_weight.size
        WeightedHeap above = WeightedHeap(n, True)
        WeightedHeap below = WeightedHeap(n, False)
        intp_t k = 0
        intp_t start = 0
        intp_t end = n
        float64_t[::1] abs_errors = np.zeros(n, dtype=np.float64)
        float64_t[::1] medians = np.zeros(n, dtype=np.float64)

    if suffix:
        start = n - 1
        end = -1

    precompute_absolute_errors(
        ys, sample_weight, sample_indices, above, below,
        k, start, end, abs_errors, medians
    )
    return np.asarray(abs_errors)


def _any_isnan_axis0(const float32_t[:, :] X):
    """Same as np.any(np.isnan(X), axis=0)"""
    cdef:
        intp_t i, j
        intp_t n_samples = X.shape[0]
        intp_t n_features = X.shape[1]
        uint8_t[::1] isnan_out = np.zeros(X.shape[1], dtype=np.bool_)

    with nogil:
        for i in range(n_samples):
            for j in range(n_features):
                if isnan_out[j]:
                    continue
                if isnan(X[i, j]):
                    isnan_out[j] = True
                    break
    return np.asarray(isnan_out)
