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


# =============================================================================
# WeightedHeap data structure
# =============================================================================

cdef class WeightedHeap:
    """Binary heap with per-item weights, supporting min-heap and max-heap modes.

    Values are stored sign-adjusted internally so that the ordering logic
    is always "min-heap" on the stored buffer:
      - if min_heap: store v
      - else (max-heap): store -v

    Attributes (all should be treated as readonly attributes)
    ----------
    capacity : intp_t
        Allocated capacity for the heap arrays.

    size : intp_t
        Current number of elements in the heap.

    heap : float64_t*
        Array of (possibly sign-adjusted) values that determines ordering.

    weights : float64_t*
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
        self.size = 0
        self.min_heap = min_heap
        self.total_weight = 0.0
        self.weighted_sum = 0.0
        self.heap = NULL
        self.weights = NULL
        # safe_realloc can raise MemoryError -> __cinit__ may propagate
        safe_realloc(&self.heap, capacity)
        safe_realloc(&self.weights, capacity)

    def __dealloc__(self):
        if self.heap != NULL:
            free(self.heap)
        if self.weights != NULL:
            free(self.weights)

    cdef void reset(self) noexcept nogil:
        """Reset to construction state (keeps capacity)."""
        self.size = 0
        self.total_weight = 0.0
        self.weighted_sum = 0.0

    cdef bint is_empty(self) noexcept nogil:
        return self.size == 0

    cdef void push(self, float64_t value, float64_t weight) noexcept nogil:
        """Insert a (value, weight)."""
        cdef intp_t n = self.size
        cdef float64_t stored = value if self.min_heap else -value

        assert n < self.capacity
        # ^ should never raise as capacity is set to the max possible size

        self.heap[n] = stored
        self.weights[n] = weight
        self.size = n + 1

        self.total_weight += weight
        self.weighted_sum += value * weight

        self._heapify_up(n)

    cdef void pop(self, float64_t* value, float64_t* weight) noexcept nogil:
        """Pop top element into pointers."""
        cdef intp_t n = self.size
        assert n > 0

        cdef float64_t stored = self.heap[0]
        cdef float64_t v = stored if self.min_heap else -stored
        cdef float64_t w = self.weights[0]
        value[0] = v
        weight[0] = w

        # Update aggregates
        self.total_weight -= w
        self.weighted_sum -= v * w

        # Move last to root and sift down
        n -= 1
        self.size = n
        if n > 0:
            self.heap[0] = self.heap[n]
            self.weights[0] = self.weights[n]
            self._heapify_down(0)

    cdef float64_t top_weight(self) noexcept nogil:
        assert self.size > 0
        return self.weights[0]

    cdef float64_t top(self) noexcept nogil:
        assert self.size > 0
        cdef float64_t s = self.heap[0]
        return s if self.min_heap else -s

    # Internal helpers (nogil):

    cdef inline void _swap(self, intp_t i, intp_t j) noexcept nogil:
        cdef float64_t tmp = self.heap[i]
        self.heap[i] = self.heap[j]
        self.heap[j] = tmp
        tmp = self.weights[i]
        self.weights[i] = self.weights[j]
        self.weights[j] = tmp

    cdef inline void _heapify_up(self, intp_t i) noexcept nogil:
        """Move up the element at index i until heap invariant is restored."""
        cdef intp_t p
        while i > 0:
            p = (i - 1) >> 1
            if self.heap[i] < self.heap[p]:
                self._swap(i, p)
                i = p
            else:
                break

    cdef inline void _heapify_down(self, intp_t i) noexcept nogil:
        """Move down the element at index i until heap invariant is restored."""
        cdef intp_t n = self.size
        cdef intp_t left, right, mc
        while True:
            left = (i << 1) + 1
            right = left + 1
            if left >= n:
                return
            mc = left
            if right < n and self.heap[right] < self.heap[left]:
                mc = right
            if self.heap[i] > self.heap[mc]:
                self._swap(i, mc)
                i = mc
            else:
                return


cdef class PytestWeightedHeap(WeightedHeap):
    """Used for testing only"""

    def py_push(self, double value, double weight):
        self.push(value, weight)

    def py_pop(self):
        cdef double v, w
        self.pop(&v, &w)
        return v, w


# =============================================================================
# MAE split precomputations algorithm
# =============================================================================

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
    """
    Fill `abs_errors` and `medians`.

    If start < end:
        Computes the "prefix" AEs/medians, i.e the AEs for each set of indices
        sample_indices[start:start + i] with i in {1, ..., n}
        where n = end - start
    Else:
        Computes the "suffix" AEs/medians, i.e the AEs for each set of indices
        sample_indices[i:] with i in {0, ..., n-1}

    Parameters
    ----------
    ys : const float64_t[:, ::1]
        Target values. Shape: (n_samples, n_outputs).
    sample_weight : const float64_t[:]
        Shape: (n_samples,)
    sample_indices : const intp_t[:]
        indices indicating which samples to use. Shape: (n_samples,)
    above : WeightedHeap
    below : WeightedHeap
    k : intp_t
        Dimension to consider in y. In [0, n_outputs - 1].
    start : intp_t
        Start index in `sample_indices`
    end : intp_t
        End index (exclusive) in `sample_indices`
    abs_errors : float64_t[::1]
        array to store (increment) the computed absolute errors. Shape: (n,)
        with n := end - start
    medians : float64_t[::1]
        array to store (overwrite) the computed medians. Shape: (n,)

    Complexity: O(n log n)
    This algorithm is an adaptation of the two heaps solution of
    the "find median from a data stream" problem
    See for instance: https://www.geeksforgeeks.org/dsa/median-of-stream-of-integers-running-integers/

    But here, it's the weighted median and we also need to compute the AE, so:
    - instead of balancing the heaps based on their number of elements,
      rebalance them based on the summed weights of their elements
    - rewrite the AE computation by splitting the sum between elements
      above and below the median, which allow to express it as a simple
      O(1) computation.
      See the maths in the PR desc:
      https://github.com/scikit-learn/scikit-learn/pull/32100
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
    cdef float64_t top_val, top_weight
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
        elif y > below.top():
            above.push(y, w)
        else:
            below.push(y, w)

        half_weight = (above.total_weight + below.total_weight) / 2.0

        # Rebalance heaps
        while above.total_weight < half_weight and not below.is_empty():
            below.pop(&top_val, &top_weight)
            above.push(top_val, top_weight)
        while (
            not above.is_empty()
            and (above.total_weight - above.top_weight()) >= half_weight
        ):
            above.pop(&top_val, &top_weight)
            below.push(top_val, top_weight)

        # Current median
        if above.total_weight > half_weight + 1e-5 * fabs(half_weight):
            median = above.top()
        else:  # above and below weight are almost exactly equals
            median = (above.top() + below.top()) / 2.
        medians[j] = median
        abs_errors[j] += (
            (below.total_weight - above.total_weight) * median
            - below.weighted_sum
            + above.weighted_sum
        )
        p += step
        j += step


def _py_precompute_absolute_errors(
    const float64_t[:, ::1] ys,
    const float64_t[:] sample_weight,
    const intp_t[:] sample_indices,
    bint suffix=False
):
    """
    Used for testing precompute_absolute_errors.
    - If `suffix` is False:
        Computes the "prefix" AEs, i.e the AEs for each set of indices
        sample_indices[:i] with i in {1, ..., n}
    - If `suffix` is True:
        Computes the "suffix" AEs, i.e the AEs for each set of indices
        sample_indices[i:] with i in {0, ..., n-1}
    """
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
