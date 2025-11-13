# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.stdlib cimport free
from libc.stdlib cimport realloc
from libc.math cimport log as ln
from libc.math cimport isnan
from libc.string cimport memset

import numpy as np
cimport numpy as cnp
cnp.import_array()

from sklearn.utils._random cimport our_rand_r

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


cdef class WeightedFenwickTree:
    """
    Fenwick tree (Binary Indexed Tree) specialized for maintaining:
      - prefix sums of weights
      - prefix sums of weight * target (y)

    Notes:
      - Implementation uses 1-based indexing internally for the Fenwick tree
        arrays, hence the +1 sized buffers. 1-based indexing is customary for this
        data structure and makes the some index handling slightly more efficient and
        natural.
      - Memory ownership: this class allocates and frees the underlying C buffers.
      - Typical operations:
          add(rank, y, w) -> O(log n)
          search(t)       -> O(log n), finds the smallest rank with
                             cumulative weight > t (see search for details).
    """

    def __cinit__(self, intp_t capacity):
        self.tree_w = NULL
        self.tree_wy = NULL

        # Allocate arrays of length (capacity + 1) because indices are 1-based.
        safe_realloc(&self.tree_w, capacity + 1)
        safe_realloc(&self.tree_wy, capacity + 1)

    cdef void reset(self, intp_t size) noexcept nogil:
        """
        Reset the tree to hold 'size' elements and clear all aggregates.
        """
        cdef intp_t p
        cdef intp_t n_bytes = (size + 1) * sizeof(float64_t)  # +1 for 1-based storage

        # Public size and zeroed aggregates.
        self.size = size
        memset(self.tree_w, 0, n_bytes)
        memset(self.tree_wy, 0, n_bytes)
        self.total_w = 0.0
        self.total_wy = 0.0

        # highest power of two <= size
        p = 1
        while p <= size:
            p <<= 1
        self.max_pow2 = p >> 1

    def __dealloc__(self):
        if self.tree_w != NULL:
            free(self.tree_w)
        if self.tree_wy != NULL:
            free(self.tree_wy)

    cdef void add(self, intp_t idx, float64_t y_value, float64_t weight) noexcept nogil:
        """
        Add a weighted observation to the Fenwick tree.

        Parameters
        ----------
        idx : intp_t
            The 0-based index where to add the observation
        y_value : float64_t
            The target value (y) of the observation
        weight : float64_t
            The sample weight

        Notes
        -----
        Updates both weight sums and weighted target sums in O(log n) time.
        """
        cdef float64_t weighted_y = weight * y_value
        cdef intp_t fenwick_idx = idx + 1  # Convert to 1-based indexing

        # Update Fenwick tree nodes by traversing up the tree
        while fenwick_idx <= self.size:
            self.tree_w[fenwick_idx] += weight
            self.tree_wy[fenwick_idx] += weighted_y
            # Move to next node using bit manipulation: add lowest set bit
            fenwick_idx += fenwick_idx & -fenwick_idx

        # Update global totals
        self.total_w += weight
        self.total_wy += weighted_y

    cdef intp_t search(
        self,
        float64_t target_weight,
        float64_t* cumul_weight_out,
        float64_t* cumul_weighted_y_out,
        intp_t* prev_idx_out,
    ) noexcept nogil:
        """
        Binary search to find the position where cumulative weight reaches target.

        This method performs a binary search on the Fenwick tree to find indices
        such that the cumulative weight at 'prev_idx' is < target_weight and
        the cumulative weight at the returned index is >= target_weight.

        Parameters
        ----------
        target_weight : float64_t
            The target cumulative weight to search for
        cumul_weight_out : float64_t*
            Output pointer for cumulative weight up to returned index (exclusive)
        cumul_weighted_y_out : float64_t*
            Output pointer for cumulative weighted y-sum up to returned index (exclusive)
        prev_idx_out : intp_t*
            Output pointer for the previous index (largest index with cumul_weight < target)

        Returns
        -------
        intp_t
            The index where cumulative weight first reaches or exceeds target_weight

        Notes
        -----
        - O(log n) complexity
        - Ignores nodes with zero weights (corresponding to uninserted y-values)
        - Assumes at least one active (positive-weight) item exists
        - Assumes 0 <= target_weight <= total_weight
        """
        cdef:
            intp_t current_idx = 0
            intp_t next_idx, prev_idx, equal_bit
            float64_t cumul_weight = 0.0
            float64_t cumul_weighted_y = 0.0
            intp_t search_bit = self.max_pow2  # Start from highest power of 2
            float64_t node_weight, equal_target

        # Phase 1: Standard Fenwick binary search with prefix accumulation
        # Traverse down the tree, moving right when we can consume more weight
        while search_bit != 0:
            next_idx = current_idx + search_bit
            if next_idx <= self.size:
                node_weight = self.tree_w[next_idx]
                if target_weight == node_weight:
                    # Exact match found - store state for later processing
                    equal_target = target_weight
                    equal_bit = search_bit
                    break
                elif target_weight > node_weight:
                    # We can consume this node's weight - move right and accumulate
                    target_weight -= node_weight
                    current_idx = next_idx
                    cumul_weight += node_weight
                    cumul_weighted_y += self.tree_wy[next_idx]
            search_bit >>= 1

        # If no exact match, we're done with standard search
        if search_bit == 0:
            cumul_weight_out[0] = cumul_weight
            cumul_weighted_y_out[0] = cumul_weighted_y
            prev_idx_out[0] = current_idx
            return current_idx

        # Phase 2: Handle exact match case - find prev_idx
        # Search for the largest index with cumulative weight < original target
        prev_idx = current_idx
        while search_bit != 0:
            next_idx = prev_idx + search_bit
            if next_idx <= self.size:
                node_weight = self.tree_w[next_idx]
                if target_weight > node_weight:
                    target_weight -= node_weight
                    prev_idx = next_idx
            search_bit >>= 1

        # Phase 3: Complete the exact match search
        # Restore state and search for the largest index with
        # cumulative weight <= original target (and this is case, we know we have ==)
        search_bit = equal_bit
        target_weight = equal_target
        while search_bit != 0:
            next_idx = current_idx + search_bit
            if next_idx <= self.size:
                node_weight = self.tree_w[next_idx]
                if target_weight >= node_weight:
                    target_weight -= node_weight
                    current_idx = next_idx
                    cumul_weight += node_weight
                    cumul_weighted_y += self.tree_wy[next_idx]
            search_bit >>= 1

        # Output results
        cumul_weight_out[0] = cumul_weight
        cumul_weighted_y_out[0] = cumul_weighted_y
        prev_idx_out[0] = prev_idx
        return current_idx


cdef class PytestWeightedFenwickTree(WeightedFenwickTree):
    """Used for testing only"""

    def py_reset(self, intp_t n):
        self.reset(n)

    def py_add(self, intp_t idx, float64_t y, float64_t w):
        self.add(idx, y, w)

    def py_search(self, float64_t t):
        cdef float64_t w, wy
        cdef intp_t prev_idx
        idx = self.search(t, &w, &wy, &prev_idx)
        return prev_idx, idx, w, wy
