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


# =============================================================================
# WeightedFenwickTree data structure
# =============================================================================

cdef class WeightedFenwickTree:
    """
    Fenwick tree (Binary Indexed Tree) for maintaining:
      - prefix sums of weights
      - prefix sums of weight*value (targets)

    Supports:
      - add(rank, w, y): point update at 'rank'
      - search(t): find the smallest rank with cumulative weight > t (or >= t),
                   also returns prefix aggregates excluding that rank.
    """

    def __cinit__(self, intp_t capacity):
        self.tree_w = NULL
        self.tree_wy = NULL
        # safe_realloc can raise MemoryError -> __cinit__ may propagate
        safe_realloc(&self.tree_w, capacity + 1)
        safe_realloc(&self.tree_wy, capacity + 1)

    cdef void reset(self, intp_t size) noexcept nogil:
        cdef intp_t p
        cdef intp_t n_bytes = (size + 1) * sizeof(float64_t)
        # +1 because 1-based

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

    cdef void add(self, intp_t idx, float64_t y, float64_t w) noexcept nogil:
        cdef float64_t wy = w * y
        idx += 1  # 1-based

        while idx <= self.size:
            self.tree_w[idx] += w
            self.tree_wy[idx] += wy
            idx += idx & -idx

        self.total_w += w
        self.total_wy += wy

    cdef intp_t search(
        self,
        float64_t t,
        float64_t* cw_out,
        float64_t* cwy_out,
        bint inclusive
    ) noexcept nogil:
        """
        Find the smallest index such that
          prefix_weight <= t < prefix_weight if inclusive

        and return:
          - idx:
          - cw (write in cw_out): prefix weight up to idx exclusive
          - cwy (write in cwy_out): prefix weighted sum to idx exclusive

        Notes:
          * Assumes there is at least one active (positive-weight) item.
          * If t >= total weight (can happen with alpha ~ 1), we clamp t slightly.
        """
        cdef:
            intp_t idx = 0
            float64_t cw = 0.0
            float64_t cwy = 0.0
            intp_t bit = self.max_pow2
            float64_t w

        # Standard Fenwick lower-bound with simultaneous prefix accumulation
        while bit != 0:
            next_idx = idx + bit
            if next_idx <= self.size:
                w = self.tree_w[next_idx]
                if (t > w) or (inclusive and t >= w):
                    t -= w
                    idx = next_idx
                    cw += w
                    cwy += self.tree_wy[next_idx]
            bit >>= 1

        cw_out[0] = cw
        cwy_out[0] = cwy

        return idx


cdef class PytestWeightedFenwickTree(WeightedFenwickTree):
    """Used for testing only"""

    def py_reset(self, intp_t n):
        self.reset(n)

    def py_add(self, intp_t idx, float64_t y, float64_t w):
        self.add(idx, y, w)

    def py_search(self, float64_t t, inclusive=True):
        cdef float64_t w, wy
        idx = self.search(t, &w, &wy, inclusive)
        return idx, w, wy
