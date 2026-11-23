# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython.parallel import prange
import numpy as np

from sklearn.ensemble._hist_gradient_boosting.common import Y_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common cimport Y_DTYPE_C


def _update_raw_predictions(
        Y_DTYPE_C [::1] raw_predictions,  # OUT
        grower,
        n_threads,
):
    """Update raw_predictions with the predictions of the newest tree.

    This is equivalent to (and much faster than):
        raw_predictions += last_estimator.predict(X_train)

    It's only possible for data X_train that is used to train the trees (it
    isn't usable for e.g. X_val).
    """
    cdef:
        unsigned int [::1] starts  # start of each leaf in partition
        unsigned int [::1] stops  # end of each leaf in partition
        Y_DTYPE_C [::1] values  # value of each leaf
        const unsigned int [::1] partition = grower.splitter.partition
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.partition_start for leaf in leaves],
                      dtype=np.uint32)
    stops = np.array([leaf.partition_stop for leaf in leaves],
                     dtype=np.uint32)
    values = np.array([leaf.value for leaf in leaves], dtype=Y_DTYPE)

    _update_raw_predictions_helper(raw_predictions, starts, stops, partition,
                                   values, n_threads)


cdef inline void _update_raw_predictions_helper(
        Y_DTYPE_C [::1] raw_predictions,  # OUT
        const unsigned int [::1] starts,
        const unsigned int [::1] stops,
        const unsigned int [::1] partition,
        const Y_DTYPE_C [::1] values,
        int n_threads,
):

    cdef:
        unsigned int position
        int leaf_idx
        int n_leaves = starts.shape[0]

    for leaf_idx in prange(n_leaves, schedule='static', nogil=True,
                           num_threads=n_threads):
        for position in range(starts[leaf_idx], stops[leaf_idx]):
            raw_predictions[partition[position]] += values[leaf_idx]


def _update_leaves_values_quantile(
        grower,
        const Y_DTYPE_C [:] y_true,
        const Y_DTYPE_C [:] raw_prediction,
        double quantile,
        bint is_median,
        int n_threads,
):
    """Set each leaf value to the quantile of its residuals, times shrinkage.

    For every finalized leaf of `grower`, this computes the quantile of
    `y_true - raw_prediction` over the samples in that leaf and stores
    `grower.shrinkage * quantile` in `leaf.value`. Leaves are processed in
    parallel.

    The result is identical to `np.median` (if `is_median`) or to
    `np.percentile(..., 100 * quantile)` with the default "linear" method,
    which is what `AbsoluteError.fit_intercept_only` and
    `PinballLoss.fit_intercept_only` compute without sample weights.
    """
    cdef:
        list leaves = grower.finalized_leaves
        const unsigned int [::1] partition = grower.splitter.partition
        unsigned int [::1] starts = np.array(
            [leaf.partition_start for leaf in leaves], dtype=np.uint32
        )
        unsigned int [::1] stops = np.array(
            [leaf.partition_stop for leaf in leaves], dtype=np.uint32
        )
        Y_DTYPE_C [::1] values = np.empty(len(leaves), dtype=Y_DTYPE)
        # Each leaf owns the disjoint slice [start, stop) of this buffer.
        Y_DTYPE_C [::1] buffer = np.empty(partition.shape[0], dtype=Y_DTYPE)
        int leaf_idx
        int n_leaves = len(leaves)
        unsigned int position, start, n
        # Same rounding as np.percentile(a, 100 * quantile), which divides by 100.
        double q = (100.0 * quantile) / 100.0

    for leaf_idx in prange(n_leaves, schedule='dynamic', nogil=True,
                           num_threads=n_threads):
        start = starts[leaf_idx]
        n = stops[leaf_idx] - start
        for position in range(start, stops[leaf_idx]):
            buffer[position] = (
                y_true[partition[position]] - raw_prediction[partition[position]]
            )
        if is_median:
            values[leaf_idx] = _median(&buffer[start], n)
        else:
            values[leaf_idx] = _percentile_linear(&buffer[start], n, q)

    for leaf_idx in range(n_leaves):
        leaves[leaf_idx].value = grower.shrinkage * values[leaf_idx]


cdef inline void _swap(Y_DTYPE_C* a, Y_DTYPE_C* b) noexcept nogil:
    cdef Y_DTYPE_C tmp = a[0]
    a[0] = b[0]
    b[0] = tmp


cdef Y_DTYPE_C _select(Y_DTYPE_C* arr, Py_ssize_t n, Py_ssize_t k) noexcept nogil:
    """Return the k-th smallest element of arr[:n], partially reordering arr.

    Quickselect with median-of-three pivot (Hoare partitioning), which copes
    well with many tied values. On return, arr[:k] <= arr[k] <= arr[k+1:n].
    """
    cdef:
        Py_ssize_t left = 0, right = n - 1, mid, i, j
        Y_DTYPE_C pivot

    while True:
        if right <= left + 1:
            if right == left + 1 and arr[right] < arr[left]:
                _swap(&arr[left], &arr[right])
            return arr[k]
        mid = (left + right) >> 1
        _swap(&arr[mid], &arr[left + 1])
        if arr[left] > arr[right]:
            _swap(&arr[left], &arr[right])
        if arr[left + 1] > arr[right]:
            _swap(&arr[left + 1], &arr[right])
        if arr[left] > arr[left + 1]:
            _swap(&arr[left], &arr[left + 1])
        i = left + 1
        j = right
        pivot = arr[left + 1]
        while True:
            i = i + 1
            while arr[i] < pivot:
                i = i + 1
            j = j - 1
            while arr[j] > pivot:
                j = j - 1
            if j < i:
                break
            _swap(&arr[i], &arr[j])
        arr[left + 1] = arr[j]
        arr[j] = pivot
        if j >= k:
            right = j - 1
        if j <= k:
            left = i


cdef inline Y_DTYPE_C _min(const Y_DTYPE_C* arr, Py_ssize_t n) noexcept nogil:
    cdef:
        Py_ssize_t i
        Y_DTYPE_C res = arr[0]
    for i in range(1, n):
        if arr[i] < res:
            res = arr[i]
    return res


cdef Y_DTYPE_C _median(Y_DTYPE_C* arr, Py_ssize_t n) noexcept nogil:
    """Same result as np.median(arr[:n])."""
    cdef:
        Py_ssize_t half = n // 2
        Y_DTYPE_C low
    if n % 2 == 1:
        return _select(arr, n, half)
    low = _select(arr, n, half - 1)
    # After selection, the next order statistic is the minimum above index k.
    return (low + _min(&arr[half], n - half)) / 2


cdef Y_DTYPE_C _percentile_linear(
    Y_DTYPE_C* arr, Py_ssize_t n, double q
) noexcept nogil:
    """Same result as np.quantile(arr[:n], q, method="linear")."""
    cdef:
        double virtual_index = (n - 1) * q
        double gamma
        Py_ssize_t k
        Y_DTYPE_C a, b, diff_b_a

    if virtual_index >= n - 1:
        return _select(arr, n, n - 1)
    k = <Py_ssize_t> virtual_index  # floor, as virtual_index >= 0
    gamma = virtual_index - k
    a = _select(arr, n, k)
    b = _min(&arr[k + 1], n - k - 1)
    # Same arithmetic as numpy's _lerp to get bitwise identical results.
    diff_b_a = b - a
    if gamma >= 0.5:
        return b - diff_b_a * (1 - gamma)
    return a + diff_b_a * gamma
