# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython.parallel import prange
from libc.math cimport floor
import numpy as np

from sklearn.ensemble._hist_gradient_boosting.common import Y_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common cimport Y_DTYPE_C


cdef extern from "<algorithm>" namespace "std" nogil:
    void nth_element[Iter](Iter first, Iter nth, Iter last)


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


def _update_raw_predictions_vector(
        Y_DTYPE_C [::1, :] raw_predictions,  # OUT
        grower,
        n_threads,
):
    """Update raw_predictions with the predictions of the newest vector tree.

    Vector leaves counterpart of `_update_raw_predictions`: raw_predictions
    has shape (n_samples, n_outputs) and each leaf value is an array of shape
    (n_outputs,).
    """
    cdef:
        unsigned int [::1] starts  # start of each leaf in partition
        unsigned int [::1] stops  # end of each leaf in partition
        Y_DTYPE_C [:, ::1] values  # values of each leaf
        const unsigned int [::1] partition = grower.splitter.partition
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.partition_start for leaf in leaves],
                      dtype=np.uint32)
    stops = np.array([leaf.partition_stop for leaf in leaves],
                     dtype=np.uint32)
    values = np.array([leaf.value for leaf in leaves], dtype=Y_DTYPE)
    _update_raw_predictions_vector_helper(raw_predictions, starts, stops,
                                          partition, values, n_threads)


cdef inline void _update_raw_predictions_vector_helper(
        Y_DTYPE_C [::1, :] raw_predictions,  # OUT
        const unsigned int [::1] starts,
        const unsigned int [::1] stops,
        const unsigned int [::1] partition,
        const Y_DTYPE_C [:, ::1] values,
        int n_threads,
):

    cdef:
        unsigned int position
        unsigned int sample_idx
        unsigned int k
        unsigned int n_outputs = values.shape[1]
        int leaf_idx
        int n_leaves = starts.shape[0]

    for leaf_idx in prange(n_leaves, schedule='static', nogil=True,
                           num_threads=n_threads):
        for position in range(starts[leaf_idx], stops[leaf_idx]):
            sample_idx = partition[position]
            for k in range(n_outputs):
                raw_predictions[sample_idx, k] += values[leaf_idx, k]


def _compute_quantile_leaves_values(
        const Y_DTYPE_C [::1, :] y_true,
        const Y_DTYPE_C [::1, :] raw_predictions,
        grower,
        const Y_DTYPE_C [::1] quantiles,
        n_threads,
):
    """Compute the unweighted quantiles of the residuals in each leaf.

    For each leaf and each output k, this computes the quantile of level
    quantiles[k] of y_true[:, k] - raw_predictions[:, k] over the samples of
    the leaf, with the same (linear) interpolation as np.percentile.

    y_true has either one column per output or a single column shared by all
    the outputs. Returns an array of shape (n_leaves, n_outputs).
    """
    cdef:
        unsigned int [::1] starts
        unsigned int [::1] stops
        const unsigned int [::1] partition = grower.splitter.partition
        int n_outputs = raw_predictions.shape[1]
        Y_DTYPE_C [:, ::1] buffer
        Y_DTYPE_C [:, ::1] values
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.partition_start for leaf in leaves],
                      dtype=np.uint32)
    stops = np.array([leaf.partition_stop for leaf in leaves],
                     dtype=np.uint32)
    # Leaves cover disjoint ranges of the partition: each (leaf, output) task
    # writes its residuals to buffer[output, start:stop] without conflict.
    buffer = np.empty((n_outputs, partition.shape[0]), dtype=Y_DTYPE)
    values = np.zeros((len(leaves), n_outputs), dtype=Y_DTYPE)
    _compute_quantile_leaves_values_helper(
        y_true, raw_predictions, starts, stops, partition, quantiles, buffer,
        values, n_threads
    )
    return np.asarray(values)


cdef void _compute_quantile_leaves_values_helper(
        const Y_DTYPE_C [::1, :] y_true,
        const Y_DTYPE_C [::1, :] raw_predictions,
        const unsigned int [::1] starts,
        const unsigned int [::1] stops,
        const unsigned int [::1] partition,
        const Y_DTYPE_C [::1] quantiles,
        Y_DTYPE_C [:, ::1] buffer,
        Y_DTYPE_C [:, ::1] values,  # OUT
        int n_threads,
) noexcept nogil:
    cdef:
        int n_outputs = raw_predictions.shape[1]
        int n_leaves = starts.shape[0]
        bint shared_target = y_true.shape[1] == 1
        int task_idx
        int leaf_idx
        int k
        int k_y
        unsigned int position
        unsigned int sample_idx
        Py_ssize_t n
        Py_ssize_t i
        Py_ssize_t lo
        Y_DTYPE_C *b
        Y_DTYPE_C q
        Y_DTYPE_C virtual_idx
        Y_DTYPE_C gamma
        Y_DTYPE_C a_val
        Y_DTYPE_C b_val
        Y_DTYPE_C diff

    for task_idx in prange(n_leaves * n_outputs, schedule='dynamic',
                           num_threads=n_threads):
        leaf_idx = task_idx // n_outputs
        k = task_idx % n_outputs
        k_y = 0 if shared_target else k
        n = stops[leaf_idx] - starts[leaf_idx]
        if n == 0:
            continue
        b = &buffer[k, starts[leaf_idx]]
        for position in range(starts[leaf_idx], stops[leaf_idx]):
            sample_idx = partition[position]
            b[position - starts[leaf_idx]] = (
                y_true[sample_idx, k_y] - raw_predictions[sample_idx, k]
            )

        # Same virtual index and interpolation as np.percentile with
        # method="linear".
        q = quantiles[k]
        virtual_idx = (n * q + (1. - q)) - 1.
        lo = <Py_ssize_t> floor(virtual_idx)
        if lo < 0:
            lo = 0
        if lo > n - 1:
            lo = n - 1
        gamma = virtual_idx - lo
        nth_element(b, b + lo, b + n)
        a_val = b[lo]
        if lo + 1 < n:
            # After nth_element, the next order statistic is the minimum of
            # the upper part.
            b_val = b[lo + 1]
            for i in range(lo + 2, n):
                if b[i] < b_val:
                    b_val = b[i]
        else:
            b_val = a_val
        diff = b_val - a_val
        if gamma >= 0.5:
            values[leaf_idx, k] = b_val - diff * (1. - gamma)
        else:
            values[leaf_idx, k] = a_val + diff * gamma
