# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""This module contains njitted routines for building histograms.

A histogram is an array with n_bins entry of type HISTOGRAM_DTYPE. Each
feature has its own histogram. A histogram contains the sum of gradients and
hessians of all the samples belonging to each bin.
"""
cimport cython

import numpy as np
cimport numpy as np

from .types import HISTOGRAM_DTYPE


ctypedef np.npy_uint8 NPY_X_BINNED_DTYPE
ctypedef np.npy_float32 NPY_Y_DTYPE

cdef packed struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count


cpdef void _build_histogram_naive(unsigned int n_bins,
                                  unsigned int [:] sample_indices,
                                  NPY_X_BINNED_DTYPE [:] binned_feature,
                                  NPY_Y_DTYPE [:] ordered_gradients,
                                  NPY_Y_DTYPE [:] ordered_hessians,
                                  hist_struct [:] out) nogil:
    """Build histogram in a naive way, without optimizing for cache hit."""
    cdef:
        unsigned int i
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int sample_idx
        unsigned int bin_idx

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        bin_idx = binned_feature[sample_idx]
        out[bin_idx].sum_gradients += ordered_gradients[i]
        out[bin_idx].sum_hessians += ordered_hessians[i]
        out[bin_idx].count += 1


cpdef void _subtract_histograms(unsigned int n_bins,
                                hist_struct [:] hist_a,
                                hist_struct [:] hist_b,
                                hist_struct [:] out) nogil:
    """Return hist_a - hist_b"""

    cdef:
        unsigned int i = 0
    for i in range(n_bins):
        out[i].sum_gradients = hist_a[i].sum_gradients - hist_b[i].sum_gradients
        out[i].sum_hessians = hist_a[i].sum_hessians - hist_b[i].sum_hessians
        out[i].count = hist_a[i].count - hist_b[i].count


cpdef void _build_histogram(unsigned int n_bins,
                            unsigned int [:] sample_indices,
                            NPY_X_BINNED_DTYPE [:] binned_feature,
                            NPY_Y_DTYPE [:] ordered_gradients,
                            NPY_Y_DTYPE [:] ordered_hessians,
                            hist_struct [:] out) nogil:
    """Return histogram for a given feature."""
    cdef:
        unsigned int i = 0
        unsigned int n_node_samples = sample_indices.shape[0]
        unsigned int unrolled_upper = (n_node_samples // 4) * 4

        unsigned int bin_0
        unsigned int bin_1
        unsigned int bin_2
        unsigned int bin_3
        unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        out[bin_0].sum_gradients += ordered_gradients[i]
        out[bin_1].sum_gradients += ordered_gradients[i + 1]
        out[bin_2].sum_gradients += ordered_gradients[i + 2]
        out[bin_3].sum_gradients += ordered_gradients[i + 3]

        out[bin_0].sum_hessians += ordered_hessians[i]
        out[bin_1].sum_hessians += ordered_hessians[i + 1]
        out[bin_2].sum_hessians += ordered_hessians[i + 2]
        out[bin_3].sum_hessians += ordered_hessians[i + 3]

        out[bin_0].count += 1
        out[bin_1].count += 1
        out[bin_2].count += 1
        out[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        out[bin_idx].sum_gradients += ordered_gradients[i]
        out[bin_idx].sum_hessians += ordered_hessians[i]
        out[bin_idx].count += 1


cpdef void _build_histogram_no_hessian(
    unsigned int n_bins,
    unsigned int [:] sample_indices,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] ordered_gradients,
    hist_struct [:] out) nogil:
    """Return histogram for a given feature."""
    cdef:
        unsigned int i = 0
        unsigned int n_node_samples = sample_indices.shape[0]
        unsigned int unrolled_upper = (n_node_samples // 4) * 4

        unsigned int bin_0
        unsigned int bin_1
        unsigned int bin_2
        unsigned int bin_3
        unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        out[bin_0].sum_gradients += ordered_gradients[i]
        out[bin_1].sum_gradients += ordered_gradients[i + 1]
        out[bin_2].sum_gradients += ordered_gradients[i + 2]
        out[bin_3].sum_gradients += ordered_gradients[i + 3]

        out[bin_0].count += 1
        out[bin_1].count += 1
        out[bin_2].count += 1
        out[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        out[bin_idx].sum_gradients += ordered_gradients[i]
        out[bin_idx].count += 1


cpdef void _build_histogram_root_no_hessian(
    unsigned int n_bins,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] all_gradients,
    hist_struct [:] out) nogil:
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients already have a consistent
    ordering.

    Hessians are not updated (used when hessians are constant)
    """
    cdef:
        unsigned int i = 0
        unsigned int n_samples = binned_feature.shape[0]
        unsigned int unrolled_upper = (n_samples // 4) * 4

        unsigned int bin_0
        unsigned int bin_1
        unsigned int bin_2
        unsigned int bin_3
        unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        out[bin_0].sum_gradients += all_gradients[i]
        out[bin_1].sum_gradients += all_gradients[i + 1]
        out[bin_2].sum_gradients += all_gradients[i + 2]
        out[bin_3].sum_gradients += all_gradients[i + 3]

        out[bin_0].count += 1
        out[bin_1].count += 1
        out[bin_2].count += 1
        out[bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        out[bin_idx].sum_gradients += all_gradients[i]
        out[bin_idx].count += 1


cpdef void _build_histogram_root(
    unsigned int n_bins,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] all_gradients,
    NPY_Y_DTYPE [:] all_hessians,
    hist_struct [:] out) nogil:
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients and all_hessians already
    have a consistent ordering.
    """
    cdef:
        unsigned int i = 0
        unsigned int n_samples = binned_feature.shape[0]
        unsigned int unrolled_upper = (n_samples // 4) * 4

        unsigned int bin_0
        unsigned int bin_1
        unsigned int bin_2
        unsigned int bin_3
        unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):

        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        out[bin_0].sum_gradients += all_gradients[i]
        out[bin_1].sum_gradients += all_gradients[i + 1]
        out[bin_2].sum_gradients += all_gradients[i + 2]
        out[bin_3].sum_gradients += all_gradients[i + 3]

        out[bin_0].sum_hessians += all_hessians[i]
        out[bin_1].sum_hessians += all_hessians[i + 1]
        out[bin_2].sum_hessians += all_hessians[i + 2]
        out[bin_3].sum_hessians += all_hessians[i + 3]

        out[bin_0].count += 1
        out[bin_1].count += 1
        out[bin_2].count += 1
        out[bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        out[bin_idx].sum_gradients += all_gradients[i]
        out[bin_idx].sum_hessians += all_hessians[i]
        out[bin_idx].count += 1
