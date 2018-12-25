"""This module contains njitted routines for building histograms.

A histogram is an array with n_bins entry of type HISTOGRAM_DTYPE. Each
feature has its own histogram. A histogram contains the sum of gradients and
hessians of all the samples belonging to each bin.
"""
cimport cython

import numpy as np
cimport numpy as np


HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', np.float32),
    ('sum_hessians', np.float32),
    ('count', np.uint32),
])


from libc.stdlib cimport malloc, free

cdef struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _build_histogram_naive(unsigned int n_bins, unsigned int [:]
                                sample_indices, unsigned char [:]
                                binned_feature, float [:] ordered_gradients,
                                float[:] ordered_hessians):
    """Build histogram in a naive way, without optimizing for cache hit."""
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef:
        hist_struct [:] view = histogram
        unsigned int i
        unsigned int sample_idx
        unsigned char bin_idx

    for i, sample_idx in enumerate(sample_indices):
        bin_idx = binned_feature[sample_idx]
        view[bin_idx].sum_gradients += ordered_gradients[i]
        view[bin_idx].sum_hessians += ordered_hessians[i]
        view[bin_idx].count += 1
    return histogram


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _subtract_histograms(unsigned int n_bins, np.ndarray hist_a, np.ndarray hist_b):
    """Return hist_a - hist_b"""
    # print('subtract_hist')

    cdef unsigned int i = 0
    cdef np.ndarray histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef hist_struct [:] view = histogram
    cdef hist_struct [:] view_a = hist_a
    cdef hist_struct [:] view_b = hist_b

    for i in range(n_bins):
        view[i].sum_gradients = view_a[i].sum_gradients - view_b[i].sum_gradients
        view[i].sum_hessians = view_a[i].sum_hessians - view_b[i].sum_hessians
        view[i].count = view_a[i].count - view_b[i].count

    return histogram


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _build_histogram(unsigned int n_bins, unsigned int [:]
                                sample_indices, unsigned char [:]
                                binned_feature, float [:] ordered_gradients,
                                float[:] ordered_hessians):
    """Return histogram for a given feature."""
    cdef np.ndarray histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef hist_struct [:] view = histogram
    cdef int i = 0

    cdef float [:] ordered_gradients_view = ordered_gradients
    cdef float [:] ordered_hessians_view = ordered_hessians

    cdef int n_node_samples = sample_indices.shape[0]
    cdef int unrolled_upper = (n_node_samples // 4) * 4

    cdef unsigned int bin_0
    cdef unsigned int bin_1
    cdef unsigned int bin_2
    cdef unsigned int bin_3
    cdef unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        view[bin_0].sum_gradients += ordered_gradients_view[i]
        view[bin_1].sum_gradients += ordered_gradients_view[i + 1]
        view[bin_2].sum_gradients += ordered_gradients_view[i + 2]
        view[bin_3].sum_gradients += ordered_gradients_view[i + 3]

        view[bin_0].sum_hessians += ordered_hessians_view[i]
        view[bin_1].sum_hessians += ordered_hessians_view[i + 1]
        view[bin_2].sum_hessians += ordered_hessians_view[i + 2]
        view[bin_3].sum_hessians += ordered_hessians_view[i + 3]

        view[bin_0].count += 1
        view[bin_1].count += 1
        view[bin_2].count += 1
        view[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        view[bin_idx].sum_gradients += ordered_gradients_view[i]
        view[bin_idx].sum_hessians += ordered_hessians_view[i]
        view[bin_idx].count += 1

    return histogram


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _build_histogram_no_hessian(unsigned int n_bins, unsigned int [:]
                                sample_indices, unsigned char [:]
                                binned_feature, float [:] ordered_gradients):
    """Return histogram for a given feature.

    Hessians are not updated (used when hessians are constant).
    """
    # print('build_hist_no_hessian')
    cdef np.ndarray histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef hist_struct [:] view = histogram
    cdef unsigned int i = 0

    cdef float [:] ordered_gradients_view = ordered_gradients
    cdef unsigned char [:] binned_feature_view = binned_feature
    cdef unsigned int [:] sample_indices_view = sample_indices

    cdef unsigned int n_node_samples = sample_indices.shape[0]
    cdef unsigned int unrolled_upper = (n_node_samples // 4) * 4

    cdef unsigned int bin_0
    cdef unsigned int bin_1
    cdef unsigned int bin_2
    cdef unsigned int bin_3
    cdef unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature_view[sample_indices_view[i]]
        bin_1 = binned_feature_view[sample_indices_view[i + 1]]
        bin_2 = binned_feature_view[sample_indices_view[i + 2]]
        bin_3 = binned_feature_view[sample_indices_view[i + 3]]

        view[bin_0].sum_gradients += ordered_gradients_view[i]
        view[bin_1].sum_gradients += ordered_gradients_view[i + 1]
        view[bin_2].sum_gradients += ordered_gradients_view[i + 2]
        view[bin_3].sum_gradients += ordered_gradients_view[i + 3]

        view[bin_0].count += 1
        view[bin_1].count += 1
        view[bin_2].count += 1
        view[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature_view[sample_indices_view[i]]
        view[bin_idx].sum_gradients += ordered_gradients_view[i]
        view[bin_idx].count += 1

    return histogram



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _build_histogram_root_no_hessian(unsigned int n_bins, unsigned char [:]
                                     binned_feature, float [:]all_gradients):
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients already have a consistent
    ordering.

    Hessians are not updated (used when hessians are constant)
    """
    # print('build_hist_root_no_hessian')

    cdef np.ndarray histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef hist_struct [:] view = histogram
    cdef unsigned int i = 0

    cdef float [:] all_gradients_view = all_gradients
    cdef unsigned char [:] binned_feature_view = binned_feature

    cdef unsigned int n_node_samples = binned_feature.shape[0]
    cdef unsigned int unrolled_upper = (n_node_samples // 4) * 4

    cdef unsigned int bin_0
    cdef unsigned int bin_1
    cdef unsigned int bin_2
    cdef unsigned int bin_3
    cdef unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature_view[i]
        bin_1 = binned_feature_view[i + 1]
        bin_2 = binned_feature_view[i + 2]
        bin_3 = binned_feature_view[i + 3]

        view[bin_0].sum_gradients += all_gradients_view[i]
        view[bin_1].sum_gradients += all_gradients_view[i + 1]
        view[bin_2].sum_gradients += all_gradients_view[i + 2]
        view[bin_3].sum_gradients += all_gradients_view[i + 3]

        view[bin_0].count += 1
        view[bin_1].count += 1
        view[bin_2].count += 1
        view[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature_view[i]
        view[bin_idx].sum_gradients += all_gradients_view[i]
        view[bin_idx].count += 1

    return histogram


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _build_histogram_root(unsigned int n_bins, unsigned char [:]
                          binned_feature, float [:] all_gradients,
                          float[:] all_hessians):
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients and all_hessians already
    have a consistent ordering.
    """
    cdef np.ndarray histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    cdef hist_struct [:] view = histogram
    cdef int i = 0

    cdef unsigned int n_node_samples = binned_feature.shape[0]
    cdef unsigned int unrolled_upper = (n_node_samples // 4) * 4

    cdef unsigned int bin_0
    cdef unsigned int bin_1
    cdef unsigned int bin_2
    cdef unsigned int bin_3
    cdef unsigned int bin_idx

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        view[bin_0].sum_gradients += all_gradients[i]
        view[bin_1].sum_gradients += all_gradients[i + 1]
        view[bin_2].sum_gradients += all_gradients[i + 2]
        view[bin_3].sum_gradients += all_gradients[i + 3]

        view[bin_0].sum_hessians += all_hessians[i]
        view[bin_1].sum_hessians += all_hessians[i + 1]
        view[bin_2].sum_hessians += all_hessians[i + 2]
        view[bin_3].sum_hessians += all_hessians[i + 3]

        view[bin_0].count += 1
        view[bin_1].count += 1
        view[bin_2].count += 1
        view[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[i]
        view[bin_idx].sum_gradients += all_gradients[i]
        view[bin_idx].sum_hessians += all_hessians[i]
        view[bin_idx].count += 1

    return histogram
