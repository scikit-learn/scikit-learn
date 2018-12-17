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


def _build_histogram_naive(n_bins, sample_indices, binned_feature,
                           ordered_gradients, ordered_hessians):
    """Build histogram in a naive way, without optimizing for cache hit."""
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    for i, sample_idx in enumerate(sample_indices):
        bin_idx = binned_feature[sample_idx]
        histogram[bin_idx]['sum_gradients'] += ordered_gradients[i]
        histogram[bin_idx]['sum_hessians'] += ordered_hessians[i]
        histogram[bin_idx]['count'] += 1
    return histogram


def _subtract_histograms(n_bins, hist_a, hist_b):
    """Return hist_a - hist_b"""

    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)

    sg = 'sum_gradients'
    sh = 'sum_hessians'
    c = 'count'

    for i in range(n_bins):
        histogram[i][sg] = hist_a[i][sg] - hist_b[i][sg]
        histogram[i][sh] = hist_a[i][sh] - hist_b[i][sh]
        histogram[i][c] = hist_a[i][c] - hist_b[i][c]

    return histogram


def _build_histogram(n_bins, sample_indices, binned_feature, ordered_gradients,
                     ordered_hessians):
    """Return histogram for a given feature."""
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    n_node_samples = sample_indices.shape[0]
    unrolled_upper = (n_node_samples // 4) * 4

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        histogram[bin_0]['sum_gradients'] += ordered_gradients[i]
        histogram[bin_1]['sum_gradients'] += ordered_gradients[i + 1]
        histogram[bin_2]['sum_gradients'] += ordered_gradients[i + 2]
        histogram[bin_3]['sum_gradients'] += ordered_gradients[i + 3]

        histogram[bin_0]['sum_hessians'] += ordered_hessians[i]
        histogram[bin_1]['sum_hessians'] += ordered_hessians[i + 1]
        histogram[bin_2]['sum_hessians'] += ordered_hessians[i + 2]
        histogram[bin_3]['sum_hessians'] += ordered_hessians[i + 3]

        histogram[bin_0]['count'] += 1
        histogram[bin_1]['count'] += 1
        histogram[bin_2]['count'] += 1
        histogram[bin_3]['count'] += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        histogram[bin_idx]['sum_gradients'] += ordered_gradients[i]
        histogram[bin_idx]['sum_hessians'] += ordered_hessians[i]
        histogram[bin_idx]['count'] += 1

    return histogram


def _build_histogram_no_hessian(n_bins, sample_indices, binned_feature,
                                ordered_gradients):
    """Return histogram for a given feature.

    Hessians are not updated (used when hessians are constant).
    """
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    n_node_samples = sample_indices.shape[0]
    unrolled_upper = (n_node_samples // 4) * 4

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        histogram[bin_0]['sum_gradients'] += ordered_gradients[i]
        histogram[bin_1]['sum_gradients'] += ordered_gradients[i + 1]
        histogram[bin_2]['sum_gradients'] += ordered_gradients[i + 2]
        histogram[bin_3]['sum_gradients'] += ordered_gradients[i + 3]

        histogram[bin_0]['count'] += 1
        histogram[bin_1]['count'] += 1
        histogram[bin_2]['count'] += 1
        histogram[bin_3]['count'] += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        histogram[bin_idx]['sum_gradients'] += ordered_gradients[i]
        histogram[bin_idx]['count'] += 1

    return histogram


def _build_histogram_root_no_hessian(n_bins, binned_feature, all_gradients):
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients already have a consistent
    ordering.

    Hessians are not updated (used when hessians are constant)
    """
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    n_node_samples = binned_feature.shape[0]
    unrolled_upper = (n_node_samples // 4) * 4

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        histogram[bin_0]['sum_gradients'] += all_gradients[i]
        histogram[bin_1]['sum_gradients'] += all_gradients[i + 1]
        histogram[bin_2]['sum_gradients'] += all_gradients[i + 2]
        histogram[bin_3]['sum_gradients'] += all_gradients[i + 3]

        histogram[bin_0]['count'] += 1
        histogram[bin_1]['count'] += 1
        histogram[bin_2]['count'] += 1
        histogram[bin_3]['count'] += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[i]
        histogram[bin_idx]['sum_gradients'] += all_gradients[i]
        histogram[bin_idx]['count'] += 1

    return histogram


def _build_histogram_root(n_bins, binned_feature, all_gradients,
                          all_hessians):
    """Special case for the root node

    The root node has to find the split among all the samples from the
    training set. binned_feature and all_gradients and all_hessians already
    have a consistent ordering.
    """
    histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
    n_node_samples = binned_feature.shape[0]
    unrolled_upper = (n_node_samples // 4) * 4

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        histogram[bin_0]['sum_gradients'] += all_gradients[i]
        histogram[bin_1]['sum_gradients'] += all_gradients[i + 1]
        histogram[bin_2]['sum_gradients'] += all_gradients[i + 2]
        histogram[bin_3]['sum_gradients'] += all_gradients[i + 3]

        histogram[bin_0]['sum_hessians'] += all_hessians[i]
        histogram[bin_1]['sum_hessians'] += all_hessians[i + 1]
        histogram[bin_2]['sum_hessians'] += all_hessians[i + 2]
        histogram[bin_3]['sum_hessians'] += all_hessians[i + 3]

        histogram[bin_0]['count'] += 1
        histogram[bin_1]['count'] += 1
        histogram[bin_2]['count'] += 1
        histogram[bin_3]['count'] += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[i]
        histogram[bin_idx]['sum_gradients'] += all_gradients[i]
        histogram[bin_idx]['sum_hessians'] += all_hessians[i]
        histogram[bin_idx]['count'] += 1

    return histogram
