"""This module contains routines for building histograms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

cimport cython
from cython.parallel import prange
from libc.stdint cimport int32_t, int64_t
from libc.string cimport memset

import numpy as np

from sklearn.ensemble._hist_gradient_boosting.common import HISTOGRAM_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common cimport hist_struct
from sklearn.ensemble._hist_gradient_boosting.common cimport X_BINNED_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport G_H_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport Y_DTYPE_C
from sklearn.utils._typedefs cimport uint8_t


ctypedef fused sparse_index_t:
    int32_t
    int64_t


cdef void _build_sparse_histograms_impl(
    const X_BINNED_DTYPE_C [::1] data,
    const sparse_index_t [::1] indices,
    const sparse_index_t [::1] indptr,
    const uint8_t [::1] zero_bin_idx,
    const uint8_t [::1] feature_mask,
    const unsigned int [::1] sample_indices,
    const G_H_DTYPE_C [::1] gradients,
    const G_H_DTYPE_C [::1] hessians,
    const uint8_t hessians_are_constant,
    const uint8_t missing_values_bin_idx,
    const unsigned int n_features,
    const unsigned int n_bins,
    const int n_threads,
    hist_struct [:, ::1] histograms,
) noexcept nogil:
    cdef:
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int i
        unsigned int sample_idx
        unsigned int feature_idx
        unsigned int bin_idx
        sparse_index_t j
        sparse_index_t start
        sparse_index_t stop
        unsigned int zero_bin
        unsigned int other_count
        unsigned int bin
        X_BINNED_DTYPE_C encoded_bin
        Y_DTYPE_C total_sum_gradients = 0
        Y_DTYPE_C total_sum_hessians = 0
        Y_DTYPE_C other_sum_gradients
        Y_DTYPE_C other_sum_hessians

    memset(&histograms[0, 0], 0, n_features * n_bins * sizeof(hist_struct))

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        total_sum_gradients += gradients[sample_idx]
        if hessians_are_constant:
            total_sum_hessians += hessians[0]
        else:
            total_sum_hessians += hessians[sample_idx]

        start = indptr[sample_idx]
        stop = indptr[sample_idx + 1]
        for j in range(start, stop):
            feature_idx = indices[j]
            if not feature_mask[feature_idx]:
                continue
            encoded_bin = data[j]
            if encoded_bin == missing_values_bin_idx:
                bin_idx = missing_values_bin_idx
            elif encoded_bin == 0:
                bin_idx = zero_bin_idx[feature_idx]
            elif encoded_bin <= zero_bin_idx[feature_idx]:
                bin_idx = encoded_bin - 1
            else:
                bin_idx = encoded_bin
            histograms[feature_idx, bin_idx].sum_gradients += gradients[sample_idx]
            if hessians_are_constant:
                histograms[feature_idx, bin_idx].sum_hessians += hessians[0]
            else:
                histograms[feature_idx, bin_idx].sum_hessians += hessians[sample_idx]
            histograms[feature_idx, bin_idx].count += 1

    for feature_idx in range(n_features):
        if not feature_mask[feature_idx]:
            continue
        zero_bin = zero_bin_idx[feature_idx]
        other_sum_gradients = 0
        other_sum_hessians = 0
        other_count = 0
        for bin in range(n_bins):
            if bin != zero_bin:
                other_sum_gradients += histograms[feature_idx, bin].sum_gradients
                other_sum_hessians += histograms[feature_idx, bin].sum_hessians
                other_count += histograms[feature_idx, bin].count
        histograms[feature_idx, zero_bin].sum_gradients = (
            total_sum_gradients - other_sum_gradients
        )
        histograms[feature_idx, zero_bin].sum_hessians = (
            total_sum_hessians - other_sum_hessians
        )
        histograms[feature_idx, zero_bin].count = n_samples - other_count


cdef void _build_sparse_histograms32(
    const X_BINNED_DTYPE_C [::1] data,
    const int32_t [::1] indices,
    const int32_t [::1] indptr,
    const uint8_t [::1] zero_bin_idx,
    const uint8_t [::1] feature_mask,
    const unsigned int [::1] sample_indices,
    const G_H_DTYPE_C [::1] gradients,
    const G_H_DTYPE_C [::1] hessians,
    const uint8_t hessians_are_constant,
    const uint8_t missing_values_bin_idx,
    const unsigned int n_features,
    const unsigned int n_bins,
    const int n_threads,
    hist_struct [:, ::1] histograms,
) noexcept nogil:
    _build_sparse_histograms_impl[int32_t](
        data, indices, indptr, zero_bin_idx, feature_mask, sample_indices,
        gradients, hessians, hessians_are_constant, missing_values_bin_idx,
        n_features, n_bins, n_threads, histograms
    )


cdef void _build_sparse_histograms64(
    const X_BINNED_DTYPE_C [::1] data,
    const int64_t [::1] indices,
    const int64_t [::1] indptr,
    const uint8_t [::1] zero_bin_idx,
    const uint8_t [::1] feature_mask,
    const unsigned int [::1] sample_indices,
    const G_H_DTYPE_C [::1] gradients,
    const G_H_DTYPE_C [::1] hessians,
    const uint8_t hessians_are_constant,
    const uint8_t missing_values_bin_idx,
    const unsigned int n_features,
    const unsigned int n_bins,
    const int n_threads,
    hist_struct [:, ::1] histograms,
) noexcept nogil:
    _build_sparse_histograms_impl[int64_t](
        data, indices, indptr, zero_bin_idx, feature_mask, sample_indices,
        gradients, hessians, hessians_are_constant, missing_values_bin_idx,
        n_features, n_bins, n_threads, histograms
    )


# Notes:
# - IN views are read-only, OUT views are write-only
# - In a lot of functions here, we pass feature_idx and the whole 2d
#   histograms arrays instead of just histograms[feature_idx]. This is because
#   Cython generated C code will have strange Python interactions (likely
#   related to the GIL release and the custom histogram dtype) when using 1d
#   histogram arrays that come from 2d arrays.
# - The for loops are un-wrapped, for example:
#
#   for i in range(n):
#       array[i] = i
#
#   will become
#
#   for i in range(n // 4):
#       array[i] = i
#       array[i + 1] = i + 1
#       array[i + 2] = i + 2
#       array[i + 3] = i + 3
#
#   This is to hint gcc that it can auto-vectorize these 4 operations and
#   perform them all at once.


@cython.final
cdef class HistogramBuilder:
    """A Histogram builder... used to build histograms.

    A histogram is an array with n_bins entries of type HISTOGRAM_DTYPE. Each
    feature has its own histogram. A histogram contains the sum of gradients
    and hessians of all the samples belonging to each bin.

    There are different ways to build a histogram:
    - by subtraction: hist(child) = hist(parent) - hist(sibling)
    - from scratch. In this case we have routines that update the hessians
      or not (not useful when hessians are constant for some losses e.g.
      least squares). Also, there's a special case for the root which
      contains all the samples, leading to some possible optimizations.
      Overall all the implementations look the same, and are optimized for
      cache hit.

    Parameters
    ----------
    X_binned : ndarray or CSR sparse matrix/array of int
        The binned input samples. Dense arrays must be Fortran-aligned. Sparse
        inputs use the internal CSR encoding where stored value 0 represents
        the bin containing the implicit numeric zero.
    n_bins : int
        The total number of bins, including the bin for missing values. Used
        to define the shape of the histograms.
    gradients : ndarray, shape (n_samples,)
        The gradients of each training sample. Those are the gradients of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    hessians : ndarray, shape (n_samples,)
        The hessians of each training sample. Those are the hessians of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    hessians_are_constant : bool
        Whether hessians are constant.
    """
    cdef public:
        object X_binned
        const X_BINNED_DTYPE_C [::1, :] X_binned_dense
        const X_BINNED_DTYPE_C [::1] X_binned_sparse_data
        const uint8_t [::1] zero_bin_idx
        unsigned int n_features
        unsigned int n_bins
        G_H_DTYPE_C [::1] gradients
        G_H_DTYPE_C [::1] hessians
        G_H_DTYPE_C [::1] ordered_gradients
        G_H_DTYPE_C [::1] ordered_hessians
        uint8_t hessians_are_constant
        uint8_t is_sparse
        int n_threads

    def __init__(self, X_binned,
                 unsigned int n_bins, G_H_DTYPE_C [::1] gradients,
                 G_H_DTYPE_C [::1] hessians,
                 uint8_t hessians_are_constant,
                 int n_threads,
                 zero_bin_idx=None):

        self.X_binned = X_binned
        self.is_sparse = hasattr(X_binned, "format") and X_binned.format == "csr"
        if self.is_sparse:
            self.X_binned_sparse_data = X_binned.data
            self.zero_bin_idx = zero_bin_idx
        else:
            self.X_binned_dense = X_binned
        self.n_features = X_binned.shape[1]
        # Note: all histograms will have <n_bins> bins, but some of the
        # bins may be unused if a feature has a small number of unique values.
        self.n_bins = n_bins
        self.gradients = gradients
        self.hessians = hessians
        # for root node, gradients and hessians are already ordered
        self.ordered_gradients = gradients.copy()
        self.ordered_hessians = hessians.copy()
        self.hessians_are_constant = hessians_are_constant
        self.n_threads = n_threads

    def compute_histograms_brute(
        HistogramBuilder self,
        const unsigned int [::1] sample_indices,       # IN
        const unsigned int [:] allowed_features=None,  # IN
    ):
        """Compute the histograms of the node by scanning through all the data.

        For a given feature, the complexity is O(n_samples)

        Parameters
        ----------
        sample_indices : array of int, shape (n_samples_at_node,)
            The indices of the samples at the node to split.

        allowed_features : None or ndarray, dtype=np.uint32
            Indices of the features that are allowed by interaction constraints to be
            split.

        Returns
        -------
        histograms : ndarray of HISTOGRAM_DTYPE, shape (n_features, n_bins)
            The computed histograms of the current node.
        """
        cdef:
            int n_samples
            int feature_idx
            int f_idx
            int i
            object feature_mask_array
            const uint8_t [::1] feature_mask
            # need local views to avoid python interactions
            uint8_t hessians_are_constant = self.hessians_are_constant
            int n_allowed_features = self.n_features
            G_H_DTYPE_C [::1] ordered_gradients = self.ordered_gradients
            G_H_DTYPE_C [::1] gradients = self.gradients
            G_H_DTYPE_C [::1] ordered_hessians = self.ordered_hessians
            G_H_DTYPE_C [::1] hessians = self.hessians
            # Histograms will be initialized to zero later within a prange
            hist_struct [:, ::1] histograms = np.empty(
                shape=(self.n_features, self.n_bins),
                dtype=HISTOGRAM_DTYPE
            )
            bint has_interaction_cst = allowed_features is not None
            int n_threads = self.n_threads

        if self.is_sparse:
            feature_mask_array = np.zeros(self.n_features, dtype=np.uint8)
            if has_interaction_cst:
                feature_mask_array[allowed_features] = 1
            else:
                feature_mask_array[:] = 1
            feature_mask = feature_mask_array
            if self.X_binned.indices.dtype == np.int32:
                _build_sparse_histograms32(
                    self.X_binned_sparse_data,
                    self.X_binned.indices,
                    self.X_binned.indptr,
                    self.zero_bin_idx,
                    feature_mask,
                    sample_indices,
                    gradients,
                    hessians,
                    hessians_are_constant,
                    self.n_bins - 1,
                    self.n_features,
                    self.n_bins,
                    n_threads,
                    histograms,
                )
            else:
                _build_sparse_histograms64(
                    self.X_binned_sparse_data,
                    self.X_binned.indices,
                    self.X_binned.indptr,
                    self.zero_bin_idx,
                    feature_mask,
                    sample_indices,
                    gradients,
                    hessians,
                    hessians_are_constant,
                    self.n_bins - 1,
                    self.n_features,
                    self.n_bins,
                    n_threads,
                    histograms,
                )
            return histograms

        if has_interaction_cst:
            n_allowed_features = allowed_features.shape[0]

        with nogil:
            n_samples = sample_indices.shape[0]

            # Populate ordered_gradients and ordered_hessians. (Already done
            # for root) Ordering the gradients and hessians helps to improve
            # cache hit.
            if sample_indices.shape[0] != gradients.shape[0]:
                if hessians_are_constant:
                    for i in prange(n_samples, schedule='static',
                                    num_threads=n_threads):
                        ordered_gradients[i] = gradients[sample_indices[i]]
                else:
                    for i in prange(n_samples, schedule='static',
                                    num_threads=n_threads):
                        ordered_gradients[i] = gradients[sample_indices[i]]
                        ordered_hessians[i] = hessians[sample_indices[i]]

            # Compute histogram of each feature
            for f_idx in prange(
                n_allowed_features, schedule='static', num_threads=n_threads
            ):
                if has_interaction_cst:
                    feature_idx = allowed_features[f_idx]
                else:
                    feature_idx = f_idx

                self._compute_histogram_brute_single_feature(
                    feature_idx, sample_indices, histograms
                )

        return histograms

    cdef void _compute_histogram_brute_single_feature(
            HistogramBuilder self,
            const int feature_idx,
            const unsigned int [::1] sample_indices,  # IN
            hist_struct [:, ::1] histograms) noexcept nogil:  # OUT
        """Compute the histogram for a given feature."""

        cdef:
            unsigned int n_samples = sample_indices.shape[0]
            const X_BINNED_DTYPE_C [::1] X_binned = \
                self.X_binned_dense[:, feature_idx]
            unsigned int root_node = X_binned.shape[0] == n_samples
            G_H_DTYPE_C [::1] ordered_gradients = \
                self.ordered_gradients[:n_samples]
            G_H_DTYPE_C [::1] ordered_hessians = \
                self.ordered_hessians[:n_samples]
            uint8_t hessians_are_constant = \
                self.hessians_are_constant

        # Set histograms to zero.
        memset(&histograms[feature_idx, 0], 0, self.n_bins * sizeof(hist_struct))

        if root_node:
            if hessians_are_constant:
                _build_histogram_root_no_hessian(feature_idx, X_binned,
                                                 ordered_gradients,
                                                 histograms)
            else:
                _build_histogram_root(feature_idx, X_binned,
                                      ordered_gradients, ordered_hessians,
                                      histograms)
        else:
            if hessians_are_constant:
                _build_histogram_no_hessian(feature_idx,
                                            sample_indices, X_binned,
                                            ordered_gradients, histograms)
            else:
                _build_histogram(feature_idx, sample_indices,
                                 X_binned, ordered_gradients,
                                 ordered_hessians, histograms)

    def compute_histograms_subtraction(
        HistogramBuilder self,
        hist_struct [:, ::1] parent_histograms,        # IN and OUT
        hist_struct [:, ::1] sibling_histograms,       # IN
        const unsigned int [:] allowed_features=None,  # IN
    ):
        """Compute the histograms of the node using the subtraction trick.

        hist(parent) = hist(left_child) + hist(right_child)

        For a given feature, the complexity is O(n_bins). This is much more
        efficient than compute_histograms_brute, but it's only possible for one
        of the siblings.

        Parameters
        ----------
        parent_histograms : ndarray of HISTOGRAM_DTYPE, \
                shape (n_features, n_bins)
            The histograms of the parent.
        sibling_histograms : ndarray of HISTOGRAM_DTYPE, \
                shape (n_features, n_bins)
            The histograms of the sibling.
        allowed_features : None or ndarray, dtype=np.uint32
            Indices of the features that are allowed by interaction constraints to be
            split.

        Returns
        -------
        histograms : ndarray of HISTOGRAM_DTYPE, shape(n_features, n_bins)
            The computed histograms of the current node.
            We repurpose parent_histograms for this and don't need to allocate new
            memory.
        """

        cdef:
            int feature_idx
            int f_idx
            int n_allowed_features = self.n_features
            bint has_interaction_cst = allowed_features is not None
            int n_threads = self.n_threads

        if has_interaction_cst:
            n_allowed_features = allowed_features.shape[0]

        # Compute histogram of each feature
        for f_idx in prange(n_allowed_features, schedule='static', nogil=True,
                            num_threads=n_threads):
            if has_interaction_cst:
                feature_idx = allowed_features[f_idx]
            else:
                feature_idx = f_idx

            _subtract_histograms(
                feature_idx,
                self.n_bins,
                parent_histograms,
                sibling_histograms,
            )
        return parent_histograms


cpdef void _build_histogram_naive(
        const int feature_idx,
        unsigned int [:] sample_indices,  # IN
        X_BINNED_DTYPE_C [:] binned_feature,  # IN
        G_H_DTYPE_C [:] ordered_gradients,  # IN
        G_H_DTYPE_C [:] ordered_hessians,  # IN
        hist_struct [:, :] out) noexcept nogil:  # OUT
    """Build histogram in a naive way, without optimizing for cache hit.

    Used in tests to compare with the optimized version."""
    cdef:
        unsigned int i
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int sample_idx
        unsigned int bin_idx

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        bin_idx = binned_feature[sample_idx]
        out[feature_idx, bin_idx].sum_gradients += ordered_gradients[i]
        out[feature_idx, bin_idx].sum_hessians += ordered_hessians[i]
        out[feature_idx, bin_idx].count += 1


cpdef void _subtract_histograms(
        const int feature_idx,
        unsigned int n_bins,
        hist_struct [:, ::1] hist_a,  # IN and OUT
        hist_struct [:, ::1] hist_b,  # IN
) noexcept nogil:  # OUT
    """compute hist_a = hist_a - hist_b"""
    # Note that subtraction of large sums of floating point numbers, as we have here,
    # can exhibit catastrophic cancallation. This is in particular true for gradients
    # as they can be positive and negative, while hessians are non-negative.
    # Remember that gradients and hessians are originally computed in
    # G_H_DTYPE_C = float32 precision. Therefore, if sum_gradients and sum_hessians are
    # float64, we don't loose precision. But if we also used float32 for summation, we
    # would need to take care of floating point errors.
    #
    # Note that we could protect for negative hessians by setting:
    #     sum_hessians = max(0, sum_hessians)
    # But as we use float64 for summing float32, that's veeeery unlikely.
    cdef:
        unsigned int i = 0
    for i in range(n_bins):
        hist_a[feature_idx, i].sum_gradients -= hist_b[feature_idx, i].sum_gradients
        hist_a[feature_idx, i].sum_hessians -= hist_b[feature_idx, i].sum_hessians
        hist_a[feature_idx, i].count -= hist_b[feature_idx, i].count


cpdef void _build_histogram(
        const int feature_idx,
        const unsigned int [::1] sample_indices,  # IN
        const X_BINNED_DTYPE_C [::1] binned_feature,  # IN
        const G_H_DTYPE_C [::1] ordered_gradients,  # IN
        const G_H_DTYPE_C [::1] ordered_hessians,  # IN
        hist_struct [:, ::1] out) noexcept nogil:  # OUT
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

        out[feature_idx, bin_0].sum_gradients += ordered_gradients[i]
        out[feature_idx, bin_1].sum_gradients += ordered_gradients[i + 1]
        out[feature_idx, bin_2].sum_gradients += ordered_gradients[i + 2]
        out[feature_idx, bin_3].sum_gradients += ordered_gradients[i + 3]

        out[feature_idx, bin_0].sum_hessians += ordered_hessians[i]
        out[feature_idx, bin_1].sum_hessians += ordered_hessians[i + 1]
        out[feature_idx, bin_2].sum_hessians += ordered_hessians[i + 2]
        out[feature_idx, bin_3].sum_hessians += ordered_hessians[i + 3]

        out[feature_idx, bin_0].count += 1
        out[feature_idx, bin_1].count += 1
        out[feature_idx, bin_2].count += 1
        out[feature_idx, bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        out[feature_idx, bin_idx].sum_gradients += ordered_gradients[i]
        out[feature_idx, bin_idx].sum_hessians += ordered_hessians[i]
        out[feature_idx, bin_idx].count += 1


cpdef void _build_histogram_no_hessian(
        const int feature_idx,
        const unsigned int [::1] sample_indices,  # IN
        const X_BINNED_DTYPE_C [::1] binned_feature,  # IN
        const G_H_DTYPE_C [::1] ordered_gradients,  # IN
        hist_struct [:, ::1] out) noexcept nogil:  # OUT
    """Return histogram for a given feature, not updating hessians.

    Used when the hessians of the loss are constant (typically LS loss).
    """

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

        out[feature_idx, bin_0].sum_gradients += ordered_gradients[i]
        out[feature_idx, bin_1].sum_gradients += ordered_gradients[i + 1]
        out[feature_idx, bin_2].sum_gradients += ordered_gradients[i + 2]
        out[feature_idx, bin_3].sum_gradients += ordered_gradients[i + 3]

        out[feature_idx, bin_0].count += 1
        out[feature_idx, bin_1].count += 1
        out[feature_idx, bin_2].count += 1
        out[feature_idx, bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        out[feature_idx, bin_idx].sum_gradients += ordered_gradients[i]
        out[feature_idx, bin_idx].count += 1


cpdef void _build_histogram_root(
        const int feature_idx,
        const X_BINNED_DTYPE_C [::1] binned_feature,  # IN
        const G_H_DTYPE_C [::1] all_gradients,  # IN
        const G_H_DTYPE_C [::1] all_hessians,  # IN
        hist_struct [:, ::1] out) noexcept nogil:  # OUT
    """Compute histogram of the root node.

    Unlike other nodes, the root node has to find the split among *all* the
    samples from the training set. binned_feature and all_gradients /
    all_hessians already have a consistent ordering.
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

        out[feature_idx, bin_0].sum_gradients += all_gradients[i]
        out[feature_idx, bin_1].sum_gradients += all_gradients[i + 1]
        out[feature_idx, bin_2].sum_gradients += all_gradients[i + 2]
        out[feature_idx, bin_3].sum_gradients += all_gradients[i + 3]

        out[feature_idx, bin_0].sum_hessians += all_hessians[i]
        out[feature_idx, bin_1].sum_hessians += all_hessians[i + 1]
        out[feature_idx, bin_2].sum_hessians += all_hessians[i + 2]
        out[feature_idx, bin_3].sum_hessians += all_hessians[i + 3]

        out[feature_idx, bin_0].count += 1
        out[feature_idx, bin_1].count += 1
        out[feature_idx, bin_2].count += 1
        out[feature_idx, bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        out[feature_idx, bin_idx].sum_gradients += all_gradients[i]
        out[feature_idx, bin_idx].sum_hessians += all_hessians[i]
        out[feature_idx, bin_idx].count += 1


cpdef void _build_histogram_root_no_hessian(
        const int feature_idx,
        const X_BINNED_DTYPE_C [::1] binned_feature,  # IN
        const G_H_DTYPE_C [::1] all_gradients,  # IN
        hist_struct [:, ::1] out) noexcept nogil:  # OUT
    """Compute histogram of the root node, not updating hessians.

    Used when the hessians of the loss are constant (typically LS loss).
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

        out[feature_idx, bin_0].sum_gradients += all_gradients[i]
        out[feature_idx, bin_1].sum_gradients += all_gradients[i + 1]
        out[feature_idx, bin_2].sum_gradients += all_gradients[i + 2]
        out[feature_idx, bin_3].sum_gradients += all_gradients[i + 3]

        out[feature_idx, bin_0].count += 1
        out[feature_idx, bin_1].count += 1
        out[feature_idx, bin_2].count += 1
        out[feature_idx, bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        out[feature_idx, bin_idx].sum_gradients += all_gradients[i]
        out[feature_idx, bin_idx].count += 1
