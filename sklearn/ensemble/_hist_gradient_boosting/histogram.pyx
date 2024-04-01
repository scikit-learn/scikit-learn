"""This module contains routines for building histograms."""

# Author: Nicolas Hug

cimport cython
from cython.parallel import prange
from libc.string cimport memset

import numpy as np

from ...utils._typedefs cimport uint8_t, uint16_t, uint32_t
from .common cimport BinnedData
from .common cimport hist_struct
from .common cimport Histograms
from .common cimport X_BINNED_DTYPE_FUSED_C
from .common cimport G_H_DTYPE_C


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
    X_binned : BinnedData of shape (n_samples, n_features)
        The binned input samples. Underlying arrays are Fortran-aligned.
    n_bins : int or ndarray of shape (n_features,), dtype=np.uint32
        The total number of bins for each feature, always including the bin for missing
        values as the last bin. Used to define the shape of the histograms.
        If an integer is passed, it is considered as `np.array([n_bins] * n_features)`.
    gradients : ndarray, shape (n_samples,)
        The gradients of each training sample. Those are the gradients of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    hessians : ndarray, shape (n_samples,)
        The hessians of each training sample. Those are the hessians of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    hessians_are_constant : bool
        Whether hessians are constant.

    Attributes
    ----------
    bin_offsets : ndarray of shape (n_features + 1), dtype=np.uint32
        The bin offsets specify which partition of the histograms ndarray belongs
        to which features: feature j goes from `histograms[bin_offsets[j]]` until
        `histograms[bin_offsets[j + 1] - 1]`. `bin_offsets[n_features + 1]` gives
        the total number of bins over all features.
    """
    cdef public:
        BinnedData X_binned
        unsigned int n_features
        uint32_t [::1] n_bins
        uint32_t [::1] bin_offsets
        G_H_DTYPE_C [::1] gradients
        G_H_DTYPE_C [::1] hessians
        G_H_DTYPE_C [::1] ordered_gradients
        G_H_DTYPE_C [::1] ordered_hessians
        unsigned char hessians_are_constant
        int n_threads

    def __init__(
        self,
        BinnedData X_binned,
        object n_bins,
        G_H_DTYPE_C [::1] gradients,
        G_H_DTYPE_C [::1] hessians,
        unsigned char hessians_are_constant,
        int n_threads,
    ):
        self.X_binned = X_binned
        self.n_features = X_binned.shape[1]
        if isinstance(n_bins, int):
            # Note: all histograms will have <n_bins> bins, but some of the
            # bins may be unused if a feature has a small number of unique values.
            self.n_bins = np.full(
                shape=self.n_features, fill_value=n_bins, dtype=np.uint32,
            )
        else:
            self.n_bins = n_bins
        self.gradients = gradients
        self.hessians = hessians
        # for root node, gradients and hessians are already ordered
        self.ordered_gradients = gradients.copy()
        self.ordered_hessians = hessians.copy()
        self.hessians_are_constant = hessians_are_constant
        self.n_threads = n_threads
        # bin_offsets[j] is the start of the bins of feature j,
        # bin_offsets[n_features + 1] gives the total number of bins
        bin_offsets = np.zeros(shape=self.n_features + 1, dtype=np.uint32)
        bin_offsets[1:] = np.cumsum(self.n_bins)
        self.bin_offsets = bin_offsets

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
            # need local views to avoid python interactions
            unsigned char hessians_are_constant = self.hessians_are_constant
            int n_allowed_features = self.n_features
            G_H_DTYPE_C [::1] ordered_gradients = self.ordered_gradients
            G_H_DTYPE_C [::1] gradients = self.gradients
            G_H_DTYPE_C [::1] ordered_hessians = self.ordered_hessians
            G_H_DTYPE_C [::1] hessians = self.hessians
            # Histograms will be initialized to zero later within a prange.
            # Here we just allocate the array, i.e. __init__ calls:
            # np.empty(shape=self.bin_offsets[-1], dtype=HISTOGRAM_DTYPE)
            Histograms histograms = Histograms(
                n_features=self.n_features, bin_offsets=self.bin_offsets
            )
            bint has_interaction_cst = allowed_features is not None
            int n_threads = self.n_threads

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

                if self.X_binned.feature_is_8bit_view[feature_idx]:
                    self._compute_histogram_brute_single_feature[uint8_t](
                        X_binned=self.X_binned.get_feature_view8(feature_idx),
                        feature_idx=feature_idx,
                        sample_indices=sample_indices,
                        histograms=histograms,
                    )
                else:
                    self._compute_histogram_brute_single_feature[uint16_t](
                        X_binned=self.X_binned.get_feature_view16(feature_idx),
                        feature_idx=feature_idx,
                        sample_indices=sample_indices,
                        histograms=histograms,
                    )

        return histograms

    cdef void _compute_histogram_brute_single_feature(
        HistogramBuilder self,
        const X_BINNED_DTYPE_FUSED_C [::1] X_binned,
        const int feature_idx,
        const unsigned int [::1] sample_indices,  # IN
        Histograms histograms,                    # OUT
    ) noexcept nogil:
        """Compute the histogram for a given feature."""

        cdef:
            unsigned int n_samples = sample_indices.shape[0]
            unsigned int root_node = X_binned.shape[0] == n_samples
            G_H_DTYPE_C [::1] ordered_gradients = \
                self.ordered_gradients[:n_samples]
            G_H_DTYPE_C [::1] ordered_hessians = \
                self.ordered_hessians[:n_samples]
            unsigned char hessians_are_constant = \
                self.hessians_are_constant
            uint32_t n_bins = histograms.n_bins(feature_idx)
            hist_struct * hist = histograms.at(feature_idx, 0)

        # Set histograms to zero.
        memset(hist, 0, n_bins * sizeof(hist_struct))

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
        Histograms parent_histograms,        # IN and OUT
        Histograms sibling_histograms,       # IN
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
                parent_histograms,
                sibling_histograms,
            )
        return parent_histograms


cpdef void _build_histogram_naive(
    const int feature_idx,
    unsigned int [:] sample_indices,      # IN
    X_BINNED_DTYPE_FUSED_C [:] binned_feature,  # IN
    G_H_DTYPE_C [:] ordered_gradients,    # IN
    G_H_DTYPE_C [:] ordered_hessians,     # IN
    Histograms out,                       # OUT
) noexcept nogil:
    """Build histogram in a naive way, without optimizing for cache hit.

    Used in tests to compare with the optimized version."""
    cdef:
        unsigned int i
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int sample_idx
        uint32_t bin_idx

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        bin_idx = binned_feature[sample_idx]
        out.at(feature_idx, bin_idx).sum_gradients += ordered_gradients[i]
        out.at(feature_idx, bin_idx).sum_hessians += ordered_hessians[i]
        out.at(feature_idx, bin_idx).count += 1


cpdef void _subtract_histograms(
    const int feature_idx,
    Histograms hist_a,  # IN and OUT
    Histograms hist_b,  # IN
) noexcept nogil:
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
        uint32_t i = 0
        uint32_t n_bins = hist_a.n_bins(feature_idx)
        hist_struct * ha = hist_a.at(feature_idx, 0)
        hist_struct * hb = hist_b.at(feature_idx, 0)
    for i in range(n_bins):
        ha[i].sum_gradients -= hb[i].sum_gradients  # no-cython-lint
        ha[i].sum_hessians  -= hb[i].sum_hessians   # no-cython-lint
        ha[i].count         -= hb[i].count          # no-cython-lint


cpdef void _build_histogram(
    const int feature_idx,
    const unsigned int [::1] sample_indices,      # IN
    const X_BINNED_DTYPE_FUSED_C [::1] binned_feature,  # IN
    const G_H_DTYPE_C [::1] ordered_gradients,    # IN
    const G_H_DTYPE_C [::1] ordered_hessians,     # IN
    Histograms out,                               # OUT
) noexcept nogil:
    """Return histogram for a given feature."""
    cdef:
        unsigned int i = 0
        unsigned int n_node_samples = sample_indices.shape[0]
        unsigned int unrolled_upper = (n_node_samples // 4) * 4

        uint32_t bin_0
        uint32_t bin_1
        uint32_t bin_2
        uint32_t bin_3
        uint32_t bin_idx
        hist_struct * hist = out.at(feature_idx, 0)

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        hist[bin_0].sum_gradients += ordered_gradients[i]
        hist[bin_1].sum_gradients += ordered_gradients[i + 1]
        hist[bin_2].sum_gradients += ordered_gradients[i + 2]
        hist[bin_3].sum_gradients += ordered_gradients[i + 3]

        hist[bin_0].sum_hessians += ordered_hessians[i]
        hist[bin_1].sum_hessians += ordered_hessians[i + 1]
        hist[bin_2].sum_hessians += ordered_hessians[i + 2]
        hist[bin_3].sum_hessians += ordered_hessians[i + 3]

        hist[bin_0].count += 1
        hist[bin_1].count += 1
        hist[bin_2].count += 1
        hist[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        hist[bin_idx].sum_gradients += ordered_gradients[i]
        hist[bin_idx].sum_hessians += ordered_hessians[i]
        hist[bin_idx].count += 1


cpdef void _build_histogram_no_hessian(
    const int feature_idx,
    const unsigned int [::1] sample_indices,      # IN
    const X_BINNED_DTYPE_FUSED_C [::1] binned_feature,  # IN
    const G_H_DTYPE_C [::1] ordered_gradients,    # IN
    Histograms out,                     # OUT
) noexcept nogil:
    """Return histogram for a given feature, not updating hessians.

    Used when the hessians of the loss are constant (typically LS loss).
    """
    cdef:
        unsigned int i = 0
        unsigned int n_node_samples = sample_indices.shape[0]
        unsigned int unrolled_upper = (n_node_samples // 4) * 4

        uint32_t bin_0
        uint32_t bin_1
        uint32_t bin_2
        uint32_t bin_3
        uint32_t bin_idx
        hist_struct * hist = out.at(feature_idx, 0)

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[sample_indices[i]]
        bin_1 = binned_feature[sample_indices[i + 1]]
        bin_2 = binned_feature[sample_indices[i + 2]]
        bin_3 = binned_feature[sample_indices[i + 3]]

        hist[bin_0].sum_gradients += ordered_gradients[i]
        hist[bin_1].sum_gradients += ordered_gradients[i + 1]
        hist[bin_2].sum_gradients += ordered_gradients[i + 2]
        hist[bin_3].sum_gradients += ordered_gradients[i + 3]

        hist[bin_0].count += 1
        hist[bin_1].count += 1
        hist[bin_2].count += 1
        hist[bin_3].count += 1

    for i in range(unrolled_upper, n_node_samples):
        bin_idx = binned_feature[sample_indices[i]]
        hist[bin_idx].sum_gradients += ordered_gradients[i]
        hist[bin_idx].count += 1


cpdef void _build_histogram_root(
    const int feature_idx,
    const X_BINNED_DTYPE_FUSED_C [::1] binned_feature,  # IN
    const G_H_DTYPE_C [::1] all_gradients,        # IN
    const G_H_DTYPE_C [::1] all_hessians,         # IN
    Histograms out,                     # OUT
) noexcept nogil:
    """Compute histogram of the root node.

    Unlike other nodes, the root node has to find the split among *all* the
    samples from the training set. binned_feature and all_gradients /
    all_hessians already have a consistent ordering.
    """

    cdef:
        unsigned int i = 0
        unsigned int n_samples = binned_feature.shape[0]
        unsigned int unrolled_upper = (n_samples // 4) * 4

        uint32_t bin_0
        uint32_t bin_1
        uint32_t bin_2
        uint32_t bin_3
        uint32_t bin_idx
        hist_struct * hist = out.at(feature_idx, 0)

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        hist[bin_0].sum_gradients += all_gradients[i]
        hist[bin_1].sum_gradients += all_gradients[i + 1]
        hist[bin_2].sum_gradients += all_gradients[i + 2]
        hist[bin_3].sum_gradients += all_gradients[i + 3]

        hist[bin_0].sum_hessians += all_hessians[i]
        hist[bin_1].sum_hessians += all_hessians[i + 1]
        hist[bin_2].sum_hessians += all_hessians[i + 2]
        hist[bin_3].sum_hessians += all_hessians[i + 3]

        hist[bin_0].count += 1
        hist[bin_1].count += 1
        hist[bin_2].count += 1
        hist[bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        hist[bin_idx].sum_gradients += all_gradients[i]
        hist[bin_idx].sum_hessians += all_hessians[i]
        hist[bin_idx].count += 1


cpdef void _build_histogram_root_no_hessian(
    const int feature_idx,
    const X_BINNED_DTYPE_FUSED_C [::1] binned_feature,  # IN
    const G_H_DTYPE_C [::1] all_gradients,        # IN
    Histograms out,                               # OUT
) noexcept nogil:
    """Compute histogram of the root node, not updating hessians.

    Used when the hessians of the loss are constant (typically LS loss).
    """

    cdef:
        unsigned int i = 0
        unsigned int n_samples = binned_feature.shape[0]
        unsigned int unrolled_upper = (n_samples // 4) * 4

        uint32_t bin_0
        uint32_t bin_1
        uint32_t bin_2
        uint32_t bin_3
        uint32_t bin_idx
        hist_struct * hist = out.at(feature_idx, 0)

    for i in range(0, unrolled_upper, 4):
        bin_0 = binned_feature[i]
        bin_1 = binned_feature[i + 1]
        bin_2 = binned_feature[i + 2]
        bin_3 = binned_feature[i + 3]

        hist[bin_0].sum_gradients += all_gradients[i]
        hist[bin_1].sum_gradients += all_gradients[i + 1]
        hist[bin_2].sum_gradients += all_gradients[i + 2]
        hist[bin_3].sum_gradients += all_gradients[i + 3]

        hist[bin_0].count += 1
        hist[bin_1].count += 1
        hist[bin_2].count += 1
        hist[bin_3].count += 1

    for i in range(unrolled_upper, n_samples):
        bin_idx = binned_feature[i]
        hist[bin_idx].sum_gradients += all_gradients[i]
        hist[bin_idx].count += 1
