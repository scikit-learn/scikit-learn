"""This module contains routines and data structures to:

- Find the best possible split of a node. For a given node, a split is
  characterized by a feature and a bin.
- Apply a split to a node, i.e. split the indices of the samples at the node
  into the newly created left and right children.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

cimport cython
from cython.parallel import prange
import numpy as np
from libc.math cimport INFINITY, ceil
from libc.stdlib cimport malloc, free, qsort
from libc.string cimport memcpy

from ...utils._typedefs cimport uint8_t
from .common cimport X_BINNED_DTYPE_C
from .common cimport Y_DTYPE_C
from .common cimport hist_struct
from .common cimport BITSET_INNER_DTYPE_C
from .common cimport BITSET_DTYPE_C
from .common cimport MonotonicConstraint
from ._bitset cimport init_bitset
from ._bitset cimport set_bitset
from ._bitset cimport in_bitset


cdef struct split_info_struct:
    # Same as the SplitInfo class, but we need a C struct to use it in the
    # nogil sections and to use in arrays.
    Y_DTYPE_C gain
    int feature_idx
    unsigned int bin_idx
    uint8_t missing_go_to_left
    Y_DTYPE_C sum_gradient_left
    Y_DTYPE_C sum_gradient_right
    Y_DTYPE_C sum_hessian_left
    Y_DTYPE_C sum_hessian_right
    unsigned int n_samples_left
    unsigned int n_samples_right
    Y_DTYPE_C value_left
    Y_DTYPE_C value_right
    uint8_t is_categorical
    BITSET_DTYPE_C left_cat_bitset


# used in categorical splits for sorting categories by increasing values of
# sum_gradients / sum_hessians
cdef struct categorical_info:
    X_BINNED_DTYPE_C bin_idx
    Y_DTYPE_C value


class SplitInfo:
    """Pure data class to store information about a potential split.

    Parameters
    ----------
    gain : float
        The gain of the split.
    feature_idx : int
        The index of the feature to be split.
    bin_idx : int
        The index of the bin on which the split is made. Should be ignored if
        `is_categorical` is True: `left_cat_bitset` will be used to determine
        the split.
    missing_go_to_left : bool
        Whether missing values should go to the left child. This is used
        whether the split is categorical or not.
    sum_gradient_left : float
        The sum of the gradients of all the samples in the left child.
    sum_hessian_left : float
        The sum of the hessians of all the samples in the left child.
    sum_gradient_right : float
        The sum of the gradients of all the samples in the right child.
    sum_hessian_right : float
        The sum of the hessians of all the samples in the right child.
    n_samples_left : int, default=0
        The number of samples in the left child.
    n_samples_right : int
        The number of samples in the right child.
    is_categorical : bool
        Whether the split is done on a categorical feature.
    left_cat_bitset : ndarray of shape=(8,), dtype=uint32 or None
        Bitset representing the categories that go to the left. This is used
        only when `is_categorical` is True.
        Note that missing values are part of that bitset if there are missing
        values in the training data. For missing values, we rely on that
        bitset for splitting, but at prediction time, we rely on
        missing_go_to_left.
    """
    def __init__(self, gain, feature_idx, bin_idx,
                 missing_go_to_left, sum_gradient_left, sum_hessian_left,
                 sum_gradient_right, sum_hessian_right, n_samples_left,
                 n_samples_right, value_left, value_right,
                 is_categorical, left_cat_bitset):
        self.gain = gain
        self.feature_idx = feature_idx
        self.bin_idx = bin_idx
        self.missing_go_to_left = missing_go_to_left
        self.sum_gradient_left = sum_gradient_left
        self.sum_hessian_left = sum_hessian_left
        self.sum_gradient_right = sum_gradient_right
        self.sum_hessian_right = sum_hessian_right
        self.n_samples_left = n_samples_left
        self.n_samples_right = n_samples_right
        self.value_left = value_left
        self.value_right = value_right
        self.is_categorical = is_categorical
        self.left_cat_bitset = left_cat_bitset


@cython.final
cdef class Splitter:
    """Splitter used to find the best possible split at each node.

    A split (see SplitInfo) is characterized by a feature and a bin.

    The Splitter is also responsible for partitioning the samples among the
    leaves of the tree (see split_indices() and the partition attribute).

    Parameters
    ----------
    X_binned : ndarray of int, shape (n_samples, n_features)
        The binned input samples. Must be Fortran-aligned.
    n_bins_non_missing : ndarray, shape (n_features,)
        For each feature, gives the number of bins actually used for
        non-missing values.
    missing_values_bin_idx : uint8
        Index of the bin that is used for missing values. This is the index of
        the last bin and is always equal to max_bins (as passed to the GBDT
        classes), or equivalently to n_bins - 1.
    has_missing_values : ndarray, shape (n_features,)
        Whether missing values were observed in the training data, for each
        feature.
    is_categorical : ndarray of bool of shape (n_features,)
        Indicates categorical features.
    monotonic_cst : ndarray of int of shape (n_features,), dtype=int
        Indicates the monotonic constraint to enforce on each feature.
          - 1: monotonic increase
          - 0: no constraint
          - -1: monotonic decrease

        Read more in the :ref:`User Guide <monotonic_cst_gbdt>`.
    l2_regularization : float
        The L2 regularization parameter.
    min_hessian_to_split : float, default=1e-3
        The minimum sum of hessians needed in each node. Splits that result in
        at least one child having a sum of hessians less than
        min_hessian_to_split are discarded.
    min_samples_leaf : int, default=20
        The minimum number of samples per leaf.
    min_gain_to_split : float, default=0.0
        The minimum gain needed to split a node. Splits with lower gain will
        be ignored.
    hessians_are_constant: bool, default is False
        Whether hessians are constant.
    feature_fraction_per_split : float, default=1
        Proportion of randomly chosen features in each and every node split.
        This is a form of regularization, smaller values make the trees weaker
        learners and might prevent overfitting.
    rng : Generator
    n_threads : int, default=1
        Number of OpenMP threads to use.
    """
    cdef public:
        const X_BINNED_DTYPE_C [::1, :] X_binned
        unsigned int n_features
        const unsigned int [::1] n_bins_non_missing
        uint8_t missing_values_bin_idx
        const uint8_t [::1] has_missing_values
        const uint8_t [::1] is_categorical
        const signed char [::1] monotonic_cst
        uint8_t hessians_are_constant
        Y_DTYPE_C l2_regularization
        Y_DTYPE_C min_hessian_to_split
        unsigned int min_samples_leaf
        Y_DTYPE_C min_gain_to_split
        Y_DTYPE_C feature_fraction_per_split
        rng

        unsigned int [::1] partition
        unsigned int [::1] left_indices_buffer
        unsigned int [::1] right_indices_buffer
        int n_threads

    def __init__(self,
                 const X_BINNED_DTYPE_C [::1, :] X_binned,
                 const unsigned int [::1] n_bins_non_missing,
                 const uint8_t missing_values_bin_idx,
                 const uint8_t [::1] has_missing_values,
                 const uint8_t [::1] is_categorical,
                 const signed char [::1] monotonic_cst,
                 Y_DTYPE_C l2_regularization,
                 Y_DTYPE_C min_hessian_to_split=1e-3,
                 unsigned int min_samples_leaf=20,
                 Y_DTYPE_C min_gain_to_split=0.,
                 uint8_t hessians_are_constant=False,
                 Y_DTYPE_C feature_fraction_per_split=1.0,
                 rng=np.random.RandomState(),
                 unsigned int n_threads=1):

        self.X_binned = X_binned
        self.n_features = X_binned.shape[1]
        self.n_bins_non_missing = n_bins_non_missing
        self.missing_values_bin_idx = missing_values_bin_idx
        self.has_missing_values = has_missing_values
        self.is_categorical = is_categorical
        self.monotonic_cst = monotonic_cst
        self.l2_regularization = l2_regularization
        self.min_hessian_to_split = min_hessian_to_split
        self.min_samples_leaf = min_samples_leaf
        self.min_gain_to_split = min_gain_to_split
        self.hessians_are_constant = hessians_are_constant
        self.feature_fraction_per_split = feature_fraction_per_split
        self.rng = rng
        self.n_threads = n_threads

        # The partition array maps each sample index into the leaves of the
        # tree (a leaf in this context is a node that isn't split yet, not
        # necessarily a 'finalized' leaf). Initially, the root contains all
        # the indices, e.g.:
        # partition = [abcdefghijkl]
        # After a call to split_indices, it may look e.g. like this:
        # partition = [cef|abdghijkl]
        # we have 2 leaves, the left one is at position 0 and the second one at
        # position 3. The order of the samples is irrelevant.
        self.partition = np.arange(X_binned.shape[0], dtype=np.uint32)
        # buffers used in split_indices to support parallel splitting.
        self.left_indices_buffer = np.empty_like(self.partition)
        self.right_indices_buffer = np.empty_like(self.partition)

    def split_indices(Splitter self, split_info, unsigned int [::1]
                      sample_indices):
        """Split samples into left and right arrays.

        The split is performed according to the best possible split
        (split_info).

        Ultimately, this is nothing but a partition of the sample_indices
        array with a given pivot, exactly like a quicksort subroutine.

        Parameters
        ----------
        split_info : SplitInfo
            The SplitInfo of the node to split.
        sample_indices : ndarray of unsigned int, shape (n_samples_at_node,)
            The indices of the samples at the node to split. This is a view
            on self.partition, and it is modified inplace by placing the
            indices of the left child at the beginning, and the indices of
            the right child at the end.

        Returns
        -------
        left_indices : ndarray of int, shape (n_left_samples,)
            The indices of the samples in the left child. This is a view on
            self.partition.
        right_indices : ndarray of int, shape (n_right_samples,)
            The indices of the samples in the right child. This is a view on
            self.partition.
        right_child_position : int
            The position of the right child in ``sample_indices``.
        """
        # This is a multi-threaded implementation inspired by lightgbm. Here
        # is a quick break down. Let's suppose we want to split a node with 24
        # samples named from a to x. self.partition looks like this (the * are
        # indices in other leaves that we don't care about):
        # partition = [*************abcdefghijklmnopqrstuvwx****************]
        #                           ^                       ^
        #                     node_position     node_position + node.n_samples

        # Ultimately, we want to reorder the samples inside the boundaries of
        # the leaf (which becomes a node) to now represent the samples in its
        # left and right child. For example:
        # partition = [*************abefilmnopqrtuxcdghjksvw*****************]
        #                           ^              ^
        #                   left_child_pos     right_child_pos
        # Note that left_child_pos always takes the value of node_position,
        # and right_child_pos = left_child_pos + left_child.n_samples. The
        # order of the samples inside a leaf is irrelevant.

        # 1. sample_indices is a view on this region a..x. We conceptually
        #    divide it into n_threads regions. Each thread will be responsible
        #    for its own region. Here is an example with 4 threads:
        #    sample_indices = [abcdef|ghijkl|mnopqr|stuvwx]
        # 2. Each thread processes 6 = 24 // 4 entries and maps them into
        #    left_indices_buffer or right_indices_buffer. For example, we could
        #    have the following mapping ('.' denotes an undefined entry):
        #    - left_indices_buffer =  [abef..|il....|mnopqr|tux...]
        #    - right_indices_buffer = [cd....|ghjk..|......|svw...]
        # 3. We keep track of the start positions of the regions (the '|') in
        #    ``offset_in_buffers`` as well as the size of each region. We also
        #    keep track of the number of samples put into the left/right child
        #    by each thread. Concretely:
        #    - left_counts =  [4, 2, 6, 3]
        #    - right_counts = [2, 4, 0, 3]
        # 4. Finally, we put left/right_indices_buffer back into the
        #    sample_indices, without any undefined entries and the partition
        #    looks as expected
        #    partition = [*************abefilmnopqrtuxcdghjksvw***************]

        # Note: We here show left/right_indices_buffer as being the same size
        # as sample_indices for simplicity, but in reality they are of the
        # same size as partition.

        cdef:
            int n_samples = sample_indices.shape[0]
            X_BINNED_DTYPE_C bin_idx = split_info.bin_idx
            uint8_t missing_go_to_left = split_info.missing_go_to_left
            uint8_t missing_values_bin_idx = self.missing_values_bin_idx
            int feature_idx = split_info.feature_idx
            const X_BINNED_DTYPE_C [::1] X_binned = \
                self.X_binned[:, feature_idx]
            unsigned int [::1] left_indices_buffer = self.left_indices_buffer
            unsigned int [::1] right_indices_buffer = self.right_indices_buffer
            uint8_t is_categorical = split_info.is_categorical
            # Cython is unhappy if we set left_cat_bitset to
            # split_info.left_cat_bitset directly, so we need a tmp var
            BITSET_INNER_DTYPE_C [:] cat_bitset_tmp = split_info.left_cat_bitset
            BITSET_DTYPE_C left_cat_bitset
            int n_threads = self.n_threads

            int [:] sizes = np.full(n_threads, n_samples // n_threads,
                                    dtype=np.int32)
            int [:] offset_in_buffers = np.zeros(n_threads, dtype=np.int32)
            int [:] left_counts = np.empty(n_threads, dtype=np.int32)
            int [:] right_counts = np.empty(n_threads, dtype=np.int32)
            int left_count
            int right_count
            int start
            int stop
            int i
            int thread_idx
            int sample_idx
            int right_child_position
            uint8_t turn_left
            int [:] left_offset = np.zeros(n_threads, dtype=np.int32)
            int [:] right_offset = np.zeros(n_threads, dtype=np.int32)

        # only set left_cat_bitset when is_categorical is True
        if is_categorical:
            left_cat_bitset = &cat_bitset_tmp[0]

        with nogil:
            for thread_idx in range(n_samples % n_threads):
                sizes[thread_idx] += 1

            for thread_idx in range(1, n_threads):
                offset_in_buffers[thread_idx] = \
                    offset_in_buffers[thread_idx - 1] + sizes[thread_idx - 1]

            # map indices from sample_indices to left/right_indices_buffer
            for thread_idx in prange(n_threads, schedule='static',
                                     chunksize=1, num_threads=n_threads):
                left_count = 0
                right_count = 0

                start = offset_in_buffers[thread_idx]
                stop = start + sizes[thread_idx]
                for i in range(start, stop):
                    sample_idx = sample_indices[i]
                    turn_left = sample_goes_left(
                        missing_go_to_left,
                        missing_values_bin_idx, bin_idx,
                        X_binned[sample_idx], is_categorical,
                        left_cat_bitset)

                    if turn_left:
                        left_indices_buffer[start + left_count] = sample_idx
                        left_count = left_count + 1
                    else:
                        right_indices_buffer[start + right_count] = sample_idx
                        right_count = right_count + 1

                left_counts[thread_idx] = left_count
                right_counts[thread_idx] = right_count

            # position of right child = just after the left child
            right_child_position = 0
            for thread_idx in range(n_threads):
                right_child_position += left_counts[thread_idx]

            # offset of each thread in sample_indices for left and right
            # child, i.e. where each thread will start to write.
            right_offset[0] = right_child_position
            for thread_idx in range(1, n_threads):
                left_offset[thread_idx] = \
                    left_offset[thread_idx - 1] + left_counts[thread_idx - 1]
                right_offset[thread_idx] = \
                    right_offset[thread_idx - 1] + right_counts[thread_idx - 1]

            # map indices in left/right_indices_buffer back into
            # sample_indices. This also updates self.partition since
            # sample_indices is a view.
            for thread_idx in prange(n_threads, schedule='static',
                                     chunksize=1, num_threads=n_threads):
                memcpy(
                    &sample_indices[left_offset[thread_idx]],
                    &left_indices_buffer[offset_in_buffers[thread_idx]],
                    sizeof(unsigned int) * left_counts[thread_idx]
                )
                if right_counts[thread_idx] > 0:
                    # If we're splitting the rightmost node of the tree, i.e. the
                    # rightmost node in the partition array, and if n_threads >= 2, one
                    # might have right_counts[-1] = 0 and right_offset[-1] = len(sample_indices)
                    # leading to evaluating
                    #
                    #    &sample_indices[right_offset[-1]] = &samples_indices[n_samples_at_node]
                    #                                      = &partition[n_samples_in_tree]
                    #
                    # which is an out-of-bounds read access that can cause a segmentation fault.
                    # When boundscheck=True, removing this check produces this exception:
                    #
                    #    IndexError: Out of bounds on buffer access
                    #
                    memcpy(
                        &sample_indices[right_offset[thread_idx]],
                        &right_indices_buffer[offset_in_buffers[thread_idx]],
                        sizeof(unsigned int) * right_counts[thread_idx]
                    )

        return (sample_indices[:right_child_position],
                sample_indices[right_child_position:],
                right_child_position)

    def find_node_split(
            Splitter self,
            unsigned int n_samples,
            hist_struct [:, ::1] histograms,  # IN
            const Y_DTYPE_C sum_gradients,
            const Y_DTYPE_C sum_hessians,
            const Y_DTYPE_C value,
            const Y_DTYPE_C lower_bound=-INFINITY,
            const Y_DTYPE_C upper_bound=INFINITY,
            const unsigned int [:] allowed_features=None,
            ):
        """For each feature, find the best bin to split on at a given node.

        Return the best split info among all features.

        Parameters
        ----------
        n_samples : int
            The number of samples at the node.
        histograms : ndarray of HISTOGRAM_DTYPE of \
                shape (n_features, max_bins)
            The histograms of the current node.
        sum_gradients : float
            The sum of the gradients for each sample at the node.
        sum_hessians : float
            The sum of the hessians for each sample at the node.
        value : float
            The bounded value of the current node. We directly pass the value
            instead of re-computing it from sum_gradients and sum_hessians,
            because we need to compute the loss and the gain based on the
            *bounded* value: computing the value from
            sum_gradients / sum_hessians would give the unbounded value, and
            the interaction with min_gain_to_split would not be correct
            anymore. Side note: we can't use the lower_bound / upper_bound
            parameters either because these refer to the bounds of the
            children, not the bounds of the current node.
        lower_bound : float
            Lower bound for the children values for respecting the monotonic
            constraints.
        upper_bound : float
            Upper bound for the children values for respecting the monotonic
            constraints.
        allowed_features : None or ndarray, dtype=np.uint32
            Indices of the features that are allowed by interaction constraints to be
            split.

        Returns
        -------
        best_split_info : SplitInfo
            The info about the best possible split among all features.
        """
        cdef:
            int feature_idx
            int split_info_idx
            int best_split_info_idx
            int n_allowed_features
            split_info_struct split_info
            split_info_struct * split_infos
            const uint8_t [::1] has_missing_values = self.has_missing_values
            const uint8_t [::1] is_categorical = self.is_categorical
            const signed char [::1] monotonic_cst = self.monotonic_cst
            int n_threads = self.n_threads
            bint has_interaction_cst = False
            Y_DTYPE_C feature_fraction_per_split = self.feature_fraction_per_split
            uint8_t [:] subsample_mask  # same as npy_bool
            int n_subsampled_features

        has_interaction_cst = allowed_features is not None
        if has_interaction_cst:
            n_allowed_features = allowed_features.shape[0]
        else:
            n_allowed_features = self.n_features

        if feature_fraction_per_split < 1.0:
            # We do all random sampling before the nogil and make sure that we sample
            # exactly n_subsampled_features >= 1 features.
            n_subsampled_features = max(
                1,
                int(ceil(feature_fraction_per_split * n_allowed_features)),
            )
            subsample_mask_arr = np.full(n_allowed_features, False)
            subsample_mask_arr[:n_subsampled_features] = True
            self.rng.shuffle(subsample_mask_arr)
            # https://github.com/numpy/numpy/issues/18273
            subsample_mask = subsample_mask_arr

        with nogil:

            split_infos = <split_info_struct *> malloc(
                n_allowed_features * sizeof(split_info_struct))

            # split_info_idx is index of split_infos of size n_allowed_features.
            # features_idx is the index of the feature column in X.
            for split_info_idx in prange(n_allowed_features, schedule='static',
                                         num_threads=n_threads):
                if has_interaction_cst:
                    feature_idx = allowed_features[split_info_idx]
                else:
                    feature_idx = split_info_idx

                split_infos[split_info_idx].feature_idx = feature_idx

                # For each feature, find best bin to split on
                # Start with a gain of -1 if no better split is found, that
                # means one of the constraints isn't respected
                # (min_samples_leaf, etc.) and the grower will later turn the
                # node into a leaf.
                split_infos[split_info_idx].gain = -1
                split_infos[split_info_idx].is_categorical = is_categorical[feature_idx]

                # Note that subsample_mask is indexed by split_info_idx and not by
                # feature_idx because we only need to exclude the same features again
                # and again. We do NOT need to access the features directly by using
                # allowed_features.
                if feature_fraction_per_split < 1.0 and not subsample_mask[split_info_idx]:
                    continue

                if is_categorical[feature_idx]:
                    self._find_best_bin_to_split_category(
                        feature_idx, has_missing_values[feature_idx],
                        histograms, n_samples, sum_gradients, sum_hessians,
                        value, monotonic_cst[feature_idx], lower_bound,
                        upper_bound, &split_infos[split_info_idx])
                else:
                    # We will scan bins from left to right (in all cases), and
                    # if there are any missing values, we will also scan bins
                    # from right to left. This way, we can consider whichever
                    # case yields the best gain: either missing values go to
                    # the right (left to right scan) or to the left (right to
                    # left case). See algo 3 from the XGBoost paper
                    # https://arxiv.org/abs/1603.02754
                    # Note: for the categorical features above, this isn't
                    # needed since missing values are considered a native
                    # category.
                    self._find_best_bin_to_split_left_to_right(
                        feature_idx, has_missing_values[feature_idx],
                        histograms, n_samples, sum_gradients, sum_hessians,
                        value, monotonic_cst[feature_idx],
                        lower_bound, upper_bound, &split_infos[split_info_idx])

                    if has_missing_values[feature_idx]:
                        # We need to explore both directions to check whether
                        # sending the nans to the left child would lead to a higher
                        # gain
                        self._find_best_bin_to_split_right_to_left(
                            feature_idx, histograms, n_samples,
                            sum_gradients, sum_hessians,
                            value, monotonic_cst[feature_idx],
                            lower_bound, upper_bound, &split_infos[split_info_idx])

            # then compute best possible split among all features
            # split_info is set to the best of split_infos
            best_split_info_idx = self._find_best_feature_to_split_helper(
                split_infos, n_allowed_features
            )
            split_info = split_infos[best_split_info_idx]

        out = SplitInfo(
            split_info.gain,
            split_info.feature_idx,
            split_info.bin_idx,
            split_info.missing_go_to_left,
            split_info.sum_gradient_left,
            split_info.sum_hessian_left,
            split_info.sum_gradient_right,
            split_info.sum_hessian_right,
            split_info.n_samples_left,
            split_info.n_samples_right,
            split_info.value_left,
            split_info.value_right,
            split_info.is_categorical,
            None,  # left_cat_bitset will only be set if the split is categorical
        )
        # Only set bitset if the split is categorical
        if split_info.is_categorical:
            out.left_cat_bitset = np.asarray(split_info.left_cat_bitset, dtype=np.uint32)

        free(split_infos)
        return out

    cdef int _find_best_feature_to_split_helper(
        self,
        split_info_struct * split_infos,  # IN
        int n_allowed_features,
    ) noexcept nogil:
        """Return the index of split_infos with the best feature split."""
        cdef:
            int split_info_idx
            int best_split_info_idx = 0

        for split_info_idx in range(1, n_allowed_features):
            if (split_infos[split_info_idx].gain > split_infos[best_split_info_idx].gain):
                best_split_info_idx = split_info_idx
        return best_split_info_idx

    cdef void _find_best_bin_to_split_left_to_right(
            Splitter self,
            unsigned int feature_idx,
            uint8_t has_missing_values,
            const hist_struct [:, ::1] histograms,  # IN
            unsigned int n_samples,
            Y_DTYPE_C sum_gradients,
            Y_DTYPE_C sum_hessians,
            Y_DTYPE_C value,
            signed char monotonic_cst,
            Y_DTYPE_C lower_bound,
            Y_DTYPE_C upper_bound,
            split_info_struct * split_info) noexcept nogil:  # OUT
        """Find best bin to split on for a given feature.

        Splits that do not satisfy the splitting constraints
        (min_gain_to_split, etc.) are discarded here.

        We scan node from left to right. This version is called whether there
        are missing values or not. If any, missing values are assigned to the
        right node.
        """
        cdef:
            unsigned int bin_idx
            unsigned int n_samples_left
            unsigned int n_samples_right
            unsigned int n_samples_ = n_samples
            # We set the 'end' variable such that the last non-missing-values
            # bin never goes to the left child (which would result in and
            # empty right child), unless there are missing values, since these
            # would go to the right child.
            unsigned int end = \
                self.n_bins_non_missing[feature_idx] - 1 + has_missing_values
            Y_DTYPE_C sum_hessian_left
            Y_DTYPE_C sum_hessian_right
            Y_DTYPE_C sum_gradient_left
            Y_DTYPE_C sum_gradient_right
            Y_DTYPE_C loss_current_node
            Y_DTYPE_C gain
            uint8_t found_better_split = False

            Y_DTYPE_C best_sum_hessian_left
            Y_DTYPE_C best_sum_gradient_left
            unsigned int best_bin_idx
            unsigned int best_n_samples_left
            Y_DTYPE_C best_gain = -1
            hist_struct hist

        sum_gradient_left, sum_hessian_left = 0., 0.
        n_samples_left = 0

        loss_current_node = _loss_from_value(value, sum_gradients)

        for bin_idx in range(end):
            hist = histograms[feature_idx, bin_idx]
            n_samples_left += hist.count
            n_samples_right = n_samples_ - n_samples_left

            if self.hessians_are_constant:
                sum_hessian_left += hist.count
            else:
                sum_hessian_left += \
                    hist.sum_hessians
            sum_hessian_right = sum_hessians - sum_hessian_left

            sum_gradient_left += hist.sum_gradients
            sum_gradient_right = sum_gradients - sum_gradient_left

            if n_samples_left < self.min_samples_leaf:
                continue
            if n_samples_right < self.min_samples_leaf:
                # won't get any better
                break

            if sum_hessian_left < self.min_hessian_to_split:
                continue
            if sum_hessian_right < self.min_hessian_to_split:
                # won't get any better (hessians are > 0 since loss is convex)
                break

            gain = _split_gain(sum_gradient_left, sum_hessian_left,
                               sum_gradient_right, sum_hessian_right,
                               loss_current_node,
                               monotonic_cst,
                               lower_bound,
                               upper_bound,
                               self.l2_regularization)

            if gain > best_gain and gain > self.min_gain_to_split:
                found_better_split = True
                best_gain = gain
                best_bin_idx = bin_idx
                best_sum_gradient_left = sum_gradient_left
                best_sum_hessian_left = sum_hessian_left
                best_n_samples_left = n_samples_left

        if found_better_split:
            split_info.gain = best_gain
            split_info.bin_idx = best_bin_idx
            # we scan from left to right so missing values go to the right
            split_info.missing_go_to_left = False
            split_info.sum_gradient_left = best_sum_gradient_left
            split_info.sum_gradient_right = sum_gradients - best_sum_gradient_left
            split_info.sum_hessian_left = best_sum_hessian_left
            split_info.sum_hessian_right = sum_hessians - best_sum_hessian_left
            split_info.n_samples_left = best_n_samples_left
            split_info.n_samples_right = n_samples - best_n_samples_left

            # We recompute best values here but it's cheap
            split_info.value_left = compute_node_value(
                split_info.sum_gradient_left, split_info.sum_hessian_left,
                lower_bound, upper_bound, self.l2_regularization)

            split_info.value_right = compute_node_value(
                split_info.sum_gradient_right, split_info.sum_hessian_right,
                lower_bound, upper_bound, self.l2_regularization)

    cdef void _find_best_bin_to_split_right_to_left(
            self,
            unsigned int feature_idx,
            const hist_struct [:, ::1] histograms,  # IN
            unsigned int n_samples,
            Y_DTYPE_C sum_gradients,
            Y_DTYPE_C sum_hessians,
            Y_DTYPE_C value,
            signed char monotonic_cst,
            Y_DTYPE_C lower_bound,
            Y_DTYPE_C upper_bound,
            split_info_struct * split_info) noexcept nogil:  # OUT
        """Find best bin to split on for a given feature.

        Splits that do not satisfy the splitting constraints
        (min_gain_to_split, etc.) are discarded here.

        We scan node from right to left. This version is only called when
        there are missing values. Missing values are assigned to the left
        child.

        If no missing value are present in the data this method isn't called
        since only calling _find_best_bin_to_split_left_to_right is enough.
        """

        cdef:
            unsigned int bin_idx
            unsigned int n_samples_left
            unsigned int n_samples_right
            unsigned int n_samples_ = n_samples
            Y_DTYPE_C sum_hessian_left
            Y_DTYPE_C sum_hessian_right
            Y_DTYPE_C sum_gradient_left
            Y_DTYPE_C sum_gradient_right
            Y_DTYPE_C loss_current_node
            Y_DTYPE_C gain
            unsigned int start = self.n_bins_non_missing[feature_idx] - 2
            uint8_t found_better_split = False

            Y_DTYPE_C best_sum_hessian_left
            Y_DTYPE_C best_sum_gradient_left
            unsigned int best_bin_idx
            unsigned int best_n_samples_left
            Y_DTYPE_C best_gain = split_info.gain  # computed during previous scan
            hist_struct hist

        sum_gradient_right, sum_hessian_right = 0., 0.
        n_samples_right = 0

        loss_current_node = _loss_from_value(value, sum_gradients)

        for bin_idx in range(start, -1, -1):
            hist = histograms[feature_idx, bin_idx + 1]
            n_samples_right += hist.count
            n_samples_left = n_samples_ - n_samples_right

            if self.hessians_are_constant:
                sum_hessian_right += hist.count
            else:
                sum_hessian_right += \
                    hist.sum_hessians
            sum_hessian_left = sum_hessians - sum_hessian_right

            sum_gradient_right += \
                hist.sum_gradients
            sum_gradient_left = sum_gradients - sum_gradient_right

            if n_samples_right < self.min_samples_leaf:
                continue
            if n_samples_left < self.min_samples_leaf:
                # won't get any better
                break

            if sum_hessian_right < self.min_hessian_to_split:
                continue
            if sum_hessian_left < self.min_hessian_to_split:
                # won't get any better (hessians are > 0 since loss is convex)
                break

            gain = _split_gain(sum_gradient_left, sum_hessian_left,
                               sum_gradient_right, sum_hessian_right,
                               loss_current_node,
                               monotonic_cst,
                               lower_bound,
                               upper_bound,
                               self.l2_regularization)

            if gain > best_gain and gain > self.min_gain_to_split:
                found_better_split = True
                best_gain = gain
                best_bin_idx = bin_idx
                best_sum_gradient_left = sum_gradient_left
                best_sum_hessian_left = sum_hessian_left
                best_n_samples_left = n_samples_left

        if found_better_split:
            split_info.gain = best_gain
            split_info.bin_idx = best_bin_idx
            # we scan from right to left so missing values go to the left
            split_info.missing_go_to_left = True
            split_info.sum_gradient_left = best_sum_gradient_left
            split_info.sum_gradient_right = sum_gradients - best_sum_gradient_left
            split_info.sum_hessian_left = best_sum_hessian_left
            split_info.sum_hessian_right = sum_hessians - best_sum_hessian_left
            split_info.n_samples_left = best_n_samples_left
            split_info.n_samples_right = n_samples - best_n_samples_left

            # We recompute best values here but it's cheap
            split_info.value_left = compute_node_value(
                split_info.sum_gradient_left, split_info.sum_hessian_left,
                lower_bound, upper_bound, self.l2_regularization)

            split_info.value_right = compute_node_value(
                split_info.sum_gradient_right, split_info.sum_hessian_right,
                lower_bound, upper_bound, self.l2_regularization)

    cdef void _find_best_bin_to_split_category(
            self,
            unsigned int feature_idx,
            uint8_t has_missing_values,
            const hist_struct [:, ::1] histograms,  # IN
            unsigned int n_samples,
            Y_DTYPE_C sum_gradients,
            Y_DTYPE_C sum_hessians,
            Y_DTYPE_C value,
            char monotonic_cst,
            Y_DTYPE_C lower_bound,
            Y_DTYPE_C upper_bound,
            split_info_struct * split_info) noexcept nogil:  # OUT
        """Find best split for categorical features.

        Categories are first sorted according to their variance, and then
        a scan is performed as if categories were ordered quantities.

        Ref: "On Grouping for Maximum Homogeneity", Walter D. Fisher
        """

        cdef:
            unsigned int bin_idx
            unsigned int n_bins_non_missing = self.n_bins_non_missing[feature_idx]
            unsigned int missing_values_bin_idx = self.missing_values_bin_idx
            categorical_info * cat_infos
            unsigned int sorted_cat_idx
            unsigned int n_used_bins = 0
            int [2] scan_direction
            int direction = 0
            int best_direction = 0
            unsigned int middle
            unsigned int i
            const hist_struct[::1] feature_hist = histograms[feature_idx, :]
            hist_struct hist
            Y_DTYPE_C sum_gradients_bin
            Y_DTYPE_C sum_hessians_bin
            Y_DTYPE_C loss_current_node
            Y_DTYPE_C sum_gradient_left, sum_hessian_left
            Y_DTYPE_C sum_gradient_right, sum_hessian_right
            unsigned int n_samples_left, n_samples_right
            Y_DTYPE_C gain
            Y_DTYPE_C best_gain = -1.0
            uint8_t found_better_split = False
            Y_DTYPE_C best_sum_hessian_left
            Y_DTYPE_C best_sum_gradient_left
            unsigned int best_n_samples_left
            unsigned int best_cat_infos_thresh
            # Reduces the effect of noises in categorical features,
            # especially for categories with few data. Called cat_smooth in
            # LightGBM. TODO: Make this user adjustable?
            Y_DTYPE_C MIN_CAT_SUPPORT = 10.
            # this is equal to 1 for losses where hessians are constant
            Y_DTYPE_C support_factor = n_samples / sum_hessians

        # Details on the split finding:
        # We first order categories by their sum_gradients / sum_hessians
        # values, and we exclude categories that don't respect MIN_CAT_SUPPORT
        # from this sorted array. Missing values are treated just like any
        # other category. The low-support categories will always be mapped to
        # the right child. We scan the sorted categories array from left to
        # right and from right to left, and we stop at the middle.

        # Considering ordered categories A B C D, with E being a low-support
        # category: A B C D
        #              ^
        #           midpoint
        # The scans will consider the following split-points:
        # * left to right:
        #   A - B C D E
        #   A B - C D E
        # * right to left:
        #   D - A B C E
        #   C D - A B E

        # Note that since we stop at the middle and since low-support
        # categories (E) are always mapped to the right, the following splits
        # aren't considered:
        # A E - B C D
        # D E - A B C
        # Basically, we're forcing E to always be mapped to the child that has
        # *at least half of the categories* (and this child is always the right
        # child, by convention).

        # Also note that if we scanned in only one direction (e.g. left to
        # right), we would only consider the following splits:
        # A - B C D E
        # A B - C D E
        # A B C - D E
        # and thus we would be missing on D - A B C E and on C D - A B E

        cat_infos = <categorical_info *> malloc(
            (n_bins_non_missing + has_missing_values) * sizeof(categorical_info))

        # fill cat_infos while filtering out categories based on MIN_CAT_SUPPORT
        for bin_idx in range(n_bins_non_missing):
            hist = feature_hist[bin_idx]
            if self.hessians_are_constant:
                sum_hessians_bin = hist.count
            else:
                sum_hessians_bin = hist.sum_hessians
            if sum_hessians_bin * support_factor >= MIN_CAT_SUPPORT:
                cat_infos[n_used_bins].bin_idx = bin_idx
                sum_gradients_bin = hist.sum_gradients

                cat_infos[n_used_bins].value = (
                    sum_gradients_bin / (sum_hessians_bin + MIN_CAT_SUPPORT)
                )
                n_used_bins += 1

        # Also add missing values bin so that nans are considered as a category
        if has_missing_values:
            hist = feature_hist[missing_values_bin_idx]
            if self.hessians_are_constant:
                sum_hessians_bin = hist.count
            else:
                sum_hessians_bin = hist.sum_hessians
            if sum_hessians_bin * support_factor >= MIN_CAT_SUPPORT:
                cat_infos[n_used_bins].bin_idx = missing_values_bin_idx
                sum_gradients_bin = (
                    hist.sum_gradients
                )

                cat_infos[n_used_bins].value = (
                    sum_gradients_bin / (sum_hessians_bin + MIN_CAT_SUPPORT)
                )
                n_used_bins += 1

        # not enough categories to form a split
        if n_used_bins <= 1:
            free(cat_infos)
            return

        qsort(cat_infos, n_used_bins, sizeof(categorical_info),
              compare_cat_infos)

        loss_current_node = _loss_from_value(value, sum_gradients)

        scan_direction[0], scan_direction[1] = 1, -1
        for direction in scan_direction:
            if direction == 1:
                middle = (n_used_bins + 1) // 2
            else:
                middle = (n_used_bins + 1) // 2 - 1

            # The categories we'll consider will go to the left child
            sum_gradient_left, sum_hessian_left = 0., 0.
            n_samples_left = 0

            for i in range(middle):
                sorted_cat_idx = i if direction == 1 else n_used_bins - 1 - i
                bin_idx = cat_infos[sorted_cat_idx].bin_idx
                hist = feature_hist[bin_idx]

                n_samples_left += hist.count
                n_samples_right = n_samples - n_samples_left

                if self.hessians_are_constant:
                    sum_hessian_left += hist.count
                else:
                    sum_hessian_left += hist.sum_hessians
                sum_hessian_right = sum_hessians - sum_hessian_left

                sum_gradient_left += hist.sum_gradients
                sum_gradient_right = sum_gradients - sum_gradient_left

                if (
                    n_samples_left < self.min_samples_leaf or
                    sum_hessian_left < self.min_hessian_to_split
                ):
                    continue
                if (
                    n_samples_right < self.min_samples_leaf or
                    sum_hessian_right < self.min_hessian_to_split
                ):
                    break

                gain = _split_gain(sum_gradient_left, sum_hessian_left,
                                   sum_gradient_right, sum_hessian_right,
                                   loss_current_node, monotonic_cst,
                                   lower_bound, upper_bound,
                                   self.l2_regularization)
                if gain > best_gain and gain > self.min_gain_to_split:
                    found_better_split = True
                    best_gain = gain
                    best_cat_infos_thresh = sorted_cat_idx
                    best_sum_gradient_left = sum_gradient_left
                    best_sum_hessian_left = sum_hessian_left
                    best_n_samples_left = n_samples_left
                    best_direction = direction

        if found_better_split:
            split_info.gain = best_gain

            # split_info.bin_idx is unused for categorical splits: left_cat_bitset
            # is used instead and set below
            split_info.bin_idx = 0

            split_info.sum_gradient_left = best_sum_gradient_left
            split_info.sum_gradient_right = sum_gradients - best_sum_gradient_left
            split_info.sum_hessian_left = best_sum_hessian_left
            split_info.sum_hessian_right = sum_hessians - best_sum_hessian_left
            split_info.n_samples_left = best_n_samples_left
            split_info.n_samples_right = n_samples - best_n_samples_left

            # We recompute best values here but it's cheap
            split_info.value_left = compute_node_value(
                split_info.sum_gradient_left, split_info.sum_hessian_left,
                lower_bound, upper_bound, self.l2_regularization)

            split_info.value_right = compute_node_value(
                split_info.sum_gradient_right, split_info.sum_hessian_right,
                lower_bound, upper_bound, self.l2_regularization)

            # create bitset with values from best_cat_infos_thresh
            init_bitset(split_info.left_cat_bitset)
            if best_direction == 1:
                for sorted_cat_idx in range(best_cat_infos_thresh + 1):
                    bin_idx = cat_infos[sorted_cat_idx].bin_idx
                    set_bitset(split_info.left_cat_bitset, bin_idx)
            else:
                for sorted_cat_idx in range(n_used_bins - 1, best_cat_infos_thresh - 1, -1):
                    bin_idx = cat_infos[sorted_cat_idx].bin_idx
                    set_bitset(split_info.left_cat_bitset, bin_idx)

            if has_missing_values:
                split_info.missing_go_to_left = in_bitset(
                    split_info.left_cat_bitset, missing_values_bin_idx)

        free(cat_infos)


cdef int compare_cat_infos(const void * a, const void * b) noexcept nogil:
    return -1 if (<categorical_info *>a).value < (<categorical_info *>b).value else 1

cdef inline Y_DTYPE_C _split_gain(
        Y_DTYPE_C sum_gradient_left,
        Y_DTYPE_C sum_hessian_left,
        Y_DTYPE_C sum_gradient_right,
        Y_DTYPE_C sum_hessian_right,
        Y_DTYPE_C loss_current_node,
        signed char monotonic_cst,
        Y_DTYPE_C lower_bound,
        Y_DTYPE_C upper_bound,
        Y_DTYPE_C l2_regularization) noexcept nogil:
    """Loss reduction

    Compute the reduction in loss after taking a split, compared to keeping
    the node a leaf of the tree.

    See Equation 7 of:
    :arxiv:`T. Chen, C. Guestrin, (2016) XGBoost: A Scalable Tree Boosting System,
    <1603.02754>.`
    """
    cdef:
        Y_DTYPE_C gain
        Y_DTYPE_C value_left
        Y_DTYPE_C value_right

    # Compute values of potential left and right children
    value_left = compute_node_value(sum_gradient_left, sum_hessian_left,
                                    lower_bound, upper_bound,
                                    l2_regularization)
    value_right = compute_node_value(sum_gradient_right, sum_hessian_right,
                                     lower_bound, upper_bound,
                                     l2_regularization)

    if ((monotonic_cst == MonotonicConstraint.POS and value_left > value_right) or
            (monotonic_cst == MonotonicConstraint.NEG and value_left < value_right)):
        # don't consider this split since it does not respect the monotonic
        # constraints. Note that these comparisons need to be done on values
        # that have already been clipped to take the monotonic constraints into
        # account (if any).
        return -1

    gain = loss_current_node
    gain -= _loss_from_value(value_left, sum_gradient_left)
    gain -= _loss_from_value(value_right, sum_gradient_right)
    # Note that for the gain to be correct (and for min_gain_to_split to work
    # as expected), we need all values to be bounded (current node, left child
    # and right child).

    return gain

cdef inline Y_DTYPE_C _loss_from_value(
        Y_DTYPE_C value,
        Y_DTYPE_C sum_gradient) noexcept nogil:
    """Return loss of a node from its (bounded) value

    See Equation 6 of:
    :arxiv:`T. Chen, C. Guestrin, (2016) XGBoost: A Scalable Tree Boosting System,
    <1603.02754>.`
    """
    return sum_gradient * value

cdef inline uint8_t sample_goes_left(
        uint8_t missing_go_to_left,
        uint8_t missing_values_bin_idx,
        X_BINNED_DTYPE_C split_bin_idx,
        X_BINNED_DTYPE_C bin_value,
        uint8_t is_categorical,
        BITSET_DTYPE_C left_cat_bitset) noexcept nogil:
    """Helper to decide whether sample should go to left or right child."""

    if is_categorical:
        # note: if any, missing values are encoded in left_cat_bitset
        return in_bitset(left_cat_bitset, bin_value)
    else:
        return (
            (
                missing_go_to_left and
                bin_value == missing_values_bin_idx
            )
            or (
                bin_value <= split_bin_idx
            ))


cpdef inline Y_DTYPE_C compute_node_value(
        Y_DTYPE_C sum_gradient,
        Y_DTYPE_C sum_hessian,
        Y_DTYPE_C lower_bound,
        Y_DTYPE_C upper_bound,
        Y_DTYPE_C l2_regularization) noexcept nogil:
    """Compute a node's value.

    The value is capped in the [lower_bound, upper_bound] interval to respect
    monotonic constraints. Shrinkage is ignored.

    See Equation 5 of:
    :arxiv:`T. Chen, C. Guestrin, (2016) XGBoost: A Scalable Tree Boosting System,
    <1603.02754>.`
    """

    cdef:
        Y_DTYPE_C value

    value = -sum_gradient / (sum_hessian + l2_regularization + 1e-15)

    if value < lower_bound:
        value = lower_bound
    elif value > upper_bound:
        value = upper_bound

    return value
