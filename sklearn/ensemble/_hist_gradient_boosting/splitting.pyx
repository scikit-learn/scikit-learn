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
from libc.string cimport memcpy, memset

from sklearn.utils._bitset cimport BITSET_DTYPE_C, BITSET_INNER_DTYPE_C
from sklearn.utils._bitset cimport in_bitset, init_bitset, set_bitset
from sklearn.utils._typedefs cimport uint8_t
from sklearn.ensemble._hist_gradient_boosting.common cimport X_BINNED_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport Y_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport hist_struct
from sklearn.ensemble._hist_gradient_boosting.common cimport MonotonicConstraint


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


cdef struct vector_split_info_struct:
    # Per-feature best split for trees with vector leaves. The per-output sums
    # of the left child are stored in a separate buffer.
    Y_DTYPE_C gain
    int feature_idx
    unsigned int bin_idx
    uint8_t missing_go_to_left
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
        unsigned int [::1] all_features
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
        self.all_features = np.arange(self.n_features, dtype=np.uint32)
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
            int n_split_candidates
            split_info_struct split_info
            split_info_struct * split_infos
            const uint8_t [::1] has_missing_values = self.has_missing_values
            const uint8_t [::1] is_categorical = self.is_categorical
            const signed char [::1] monotonic_cst = self.monotonic_cst
            int n_threads = self.n_threads
            Y_DTYPE_C feature_fraction_per_split = self.feature_fraction_per_split
            const unsigned int [:] split_features
            uint8_t missing_go_to_left

        if allowed_features is None:
            split_features = self.all_features
        else:
            split_features = allowed_features

        n_split_candidates = split_features.shape[0]

        if feature_fraction_per_split < 1.0:
            # We do all random sampling before the nogil and make sure that we sample
            # exactly n_split_candidates >= 1 features.
            n_split_candidates = max(
                1,
                int(ceil(feature_fraction_per_split * n_split_candidates)),
            )
            split_features = self.rng.choice(
                split_features, n_split_candidates, replace=False,
            )
            split_features = np.sort(split_features)

        with nogil:

            split_infos = <split_info_struct *> malloc(
                n_split_candidates * sizeof(split_info_struct))

            # split_info_idx is index of split_infos of size n_split_candidates.
            # features_idx is the index of the feature column in X.
            for split_info_idx in prange(n_split_candidates, schedule='static',
                                         num_threads=n_threads):
                feature_idx = split_features[split_info_idx]

                split_infos[split_info_idx].feature_idx = feature_idx

                # For each feature, find best bin to split on
                # Start with a gain of -1 if no better split is found, that
                # means one of the constraints isn't respected
                # (min_samples_leaf, etc.) and the grower will later turn the
                # node into a leaf.
                split_infos[split_info_idx].gain = -1
                split_infos[split_info_idx].is_categorical = is_categorical[feature_idx]

                if is_categorical[feature_idx]:
                    self._find_best_bin_to_split_category(
                        feature_idx, has_missing_values[feature_idx],
                        histograms, n_samples, sum_gradients, sum_hessians,
                        value, monotonic_cst[feature_idx], lower_bound,
                        upper_bound, &split_infos[split_info_idx])
                else:
                    # We scan bins from left to right and, if there are any
                    # missing values, we scan a second time with missing
                    # values in the left child. This way, we can consider
                    # whichever case yields the best gain: either missing
                    # values go to the right or to the left. See algo 3 from
                    # the XGBoost paper
                    # https://arxiv.org/abs/1603.02754
                    # Note: for the categorical features above, this isn't
                    # needed since missing values are considered a native
                    # category.
                    for missing_go_to_left in range(has_missing_values[feature_idx] + 1):

                        self._find_best_bin_to_split_numerical(
                            feature_idx, has_missing_values[feature_idx],
                            histograms, n_samples, sum_gradients, sum_hessians,
                            value, monotonic_cst[feature_idx],
                            lower_bound, upper_bound, missing_go_to_left,
                            &split_infos[split_info_idx])

            # then compute best possible split among all features
            # split_info is set to the best of split_infos
            best_split_info_idx = self._find_best_feature_to_split_helper(
                split_infos, n_split_candidates
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
        int n_split_candidates,
    ) noexcept nogil:
        """Return the index of split_infos with the best feature split."""
        cdef:
            int split_info_idx
            int best_split_info_idx = 0

        for split_info_idx in range(1, n_split_candidates):
            if (split_infos[split_info_idx].gain > split_infos[best_split_info_idx].gain):
                best_split_info_idx = split_info_idx
        return best_split_info_idx

    cdef void _find_best_bin_to_split_numerical(
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
            uint8_t missing_go_to_left,
            split_info_struct * split_info) noexcept nogil:  # OUT
        """Find best bin to split on for a given feature.

        Splits that do not satisfy the splitting constraints
        (min_gain_to_split, etc.) are discarded here.

        We scan node from left to right. When missing_go_to_left is True,
        the scan starts by adding the missing-value bin to the left child
        before scanning the non-missing bins.
        """
        cdef:
            int bin_idx
            unsigned int hist_bin_idx
            int n_bins_non_missing = self.n_bins_non_missing[feature_idx]
            unsigned int n_samples_left
            unsigned int n_samples_right
            unsigned int n_samples_ = n_samples
            int start
            int end
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
            Y_DTYPE_C best_gain = split_info.gain
            hist_struct hist

        sum_gradient_left, sum_hessian_left = 0., 0.
        n_samples_left = 0

        loss_current_node = _loss_from_value(value, sum_gradients)

        # When missing_go_to_left is true, bin_idx == -1 maps to the missing
        # bin. Missing values are always in the last bin.
        start = -missing_go_to_left
        # Do not put the last non-missing bin in the left child unless
        # missing values can still go to the right child.
        end = n_bins_non_missing - 1 + has_missing_values * (1 - missing_go_to_left)

        for bin_idx in range(start, end):
            if bin_idx == -1:
                hist_bin_idx = self.missing_values_bin_idx
            else:
                hist_bin_idx = bin_idx

            hist = histograms[feature_idx, hist_bin_idx]
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

            if bin_idx == -1:
                # This iteration only adds missing values to the left child.
                # The corresponding all-non-missing-values-on-one-side split
                # is already considered when missing_go_to_left is False.
                continue

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
            split_info.missing_go_to_left = missing_go_to_left
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

    def find_node_split_vector(
            Splitter self,
            unsigned int n_samples,
            const Y_DTYPE_C [:, :, ::1] histograms,  # IN
            const Y_DTYPE_C [::1] sum_gradients,
            const Y_DTYPE_C [::1] sum_hessians,
            const Y_DTYPE_C [::1] value,
            const unsigned int [:] allowed_features=None,
            ):
        """Find the best split of a node of a tree with vector leaves.

        The gain of a split is the sum over the outputs of the per-output
        gains (see `_split_gain`), so that a single tree structure is shared
        by all the outputs while each leaf stores one value per output.

        Monotonic constraints are not supported.

        Parameters
        ----------
        n_samples : int
            The number of samples at the node.
        histograms : ndarray of float64 of shape (n_features, n_bins, n_stats)
            The histograms of the current node, as computed by
            `VectorHistogramBuilder`.
        sum_gradients : ndarray of shape (n_outputs,)
            The per-output sum of the gradients of the samples at the node.
        sum_hessians : ndarray of shape (n_outputs,)
            The per-output sum of the hessians of the samples at the node.
        value : ndarray of shape (n_outputs,)
            The values of the current node.
        allowed_features : None or ndarray, dtype=np.uint32
            Indices of the features that are allowed by interaction constraints to be
            split.

        Returns
        -------
        best_split_info : SplitInfo
            The info about the best possible split among all features. The
            sums of gradients / hessians and the children values are arrays
            of shape (n_outputs,).
        """
        cdef:
            int feature_idx
            int split_info_idx
            int best_split_info_idx
            int n_split_candidates
            unsigned int k
            unsigned int n_outputs = sum_gradients.shape[0]
            unsigned int n_stats = histograms.shape[2]
            vector_split_info_struct split_info
            vector_split_info_struct * split_infos
            Y_DTYPE_C [:, ::1] left_stats
            Y_DTYPE_C [:, ::1] best_left_stats
            Y_DTYPE_C loss_current_node = 0.
            const uint8_t [::1] has_missing_values = self.has_missing_values
            const uint8_t [::1] is_categorical = self.is_categorical
            int n_threads = self.n_threads
            Y_DTYPE_C feature_fraction_per_split = self.feature_fraction_per_split
            const unsigned int [:] split_features
            uint8_t missing_go_to_left

        if allowed_features is None:
            split_features = self.all_features
        else:
            split_features = allowed_features

        n_split_candidates = split_features.shape[0]

        if feature_fraction_per_split < 1.0:
            n_split_candidates = max(
                1,
                int(ceil(feature_fraction_per_split * n_split_candidates)),
            )
            split_features = self.rng.choice(
                split_features, n_split_candidates, replace=False,
            )
            split_features = np.sort(split_features)

        for k in range(n_outputs):
            loss_current_node += _loss_from_value(value[k], sum_gradients[k])

        # Per-candidate scratch buffer for the running sums of the left child
        # and per-candidate statistics of the left child of the best split.
        left_stats = np.empty((n_split_candidates, n_stats), dtype=np.float64)
        best_left_stats = np.zeros((n_split_candidates, n_stats), dtype=np.float64)

        with nogil:
            split_infos = <vector_split_info_struct *> malloc(
                n_split_candidates * sizeof(vector_split_info_struct))

            for split_info_idx in prange(n_split_candidates, schedule='static',
                                         num_threads=n_threads):
                feature_idx = split_features[split_info_idx]
                split_infos[split_info_idx].feature_idx = feature_idx
                split_infos[split_info_idx].gain = -1
                split_infos[split_info_idx].is_categorical = is_categorical[feature_idx]

                if is_categorical[feature_idx]:
                    self._find_best_bin_to_split_category_vector(
                        feature_idx, has_missing_values[feature_idx],
                        histograms, n_samples, &sum_gradients[0],
                        &sum_hessians[0], n_outputs, loss_current_node,
                        &split_infos[split_info_idx],
                        &left_stats[split_info_idx, 0],
                        &best_left_stats[split_info_idx, 0])
                else:
                    for missing_go_to_left in range(has_missing_values[feature_idx] + 1):
                        self._find_best_bin_to_split_numerical_vector(
                            feature_idx, has_missing_values[feature_idx],
                            histograms, n_samples, &sum_gradients[0],
                            &sum_hessians[0], n_outputs, loss_current_node,
                            missing_go_to_left, &split_infos[split_info_idx],
                            &left_stats[split_info_idx, 0],
                            &best_left_stats[split_info_idx, 0])

            best_split_info_idx = 0
            for split_info_idx in range(1, n_split_candidates):
                if split_infos[split_info_idx].gain > split_infos[best_split_info_idx].gain:
                    best_split_info_idx = split_info_idx
            split_info = split_infos[best_split_info_idx]

        left = np.asarray(best_left_stats[best_split_info_idx])
        n_samples_left = int(left[n_stats - 1])
        sum_gradient_left = left[:n_outputs].copy()
        if self.hessians_are_constant:
            sum_hessian_left = np.full(n_outputs, left[n_stats - 1])
        else:
            sum_hessian_left = left[n_outputs:2 * n_outputs].copy()
        sum_gradient_right = np.asarray(sum_gradients) - sum_gradient_left
        sum_hessian_right = np.asarray(sum_hessians) - sum_hessian_left
        l2_regularization = self.l2_regularization
        value_left = -sum_gradient_left / (sum_hessian_left + l2_regularization + 1e-15)
        value_right = -sum_gradient_right / (
            sum_hessian_right + l2_regularization + 1e-15
        )

        out = SplitInfo(
            split_info.gain,
            split_info.feature_idx,
            split_info.bin_idx,
            split_info.missing_go_to_left,
            sum_gradient_left,
            sum_hessian_left,
            sum_gradient_right,
            sum_hessian_right,
            n_samples_left,
            n_samples - n_samples_left,
            value_left,
            value_right,
            split_info.is_categorical,
            None,  # left_cat_bitset will only be set if the split is categorical
        )
        if split_info.is_categorical:
            out.left_cat_bitset = np.asarray(split_info.left_cat_bitset, dtype=np.uint32)

        free(split_infos)
        return out

    cdef void _find_best_bin_to_split_numerical_vector(
            Splitter self,
            unsigned int feature_idx,
            uint8_t has_missing_values,
            const Y_DTYPE_C [:, :, ::1] histograms,  # IN
            unsigned int n_samples,
            const Y_DTYPE_C * sum_gradients,
            const Y_DTYPE_C * sum_hessians,
            unsigned int n_outputs,
            Y_DTYPE_C loss_current_node,
            uint8_t missing_go_to_left,
            vector_split_info_struct * split_info,  # OUT
            Y_DTYPE_C * left_stats,  # scratch buffer of size n_stats
            Y_DTYPE_C * best_left_stats) noexcept nogil:  # OUT
        """Vector leaves counterpart of _find_best_bin_to_split_numerical."""
        cdef:
            int bin_idx
            unsigned int hist_bin_idx
            unsigned int s
            unsigned int k
            unsigned int n_stats = histograms.shape[2]
            unsigned int count_idx = n_stats - 1
            int n_bins_non_missing = self.n_bins_non_missing[feature_idx]
            uint8_t hessians_are_constant = self.hessians_are_constant
            unsigned int n_samples_left
            unsigned int n_samples_right
            int start
            int end
            Y_DTYPE_C total_sum_hessians = 0.
            Y_DTYPE_C sum_hessian_left
            Y_DTYPE_C sum_hessian_right
            Y_DTYPE_C gain
            Y_DTYPE_C best_gain = split_info.gain
            uint8_t found_better_split = False
            unsigned int best_bin_idx
            const Y_DTYPE_C * hist_bin

        memset(left_stats, 0, n_stats * sizeof(Y_DTYPE_C))
        for k in range(n_outputs):
            total_sum_hessians += sum_hessians[k]

        start = -missing_go_to_left
        end = n_bins_non_missing - 1 + has_missing_values * (1 - missing_go_to_left)

        for bin_idx in range(start, end):
            if bin_idx == -1:
                hist_bin_idx = self.missing_values_bin_idx
            else:
                hist_bin_idx = bin_idx

            hist_bin = &histograms[feature_idx, hist_bin_idx, 0]
            for s in range(n_stats):
                left_stats[s] += hist_bin[s]

            if bin_idx == -1:
                continue

            n_samples_left = <unsigned int> left_stats[count_idx]
            n_samples_right = n_samples - n_samples_left

            if n_samples_left < self.min_samples_leaf:
                continue
            if n_samples_right < self.min_samples_leaf:
                break

            sum_hessian_left = _total_hessian(
                left_stats, n_outputs, hessians_are_constant)
            sum_hessian_right = total_sum_hessians - sum_hessian_left
            if sum_hessian_left < self.min_hessian_to_split:
                continue
            if sum_hessian_right < self.min_hessian_to_split:
                break

            gain = _vector_split_gain(
                left_stats, sum_gradients, sum_hessians, n_outputs,
                hessians_are_constant, loss_current_node,
                self.l2_regularization)

            if gain > best_gain and gain > self.min_gain_to_split:
                found_better_split = True
                best_gain = gain
                best_bin_idx = bin_idx
                memcpy(best_left_stats, left_stats, n_stats * sizeof(Y_DTYPE_C))

        if found_better_split:
            split_info.gain = best_gain
            split_info.bin_idx = best_bin_idx
            split_info.missing_go_to_left = missing_go_to_left

    cdef void _find_best_bin_to_split_category_vector(
            Splitter self,
            unsigned int feature_idx,
            uint8_t has_missing_values,
            const Y_DTYPE_C [:, :, ::1] histograms,  # IN
            unsigned int n_samples,
            const Y_DTYPE_C * sum_gradients,
            const Y_DTYPE_C * sum_hessians,
            unsigned int n_outputs,
            Y_DTYPE_C loss_current_node,
            vector_split_info_struct * split_info,  # OUT
            Y_DTYPE_C * left_stats,  # scratch buffer of size n_stats
            Y_DTYPE_C * best_left_stats) noexcept nogil:  # OUT
        """Vector leaves counterpart of _find_best_bin_to_split_category.

        With a single output, categories are sorted by their value as in the
        scalar case. With several outputs, there is no natural ordering of
        the categories: they are sorted by the projection of their
        (smoothed) value vectors onto the first principal direction of these
        vectors (weighted by the hessians), computed by power iteration. A
        Fisher-style scan is then performed over this ordering.
        """
        cdef:
            unsigned int bin_idx
            unsigned int n_bins_non_missing = self.n_bins_non_missing[feature_idx]
            unsigned int missing_values_bin_idx = self.missing_values_bin_idx
            unsigned int n_candidate_bins = n_bins_non_missing + has_missing_values
            unsigned int n_stats = histograms.shape[2]
            unsigned int count_idx = n_stats - 1
            uint8_t hessians_are_constant = self.hessians_are_constant
            categorical_info * cat_infos
            Y_DTYPE_C * cat_values
            Y_DTYPE_C * cat_weights
            Y_DTYPE_C * direction_vec
            Y_DTYPE_C * new_direction_vec
            unsigned int sorted_cat_idx
            unsigned int n_used_bins = 0
            int [2] scan_direction
            int direction = 0
            int best_direction = 0
            unsigned int middle
            unsigned int i
            unsigned int s
            unsigned int k
            unsigned int c
            unsigned int _
            unsigned int best_c
            const Y_DTYPE_C * hist_bin
            Y_DTYPE_C h_bin
            Y_DTYPE_C sum_hessians_bin
            Y_DTYPE_C total_sum_hessians = 0.
            Y_DTYPE_C sum_hessian_left, sum_hessian_right
            unsigned int n_samples_left, n_samples_right
            Y_DTYPE_C gain
            Y_DTYPE_C best_gain = -1.0
            Y_DTYPE_C tmp
            Y_DTYPE_C norm
            Y_DTYPE_C total_weight
            uint8_t found_better_split = False
            unsigned int best_cat_infos_thresh
            Y_DTYPE_C MIN_CAT_SUPPORT = 10.
            Y_DTYPE_C support_factor
            unsigned int N_POWER_ITERATIONS = 10

        for k in range(n_outputs):
            total_sum_hessians += sum_hessians[k]
        support_factor = n_samples / total_sum_hessians

        cat_infos = <categorical_info *> malloc(
            n_candidate_bins * sizeof(categorical_info))
        cat_values = <Y_DTYPE_C *> malloc(
            n_candidate_bins * n_outputs * sizeof(Y_DTYPE_C))
        cat_weights = <Y_DTYPE_C *> malloc(n_candidate_bins * sizeof(Y_DTYPE_C))
        direction_vec = <Y_DTYPE_C *> malloc(n_outputs * sizeof(Y_DTYPE_C))
        new_direction_vec = <Y_DTYPE_C *> malloc(n_outputs * sizeof(Y_DTYPE_C))

        # fill cat_infos while filtering out categories based on MIN_CAT_SUPPORT
        for i in range(n_candidate_bins):
            if i < n_bins_non_missing:
                bin_idx = i
            else:
                bin_idx = missing_values_bin_idx
            hist_bin = &histograms[feature_idx, bin_idx, 0]
            sum_hessians_bin = _total_hessian(
                hist_bin, n_outputs, hessians_are_constant)
            if sum_hessians_bin * support_factor >= MIN_CAT_SUPPORT:
                cat_infos[n_used_bins].bin_idx = bin_idx
                cat_weights[n_used_bins] = sum_hessians_bin
                for k in range(n_outputs):
                    if hessians_are_constant:
                        h_bin = hist_bin[count_idx]
                    else:
                        h_bin = hist_bin[n_outputs + k]
                    cat_values[n_used_bins * n_outputs + k] = (
                        hist_bin[k] / (h_bin + MIN_CAT_SUPPORT)
                    )
                n_used_bins += 1

        if n_used_bins <= 1:
            free(cat_infos)
            free(cat_values)
            free(cat_weights)
            free(direction_vec)
            free(new_direction_vec)
            return

        if n_outputs == 1:
            for c in range(n_used_bins):
                cat_infos[c].value = cat_values[c]
        else:
            # Center the category values (weighted by their hessians).
            total_weight = 0.
            for c in range(n_used_bins):
                total_weight += cat_weights[c]
            for k in range(n_outputs):
                tmp = 0.
                for c in range(n_used_bins):
                    tmp += cat_weights[c] * cat_values[c * n_outputs + k]
                tmp /= total_weight
                for c in range(n_used_bins):
                    cat_values[c * n_outputs + k] -= tmp

            # Initialize the direction with the category of largest weighted
            # squared norm and run a few power iterations on the weighted
            # covariance matrix of the centered values.
            best_c = 0
            norm = -1.
            for c in range(n_used_bins):
                tmp = 0.
                for k in range(n_outputs):
                    tmp += cat_values[c * n_outputs + k] ** 2
                tmp *= cat_weights[c]
                if tmp > norm:
                    norm = tmp
                    best_c = c
            for k in range(n_outputs):
                direction_vec[k] = cat_values[best_c * n_outputs + k]

            for _ in range(N_POWER_ITERATIONS):
                memset(new_direction_vec, 0, n_outputs * sizeof(Y_DTYPE_C))
                for c in range(n_used_bins):
                    tmp = 0.
                    for k in range(n_outputs):
                        tmp += cat_values[c * n_outputs + k] * direction_vec[k]
                    tmp *= cat_weights[c]
                    for k in range(n_outputs):
                        new_direction_vec[k] += tmp * cat_values[c * n_outputs + k]
                norm = 0.
                for k in range(n_outputs):
                    norm += new_direction_vec[k] ** 2
                if norm <= 0.:
                    break
                norm = norm ** 0.5
                for k in range(n_outputs):
                    direction_vec[k] = new_direction_vec[k] / norm

            for c in range(n_used_bins):
                tmp = 0.
                for k in range(n_outputs):
                    tmp += cat_values[c * n_outputs + k] * direction_vec[k]
                cat_infos[c].value = tmp

        qsort(cat_infos, n_used_bins, sizeof(categorical_info),
              compare_cat_infos)

        scan_direction[0], scan_direction[1] = 1, -1
        for direction in scan_direction:
            if direction == 1:
                middle = (n_used_bins + 1) // 2
            else:
                middle = (n_used_bins + 1) // 2 - 1

            memset(left_stats, 0, n_stats * sizeof(Y_DTYPE_C))

            for i in range(middle):
                sorted_cat_idx = i if direction == 1 else n_used_bins - 1 - i
                bin_idx = cat_infos[sorted_cat_idx].bin_idx
                hist_bin = &histograms[feature_idx, bin_idx, 0]
                for s in range(n_stats):
                    left_stats[s] += hist_bin[s]

                n_samples_left = <unsigned int> left_stats[count_idx]
                n_samples_right = n_samples - n_samples_left
                sum_hessian_left = _total_hessian(
                    left_stats, n_outputs, hessians_are_constant)
                sum_hessian_right = total_sum_hessians - sum_hessian_left

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

                gain = _vector_split_gain(
                    left_stats, sum_gradients, sum_hessians, n_outputs,
                    hessians_are_constant, loss_current_node,
                    self.l2_regularization)
                if gain > best_gain and gain > self.min_gain_to_split:
                    found_better_split = True
                    best_gain = gain
                    best_cat_infos_thresh = sorted_cat_idx
                    best_direction = direction
                    memcpy(best_left_stats, left_stats, n_stats * sizeof(Y_DTYPE_C))

        if found_better_split:
            split_info.gain = best_gain
            # bin_idx is unused for categorical splits: left_cat_bitset is used
            split_info.bin_idx = 0

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
        free(cat_values)
        free(cat_weights)
        free(direction_vec)
        free(new_direction_vec)


cdef inline Y_DTYPE_C _total_hessian(
        const Y_DTYPE_C * stats,
        unsigned int n_outputs,
        uint8_t hessians_are_constant) noexcept nogil:
    """Sum over the outputs of the hessians stored in a vector histogram bin.

    With constant hessians (equal to 1), this is n_outputs * count.
    """
    cdef:
        unsigned int k
        Y_DTYPE_C total = 0.
    if hessians_are_constant:
        return n_outputs * stats[n_outputs]
    for k in range(n_outputs):
        total += stats[n_outputs + k]
    return total


cdef inline Y_DTYPE_C _vector_split_gain(
        const Y_DTYPE_C * left_stats,
        const Y_DTYPE_C * sum_gradients,
        const Y_DTYPE_C * sum_hessians,
        unsigned int n_outputs,
        uint8_t hessians_are_constant,
        Y_DTYPE_C loss_current_node,
        Y_DTYPE_C l2_regularization) noexcept nogil:
    """Loss reduction of a split, summed over the outputs.

    With a single output, this performs exactly the same floating point
    operations as `_split_gain` without monotonic constraints.
    """
    cdef:
        unsigned int k
        Y_DTYPE_C gain = loss_current_node
        Y_DTYPE_C sum_gradient_left, sum_gradient_right
        Y_DTYPE_C sum_hessian_left, sum_hessian_right
        Y_DTYPE_C value_left, value_right

    for k in range(n_outputs):
        sum_gradient_left = left_stats[k]
        sum_gradient_right = sum_gradients[k] - sum_gradient_left
        if hessians_are_constant:
            sum_hessian_left = left_stats[n_outputs]
        else:
            sum_hessian_left = left_stats[n_outputs + k]
        sum_hessian_right = sum_hessians[k] - sum_hessian_left
        value_left = compute_node_value(sum_gradient_left, sum_hessian_left,
                                        -INFINITY, INFINITY, l2_regularization)
        value_right = compute_node_value(sum_gradient_right, sum_hessian_right,
                                         -INFINITY, INFINITY, l2_regularization)
        gain -= _loss_from_value(value_left, sum_gradient_left)
        gain -= _loss_from_value(value_right, sum_gradient_right)
    return gain


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
