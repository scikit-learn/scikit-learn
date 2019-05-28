# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""This module contains routines and data structures to:

- Find the best possible split of a node. For a given node, a split is
  characterized by a feature and a bin.
- Apply a split to a node, i.e. split the indices of the samples at the node
  into the newly created left and right childs.
"""
# Author: Nicolas Hug

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np
IF SKLEARN_OPENMP_SUPPORTED:
    from openmp cimport omp_get_max_threads
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

from .types cimport X_BINNED_DTYPE_C
from .types cimport Y_DTYPE_C
from .types cimport hist_struct
from .types import HISTOGRAM_DTYPE


cdef struct split_info_struct:
    # Same as the SplitInfo class, but we need a C struct to use it in the
    # nogil sections and to use in arrays.
    Y_DTYPE_C gain
    int feature_idx
    unsigned int bin_idx
    Y_DTYPE_C sum_gradient_left
    Y_DTYPE_C sum_gradient_right
    Y_DTYPE_C sum_hessian_left
    Y_DTYPE_C sum_hessian_right
    unsigned int n_samples_left
    unsigned int n_samples_right


class SplitInfo:
    """Pure data class to store information about a potential split.

    Parameters
    ----------
    gain : float
        The gain of the split.
    feature_idx : int
        The index of the feature to be split.
    bin_idx : int
        The index of the bin on which the split is made.
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
    """
    def __init__(self, gain, feature_idx, bin_idx, sum_gradient_left,
                 sum_hessian_left, sum_gradient_right, sum_hessian_right,
                 n_samples_left, n_samples_right):
        self.gain = gain
        self.feature_idx = feature_idx
        self.bin_idx = bin_idx
        self.sum_gradient_left = sum_gradient_left
        self.sum_hessian_left = sum_hessian_left
        self.sum_gradient_right = sum_gradient_right
        self.sum_hessian_right = sum_hessian_right
        self.n_samples_left = n_samples_left
        self.n_samples_right = n_samples_right


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
    actual_n_bins : ndarray, shape (n_features,)
        The actual number of bins needed for each feature, which is lower or
        equal to max_bins.
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
    """
    cdef public:
        const X_BINNED_DTYPE_C [::1, :] X_binned
        unsigned int n_features
        unsigned int [::1] actual_n_bins
        unsigned char hessians_are_constant
        Y_DTYPE_C l2_regularization
        Y_DTYPE_C min_hessian_to_split
        unsigned int min_samples_leaf
        Y_DTYPE_C min_gain_to_split

        unsigned int [::1] partition
        unsigned int [::1] left_indices_buffer
        unsigned int [::1] right_indices_buffer

    def __init__(self, const X_BINNED_DTYPE_C [::1, :] X_binned,
                 np.ndarray[np.uint32_t] actual_n_bins,
                 Y_DTYPE_C l2_regularization, Y_DTYPE_C
                 min_hessian_to_split=1e-3, unsigned int
                 min_samples_leaf=20, Y_DTYPE_C min_gain_to_split=0.,
                 unsigned char hessians_are_constant=False):

        self.X_binned = X_binned
        self.n_features = X_binned.shape[1]
        self.actual_n_bins = actual_n_bins
        self.l2_regularization = l2_regularization
        self.min_hessian_to_split = min_hessian_to_split
        self.min_samples_leaf = min_samples_leaf
        self.min_gain_to_split = min_gain_to_split
        self.hessians_are_constant = hessians_are_constant

        # The partition array maps each sample index into the leaves of the
        # tree (a leaf in this context is a node that isn't splitted yet, not
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
            int feature_idx = split_info.feature_idx
            const X_BINNED_DTYPE_C [::1] X_binned = \
                self.X_binned[:, feature_idx]
            unsigned int [::1] left_indices_buffer = self.left_indices_buffer
            unsigned int [::1] right_indices_buffer = self.right_indices_buffer

            IF SKLEARN_OPENMP_SUPPORTED:
                int n_threads = omp_get_max_threads()
            ELSE:
                int n_threads = 1

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
            int [:] left_offset = np.zeros(n_threads, dtype=np.int32)
            int [:] right_offset = np.zeros(n_threads, dtype=np.int32)

        with nogil:
            for thread_idx in range(n_samples % n_threads):
                sizes[thread_idx] += 1

            for thread_idx in range(1, n_threads):
                offset_in_buffers[thread_idx] = \
                    offset_in_buffers[thread_idx - 1] + sizes[thread_idx - 1]

            # map indices from sample_indices to left/right_indices_buffer
            for thread_idx in prange(n_threads, schedule='static',
                                     chunksize=1):
                left_count = 0
                right_count = 0

                start = offset_in_buffers[thread_idx]
                stop = start + sizes[thread_idx]
                for i in range(start, stop):
                    sample_idx = sample_indices[i]
                    if X_binned[sample_idx] <= bin_idx:
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
                                     chunksize=1):
                memcpy(
                    &sample_indices[left_offset[thread_idx]],
                    &left_indices_buffer[offset_in_buffers[thread_idx]],
                    sizeof(unsigned int) * left_counts[thread_idx]
                )
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
            const unsigned int [::1] sample_indices,  # IN
            hist_struct [:, ::1] histograms,  # IN
            const Y_DTYPE_C sum_gradients,
            const Y_DTYPE_C sum_hessians):
        """For each feature, find the best bin to split on at a given node.

        Return the best split info among all features.

        Parameters
        ----------
        sample_indices : ndarray of unsigned int, shape (n_samples_at_node,)
            The indices of the samples at the node to split.
        histograms : ndarray of HISTOGRAM_DTYPE of \
                shape (n_features, max_bins)
            The histograms of the current node.
        sum_gradients : float
            The sum of the gradients for each sample at the node.
        sum_hessians : float
            The sum of the hessians for each sample at the node.

        Returns
        -------
        best_split_info : SplitInfo
            The info about the best possible split among all features.
        """
        cdef:
            int n_samples
            int feature_idx
            int best_feature_idx
            int n_features = self.n_features
            split_info_struct split_info
            split_info_struct * split_infos

        with nogil:
            n_samples = sample_indices.shape[0]

            split_infos = <split_info_struct *> malloc(
                self.n_features * sizeof(split_info_struct))

            for feature_idx in prange(n_features, schedule='static'):
                # For each feature, find best bin to split on
                split_info = self._find_best_bin_to_split_helper(
                    feature_idx, histograms, n_samples,
                    sum_gradients, sum_hessians)
                split_infos[feature_idx] = split_info

            # then compute best possible split among all features
            best_feature_idx = self._find_best_feature_to_split_helper(
                split_infos)
            split_info = split_infos[best_feature_idx]

        out = SplitInfo(
            split_info.gain,
            split_info.feature_idx,
            split_info.bin_idx,
            split_info.sum_gradient_left,
            split_info.sum_hessian_left,
            split_info.sum_gradient_right,
            split_info.sum_hessian_right,
            split_info.n_samples_left,
            split_info.n_samples_right,
        )
        free(split_infos)
        return out

    cdef int _find_best_feature_to_split_helper(
            self,
            split_info_struct * split_infos) nogil:  # IN
        """Returns the best feature among those in splits_infos."""
        cdef:
            int feature_idx
            int best_feature_idx = 0

        for feature_idx in range(1, self.n_features):
            if (split_infos[feature_idx].gain >
                    split_infos[best_feature_idx].gain):
                best_feature_idx = feature_idx
        return best_feature_idx

    cdef split_info_struct _find_best_bin_to_split_helper(
            self,
            unsigned int feature_idx,
            const hist_struct [:, ::1] histograms,  # IN
            unsigned int n_samples,
            Y_DTYPE_C sum_gradients,
            Y_DTYPE_C sum_hessians) nogil:
        """Find best bin to split on for a given feature.

        Splits that do not satisfy the splitting constraints
        (min_gain_to_split, etc.) are discarded here. If no split can
        satisfy the constraints, a SplitInfo with a gain of -1 is returned.
        If for a given node the best SplitInfo has a gain of -1, it is
        finalized into a leaf in the grower.
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
            Y_DTYPE_C negative_loss_current_node
            Y_DTYPE_C gain
            split_info_struct best_split

        best_split.gain = -1.
        sum_gradient_left, sum_hessian_left = 0., 0.
        n_samples_left = 0
        negative_loss_current_node = negative_loss(sum_gradients,
            sum_hessians, self.l2_regularization)

        for bin_idx in range(self.actual_n_bins[feature_idx]):
            n_samples_left += histograms[feature_idx, bin_idx].count
            n_samples_right = n_samples_ - n_samples_left

            if self.hessians_are_constant:
                sum_hessian_left += histograms[feature_idx, bin_idx].count
            else:
                sum_hessian_left += \
                    histograms[feature_idx, bin_idx].sum_hessians
            sum_hessian_right = sum_hessians - sum_hessian_left

            sum_gradient_left += histograms[feature_idx, bin_idx].sum_gradients
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
                               negative_loss_current_node,
                               self.l2_regularization)

            if gain > best_split.gain and gain > self.min_gain_to_split:
                best_split.gain = gain
                best_split.feature_idx = feature_idx
                best_split.bin_idx = bin_idx
                best_split.sum_gradient_left = sum_gradient_left
                best_split.sum_gradient_right = sum_gradient_right
                best_split.sum_hessian_left = sum_hessian_left
                best_split.sum_hessian_right = sum_hessian_right
                best_split.n_samples_left = n_samples_left
                best_split.n_samples_right = n_samples_right

        return best_split


cdef inline Y_DTYPE_C _split_gain(
        Y_DTYPE_C sum_gradient_left,
        Y_DTYPE_C sum_hessian_left,
        Y_DTYPE_C sum_gradient_right,
        Y_DTYPE_C sum_hessian_right,
        Y_DTYPE_C negative_loss_current_node,
        Y_DTYPE_C l2_regularization) nogil:
    """Loss reduction

    Compute the reduction in loss after taking a split, compared to keeping
    the node a leaf of the tree.

    See Equation 7 of:
    XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
    https://arxiv.org/abs/1603.02754
    """
    cdef:
        Y_DTYPE_C gain
    gain = negative_loss(sum_gradient_left, sum_hessian_left,
                         l2_regularization)
    gain += negative_loss(sum_gradient_right, sum_hessian_right,
                          l2_regularization)
    gain -= negative_loss_current_node
    return gain

cdef inline Y_DTYPE_C negative_loss(
        Y_DTYPE_C gradient,
        Y_DTYPE_C hessian,
        Y_DTYPE_C l2_regularization) nogil:
    return (gradient * gradient) / (hessian + l2_regularization)
