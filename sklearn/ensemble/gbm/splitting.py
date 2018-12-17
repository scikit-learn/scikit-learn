"""This module contains njitted routines and data structures to:

- Find the best possible split of a node. For a given node, a split is
  characterized by a feature and a bin.
- Apply a split to a node, i.e. split the indices of the samples at the node
  into the newly created left and right childs.
"""
import numpy as np

from .histogram import _build_histogram
from .histogram import _subtract_histograms
from .histogram import _build_histogram_no_hessian
from .histogram import _build_histogram_root
from .histogram import _build_histogram_root_no_hessian
from .histogram import HISTOGRAM_DTYPE
from .utils import get_threads_chunks


class SplitInfo:
    """Pure data class to store information about a potential split.

    Parameters
    ----------
    gain : float32
        The gain of the split
    feature_idx : int
        The index of the feature to be split
    bin_idx : int
        The index of the bin on which the split is made
    gradient_left : float32
        The sum of the gradients of all the samples in the left child
    hessian_left : float32
        The sum of the hessians of all the samples in the left child
    gradient_right : float32
        The sum of the gradients of all the samples in the right child
    hessian_right : float32
        The sum of the hessians of all the samples in the right child
    n_samples_left : int
        The number of samples in the left child
    n_samples_right : int
        The number of samples in the right child
    """
    def __init__(self, gain=-1., feature_idx=0, bin_idx=0,
                 gradient_left=0., hessian_left=0.,
                 gradient_right=0., hessian_right=0.,
                 n_samples_left=0, n_samples_right=0):
        self.gain = gain
        self.feature_idx = feature_idx
        self.bin_idx = bin_idx
        self.gradient_left = gradient_left
        self.hessian_left = hessian_left
        self.gradient_right = gradient_right
        self.hessian_right = hessian_right
        self.n_samples_left = n_samples_left
        self.n_samples_right = n_samples_right


class SplittingContext:
    """Pure data class defining a splitting context.

    Ideally it would also have methods but numba does not support annotating
    jitclasses (so we can't use parallel=True). This structure is
    instanciated in the grower and stores all the required information to
    compute the SplitInfo and histograms of each node.

    Parameters
    ----------
    X_binned : array of int
        The binned input samples. Must be Fortran-aligned.
    max_bins : int, optional(default=256)
        The maximum number of bins. Used to define the shape of the
        histograms.
    n_bins_per_feature : array-like of int
        The actual number of bins needed for each feature, which is lower or
        equal to max_bins.
    gradients : array-like, shape=(n_samples,)
        The gradients of each training sample. Those are the gradients of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    hessians : array-like, shape=(n_samples,)
        The hessians of each training sample. Those are the hessians of the
        loss w.r.t the predictions, evaluated at iteration i - 1.
    l2_regularization : float
        The L2 regularization parameter.
    min_hessian_to_split : float
        The minimum sum of hessians needed in each node. Splits that result in
        at least one child having a sum of hessians less than
        min_hessian_to_split are discarded.
    min_samples_leaf : int
        The minimum number of samples per leaf.
    min_gain_to_split : float, optional(default=0.)
        The minimum gain needed to split a node. Splits with lower gain will
        be ignored.
    """
    def __init__(self, X_binned, max_bins, n_bins_per_feature,
                 gradients, hessians, l2_regularization,
                 min_hessian_to_split=1e-3, min_samples_leaf=20,
                 min_gain_to_split=0.):

        self.X_binned = X_binned
        self.n_features = X_binned.shape[1]
        # Note: all histograms will have <max_bins> bins, but some of the
        # last bins may be unused if n_bins_per_feature[f] < max_bins
        self.max_bins = max_bins
        self.n_bins_per_feature = n_bins_per_feature
        self.gradients = gradients
        self.hessians = hessians
        # for root node, gradients and hessians are already ordered
        self.ordered_gradients = gradients.copy()
        self.ordered_hessians = hessians.copy()
        self.sum_gradients = self.gradients.sum()
        self.sum_hessians = self.hessians.sum()
        self.constant_hessian = hessians.shape[0] == 1
        self.l2_regularization = l2_regularization
        self.min_hessian_to_split = min_hessian_to_split
        self.min_samples_leaf = min_samples_leaf
        self.min_gain_to_split = min_gain_to_split
        if self.constant_hessian:
            self.constant_hessian_value = self.hessians[0]  # 1 scalar
        else:
            self.constant_hessian_value = np.float32(1.)  # won't be used anyway

        # The partition array maps each sample index into the leaves of the
        # tree (a leaf in this context is a node that isn't splitted yet, not
        # necessarily a 'finalized' leaf). Initially, the root contains all
        # the indices, e.g.:
        # partition = [abcdefghijkl]
        # After a call to split_indices, it may look e.g. like this:
        # partition = [cef|abdghijkl]
        # we have 2 leaves, the left one is at position 0 and the second one at
        # position 3. The order of the samples is irrelevant.
        self.partition = np.arange(0, X_binned.shape[0], 1, np.uint32)
        # buffers used in split_indices to support parallel splitting.
        self.left_indices_buffer = np.empty_like(self.partition)
        self.right_indices_buffer = np.empty_like(self.partition)


def split_indices(context, split_info, sample_indices):
    """Split samples into left and right arrays.

    Parameters
    ----------
    context : SplittingContext
        The splitting context
    split_ingo : SplitInfo
        The SplitInfo of the node to split
    sample_indices : array of int
        The indices of the samples at the node to split. This is a view on
        context.partition, and it is modified inplace by placing the indices
        of the left child at the beginning, and the indices of the right child
        at the end.

    Returns
    -------
    left_indices : array of int
        The indices of the samples in the left child. This is a view on
        context.partition.
    right_indices : array of int
        The indices of the samples in the right child. This is a view on
        context.partition.
    """
    # This is a multi-threaded implementation inspired by lightgbm.
    # Here is a quick break down. Let's suppose we want to split a node with
    # 24 samples named from a to x. context.partition looks like this (the *
    # are indices in other leaves that we don't care about):
    # partition = [*************abcdefghijklmnopqrstuvwx****************]
    #                           ^                       ^
    #                     node_position     node_position + node.n_samples

    # Ultimately, we want to reorder the samples inside the boundaries of the
    # leaf (which becomes a node) to now represent the samples in its left and
    # right child. For example:
    # partition = [*************abefilmnopqrtuxcdghjksvw*****************]
    #                           ^              ^
    #                   left_child_pos     right_child_pos
    # Note that left_child_pos always takes the value of node_position, and
    # right_child_pos = left_child_pos + left_child.n_samples. The order of
    # the samples inside a leaf is irrelevant.

    # 1. samples_indices is a view on this region a..x. We conceptually
    #    divide it into n_threads regions. Each thread will be responsible for
    #    its own region. Here is an example with 4 threads:
    #    samples_indices = [abcdef|ghijkl|mnopqr|stuvwx]
    # 2. Each thread processes 6 = 24 // 4 entries and maps them into
    #    left_indices_buffer or right_indices_buffer. For example, we could
    #    have the following mapping ('.' denotes an undefined entry):
    #    - left_indices_buffer =  [abef..|il....|mnopqr|tux...]
    #    - right_indices_buffer = [cd....|ghjk..|......|svw...]
    # 3. We keep track of the start positions of the regions (the '|') in
    #    ``offset_in_buffers`` as well as the size of each region. We also keep
    #    track of the number of samples put into the left/right child by each
    #    thread. Concretely:
    #    - left_counts =  [4, 2, 6, 3]
    #    - right_counts = [2, 4, 0, 3]
    # 4. Finally, we put left/right_indices_buffer back into the
    #    samples_indices, without any undefined entries and the partition looks
    #    as expected
    #    partition = [*************abefilmnopqrtuxcdghjksvw*****************]

    # Note: We here show left/right_indices_buffer as being the same size as
    # sample_indices for simplicity, but in reality they are of the same size
    # as partition.

    X_binned = context.X_binned.T[split_info.feature_idx]

    n_threads = 4  # TODO: change this
    n_samples = sample_indices.shape[0]

    # Note: we could probably allocate all the arrays of size n_threads in the
    # splitting context as well, but gains are probably going to be minimal
    sizes = np.full(n_threads, n_samples // n_threads, dtype=np.int32)
    if n_samples % n_threads > 0:
        # array[:0] will cause a bug in numba 0.41 so we need the if. Remove
        # once issue numba 3554 is fixed.
        sizes[:n_samples % n_threads] += 1
    offset_in_buffers = np.zeros(n_threads, dtype=np.int32)
    offset_in_buffers[1:] = np.cumsum(sizes[:-1])

    left_counts = np.empty(n_threads, dtype=np.int32)
    right_counts = np.empty(n_threads, dtype=np.int32)

    # Need to declare local variables, else they're not updated :/
    # (see numba issue 3459)
    left_indices_buffer = context.left_indices_buffer
    right_indices_buffer = context.right_indices_buffer

    # map indices from samples_indices to left/right_indices_buffer
    for thread_idx in range(n_threads):
        left_count = 0
        right_count = 0

        start = offset_in_buffers[thread_idx]
        stop = start + sizes[thread_idx]
        for i in range(start, stop):
            sample_idx = sample_indices[i]
            if X_binned[sample_idx] <= split_info.bin_idx:
                left_indices_buffer[start + left_count] = sample_idx
                left_count += 1
            else:
                right_indices_buffer[start + right_count] = sample_idx
                right_count += 1

        left_counts[thread_idx] = left_count
        right_counts[thread_idx] = right_count

    # position of right child = just after the left child
    right_child_position = left_counts.sum()

    # offset of each thread in samples_indices for left and right child, i.e.
    # where each thread will start to write.
    left_offset = np.zeros(n_threads, dtype=np.int32)
    left_offset[1:] = np.cumsum(left_counts[:-1])
    right_offset = np.full(n_threads, right_child_position, dtype=np.int32)
    right_offset[1:] += np.cumsum(right_counts[:-1])

    # map indices in left/right_indices_buffer back into samples_indices. This
    # also updates context.partition since samples_indice is a view.
    for thread_idx in range(n_threads):

        for i in range(left_counts[thread_idx]):
            sample_indices[left_offset[thread_idx] + i] = \
                left_indices_buffer[offset_in_buffers[thread_idx] + i]
        for i in range(right_counts[thread_idx]):
            sample_indices[right_offset[thread_idx] + i] = \
                right_indices_buffer[offset_in_buffers[thread_idx] + i]

    return (sample_indices[:right_child_position],
            sample_indices[right_child_position:])


def find_node_split(context, sample_indices):
    """For each feature, find the best bin to split on at a given node.

    Returns the best split info among all features, and the histograms of
    all the features. The histograms are computed by scanning the whole
    data.

    Parameters
    ----------
    context : SplittingContext
        The splitting context
    sample_indices : array of int
        The indices of the samples at the node to split.

    Returns
    -------
    best_split_info : SplitInfo
        The info about the best possible split among all features.
    histograms : array of HISTOGRAM_DTYPE, shape=(n_features, max_bins)
        The histograms of each feature. A histogram is an array of
        HISTOGRAM_DTYPE of size ``max_bins`` (only
        ``n_bins_per_features[feature]`` entries are relevant).
    """

    ctx = context  # shorter name to avoid various line breaks
    n_samples = sample_indices.shape[0]

    # Need to declare local variables, else they're not updated
    # (see numba issue 3459)
    ordered_gradients = ctx.ordered_gradients
    ordered_hessians = ctx.ordered_hessians

    # Populate ordered_gradients and ordered_hessians. (Already done for root)
    # Ordering the gradients and hessians helps to improve cache hit.
    # This is a parallelized version of the following vanilla code:
    # for i range(n_samples):
    #     ctx.ordered_gradients[i] = ctx.gradients[samples_indices[i]]
    if sample_indices.shape[0] != ctx.gradients.shape[0]:
        starts, ends, n_threads = get_threads_chunks(n_samples)
        if ctx.constant_hessian:
            for thread_idx in range(n_threads):
                for i in range(starts[thread_idx], ends[thread_idx]):
                    ordered_gradients[i] = ctx.gradients[sample_indices[i]]
        else:
            for thread_idx in range(n_threads):
                for i in range(starts[thread_idx], ends[thread_idx]):
                    ordered_gradients[i] = ctx.gradients[sample_indices[i]]
                    ordered_hessians[i] = ctx.hessians[sample_indices[i]]

    ctx.sum_gradients = ctx.ordered_gradients[:n_samples].sum()
    if ctx.constant_hessian:
        ctx.sum_hessians = ctx.constant_hessian_value * float32(n_samples)
    else:
        ctx.sum_hessians = ctx.ordered_hessians[:n_samples].sum()

    # Pre-allocate the results datastructure to be able to use prange:
    # numba jitclass do not seem to properly support default values for kwargs.
    split_infos = [SplitInfo(-1., 0, 0, 0., 0., 0., 0., 0, 0)
                   for i in range(context.n_features)]
    histograms = np.empty(
        shape=(np.int64(context.n_features), np.int64(context.max_bins)),
        dtype=HISTOGRAM_DTYPE
    )
    for feature_idx in range(context.n_features):
        split_info, histogram = _find_histogram_split(
            context, feature_idx, sample_indices)
        split_infos[feature_idx] = split_info
        histograms[feature_idx, :] = histogram

    split_info = _find_best_feature_to_split_helper(split_infos)
    return split_info, histograms


def find_node_split_subtraction(context, sample_indices, parent_histograms,
                                sibling_histograms):
    """For each feature, find the best bin to split on at a given node.

    Returns the best split info among all features, and the histograms of
    all the features.

    This does the same job as ``find_node_split()`` but uses the histograms
    of the parent and sibling of the node to split. This allows to use the
    identity: ``histogram(parent) = histogram(node) - histogram(sibling)``,
    which is significantly faster than computing the histograms from data.

    Returns the best SplitInfo among all features, along with all the feature
    histograms that can be latter used to compute the sibling or children
    histograms by substraction.

    Parameters
    ----------
    context : SplittingContext
        The splitting context
    sample_indices : array of int
        The indices of the samples at the node to split.
    parent_histograms : array of HISTOGRAM_DTYPE of shape(n_features, max_bins)
        The histograms of the parent
    sibling_histograms : array of HISTOGRAM_DTYPE of \
        shape(n_features, max_bins)
        The histograms of the sibling

    Returns
    -------
    best_split_info : SplitInfo
        The info about the best possible split among all features.
    histograms : array of HISTOGRAM_DTYPE, shape=(n_features, max_bins)
        The histograms of each feature. A histogram is an array of
        HISTOGRAM_DTYPE of size ``max_bins`` (only
        ``n_bins_per_features[feature]`` entries are relevant).
    """

    # We can pick any feature (here the first) in the histograms to
    # compute the gradients: they must be the same across all features
    # anyway, we have tests ensuring this. Maybe a more robust way would
    # be to compute an average but it's probably not worth it.
    context.sum_gradients = (parent_histograms[0]['sum_gradients'].sum() -
                             sibling_histograms[0]['sum_gradients'].sum())

    n_samples = sample_indices.shape[0]
    if context.constant_hessian:
        context.sum_hessians = \
            context.constant_hessian_value * np.float32(n_samples)
    else:
        context.sum_hessians = (parent_histograms[0]['sum_hessians'].sum() -
                                sibling_histograms[0]['sum_hessians'].sum())

    # Pre-allocate the results datastructure to be able to use prange
    split_infos = [SplitInfo(-1., 0, 0, 0., 0., 0., 0., 0, 0)
                   for i in range(context.n_features)]
    histograms = np.empty(
        shape=(np.int64(context.n_features), np.int64(context.max_bins)),
        dtype=HISTOGRAM_DTYPE
    )
    for feature_idx in range(context.n_features):
        split_info, histogram = _find_histogram_split_subtraction(
            context, feature_idx, parent_histograms,
            sibling_histograms, n_samples)
        split_infos[feature_idx] = split_info
        histograms[feature_idx, :] = histogram

    split_info = _find_best_feature_to_split_helper(split_infos)
    return split_info, histograms


def _find_best_feature_to_split_helper(split_infos):
    best_gain = None
    for i, split_info in enumerate(split_infos):
        gain = split_info.gain
        if best_gain is None or gain > best_gain:
            best_gain = gain
            best_split_info = split_info
    return best_split_info


def _find_histogram_split(context, feature_idx, sample_indices):
    """Compute the histogram for a given feature

    Returns the best SplitInfo among all the possible bins of the feature.
    """
    n_samples = sample_indices.shape[0]
    X_binned = context.X_binned.T[feature_idx]

    root_node = X_binned.shape[0] == n_samples
    ordered_gradients = context.ordered_gradients[:n_samples]
    ordered_hessians = context.ordered_hessians[:n_samples]

    if root_node:
        if context.constant_hessian:
            histogram = _build_histogram_root_no_hessian(
                context.max_bins, X_binned, ordered_gradients)
        else:
            histogram = _build_histogram_root(
                context.max_bins, X_binned, ordered_gradients,
                context.ordered_hessians)
    else:
        if context.constant_hessian:
            histogram = _build_histogram_no_hessian(
                context.max_bins, sample_indices, X_binned,
                ordered_gradients)
        else:
            histogram = _build_histogram(
                context.max_bins, sample_indices, X_binned,
                ordered_gradients, ordered_hessians)

    return _find_best_bin_to_split_helper(context, feature_idx, histogram,
                                          n_samples)


def _find_histogram_split_subtraction(context, feature_idx,
                                      parent_histograms, sibling_histograms,
                                      n_samples):
    """Compute the histogram by substraction of parent and sibling

    Uses the identity: hist(parent) = hist(left) + hist(right).
    Returns the best SplitInfo among all the possible bins of the feature.
    """
    histogram = _subtract_histograms(
        context.max_bins,
        parent_histograms[feature_idx], sibling_histograms[feature_idx])

    return _find_best_bin_to_split_helper(context, feature_idx, histogram,
                                          n_samples)


def _find_best_bin_to_split_helper(context, feature_idx, histogram, n_samples):
    """Find best bin to split on, and return the corresponding SplitInfo.

    Splits that do not satisfy the splitting constraints (min_gain_to_split,
    etc.) are discarded here. If no split can satisfy the constraints, a
    SplitInfo with a gain of -1 is returned. If for a given node the best
    SplitInfo has a gain of -1, it is finalized into a leaf.
    """
    # Allocate the structure for the best split information. It can be
    # returned as such (with a negative gain) if the min_hessian_to_split
    # condition is not satisfied. Such invalid splits are later discarded by
    # the TreeGrower.
    best_split = SplitInfo(-1., 0, 0, 0., 0., 0., 0., 0, 0)
    gradient_left, hessian_left = 0., 0.
    n_samples_left = 0

    for bin_idx in range(context.n_bins_per_feature[feature_idx]):
        n_samples_left += histogram[bin_idx]['count']
        n_samples_right = n_samples - n_samples_left

        if context.constant_hessian:
            hessian_left += (histogram[bin_idx]['count']
                             * context.constant_hessian_value)
        else:
            hessian_left += histogram[bin_idx]['sum_hessians']
        hessian_right = context.sum_hessians - hessian_left

        gradient_left += histogram[bin_idx]['sum_gradients']
        gradient_right = context.sum_gradients - gradient_left

        if n_samples_left < context.min_samples_leaf:
            continue
        if n_samples_right < context.min_samples_leaf:
            # won't get any better
            break

        if hessian_left < context.min_hessian_to_split:
            continue
        if hessian_right < context.min_hessian_to_split:
            # won't get any better (hessians are > 0 since loss is convex)
            break

        gain = _split_gain(gradient_left, hessian_left,
                           gradient_right, hessian_right,
                           context.sum_gradients, context.sum_hessians,
                           context.l2_regularization)

        if gain > best_split.gain and gain > context.min_gain_to_split:
            best_split.gain = gain
            best_split.feature_idx = feature_idx
            best_split.bin_idx = bin_idx
            best_split.gradient_left = gradient_left
            best_split.hessian_left = hessian_left
            best_split.n_samples_left = n_samples_left
            best_split.gradient_right = gradient_right
            best_split.hessian_right = hessian_right
            best_split.n_samples_right = n_samples_right

    return best_split, histogram


def _split_gain(gradient_left, hessian_left, gradient_right, hessian_right,
                sum_gradients, sum_hessians, l2_regularization):
    """Loss reduction

    Compute the reduction in loss after taking a split, compared to keeping
    the node a leaf of the tree.

    See Equation 7 of:
    XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
    https://arxiv.org/abs/1603.02754
    """
    def negative_loss(gradient, hessian):
        return (gradient ** 2) / (hessian + l2_regularization)

    gain = negative_loss(gradient_left, hessian_left)
    gain += negative_loss(gradient_right, hessian_right)
    gain -= negative_loss(sum_gradients, sum_hessians)
    return gain
