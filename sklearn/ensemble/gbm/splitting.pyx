"""This module contains njitted routines and data structures to:

- Find the best possible split of a node. For a given node, a split is
  characterized by a feature and a bin.
- Apply a split to a node, i.e. split the indices of the samples at the node
  into the newly created left and right childs.
"""
cimport cython

import numpy as np
cimport numpy as np

from .histogram import _build_histogram
from .histogram import _subtract_histograms
from .histogram import _build_histogram_no_hessian
from .histogram import _build_histogram_root
from .histogram import _build_histogram_root_no_hessian
from .histogram import HISTOGRAM_DTYPE

cdef struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count


cdef get_threads_chunks(unsigned int total_size):
    """Get start and end indices of threads in an array of size total_size.

    The interval [0, total_size - 1] is divided into n_threads contiguous
    regions, and the starts and ends of each region are returned. Used to
    simulate a 'static' scheduling.
    """
    cdef:
        np.ndarray[np.uint32_t] sizes
        np.ndarray[np.uint32_t] starts
        np.ndarray[np.uint32_t] ends
        unsigned int n_threads

    n_threads = 1  # TODO: change this
    sizes = np.full(n_threads, total_size // n_threads, dtype=np.uint32)
    sizes[:total_size % n_threads] += 1
    starts = np.zeros(n_threads, dtype=np.uint32)
    starts[1:] = np.cumsum(sizes[:-1])
    ends = starts + sizes

    return starts, ends, n_threads

@cython.freelist(100)
cdef class SplitInfo:
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
    cdef public:
        float gain
        unsigned int feature_idx
        unsigned int bin_idx
        float gradient_left
        float gradient_right
        float hessian_left
        float hessian_right
        unsigned int n_samples_left
        unsigned int n_samples_right

    def __cinit__(self, float gain=-1., unsigned int feature_idx=0, unsigned
                  int bin_idx=0,
                 float gradient_left=0., float hessian_left=0.,
                 float gradient_right=0., float hessian_right=0.,
                 unsigned int n_samples_left=0, unsigned int n_samples_right=0):
        self.gain = gain
        self.feature_idx = feature_idx
        self.bin_idx = bin_idx
        self.gradient_left = gradient_left
        self.hessian_left = hessian_left
        self.gradient_right = gradient_right
        self.hessian_right = hessian_right
        self.n_samples_left = n_samples_left
        self.n_samples_right = n_samples_right


cdef class SplittingContext:
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
    cdef public:
        unsigned char [:, :] X_binned
        unsigned int n_features
        unsigned int max_bins
        unsigned int [:] n_bins_per_feature
        float [:] gradients
        float [:] hessians
        float [:] ordered_gradients
        float [:] ordered_hessians
        float sum_gradients
        float sum_hessians
        unsigned char constant_hessian
        float constant_hessian_value
        float l2_regularization
        float min_hessian_to_split
        unsigned int min_samples_leaf
        float min_gain_to_split

        unsigned int [:] partition
        unsigned int [:] left_indices_buffer
        unsigned int [:] right_indices_buffer

    def __cinit__(self, np.ndarray[np.uint8_t, ndim=2] X_binned, unsigned int max_bins,
                 np.ndarray[np.uint32_t] n_bins_per_feature,
                 np.ndarray [np.float32_t] gradients, np.ndarray[np.float32_t] hessians, float l2_regularization,
                 float min_hessian_to_split=1e-3, unsigned int min_samples_leaf=20,
                 float min_gain_to_split=0.):

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
        self.sum_gradients = gradients.sum()
        self.sum_hessians = hessians.sum()
        self.constant_hessian = hessians.shape[0] == 1
        self.l2_regularization = l2_regularization
        self.min_hessian_to_split = min_hessian_to_split
        self.min_samples_leaf = min_samples_leaf
        self.min_gain_to_split = min_gain_to_split
        if self.constant_hessian:
            self.constant_hessian_value = hessians[0]  # 1 scalar
        else:
            self.constant_hessian_value = 1.  # won't be used anyway

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


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def split_indices(SplittingContext context, SplitInfo split_info, unsigned int [:] sample_indices):
    cdef:
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int i = 0
        unsigned int j = n_samples - 1
        unsigned char pivot = split_info.bin_idx
        unsigned int [:] view = sample_indices
        unsigned char [:] binned_feature = context.X_binned.T[split_info.feature_idx]

    while i != j:
        # continue until we find an element that should be on right
        while binned_feature[view[i]] <= pivot and i < n_samples:
            i += 1
        # same, but now an element that should be on the left
        while binned_feature[view[j]] > pivot and j >= 0:
            j -= 1
        if i >= j:  # j can become smaller than j!
            break
        else:
            # swap
            view[i], view[j] = view[j], view[i]
            i += 1
            j -= 1

    return sample_indices[:i], sample_indices[i:]


def find_node_split(SplittingContext context, unsigned int [:] sample_indices):
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
    cdef:
        unsigned int n_samples
        hist_struct [:, :] view
        hist_struct [:] histogram
        unsigned int feature_idx
        unsigned int i
        unsigned int thread_idx
        SplittingContext ctx
        unsigned int [:] starts
        unsigned int [:] ends
        unsigned int n_threads
        SplitInfo split_info
        list split_infos

    ctx = context  # shorter name to avoid various line breaks
    n_samples = sample_indices.shape[0]

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
                    ctx.ordered_gradients[i] = ctx.gradients[sample_indices[i]]
        else:
            for thread_idx in range(n_threads):
                for i in range(starts[thread_idx], ends[thread_idx]):
                    ctx.ordered_gradients[i] = ctx.gradients[sample_indices[i]]
                    ctx.ordered_hessians[i] = ctx.hessians[sample_indices[i]]

    # ctx.sum_gradients = ctx.ordered_gradients[:n_samples].sum()
    ctx.sum_gradients = np.sum(ctx.ordered_gradients[:n_samples])
    if ctx.constant_hessian:
        ctx.sum_hessians = ctx.constant_hessian_value * np.float32(n_samples)
    else:
        # ctx.sum_hessians = ctx.ordered_hessians[:n_samples].sum()
        ctx.sum_hessians = np.sum(ctx.ordered_hessians[:n_samples])

    split_infos = [SplitInfo(-1., 0, 0, 0., 0., 0., 0., 0, 0)
                   for i in range(context.n_features)]
    histograms = np.empty(
        shape=(np.int64(context.n_features), np.int64(context.max_bins)),
        dtype=HISTOGRAM_DTYPE
    )
    view = histograms
    for feature_idx in range(context.n_features):
        split_info, histogram = _find_histogram_split(
            context, feature_idx, sample_indices)
        split_infos[feature_idx] = split_info
        view[feature_idx, :] = histogram

    split_info = _find_best_feature_to_split_helper(split_infos)
    return split_info, histograms


def find_node_split_subtraction(SplittingContext context, unsigned int [:]
                                sample_indices, np.ndarray parent_histograms,
                                np.ndarray sibling_histograms):
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

    cdef:
        hist_struct [:, :] view
        hist_struct [:] histogram
        unsigned int feature_idx
        unsigned int n_samples
        SplitInfo split_info
        list split_infos

    # We can pick any feature (here the first) in the histograms to
    # compute the gradients: they must be the same across all features
    # anyway, we have tests ensuring this. Maybe a more robust way would
    # be to compute an average but it's probably not worth it.
    context.sum_gradients = (parent_histograms[0]['sum_gradients'].sum() -
                             sibling_histograms[0]['sum_gradients'].sum())

    n_samples = sample_indices.shape[0]
    if context.constant_hessian:
        context.sum_hessians = \
            context.constant_hessian_value * float(n_samples)
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
    view = histograms
    for feature_idx in range(context.n_features):
        split_info, histogram = _find_histogram_split_subtraction(
            context, feature_idx, parent_histograms,
            sibling_histograms, n_samples)
        split_infos[feature_idx] = split_info
        view[feature_idx, :] = histogram

    split_info = _find_best_feature_to_split_helper(split_infos)
    return split_info, histograms


cdef SplitInfo _find_best_feature_to_split_helper(list split_infos):
    cdef:
        float gain
        float best_gain
        SplitInfo split_info
        SplitInfo best_split_info
        unsigned int i

    best_gain = -1.
    for i, split_info in enumerate(split_infos):
        gain = split_info.gain
        if best_gain == -1 or gain > best_gain:
            best_gain = gain
            best_split_info = split_info
    return best_split_info


cdef _find_histogram_split(SplittingContext context, unsigned int feature_idx,
                          unsigned int [:] sample_indices):
    """Compute the histogram for a given feature

    Returns the best SplitInfo among all the possible bins of the feature.
    """

    cdef:
        unsigned int n_samples = sample_indices.shape[0]
        unsigned char [:] X_binned = context.X_binned.T[feature_idx]
        unsigned int root_node = X_binned.shape[0] == n_samples
        float [:] ordered_gradients = context.ordered_gradients[:n_samples]
        float [:] ordered_hessians = context.ordered_hessians[:n_samples]
        np.ndarray histogram

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

cdef _find_histogram_split_subtraction(SplittingContext context, unsigned int feature_idx,
                                      np.ndarray parent_histograms, np.ndarray sibling_histograms,
                                      unsigned int n_samples):
    """Compute the histogram by substraction of parent and sibling

    Uses the identity: hist(parent) = hist(left) + hist(right).
    Returns the best SplitInfo among all the possible bins of the feature.
    """
    cdef:
        np.ndarray histogram

    histogram = _subtract_histograms(
        context.max_bins,
        parent_histograms[feature_idx], sibling_histograms[feature_idx])

    return _find_best_bin_to_split_helper(context, feature_idx, histogram,
                                          n_samples)


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef _find_best_bin_to_split_helper(SplittingContext context, unsigned int feature_idx,
                                    hist_struct [:] histogram, unsigned int
                                    n_samples):
    """Find best bin to split on, and return the corresponding SplitInfo.

    Splits that do not satisfy the splitting constraints (min_gain_to_split,
    etc.) are discarded here. If no split can satisfy the constraints, a
    SplitInfo with a gain of -1 is returned. If for a given node the best
    SplitInfo has a gain of -1, it is finalized into a leaf.
    """
    cdef:
        unsigned int bin_idx
        unsigned int n_samples_left
        unsigned int n_samples_right
        unsigned int n_samples_ = n_samples
        float hessian_left
        float hessian_right
        float gradient_left
        float gradient_right
        float gain
        SplitInfo best_split

        hist_struct [:] view = histogram

    best_split = SplitInfo.__new__(SplitInfo)
    gradient_left, hessian_left = 0., 0.
    n_samples_left = 0

    for bin_idx in range(context.n_bins_per_feature[feature_idx]):
        n_samples_left += view[bin_idx].count
        n_samples_right = n_samples_ - n_samples_left

        if context.constant_hessian:
            hessian_left += (<float> view[bin_idx].count
                             * context.constant_hessian_value)
        else:
            hessian_left += view[bin_idx].sum_hessians
        hessian_right = context.sum_hessians - hessian_left

        gradient_left += view[bin_idx].sum_gradients
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
            best_split.gradient_right = gradient_right
            best_split.hessian_left = hessian_left
            best_split.hessian_right = hessian_right
            best_split.n_samples_left = n_samples_left
            best_split.n_samples_right = n_samples_right
            """
            best_split = SplitInfo(
                gain,
                feature_idx,
                bin_idx,
                gradient_left,
                gradient_right,
                hessian_left,
                hessian_right,
                n_samples_left,
                n_samples_right,
            )
            """

    return best_split, histogram


cdef inline float _split_gain(float gradient_left, float hessian_left, float gradient_right,
                 float hessian_right, float sum_gradients, float
                 sum_hessians, float l2_regularization) nogil:
    """Loss reduction

    Compute the reduction in loss after taking a split, compared to keeping
    the node a leaf of the tree.

    See Equation 7 of:
    XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
    https://arxiv.org/abs/1603.02754
    """
    cdef float gain
    gain = negative_loss(gradient_left, hessian_left, l2_regularization)
    gain += negative_loss(gradient_right, hessian_right, l2_regularization)
    gain -= negative_loss(sum_gradients, sum_hessians, l2_regularization)
    return gain

@cython.cdivision(True)
cdef inline float negative_loss(float gradient, float hessian, float
l2_regularization) nogil:
    return (gradient * gradient) / (hessian + l2_regularization)
