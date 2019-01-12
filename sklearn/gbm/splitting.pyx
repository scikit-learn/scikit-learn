# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""This module contains njitted routines and data structures to:

- Find the best possible split of a node. For a given node, a split is
  characterized by a feature and a bin.
- Apply a split to a node, i.e. split the indices of the samples at the node
  into the newly created left and right childs.
"""
cimport cython
from cython.parallel import prange

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np

from .histogram cimport _build_histogram
from .histogram cimport _build_histogram_no_hessian
from .histogram cimport _build_histogram_root
from .histogram cimport _build_histogram_root_no_hessian
from .histogram cimport _subtract_histograms
from .histogram cimport NPY_X_BINNED_DTYPE
from .histogram cimport NPY_Y_DTYPE

from .types import HISTOGRAM_DTYPE

cdef struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count


cdef struct split_info_struct:
    float gain
    unsigned int feature_idx
    unsigned int bin_idx
    float gradient_left
    float gradient_right
    float hessian_left
    float hessian_right
    unsigned int n_samples_left
    unsigned int n_samples_right

@cython.freelist(100)
@cython.final
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

    def __init__(self, float gain=-1., unsigned int feature_idx=0, unsigned
                 int bin_idx=0, float gradient_left=0., float hessian_left=0.,
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


@cython.final
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
        NPY_X_BINNED_DTYPE [:, :] X_binned
        unsigned int n_features
        unsigned int max_bins
        unsigned int [:] n_bins_per_feature
        NPY_Y_DTYPE [:] gradients
        NPY_Y_DTYPE [:] hessians
        NPY_Y_DTYPE [:] ordered_gradients
        NPY_Y_DTYPE [:] ordered_hessians
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

    def __cinit__(self, NPY_X_BINNED_DTYPE [:, :] X_binned, unsigned int max_bins,
                 np.ndarray[np.uint32_t] n_bins_per_feature,
                 NPY_Y_DTYPE [:] gradients, NPY_Y_DTYPE [:] hessians, float l2_regularization,
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
        self.sum_gradients = np.sum(gradients)
        self.sum_hessians = np.sum(hessians)
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


def split_indices(
    SplittingContext context,
    SplitInfo split_info,
    unsigned int [:] sample_indices):
    cdef:
        unsigned int n_samples = sample_indices.shape[0]
        unsigned int i = 0
        unsigned int j = n_samples - 1
        unsigned char pivot = split_info.bin_idx
        unsigned int [:] view = sample_indices
        NPY_X_BINNED_DTYPE [:] binned_feature = context.X_binned.T[split_info.feature_idx]

    with nogil:
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

    return sample_indices[:i], sample_indices[i:], i


def find_node_split(
    SplittingContext context,
    unsigned int [:] sample_indices,
    hist_struct [:, :] histograms):
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
        int feature_idx
        int i
        unsigned int thread_idx
        unsigned int [:] starts
        unsigned int [:] ends
        unsigned int n_threads
        split_info_struct split_info
        split_info_struct * split_infos

    with nogil:
        n_samples = sample_indices.shape[0]

        # Populate ordered_gradients and ordered_hessians. (Already done for root)
        # Ordering the gradients and hessians helps to improve cache hit.
        if sample_indices.shape[0] != context.gradients.shape[0]:
            if context.constant_hessian:
                for i in prange(n_samples, schedule='static'):
                    context.ordered_gradients[i] = context.gradients[sample_indices[i]]
            else:
                for i in prange(n_samples, schedule='static'):
                    context.ordered_gradients[i] = context.gradients[sample_indices[i]]
                    context.ordered_hessians[i] = context.hessians[sample_indices[i]]

        context.sum_gradients = 0.
        for i in range(n_samples):
            context.sum_gradients += context.ordered_gradients[i]

        if context.constant_hessian:
            context.sum_hessians = context.constant_hessian_value * <float> (n_samples)
        else:
            context.sum_hessians = 0.
            for i in range(n_samples):
                context.sum_hessians += context.ordered_hessians[i]

        # TODO: this needs to be freed at some point
        split_infos = <split_info_struct *> malloc(context.n_features * sizeof(split_info_struct))
        for feature_idx in prange(context.n_features):
            split_info = _find_histogram_split(
                context, feature_idx, sample_indices, histograms[feature_idx])
            split_infos[feature_idx] = split_info

        split_info = _find_best_feature_to_split_helper(context, split_infos)

    return SplitInfo(
        split_info.gain,
        split_info.feature_idx,
        split_info.bin_idx,
        split_info.gradient_left,
        split_info.hessian_left,
        split_info.gradient_right,
        split_info.hessian_right,
        split_info.n_samples_left,
        split_info.n_samples_right,
    )


def find_node_split_subtraction(
    SplittingContext context,
    unsigned int [:] sample_indices,
    hist_struct [:, :] parent_histograms,
    hist_struct [:, :] sibling_histograms,
    hist_struct [:, :] histograms):
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
        int feature_idx
        unsigned int n_samples
        split_info_struct split_info
        split_info_struct * split_infos
        int i

    with nogil:
        n_samples = sample_indices.shape[0]

        # TODO: maybe change this computation... we could probably store sum_g/h in
        # the SplitInfo for a speed gain
        # Compute sum_hessians and sum_gradients.
        # We can pick any feature (here the first) in the histograms to
        # compute the gradients: they must be the same across all features
        # anyway, we have tests ensuring this. Maybe a more robust way would
        # be to compute an average but it's probably not worth it.
        context.sum_gradients = 0.
        for i in range(context.max_bins):
            context.sum_gradients += parent_histograms[0, i].sum_gradients - sibling_histograms[0, i].sum_gradients

        if context.constant_hessian:
            context.sum_hessians = \
                context.constant_hessian_value * float(n_samples)
        else:
            context.sum_hessians = 0.
            for i in range(context.max_bins):
                context.sum_hessians += parent_histograms[0, i].sum_hessians - sibling_histograms[0, i].sum_hessians

        # TODO: this needs to be freed at some point
        split_infos = <split_info_struct *> malloc(context.n_features * sizeof(split_info_struct))
        for feature_idx in prange(context.n_features):
            split_info = _find_histogram_split_subtraction(
                context, feature_idx, parent_histograms[feature_idx],
                sibling_histograms[feature_idx], histograms[feature_idx],
                n_samples)
            split_infos[feature_idx] = split_info

        split_info = _find_best_feature_to_split_helper(context, split_infos)

    return SplitInfo(
        split_info.gain,
        split_info.feature_idx,
        split_info.bin_idx,
        split_info.gradient_left,
        split_info.hessian_left,
        split_info.gradient_right,
        split_info.hessian_right,
        split_info.n_samples_left,
        split_info.n_samples_right,
    )


cdef split_info_struct _find_best_feature_to_split_helper(
    SplittingContext context,
    split_info_struct * split_infos) nogil:
    cdef:
        float gain
        float best_gain
        split_info_struct split_info
        split_info_struct best_split_info
        unsigned int feature_idx

    best_gain = -1.
    for feature_idx in range(context.n_features):
        split_info = split_infos[feature_idx]
        gain = split_info.gain
        if best_gain == -1 or gain > best_gain:
            best_gain = gain
            best_split_info = split_info
    return best_split_info

cdef split_info_struct _find_histogram_split(
    SplittingContext context,
    unsigned int feature_idx,
    unsigned int [:] sample_indices,
    hist_struct [:] histogram) nogil:
    """Compute the histogram for a given feature

    Returns the best SplitInfo among all the possible bins of the feature.
    """

    cdef:
        unsigned int n_samples = sample_indices.shape[0]
        NPY_X_BINNED_DTYPE [:] X_binned = context.X_binned.T[feature_idx]
        unsigned int root_node = X_binned.shape[0] == n_samples
        NPY_Y_DTYPE [:] ordered_gradients = context.ordered_gradients[:n_samples]
        NPY_Y_DTYPE [:] ordered_hessians = context.ordered_hessians[:n_samples]

    if root_node:
        if context.constant_hessian:
            _build_histogram_root_no_hessian(context.max_bins, X_binned,
                                             ordered_gradients, histogram)
        else:
            _build_histogram_root(context.max_bins, X_binned,
                                  ordered_gradients,
                                  ordered_hessians, histogram)
    else:
        if context.constant_hessian:
            _build_histogram_no_hessian(context.max_bins, sample_indices,
                                        X_binned, ordered_gradients, histogram)
        else:
            _build_histogram(context.max_bins, sample_indices, X_binned,
                             ordered_gradients, ordered_hessians, histogram)

    return _find_best_bin_to_split_helper(context, feature_idx, histogram,
                                          n_samples)

cdef split_info_struct _find_histogram_split_subtraction(
    SplittingContext context,
    unsigned int feature_idx,
    hist_struct [:] parent_histogram,
    hist_struct [:] sibling_histogram,
    hist_struct [:] histogram,
    unsigned int n_samples) nogil:
    """Compute the histogram by substraction of parent and sibling

    Uses the identity: hist(parent) = hist(left) + hist(right).
    Returns the best SplitInfo among all the possible bins of the feature.
    """

    _subtract_histograms(context.max_bins, parent_histogram,
                         sibling_histogram, histogram)

    return _find_best_bin_to_split_helper(context, feature_idx, histogram,
                                          n_samples)


cdef split_info_struct _find_best_bin_to_split_helper(
    SplittingContext context,
    unsigned int feature_idx,
    hist_struct [:] histogram,
    unsigned int n_samples) nogil:
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
        split_info_struct best_split

    best_split.gain = -1.
    gradient_left, hessian_left = 0., 0.
    n_samples_left = 0

    for bin_idx in range(context.n_bins_per_feature[feature_idx]):
        n_samples_left += histogram[bin_idx].count
        n_samples_right = n_samples_ - n_samples_left

        if context.constant_hessian:
            hessian_left += (<float> histogram[bin_idx].count
                             * context.constant_hessian_value)
        else:
            hessian_left += histogram[bin_idx].sum_hessians
        hessian_right = context.sum_hessians - hessian_left

        gradient_left += histogram[bin_idx].sum_gradients
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

    return best_split


cdef inline float _split_gain(
    float gradient_left,
    float hessian_left,
    float gradient_right,
    float hessian_right,
    float sum_gradients,
    float sum_hessians,
    float l2_regularization) nogil:
    """Loss reduction

    Compute the reduction in loss after taking a split, compared to keeping
    the node a leaf of the tree.

    See Equation 7 of:
    XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
    https://arxiv.org/abs/1603.02754
    """
    cdef:
        float gain
    gain = negative_loss(gradient_left, hessian_left, l2_regularization)
    gain += negative_loss(gradient_right, hessian_right, l2_regularization)
    gain -= negative_loss(sum_gradients, sum_hessians, l2_regularization)
    return gain

cdef inline float negative_loss(
    float gradient,
    float hessian,
    float l2_regularization) nogil:
    return (gradient * gradient) / (hessian + l2_regularization)

# Only used for tests... not sure how to do it
def _find_histogram_split_wrapper(
    SplittingContext context,
    unsigned int feature_idx,
    unsigned int [:] sample_indices,
    hist_struct [:] histogram):

    split_info = _find_histogram_split(context, feature_idx, sample_indices,
                                       histogram)
    return SplitInfo(
        split_info.gain,
        split_info.feature_idx,
        split_info.bin_idx,
        split_info.gradient_left,
        split_info.hessian_left,
        split_info.gradient_right,
        split_info.hessian_right,
        split_info.n_samples_left,
        split_info.n_samples_right,
    )
