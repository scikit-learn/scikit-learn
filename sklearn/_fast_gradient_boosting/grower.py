"""
This module contains the TreeGrower class.

TreeGrowee builds a regression tree fitting a Newton-Raphson step, based on
the gradients and hessians of the training data.
"""
from heapq import heappush, heappop
import numpy as np
from time import time

from .splitting import Splitter
from .predictor import TreePredictor, PREDICTOR_RECORD_DTYPE
from .utils import sum_parallel
from .types import HISTOGRAM_DTYPE


class TreeNode:
    """Tree Node class used in TreeGrower.

    This isn't used for prediction purposes, only for training (see
    TreePredictor).

    Parameters
    ----------
    depth : int
        The depth of the node, i.e. its distance from the root
    sample_indices : array of int
        The indices of the samples at the node
    sum_gradients : float
        The sum of the gradients of the samples at the node
    sum_hessians : float
        The sum of the hessians of the samples at the node
    parent : TreeNode or None, optional(default=None)
        The parent of the node. None for root.

    Attributes
    ----------
    depth : int
        The depth of the node, i.e. its distance from the root
    sample_indices : array of int
        The indices of the samples at the node
    sum_gradients : float
        The sum of the gradients of the samples at the node
    sum_hessians : float
        The sum of the hessians of the samples at the node
    parent : TreeNode or None, optional(default=None)
        The parent of the node. None for root.
    split_info : SplitInfo or None
        The result of the split evaluation
    left_child : TreeNode or None
        The left child of the node. None for leaves.
    right_child : TreeNode or None
        The right child of the node. None for leaves.
    value : float or None
        The value of the leaf, as computed in finalize_leaf(). None for
        non-leaf nodes
    find_split_time : float
        The total time spent computing the histogram and finding the best
        split at the node.
    apply_split_time : float
        The total time spent actually splitting the node, e.g. splitting
        sample_indices into left and right child.
    hist_subtraction : bool
        Wheter the subtraction method was used for computing the histograms.
    partition_start : int
        start position of the node's sample_indices in splitter.partition
    partition_stop : int
        stop position of the node's sample_indices in splitter.partition
    """

    split_info = None
    left_child = None
    right_child = None
    value = None
    histograms = None
    sibling = None
    parent = None
    find_split_time = 0.
    apply_split_time = 0.
    hist_subtraction = False

    # start and stop indices of the node in the splitter.partition
    # array. Concretely,
    # self.sample_indices = view(self.splitter.partition[start:stop])
    # Only used in _update_raw_prediction, because we need to iterate over the
    # leaves and I don't know how to efficiently store the sample_indices
    # views because they're all of different sizes.
    partition_start = 0
    partition_stop = 0

    def __init__(self, depth, sample_indices, sum_gradients,
                 sum_hessians, parent=None):
        self.depth = depth
        self.sample_indices = sample_indices
        self.n_samples = sample_indices.shape[0]
        self.sum_gradients = sum_gradients
        self.sum_hessians = sum_hessians
        self.parent = parent

    def __repr__(self):
        # To help with debugging
        out = "TreeNode: depth={}, ".format(self.depth)
        out += "samples={}".format(len(self.sample_indices))
        if self.split_info is not None:
            out += ", feature_idx={}".format(self.split_info.feature_idx)
            out += ", bin_idx={}".format(self.split_info.bin_idx)
        return out

    def __lt__(self, other_node):
        """Comparison for priority queue.

        Nodes with high gain are higher priority than nodes with low gain.

        heapq.heappush only need the '<' operator.
        heapq.heappop take the smallest item first (smaller is higher
        priority).

        Parameters
        -----------
        other_node : TreeNode
            The node to compare with.
        """
        if self.split_info is None or other_node.split_info is None:
            raise ValueError("Cannot compare nodes with split_info")
        return self.split_info.gain > other_node.split_info.gain


class TreeGrower:
    """Tree grower class used to build a tree.

    The tree is fitted to predict the values of a Newton-Raphson step. The
    splits are considered in a best-first fashion, and the quality of a
    split is defined in splitting._split_gain.

    Parameters
    ----------
    X_binned : array-like of int, shape=(n_samples, n_features)
        The binned input samples. Must be Fortran-aligned.
    gradients : array-like, shape=(n_samples,)
        The gradients of each training sample. Those are the gradients of the
        loss w.r.t the predictions, evaluated at iteration ``i - 1``.
    hessians : array-like, shape=(n_samples,)
        The hessians of each training sample. Those are the hessians of the
        loss w.r.t the predictions, evaluated at iteration ``i - 1``.
    max_leaf_nodes : int or None, optional(default=None)
        The maximum number of leaves for each tree. If None, there is no
        maximum limit.
    max_depth : int or None, optional(default=None)
        The maximum depth of each tree. The depth of a tree is the number of
        nodes to go from the root to the deepest leaf.
    min_samples_leaf : int, optional(default=20)
        The minimum number of samples per leaf.
    min_gain_to_split : float, optional(default=0.)
        The minimum gain needed to split a node. Splits with lower gain will
        be ignored.
    max_bins : int, optional(default=256)
        The maximum number of bins. Used to define the shape of the
        histograms.
    n_bins_per_feature : array-like of int or int, optional(default=None)
        The actual number of bins needed for each feature, which is lower or
        equal to ``max_bins``. If it's an int, all features are considered to
        have the same number of bins. If None, all features are considered to
        have ``max_bins`` bins.
    l2_regularization : float, optional(default=0)
        The L2 regularization parameter.
    min_hessian_to_split : float, optional(default=1e-3)
        The minimum sum of hessians needed in each node. Splits that result in
        at least one child having a sum of hessians less than
        min_hessian_to_split are discarded.
    shrinkage : float, optional(default=1)
        The shrinkage parameter to apply to the leaves values, also known as
        learning rate.
    """
    def __init__(self, X_binned, gradients, hessians, max_leaf_nodes=None,
                 max_depth=None, min_samples_leaf=20, min_gain_to_split=0.,
                 max_bins=256, n_bins_per_feature=None, l2_regularization=0.,
                 min_hessian_to_split=1e-3, shrinkage=1.):

        self._validate_parameters(X_binned, max_leaf_nodes, max_depth,
                                  min_samples_leaf, min_gain_to_split,
                                  l2_regularization, min_hessian_to_split)

        if n_bins_per_feature is None:
            n_bins_per_feature = max_bins

        if isinstance(n_bins_per_feature, int):
            n_bins_per_feature = np.array(
                [n_bins_per_feature] * X_binned.shape[1],
                dtype=np.uint32)

        self.splitter = Splitter(
            X_binned, max_bins, n_bins_per_feature, gradients,
            hessians, l2_regularization, min_hessian_to_split,
            min_samples_leaf, min_gain_to_split)
        self.max_leaf_nodes = max_leaf_nodes
        self.max_bins = max_bins
        self.n_features = X_binned.shape[1]
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.X_binned = X_binned
        self.min_gain_to_split = min_gain_to_split
        self.shrinkage = shrinkage
        self.splittable_nodes = []
        self.finalized_leaves = []
        self.total_find_split_time = 0.  # time spent finding the best splits
        self.total_apply_split_time = 0.  # time spent splitting nodes
        self._intilialize_root()
        self.n_nodes = 1

    def _validate_parameters(self, X_binned, max_leaf_nodes, max_depth,
                             min_samples_leaf, min_gain_to_split,
                             l2_regularization, min_hessian_to_split):
        """Validate parameters passed to __init__.

        Also validate parameters passed to splitter.
        """
        if X_binned.dtype != np.uint8:
            raise NotImplementedError(
                "Explicit feature binning required for now")
        if not X_binned.flags.f_contiguous:
            raise ValueError(
                "X_binned should be passed as Fortran contiguous "
                "array for maximum efficiency.")
        if max_leaf_nodes is not None and max_leaf_nodes < 1:
            raise ValueError('max_leaf_nodes={} should not be'
                             ' smaller than 1'.format(max_leaf_nodes))
        if max_depth is not None and max_depth < 1:
            raise ValueError('max_depth={} should not be'
                             ' smaller than 1'.format(max_depth))
        if min_samples_leaf < 1:
            raise ValueError('min_samples_leaf={} should '
                             'not be smaller than 1'.format(min_samples_leaf))
        if min_gain_to_split < 0:
            raise ValueError('min_gain_to_split={} '
                             'must be positive.'.format(min_gain_to_split))
        if l2_regularization < 0:
            raise ValueError('l2_regularization={} must be '
                             'positive.'.format(l2_regularization))
        if min_hessian_to_split < 0:
            raise ValueError('min_hessian_to_split={} '
                             'must be positive.'.format(min_hessian_to_split))

    def grow(self):
        """Grow the tree, from root to leaves."""
        while self.can_split_further():
            self.split_next()

    def _intilialize_root(self):
        """Initialize root node and finalize it if needed."""
        n_samples = self.X_binned.shape[0]
        depth = 0
        # sum_gradients = np.sum(self.splitter.gradients)
        sum_gradients = sum_parallel(self.splitter.gradients)
        if self.splitter.hessians_are_constant:
            sum_hessians = self.splitter.hessians[0] * n_samples
        else:
            sum_hessians = np.sum(self.splitter.hessians)
        self.root = TreeNode(
            depth=depth,
            sample_indices=self.splitter.partition,
            sum_gradients=sum_gradients,
            sum_hessians=sum_hessians
        )

        self.root.partition_start = 0
        self.root.partition_stop = n_samples

        if (self.max_leaf_nodes is not None and self.max_leaf_nodes == 1):
            self._finalize_leaf(self.root)
            return
        if self.root.n_samples < 2 * self.min_samples_leaf:
            # Do not even bother computing any splitting statistics.
            self._finalize_leaf(self.root)
            return
        if sum_hessians < self.splitter.min_hessian_to_split:
            self._finalize_leaf(self.root)
            return

        self._compute_spittability(self.root)

    def _compute_spittability(self, node, only_hist=False):
        """Compute histograms and best possible split of a node.

        If the best possible gain is 0 or if the constraints aren't met
        (min_samples_leaf, min_hessian_to_split, min_gain_to_split) then the
        node is finalized (transformed into a leaf), else it is pushed on
        the splittable node heap.

        Parameters
        ----------
        node : TreeNode
            The node to evaluate.
        only_hist : bool, optional (default=False)
            Whether to only compute the histograms and the SplitInfo. It is
            set to ``True`` when ``_compute_spittability`` was called by a
            sibling node: we only want to compute the histograms (which also
            computes the ``SplitInfo``), not finalize or push the node. If
            ``_compute_spittability`` is called again by the grower on this
            same node, the histograms won't be computed again.
        """
        # Compute split_info and histograms if not already done
        if node.split_info is None and node.histograms is None:
            # If the sibling has less samples, compute its hist first (with
            # the regular method) and use the subtraction method for the
            # current node
            if node.sibling is not None:  # root has no sibling
                if node.sibling.n_samples < node.n_samples:
                    self._compute_spittability(node.sibling, only_hist=True)
                    # As hist of sibling is now computed we'll use the hist
                    # subtraction method for the current node.
                    node.hist_subtraction = True

            tic = time()
            histograms = np.zeros(shape=(self.n_features, self.max_bins),
                                  dtype=HISTOGRAM_DTYPE)
            if node.hist_subtraction:
                if node is node.parent.right_child:
                    sum_gradients = node.parent.split_info.sum_gradient_right
                    sum_hessians = node.parent.split_info.sum_hessian_right
                else:
                    sum_gradients = node.parent.split_info.sum_gradient_left
                    sum_hessians = node.parent.split_info.sum_hessian_left
                split_info = self.splitter.find_node_split_subtraction(
                    node.sample_indices,
                    sum_gradients, sum_hessians, node.parent.histograms,
                    node.sibling.histograms, histograms)
            else:
                split_info = self.splitter.find_node_split(
                    node.sample_indices, histograms)
            toc = time()
            node.find_split_time = toc - tic
            self.total_find_split_time += node.find_split_time
            node.split_info = split_info
            node.histograms = histograms

        if only_hist:
            # _compute_spittability was called by a sibling. We only needed to
            # compute the histogram.
            return

        if node.split_info.gain <= 0:  # no valid split
            # Note: this condition is reached if either all the leaves are
            # pure (best gain = 0), or if no split would satisfy the
            # constraints, (min_hessians_to_split, min_gain_to_split,
            # min_samples_leaf)
            self._finalize_leaf(node)

        else:
            heappush(self.splittable_nodes, node)

    def split_next(self):
        """Split the node with highest potential gain.

        Returns
        -------
        left : TreeNode
            The resulting left child.
        right : TreeNode
            The resulting right child.
        """
        if len(self.splittable_nodes) == 0:
            raise StopIteration("No more splittable nodes")

        # Consider the node with the highest loss reduction (a.k.a. gain)
        node = heappop(self.splittable_nodes)

        tic = time()
        (sample_indices_left,
         sample_indices_right,
         right_child_pos) = self.splitter.split_indices(node.split_info,
                                                        node.sample_indices)
        toc = time()
        node.apply_split_time = toc - tic
        self.total_apply_split_time += node.apply_split_time

        depth = node.depth + 1
        n_leaf_nodes = len(self.finalized_leaves) + len(self.splittable_nodes)
        n_leaf_nodes += 2

        left_child_node = TreeNode(depth,
                                   sample_indices_left,
                                   node.split_info.sum_gradient_left,
                                   node.split_info.sum_hessian_left,
                                   parent=node)
        right_child_node = TreeNode(depth,
                                    sample_indices_right,
                                    node.split_info.sum_gradient_right,
                                    node.split_info.sum_hessian_right,
                                    parent=node)
        left_child_node.sibling = right_child_node
        right_child_node.sibling = left_child_node
        node.right_child = right_child_node
        node.left_child = left_child_node

        # set start and stop indices
        left_child_node.partition_start = node.partition_start
        left_child_node.partition_stop = node.partition_start + right_child_pos
        right_child_node.partition_start = left_child_node.partition_stop
        right_child_node.partition_stop = node.partition_stop

        self.n_nodes += 2

        if self.max_depth is not None and depth == self.max_depth:
            self._finalize_leaf(left_child_node)
            self._finalize_leaf(right_child_node)
            return left_child_node, right_child_node

        if (self.max_leaf_nodes is not None
                and n_leaf_nodes == self.max_leaf_nodes):
            self._finalize_leaf(left_child_node)
            self._finalize_leaf(right_child_node)
            self._finalize_splittable_nodes()
            return left_child_node, right_child_node

        if left_child_node.n_samples < self.min_samples_leaf * 2:
            self._finalize_leaf(left_child_node)
        else:
            self._compute_spittability(left_child_node)

        if right_child_node.n_samples < self.min_samples_leaf * 2:
            self._finalize_leaf(right_child_node)
        else:
            self._compute_spittability(right_child_node)

        return left_child_node, right_child_node

    def can_split_further(self):
        """Return True if there are still nodes to split."""
        return len(self.splittable_nodes) >= 1

    def _finalize_leaf(self, node):
        """Compute the prediction value that minimizes the objective function.

        This sets the node.value attribute (node is a leaf iff node.value is
        not None).

        See Equation 5 of:
        XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
        https://arxiv.org/abs/1603.02754
        """
        node.value = -self.shrinkage * node.sum_gradients / (
            node.sum_hessians + self.splitter.l2_regularization)
        self.finalized_leaves.append(node)

    def _finalize_splittable_nodes(self):
        """Transform all splittable nodes into leaves.

        Used when some constraint is met e.g. maximum number of leaves or
        maximum depth."""
        while len(self.splittable_nodes) > 0:
            node = self.splittable_nodes.pop()
            self._finalize_leaf(node)

    def make_predictor(self, bin_thresholds=None):
        """Make a TreePredictor object out of the current tree.

        Parameters
        ----------
        bin_thresholds : array-like of floats, optional (default=None)
            The actual thresholds values of each bin.

        Returns
        -------
        A TreePredictor object.
        """
        predictor_nodes = np.zeros(self.n_nodes, dtype=PREDICTOR_RECORD_DTYPE)
        _fill_predictor_node_array(predictor_nodes, self.root,
                                   bin_thresholds=bin_thresholds)
        return TreePredictor(predictor_nodes)


def _fill_predictor_node_array(predictor_nodes, grower_node,
                               bin_thresholds=None, next_free_idx=0):
    """Helper used in make_predictor to set the TreePredictor fields."""
    node = predictor_nodes[next_free_idx]
    node['count'] = grower_node.n_samples
    node['depth'] = grower_node.depth
    if grower_node.split_info is not None:
        node['gain'] = grower_node.split_info.gain
    else:
        node['gain'] = -1

    if grower_node.value is not None:
        # Leaf node
        node['is_leaf'] = True
        node['value'] = grower_node.value
        return next_free_idx + 1
    else:
        # Decision node
        split_info = grower_node.split_info
        feature_idx, bin_idx = split_info.feature_idx, split_info.bin_idx
        node['feature_idx'] = feature_idx
        node['bin_threshold'] = bin_idx
        if bin_thresholds is not None:
            threshold = bin_thresholds[feature_idx][bin_idx]
            node['threshold'] = threshold
        next_free_idx += 1

        node['left'] = next_free_idx
        next_free_idx = _fill_predictor_node_array(
            predictor_nodes, grower_node.left_child,
            bin_thresholds=bin_thresholds, next_free_idx=next_free_idx)

        node['right'] = next_free_idx
        return _fill_predictor_node_array(
            predictor_nodes, grower_node.right_child,
            bin_thresholds=bin_thresholds, next_free_idx=next_free_idx)
