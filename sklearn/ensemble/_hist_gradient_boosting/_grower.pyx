# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
import numpy as np
cimport numpy as np

from .common cimport node_struct
from ._predictor cimport TreePredictor

np.import_array()

def _fill_predictor_node_array(TreePredictor predictor_nodes,
                               grower_node,  # TreeNode
                               list bin_thresholds,
                               const np.npy_uint32[:] n_bins_non_missing,
                               int next_free_idx=0):
    """Helper used in make_predictor to set the TreePredictor fields."""
    cdef:
        node_struct * node = predictor_nodes.get(next_free_idx)
        int feature_idx
        int bin_idx


    node.count = grower_node.n_samples
    node.depth = grower_node.depth
    if grower_node.split_info is not None:
        node.gain = grower_node.split_info.gain
    else:
        node.gain = -1

    node.value = grower_node.value

    if grower_node.is_leaf:
        # Leaf node
        node.is_leaf = True
        return next_free_idx + 1
    else:
        # Decision node
        split_info = grower_node.split_info
        feature_idx, bin_idx = split_info.feature_idx, split_info.bin_idx
        node.feature_idx = feature_idx
        node.bin_threshold = bin_idx
        node.missing_go_to_left = split_info.missing_go_to_left

        if split_info.bin_idx == n_bins_non_missing[feature_idx] - 1:
            # Split is on the last non-missing bin: it's a "split on nans". All
            # nans go to the right, the rest go to the left.
            node.threshold = np.inf
        elif bin_thresholds is not None:
            node.threshold = bin_thresholds[feature_idx][bin_idx]

        next_free_idx += 1
        node.left = next_free_idx
        next_free_idx = _fill_predictor_node_array(
            predictor_nodes, grower_node.left_child,
            bin_thresholds=bin_thresholds,
            n_bins_non_missing=n_bins_non_missing,
            next_free_idx=next_free_idx)

        node.right = next_free_idx
        return _fill_predictor_node_array(
            predictor_nodes, grower_node.right_child,
            bin_thresholds=bin_thresholds,
            n_bins_non_missing=n_bins_non_missing,
            next_free_idx=next_free_idx)
